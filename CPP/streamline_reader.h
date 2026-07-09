#ifndef NPSAT_URF_STREAMLINE_READER_H
#define NPSAT_URF_STREAMLINE_READER_H

#include "my_structures.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <istream>
#include <sstream>
#include <stdexcept>
#include <string>

inline std::uint64_t streamlineId(const double value) {
    return static_cast<std::uint64_t>(std::llround(value));
}

inline int parseEndReasonToken(const std::string& token) {
    char* endPtr = NULL;
    const long value = std::strtol(token.c_str(), &endPtr, 10);
    if (endPtr != token.c_str() && *endPtr == '\0') {
        return static_cast<int>(value);
    }
    return parseExitReason(token);
}

inline void writeDiscardedStreamline(std::ofstream& discardFile,
                                     const StreamlineTrajectory& trajectory,
                                     const std::string& reason) {
    discardFile << static_cast<unsigned long long>(trajectory.Eid) << ", "
                << static_cast<unsigned long long>(trajectory.Sid) << ", "
                << trajectory.samples.size() << ", "
                << (trajectory.has_termination ? 1 : 0) << ", "
                << static_cast<unsigned long long>(trajectory.termination_pid) << ", "
                << trajectory.end_reason << ", "
                << reason << std::endl;
}

inline void appendSampleOrDiscard(StreamlineTrajectory& trajectory,
                                  const StreamlinePoint& sample,
                                  std::ofstream& discardFile) {
    if (!trajectory.samples.empty() &&
        (trajectory.Eid != sample.Eid || trajectory.Sid != sample.Sid)) {
        writeDiscardedStreamline(discardFile, trajectory, "eid_sid_changed_before_termination");
        trajectory.clear();
    }

    if (trajectory.samples.empty()) {
        trajectory.Eid = sample.Eid;
        trajectory.Sid = sample.Sid;
    }

    trajectory.samples.push_back(sample);
}

inline bool completeTrajectory(StreamlineTrajectory& trajectory,
                               const std::uint64_t terminationPid,
                               const std::uint64_t Eid,
                               const std::uint64_t Sid,
                               const int endReason,
                               std::ofstream& discardFile) {
    if (trajectory.samples.empty()) {
        trajectory.Eid = Eid;
        trajectory.Sid = Sid;
        trajectory.termination_pid = terminationPid;
        trajectory.end_reason = endReason;
        trajectory.has_termination = true;
        writeDiscardedStreamline(discardFile, trajectory, "termination_without_samples");
        trajectory.clear();
        return false;
    }

    if (trajectory.Eid != Eid || trajectory.Sid != Sid) {
        writeDiscardedStreamline(discardFile, trajectory, "termination_eid_sid_mismatch");
        trajectory.clear();
        return false;
    }

    trajectory.termination_pid = terminationPid;
    trajectory.end_reason = endReason;
    trajectory.has_termination = true;
    return true;
}

inline bool readNextAsciiStreamline(std::istream& in,
                                    StreamlineTrajectory& trajectory,
                                    std::ofstream& discardFile) {
    trajectory.clear();

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        std::istringstream inp(line);
        double marker = 0.0;
        if (!(inp >> marker)) {
            continue;
        }

        if (marker == -1.0) {
            double terminationPid = 0.0;
            double Eid = 0.0;
            double Sid = 0.0;
            std::string endReason;
            if (!(inp >> terminationPid >> Eid >> Sid >> endReason)) {
                writeDiscardedStreamline(discardFile, trajectory, "malformed_termination_record");
                trajectory.clear();
                continue;
            }
            if (completeTrajectory(trajectory, streamlineId(terminationPid), streamlineId(Eid),
                                   streamlineId(Sid), parseEndReasonToken(endReason),
                                   discardFile)) {
                return true;
            }
            continue;
        }

        if (marker == -9.0) {
            double Eid = 0.0;
            double Sid = 0.0;
            std::string endReason;
            if (!(inp >> Eid >> Sid >> endReason)) {
                writeDiscardedStreamline(discardFile, trajectory, "malformed_legacy_termination_record");
                trajectory.clear();
                continue;
            }

            const std::uint64_t terminationPid = trajectory.samples.empty()
                                                 ? 0
                                                 : trajectory.samples.back().pid;
            if (completeTrajectory(trajectory, terminationPid, streamlineId(Eid),
                                   streamlineId(Sid), parseExitReason(endReason),
                                   discardFile)) {
                return true;
            }
            continue;
        }

        StreamlinePoint sample;
        double Eid = 0.0;
        double Sid = 0.0;
        sample.pid = streamlineId(marker);
        if (!(inp >> Eid >> Sid >> sample.x >> sample.y >> sample.z >> sample.vmag)) {
            writeDiscardedStreamline(discardFile, trajectory, "malformed_sample_record");
            trajectory.clear();
            continue;
        }
        sample.Eid = streamlineId(Eid);
        sample.Sid = streamlineId(Sid);
        appendSampleOrDiscard(trajectory, sample, discardFile);
    }

    if (!trajectory.samples.empty()) {
        writeDiscardedStreamline(discardFile, trajectory, "end_of_file_before_termination");
        trajectory.clear();
    }

    return false;
}

inline bool readBinaryStreamlineRow(std::istream& in,
                                    std::array<double, 7>& row,
                                    const std::string& filename) {
    in.read(reinterpret_cast<char *>(row.data()),
            static_cast<std::streamsize>(row.size() * sizeof(double)));

    if (in) {
        return true;
    }

    if (in.eof() && in.gcount() == 0) {
        return false;
    }

    throw std::runtime_error("Incomplete binary streamline record in: " + filename);
}

inline bool readNextBinaryStreamline(std::istream& in,
                                     StreamlineTrajectory& trajectory,
                                     const std::string& filename,
                                     std::ofstream& discardFile) {
    trajectory.clear();

    std::array<double, 7> row;
    while (readBinaryStreamlineRow(in, row, filename)) {
        if (row[0] == -1.0) {
            if (completeTrajectory(trajectory, streamlineId(row[1]), streamlineId(row[2]),
                                   streamlineId(row[3]), static_cast<int>(std::llround(row[4])),
                                   discardFile)) {
                return true;
            }
            continue;
        }

        StreamlinePoint sample;
        sample.pid = streamlineId(row[0]);
        sample.Eid = streamlineId(row[1]);
        sample.Sid = streamlineId(row[2]);
        sample.x = row[3];
        sample.y = row[4];
        sample.z = row[5];
        sample.vmag = row[6];
        appendSampleOrDiscard(trajectory, sample, discardFile);
    }

    if (!trajectory.samples.empty()) {
        writeDiscardedStreamline(discardFile, trajectory, "end_of_file_before_termination");
        trajectory.clear();
    }

    return false;
}

#endif //NPSAT_URF_STREAMLINE_READER_H
