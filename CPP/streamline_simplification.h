#ifndef NPSAT_URF_STREAMLINE_SIMPLIFICATION_H
#define NPSAT_URF_STREAMLINE_SIMPLIFICATION_H

#include "my_structures.h"
#include "streamline_reader.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

struct SimplifiedPoint{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double v = 0.0;
    double age = 0.0;
};

struct SimplifiedVtkRecord{
    unsigned long long Eid = 0;
    unsigned long long Sid = 0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double v = 0.0;
    double age = 0.0;
};

struct SimplifiedVtkLine{
    unsigned long long Eid = 0;
    unsigned long long Sid = 0;
    unsigned long long start = 0;
    unsigned long long count = 0;
};

inline bool parseSimplifiedStreamlineRow(const std::string& line,
                                         SimplifiedVtkRecord& record) {
    std::string row = line;
    std::replace(row.begin(), row.end(), ',', ' ');
    std::istringstream inp(row);
    return static_cast<bool>(inp >> record.Eid >> record.Sid
                                 >> record.x >> record.y >> record.z
                                 >> record.v >> record.age);
}

inline bool writeSimplifiedStreamlinesVtk(const std::string& simplifiedFilename,
                                          const std::string& vtkFilename) {
    std::ifstream simplifiedFile(simplifiedFilename.c_str());
    if (!simplifiedFile.good()) {
        return false;
    }

    std::vector<SimplifiedVtkRecord> records;
    std::vector<SimplifiedVtkLine> lines;
    std::string line;
    bool firstDataRow = true;

    while (std::getline(simplifiedFile, line)) {
        if (line.empty()) {
            continue;
        }

        SimplifiedVtkRecord record;
        if (!parseSimplifiedStreamlineRow(line, record)) {
            if (firstDataRow) {
                firstDataRow = false;
                continue;
            }
            return false;
        }
        firstDataRow = false;

        if (lines.empty() ||
            lines.back().Eid != record.Eid ||
            lines.back().Sid != record.Sid) {
            SimplifiedVtkLine vtkLine;
            vtkLine.Eid = record.Eid;
            vtkLine.Sid = record.Sid;
            vtkLine.start = static_cast<unsigned long long>(records.size());
            vtkLine.count = 0;
            lines.push_back(vtkLine);
        }

        records.push_back(record);
        lines.back().count++;
    }

    std::ofstream vtkFile(vtkFilename.c_str());
    if (!vtkFile.good()) {
        return false;
    }

    unsigned long long lineListSize = 0;
    for (unsigned int i = 0; i < lines.size(); ++i) {
        lineListSize += lines[i].count + 1;
    }

    vtkFile << "# vtk DataFile Version 3.0" << std::endl;
    vtkFile << "NPSAT simplified streamlines" << std::endl;
    vtkFile << "ASCII" << std::endl;
    vtkFile << "DATASET POLYDATA" << std::endl;
    vtkFile << "POINTS " << records.size() << " double" << std::endl;
    vtkFile << std::setprecision(10) << std::fixed;
    for (unsigned int i = 0; i < records.size(); ++i) {
        vtkFile << records[i].x << " "
                << records[i].y << " "
                << records[i].z << std::endl;
    }

    vtkFile << "LINES " << lines.size() << " " << lineListSize << std::endl;
    for (unsigned int i = 0; i < lines.size(); ++i) {
        vtkFile << lines[i].count;
        for (unsigned long long j = 0; j < lines[i].count; ++j) {
            vtkFile << " " << lines[i].start + j;
        }
        vtkFile << std::endl;
    }

    vtkFile << "POINT_DATA " << records.size() << std::endl;
    vtkFile << "SCALARS velocity double 1" << std::endl;
    vtkFile << "LOOKUP_TABLE default" << std::endl;
    for (unsigned int i = 0; i < records.size(); ++i) {
        vtkFile << records[i].v << std::endl;
    }

    vtkFile << "SCALARS age double 1" << std::endl;
    vtkFile << "LOOKUP_TABLE default" << std::endl;
    for (unsigned int i = 0; i < records.size(); ++i) {
        vtkFile << records[i].age << std::endl;
    }

    vtkFile << std::setprecision(0) << std::fixed;
    vtkFile << "SCALARS Eid double 1" << std::endl;
    vtkFile << "LOOKUP_TABLE default" << std::endl;
    for (unsigned int i = 0; i < records.size(); ++i) {
        vtkFile << static_cast<double>(records[i].Eid) << std::endl;
    }

    vtkFile << "SCALARS Sid double 1" << std::endl;
    vtkFile << "LOOKUP_TABLE default" << std::endl;
    for (unsigned int i = 0; i < records.size(); ++i) {
        vtkFile << static_cast<double>(records[i].Sid) << std::endl;
    }

    return true;
}

inline double pointSegmentDistance(const StreamlinePoint& p,
                                   const StreamlinePoint& a,
                                   const StreamlinePoint& b) {
    const double abx = b.x - a.x;
    const double aby = b.y - a.y;
    const double abz = b.z - a.z;
    const double apx = p.x - a.x;
    const double apy = p.y - a.y;
    const double apz = p.z - a.z;
    const double ab2 = abx*abx + aby*aby + abz*abz;

    double t = 0.0;
    if (ab2 > 0.0) {
        t = (apx*abx + apy*aby + apz*abz) / ab2;
        if (t < 0.0) {
            t = 0.0;
        }
        else if (t > 1.0) {
            t = 1.0;
        }
    }

    const double cx = a.x + t*abx;
    const double cy = a.y + t*aby;
    const double cz = a.z + t*abz;
    const double dx = p.x - cx;
    const double dy = p.y - cy;
    const double dz = p.z - cz;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

inline void markDouglasPeuckerPoints(const std::vector<StreamlinePoint>& points,
                                     const unsigned int first,
                                     const unsigned int last,
                                     const double tolerance,
                                     std::vector<int>& keep) {
    if (last <= first + 1) {
        return;
    }

    double maxDistance = -1.0;
    unsigned int splitIndex = first;
    for (unsigned int i = first + 1; i < last; ++i) {
        const double distance = pointSegmentDistance(points[i], points[first], points[last]);
        if (distance > maxDistance) {
            maxDistance = distance;
            splitIndex = i;
        }
    }

    if (maxDistance > tolerance) {
        keep[splitIndex] = 1;
        markDouglasPeuckerPoints(points, first, splitIndex, tolerance, keep);
        markDouglasPeuckerPoints(points, splitIndex, last, tolerance, keep);
    }
}

inline bool calculateDetailedAges(const StreamlineTrajectory& trajectory,
                                  std::vector<double>& ages) {
    if (trajectory.samples.empty()) {
        return false;
    }

    ages.assign(trajectory.samples.size(), 0.0);
    for (unsigned int i = 1; i < trajectory.samples.size(); ++i) {
        const StreamlinePoint& a = trajectory.samples[i - 1];
        const StreamlinePoint& b = trajectory.samples[i];
        const double dx = b.x - a.x;
        const double dy = b.y - a.y;
        const double dz = b.z - a.z;
        const double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (b.vmag <= 0.0) {
            return false;
        }
        ages[i] = ages[i - 1] + distance / b.vmag;
    }

    return true;
}

inline void writeSimplifiedStreamline(const StreamlineTrajectory& trajectory,
                                      const URFoptions& opt,
                                      std::ofstream& simplifiedFile,
                                      std::ofstream& discardFile) {
    if (!opt.simplifyStreamline) {
        return;
    }

    if (trajectory.samples.size() < 2) {
        writeDiscardedStreamline(discardFile, trajectory, "simplification_not_enough_samples");
        return;
    }

    std::vector<double> detailedAges;
    if (!calculateDetailedAges(trajectory, detailedAges)) {
        writeDiscardedStreamline(discardFile, trajectory, "simplification_invalid_velocity");
        return;
    }

    std::vector<int> keep(trajectory.samples.size(), 0);
    keep.front() = 1;
    keep.back() = 1;
    markDouglasPeuckerPoints(trajectory.samples, 0,
                             static_cast<unsigned int>(trajectory.samples.size() - 1),
                             opt.simplifyTolerance, keep);

    std::vector<unsigned int> indices;
    for (unsigned int i = 0; i < keep.size(); ++i) {
        if (keep[i]) {
            indices.push_back(i);
        }
    }

    std::vector<SimplifiedPoint> simplified(indices.size());
    for (unsigned int i = 0; i < indices.size(); ++i) {
        const StreamlinePoint& point = trajectory.samples[indices[i]];
        simplified[i].x = point.x;
        simplified[i].y = point.y;
        simplified[i].z = point.z;
        simplified[i].age = detailedAges[indices[i]];
    }

    for (unsigned int i = 0; i + 1 < simplified.size(); ++i) {
        const double dx = simplified[i + 1].x - simplified[i].x;
        const double dy = simplified[i + 1].y - simplified[i].y;
        const double dz = simplified[i + 1].z - simplified[i].z;
        const double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
        const double dt = simplified[i + 1].age - simplified[i].age;
        simplified[i].v = dt > 0.0 ? distance / dt : 0.0;
    }
    simplified.back().v = simplified.size() > 1 ? simplified[simplified.size() - 2].v : 0.0;

    for (unsigned int i = 0; i < simplified.size(); ++i) {
        simplifiedFile << static_cast<unsigned long long>(trajectory.Eid) << ","
                       << static_cast<unsigned long long>(trajectory.Sid) << ","
                       << std::setprecision(10) << std::fixed
                       << simplified[i].x << ","
                       << simplified[i].y << ","
                       << simplified[i].z << ","
                       << simplified[i].v << ","
                       << simplified[i].age << std::endl;
    }
}

#endif //NPSAT_URF_STREAMLINE_SIMPLIFICATION_H
