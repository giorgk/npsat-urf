#include <chrono>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <mpi.h>

#include "my_structures.h"
#include "NPSAT_URF_main.h"
#include "streamline_reader.h"
#include "streamline_simplification.h"

void writeOutputHeader(std::ofstream& ofile, const URFoptions& opt) {
    ofile << "Eid, Sid, ER, p_cdsX, p_cdsY, p_cdsZ, v_cds, p_lndX, p_lndY, Len";
    for (int i = opt.por.startValue; i <= opt.por.endValue; i = i + opt.por.interval) {
        ofile << ", Age" << i << ", mean" << i << ", std" << i << ", err" << i;
        if (opt.calcDecay) {
            ofile << ", meanDc" << i << ", stdDc" << i << ", ScaleDc" << i << ", errDc" << i;
        }
        if (opt.calcDiff) {
            ofile << ", meanDf" << i << ", stdDf" << i << ", ScaleDf" << i << ", errDf" << i;
        }
        if (i > 200) {
            break;
        }
    }
    ofile << std::endl;
}

bool shouldFitEndReason(const int endReason, const URFoptions& opt) {
    for (unsigned int i = 0; i < opt.er_to_run.size(); ++i) {
        if (opt.er_to_run[i] < 0 || opt.er_to_run[i] == endReason) {
            return true;
        }
    }
    return false;
}

std::string buildArrayInputFilename(const URFoptions& opt,
                                    const int rankId,
                                    const int iterId) {
    return opt.prefixInput +
           "rank_" +
           num2Padstr(rankId, opt.paddingZeros) +
           "_iter_" +
           num2Padstr(iterId, opt.iterPaddingZeros) +
           "." + opt.suffixInput;
}

bool buildStreamlineSegments(const StreamlineTrajectory& trajectory,
                             const URFoptions& opt,
                             std::vector<segInfo>& strmlnSeg,
                             double& streamlineLength,
                             double& pCdsX,
                             double& pCdsY,
                             double& pCdsZ,
                             double& pLndX,
                             double& pLndY) {
    if (trajectory.samples.size() < 2) {
        return false;
    }

    strmlnSeg.clear();
    streamlineLength = 0.0;
    double leftOverLen = 0.0;

    pCdsX = trajectory.samples.front().x;
    pCdsY = trajectory.samples.front().y;
    pCdsZ = trajectory.samples.front().z;
    pLndX = trajectory.samples.back().x;
    pLndY = trajectory.samples.back().y;

    for (unsigned int i = 1; i < trajectory.samples.size(); ++i) {
        const StreamlinePoint& a = trajectory.samples[i - 1];
        const StreamlinePoint& b = trajectory.samples[i];
        const double dx = b.x - a.x;
        const double dy = b.y - a.y;
        const double dz = b.z - a.z;
        const double segLen = std::sqrt(dx*dx + dy*dy + dz*dz) + leftOverLen;

        if (segLen < opt.minElemSize) {
            leftOverLen = leftOverLen + segLen;
        }
        else if (segLen > opt.maxElemSize) {
            const double nSplits = std::ceil(segLen / opt.maxElemSize);
            const int inSplits = static_cast<int>(nSplits);
            const double lenSplit = segLen / nSplits;
            for (int j = 0; j < inSplits; ++j) {
                strmlnSeg.emplace_back(b.vmag, lenSplit);
                streamlineLength = streamlineLength + lenSplit;
            }
            leftOverLen = 0.0;
        }
        else {
            strmlnSeg.emplace_back(b.vmag, segLen);
            streamlineLength = streamlineLength + segLen;
            leftOverLen = 0.0;
        }
    }

    if (!strmlnSeg.empty()) {
        strmlnSeg[strmlnSeg.size() - 1].l = strmlnSeg[strmlnSeg.size() - 1].l + leftOverLen;
        streamlineLength = streamlineLength + leftOverLen;
    }

    return !strmlnSeg.empty();
}

bool processCompleteStreamline(const StreamlineTrajectory& trajectory,
                               URFoptions& opt,
                               std::ofstream& ofile,
                               std::ofstream& discardFile,
                               int& cntStrml,
                               std::chrono::steady_clock::time_point& beginTime) {
    std::vector<segInfo> strmlnSeg;
    double streamlineLength = 0.0;
    double pLndX = 0.0;
    double pLndY = 0.0;
    double pCdsX = 0.0;
    double pCdsY = 0.0;
    double pCdsZ = 0.0;

    if (!buildStreamlineSegments(trajectory, opt, strmlnSeg, streamlineLength,
                                 pCdsX, pCdsY, pCdsZ, pLndX, pLndY)) {
        writeDiscardedStreamline(discardFile, trajectory, "not_enough_samples");
        return false;
    }

    const std::chrono::steady_clock::time_point endTime1 = std::chrono::steady_clock::now();
    std::cout << static_cast<unsigned long long>(trajectory.Eid) << " "
              << static_cast<unsigned long long>(trajectory.Sid) << " "
              << ++cntStrml << " ["
              << std::chrono::duration_cast<std::chrono::microseconds>(endTime1 - beginTime).count()/1000000.0;

    FittedParam fp;
    ofile << static_cast<unsigned long long>(trajectory.Eid) << ", "
          << static_cast<unsigned long long>(trajectory.Sid) << ", "
          << trajectory.end_reason << ", "
          << std::setprecision(2) << std::fixed
          << pCdsX << ", " << pCdsY << ", " << pCdsZ << ", "
          << std::setprecision(5) << strmlnSeg[0].v << ", "
          << std::setprecision(2) << pLndX << ", " << pLndY << ", "
          << streamlineLength << ", ";

    for (int i = opt.por.startValue; i <= opt.por.endValue; i = i + opt.por.interval) {
        const double velMult = static_cast<double>(i) / 10.0;
        fp.reset();
        if (shouldFitEndReason(trajectory.end_reason, opt)) {
            NPSATurf(strmlnSeg, streamlineLength, velMult, opt, fp);
        }
        else {
            fp.setVal(0.0);
        }

        ofile << std::setprecision(2) << std::fixed << fp.Age
              << std::setprecision(6) << std::scientific
              << ", " << fp.urf.m << ", " << fp.urf.s << ", " << fp.urf.err;
        if (opt.calcDecay) {
            ofile << ", " << fp.Decay.m << ", " << fp.Decay.s << ", "
                  << fp.Decay.sc << ", " << fp.Decay.err;
        }
        if (opt.calcDiff) {
            ofile << ", " << fp.Diff.m << ", " << fp.Diff.s << ", "
                  << fp.Diff.sc << ", " << fp.Diff.err;
        }
        if (i != opt.por.endValue) {
            ofile << ", ";
        }
    }
    ofile << std::endl;

    const std::chrono::steady_clock::time_point endTime2 = std::chrono::steady_clock::now();
    std::cout << ", "
              << std::chrono::duration_cast<std::chrono::microseconds>(endTime2 - beginTime).count()/1000000.0
              << "]" << std::endl;

    beginTime = std::chrono::steady_clock::now();
    return true;
}

struct OutputFiles {
    std::ofstream result;
    std::ofstream discard;
    std::ofstream simplified;
    std::string simplifiedFilename;
    std::string simplifiedVtkFilename;
};

std::string buildOutputBase(const std::string& prefix,
                            const int mpiRank,
                            const int partId) {
    std::string base = prefix + "_rank_" + std::to_string(mpiRank);
    if (partId > 0) {
        base += "_" + std::to_string(partId);
    }
    return base;
}

bool openOutputFiles(const URFoptions& opt,
                     const int mpiRank,
                     const int partId,
                     OutputFiles& files) {
    const std::string resultFilename = buildOutputBase(opt.prefixOutput, mpiRank, partId) + ".dat";
    const std::string discardFilename = buildOutputBase(opt.prefixDiscard, mpiRank, partId) + ".dat";
    files.simplifiedFilename = buildOutputBase(opt.prefixSimplified, mpiRank, partId) + ".dat";
    files.simplifiedVtkFilename = buildOutputBase(opt.prefixSimplified, mpiRank, partId) + ".vtk";

    files.result.open(resultFilename.c_str());
    if (!files.result.good()) {
        std::cerr << "MPI rank " << mpiRank << ": can't open output file "
                  << resultFilename << std::endl;
        return false;
    }
    writeOutputHeader(files.result, opt);

    files.discard.open(discardFilename.c_str());
    if (!files.discard.good()) {
        std::cerr << "MPI rank " << mpiRank << ": can't open discard file "
                  << discardFilename << std::endl;
        return false;
    }
    files.discard << "Eid, Sid, samples, has_termination, termination_pid, ER, reason" << std::endl;

    if (opt.simplifyStreamline) {
        files.simplified.open(files.simplifiedFilename.c_str());
        if (!files.simplified.good()) {
            std::cerr << "MPI rank " << mpiRank << ": can't open simplified file "
                      << files.simplifiedFilename << std::endl;
            return false;
        }
        files.simplified << "Eid,Sid,x,y,z,v,a" << std::endl;
    }
    return true;
}

bool closeOutputFiles(const URFoptions& opt,
                      const int mpiRank,
                      OutputFiles& files) {
    files.result.close();
    files.discard.close();
    if (opt.simplifyStreamline) {
        files.simplified.close();
    }
    if (opt.writeSimplifiedVtk &&
        !writeSimplifiedStreamlinesVtk(files.simplifiedFilename,
                                       files.simplifiedVtkFilename)) {
        std::cerr << "MPI rank " << mpiRank << ": can't write simplified VTK file "
                  << files.simplifiedVtkFilename << std::endl;
        return false;
    }
    return true;
}

bool processInputFile(const int workId,
                      URFoptions& opt,
                      OutputFiles& files,
                      const int mpiRank) {
    const int inputRankId = workId % opt.nproc;
    const int iterId = workId / opt.nproc;
    const std::string filename = buildArrayInputFilename(opt, inputRankId, iterId);

    std::cout << "MPI rank " << mpiRank << " reading work " << workId
              << " (input rank " << inputRankId << ", iter " << iterId
              << "): " << filename << std::endl;

    std::ifstream datafile;
    if (opt.fileType.compare("npsat_bin") == 0) {
        datafile.open(filename.c_str(), std::ios::binary);
    }
    else {
        datafile.open(filename.c_str());
    }
    if (!datafile.good()) {
        std::cerr << "MPI rank " << mpiRank << ": can't open input file "
                  << filename << std::endl;
        return false;
    }

    std::chrono::steady_clock::time_point beginTime = std::chrono::steady_clock::now();
    int cntStrml = 0;
    StreamlineTrajectory trajectory;
    try {
        while (true) {
            bool found = false;
            if (opt.fileType.compare("npsat_bin") == 0) {
                found = readNextBinaryStreamline(datafile, trajectory, filename, files.discard);
            }
            else {
                found = readNextAsciiStreamline(datafile, trajectory, files.discard);
            }
            if (!found) {
                break;
            }
            writeSimplifiedStreamline(trajectory, opt, files.simplified, files.discard);
            processCompleteStreamline(trajectory, opt, files.result, files.discard,
                                      cntStrml, beginTime);
        }
    }
    catch (const std::exception& e) {
        std::cerr << "MPI rank " << mpiRank << ", work " << workId
                  << ": " << e.what() << std::endl;
        return false;
    }
    return true;
}

void reportProgress(const int completed,
                    const int total,
                    const int interval,
                    int& nextPercent) {
    const int percent = static_cast<int>((100LL * completed) / total);
    while (nextPercent <= 100 && percent >= nextPercent) {
        std::cout << "Progress: " << nextPercent << "% (" << completed
                  << "/" << total << " input files completed)" << std::endl;
        nextPercent += interval;
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int mpiRank = 0;
    int mpiSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

    if (argc != 2) {
        if (mpiRank == 0) {
            std::cout << "Usage: NPSAT_URF <options_file> or NPSAT_URF -v" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    const std::string inputArg(argv[1]);
    if (inputArg.compare("-v") == 0) {
        if (mpiRank == 0) {
            std::cout << "version 2.0.0" << std::endl;
        }
        MPI_Finalize();
        return 0;
    }

    URFoptions opt;
    int localValid = readOptionFile(inputArg, opt) ? 1 : 0;
    if (localValid && (opt.nproc <= 0 || opt.niter <= 0 ||
                       opt.outputFilesPerPart < 0 ||
                       opt.progressPercent < 1 || opt.progressPercent > 100 ||
                       opt.nproc > std::numeric_limits<int>::max() / opt.niter)) {
        if (mpiRank == 0) {
            std::cerr << "Invalid options: nproc and niter must be positive, "
                      << "output_files_per_part must be nonnegative, and "
                      << "progress_percent must be between 1 and 100." << std::endl;
        }
        localValid = 0;
    }
    if (localValid && opt.fileType.compare("npsat_ascii") != 0 &&
        opt.fileType.compare("npsat_bin") != 0) {
        if (mpiRank == 0) {
            std::cerr << "Unsupported file_type: " << opt.fileType << std::endl;
        }
        localValid = 0;
    }

    int allValid = 0;
    MPI_Allreduce(&localValid, &allValid, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (!allValid) {
        MPI_Finalize();
        return 1;
    }

    const int totalFiles = opt.nproc * opt.niter;
    if (mpiRank == 0) {
        std::cout << "Starting " << totalFiles << " input files on " << mpiSize
                  << " MPI processes; progress interval " << opt.progressPercent
                  << "%" << std::endl;
    }

    const std::chrono::steady_clock::time_point beginTime = std::chrono::steady_clock::now();
    const int progressTag = 101;
    int completed = 0;
    int nextPercent = opt.progressPercent;
    int localFailures = 0;
    int filesInPart = 0;
    int partId = 0;
    bool outputOpen = false;
    OutputFiles outputs;

    for (int workId = mpiRank; workId < totalFiles; workId += mpiSize) {
        if (!outputOpen || (opt.outputFilesPerPart > 0 &&
                            filesInPart == opt.outputFilesPerPart)) {
            if (outputOpen && !closeOutputFiles(opt, mpiRank, outputs)) {
                ++localFailures;
            }
            if (outputOpen) {
                ++partId;
            }
            outputs = OutputFiles();
            outputOpen = openOutputFiles(opt, mpiRank, partId, outputs);
            filesInPart = 0;
        }

        bool success = outputOpen && processInputFile(workId, opt, outputs, mpiRank);
        if (!success) {
            ++localFailures;
        }
        ++filesInPart;

        const int status = success ? 1 : 0;
        if (mpiRank == 0) {
            ++completed;
            int messageWaiting = 0;
            MPI_Status messageStatus;
            do {
                MPI_Iprobe(MPI_ANY_SOURCE, progressTag, MPI_COMM_WORLD,
                           &messageWaiting, &messageStatus);
                if (messageWaiting) {
                    int workerStatus = 0;
                    MPI_Recv(&workerStatus, 1, MPI_INT, messageStatus.MPI_SOURCE,
                             progressTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    ++completed;
                }
            } while (messageWaiting);
            reportProgress(completed, totalFiles, opt.progressPercent, nextPercent);
        }
        else {
            MPI_Send(&status, 1, MPI_INT, 0, progressTag, MPI_COMM_WORLD);
        }
    }

    if (outputOpen && !closeOutputFiles(opt, mpiRank, outputs)) {
        ++localFailures;
    }

    if (mpiRank == 0) {
        while (completed < totalFiles) {
            int workerStatus = 0;
            MPI_Recv(&workerStatus, 1, MPI_INT, MPI_ANY_SOURCE, progressTag,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            ++completed;
            reportProgress(completed, totalFiles, opt.progressPercent, nextPercent);
        }
        if (100 % opt.progressPercent != 0) {
            std::cout << "Progress: 100% (" << completed << "/" << totalFiles
                      << " input files completed)" << std::endl;
        }
    }

    int totalFailures = 0;
    MPI_Reduce(&localFailures, &totalFailures, 1, MPI_INT, MPI_SUM, 0,
               MPI_COMM_WORLD);

    if (mpiRank == 0) {
        const std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
        std::cout << "Done in "
                  << std::chrono::duration_cast<std::chrono::microseconds>(endTime - beginTime).count()/1000000.0/60.0
                  << " minutes with " << totalFailures << " failure(s)" << std::endl;
    }

    MPI_Finalize();
    return totalFailures == 0 ? 0 : 1;
}
