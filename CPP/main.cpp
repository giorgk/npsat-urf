#include <chrono>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

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
        if (opt.er_to_run < 0) {
            NPSATurf(strmlnSeg, streamlineLength, velMult, opt, fp);
        }
        else {
            if (trajectory.end_reason == opt.er_to_run) {
                NPSATurf(strmlnSeg, streamlineLength, velMult, opt, fp);
            }
            else {
                fp.setVal(0.0);
            }
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

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cout << "Usage: NPSAT_URF <process_id> or NPSAT_URF -v" << std::endl;
        return 1;
    }

    std::string inputArg(argv[1]);
    if (inputArg.compare("-v") == 0) {
        std::cout << "version 1.3.0" << std::endl;
        return 0;
    }

    URFoptions opt;
    if (!readOptionFile(opt)) {
        return 1;
    }
    opt.ProcId = std::atoi(argv[1]);

    if (opt.fileType.compare("modpath") == 0) {
        std::cout << "file_type modpath is not available yet." << std::endl;
        return 1;
    }

    if (opt.fileType.compare("npsat_ascii") != 0 &&
        opt.fileType.compare("npsat_bin") != 0) {
        std::cout << "Unsupported file_type: " << opt.fileType << std::endl;
        return 1;
    }

    const std::string filename = opt.prefixInput + num2Padstr(opt.ProcId, opt.paddingZeros) + "." + opt.suffixInput;
    std::cout << "Reading: " << filename << std::endl;

    const std::string outfile = opt.prefixOutput + "_" + std::to_string(opt.ProcId) + ".dat";
    const std::string discardFilename = opt.prefixDiscard + "_" + std::to_string(opt.ProcId) + ".dat";
    const std::string simplifiedFilename = opt.prefixSimplified + "_" + std::to_string(opt.ProcId) + ".dat";
    std::cout << "Output file: " << outfile << std::endl;
    std::cout << "Discard file: " << discardFilename << std::endl;
    if (opt.simplifyStreamline) {
        std::cout << "Simplified streamline file: " << simplifiedFilename << std::endl;
    }

    std::ofstream ofile(outfile.c_str());
    if (!ofile.good()) {
        std::cout << "Can't open the output file " << outfile << std::endl;
        return 1;
    }
    writeOutputHeader(ofile, opt);
    std::cout << "Output file header prepared" << std::endl;

    std::ofstream discardFile(discardFilename.c_str());
    if (!discardFile.good()) {
        std::cout << "Can't open the discard file " << discardFilename << std::endl;
        return 1;
    }
    discardFile << "Eid, Sid, samples, has_termination, termination_pid, ER, reason" << std::endl;

    std::ofstream simplifiedFile;
    if (opt.simplifyStreamline) {
        simplifiedFile.open(simplifiedFilename.c_str());
        if (!simplifiedFile.good()) {
            std::cout << "Can't open the simplified streamline file " << simplifiedFilename << std::endl;
            return 1;
        }
        simplifiedFile << "Eid,Sid,x,y,z,v,a" << std::endl;
    }

    std::ifstream datafile;
    if (opt.fileType.compare("npsat_bin") == 0) {
        datafile.open(filename.c_str(), std::ios::binary);
    }
    else {
        datafile.open(filename.c_str());
    }

    if (!datafile.good()) {
        std::cout << "Can't open the file " << filename << std::endl;
        return 1;
    }

    std::chrono::steady_clock::time_point beginTimeALL = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point beginTime = std::chrono::steady_clock::now();
    int cntStrml = 0;
    StreamlineTrajectory trajectory;

    try {
        while (true) {
            bool found = false;
            if (opt.fileType.compare("npsat_bin") == 0) {
                found = readNextBinaryStreamline(datafile, trajectory, filename, discardFile);
            }
            else {
                found = readNextAsciiStreamline(datafile, trajectory, discardFile);
            }

            if (!found) {
                break;
            }

            writeSimplifiedStreamline(trajectory, opt, simplifiedFile, discardFile);
            processCompleteStreamline(trajectory, opt, ofile, discardFile, cntStrml, beginTime);
        }
    }
    catch (const std::exception& e) {
        std::cout << e.what() << std::endl;
        return 1;
    }

    const std::chrono::steady_clock::time_point endTimeALL = std::chrono::steady_clock::now();
    std::cout << "Done in "
              << std::chrono::duration_cast<std::chrono::microseconds>(endTimeALL - beginTimeALL).count()/1000000.0/60.0
              << std::endl;

    return 0;
}
