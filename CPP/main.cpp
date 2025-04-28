#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <chrono>

#include "my_structures.h"
#include "NPSAT_URF_main.h"

int main(int argc, char *argv[]) {

    std::string input_arg(argv[1]);
    if (input_arg.compare("-v") == 0){
        std::cout << "version 1.2.1" << std::endl;
        return 0;
    }

    URFoptions opt;
    readOptionFile(opt);
    opt.ProcId = atoi(argv[1]);


    // Start with reading the file
    std::string filename = opt.prefixInput + num2Padstr(opt.ProcId, opt.paddingZeros) + "." + opt.suffixInput;
    std::cout << "Reading: " << filename << std::endl;
    std::ifstream datafile(filename.c_str());
    if (!datafile.good()) {
        std::cout << "Can't open the file " << filename << std::endl;
        return false;
    }
    else{
        std::string outfile;
        outfile = opt.prefixOutput + "_" + std::to_string(opt.ProcId) + ".dat";
        std::cout << "Output file: " << outfile << std::endl;
        std::ofstream ofile(outfile.c_str());
        ofile << "Eid, Sid, ER, p_cdsX, p_cdsY, p_cdsZ, v_cds, p_lndX, p_lndY, Len";
        for (int i = opt.por.startValue; i <= opt.por.endValue; ++opt.por.interval){
            ofile << ", Age" << i << ", mean" << i << ", std" << i << ", err" << i;
            if (opt.calcDecay){
                ofile << ", meanDc" << i << ", stdDc" << i << ", ScaleDc" << i << ", errDc" << i;
            }
            if (opt.calcDiff){
                ofile << ", meanDf" << i << ", stdDf" << i << ", ScaleDf" << i << ", errDf" << i;
            }
        }
        ofile << std::endl;

        std::chrono::steady_clock::time_point beginTimeALL = std::chrono::steady_clock::now();
        std::chrono::steady_clock::time_point beginTime = std::chrono::steady_clock::now();
        std::chrono::steady_clock::time_point endTime1, endTime2;

        std::string line;
        int tmpInt, Eid, Sid, iER;
        int cntStrml = 0;
        std::string er;
        trajP p;
        //std::vector<trajP> strmlnSeg;
        std::vector<segInfo> strmlnSeg;
        segInfo si;
        double StreamlineLength = 0;
        double ax, ay, az, bx, by, bz, v, segLen;
        double p_lndX, p_lndY, p_cdsX, p_cdsY, p_cdsZ;
        bool CompleteStreamlineFound = false;
        double leftOverLen = 0.0;
        bool bCalcLen = false;
        while (getline(datafile, line)){
            if (line.size() > 1){
                std::istringstream inp(line);
                inp >> tmpInt;
                if (tmpInt == -9){
                    inp >> Eid;
                    inp >> Sid;
                    inp >> er;
                    iER = parseExitReason(er);
                    strmlnSeg[strmlnSeg.size() - 1].l = strmlnSeg[strmlnSeg.size() - 1].l + leftOverLen;
                    StreamlineLength = StreamlineLength + leftOverLen;
                    CompleteStreamlineFound = true;
                    bCalcLen = false;
                    leftOverLen = 0;
                    p_lndX = ax;
                    p_lndY = ay;
                }
                else{
                    if (!opt.bIsGather){
                        inp >> Eid;
                    }
                    else{
                        Eid = tmpInt;
                    }
                    inp >> Sid;
                    inp >> bx;
                    inp >> by;
                    inp >> bz;
                    if (!opt.bIsGather){
                        double vx, vy, vz;
                        inp >> vx;
                        inp >> vy;
                        inp >> vz;
                        v = std::sqrt(vx*vx + vy*vy + vz*vz);
                    }
                    else{
                        inp >> v;
                    }

                    if (bCalcLen){
                        segLen = std::sqrt((bx - ax)*(bx - ax) + (by - ay)*(by - ay) + (bz - az)*(bz - az)) + leftOverLen;
                        if (segLen < opt.minElemSize){
                            leftOverLen = leftOverLen + segLen;
                        }
                        else if (segLen > opt.maxElemSize){
                            double nSplits = std::ceil(segLen/opt.maxElemSize);
                            int inSplits = static_cast<int>(nSplits);
                            double lenSplit = segLen/nSplits;
                            for (int i = 0; i < inSplits; ++i){
                                strmlnSeg.emplace_back(v, lenSplit);
                                StreamlineLength = StreamlineLength + lenSplit;
                            }
                            leftOverLen = 0.0;
                        }
                        else{
                            strmlnSeg.emplace_back(v, segLen);
                            StreamlineLength = StreamlineLength + segLen;
                            leftOverLen = 0.0;
                        }
                    }
                    else{
                        bCalcLen = true;
                        p_cdsX = bx;
                        p_cdsY = by;
                        p_cdsZ = bz;
                    }
                    ax = bx;
                    ay = by;
                    az = bz;
                }
            }
            endTime1 = std::chrono::steady_clock::now();

            if (CompleteStreamlineFound){
                std::cout << Eid << " " << Sid << " " << ++cntStrml << " [" << std::chrono::duration_cast<std::chrono::microseconds>(endTime1 - beginTime).count()/1000000.0;
                //std::vector<FittedParam> AllPorFP;
                FittedParam fp;

                ofile << Eid << ", " << Sid << ", " << iER << ", "
                      << std::setprecision(2) << std::fixed
                      << p_cdsX << ", " << p_cdsY << ", " << p_cdsZ << ", "
                      << std::setprecision(5) << strmlnSeg[0].v << ", " << std::setprecision(2)
                      << p_lndX << ", " << p_lndY << ", " << StreamlineLength << ", " ;
                for (int i = opt.por.startValue; i <= opt.por.endValue; ++opt.por.interval){
                    double velMult = static_cast<double>(i)/10.0;
                    fp.reset();
                    if (iER == 1){
                        bool tf = NPSATurf(strmlnSeg, StreamlineLength, velMult, opt, fp);
                    }
                    else{
                        fp.setVal(0.0);
                    }

                    ofile << std::setprecision(2) << std::fixed << fp.Age << std::setprecision(6) << std::scientific
                          << ", " << fp.urf.m << ", " << fp.urf.s << ", " << fp.urf.err;
                    if (opt.calcDecay){
                        ofile << ", " << fp.Decay.m << ", " << fp.Decay.s << ", " << fp.Decay.sc << ", " << fp.Decay.err;
                    }
                    if (opt.calcDiff){
                        ofile << ", " << fp.Diff.m << ", " << fp.Diff.s << ", " << fp.Diff.sc << ", " << fp.Diff.err;
                    }
                    if (i != opt.por.endValue){
                        ofile << ", ";
                    }
                    //AllPorFP.push_back(fp);
                    //velMult = velMult + 1.0;
                }
                ofile << std::endl;

                endTime2 = std::chrono::steady_clock::now();

                std::cout << ", " << std::chrono::duration_cast<std::chrono::microseconds>(endTime2 - beginTime).count()/1000000.0 << "]" << std::endl;

                beginTime = std::chrono::steady_clock::now();

                CompleteStreamlineFound = false;
                StreamlineLength = 0.0;
                strmlnSeg.clear();
            }
        }
        std::chrono::steady_clock::time_point endTimeALL = std::chrono::steady_clock::now();
        std::cout << "Done in " << std::chrono::duration_cast<std::chrono::microseconds>(endTimeALL - beginTimeALL).count()/1000000.0/60.0 << std::endl;
        ofile.close();
    }


    return 0;
}
