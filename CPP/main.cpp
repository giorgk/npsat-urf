#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

#include "my_structures.h"
#include "NPSAT_URF_main.h"

int main(int argc, char *argv[]) {
    std::cout << argc << std::endl;
    if (argc < 3){
        std::cout << "Two inputs are needed" << std::endl;
        return 0;
    }
    std::cout << argv[0] << std::endl;
    std::cout << argv[1] << std::endl;
    std::cout << argv[2] << std::endl;
    std::cout << "Hello, World!" << std::endl;

    URFoptions opt;


    // Start with reading the file
    std::string filename = "teststrmlinfit.traj";
    std::ifstream datafile(filename.c_str());
    if (!datafile.good()) {
        std::cout << "Can't open the file " << filename << std::endl;
        return false;
    }
    else{
        std::string outfile(argv[2]);
        outfile = outfile + "_" + std::to_string(atoi(argv[1])) + ".dat";
        std::ofstream ofile(outfile.c_str());
        ofile << "Eid, Sid, Len, ER";
        for (int i = 1; i < 6;++i){
            ofile << ", Age" << i << ", mean" << i << ", std" << i
                  << ", meanDc" << i << ", stdDc" << i << ", ScaleDc" << i
                  << ", meanDf" << i << ", stdDf" << i << ", ScaleDf" << i;
        }
        ofile << std::endl;


        std::string line;
        int tmpInt, Eid, Sid, iER;
        int cntStrml = 0;
        std::string er;
        trajP p;
        //std::vector<trajP> S;
        std::vector<segInfo> S;
        segInfo si;
        double StreamlineLength = 0;
        double ax, ay, az, bx, by, bz, v, segLen;
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
                    S[S.size()-1].l = S[S.size()-1].l + leftOverLen;
                    StreamlineLength = StreamlineLength + leftOverLen;
                    CompleteStreamlineFound = true;
                    bCalcLen = false;
                    leftOverLen = 0;
                }
                else{
                    Eid = tmpInt;
                    inp >> Sid;
                    inp >> bx;
                    inp >> by;
                    inp >> bz;
                    inp >> v;
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
                                S.emplace_back(v, lenSplit);
                                StreamlineLength = StreamlineLength + lenSplit;
                            }
                            leftOverLen = 0;
                        }
                        else{
                            S.emplace_back(v, segLen);
                            StreamlineLength = StreamlineLength + segLen;
                            leftOverLen = 0;
                        }
                    }
                    else{
                        bCalcLen = true;
                    }
                    ax = bx;
                    ay = by;
                    az = bz;
                }
            }

            if (CompleteStreamlineFound){
                std::cout << Eid << " " << Sid << " " << ++cntStrml << std::endl;
                //std::vector<FittedParam> AllPorFP;
                FittedParam fp;
                double velMult = 1.0;
                ofile << Eid << ", " << Sid << ", " << StreamlineLength << ", " << iER << ", ";
                for (int i = 0; i < 5; ++i){
                    bool tf = NPSATurf(S, StreamlineLength, velMult,  opt, fp);
                    ofile << fp.Age << ", " << fp.urf.m << ", " << fp.urf.s << ", ";
                    ofile << fp.Decay.m << ", " << fp.Decay.s << ", " << fp.Decay.sc << ", ";
                    ofile << fp.Diff.m << ", " << fp.Diff.s << ", " << fp.Diff.sc << ", ";
                    //AllPorFP.push_back(fp);
                    velMult = velMult + 1.0;
                }
                ofile << std::endl;

                CompleteStreamlineFound = false;
            }
        }
        ofile.close();
    }


    return 0;
}
