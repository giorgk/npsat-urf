#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

#include "my_structures.h"
#include "NPSAT_URF_main.h"

int main(int argc, char *argv[]) {
    std::cout << argc << std::endl;
    std::cout << argv[0] << std::endl;
    std::cout << argv[1] << std::endl;
    std::cout << argv[2] << std::endl;
    std::cout << "Hello, World!" << std::endl;

    URFoptions opt;


    // Start with reading the file
    std::string filename = "d:\\UCDAVIS\\npsat-urf\\test_data\\teststrmlinfit.traj";
    std::ifstream datafile(filename.c_str());
    if (!datafile.good()) {
        std::cout << "Can't open the file " << filename << std::endl;
        return false;
    }
    else{
        std::string line;
        int tmpInt, Eid, Sid, iER;
        int cntStrml = 0;
        std::string er;
        trajP p;
        std::vector<trajP> S;
        double StreamlineLength = 0;
        bool CompleteStreamlineFound = false;
        while (getline(datafile, line)){
            if (line.size() > 1){
                std::istringstream inp(line.c_str());
                inp >> tmpInt;
                if (tmpInt == -9){
                    inp >> Eid;
                    inp >> Sid;
                    inp >> er;
                    iER = parseExitReason(er);
                    CompleteStreamlineFound = true;
                }
                else{
                    Eid = tmpInt;
                    inp >> Sid;
                    inp >> p.x;
                    inp >> p.y;
                    inp >> p.z;
                    inp >> p.v;
                    if (S.size() > 1){
                        double ddx = p.x - S[S.size()-1].x;
                        double ddy = p.y - S[S.size()-1].y;
                        double ddz = p.z - S[S.size()-1].z;
                        StreamlineLength += std::sqrt(ddx*ddx + ddy*ddy + ddz*ddz);
                    }
                    S.push_back(p);

                }
            }

            if (CompleteStreamlineFound){
                std::cout << Eid << " " << Sid << " " << ++cntStrml << std::endl;
                bool tf = NPSATurf(S, StreamlineLength, opt);

            }
        }
    }


    return 0;
}
