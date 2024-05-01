//
// Created by giorgk on 4/30/2024.
//

#ifndef NPSAT_URF_MY_STRUCTURES_H
#define NPSAT_URF_MY_STRUCTURES_H

#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<double> eigenMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> eigenTriplet;

struct trajP{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double v = 0.0;
    double len = 0.0;
};

struct segInfo{
    segInfo(){
        v = 0;
        l = 0;
    }
    segInfo(double v_in, double l_in){
        v = v_in;
        l = l_in;
    }
    double v = 0;
    double l = 0;
};

struct URFoptions{
    double alpha = 0.32;
    double beta = 0.83;
    double Dm =1.1578e-4;
    double lambda = 0;
    double K_d = 0;
    double rho_b = 1.0;

    double minElemSize = 0.01;
    double maxElemSize = 300;
    double TimeStep = 365.0; // in years
    double maxTotalTime = 1000; // in years
    double wmega = 0.5;
};


int parseExitReason(std::string t){
    if (t.compare("EXIT_TOP") == 0){
        return 1;
    }
    else if (t.compare("EXIT_SIDE") == 0){
        return 2;
    }
    else if (t.compare("EXIT_BOTTOM") == 0){
        return 3;
    }
    else if (t.compare("MAX_INNER_ITER") == 0){
        return 4;
    }
    else if (t.compare("STUCK") == 0){
        return 5;
    }
    else if (t.compare("MAX_AGE") == 0){
        return 6;
    }
    else if (t.compare("ATTRACT") == 0){
        return 7;
    }
    else if (t.compare("INIT_OUT") == 0){
        return 8;
    }
    else if (t.compare("FIRST_POINT_GHOST") == 0){
        return 9;
    }
    else if (t.compare("FAR_AWAY") == 0){
        return 10;
    }
    else if (t.compare("NOT_IN_SUBDOMAIN") == 0){
        return 11;
    }
    else if (t.compare("CHANGE_PROCESSOR") == 0){
        return 12;
    }
    else if (t.compare("NO_EXIT") == 0){
        return 13;
    }
    else if (t.compare("NO_REASON") == 0){
        return 14;
    }
    else
        return 999;
}

#endif //NPSAT_URF_MY_STRUCTURES_H
