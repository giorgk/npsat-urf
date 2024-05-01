//
// Created by giorgk on 4/30/2024.
//

#ifndef NPSAT_URF_NPSAT_URF_MAIN_H
#define NPSAT_URF_NPSAT_URF_MAIN_H

//#include "helper_func.h"

bool NPSATurf(std::vector<segInfo>& S, double SLen, URFoptions& opt){

    // This is the discretized streamline
    int nel = S.size();

    double aL = opt.alpha*std::pow(SLen, opt.beta);

    //Sdx.push_back(S[np-1])
    double Lel, vel, Del;
    double ddx, ddy, ddz, rho_term;

    std::vector<eigenTriplet> DgloTri, MgloTri;
    eigenMat Dglo(nel+1,nel+1);
    eigenMat Mglo(nel+1,nel+1);

    //createTriplets(np, T);

    //Dglo.setFromTriplets(T.begin(), T.end());
    //Mglo.setFromTriplets(T.begin(), T.end());


    double DdiagPrev = 0;
    double MdiagPrev = 0;

    std::vector<segInfo>::reverse_iterator  rit;
    int ii = 0;// The index of the diagonal element
    for (rit = S.rbegin(); rit != S.rend(); ++rit){
        std::cout << ii << std::endl;
        vel = rit->v;

        Lel = rit->l;
        Del = aL * vel + opt.Dm;

        /*
              Dx   | 1  -1|    Vx  |-1  1|                         L |2  1|
           --------|      | + ---- |     | + lambda(1 + pho*Kd/n)*---|    |
              L    |-1   1|    2n  |-1  1|                         6 |1  2|
        */

        // Since we use Kd = 0 there is no need for that term
        rho_term = 1.0;// +opt.rho_b*opt.K_d/n
        double D11 =  Del/Lel  -  vel/2 + 2*opt.lambda*rho_term*Lel/6;
        double D12 = -Del/Lel  +  vel/2 + 1*opt.lambda*rho_term*Lel/6;
        double D21 = -Del/Lel  -  vel/2 + 1*opt.lambda*rho_term*Lel/6;
        double D22 =  Del/Lel  +  vel/2 + 2*opt.lambda*rho_term*Lel/6;

        /*
                            L  |2  1|
             (1+ pho*Kd/n) --- |    |
                            6  |1  2|
        */
        double A11 = 2*rho_term*Lel/6;
        double A12 = 1*rho_term*Lel/6;
        double A21 = 1*rho_term*Lel/6;
        double A22 = 2*rho_term*Lel/6;

        if (ii == 0){
            DgloTri.emplace_back(ii, ii, D11);
            DdiagPrev = D22;

            MgloTri.emplace_back(ii, ii, A11);
            MdiagPrev = A22;
        }
        else if (ii == nel-1){
            DgloTri.emplace_back(ii, ii, D11 + DdiagPrev);
            DgloTri.emplace_back(ii+1, ii+1, D22);

            MgloTri.emplace_back(ii, ii, A11 + MdiagPrev);
            MgloTri.emplace_back(ii+1, ii+1, A22);
        }
        else{
            DgloTri.emplace_back(ii, ii, D11 + DdiagPrev);
            DdiagPrev = D22;

            MgloTri.emplace_back(ii, ii, A11 + MdiagPrev);
            MdiagPrev = A22;
        }
        DgloTri.emplace_back(ii, ii+1, D12);
        DgloTri.emplace_back(ii+1, ii, D21);

        MgloTri.emplace_back(ii, ii+1, A12);
        MgloTri.emplace_back(ii+1, ii, A21);

        ii = ii + 1;
    }


    Dglo.setFromTriplets(DgloTri.begin(), DgloTri.end());
    Mglo.setFromTriplets(MgloTri.begin(), MgloTri.end());
    //std::cout << Dglo << std::endl;
    //std::cout << Mglo << std::endl;

    eigenMat Aglo, KK, GG, Bglo, BgloRed, RHS1;
    Aglo = Mglo+opt.wmega*opt.TimeStep*Dglo;
    std::cout << Aglo << std::endl;
    Bglo = (Mglo-(1-opt.wmega)*opt.TimeStep*Dglo);
    std::cout << Bglo << std::endl;
    BgloRed = Bglo.block(1,0,nel,nel+1);
    std::cout << BgloRed << std::endl;

    KK = Aglo.block(1,1,nel,nel);
    std::cout << KK << std::endl;
    GG = Aglo.block(1,0,nel,1);
    std::cout << GG << std::endl;

    Eigen::VectorXd Cinit(nel+1);
    Cinit.setZero();
    std::cout << Cinit << std::endl;
    Cinit(0) = 1.0;
    std::cout << Cinit << std::endl;
    //eigenMat
    //Cinit




    int cnt = 0;
    while(true){
        cnt = cnt + 1;
        Eigen::VectorXd RHS2 = BgloRed*Cinit;
        std::cout << RHS2 << std::endl;



        if (cnt > 500){
            break;
        }
    }

    return true;
}


#endif //NPSAT_URF_NPSAT_URF_MAIN_H
