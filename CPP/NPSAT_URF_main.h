//
// Created by giorgk on 4/30/2024.
//

#ifndef NPSAT_URF_NPSAT_URF_MAIN_H
#define NPSAT_URF_NPSAT_URF_MAIN_H

#include "helper_func.h"

bool NPSATurf(std::vector<segInfo>& S, double SLen, double velMult, URFoptions& opt, FittedParam& fp){

    double lambda = halfTime(12.32);
    // This is the discretized streamline
    int nel = S.size();

    double aL = opt.alpha*std::pow(SLen, opt.beta);

    //Sdx.push_back(S[np-1])
    double Lel, vel, Del;
    double ddx, ddy, ddz, rho_term;

    std::vector<eigenTriplet> DgloTri, MgloTri, DgloTriDecay;
    eigenMat Dglo(nel+1,nel+1);
    eigenMat Mglo(nel+1,nel+1);
    eigenMat DgloDecay(nel+1,nel+1);


    double DdiagPrev = 0;
    double MdiagPrev = 0;
    double DdiagPrevDecay = 0;
    double MdiagPrevDecay = 0;

    std::vector<segInfo>::reverse_iterator  rit;
    int ii = 0;// The index of the diagonal element
    double age = 0;
    double len = 0;
    for (rit = S.rbegin(); rit != S.rend(); ++rit){
        //std::cout << ii << std::endl;
        vel = rit->v / velMult;
        Lel = rit->l;
        age = age + Lel / vel;
        len = len + Lel;

        Del = aL * vel + opt.Dm;

        /*
              Dx   | 1  -1|    Vx  |-1  1|                         L |2  1|
           --------|      | + ---- |     | + lambda(1 + pho*Kd/n)*---|    |
              L    |-1   1|    2n  |-1  1|                         6 |1  2|
        */

        // Since we use Kd = 0 there is no need for that term
        //rho_term = 1.0;// +opt.rho_b*opt.K_d/n
        // Terms without decay
        double D11 =  Del/Lel  -  vel/2;
        double D12 = -Del/Lel  +  vel/2;
        double D21 = -Del/Lel  -  vel/2;
        double D22 =  Del/Lel  +  vel/2;

        // Terms with decay
        double D11dec =  Del/Lel  -  vel/2 + 2*lambda*Lel/6;
        double D12dec = -Del/Lel  +  vel/2 + 1*lambda*Lel/6;
        double D21dec = -Del/Lel  -  vel/2 + 1*lambda*Lel/6;
        double D22dec =  Del/Lel  +  vel/2 + 2*lambda*Lel/6;

        /*
                            L  |2  1|
             (1+ pho*Kd/n) --- |    |
                            6  |1  2|
        */
        double A11 = 2*Lel/6;
        double A12 = 1*Lel/6;
        double A21 = 1*Lel/6;
        double A22 = 2*Lel/6;

        if (ii == 0){
            DgloTri.emplace_back(ii, ii, D11);
            DdiagPrev = D22;

            DgloTriDecay.emplace_back(ii, ii, D11dec);
            DdiagPrevDecay = D22dec;

            MgloTri.emplace_back(ii, ii, A11);
            MdiagPrev = A22;
        }
        else if (ii == nel-1){
            DgloTri.emplace_back(ii, ii, D11 + DdiagPrev);
            DgloTri.emplace_back(ii+1, ii+1, D22);

            DgloTriDecay.emplace_back(ii, ii, D11dec + DdiagPrevDecay);
            DgloTriDecay.emplace_back(ii+1, ii+1, D22dec);

            MgloTri.emplace_back(ii, ii, A11 + MdiagPrev);
            MgloTri.emplace_back(ii+1, ii+1, A22);
        }
        else{
            DgloTri.emplace_back(ii, ii, D11 + DdiagPrev);
            DdiagPrev = D22;

            DgloTriDecay.emplace_back(ii, ii, D11 + DdiagPrevDecay);
            DdiagPrevDecay = D22dec;

            MgloTri.emplace_back(ii, ii, A11 + MdiagPrev);
            MdiagPrev = A22;
        }
        DgloTri.emplace_back(ii, ii+1, D12);
        DgloTri.emplace_back(ii+1, ii, D21);

        DgloTriDecay.emplace_back(ii, ii+1, D12dec);
        DgloTriDecay.emplace_back(ii+1, ii, D21dec);

        MgloTri.emplace_back(ii, ii+1, A12);
        MgloTri.emplace_back(ii+1, ii, A21);

        ii = ii + 1;
    }
    fp.Age = age/365.0;
    fp.Len = len;
    if (fp.Age < opt.skipAge){
        for (int i = 0; i < opt.skipAge; ++i){
            double d =static_cast<double>(i);
            if (fp.Age <= d){
                fp.setVal(-d);
                return true;
            }
        }
    }



    Dglo.setFromTriplets(DgloTri.begin(), DgloTri.end());
    DgloDecay.setFromTriplets(DgloTriDecay.begin(), DgloTriDecay.end());
    Mglo.setFromTriplets(MgloTri.begin(), MgloTri.end());
    //std::cout << Dglo << std::endl;
    //std::cout << Mglo << std::endl;

    eigenMat Aglo, KK, GG, Bglo, BgloRed, RHS1;
    eigenMat AgloDecay, BgloDecay, BgloRedDecay, KKDecay, GGDecay;

    Aglo = Mglo+opt.wmega*opt.TimeStep*Dglo;
    AgloDecay = Mglo+opt.wmega*opt.TimeStep*DgloDecay;

    //std::cout << Aglo << std::endl;
    Bglo = (Mglo-(1-opt.wmega)*opt.TimeStep*Dglo);
    BgloDecay = (Mglo-(1-opt.wmega)*opt.TimeStep*DgloDecay);

    //std::cout << Bglo << std::endl;
    BgloRed = Bglo.block(1,0,nel,nel+1);
    BgloRedDecay = BgloDecay.block(1,0,nel,nel+1);
    //std::cout << BgloRed << std::endl;

    KK = Aglo.block(1,1,nel,nel);
    KKDecay = AgloDecay.block(1,1,nel,nel);
    //std::cout << KK << std::endl;

    GG = Aglo.block(1,0,nel,1);
    GGDecay = AgloDecay.block(1,0,nel,1);
    //std::cout << GG << std::endl;

    Eigen::VectorXd Cprev(nel+1);
    Eigen::VectorXd CprevDecay(nel+1);
    Cprev.setZero();
    CprevDecay.setZero();
    //std::cout << Cprev << std::endl;
    Cprev(0) = 1.0;
    CprevDecay(0) = 1.0;
    //std::cout << Cprev << std::endl;
    //eigenMat
    //Cinit

    //Eigen::BiCGSTAB<eigenMat > solver;
    Eigen::SparseLU<eigenMat> solver;
    Eigen::SparseLU<eigenMat> solverDecay;
    solver.compute(KK);
    if (solver.info() != Eigen::Success){
        fp.setVal(-77);
        return false;
    }
    solverDecay.compute(KKDecay);
    if (solverDecay.info() != Eigen::Success){
        fp.setVal(-77);
        return false;
    }

    Eigen::VectorXd C(nel);
    Eigen::VectorXd CDecay(nel);
    std::vector<double> unitBTC;
    std::vector<double> unitBTCDecay;

    Eigen::VectorXd RHS2(nel);
    Eigen::VectorXd RHS2Decay(nel);
    Eigen::VectorXd RHS(nel);
    Eigen::VectorXd RHSDecay(nel);

    int cnt = 0;
    bool bExceedMaxTime = false;
    while(true){

        RHS2 = BgloRed*Cprev;
        RHS2Decay = BgloRedDecay*CprevDecay;
        //std::cout << RHS2 << std::endl;
        RHS = RHS2 - GG;
        RHSDecay = RHS2Decay - GGDecay;

        C = solver.solve(RHS);
        if (solver.info() != Eigen::Success){
            fp.setVal(-99);
            return false;
        }
        CDecay = solverDecay.solve(RHSDecay);
        if (solverDecay.info() != Eigen::Success){
            fp.setVal(-99);
            return false;
        }
        //std::cout << "C:" << std::endl;
        //std::cout << C << std::endl;
        Cprev.segment(1,nel) = C;
        CprevDecay.segment(1,nel) = CDecay;
        //std::cout << "Cprev:" << std::endl;
        //std::cout << Cprev << std::endl;
        //std::cout << Cprev.size() << std::endl;
        //std::cout << Cprev(nel) << std::endl;
        //std::cout << URF.size() << std::endl;

        unitBTC.push_back(Cprev(nel));
        unitBTCDecay.push_back(CprevDecay(nel));

        //std::cout << Cprev(nel) << " " << CprevDecay(nel) << std::endl;
        if (Cprev(nel) > opt.URFtol){
            break;
        }
        cnt = cnt + 1;
        if (cnt > opt.maxTotalTime){
            bExceedMaxTime = true;
            break;
        }
    }
    if (bExceedMaxTime){
        fp.setVal(-99);
        return true;
    }

    // Fitting URFs
    std::vector<double> URFDiff;
    data_samples DS_urf, DS_urfDecay, DS_urfDiff;
    input_vector input;
    double x = 1;
    double urf_val, urf_dec_val, urf_diff_val;
    double scaleDecay = 1/unitBTCDecay.back();
    double scaleDiff = 1/(unitBTC.back() - unitBTCDecay.back());
    double maxYurf = 0.0;
    double maxYurfDec = 0.0;
    double maxYurfDiff = 0.0;
    double maxYPosurf = 0.0;
    double maxYPosurfDec = 0.0;
    double maxYPosurfDiff = 0.0;
    for (unsigned int i = 0; i < unitBTC.size(); ++i){
        input(0) = x;
        if (i == 0){
            urf_val = unitBTC[0];
            urf_dec_val = unitBTCDecay[0];
            urf_diff_val = urf_val - urf_dec_val;
        }
        else{
            urf_val = unitBTC[i] - unitBTC[i-1];
            urf_dec_val = unitBTCDecay[i] - unitBTCDecay[i-1];
            urf_diff_val = urf_val - urf_dec_val;
        }
        if (urf_val > maxYurf){
            maxYurf = urf_val;
            maxYPosurf = x;
        }
        if (urf_dec_val > maxYurfDec){
            maxYurfDec = urf_dec_val;
            maxYPosurfDec = x;
        }
        if (urf_diff_val > maxYurfDiff){
            maxYurfDiff = urf_diff_val;
            maxYPosurfDiff = x;
        }

        //std::cout << urf_val << " " << urf_dec_val << " " << urf_dec_val*scaleDecay << " "
        //          << urf_diff_val << " " << urf_diff_val*scaleDiff << std::endl;

        DS_urf.emplace_back(input, urf_val);
        DS_urfDecay.emplace_back(input, urf_dec_val*scaleDecay);
        DS_urfDiff.emplace_back(input, urf_diff_val*scaleDiff);

        x = x + 1.0;
    }

    parameter_vector Xurf, XurfDec, XurfDiff;
    bool tf;
    tf = fitLgnrm(DS_urf,maxYPosurf, Xurf);
    if (tf){
        fp.urf.err = fitError(DS_urf, Xurf);
        fp.urf.m = Xurf(0,0);
        fp.urf.s = Xurf(1,0);
        fp.urf.sc = 1.0;
    }
    else{
        fp.urf.m = -88.0;
        fp.urf.s = -88.0;
        fp.urf.sc = 1.0;
        fp.urf.err = -88.0;
    }

    tf = fitLgnrm(DS_urfDecay,maxYPosurfDec, XurfDec);
    if (tf){
        fp.Decay.err = fitError(DS_urfDecay, XurfDec);
        fp.Decay.m = XurfDec(0,0);
        fp.Decay.s = XurfDec(1,0);
        fp.Decay.sc = scaleDecay;
    }
    else{
        fp.Decay.m = -88.0;
        fp.Decay.s = -88.0;
        fp.Decay.sc = scaleDecay;
        fp.Decay.err = -88.0;
    }

    tf = fitLgnrm(DS_urfDiff,maxYPosurfDiff, XurfDiff);
    if (tf){
        fp.Diff.err = fitError(DS_urfDiff, XurfDiff);
        fp.Diff.m = XurfDiff(0,0);
        fp.Diff.s = XurfDiff(1,0);
        fp.Diff.sc = scaleDiff;
    }
    else{
        fp.Diff.m = -88.0;
        fp.Diff.s = -88.0;
        fp.Diff.sc = scaleDiff;
        fp.Diff.err = -88.0;
    }

    return true;
}


#endif //NPSAT_URF_NPSAT_URF_MAIN_H
