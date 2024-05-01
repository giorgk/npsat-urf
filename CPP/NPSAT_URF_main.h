//
// Created by giorgk on 4/30/2024.
//

#ifndef NPSAT_URF_NPSAT_URF_MAIN_H
#define NPSAT_URF_NPSAT_URF_MAIN_H

bool NPSATurf(std::vector<trajP>& S, double SLen, URFoptions& opt){

    // This is the discretized streamline
    std::vector<trajP>Sdx;
    int np = S.size();

    double aL = std::pow(opt.alpha*SLen, opt.beta);

    //Sdx.push_back(S[np-1])
    double Lel, vel, Del;
    double ddx, ddy, ddz, rho_term;
    for (int i = np-1; i > 0; i = i - 1 ){
        vel = S[i].v;
        ddx = S[i].x - S[i-1].x;
        ddy = S[i].y - S[i-1].y;
        ddz = S[i].z - S[i-1].z;
        Lel = std::sqrt(ddx*ddx + ddy*ddy + ddz*ddz);

        Del = aL * vel + opt.Dm;
        
        //rho_term = 1.0 +opt.rho_b*opt.K_d/n



        std::cout << i << std::endl;
    }


    return true;
}


#endif //NPSAT_URF_NPSAT_URF_MAIN_H
