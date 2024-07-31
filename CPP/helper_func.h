//
// Created by giorg on 5/1/2024.
//

#ifndef NPSAT_URF_HELPER_FUNC_H
#define NPSAT_URF_HELPER_FUNC_H

#include "my_structures.h"


bool readOptionFile(URFoptions& opt){
    const std::string optfile = "npsat_urf.opt";
    std::ifstream datafile(optfile.c_str());
    if (!datafile.good()) {
        std::cout << "Can't open the file " << optfile << std::endl;
        return false;
    }
    else{
        std::string line, propname;
        int count = 0;
        while (getline(datafile, line)){
            std::istringstream inp(line);
            inp >> propname;
            if (propname.compare("prefix") == 0){
                inp >> opt.prefixInput;
                continue;
            }
            if (propname.compare("suffix") == 0){
                inp >> opt.suffixInput;
                continue;
            }
            if (propname.compare("output_prefix") == 0){
                inp >> opt.prefixOutput;
                continue;
            }
            if (propname.compare("alpha") == 0){
                inp >> opt.alpha;
                continue;
            }
            if (propname.compare("beta") == 0){
                inp >> opt.beta;
                continue;
            }
            if (propname.compare("Dm") == 0){
                inp >> opt.Dm;
                continue;
            }
            if (propname.compare("minElemSize") == 0){
                inp >> opt.minElemSize;
                continue;
            }
            if (propname.compare("maxElemSize") == 0){
                inp >> opt.maxElemSize;
                continue;
            }
            if (propname.compare("TimeStep") == 0){
                inp >> opt.TimeStep;
                continue;
            }
            if (propname.compare("maxTotalTime") == 0){
                inp >> opt.maxTotalTime;
                continue;
            }
            if (propname.compare("URFtol") == 0){
                inp >> opt.URFtol;
                continue;
            }
            if (propname.compare("skipAge") == 0){
                inp >> opt.skipAge;
                continue;
            }
            if (propname.compare("paddingZeros") == 0){
                inp >> opt.paddingZeros;
                continue;
            }
            if (propname.compare("startPor") == 0){
                inp >> opt.startPor;
                continue;
            }
            if (propname.compare("endPor") == 0){
                inp >> opt.endPor;
                continue;
            }
            count++;
            if (count > 20){
                break;
            }
        }
        return true;
    }

}

double halfTime(double val){
    return std::log(2.0)/(val*365.0);
}


// We will use this function to generate data.  It represents a function of 2 variables
// and 3 parameters.   The least squares procedure will be used to infer the values of
// the 3 parameters based on a set of input/output pairs.
double lgnrmlFunction(const input_vector& input, const parameter_vector& params){
    const double a = params(0);
    const double b = params(1);

    const double x = input(0);

    /*
               /               2 \
               |   (a - log(x))  |
    sqrt(2) exp| - ------------- |
               |           2     |
               \        2 b      /
    ------------------------------
            b x sqrt(pi) 2
     */

    const double logx = std::log(x);
    const double nominator = dlib::sqrt_2 * exp( -((a - logx)*(a - logx))/(2.0*b*b) );
    const double denominator = b * x * sqrtpi * 2.0;
    double out = nominator/denominator;
    return out;
}

// This function is the "residual" for a least squares problem.   It takes an input/output
// pair and compares it to the output of our model and returns the amount of error.  The idea
// is to find the set of parameters which makes the residual small on all the data pairs.
double residual(const std::pair<input_vector, double>& data, const parameter_vector& params){
    return lgnrmlFunction(data.first, params) - data.second;
}

// This function is the derivative of the residual() function with respect to the parameters.
parameter_vector residual_derivative(const std::pair<input_vector, double>& data, const parameter_vector& params){
    /*
    wrt a (mean)
                 /               2 \
                 |   (a - log(x))  |
      sqrt(2) exp| - ------------- | (2 a - 2 log(x))
                 |           2     |
                 \        2 b      /
    - -----------------------------------------------
                       3
                      b  x sqrt(pi) 4

    wrt b (std)
               /               2 \                            /               2 \
               |   (a - log(x))  |             2              |   (a - log(x))  |
    sqrt(2) exp| - ------------- | (a - log(x))    sqrt(2) exp| - ------------- |
               |           2     |                            |           2     |
               \        2 b      /                            \        2 b      /
    -------------------------------------------- - ------------------------------
                    4                                       2
                   b  x sqrt(pi) 2                         b  x sqrt(pi) 2
        */
    parameter_vector der;
    const double a = params(0);
    const double b = params(1);

    const double x = data.first(0);
    const double logx = std::log(x);
    const double a_logxsq = (a - logx)*(a - logx);
    const double sqb = b*b;
    const double exp_fraction = dlib::sqrt_2 * std::exp( -a_logxsq/(2*sqb));

    const double numeratorA = exp_fraction * (2 * a - 2 * logx);
    const double denominatorA = b*sqb * x * sqrtpi * 4;

    const double numeratorB = exp_fraction * a_logxsq;
    const double denominatorB = sqb*sqb * x * sqrtpi * 2;

    der(0) = -numeratorA/denominatorA;
    der(1) = numeratorB/denominatorB - exp_fraction/(sqb * x * sqrtpi * 2);

    return der;
}

bool fitLgnrm(data_samples& DS, double& maxYpos, parameter_vector& x){
    parameter_vector params;
    params(0) = log(maxYpos);
    params(1) = 0.3;
    try{
        x = params;
        dlib::solve_least_squares_lm(dlib::objective_delta_stop_strategy(1e-7),
                               residual,
                               residual_derivative,
                               DS,
                               x);
        //std::cout << dlib::trans(x) << "| " << length(x - params) << std::endl;

        return true;
    }
    catch (std::exception& e)
    {
        std::cout << e.what() << std::endl;
        return false;
    }
}

double fitError(data_samples& DS, parameter_vector& params){
    double urfx;
    double sumSqErr = 0.0;
    double cumUrf = 0.0;
    double cntVal = 0.0;
    double err;
    for (unsigned int i = 0; i < DS.size(); ++i){
        urfx = lgnrmlFunction(DS[i].first, params);
        cumUrf = cumUrf + urfx;
        cntVal = cntVal +1.0;
        err = DS[i].second - urfx;
        sumSqErr = sumSqErr + err*err;
        if (cumUrf > 0.95){
            break;
        }
    }
    return std::sqrt(sumSqErr/cntVal);
}

std::string num2Padstr(int i, int n) {
    std::stringstream ss;
    ss << std::setw(n) << std::setfill('0') << i;
    return ss.str();
}


/*
void createTriplets(int n, std::vector<eigenTriplet>& T){
    eigenMat A(n,n);

    eigenTriplet tt;
    tt.value()
    for (int i = 0; i < n; ++i){
        if (i == 0){
            T.push_back(eigenTriplet(i,i,0.0));
            T.push_back(eigenTriplet(i,i+1, 0.0));
        }
        else if (i == n-1){
            T.push_back(eigenTriplet(i,i-1, 0.0));
            T.push_back(eigenTriplet(i,i, 0.0));
        }
        else{
            T.push_back(eigenTriplet(i,i-1, 0.0));
            T.push_back(eigenTriplet(i,i, 0.0));
            T.push_back(eigenTriplet(i,i+1, 0.0));
        }
    }
}
*/

#endif //NPSAT_URF_HELPER_FUNC_H
