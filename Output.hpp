//
//  Output.hpp
//  SimMethod
//
//  Created by Portenard Baptiste on 17/03/2017.
//  Copyright © 2017 Portenard Baptiste. All rights reserved.
//

#ifndef Output_hpp
#define Output_hpp

#include <stdio.h>
#include "basic.h"
#include "heston.h"
#include "analytic.h"
#include <vector>

enum SolvingType{
    Analytic,
    MC //Monte Carlo
};
enum Method{
    Direct,
    Antithetic,
    ControlVariate,
    None
};
enum X_axis{
    Spot,
    Volatility,
    Maturity
};

class Output{
public:
    Output(X_axis X_type, double X_min, double X_max, int N): m_X_type(X_type), m_X_min(X_min), m_X_max(X_max), m_N(N) {}
    virtual ~Output(){}
    virtual void set_values(SolvingType RES, Method MTD, double &S, double &K, double &r, double &sigma, double &T, int &M, double d_S){};
    virtual void set_values(SolvingType RES, Method MTD, double &S, double &K, double &r, double &kappa, double &theta, double &sigma, double &rho, double &T, int &M, double &delta_S, double &delta_sigma){};
    virtual void show() = 0;
protected:
    std::vector<std::pair<double, double> > m_pts_price;
    std::vector<std::pair<double, double> > m_pts_delta;
    std::vector<std::pair<double, double> > m_pts_gamma;
    std::vector<std::pair<double, double> > m_pts_vega;
    X_axis                                  m_X_type;
    double                                  m_X_min;
    double                                  m_X_max;
    int                                     m_N;
};

class Output_GBM : public Output {
public:
    Output_GBM(X_axis X_type, double X_min, double X_max, int N): Output(X_type, X_min, X_max, N) {}
    ~Output_GBM(){}
    void set_values(SolvingType RES, Method MTD, double &S, double &K, double &r, double &sigma, double &T, int &M, double d_S);
    void show();
};

class Output_Heston : public Output {
public:
    Output_Heston(X_axis X_type, double X_min, double X_max, int N): Output(X_type, X_min, X_max, N) {}
    ~Output_Heston(){}
    void set_values(SolvingType RES, Method MTD, double &S, double &K, double &r, double &kappa, double &theta, double &sigma, double &rho, double &T, int &M, double &delta_S, double &delta_sigma);
    void show();
};

#endif /* Output_hpp */
