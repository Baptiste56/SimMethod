#include <iostream>
#include <iomanip>
#include "basic.h"
#include "analytic.h"
#include "heston.h"
#include "plot.h"

int main() {

    double S_0 = 100.;
    double K = 100. ;
    double r = 0.05;
    double sigma = 0.4;
    double T = 1;

    //compute the greeks
    double delta_S = 0.01;
    double delta_sigma = 0.01;

    //number of simulations
    int M [3] = {100, 10000, 100000};

    std::cout << "\n" << std::endl;
    std::cout << std::setw(100) << "****************************************************************" << std::endl;
    std::cout << std::setw(100) << "*************************  Basic Task  *************************" << std::endl;
    std::cout << std::setw(100) << "****************************************************************" << std::endl;
    std::cout << "-- Analytic --" << std::endl;
    double a_price = GBM::Analytic::call_price(S_0, K, r, sigma, T);
    double a_delta = GBM::Analytic::call_delta(S_0, K, r, sigma, T);
    double a_gamma = GBM::Analytic::call_gamma(S_0, K, r, sigma, T);
    double a_vega =  GBM::Analytic::call_vega(S_0, K, r, sigma, T);
    std::cout<< " price: " << a_price <<
              "\n delta: " << a_delta <<
              "\n gamma: " << a_gamma <<
              "\n vega:  " << a_vega  <<
             std::endl;

    ////////plot price against S_0/////////
    //////parameters///////
    int N = 100;
    double S_range [2] = {20., 200.};
    double V_range [2] = {0., 1.};
    double T_range [2] = {1., 10.};
    //////////////////////
    Plot *plt_1 = new Plot_GBM(Spot, S_range[0], S_range[1], N);
    plt_1->set_values(Analytic, None, S_0, K, r, sigma, T, M[1], delta_S);
    plt_1->show();
    delete plt_1;

    Plot *plt_2 = new Plot_GBM(Volatility, V_range[0], V_range[1], N);
    plt_2->set_values(Analytic, None, S_0, K, r, sigma, T, M[1], delta_S);
    plt_2->show();
    delete plt_2;

    Plot *plt_3 = new Plot_GBM(Maturity, T_range[0], T_range[1], N);
    plt_3->set_values(Analytic, None, S_0, K, r, sigma, T, M[1], delta_S);
    plt_3->show();
    delete plt_3;

    std::cout << "\n--- MC-Direct Method" << std::endl;
    std::cout << std::setw(40) << "Value" << std::setw(30) << "Relative error" << std::setw(30) << "Half-Width" << std::setw(30) << "Computation Time" << std:: endl;
    for(int i(0); i< 3; ++i)
    {
        std::vector<double> mc_price = GBM::MC::Direct::call_price(S_0, K, r, sigma, T, M[i]);
        std::vector<double> mc_delta = GBM::MC::Direct::call_delta(S_0, K, r, sigma, T, M[i], delta_S);
        std::vector<double> mc_gamma = GBM::MC::Direct::call_gamma(S_0, K, r, sigma, T, M[i], delta_S);
        std::vector<double> mc_vega  = GBM::MC::Direct::call_vega(S_0, K, r, sigma, T, M[i], delta_sigma);
        std::cout << "\n-- Simulations: " << M[i] <<
                  "\n   price: " << std::setw(30) << mc_price[0] << std::setw(30) << (mc_price[0] - a_price) / a_price << std::setw(30) << mc_price[1] << std::setw(30) << mc_price[2] <<
                  "\n   delta: " << std::setw(30) << mc_delta[0] << std::setw(30) << (mc_delta[0] - a_delta) / a_delta << std::setw(30) << mc_delta[1] << std::setw(30) << mc_delta[2] <<
                  "\n   gamma: " << std::setw(30) << mc_gamma[0] << std::setw(30) << (mc_gamma[0] - a_gamma) / a_gamma << std::setw(30) << mc_gamma[1] << std::setw(30) << mc_gamma[2] <<
                  "\n   vega:  " << std::setw(30) << mc_vega[0]  << std::setw(30) << (mc_vega[0] - a_vega) / a_vega    << std::setw(30) << mc_vega[1]  << std::setw(30) << mc_vega[2]  <<
                  std::endl;
    }


    std::cout << "\n--- MC-Antithetic Method" << std::endl;
    std::cout << std::setw(40) << "Value" << std::setw(30) << "Relative error" << std::setw(30) << "Half-Width" << std::setw(30) << "Computation Time" << std:: endl;
    for(int i(0); i< 3; ++i)
    {
        std::vector<double> mc_price = GBM::MC::Antithetic::call_price(S_0, K, r, sigma, T, M[i]);
        std::vector<double> mc_delta = GBM::MC::Antithetic::call_delta(S_0, K, r, sigma, T, M[i], delta_S);
        std::vector<double> mc_gamma = GBM::MC::Antithetic::call_gamma(S_0, K, r, sigma, T, M[i], delta_S);
        std::vector<double> mc_vega  = GBM::MC::Antithetic::call_vega(S_0, K, r, sigma, T, M[i], delta_sigma);
        std::cout << "\n-- Simulations: " << M[i] <<
                  "\n   price: " << std::setw(30) << mc_price[0] << std::setw(30) << (mc_price[0] - a_price) / a_price << std::setw(30) << mc_price[1] << std::setw(30) << mc_price[2] <<
                  "\n   delta: " << std::setw(30) << mc_delta[0] << std::setw(30) << (mc_delta[0] - a_delta) / a_delta << std::setw(30) << mc_delta[1] << std::setw(30) << mc_delta[2] <<
                  "\n   gamma: " << std::setw(30) << mc_gamma[0] << std::setw(30) << (mc_gamma[0] - a_gamma) / a_gamma << std::setw(30) << mc_gamma[1] << std::setw(30) << mc_gamma[2] <<
                  "\n   vega:  " << std::setw(30) << mc_vega[0]  << std::setw(30) << (mc_vega[0] - a_vega) / a_vega    << std::setw(30) << mc_vega[1]  << std::setw(30) << mc_vega[2]  <<
                  std::endl;
    }

    Plot *plt_4 = new Plot_GBM(Spot, S_range[0], S_range[1], N);
    plt_4->set_values(MC, Antithetic, S_0, K, r, sigma, T, M[1], delta_S);
    plt_4->show();
    delete plt_4;

    Plot *plt_5 = new Plot_GBM(Volatility, V_range[0], V_range[1], N);
    plt_5->set_values(MC, Antithetic, S_0, K, r, sigma, T, M[1], delta_S);
    plt_5->show();
    delete plt_5;

    Plot *plt_6 = new Plot_GBM(Maturity, T_range[0], T_range[1], N);
    plt_6->set_values(MC, Antithetic, S_0, K, r, sigma, T, M[1], delta_S);
    plt_6->show();
    delete plt_6;

    std::cout << "\n" << std::endl;
    std::cout << std::setw(100) << "****************************************************************" << std::endl;
    std::cout << std::setw(100) << "************************* Heston Model *************************" << std::endl;
    std::cout << std::setw(100) << "****************************************************************" << std::endl;

    //////// parameters ///////
    double NM_price = 8.894869;
    S_0 = 100.;
    K = 100. ;
    r = 0.025;
    sigma = 0.3;
    T = 1.;
    double kappa = 1.5;
    double theta = 0.04;
    double rho = -0.9;

    delta_S = 0.01;
    delta_sigma = 0.01;

    N = 50;
    //////////////////////////

    std::cout << "\n--- MC-Direct Method" << std::endl;
    std::cout << std::setw(40) << "Value" << std::setw(30) << "Relative error" << std::setw(30) << "Half-Width" << std::setw(30) << "Computation Time" << std:: endl;
    for(int i(0); i< 3; ++i)
    {
        std::vector<double> mc_price = Heston::MC::Direct::call_price(S_0, K, r, kappa, theta, sigma, rho, T, M[i]);
        std::vector<double> mc_delta = Heston::MC::Direct::call_delta(S_0, K, r, kappa, theta, sigma, rho, T, M[i], delta_S);
        std::vector<double> mc_gamma = Heston::MC::Direct::call_gamma(S_0, K, r, kappa, theta, sigma, rho, T, M[i], delta_S);
        std::vector<double> mc_vega  = Heston::MC::Direct::call_vega(S_0, K, r, kappa, theta, sigma, rho, T, M[i], delta_sigma);
        std::cout << "\n-- Simulations: " << M[i] <<
                  "\n   price: " << std::setw(30) << mc_price[0] << std::setw(30) << (mc_price[0] - NM_price) / NM_price << std::setw(30) << mc_price[1] << std::setw(30) << mc_price[2] <<
                  "\n   delta: " << std::setw(30) << mc_delta[0] << std::setw(30) << "N/A"                             << std::setw(30) << mc_delta[1] << std::setw(30) << mc_delta[2] <<
                  "\n   gamma: " << std::setw(30) << mc_gamma[0] << std::setw(30) << "N/A"                             << std::setw(30) << mc_gamma[1] << std::setw(30) << mc_gamma[2] <<
                  "\n   vega:  " << std::setw(30) << mc_vega[0]  << std::setw(30) << "N/A"                             << std::setw(30) << mc_vega[1]  << std::setw(30) << mc_vega[2]  <<
                  std::endl;
    }

    std::cout << "\n--- MC-Antithetic Method" << std::endl;
    std::cout << std::setw(40) << "Value" << std::setw(30) << "Relative error" << std::setw(30) << "Half-Width" << std::setw(30) << "Computation Time" << std:: endl;
    for(int i(0); i< 3; ++i)
    {
        std::vector<double> mc_price = Heston::MC::Antithetic::call_price(S_0, K, r, kappa, theta, sigma, rho, T, M[i]);
        std::vector<double> mc_delta = Heston::MC::Antithetic::call_delta(S_0, K, r, kappa, theta, sigma, rho, T, M[i], delta_S);
        std::vector<double> mc_gamma = Heston::MC::Antithetic::call_gamma(S_0, K, r, kappa, theta, sigma, rho, T, M[i], delta_S);
        std::vector<double> mc_vega  = Heston::MC::Antithetic::call_vega(S_0, K, r, kappa, theta, sigma, rho, T, M[i], delta_sigma);
        std::cout << "\n-- Simulations: " << M[i] <<
                  "\n   price: " << std::setw(30) << mc_price[0] << std::setw(30) << (mc_price[0] - NM_price) / NM_price << std::setw(30) << mc_price[1] << std::setw(30) << mc_price[2] <<
                  "\n   delta: " << std::setw(30) << mc_delta[0] << std::setw(30) << "N/A"                             << std::setw(30) << mc_delta[1] << std::setw(30) << mc_delta[2] <<
                  "\n   gamma: " << std::setw(30) << mc_gamma[0] << std::setw(30) << "N/A"                             << std::setw(30) << mc_gamma[1] << std::setw(30) << mc_gamma[2] <<
                  "\n   vega:  " << std::setw(30) << mc_vega[0]  << std::setw(30) << "N/A"                             << std::setw(30) << mc_vega[1]  << std::setw(30) << mc_vega[2]  <<
                  std::endl;
    }

    Plot *plt_7 = new Plot_Heston(Spot, S_range[0], S_range[1], N);
    plt_7->set_values(MC, Antithetic, S_0, K, r, kappa, theta, sigma, rho, T, M[1], delta_S, delta_sigma);
    plt_7->show();
    delete plt_7;

    Plot *plt_8 = new Plot_Heston(Maturity, T_range[0], T_range[1], N);
    plt_8->set_values(MC, Antithetic, S_0, K, r, kappa, theta, sigma, rho, T, M[1], delta_S, delta_sigma);
    plt_8->show();
    delete plt_8;

    return 0;
}
