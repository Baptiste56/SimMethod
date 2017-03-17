#include "basic.h"

namespace GBM {
    namespace MC {
        namespace Direct {
            std::vector<double> call_price(double &S_0, double &K, double &r, double &sigma, double &T, int &M) {
                std::clock_t t = clock();
                double sum_price = 0.;
                double sum_2_price = 0.;
                double D = exp(-r * T);
                for (int i = 0; i < M; ++i) {
                    double X_1 = (double) rand() / RAND_MAX;
                    double X_2 = (double) rand() / RAND_MAX;
                    double W = sqrt(-2 * log(X_1)) * sin(2 * M_PI * X_2);
                    double S_T = S_0 * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    sum_price += D * std::max(S_T - K, 0.);
                    sum_2_price += std::pow(D * std::max(S_T - K, 0.), 2);
                }
                double m = sum_price / (double) M;
                double v = ((1 / (double) M) * sum_2_price - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }

            std::vector<double>
            call_delta(double &S_0, double &K, double &r, double &sigma, double &T, int &M, double d_S) {
                std::clock_t t = clock();
                double sum_delta = 0.;
                double sum_2_delta = 0.;
                double D = exp(-r * T);
                for (int i = 0; i < M; ++i) {
                    double X_1 = (double) rand() / RAND_MAX;
                    double X_2 = (double) rand() / RAND_MAX;
                    double W = sqrt(-2 * log(X_1)) * sin(2 * M_PI * X_2);
                    double S_T_m = (S_0 - d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double S_T_h = (S_0 + d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double price_m = D * std::max(S_T_m - K, 0.);
                    double price_h = D * std::max(S_T_h - K, 0.);
                    sum_delta += (price_h - price_m) / (2 * d_S);
                    sum_2_delta += std::pow((price_h - price_m) / (2 * d_S), 2);
                }
                double m = (sum_delta) / (double) M;
                double v = ((1 / (double) M) * sum_2_delta - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }

            std::vector<double>
            call_gamma(double &S_0, double &K, double &r, double &sigma, double &T, int &M, double d_S) {
                std::clock_t t = clock();
                double sum_gamma = 0.;
                double sum_2_gamma = 0.;
                double D = exp(-r * T);
                for (int i = 0; i < M; ++i) {
                    double X_1 = (double) rand() / RAND_MAX;
                    double X_2 = (double) rand() / RAND_MAX;
                    double W = sqrt(-2 * log(X_1)) * sin(2 * M_PI * X_2);
                    double S_T_l = (S_0 - d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double S_T_m = S_0 * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double S_T_h = (S_0 + d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double price_l = D * std::max(S_T_l - K, 0.);
                    double price_m = D * std::max(S_T_m - K, 0.);
                    double price_h = D * std::max(S_T_h - K, 0.);
                    sum_gamma += (price_h - 2 * price_m + price_l) / (d_S * d_S);
                    sum_2_gamma += std::pow((price_h - 2 * price_m + price_l) / (d_S * d_S), 2);
                }
                double m = (sum_gamma) / (double) M;
                double v = ((1 / (double) M) * sum_2_gamma - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }

            std::vector<double>
            call_vega(double &S_0, double &K, double &r, double &sigma, double &T, int &M, double d_sigma) {
                std::clock_t t = clock();
                double sum_vega = 0.;
                double sum_2_vega = 0.;
                double D = exp(-r * T);
                for (int i = 0; i < M; ++i) {
                    double X_1 = (double) rand() / RAND_MAX;
                    double X_2 = (double) rand() / RAND_MAX;
                    double W = sqrt(-2 * log(X_1)) * sin(2 * M_PI * X_2);
                    double S_T_m = S_0 * exp(T * (r - 0.5 * (sigma - d_sigma) * (sigma - d_sigma)) + (sigma - d_sigma) * sqrt(T) * W);
                    double S_T_h = S_0 * exp(T * (r - 0.5 * (sigma + d_sigma) * (sigma + d_sigma)) + (sigma + d_sigma) * sqrt(T) * W);
                    double price_m = D * std::max(S_T_m - K, 0.);
                    double price_h = D * std::max(S_T_h - K, 0.);
                    sum_vega += (price_h - price_m) / (2*d_sigma);
                    sum_2_vega += std::pow((price_h - price_m) / d_sigma, 2);
                }
                double m = (sum_vega) / (double) M;
                double v = ((1 / (double) M) * sum_2_vega - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }
        }
        namespace Antithetic {
            std::vector<double> call_price(double &S_0, double &K, double &r, double &sigma, double &T, int &M) {
                std::clock_t t = clock();
                double sum_price = 0.;
                double sum_2_price = 0.;
                double D = exp(-r * T);
                for (int i = 0; i < M; ++i) {
                    double X_1 = (double) rand() / RAND_MAX;
                    double X_2 = (double) rand() / RAND_MAX;
                    double W = sqrt(-2 * log(X_1)) * sin(2 * M_PI * X_2);
                    double S_T_1 = S_0 * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double S_T_2 = S_0 * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * (-W));
                    sum_price += D * (std::max(S_T_1 - K, 0.) + std::max(S_T_2 - K, 0.)) / 2;
                    sum_2_price += std::pow(D * (std::max(S_T_1 - K, 0.) + std::max(S_T_2 - K, 0.)) / 2, 2);
                }
                double m = sum_price / (double) M;
                double v = ((1 / (double) M) * sum_2_price - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }

            std::vector<double>
            call_delta(double &S_0, double &K, double &r, double &sigma, double &T, int &M, double d_S) {
                std::clock_t t = clock();
                double sum_delta = 0.;
                double sum_2_delta = 0.;
                double D = exp(-r * T);
                for (int i = 0; i < M; ++i) {
                    double X_1 = (double) rand() / RAND_MAX;
                    double X_2 = (double) rand() / RAND_MAX;
                    double W = sqrt(-2 * log(X_1)) * sin(2 * M_PI * X_2);
                    double S_T_m_1 = (S_0 - d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double S_T_h_1 = (S_0 + d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double S_T_m_2 = (S_0 - d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * (-W));
                    double S_T_h_2 = (S_0 + d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * (-W));
                    double price_m_1 = D * std::max(S_T_m_1 - K, 0.);
                    double price_h_1 = D * std::max(S_T_h_1 - K, 0.);
                    double price_m_2 = D * std::max(S_T_m_2 - K, 0.);
                    double price_h_2 = D * std::max(S_T_h_2 - K, 0.);
                    sum_delta += ((price_h_1 - price_m_1) / (2 * d_S) + (price_h_2 - price_m_2) / (2 * d_S)) / 2;
                    sum_2_delta += std::pow(((price_h_1 - price_m_1) / (2 * d_S) + (price_h_2 - price_m_2) / (2 * d_S)) / 2, 2);
                }
                double m = (sum_delta) / (double) M;
                double v = ((1 / (double) M) * sum_2_delta - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }

            std::vector<double>
            call_gamma(double &S_0, double &K, double &r, double &sigma, double &T, int &M, double d_S) {
                std::clock_t t = clock();
                double sum_gamma = 0.;
                double sum_2_gamma = 0.;
                double D = exp(-r * T);
                for (int i = 0; i < M; ++i) {
                    double X_1 = (double) rand() / RAND_MAX;
                    double X_2 = (double) rand() / RAND_MAX;
                    double W = sqrt(-2 * log(X_1)) * sin(2 * M_PI * X_2);
                    double S_T_l_1 = (S_0 - d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double S_T_m_1 = S_0 * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double S_T_h_1 = (S_0 + d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double S_T_l_2 = (S_0 - d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * (-W));
                    double S_T_m_2 = S_0 * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * (-W));
                    double S_T_h_2 = (S_0 + d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * (-W));
                    double price_l_1 = D * std::max(S_T_l_1 - K, 0.);
                    double price_m_1 = D * std::max(S_T_m_1 - K, 0.);
                    double price_h_1 = D * std::max(S_T_h_1 - K, 0.);
                    double price_l_2 = D * std::max(S_T_l_2 - K, 0.);
                    double price_m_2 = D * std::max(S_T_m_2 - K, 0.);
                    double price_h_2 = D * std::max(S_T_h_2 - K, 0.);
                    sum_gamma += ((price_h_1 - 2 * price_m_1 + price_l_1) / (d_S * d_S) + (price_h_2 - 2 * price_m_2 + price_l_2) / (d_S * d_S)) / 2;
                    sum_2_gamma += std::pow(((price_h_1 - 2 * price_m_1 + price_l_1) / (d_S * d_S) + (price_h_2 - 2 * price_m_2 + price_l_2) / (d_S * d_S)) / 2, 2);
                }
                double m = (sum_gamma) / (double) M;
                double v = ((1 / (double) M) * sum_2_gamma - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }

            std::vector<double>
            call_vega(double &S_0, double &K, double &r, double &sigma, double &T, int &M, double d_sigma) {
                std::clock_t t = clock();
                double sum_vega = 0.;
                double sum_2_vega = 0.;
                double D = exp(-r * T);
                for (int i = 0; i < M; ++i) {
                    double X_1 = (double) rand() / RAND_MAX;
                    double X_2 = (double) rand() / RAND_MAX;
                    double W = sqrt(-2 * log(X_1)) * sin(2 * M_PI * X_2);
                    double S_T_m_1 = S_0 * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double S_T_h_1 = S_0 * exp(T * (r - 0.5 * (sigma + d_sigma) * (sigma + d_sigma)) + (sigma + d_sigma) * sqrt(T) * W);
                    double S_T_m_2 = S_0 * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * (-W));
                    double S_T_h_2 = S_0 * exp(T * (r - 0.5 * (sigma + d_sigma) * (sigma + d_sigma)) + (sigma + d_sigma) * sqrt(T) * (-W));
                    double price_m_1 = D * std::max(S_T_m_1 - K, 0.);
                    double price_h_1 = D * std::max(S_T_h_1 - K, 0.);
                    double price_m_2 = D * std::max(S_T_m_2 - K, 0.);
                    double price_h_2 = D * std::max(S_T_h_2 - K, 0.);
                    sum_vega += ((price_h_1 - price_m_1) / d_sigma + (price_h_2 - price_m_2) / d_sigma) / 2;
                    sum_2_vega += std::pow(((price_h_1 - price_m_1) / d_sigma + (price_h_2 - price_m_2) / d_sigma) / 2, 2);
                }
                double m = (sum_vega) / (double) M;
                double v = ((1 / (double) M) * sum_2_vega - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }
        }
    }
}