#include "heston.h"

namespace Heston {
    namespace MC {
        namespace Direct {
            std::vector<double>
            call_price(double &S_0, double &K, double &r, double &kappa, double &theta, double &sigma, double &rho,
                       double &T, int &M) {
                if (2 * kappa * theta < sigma * sigma) {
                    std::cout << "Condition to ensure V_t > 0 not satisfied" << std::endl;
                    std::exit(0);
                }

                std::clock_t t = clock();
                double sum_price = 0.;
                double sum_2_price = 0.;
                double D = exp(-r * T);
                int N = 256;
                double dt = T / (double) N;
                for (int i = 0; i < M; ++i) {
                    double S = S_0;
                    double V = 0.04;
                    for (int j = 0; j < N; ++j) {
                        double X_1 = (double) rand() / RAND_MAX;
                        double X_2 = (double) rand() / RAND_MAX;
                        double Z_1 = sqrt(std::abs(2 * log(X_1))) * sin(2 * M_PI * X_2);
                        double Z_2 = sqrt(std::abs(2 * log(X_1))) * cos(2 * M_PI * X_2);
                        double dW_v = sqrt(dt) * Z_1;
                        double dW_s = sqrt(dt) * (rho * Z_1 + sqrt(1 - rho * rho) * Z_2);
                        double d_V = kappa * (theta - V) * dt + sigma * sqrt(std::abs(V)) * dW_v;
                        double d_S = r * S * dt + S * sqrt(std::abs(V)) * dW_s;
                        V += d_V;
                        S += d_S;
                    }
                    sum_price += D * std::max(S - K, 0.);
                    sum_2_price += std::pow(D * std::max(S - K, 0.), 2);
                }
                double m = sum_price / (double) M;
                double v = ((1 / (double) M) * sum_2_price - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }

            std::vector<double>
            call_delta(double &S_0, double &K, double &r, double &kappa, double &theta, double &sigma, double &rho,
                       double &T, int &M, double &delta_S) {
                if (2 * kappa * theta < sigma * sigma) {
                    std::cout << "Condition to ensure V_t > 0 not satisfied" << std::endl;
                    std::exit(0);
                }

                std::clock_t t = clock();
                double sum_delta = 0.;
                double sum_2_delta = 0.;
                double D = exp(-r * T);
                int N = 256;
                double dt = T / (double) N;
                for (int i = 0; i < M; ++i) {
                    double S = S_0;
                    double V = 0.04;
                    for (int j = 0; j < N; ++j) {
                        double X_1 = (double) rand() / RAND_MAX;
                        double X_2 = (double) rand() / RAND_MAX;
                        double Z_1 = sqrt(std::abs(2 * log(X_1))) * sin(2 * M_PI * X_2);
                        double Z_2 = sqrt(std::abs(2 * log(X_1))) * cos(2 * M_PI * X_2);
                        double dW_v = sqrt(dt) * Z_1;
                        double dW_s = sqrt(dt) * (rho * Z_1 + sqrt(1 - rho * rho) * Z_2);
                        double d_V = kappa * (theta - V) * dt + sigma * sqrt(std::abs(V)) * dW_v;
                        double d_S = r * S * dt + S * sqrt(std::abs(V)) * dW_s;
                        V += d_V;
                        S += d_S;
                    }
                    double price_h = D * std::max(S + delta_S - K, 0.);
                    double price_m = D * std::max(S - K, 0.);
                    sum_delta += (price_h - price_m) / delta_S;
                    sum_2_delta += std::pow((price_h - price_m) / delta_S, 2);
                }
                double m = sum_delta / (double) M;
                double v = ((1 / (double) M) * sum_2_delta - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }

            std::vector<double>
            call_gamma(double &S_0, double &K, double &r, double &kappa, double &theta, double &sigma, double &rho,
                       double &T, int &M, double &delta_S) {
                if (2 * kappa * theta < sigma * sigma) {
                    std::cout << "Condition to ensure V_t > 0 not satisfied" << std::endl;
                    std::exit(0);
                }

                std::clock_t t = clock();
                double sum_gamma = 0.;
                double sum_2_gamma = 0.;
                double D = exp(-r * T);
                int N = 256;
                double dt = T / (double) N;
                for (int i = 0; i < M; ++i) {
                    double S = S_0;
                    double V = 0.04;
                    for (int j = 0; j < N; ++j) {
                        double X_1 = (double) rand() / RAND_MAX;
                        double X_2 = (double) rand() / RAND_MAX;
                        double Z_1 = sqrt(std::abs(2 * log(X_1))) * sin(2 * M_PI * X_2);
                        double Z_2 = sqrt(std::abs(2 * log(X_1))) * cos(2 * M_PI * X_2);
                        double dW_v = sqrt(dt) * Z_1;
                        double dW_s = sqrt(dt) * (rho * Z_1 + sqrt(1 - rho * rho) * Z_2);
                        double d_V = kappa * (theta - V) * dt + sigma * sqrt(std::abs(V)) * dW_v;
                        double d_S = r * S * dt + S * sqrt(std::abs(V)) * dW_s;
                        V += d_V;
                        S += d_S;
                    }
                    double price_h = D * std::max(S + delta_S - K, 0.);
                    double price_m = D * std::max(S - K, 0.);
                    double price_l = D * std::max(S - delta_S - K, 0.);
                    sum_gamma += (price_h - 2 * price_m + price_l) / (delta_S * delta_S);
                    sum_2_gamma += std::pow((price_h - 2 * price_m + price_l) / (delta_S * delta_S), 2);
                }
                double m = sum_gamma / (double) M;
                double v = ((1 / (double) M) * sum_2_gamma - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }

            std::vector<double>
            call_vega(double &S_0, double &K, double &r, double &kappa, double &theta, double &sigma, double &rho,
                      double &T, int &M, double &delta_sigma) {
                if (2 * kappa * theta < sigma * sigma) {
                    std::cout << "Condition to ensure V_t > 0 not satisfied" << std::endl;
                    std::exit(0);
                }

                std::clock_t t = clock();
                double sum_gamma = 0.;
                double sum_2_gamma = 0.;
                double D = exp(-r * T);
                int N = 256;
                double dt = T / (double) N;
                for (int i = 0; i < M; ++i) {
                    double S_h = S_0;
                    double S_m = S_0;
                    double V_h = 0.04;
                    double V_m = 0.04;
                    for (int j = 0; j < N; ++j) {
                        double X_1 = (double) rand() / RAND_MAX;
                        double X_2 = (double) rand() / RAND_MAX;
                        double Z_1 = sqrt(std::abs(2 * log(X_1))) * sin(2 * M_PI * X_2);
                        double Z_2 = sqrt(std::abs(2 * log(X_1))) * cos(2 * M_PI * X_2);
                        double dW_v = sqrt(dt) * Z_1;
                        double dW_s = sqrt(dt) * (rho * Z_1 + sqrt(1 - rho * rho) * Z_2);
                        double d_V_m = kappa * (theta - V_m) * dt + sigma * sqrt(std::abs(V_m)) * dW_v;
                        double d_V_h = d_V_m + delta_sigma / N;
                        double d_S_h = r * S_h * dt + S_h * sqrt(std::abs(V_h)) * dW_s;
                        double d_S_m = r * S_m * dt + S_m * sqrt(std::abs(V_m)) * dW_s;
                        V_h += d_V_h;
                        V_m += d_V_m;
                        S_h += d_S_h;
                        S_m += d_S_m;
                    }
                    double price_h = D * std::max(S_h - K, 0.);
                    double price_m = D * std::max(S_m - K, 0.);
                    sum_gamma += (price_h - price_m) / delta_sigma;
                    sum_2_gamma += std::pow((price_h - price_m) / delta_sigma, 2);
                }
                double m = sum_gamma / (double) M;
                double v = ((1 / (double) M) * sum_2_gamma - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }
        }
        namespace Antithetic {
            std::vector<double>
            call_price(double &S_0, double &K, double &r, double &kappa, double &theta, double &sigma, double &rho,
                       double &T, int &M) {
                if (2 * kappa * theta < sigma * sigma) {
                    std::cout << "Condition to ensure V_t > 0 not satisfied" << std::endl;
                    std::exit(0);
                }

                std::clock_t t = clock();
                double sum_price = 0.;
                double sum_2_price = 0.;
                double D = exp(-r * T);
                int N = 256;
                double dt = T / (double) N;
                for (int i = 0; i < M; ++i) {
                    double S_1 = S_0;
                    double V_1 = 0.04;
                    double S_2 = S_0;
                    double V_2 = 0.04;
                    for (int j = 0; j < N; ++j) {
                        double X_1 = (double) rand() / RAND_MAX;
                        double X_2 = (double) rand() / RAND_MAX;
                        double Z_1 = sqrt(std::abs(2 * log(X_1))) * sin(2 * M_PI * X_2);
                        double Z_2 = sqrt(std::abs(2 * log(X_1))) * cos(2 * M_PI * X_2);
                        double dW_v = sqrt(dt) * Z_1;
                        double dW_s = sqrt(dt) * (rho * Z_1 + sqrt(1 - rho * rho) * Z_2);
                        double d_V_1 = kappa * (theta - V_1) * dt + sigma * sqrt(std::abs(V_1)) * dW_v;
                        double d_S_1 = r * S_1 * dt + S_1 * sqrt(std::abs(V_1)) * dW_s;
                        double d_V_2 = kappa * (theta - V_2) * dt + sigma * sqrt(std::abs(V_2)) * (-dW_v);
                        double d_S_2 = r * S_2 * dt + S_2 * sqrt(std::abs(V_2)) * (-dW_s);
                        V_1 += d_V_1;
                        S_1 += d_S_1;
                        V_2 += d_V_2;
                        S_2 += d_S_2;
                    }
                    sum_price += D * (std::max(S_1 - K, 0.) + std::max(S_2 - K, 0.)) / 2;
                    sum_2_price += std::pow(D * (std::max(S_1 - K, 0.) + std::max(S_2 - K, 0.)) / 2, 2);
                }
                double m = sum_price / (double) M;
                double v = ((1 / (double) M) * sum_2_price - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }

            std::vector<double>
            call_delta(double &S_0, double &K, double &r, double &kappa, double &theta, double &sigma, double &rho,
                       double &T, int &M, double &delta_S) {
                if (2 * kappa * theta < sigma * sigma) {
                    std::cout << "Condition to ensure V_t > 0 not satisfied" << std::endl;
                    std::exit(0);
                }

                std::clock_t t = clock();
                double sum_delta = 0.;
                double sum_2_delta = 0.;
                double D = exp(-r * T);
                int N = 256;
                double dt = T / (double) N;
                for (int i = 0; i < M; ++i) {
                    double S_1 = S_0;
                    double V_1 = 0.04;
                    double S_2 = S_0;
                    double V_2 = 0.04;
                    for (int j = 0; j < N; ++j) {
                        double X_1 = (double) rand() / RAND_MAX;
                        double X_2 = (double) rand() / RAND_MAX;
                        double Z_1 = sqrt(std::abs(2 * log(X_1))) * sin(2 * M_PI * X_2);
                        double Z_2 = sqrt(std::abs(2 * log(X_1))) * cos(2 * M_PI * X_2);
                        double dW_v = sqrt(dt) * Z_1;
                        double dW_s = sqrt(dt) * (rho * Z_1 + sqrt(1 - rho * rho) * Z_2);
                        double d_V_1 = kappa * (theta - V_1) * dt + sigma * sqrt(std::abs(V_1)) * dW_v;
                        double d_S_1 = r * S_1 * dt + S_1 * sqrt(std::abs(V_1)) * dW_s;
                        double d_V_2 = kappa * (theta - V_2) * dt + sigma * sqrt(std::abs(V_2)) * (-dW_v);
                        double d_S_2 = r * S_2 * dt + S_2 * sqrt(std::abs(V_2)) * (-dW_s);
                        V_1 += d_V_1;
                        S_1 += d_S_1;
                        V_2 += d_V_2;
                        S_2 += d_S_2;
                    }
                    double price_1_h = D * std::max(S_1 + delta_S - K, 0.);
                    double price_1_m = D * std::max(S_1 - K, 0.);
                    double price_2_h = D * std::max(S_2 + delta_S - K, 0.);
                    double price_2_m = D * std::max(S_2 - K, 0.);
                    sum_delta += ((price_1_h - price_1_m) / delta_S + (price_2_h - price_2_m) / delta_S) / 2;
                    sum_2_delta += std::pow(((price_1_h - price_1_m) / delta_S + (price_2_h - price_2_m) / delta_S) / 2, 2);
                }
                double m = sum_delta / (double) M;
                double v = ((1 / (double) M) * sum_2_delta - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }

            std::vector<double>
            call_gamma(double &S_0, double &K, double &r, double &kappa, double &theta, double &sigma, double &rho,
                       double &T, int &M, double &delta_S) {
                if (2 * kappa * theta < sigma * sigma) {
                    std::cout << "Condition to ensure V_t > 0 not satisfied" << std::endl;
                    std::exit(0);
                }

                std::clock_t t = clock();
                double sum_gamma = 0.;
                double sum_2_gamma = 0.;
                double D = exp(-r * T);
                int N = 256;
                double dt = T / (double) N;
                for (int i = 0; i < M; ++i) {
                    double S_1 = S_0;
                    double V_1 = 0.04;
                    double S_2 = S_0;
                    double V_2 = 0.04;
                    for (int j = 0; j < N; ++j) {
                        double X_1 = (double) rand() / RAND_MAX;
                        double X_2 = (double) rand() / RAND_MAX;
                        double Z_1 = sqrt(std::abs(2 * log(X_1))) * sin(2 * M_PI * X_2);
                        double Z_2 = sqrt(std::abs(2 * log(X_1))) * cos(2 * M_PI * X_2);
                        double dW_v = sqrt(dt) * Z_1;
                        double dW_s = sqrt(dt) * (rho * Z_1 + sqrt(1 - rho * rho) * Z_2);
                        double d_V_1 = kappa * (theta - V_1) * dt + sigma * sqrt(std::abs(V_1)) * dW_v;
                        double d_S_1 = r * S_1 * dt + S_1 * sqrt(std::abs(V_1)) * dW_s;
                        double d_V_2 = kappa * (theta - V_2) * dt + sigma * sqrt(std::abs(V_2)) * (-dW_v);
                        double d_S_2 = r * S_2 * dt + S_2 * sqrt(std::abs(V_2)) * (-dW_s);
                        V_1 += d_V_1;
                        S_1 += d_S_1;
                        V_2 += d_V_2;
                        S_2 += d_S_2;
                    }
                    double price_1_h = D * std::max(S_1 + delta_S - K, 0.);
                    double price_1_m = D * std::max(S_1 - K, 0.);
                    double price_1_l = D * std::max(S_1 - delta_S - K, 0.);
                    double price_2_h = D * std::max(S_2 + delta_S - K, 0.);
                    double price_2_m = D * std::max(S_2 - K, 0.);
                    double price_2_l = D * std::max(S_2 - delta_S - K, 0.);
                    sum_gamma += ((price_1_h - 2 * price_1_m + price_1_l) / (delta_S * delta_S) + (price_2_h - 2 * price_2_m + price_2_l) / (delta_S * delta_S)) / 2;
                    sum_2_gamma += std::pow(((price_1_h - 2 * price_1_m + price_1_l) / (delta_S * delta_S) + (price_2_h - 2 * price_2_m + price_2_l) / (delta_S * delta_S)) / 2, 2);
                }
                double m = sum_gamma / (double) M;
                double v = ((1 / (double) M) * sum_2_gamma - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }

            std::vector<double>
            call_vega(double &S_0, double &K, double &r, double &kappa, double &theta, double &sigma, double &rho,
                      double &T, int &M, double &delta_sigma) {
                if (2 * kappa * theta < sigma * sigma) {
                    std::cout << "Condition to ensure V_t > 0 not satisfied" << std::endl;
                    std::exit(0);
                }

                std::clock_t t = clock();
                double sum_gamma = 0.;
                double sum_2_gamma = 0.;
                double D = exp(-r * T);
                int N = 256;
                double dt = T / (double) N;
                for (int i = 0; i < M; ++i) {
                    double S_1_h = S_0;
                    double S_1_m = S_0;
                    double V_1_h = 0.04;
                    double V_1_m = 0.04;
                    double S_2_h = S_0;
                    double S_2_m = S_0;
                    double V_2_h = 0.04;
                    double V_2_m = 0.04;
                    for (int j = 0; j < N; ++j) {
                        double X_1 = (double) rand() / RAND_MAX;
                        double X_2 = (double) rand() / RAND_MAX;
                        double Z_1 = sqrt(std::abs(2 * log(X_1))) * sin(2 * M_PI * X_2);
                        double Z_2 = sqrt(std::abs(2 * log(X_1))) * cos(2 * M_PI * X_2);
                        double dW_v = sqrt(dt) * Z_1;
                        double dW_s = sqrt(dt) * (rho * Z_1 + sqrt(1 - rho * rho) * Z_2);
                        double d_V_1_m = kappa * (theta - V_1_m) * dt + sigma * sqrt(std::abs(V_1_m)) * dW_v;
                        double d_V_1_h = d_V_1_m + delta_sigma / N;
                        double d_V_2_m = kappa * (theta - V_2_m) * dt + sigma * sqrt(std::abs(V_2_m)) * (-dW_v);
                        double d_V_2_h = d_V_2_m + delta_sigma / N;
                        double d_S_1_h = r * S_1_h * dt + S_1_h * sqrt(std::abs(V_1_h)) * dW_s;
                        double d_S_1_m = r * S_1_m * dt + S_1_m * sqrt(std::abs(V_1_m)) * dW_s;
                        double d_S_2_h = r * S_2_h * dt + S_2_h * sqrt(std::abs(V_2_h)) * (-dW_s);
                        double d_S_2_m = r * S_2_m * dt + S_2_m * sqrt(std::abs(V_2_m)) * (-dW_s);
                        V_1_h += d_V_1_h;
                        V_1_m += d_V_1_m;
                        S_1_h += d_S_1_h;
                        S_1_m += d_S_1_m;
                        V_2_h += d_V_2_h;
                        V_2_m += d_V_2_m;
                        S_2_h += d_S_2_h;
                        S_2_m += d_S_2_m;
                    }
                    double price_1_h = D * std::max(S_1_h - K, 0.);
                    double price_1_m = D * std::max(S_1_m - K, 0.);
                    double price_2_h = D * std::max(S_2_h - K, 0.);
                    double price_2_m = D * std::max(S_2_m - K, 0.);
                    sum_gamma += ((price_1_h - price_1_m) / delta_sigma + (price_2_h - price_2_m) / delta_sigma) / 2;
                    sum_2_gamma += std::pow(((price_1_h - price_1_m) / delta_sigma + (price_2_h - price_2_m) / delta_sigma) / 2, 2);
                }
                double m = sum_gamma / (double) M;
                double v = ((1 / (double) M) * sum_2_gamma - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }
        }
    }
}