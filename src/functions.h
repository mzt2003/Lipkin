//    h files of the functions
#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <utility>         // For std::pair
#include <Eigen/Dense> 
extern bool debug_mode;

std::pair<Eigen::VectorXd, Eigen::MatrixXd> lipkin_tot(int N, double chi);
Eigen::VectorXd get_entropy(int n, const Eigen::MatrixXd& coefficients);
void saveDataWithParams(const std::vector<Eigen::VectorXd>& results, 
                        const std::string& filename, 
                        const std::vector<std::string>& params);
void savecoefficients(const Eigen::MatrixXd& coefficient, const std::string& filename);
void get_system(int N, double chi, int n, const std::string& filename);
void gs_S_vs_chi(int N, int n, double chi_1, double chi_2, int num_chi, const std::string& filename);
void gs_S_vs_n(int N, double chi, int n_1, int n_2, int interval, const std::string& filename);
void gs_S_vs_chi_n(int N, double chi_1, double chi_2, int num_chi, int n_1, int n_2, int interval, const std::string& filename);

#endif
