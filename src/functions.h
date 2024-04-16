//    h files of the functions
#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <utility>         // For std::pair
#include <Eigen/Dense> 

std::pair<Eigen::VectorXd, Eigen::MatrixXd> lipkin_tot(int N, double chi);
Eigen::VectorXd get_entropy(int n, const Eigen::MatrixXd& coefficients);
void saveDataWithParams(const Eigen::VectorXd& Energy, const Eigen::VectorXd& Entropy, const std::string& filename, int N, double chi, int n);
void savecoefficients(const Eigen::MatrixXd& coefficient, const std::string& filename);
void get_system(int N, double chi, int n);

#endif
