// main.cpp
#include <fstream>
#include <vector>
#include <iostream> 
#include "functions.h"
#include <Eigen/Dense>
#include <cmath>

using namespace Eigen;
using namespace std;
int main() {
    get_system(100, 1.0, 20); //int N, double chi, int n
	///////
	
	////////////////////////



    int n1 =20;
    int n2 = 80;
    double j1 = n1 / 2.0;
    double j2 = n2 / 2.0;
    int s1 = n1 + 1;
    int s2 = n2 + 1;
    double ep1 = 1, ep2 = 1;
    double V1 = 1.0 / 99.0, V2 = 1.0 / 99.0, V12 = 1.0 / 99.0;

    vector<double> m1(s1), m2(s2);
    for (int i = 0; i < s1; i++) {
        m1[i] = -j1 + i;
    }
    for (int i = 0; i < s2; i++) {
        m2[i] = -j2 + i;
    }

    MatrixXd H = MatrixXd::Zero(s1 * s2, s1 * s2);

    // Populating the Hamiltonian matrix
    for (int ii = 0; ii < s1; ++ii) {
        for (int jj = 0; jj < s2; ++jj) {
            for (int kk = 0; kk < s1; ++kk) {
                for (int ll = 0; ll < s2; ++ll) {
                    if (ii == kk && jj == ll) {
                        H(ii * s2 + jj, kk * s2 + ll) = ep1 * m1[ii] + ep2 * m2[jj];
                    }
                    if (ii + 2 == kk && jj == ll) {
                        double t1 = sqrt(j1 * (j1 + 1) - m1[ii] * (m1[ii] + 1));
                        double t2 = sqrt(j1 * (j1 + 1) - (m1[ii] + 1) * (m1[ii] + 2));
                        H(ii * s2 + jj, kk * s2 + ll) = -V1 * t1 * t2 / 2;
                        H(kk * s2 + ll, ii * s2 + jj) = -V1 * t1 * t2 / 2;
                    }
                    if (ii == kk && jj + 2 == ll) {
                        double t1 = sqrt(j2 * (j2 + 1) - m2[jj] * (m2[jj] + 1));
                        double t2 = sqrt(j2 * (j2 + 1) - (m2[jj] + 1) * (m2[jj] + 2));
                        H(ii * s2 + jj, kk * s2 + ll) += -V2 * t1 * t2 / 2;
                        H(kk * s2 + ll, ii * s2 + jj) += -V2 * t1 * t2 / 2;
                    }
                    if (ii + 1 == kk && jj + 1 == ll) {
                        double t1 = sqrt(j1 * (j1 + 1) - m1[ii] * (m1[ii] + 1));
                        double t2 = sqrt(j2 * (j2 + 1) - m2[jj] * (m2[jj] + 1));
                        H(ii * s2 + jj, kk * s2 + ll) = -V12 * t1 * t2;
                        H(kk * s2 + ll, ii * s2 + jj) = -V12 * t1 * t2;
                    }
                }
            }
        }
    }
    
    
    // write H to file
    ofstream outFile("../data/H_matrix.txt");
    if (!outFile.is_open()) {
        cerr << "Failed to open file for writing." << endl;
        return -1;
    }
    for (int i = 0; i < H.rows(); ++i) {
        for (int j = 0; j < H.cols(); ++j) {
            outFile << H(i, j);
            if (j != H.cols() - 1) outFile << " ";  
        }
        outFile << "\n";  
    }
    outFile.close();
    
    // Compute eigenvalues and eigenvectors
    SelfAdjointEigenSolver<MatrixXd> eigensolver(H);
    if (eigensolver.info() != Success) abort();
    VectorXd eigenvalues = eigensolver.eigenvalues();
    MatrixXd eigenvectors = eigensolver.eigenvectors();

    // Compute entropies
    vector<MatrixXd> reduced_rho_1(s1 * s2, MatrixXd::Zero(s1,s1));
    VectorXd entropy1 = VectorXd::Zero(s1 * s2);
    for (int ii = 0; ii < s1 * s2; ++ii) {
        Map<MatrixXd> t1(eigenvectors.col(ii).data(), s2, s1);
        MatrixXd rho = t1.transpose() * t1;
        reduced_rho_1[ii] = rho;       
         std::string filename = "../data/reduced_rho111_" + std::to_string(ii) + ".txt";
        std::ofstream file1(filename);
        if (file1.is_open()) {
         file1 << reduced_rho_1[ii] << "\n";
            file1.close();
        } else {
            std::cerr << "Unable to open file: " << filename << std::endl;
        }        
     }   
    std::cout << "Start diagonalizing rho_matrix " << std::endl;    
    for (int ind_e = 0; ind_e < s1*s2; ++ind_e) {
        VectorXd eigvals = SelfAdjointEigenSolver<MatrixXd>(reduced_rho_1[ind_e]).eigenvalues();
       // std::cout << "eigenvals "<< eigvals << std::endl; 
        for (int i = 0; i < eigvals.size(); ++i) {
        	if (eigvals[i] > 0) {
            	entropy1(ind_e) -= eigvals[i] * std::log(eigvals[i]);
            }else {
                std::cerr << "Warning:n, ind_e: " << ind_e
                          << ", i: " << i << ", eigval: " << eigvals[i] << std::endl;
            }	
        }
    }
	 std::cout << "eigenvals "<< eigenvalues << std::endl;
    // Save data for plotting
    ofstream file("../data/Energy_Entropy22.txt");
    for (int i = 0; i < eigenvalues.size(); ++i) {
        file << eigenvalues[i] << " " << entropy1[i] << endl;
    }
    file.close();

	
	///////
	
    return 0;
}




