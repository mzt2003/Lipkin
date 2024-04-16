//      define the functions
#include "functions.h"
#include "wignerSymbols.h"  //contains the CG coefficients
#include <iostream>
#include <vector>
#include <Eigen/Dense> 
#include <cmath>
#include <fstream>

//////////////////////
// get the energies and wave functions of Lipkin model with para N and chi.
//////////////////////
std::pair<Eigen::VectorXd, Eigen::MatrixXd> lipkin_tot(int N, double chi) {
    double J = N / 2.0;
    double epsilon = 1;

    int size_H = 2 * J + 1;
    std::vector<double> mm(size_H);
    for (int i = 0; i < size_H; ++i) {
        mm[i] = -J + i;
    }

    double V = chi * epsilon / (N - 1);
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(size_H, size_H);

    for (int i = 0; i < size_H; ++i) {
        H(i, i) = epsilon * mm[i];
    }

    for (int i = 0; i < size_H - 2; ++i) {
        double m = mm[i];
        double temp1 = -V / 2;
        double t2 = J * (J + 1) - m * (m + 1);
        double t3 = J * (J + 1) - (m + 1) * (m + 2);
        double tt = temp1 * std::sqrt(t2 * t3);
        H(i, i + 2) = tt;
        H(i + 2, i) = tt;
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(H);
    if (eigensolver.info() != Eigen::Success) {
        std::cerr << "Eigenvalue decomposition failed!" << std::endl;
        return std::make_pair(Eigen::VectorXd(), Eigen::MatrixXd());
    }

    Eigen::VectorXd E = eigensolver.eigenvalues();
    Eigen::MatrixXd coefficients = eigensolver.eigenvectors();

    return std::make_pair(E, coefficients);
}
/////////////////////////////////////
//entropy is determined when the parameter of the system n, N, chi is known
//get the entropy given n and the wave function under basis |J,M> 
/////////////////////////////////////
Eigen::VectorXd get_entropy(int n, const Eigen::MatrixXd& coefficients) {
    int N = coefficients.rows() - 1;
    double J = N / 2.0;
    double j1 = n / 2.0;
    double j2 = (N - n) / 2.0;
    
    std::vector<double> mm1(2 * j1 + 1), mm2(2 * j2 + 1);
    for (int i = 0; i < mm1.size(); ++i) mm1[i] = -j1 + i;
    for (int i = 0; i < mm2.size(); ++i) mm2[i] = -j2 + i;

    int size_H = 2 * (j1 + j2) + 1;
    std::vector<Eigen::MatrixXd> reduced_rho(size_H, Eigen::MatrixXd::Zero(mm1.size(), mm1.size()));
	Eigen::MatrixXd cg=Eigen::MatrixXd::Zero(mm1.size(), mm2.size());
	Eigen::MatrixXd c=Eigen::MatrixXd::Zero(mm1.size(), mm2.size());

	std::cout << "Start CG coefficients table calculation" << std::endl;
	for (int ii=0; ii<mm1.size(); ++ii){
		for (int jj=0; jj<mm2.size(); ++jj){
            cg(ii, jj) = WignerSymbols::clebschGordan(j1, j2, J, mm1[ii], mm2[jj], mm1[ii] + mm2[jj]);		
		}
	}
	std::cout << "Finish CG coefficients table calculation" << std::endl;

	std::cout << "Start rho_matrix calculation" << std::endl;
    for (int ind_e = 0; ind_e < size_H; ++ind_e) {
    	std::cout << "Process:"<< ind_e << '/'<< size_H << std::endl;
        for (int i = 0; i < mm1.size(); ++i) {
            for (int j = 0; j < mm2.size(); ++j) {
                c(i, j) = cg(i, j) * coefficients(i + j, ind_e);                
            }
        }
        reduced_rho[ind_e] = c * c.transpose();
    }
    std::cout << "Finish rho_matrix calculation" << std::endl;

	// Saving and printing reduced_rho matrices
    for (int ind_e = 0; ind_e < size_H; ++ind_e) {
        std::string filename = "../data/reduced_rho_" + std::to_string(ind_e) + ".txt";
        std::ofstream file(filename);
        if (file.is_open()) {
         file << reduced_rho[ind_e] << "\n";
            file.close();
        } else {
            std::cerr << "Unable to open file: " << filename << std::endl;
        }
        
        //  print to console
       // std::cout << "Reduced Rho Matrix " << ind_e << ":\n" << reduced_rho[ind_e] << "\n\n";
    }
    
	std::cout << "Start diagonalizing rho_matrix " << std::endl;    
    Eigen::VectorXd entropy = Eigen::VectorXd::Zero(size_H);
    for (int ind_e = 0; ind_e < size_H; ++ind_e) {
        Eigen::VectorXd eigvals = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(reduced_rho[ind_e]).eigenvalues();
        for (int i = 0; i < eigvals.size(); ++i) {
        	if (eigvals[i] > 0) {
            	entropy(ind_e) -= eigvals[i] * std::log(eigvals[i]);
            }else {
                std::cerr << "Warning: log of non-positive number, ind_e: " << ind_e
                          << ", i: " << i << ", eigval: " << eigvals[i] << std::endl;
            }	
        }
    }
	std::cout << "End diagonalizing rho_matrix " << std::endl;    
    return entropy;
}

/////////////////////
//save energy-entropy
/////////////////////
void saveDataWithParams(const Eigen::VectorXd& Energy, const Eigen::VectorXd& Entropy, const std::string& filename, int N, double chi, int n) {
    std::ofstream file(filename);
    if (file.is_open()) {
        // para
        file << N << " " << chi << " " << n << "\n";

        // results: energy,  entropy
        for (int i = 0; i < Energy.size(); ++i) {
            file << Energy(i) << " " << Entropy(i) << "\n";
        }
        file.close();
    } else {
        std::cerr << "Unable to open file";
    }
}

/////////////////////
//save coefficients
/////////////////////
void savecoefficients(const Eigen::MatrixXd& coefficient, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        file << coefficient << "\n";
        file.close();
    } else {
        std::cerr << "Unable to open file";
    }
}


//////////////////
//calculate a certain system
//////////////////
void get_system(int N, double chi, int n){
	auto [E, coefficients] = lipkin_tot(N, chi);
    auto entropy = get_entropy(n, coefficients);

    std::cout << "Energy:\n" << E << std::endl;
    //std::cout << "coefficients:\n" << coefficients << std::endl;
	std::cout << "Entropy:\n" << entropy << std::endl;
	
	saveDataWithParams(E, entropy, "../data/Energy_Entropy.txt", N, chi, n);
	savecoefficients(coefficients, "../data/coefficients.txt") ;
}





