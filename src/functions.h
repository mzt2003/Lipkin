//    h files of the functions
#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <utility>         // For std::pair
#include <Eigen/Dense> 

using namespace Eigen;
using namespace std;

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
void get_sys2(int n1, int n2, double ep1, double ep2, double V1, double V2, double V12);

class QuantumSystem {
private:
    int n1, n2, s1, s2;
    double ep1, ep2, V1, V2, V12, j1, j2;
    vector<double> m1, m2;
    MatrixXd H;
    MatrixXd H1, H2;
    VectorXd eigenvalues;
    MatrixXd eigenvectors;
    VectorXd entropy;
    VectorXd H1_expected, H2_expected;
    std::vector<Eigen::MatrixXd> reduced_rho_1;
    std::vector<Eigen::MatrixXd> reduced_rho_2;
    double gs_E, gs_S;
    bool hamiltonianComputed = false;
    bool eigenSystemSolved = false;
    bool entropyComputed = false;
    bool H1Computed = false;
    bool H2Computed = false;
    bool H1_expected_Computed = false;
    bool H2_expected_Computed = false;
    bool reduced_rho_1_Computed = false;
    bool reduced_rho_2_Computed = false;

public:
    // Constructor
    QuantumSystem(int n1, int n2, double ep1, double ep2, double V1, double V2, double V12);
    // Methods
    void computeHamiltonian();
    void solveEigenSystem();
    void computeEntropy();
    void computeH1();
	void computeH2();
	void computeH1_expected();
	void computeH2_expected();
    void outputResults(const string& filename);
    
    const MatrixXd& getHamiltonian() ;
    const VectorXd& getEigenvalues();
    const MatrixXd& getEigenvectors(); 
    const VectorXd& getEntropy() ;
    double getGroundStateEnergy();
    double getGroundStateEntropy();
    const VectorXd& getH1_expected();
    const VectorXd& getH2_expected();
    
};

#endif
