//      define the functions
#include "functions.h"
#include "wignerSymbols.h"  //contains the CG coefficients
#include <iostream>
#include <vector>
#include <Eigen/Dense> 
#include <cmath>
#include <fstream>

//////////////////////////////////////////////////////////////
// get the energies and wave functions of Lipkin model with para N and chi.
//////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////
//entropy is determined when the parameter of the system n, N, chi is known
//get the entropy given n and the wave function under basis |J,M> 
///////////////////////////////////////////////////////////////
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
    	//std::cout << "Process:"<< ind_e << '/'<< size_H << std::endl;
        for (int i = 0; i < mm1.size(); ++i) {
            for (int j = 0; j < mm2.size(); ++j) {
                c(i, j) = cg(i, j) * coefficients(i + j, ind_e);                
            }
        }
        reduced_rho[ind_e] = c * c.transpose();
    }
    std::cout << "Finish rho_matrix calculation" << std::endl;

	// Saving and printing reduced_rho matrices
	if (debug_mode){
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
    }
    
	std::cout << "Start diagonalizing rho_matrix " << std::endl;    
    Eigen::VectorXd entropy = Eigen::VectorXd::Zero(size_H);
    for (int ind_e = 0; ind_e < size_H; ++ind_e) {
        Eigen::VectorXd eigvals = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(reduced_rho[ind_e]).eigenvalues();
        for (int i = 0; i < eigvals.size(); ++i) {
        	if (eigvals[i] > 0) {
            	entropy(ind_e) -= eigvals[i] * std::log(eigvals[i]);
            }else {
                //std::cerr << "Warning: log of non-positive number, ind_e: " << ind_e
                //          << ", i: " << i << ", eigval: " << eigvals[i] << std::endl;
            }	
        }
    }
	std::cout << "End diagonalizing rho_matrix " << std::endl;    
    return entropy;
}

/////////////////////////////////////////////////
//save para and data
/////////////////////////////////////////////////
void saveDataWithParams(const std::vector<Eigen::VectorXd>& results, 
                        const std::string& filename, 
                        const std::vector<std::string>& params) {
 std::ofstream file(filename);
    if (file.is_open()) {
        for (const auto& param : params) {
            file << param << " ";
        }
        file << "\n";
        
        for (int i = 0; i < results.front().size(); ++i) {
            for (const auto& vec : results) {
                if (vec.size() > i) {  
                    file << vec(i) << " ";
                }
            }
            file << "\n";
        }
        file.close();
    } else {
        std::cerr << "Unable to open file\n";
    }
}


//////////////////////////////////////////////////
//save coefficients
//////////////////////////////////////////////////
void savecoefficients(const Eigen::MatrixXd& coefficient, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        file << coefficient << "\n";
        file.close();
    } else {
        std::cerr << "Unable to open file";
    }
}


/////////////////////////////////////////////////
//calculate a certain system
/////////////////////////////////////////////////
void get_system(int N, double chi, int n, const std::string& filename){
	auto [E, coefficients] = lipkin_tot(N, chi);
    auto entropy = get_entropy(n, coefficients);
    std::vector<Eigen::VectorXd> results = {E, entropy};
	std::vector<std::string> parameters = {std::to_string(N),std::to_string(chi),std::to_string(n)};
	saveDataWithParams(results, filename  , parameters);
    if(debug_mode){
    	std::cout << "Energy:\n" << E << std::endl;
    	//std::cout << "coefficients:\n" << coefficients << std::endl;
		std::cout << "Entropy:\n" << entropy << std::endl;
		savecoefficients(coefficients, "../data/coefficients.txt") ;
	}

}
/////////////////////////////////////////////////
//save gs_S_vs_chi
/////////////////////////////////////////////////
 void gs_S_vs_chi(int N, int n, double chi_1, double chi_2, int num_chi, const std::string& filename) {
    Eigen::VectorXd chis = Eigen::VectorXd::Zero(num_chi);
    Eigen::VectorXd gsentropy = Eigen::VectorXd::Zero(num_chi);

    for (int i = 0; i < num_chi; ++i) {
        double chi = chi_1+ (chi_2-chi_1) * i / (num_chi - 1); 
        chis(i) = chi;
        auto [E, coefficients] = lipkin_tot(N, chi);
        auto entropy = get_entropy(n, coefficients);
        gsentropy(i) = entropy(0); 
    }
    
    std::vector<Eigen::VectorXd> results = {chis, gsentropy};
	std::vector<std::string> parameters = {std::to_string(N),std::to_string(n)};
	saveDataWithParams(results,filename , parameters);

	}

/////////////////////////////////////////////////
//save gs_S_vs_n
/////////////////////////////////////////////////
 void gs_S_vs_n(int N, double chi, int n_1, int n_2, int interval, const std::string& filename) {
    int num_n= (n_2-n_1)/ interval +1;
    Eigen::VectorXd ns = Eigen::VectorXd::Zero(num_n);
    Eigen::VectorXd gsentropy = Eigen::VectorXd::Zero(num_n);

    for (int i = 0; i < num_n; ++i) {
        double n = n_1+ i*interval; 
        ns(i) = n;
        auto [E, coefficients] = lipkin_tot(N, chi);
        auto entropy = get_entropy(n, coefficients);
        gsentropy(i) = entropy(0); 
    }
    
    std::vector<Eigen::VectorXd> results = {ns, gsentropy};
	std::vector<std::string> parameters = {std::to_string(N),std::to_string(chi)};
	saveDataWithParams(results, filename , parameters);

	}

void gs_S_vs_chi_n(int N, double chi_1, double chi_2, int num_chi, int n_1, int n_2, int interval , const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Unable to open file" << std::endl;
        return;
    }
    // 生成chi和n的值
    std::vector<double> chis(num_chi);
    double chi_step = (chi_2 - chi_1) / (num_chi - 1);
    for (int i = 0; i < num_chi; ++i) {
        chis[i] = chi_1 + i * chi_step;
    }
    int num_n= (n_2-n_1)/ interval +1;
    std::vector<int> ns(num_n);
    for (int i = 0; i < num_n; ++i) {
        ns[i] = n_1 + i*interval;
    }
    // 写入文件
    file << N << std::endl;
    for (double chi : chis) {
        file << chi << " ";
    }
    file << std::endl;
    for (int n : ns) {
        file << n << " ";
    }
    file << std::endl;
    // 计算gs_S并写入文件
    for (int n : ns) {
        std::cout << "Process:"<< n  << '/'<< n_2 << std::endl;
        for (double chi : chis) {
        	auto [E, coefficients] = lipkin_tot(N, chi);
        	auto entropy = get_entropy(n, coefficients);
            file << entropy(0) << " ";
        }
        file << std::endl;
    }
    file.close();
}

void get_sys2(int n1, int n2, double ep1, double ep2, double V1, double V2, double V12){
    using namespace Eigen;
	using namespace std;
    double j1 = n1 / 2.0;
    double j2 = n2 / 2.0;
    int s1 = n1 + 1;
    int s2 = n2 + 1;

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
    if (debug_mode){
    ofstream outFile("../data/H_matrix.txt");
    if (!outFile.is_open()) {
        cerr << "Failed to open file for writing." << endl;
    }
    for (int i = 0; i < H.rows(); ++i) {
        for (int j = 0; j < H.cols(); ++j) {
            outFile << H(i, j);
            if (j != H.cols() - 1) outFile << " ";  
        }
        outFile << "\n";  
    }
    outFile.close();
    }
    
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
   
}


QuantumSystem::QuantumSystem(int n1, int n2, double ep1, double ep2, double V1, double V2, double V12)
    : n1(n1), n2(n2), ep1(ep1), ep2(ep2), V1(V1), V2(V2), V12(V12) {
    s1 = n1 + 1;
    s2 = n2 + 1;
    m1.resize(s1);
    m2.resize(s2);
    j1 = n1 / 2.0;
    j2 = n2 / 2.0;
    for (int i = 0; i < s1; i++) {
        m1[i] = -n1 / 2.0 + i;
    }
    for (int i = 0; i < s2; i++) {
        m2[i] = -n2 / 2.0 + i;
    }
    H = MatrixXd::Zero(s1 * s2, s1 * s2);
}

void QuantumSystem::computeHamiltonian() {
    // Implement your Hamiltonian calculation here
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
    if (debug_mode){
    ofstream outFile("../data/H_matrix.txt");
    if (!outFile.is_open()) {
        cerr << "Failed to open file for writing." << endl;
    }
    for (int i = 0; i < H.rows(); ++i) {
        for (int j = 0; j < H.cols(); ++j) {
            outFile << H(i, j);
            if (j != H.cols() - 1) outFile << " ";  
        }
        outFile << "\n";  
    }
    outFile.close();
    }
    hamiltonianComputed = true;
    eigenSystemSolved = false;  // Any new computation of H requires re-solving the eigen system
    entropyComputed = false; 
}

void QuantumSystem::solveEigenSystem() {
    if (!hamiltonianComputed) {
        computeHamiltonian();  // Ensure Hamiltonian is computed
    }
    // Implement your eigen system solution here
    SelfAdjointEigenSolver<MatrixXd> eigensolver(H);
    if (eigensolver.info() != Success) abort();
    eigenvalues = eigensolver.eigenvalues();
    eigenvectors = eigensolver.eigenvectors();
    gs_E=eigenvalues(0);
    eigenSystemSolved = true;
    entropyComputed = false; 
}

void QuantumSystem::computeEntropy() {
if (!eigenSystemSolved) {
        solveEigenSystem();  // Ensure eigen system is solved
    }
    // Implement your entropy computation here
    vector<MatrixXd> reduced_rho_1(s1 * s2, MatrixXd::Zero(s1,s1));
    entropy = VectorXd::Zero(s1 * s2);
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
            	entropy(ind_e) -= eigvals[i] * std::log(eigvals[i]);
            }else {
                std::cerr << "Warning:n, ind_e: " << ind_e
                          << ", i: " << i << ", eigval: " << eigvals[i] << std::endl;
            }	
        }
    }
	gs_S=entropy(0);
    entropyComputed = true;
}

double QuantumSystem::getGroundStateEntropy() {
    if (!entropyComputed) {
        computeEntropy();  // Ensure entropy is computed
    }
    return entropy[0];  // Assuming entropy for ground state is the first element
}
double QuantumSystem::getGroundStateEnergy() {
    if (!eigenSystemSolved) {
        solveEigenSystem();  
    }
    return eigenvalues[0];   
}

const MatrixXd& QuantumSystem::getHamiltonian() {
    if (!hamiltonianComputed) {
        solveEigenSystem();  
    }
    return H;   
}
const VectorXd& QuantumSystem::getEigenvalues(){
    if (!eigenSystemSolved) {
        solveEigenSystem();  
    }
    return eigenvalues;   
}
const MatrixXd& QuantumSystem::getEigenvectors() {
    if (!eigenSystemSolved) {
        solveEigenSystem();  
    }
    return eigenvectors;   
}
const VectorXd& QuantumSystem::getEntropy() {
    if (!entropyComputed) {
        solveEigenSystem();  
    }
    return entropy;   
}
void QuantumSystem::outputResults(const string& filename) {
    // Implement your result output here
    ofstream file(filename);
    for (int i = 0; i < eigenvalues.size(); ++i) {
        file << eigenvalues[i] << " " << entropy[i] << endl;
    }
    file.close();
}
























