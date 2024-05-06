//      define the functions
#include "functions.h"
#include "wignerSymbols.h"  //contains the CG coefficients
#include <iostream>
#include <vector>
#include <Eigen/Dense> 
#include <cmath>
#include <fstream>
using namespace Eigen;
using namespace std;
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

void write_sys2_gsS_v_v12(int n1,int n2,double ep1,double ep2, double v_1,double v_2, int num_v, double v12_1, double v12_2, int num_v12, const std::string& filename){
     double N=n1+n2;
       std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Unable to open file" << std::endl;
    }
    // 生成v,v12的值
    std::vector<double> vs(num_v);
    double v_step = (v_2 - v_1) / (num_v - 1);
    for (int i = 0; i < num_v; ++i) {
        vs[i] = v_1 + i * v_step;
    }
    std::vector<double> v12s(num_v12);
    double v12_step = (v12_2 - v12_1) / (num_v12 - 1);
    for (int i = 0; i < num_v12; ++i) {
        v12s[i] = v12_1 + i * v12_step;
                std::cout<<"v12s"<<v12s[i]<<std::endl;
    }

    // 写入文件
    file << N << std::endl;
    for (double v : vs) {
        file << v*N << " ";
    }
    file << std::endl;
    for (double v12 : v12s) {
        file << v12*N << " ";
    }
    file << std::endl;
    // 计算gs_S并写入文件
    for (double v12 : v12s) {
        std::cout << "Process:"<< v12  << '/'<< v12_2 << std::endl;
        for (double v : vs) {       	
        	QuantumSystem qs3(n1, n2, ep1,ep2, v,v,v12);
        	double a = qs3.getGroundStateEntropy();
            file << a << " ";
        }
        file << std::endl;
    }
    file.close();
}

void write_sys2_gsS_v12(int n1,int n2,double ep1,double ep2, double v_1,double v_2, double v12_1, double v12_2, int num_v12, const std::string& filename){
    double N=n1+n2;
    Eigen::VectorXd v12s = Eigen::VectorXd::Zero(num_v12);
    Eigen::VectorXd gsentropy = Eigen::VectorXd::Zero(num_v12);
	double v12_step = (v12_2 - v12_1) / (num_v12 - 1);
    for (int i = 0; i < num_v12; ++i) {
        double v12 = v12_1 + i * v12_step;
        v12s[i]= v12;
        QuantumSystem qs3(n1, n2, ep1,ep2, v_1,v_2,v12);
        gsentropy[i] = qs3.getGroundStateEntropy();
    }
    std::vector<Eigen::VectorXd> results = {v12s, gsentropy};
	std::vector<std::string> parameters = {std::to_string(n1),std::to_string(n2),std::to_string(ep1),std::to_string(ep2),std::to_string(v_1),std::to_string(v_2)};
	saveDataWithParams(results,filename , parameters);
}

void write_sys2_h1h2_v12(int n1,int n2,double ep1,double ep2, double v_1,double v_2, double v12_1, double v12_2, int num_v12, const std::string& filename){
    double N=n1+n2;
    Eigen::VectorXd v12s = Eigen::VectorXd::Zero(num_v12);
    Eigen::VectorXd h1 = Eigen::VectorXd::Zero(num_v12);
    Eigen::VectorXd h2 = Eigen::VectorXd::Zero(num_v12);
	double v12_step = (v12_2 - v12_1) / (num_v12 - 1);
    for (int i = 0; i < num_v12; ++i) {
        double v12 = v12_1 + i * v12_step;
        v12s[i]= v12;
        QuantumSystem qs3(n1, n2, ep1,ep2, v_1,v_2,v12);
        h1[i] = qs3.getH1_expected()[0];
        h2[i] = qs3.getH2_expected()[0];
     
    }
    std::vector<Eigen::VectorXd> results = {v12s, h1, h2};
	std::vector<std::string> parameters = {std::to_string(n1),std::to_string(n2),std::to_string(ep1),std::to_string(ep2),std::to_string(v_1),std::to_string(v_2)};
	saveDataWithParams(results,filename , parameters);
}

void write_sys2_gsS_v12_multi(int n1, int n2, double ep1, double ep2,
                              const std::vector<double>& v1s, const std::vector<double>& v2s,
                              double v12_1, double v12_2, int num_v12, const std::string& filename) {
    double N = n1 + n2;
    Eigen::VectorXd v12s = Eigen::VectorXd::Zero(num_v12);
    std::vector<Eigen::VectorXd> gsentropies(v1s.size() * v2s.size());
    double v12_step = (v12_2 - v12_1) / (num_v12 - 1);
    int count = 0;

        for (size_t k = 0; k < v2s.size(); ++k) {
            Eigen::VectorXd gsentropy = Eigen::VectorXd::Zero(num_v12);
            for (int i = 0; i < num_v12; ++i) {
                double v12 = v12_1 + i * v12_step;
                v12s[i] = v12;  
                QuantumSystem qs3(n1, n2, ep1, ep2, v1s[k], v2s[k], v12);
                gsentropy[i] = qs3.getGroundStateEntropy();
            }
            gsentropies[count++] = gsentropy;
        }
    
    std::cout<<"v12"<<v12s<<std::endl;
    std::vector<Eigen::VectorXd> results = {v12s};
    results.insert(results.end(), gsentropies.begin(), gsentropies.end());
    std::vector<std::string> parameters = {std::to_string(n1), std::to_string(n2), std::to_string(ep1), std::to_string(ep2)};
    for (double v : v1s) parameters.push_back(std::to_string(v));
    for (double v : v2s) parameters.push_back(std::to_string(v));

    saveDataWithParams(results, filename, parameters);
}

QuantumSystem::QuantumSystem(int n1, int n2, double ep1, double ep2, double V1, double V2, double V12)
    : n1(n1), n2(n2), ep1(ep1), ep2(ep2), V1(V1), V2(V2), V12(V12),s1(n1 + 1), s2(n2 + 1),  // 确保这里s1和s2先被初始化
       reduced_rho_1(s1 * s2),reduced_rho_2(s1 * s2) {
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
    H1= MatrixXd::Zero(s1 , s1 );
    H2= MatrixXd::Zero(s2,  s2);
    H1_expected=VectorXd(s1*s2);
    H2_expected=VectorXd(s1*s2);
    for (Eigen::MatrixXd& matrix : reduced_rho_1) {
        matrix = MatrixXd::Zero(s1, s1);  
    }
    for (Eigen::MatrixXd& matrix : reduced_rho_2) {
        matrix = MatrixXd::Zero(s2, s2); 
    }
    entropy = VectorXd::Zero(s1 * s2);
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
}
void QuantumSystem::computeH1(){
	for (int i = 0; i < s1; ++i) {
        H1(i, i) = ep1 * m1[i];
    }
    for (int i = 0; i < s1 - 2; ++i) {
        double m = m1[i];
        double temp1 = -V1 / 2;
        double t2 = j1 * (j1 + 1) - m * (m + 1);
        double t3 = j1 * (j1 + 1) - (m + 1) * (m + 2);
        double tt = temp1 * std::sqrt(t2 * t3);
        H1(i, i + 2) = tt;
        H1(i + 2, i) = tt;
    }
    H1Computed = true;
    H1_expected_Computed = false;
}
void QuantumSystem::computeH2(){
	for (int i = 0; i < s2; ++i) {
        H2(i, i) = ep2 * m2[i];
    }

    for (int i = 0; i < s2 - 2; ++i) {
        double m = m2[i];
        double temp1 = -V2 / 2;
        double t2 = j2 * (j2 + 1) - m * (m + 1);
        double t3 = j2 * (j2 + 1) - (m + 1) * (m + 2);
        double tt = temp1 * std::sqrt(t2 * t3);
        H2(i, i + 2) = tt;
        H2(i + 2, i) = tt;
    }
    H2Computed= true;
    H2_expected_Computed = false;
}
void QuantumSystem::computeH1_expected(){
	if (!H1Computed) {
        computeH1();  
    }
    if (!reduced_rho_1_Computed){
    	computeEntropy();
    }
    for (int ii=0; ii< s1*s2; ++ii){
        Eigen::MatrixXd product = H1 * reduced_rho_1[ii];
    	H1_expected(ii) = product.trace();
    }
    H1_expected_Computed = true ;
}
void QuantumSystem::computeH2_expected(){
	if (!H2Computed) {
        computeH2();  
    }
    if (!reduced_rho_2_Computed){
    	computeEntropy();
    }
    for (int ii=0; ii< s1*s2; ++ii){
        Eigen::MatrixXd product = H2 * reduced_rho_2[ii];
    	H2_expected(ii) = product.trace();
    }
    H2_expected_Computed = true;
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
    for (int ii = 0; ii < s1 * s2; ++ii) {
        Map<MatrixXd> t1(eigenvectors.col(ii).data(), s2, s1);
        MatrixXd rho = t1.transpose() * t1;
        reduced_rho_1[ii] = rho;  
        reduced_rho_2[ii] = t1 * t1.transpose();  
        std::string filename = "../data/reduced_rho111_" + std::to_string(ii) + ".txt";
        std::ofstream file1(filename);
        if (file1.is_open()) {
         file1 << reduced_rho_1[ii] << "\n";
            file1.close();
        } else {
            std::cerr << "Unable to open file: " << filename << std::endl;
        }        
     }  
           
    for (int ind_e = 0; ind_e < s1*s2; ++ind_e) {
        VectorXd eigvals = SelfAdjointEigenSolver<MatrixXd>(reduced_rho_1[ind_e]).eigenvalues();
       // std::cout << "eigenvals "<< eigvals << std::endl; 
        for (int i = 0; i < eigvals.size(); ++i) {
        	if (eigvals[i] > 0) {
            	entropy(ind_e) -= eigvals[i] * std::log(eigvals[i]);
            }else {
                //std::cerr << "Warning:n, ind_e: " << ind_e
                 //         << ", i: " << i << ", eigval: " << eigvals[i] << std::endl;
            }	
        }
    }
	gs_S=entropy(0);
    entropyComputed = true;
    reduced_rho_1_Computed= true ;
    reduced_rho_2_Computed= true ;
    H1_expected_Computed = false;
    H2_expected_Computed = false;

}

double QuantumSystem::getGroundStateEntropy() {
    if (!entropyComputed) {
        computeEntropy();  // Ensure entropy is computed
    }
    return gs_S;  // Assuming entropy for ground state is the first element
}
double QuantumSystem::getGroundStateEnergy() {
    if (!eigenSystemSolved) {
        solveEigenSystem();  
    }
    return eigenvalues[0];   
}
const MatrixXd& QuantumSystem::getHamiltonian() {
    if (!hamiltonianComputed) {
        computeHamiltonian();  
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
        computeEntropy();  
    }
    return entropy;   
}
const VectorXd& QuantumSystem::getH1_expected(){
	if (!H1_expected_Computed) {
        computeH1_expected();  
    }
    return H1_expected; 
}
const VectorXd& QuantumSystem::getH2_expected(){
	if (!H2_expected_Computed) {
        computeH2_expected();  
    }
    return H2_expected; 
}

void QuantumSystem::outputResults(const string& filename) {
    if (!H1_expected_Computed) {
        computeH1_expected();  
    }

    if (!H2_expected_Computed) {
        computeH2_expected();  
    }
    //  result output 
    std::vector<Eigen::VectorXd> results = {eigenvalues,entropy,H1_expected, H2_expected};
	std::vector<std::string> parameters = {std::to_string(n1),std::to_string(n2),std::to_string(ep1),std::to_string(ep2),std::to_string(V1),std::to_string(V2),std::to_string(V12)};
	saveDataWithParams(results, filename , parameters);
}
























