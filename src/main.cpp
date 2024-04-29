// main.cpp
#include <fstream>
#include <vector>
#include <iostream> 
#include "functions.h"
#include <Eigen/Dense>
#include <cmath>


using namespace Eigen;
using namespace std;
bool debug_mode=1;

int main() {
	// Calculate a certain system, save E and S
    get_system(30, 1.0, 20, "../data/Energy_Entropy.txt"); //int N, double chi, int n
	// Calculate gs_S_vs_chi, save to "../data/S_vs_chi.txt"
	//gs_S_vs_chi(100, 20, 0, 3.5, 100, "../data/S_vs_chi.txt"); //int N, int n, double chi_1, chi_2, num_chi
	//
	//gs_S_vs_n(100, 2, 10, 90, 2, "../data/S_vs_n.txt"); //int N, double chi, int n1 , n2, interval
	///////
	// input: int N double chi_1, chi_2, num_chi , int n1 , n2, interval, path
	//gs_S_vs_chi_n(100, 0, 3.5, 300, 10, 90, 5, "../data/S_vs_chi_n.txt");
	//gs_S_vs_chi_n(100, 0, 3.5, 5, 40, 60, 5, "../data/S_vs_chi_n.txt");
	
    //get_sys2(10, 20, 1, 1, 1/29.0, 1/29.0, 1/29.0);



	QuantumSystem qs(20, 10, 1, 1, 1/29.0, 1/29.0, 1/29.0);
	//qs.temp();
    //qs.computeHamiltonian();
   // qs.solveEigenSystem();
   // qs.computeEntropy();
	//double a = qs.getGroundStateEntropy();
    //std::cout<< "entropyyyyyyyy" << a <<std::endl;
	qs.computeH1();
   // std::cout<<"entropy2445"<<qs.getEntropy()<<std::endl;
  //  std::cout<<"entropy65437"<<qs.getEigenvalues()<<std::endl;

    qs.outputResults("../data/Energy_Entropy33.txt");
    
    
    
    return 0;
}




