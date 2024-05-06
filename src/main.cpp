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
    //get_system(30, 1.0, 20, "../data/sys1_Energy_Entropy.txt"); //int N, double chi, int n
    
	// Calculate gs_S_vs_chi, save chi and gsS to "../data/S_vs_chi.txt"
	//gs_S_vs_chi(30, 20, 0, 3.5, 100, "../data/sys1_S_vs_chi.txt"); //int N, int n, double chi_1, chi_2, num_chi
	
	//save n and gsS
	//gs_S_vs_n(100, 2, 10, 90, 2, "../data/sys1_S_vs_n.txt"); //int N, double chi, int n1 , n2, interval
	
	// save chi, n and gsS, input: int N double chi_1, chi_2, num_chi , int n1 , n2, interval, path
	//gs_S_vs_chi_n(100, 0, 3.5, 5, 40, 60, 5, "../data/sys1_S_vs_chi_n.txt");
	
	QuantumSystem qs2(20, 10, 1, 1, 0.01/29.0, 0.01/29.0, 10/29.0);
    qs2.outputResults("../data/sys2_Energy_Entropy_.txt");
    
     double v_1=0/30.0;
     double v_2=1/7.5 ;
     int num_v=4 ;
     double v12_1=0/30.0;
     double v12_2= 1/2.0;
     int num_v12=4 ;
     const std::string& filename="../data/sys2_S_v_v12.txt";
     const std::vector<double>& v1s={0.01/30.0,0.03/30.0,0.1/30.0,0.3/30.0,1/30.0,3/30.0};
      const std::vector<double>& v2s={0.01/30.0,0.03/30.0,0.1/30.0,0.3/30.0,1/30.0,3/30.0};
     //n1,n2,ep1,ep2, double v_1,double v_2, int num_v, double v12_1 double v12_2, int num_v12, const std::string& filename
    //write_sys2_gsS_v_v12(20, 10, 1, 1,  v_1,v_2, num_v, v12_1 ,v12_2, num_v12, filename);
    //write_sys2_gsS_v12(20, 10, 1, 1, 3/30.0, 3/30.0, v12_1, v12_2, 100, "../data/sys2_S_v12.txt");
    //write_sys2_gsS_v12(20, 10, 1, 1, 0.3/30.0, 0.3/30.0, v12_1, v12_2, 100, "../data/sys2_S_v12.txt");
    //write_sys2_gsS_v12(20, 10, 1, 1, 10/30.0, 10/30.0, v12_1, v12_2, 100, "../data/sys2_S_v12.txt");
	//write_sys2_gsS_v12_multi(20,10,1,1,v1s,v2s,v12_1, v12_2, 100, "../data/sys2_S_v12_multi.txt") ;
    //write_sys2_h1h2_v12(20, 10, 1, 1, 0.00/30.0, 0.00/30.0, v12_1, v12_2, 100, "../data/sys2_h1h2_v12.txt");

    return 0;
}




