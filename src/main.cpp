// main.cpp
#include <fstream>
#include <vector>
#include <iostream> 
#include "functions.h"
#include <Eigen/Dense>
#include <cmath>

using namespace Eigen;
using namespace std;
bool debug_mode=0;

int main() {
	// Calculate a certain system, save E and S
    //get_system(30, 2.5, 20, "../data/sys1_Energy_Entropy.txt"); //int N, double chi, int n
    
	// Calculate gs_S_vs_chi, save chi and gsS to "../data/S_vs_chi.txt"
	//gs_S_vs_chi(40, 20, 0, 3, 30, "../data/sys1_S_vs_chi.txt"); //int N, int n, double chi_1, chi_2, num_chi
	
	//save n and gsS
	//gs_S_vs_n(100, 2, 10, 90, 2, "../data/sys1_S_vs_n.txt"); //int N, double chi, int n1 , n2, interval
	
	// save chi, n and gsS, input: int N double chi_1, chi_2, num_chi , int n1 , n2, interval, path
	//gs_S_vs_chi_n(100, 0, 3.5, 5, 40, 60, 5, "../data/sys1_S_vs_chi_n.txt");
	
	QuantumSystem qs2(20, 80, 1.0, 1.0, 1/99.0, 1/99.0, 0.03/99.0);
    qs2.outputResults("../data/sys2_Energy_Entropy_.txt");
    //qs2.compute_gs_overlap1();
    qs2.outputE1soverlap1("../data/E1_overlap1_.txt");
    qs2.outputE2soverlap2("../data/E2_overlap2_.txt");
    double a, b;
    a=qs2.get_gs_temp1();
    b=qs2.get_gs_temp2();
    VectorXd SS, hh1;
    SS=qs2.getEntropy();
    hh1=qs2.getH1_expected();
    std::cout<<"a: "<<a<<"b: "<<b <<std::endl;
    std::cout<<"deltaSSSS: "<<SS(1)-SS(0)<<"deltahhhh1: "<<hh1(1)-hh1(0) <<std::endl;
    
     double v_1=0.0001/30.0;
     double v_2=1/15.0 ;
     int num_v=10 ;
     double v12_1=0/30.0;
     double v12_2= 1/2.0;
     int num_v12=40 ;
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
	 //write_sys2_temp_v(50, 50, 1, 1,  v_1, v_2,  num_v, "../data/sys2_temp_v.txt"); 

    return 0;
}




