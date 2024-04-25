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
	// Calculate a certain system, save to "../data/Energy_Entropy.txt"
    get_system(100, 1.0, 20, "../data/Energy_Entropy.txt"); //int N, double chi, int n
	// Calculate gs_S_vs_chi, save to "../data/S_vs_chi.txt"
	gs_S_vs_chi(100, 20, 0, 3, 10, "../data/S_vs_chi.txt"); //int N, int n, double chi_1, chi_2, num_chi
	//
	gs_S_vs_n(100, 2, 45, 55, 2, "../data/S_vs_n.txt"); //int N, double chi, int n1 , n2, interval
	///////
	// input: int N double chi_1, chi_2, num_chi , int n1 , n2, interval, path
	gs_S_vs_chi_n(100, 0, 3, 7, 40, 60, 2, "../data/S_vs_chi_n.txt");

    return 0;
}




