// main.cpp
#include <fstream>
#include <vector>
#include <iostream> 
#include "functions.h"
#include <fstream>

int main() {
    get_system(100, 1, 20); //int N, double chi, int n
	///////
	
	int N = 100;
    int n = 20;
	int num_chi = 301; 
    std::vector<double> chis(num_chi);
    std::vector<double> gsentropy(num_chi);

    for (int i = 0; i < num_chi; ++i) {
        double chi = 3.0 * i / (num_chi - 1); 
        chis[i] = chi;
        auto [E, coefficients] = lipkin_tot(N, chi);
        auto entropy = get_entropy(n, coefficients);
        gsentropy[i] = entropy(0); // 只取基态的entropy
    }
    std::ofstream file("../data/entropy_vs_chi.txt");
    if (file.is_open()) {
        for (int i = 0; i < num_chi; ++i) {
            file << chis[i] << " " << gsentropy[i] << std::endl;
        }
        file.close();
    } else {
        std::cerr << "Unable to open file" << std::endl;
        return -1;
    }

	
    return 0;
}




