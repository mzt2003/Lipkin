// main.cpp
#include <fstream>
#include <vector>
#include <iostream> 
#include "functions.h"
#include <fstream>

int main() {
    int N = 500;       
    double chi = 1; 
    int n= 125;

    auto [E, coefficients] = lipkin_tot(N, chi);
    auto entropy = get_entropy(n, coefficients);

    std::cout << "Energy:\n" << E << std::endl;
    //std::cout << "coefficients:\n" << coefficients << std::endl;
	std::cout << "Entropy:\n" << entropy << std::endl;
	
	saveDataWithParams(E, entropy, "../data/data.txt", N, chi, n);

    return 0;
}




