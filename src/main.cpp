#include <fstream>
#include <vector>

int main() {
    std::vector<int> x(10), y(10);
    for (int i = 0; i < 10; ++i) {
        x[i] = i;
        y[i] = i;
    }

    std::ofstream outfile("data/data.txt");
    for (int i = 0; i < 10; ++i) {
        outfile << x[i] << " " << y[i] << std::endl;
    }
    outfile.close();

    return 0;
}
