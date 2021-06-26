#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <cassert>
#include <string>
#include <sstream>
#define M_PI 3.14159265358979323846 /* pi */

struct MeshPoint {
    double x;
    double y;
};
// function that reads input parameters
struct Mesh {
    //FSMACH (mach number in infinity), gamma = 1.4, EPSE (smoothing coefficient- around 0.06), deltaT
    double FSMACH;              // = 0.9;
    double gamma;               // = 1.4;
    double EPSE;                // = 0.06;
    double deltaT;              // = 1;

    const int cols;             // = 121;
    const int rows;             // = 41;

    std::string xFilename;
    std::string yFilename;
    double* x_arr;
    double* y_arr;

    Mesh(const char* paramFilename = nullptr) : FSMACH(0.9), gamma(1.4), EPSE(0.06), deltaT(1), cols(121), rows(41), 
        xFilename("x.csv"), yFilename("y.csv") {
        if (paramFilename) {
            parseParams(paramFilename);
        }
        x_arr = parseCSV(xFilename.c_str());
        y_arr = parseCSV(yFilename.c_str());
    }

    ~Mesh() {
        delete x_arr;
        delete y_arr;
    }

    void parseParams(const char* filename) {
        std::ifstream  data(filename);
        if (!data) {
            std::cerr << "failed reading params file" << filename << " you stupid" << std::endl;
        }
        std::cout << "great! you are reading params " << filename << std::endl;
        while (!data.eof()) {
            std::string keyword;
            data >> keyword;
            if (!data.good()) {
                break;
            }
            if (keyword == "FSMACH") { data >> FSMACH; continue; }
            if (keyword == "gamma") { data >> gamma; continue; }
            if (keyword == "EPSE") { data >> EPSE; continue; }
            if (keyword == "deltaT") { data >> deltaT; continue; }
            std::cerr << "unknowm params name " << keyword << std::endl;
            exit(1);
        }
    }

    // function that reads csv files (for the grid)
    double* parseCSV(const char* filename) {
        double* arr = (double*)malloc(rows * cols * sizeof(double));

        std::ifstream  data(filename);
        std::string line;
        int row = 0;
        while (!data.bad() && std::getline(data, line))
        {
            std::stringstream lineStream(line);
            int col = 0;
            while (col < cols)
            {
                lineStream >> arr[row * cols + col];
                char c = ' ';
                do { lineStream >> c; } while (c != ',' && !lineStream.eof());
                if (col == 5 && row == 9) {
                    std::cout << "row = " << row << ", " << "col = " << col << " => " << arr[row * cols + col] << std::endl;
                }
                col++;
            }
            assert(col == cols);
            row++;
        }
        assert(row == rows);
        return arr;
    };

    double& x(int i, int j) { return x_arr[offset2D(i, j)]; }
    double& y(int i, int j) { return y_arr[offset2D(i, j)]; }

    int offset2D(int i, int j) {
        if (!(i > 0 && i <= rows)) {
            std::cout << "got bad i = " << i << std::endl;
            assert(true);
        }
        if (!(j > 0 && j <= cols)) {
            std::cout << "got bad j = " << j << std::endl;
            assert(true);
        }
        return (j - 1) * rows + (i - 1);
    }

    void solve() {

    }
};



// main
int main(int argc, const char* argv[]) {
    if (argc > 1) {
        for (int i = 1; i < argc; i++) {
            Mesh mesh(argv[i]);
            mesh.solve();
        }
    }
    else {
        Mesh mesh;
        mesh.solve();
    }
    return 0;
}