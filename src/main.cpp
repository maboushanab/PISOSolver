#include <iostream>
#include <vector>
#include "data.h"
#include "input.h"
#include "setup.h"
#include "solve.h"
#include "output.h"
#include "plots.h"
#include "interface.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>
#include <filesystem>
#include <string>
#include <ctime>


int main(int argc, char** argv) {
    auto start = std::chrono::high_resolution_clock::now();
    std::srand(std::time(0));
    Data2D data;
    const char* inputFilePath = argv[1];
    const char* setupFilePath = argv[2];

    bool input = fInput(inputFilePath, setupFilePath, data);
    bool setup = fSetup(data);
    std::string solve = fSolve(data);
    residualPlot(data, solve);
    // jsonOutput(data);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
    std::cout << "Runtime duration: "
         << duration.count() << " seconds" << std::endl;
    return 0;
}
int main(int argc, char** argv) {
    auto start = std::chrono::high_resolution_clock::now();
    std::srand(std::time(0));
    Data2D data;
    const char* inputFilePath = argv[1];
    const char* setupFilePath = argv[2];

    bool input = fInput(inputFilePath, setupFilePath, data);
    bool setup = fSetup(data);
    std::string directoryName = createDirectory();
    reconstructInterfaceLines(data);


    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
    std::cout << "Runtime duration: "
         << duration.count() << " seconds" << std::endl;
    return 0;
}

