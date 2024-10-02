#include <iostream>
#include <vector>
#include "data.h"
#include "input.h"
#include "setup.h"
#include "solve.h"
#include "output.h"
#include "plots.h"
#include "interface.h"
#include "alphaAdvection.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>
#include <filesystem>
#include <string>
#include <ctime>

int main(int argc, char** argv) {
    auto start = std::chrono::high_resolution_clock::now();
    Data2D data;
    const char* inputFilePath = argv[1];
    std::cout << "Input file path: " << inputFilePath << std::endl;
    const char* setupFilePath = argv[2];
    std::cout << "Setup file path: " << setupFilePath << std::endl;
    const char* alphaFilePath = argv[3];
    if (argc < 4) {
        std::cout << "No alpha file path provided" << std::endl;
    } else {
        const char* alphaFilePath = argv[3];
        std::cout << "Alpha file path: " << alphaFilePath << std::endl;
    }

    bool input = fInput(inputFilePath, setupFilePath, data);
    bool setup = fSetup(alphaFilePath, data);
    std::string solve = fSolve(data);
    initResidualPlot(data, solve);
    outerLoopResidualPlot(data, solve);
    // jsonOutput(data);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
    std::cout << "Runtime duration: "
         << duration.count() << " seconds" << std::endl;
    return 0;
}

// int main(int argc, char** argv) {
//     auto start = std::chrono::high_resolution_clock::now();
//     std::srand(std::time(0));
//     Data2D data;
//     const char* inputFilePath = argv[1];
//     const char* setupFilePath = argv[2];

//     bool input = fInput(inputFilePath, setupFilePath, data);
//     bool setup = fSetup(data);
//     std::string directoryName = createDirectory();
//     reconstructInterfaceLines(data);
//     unstaggerGrid(data, INITIAL);
//     fOutputVTKframe(data, directoryName, INITIAL);
//     data.timeStep++;
//     advectAlpha(data);
//     fOutputVTKframe(data, directoryName, CORRECTED_2);


//     auto stop = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
//     std::cout << "Runtime duration: "
//          << duration.count() << " seconds" << std::endl;
//     return 0;
// }

