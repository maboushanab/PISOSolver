#include <iostream>
#include <fstream>
#include <vector>
#include "data.h"
#include <filesystem>
#include <string>
#include <ctime>




std::string createDirectory() {
    // Get today's date
    std::time_t now = std::time(nullptr);
    std::tm* localTime = std::localtime(&now);
    int year = localTime->tm_year + 1900;
    int month = localTime->tm_mon + 1;
    int day = localTime->tm_mday;

    // Create directory name
    std::string directoryName = std::to_string(year) + "-" + std::to_string(month) + "-" + std::to_string(day);

    // Check if directory already exists
    int id = 1;
    std::string finalDirectoryName = directoryName;
    while (std::filesystem::exists("../out/" + finalDirectoryName)) {
        finalDirectoryName = directoryName + "_" + std::to_string(id);
        id++;
    }

    // Create the directory
    std::filesystem::create_directory("../out/" + finalDirectoryName);

    return finalDirectoryName;
}

bool fOutputVTKframe(Data2D& data, std::string finalDirectoryName, int step) {
    // Open the output file
    std::string timeStep = std::to_string(data.timeStep);
    if (!std::filesystem::exists("../out/" + finalDirectoryName + "/outParaview")){
        std::filesystem::create_directory("../out/" + finalDirectoryName + "/outParaview");
    }
    std::ofstream outputFile("../out/" + finalDirectoryName + "/outParaview/timeStep_" + timeStep + ".vtk");
    if (!outputFile) {
        std::cerr << "Error opening output file!" << std::endl;
        return 1;
    }

    // Write the header
    outputFile << "# vtk DataFile Version 2.0" << std::endl;
    outputFile << "CFD Solver Output" << std::endl;
    outputFile << "ASCII" << std::endl;
    outputFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

    // Write the points
    outputFile << "POINTS " << data.nPoints << " float" << std::endl;
    for (int i=0; i<data.nPoints; i++) {
        outputFile << data.points[i].x << " " << data.points[i].y << " " << 0 << std::endl;
    }

    // Write the cells
    outputFile << "CELLS " << data.nCells << " " << (data.nCells * 5) << std::endl;
    for (int i=0; i<data.nCells; i++) {
        outputFile << "4 " << data.cells[i].points[0]->id << " " << data.cells[i].points[1]->id << " " << data.cells[i].points[2]->id << " " << data.cells[i].points[3]->id << std::endl;
    }

    // Write the cell types
    outputFile << "CELL_TYPES " << data.nCells << std::endl;
    for (int i = 0; i < data.nCells; ++i) {
        outputFile << "10" << std::endl; // Assuming all cells are tetrahedrons
    }

    // Write the cell data (alpha)
    outputFile << "CELL_DATA " << data.nCells << std::endl;
    outputFile << "SCALARS alpha float" << std::endl;
    outputFile << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < data.nCells; ++i) {
        outputFile << data.cells[i].alpha << std::endl;
    }
    // Write the cell data (p)
    outputFile << "SCALARS p float" << std::endl;
    outputFile << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < data.nCells; ++i) {
        outputFile << data.cells[i].p[step] << std::endl;
    }
    // Write the cell data (bType_p)
    outputFile << "SCALARS bType_p int" << std::endl;
    outputFile << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < data.nCells; ++i) {
        outputFile << data.cells[i].bType_p << std::endl;
    }
    // Write the cell data (u[step])
    outputFile << "VECTORS u float" << std::endl;
    for (int i = 0; i < data.nCells; ++i) {
        outputFile << data.cells[i].u[step] << " " << data.cells[i].v[step] << " " << 0 << std::endl;
    }
    //Write the cell data (n,m)
        outputFile << "VECTORS f(n,m) float" << std::endl;
    for (int i = 0; i < data.nCells; ++i) {
        outputFile << data.cells[i].interfaceLine.n << " " << data.cells[i].interfaceLine.m << " " << 0 << std::endl;
    }
    // Close the output file
    outputFile.close();

    return true;
}
