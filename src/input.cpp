#include <iostream>
#include <fstream>
#include <limits>
#include <sstream>
#include <vector>
#include <cmath>
#include "data.h"
#include "input.h"

bool fInput(const char* inputFilePath, const char* setupFilePath, Data2D& data) { 

    std::ifstream inputFile(inputFilePath); 

    if (!inputFile) {
        std::cerr << "Failed to open the file. 1" << std::endl;
        return 1;
    }

    std::vector<double> numbers;
    std::string line;
    while (std::getline(inputFile, line)) {
        if (line[0] == '#') {
            continue;
        }
        std::istringstream iss(line);
        std::string token;
        while (iss >> token) {
            try {
                double num = std::stod(token);
                numbers.push_back(num);
            } catch (const std::exception& e) {
                std::cerr << "Error converting token to number: " << e.what() << std::endl;
            }
        }
    }

    data.dimX = numbers[0];
    data.dimY = numbers[1];
    data.nPoints = data.dimX * data.dimY;
    data.nFaces = (data.dimX - 1) * data.dimY + (data.dimY - 1) * data.dimX;
    data.nCells = (data.dimX - 1) * (data.dimY - 1);
    data.nhorizontalFaces = (data.dimX - 1) * data.dimY;
    data.nverticalFaces = (data.dimY - 1) * data.dimX;
    std::cout << "dimX: " << data.dimX << ", dimY: " << data.dimY << ", nPoints: " << data.nPoints << ", nFaces: " << data.nFaces << ", nCells: " << data.nCells << ", nhorizontalFaces: " << data.nhorizontalFaces << ", nverticalFaces: " << data.nverticalFaces << std::endl;
    data.points = new Point2D[data.nPoints];
    data.faces = new Face2D[data.nFaces];
    data.cells = new Cell2D[data.nCells];
    
    // point definition
    int j = 0; // point index
    for (int i = 2; i < data.nPoints*3; i += 3) {
        data.points[j].id = (int)numbers[i];
        data.points[j].x = numbers[i + 1];
        data.points[j].y = numbers[i + 2];
        // std::cout << "Point " << data.points[j].id << ": (" << data.points[j].x << ", " << data.points[j].y << ")" << std::endl;
        j++;
    }
    int p = 0; // cell index
    for (int i = 2 + data.nPoints*3; i < 2 + data.nPoints*3 + data.nCells*3; i+=3) {
        data.cells[p].id = (int)numbers[i];
        data.cells[p].bType_p = (int)numbers[i+1];
        if (data.cells[p].bType_p == 2){
            data.cells[p].bType_p = 0; 
        }
        // data.cells[p].bType_p = 0;
        if (data.cells[p].bType_p == DIRICHLET || data.cells[p].bType_p == INNERCELL){
            data.cells[p].p[INITIAL] = numbers[i+2];
        } else if (data.cells[p].bType_p == NEUMANN){
            data.cells[p].p[INITIAL] = numbers[i+2];
            // data.cells[p].g_p = numbers[i+3];
            data.cells[p].g_p = 0;
        }
        data.cells[p].alpha = 0;
        // std::cout << "Cell " << data.cells[p].id << ": alpha = " << data.cells[p].alpha << ", bType_sc = " << data.cells[p].bType_sc << ", bType_p = " << data.cells[p].bType_p << ", sc = " << data.cells[p].sc << ", p = " << data.cells[p].p[INITIAL] << std::endl;
        p++;
    }
    p = 0; // face index
    for (int i=2+data.nPoints*3+data.nCells*3; i < 2+data.nPoints*3+data.nCells*3+data.nFaces*4; i+=4) {
        data.faces[p].id = (int)numbers[i];
        data.faces[p].bType_u = (int)numbers[i+1];
        if (data.faces[p].bType_u == DIRICHLET || data.faces[p].bType_u == INNERCELL){
            data.faces[p].u[INITIAL] = numbers[i+2];
            data.faces[p].v[INITIAL] = numbers[i+3];
        } else if (data.faces[p].bType_u == NEUMANN){
            data.faces[p].u[INITIAL] = numbers[i+2];
            data.faces[p].v[INITIAL] = numbers[i+3];
            // data.faces[p].g_u = numbers[i+2];
            // data.faces[p].g_v = numbers[i+3];
            data.faces[p].g_u = 0;
            data.faces[p].g_v = 0;
        }
        // std::cout << "Face " << data.faces[p].id << ": bType_u = " << data.faces[p].bType_u << ", u = " << data.faces[p].u[INITIAL] << ", v = " << data.faces[p].v[INITIAL] << std::endl;
         p++;
    }
    inputFile.close();
    // Read second file in execution arguments

    std::ifstream secondInputFile(setupFilePath);

    if (!secondInputFile) {
        std::cerr << "Failed to open the file. 2" << std::endl;
        return 1;
    }

    std::vector<double> numbers_setup;
    std::string secondLine;
    while (std::getline(secondInputFile, secondLine)) {
        if (secondLine[0] == '#') {
            continue;
        }
        std::istringstream iss(secondLine);
        std::string token;
        while (iss >> token) {
            try {
                double num = std::stod(token);
                numbers_setup.push_back(num);
            } catch (const std::exception& e) {
                std::cerr << "Error converting token to number: " << e.what() << std::endl;
            }
        }
    }

    secondInputFile.close();
    data.maxTime = numbers_setup[0];
    data.dt = numbers_setup[1];
    data.maxIteration = numbers_setup[2];
    data.mode = numbers_setup[3];
    data.pecFunc = numbers_setup[4];
    data.velSolver = numbers_setup[5];
    data.fixedPressure = numbers_setup[6];
    data.presSolver = numbers_setup[7];
    std::cout << "maxTime: " << data.maxTime << std::endl;
    std::cout << "dt: " << data.dt << std::endl;
    if (data.dt == 0) {
        std::cerr << "Time step cannot be zero." << std::endl;
        return 1;
    }
    if (data.maxTime == 0) {
        std::cerr << "Maximum time cannot be zero." << std::endl;
        return 1;
    }
    double ratio = data.maxTime / data.dt;
    if (std::abs(ratio - std::round(ratio)) > std::numeric_limits<double>::epsilon()) {
        std::cerr << "maxTime is not a multiple of dt." << std::endl;
        return 1;
    }
    data.timeStep = 0;
    data.rho1 = numbers_setup[8];
    data.eta1 = numbers_setup[9];
    data.rho2 = numbers_setup[10];
    data.eta2 = numbers_setup[11];
    data.xInertiaDamper = numbers_setup[12];
    data.yInertiaDamper = numbers_setup[12];
    data.alpha_p_relax = numbers_setup[13];
    data.alpha_u_relax = numbers_setup[14];
    data.alpha_v_relax = numbers_setup[14];
    data.u_top = numbers_setup[15];
    data.u_bottom = numbers_setup[16];
    data.v_left = numbers_setup[17];
    data.v_right = numbers_setup[18];
    data.isThereGravity = numbers_setup[19];


    return 0;
}