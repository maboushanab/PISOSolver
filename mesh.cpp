#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

void processFile(std::ifstream &file, std::unordered_map<std::string, double> &keywordValues) {
    std::string line;
    std::string currentKeyword;

    while (std::getline(file, line)) {
        // Skip lines starting with #
        if (!line.empty() && line[0] == '#') {
            continue;
        }

        // Check if the line is one of the predefined keywords
        if (keywordValues.find(line) != keywordValues.end()) {
            currentKeyword = line;

            // Read the next line for the value
            if (std::getline(file, line)) {
                // Skip lines starting with #
                while (!line.empty() && line[0] == '#') {
                    if (!std::getline(file, line)) {
                        break;
                    }
                }
                try {
                    keywordValues[currentKeyword] = std::stod(line);
                } catch (const std::invalid_argument &e) {
                    std::cerr << "Invalid value for " << currentKeyword << ": " << line << std::endl;
                }
            }
        }
    }
}
//
//
// ASSUMPTION: GRID IS RECTANGULAR 
//

void cellsConstPressure(std::vector<std::string>& CellProp, int xDim, int yDim, double pressure, 
                        double left_sc, double left_p, double right_sc, double right_p, double top_sc, 
                        double top_p, double bottom_sc, double bottom_p) {
    int cellIndex = 0;
    for (int i = 0; i < (xDim - 1) * (yDim - 1); i++) {
        if (cellIndex >= (xDim - 1) * (yDim - 1) - xDim + 1) {
            CellProp.push_back(std::to_string(cellIndex) + " " + std::to_string(top_sc) + " " + std::to_string(top_p) + " " + std::to_string(pressure));
        } else if (cellIndex < xDim - 1) {
            CellProp.push_back(std::to_string(cellIndex) + " " + std::to_string(bottom_sc) + " " + std::to_string(bottom_p) + " " + std::to_string(pressure));
        } else if (cellIndex % (xDim - 1) == xDim - 2) {
            CellProp.push_back(std::to_string(cellIndex) + " " + std::to_string(right_sc) + " " + std::to_string(right_p) + " " + std::to_string(pressure));
        } else if (cellIndex % (xDim - 1) == 0) {
            CellProp.push_back(std::to_string(cellIndex) + " " + std::to_string(left_sc) + " " + std::to_string(left_p) + " " + std::to_string(pressure));
        } else {
            CellProp.push_back(std::to_string(cellIndex) + " 0 0 " + std::to_string(pressure));
        }
        cellIndex++;
    }
}

void cellsPressureGradient_y(std::vector<std::string>& cellProp, int xDim, int yDim, double pressure_up, double pressure_down, 
                            double left_sc, double left_p, double right_sc, double right_p, double top_sc, double top_p,
                            double bottom_sc, double bottom_p) {
    int cellIndex = 0;
    double p_grad_y = (pressure_down - pressure_up) / (yDim - 2);
    for (int i = 0; i < yDim - 1; i++) {
        double p_y = pressure_up + i * p_grad_y;
        for(int j = 0; j < xDim - 1; j++){
            if (cellIndex >= (xDim - 1) * (yDim - 1) - xDim + 1) {
                cellProp.push_back(std::to_string(cellIndex) + " " + std::to_string(top_sc) + " " + std::to_string(top_p) + " " + std::to_string(p_y));
            } else if (cellIndex < xDim - 1) {
                cellProp.push_back(std::to_string(cellIndex) + " " + std::to_string(bottom_sc) + " " + std::to_string(bottom_p) + " " + std::to_string(p_y));
            } else if (cellIndex % (xDim - 1) == xDim - 2) {
                cellProp.push_back(std::to_string(cellIndex) + " " + std::to_string(right_sc) + " " + std::to_string(right_p) + " " + std::to_string(p_y));

            } else if (cellIndex % (xDim - 1) == 0) {
                cellProp.push_back(std::to_string(cellIndex) + " " + std::to_string(left_sc) + " " + std::to_string(left_p) + " " + std::to_string(p_y));

            } else {
                cellProp.push_back(std::to_string(cellIndex) + " 0 0 " + std::to_string(p_y));
            }
            cellIndex++;
        }
    }

}

void cellsPressureGradient_x(std::vector<std::string>& cellProp, int xDim, int yDim, double pressure_left, double pressure_right,
                            double left_sc, double left_p, double right_sc, double right_p, double top_sc, double top_p,
                            double bottom_sc, double bottom_p) {
    int cellIndex = 0;
    double p_grad_x = (pressure_right - pressure_left) / (xDim - 2);
    for (int i = 0; i < yDim - 1; i++) {
        double p_x = pressure_left;
        for(int j = 0; j < xDim - 1; j++){
            if (cellIndex >= (xDim - 1) * (yDim - 1) - xDim + 1) {
                cellProp.push_back(std::to_string(cellIndex) + " " + std::to_string(top_sc) + " " + std::to_string(top_p) + " " + std::to_string(p_x));
            } else if (cellIndex < (xDim - 1)) {
                cellProp.push_back(std::to_string(cellIndex) + " " + std::to_string(bottom_sc) + " " + std::to_string(bottom_p) + " " + std::to_string(p_x));
            } else if (cellIndex % (xDim - 1) == xDim - 2) {
                cellProp.push_back(std::to_string(cellIndex) + " " + std::to_string(right_sc) + " " + std::to_string(right_p) + " " + std::to_string(p_x));
            } else if (cellIndex % (xDim - 1) == 0) {
                cellProp.push_back(std::to_string(cellIndex) + " " + std::to_string(left_sc) + " " + std::to_string(left_p) + " " + std::to_string(p_x));
            } else {
                cellProp.push_back(std::to_string(cellIndex) + " 0 0 " + std::to_string(p_x));
            }
            cellIndex++;
            p_x += p_grad_x;
        }
    }

}

void cellsConstanScalar(std::vector<std::string>& cellProp, int xDim, int yDim, double scalar) {
    for (int i = 0; i < cellProp.size(); i++) {
        cellProp[i] += " " + std::to_string(scalar);
    }
}

void cellsScalarGradient_y(std::vector<std::string>& cellProp, int xDim, int yDim, double scalar_up, double scalar_down) {
    int cellIndex = 0;
    double s_grad_y = (scalar_down - scalar_up) / (yDim - 2);
    for (int i = 0; i < yDim - 1; i++) {
        double s_y = scalar_up + i * s_grad_y;
        for(int j = 0; j < xDim - 1; j++){
            cellProp[cellIndex] += " " + std::to_string(s_y);
            cellIndex++;
        }
    }

}

void cellsScalarGradient_x(std::vector<std::string>& cellProp, int xDim, int yDim, double scalar_left, double scalar_right) {
    int cellIndex = 0;
    double s_grad_x = (scalar_right - scalar_left) / (xDim - 2);
    for (int i = 0; i < yDim - 1; i++) {
        double s_x = scalar_left;
        for(int j = 0; j < xDim - 1; j++){
            cellProp[cellIndex] += " " + std::to_string(s_x);
            cellIndex++;
            s_x += s_grad_x;
        }
    }

}

void createMesh(const std::unordered_map<std::string, double> &keywordValues) {
    if (keywordValues.at("dx") == 0 || keywordValues.at("dy") == 0) {
        std::cerr << "Invalid values for dx or dy." << std::endl;
        return;
    }
    // if (keywordValues.at("x") % keywordValues.at("dx") != 0 || keywordValues.at("y") % keywordValues.at("dy") != 0) {
    //     std::cerr << "x or y is not divisible by dx or dy." << std::endl;
    //     return;
    // }
    if (keywordValues.at("x") == 0 || keywordValues.at("y") == 0) {
        std::cerr << "Invalid values for x or y." << std::endl;
        return;
    }
    if (keywordValues.at("left_sc") <= 0 || keywordValues.at("left_sc") > 3) {
        std::cerr << "Invalid value for left_sc." << std::endl;
        return;
    }
        if (keywordValues.at("left_p") <= 0 || keywordValues.at("left_p") > 3) {
        std::cerr << "Invalid value for left_p." << std::endl;
        return;
    }
    if (keywordValues.at("right_sc") <= 0 || keywordValues.at("right_sc") > 3) {
        std::cerr << "Invalid value for right_sc." << std::endl;
        return;
    }
    if (keywordValues.at("right_p") <= 0 || keywordValues.at("right_p") > 3) {
        std::cerr << "Invalid value for right_p." << std::endl;
        return;
    }
    if (keywordValues.at("top_sc") <= 0 || keywordValues.at("top_sc") > 3) {
        std::cerr << "Invalid value for top_sc." << std::endl;
        return;
    }
    if (keywordValues.at("top_p") <= 0 || keywordValues.at("top_p") > 3) {
        std::cerr << "Invalid value for top_p." << std::endl;
        return;
    }
    if (keywordValues.at("bottom_sc") <= 0 || keywordValues.at("bottom_sc") > 3) {
        std::cerr << "Invalid value for bottom_sc." << std::endl;
        return;
    }
    if (keywordValues.at("bottom_p") <= 0 || keywordValues.at("bottom_p") > 3) {
        std::cerr << "Invalid value for bottom_p." << std::endl;
        return;
    }
    
    int xDim = keywordValues.at("x") / keywordValues.at("dx");
    int yDim = keywordValues.at("y") / keywordValues.at("dy");

    std::ofstream meshFile("in.mesh");

    if (meshFile.is_open()) {
        meshFile << "# point dimensions" << std::endl;
        meshFile << xDim << " " << yDim << std::endl;
        meshFile << std::endl;
        meshFile << "# points" << std::endl;
        meshFile << "# x y" << std::endl;
        int pointIndex = 0;
        for (int i = 0; i < yDim; i++) {
            for (int j = 0; j < xDim; j++) {
                meshFile << pointIndex << " " << j * keywordValues.at("dx") << " " << i * keywordValues.at("dy") << std::endl;
                pointIndex++;
            }
        }
        meshFile << std::endl;
        meshFile << "# cells" << std::endl;
        meshFile << "# boundary_sc boundary_p p sc" << std::endl;
        std::vector<std::string> cellProp;
        if (keywordValues.at("pressure_field_type") == 0) {
            cellsConstPressure(cellProp, xDim, yDim, keywordValues.at("pressure"), keywordValues.at("left_sc"), keywordValues.at("left_p"), keywordValues.at("right_sc"), keywordValues.at("right_p"), keywordValues.at("top_sc"), keywordValues.at("top_p"), keywordValues.at("bottom_sc"), keywordValues.at("bottom_p"));
        } else if (keywordValues.at("pressure_field_type") == 1){
            cellsPressureGradient_x(cellProp, xDim, yDim, keywordValues.at("pressure_left"), keywordValues.at("pressure_right"), keywordValues.at("left_sc"), keywordValues.at("left_p"), keywordValues.at("right_sc"), keywordValues.at("right_p"), keywordValues.at("top_sc"), keywordValues.at("top_p"), keywordValues.at("bottom_sc"), keywordValues.at("bottom_p"));
        } else if (keywordValues.at("pressure_field_type") == 2){
            cellsPressureGradient_y(cellProp, xDim, yDim, keywordValues.at("pressure_up"), keywordValues.at("pressure_down"), keywordValues.at("left_sc"), keywordValues.at("left_p"), keywordValues.at("right_sc"), keywordValues.at("right_p"), keywordValues.at("top_sc"), keywordValues.at("top_p"), keywordValues.at("bottom_sc"), keywordValues.at("bottom_p"));
        } else { 
            std::cerr << "Invalid value for pressure_field_type." << std::endl;
            return;
        }
        if (keywordValues.at("scalar_field_type") == 0) {
            cellsConstanScalar(cellProp, xDim, yDim, keywordValues.at("scalar"));
        } else if (keywordValues.at("scalar_field_type") == 1){
            cellsScalarGradient_x(cellProp, xDim, yDim, keywordValues.at("scalar_left"), keywordValues.at("scalar_right"));
        } else if (keywordValues.at("scalar_field_type") == 2){
            cellsScalarGradient_y(cellProp, xDim, yDim, keywordValues.at("scalar_up"), keywordValues.at("scalar_down"));
        } else { 
            std::cerr << "Invalid value for scalar_field_type." << std::endl;
            return;
        }
        for (int i = 0; i < cellProp.size(); i++) {
            meshFile << cellProp[i] << std::endl;
        }

        meshFile << std::endl;
        meshFile << "# faces" << std::endl;
        meshFile << "# boundary_vel u v" << std::endl;
        int faceIndex = 0;
        int m = 0; // row index	

        //DDUBLED BOUNDARY CONDITIONS FACES
        /*
        for (int i = 0; i < (xDim - 1) * yDim + (yDim - 1) * xDim; i++) {
            if (i > (xDim - 1) * yDim && i % xDim == 0 && i != 0) {
                m++;
            }
            if (i >= (xDim - 1) * yDim - 2 * xDim + 2 && i < (xDim - 1) * yDim) {	            //horizonal faces (bottom boundary)
                meshFile << faceIndex << " " << keywordValues.at("bottom_vel") << " " << keywordValues.at("bottom_u") << " " << keywordValues.at("bottom_v") << std::endl;

            } else if (i >= (xDim - 1) * yDim + (yDim - 1) * xDim - xDim) {                 //vertical faces (bottom boundary)
                meshFile << faceIndex << " " << keywordValues.at("bottom_vel") << " " << keywordValues.at("bottom_u") << " " << keywordValues.at("bottom_v") << std::endl;

            } else if (i < 2 * (xDim - 1)) {                                                    //horizontal faces (top boundary)  
                meshFile << faceIndex << " " << keywordValues.at("top_vel") << " " << keywordValues.at("top_u") << " " << keywordValues.at("top_v") << std::endl;

            } else if (i >= (xDim - 1) * yDim && i < (xDim - 1) * yDim + xDim) {            //vertical faces (top boundary)
                meshFile << faceIndex << " " << keywordValues.at("top_vel") << " " << keywordValues.at("top_u") << " " << keywordValues.at("top_v") << std::endl;

            } else if (i % (xDim - 1) == 0 && i <= (xDim - 1) * yDim) {                     //horizontal faces (left boundary)
                meshFile << faceIndex << " " << keywordValues.at("left_vel") << " " << keywordValues.at("left_u") << " " << keywordValues.at("left_v") << std::endl;

            } else if (i % ((xDim - 1) * yDim + xDim * m)  == 0 || i % ((xDim - 1) * yDim + xDim * m)  == 1) {                          //vertical faces (left boundary)
                meshFile << faceIndex << " " << keywordValues.at("left_vel") << " " << keywordValues.at("left_u") << " " << keywordValues.at("left_v") << std::endl;

            } else if (i % (xDim - 1) == xDim - 2 && i <= (xDim - 1) * yDim) {              //horizontal faces (right boundary)
                meshFile << faceIndex << " " << keywordValues.at("right_vel") << " " << keywordValues.at("right_u") << " " << keywordValues.at("right_v") << std::endl;

            } else if (i % ((xDim - 1) * yDim + xDim * m - 1) == 0 || i % ((xDim - 1) * yDim + xDim * m - 1) == ((xDim - 1) * yDim + xDim * m - 1) - 1) {               	    //vertical faces (right boundary)
                meshFile << faceIndex << " " << keywordValues.at("right_vel") << " " << keywordValues.at("right_u") << " " << keywordValues.at("right_v") << std::endl;

            } else {
                meshFile << faceIndex << " 0 0 0" << std::endl;
            }
            faceIndex++;
        }*/
        for (int i = 0; i < (xDim - 1) * yDim + (yDim - 1) * xDim; i++) {
            if (i > (xDim - 1) * yDim && i % xDim == 0 && i != 0) {
                m++;
            }
            if (i >= (xDim - 1) * yDim - xDim + 1 && i < (xDim - 1) * yDim) {	            //horizonal faces (bottom boundary)
                meshFile << faceIndex << " " << keywordValues.at("top_vel") << " " << keywordValues.at("top_u") << " " << keywordValues.at("top_v") << std::endl;

            } else if (i >= (xDim - 1) * yDim + (yDim - 1) * xDim - xDim) {                 //vertical faces (bottom boundary)
                meshFile << faceIndex << " " << keywordValues.at("top_vel") << " " << keywordValues.at("top_u") << " " << keywordValues.at("top_v") << std::endl;

            } else if (i < (xDim - 1)) {                                                    //horizontal faces (top boundary)  
                meshFile << faceIndex << " " << keywordValues.at("bottom_vel") << " " << keywordValues.at("bottom_u") << " " << keywordValues.at("bottom_v") << std::endl;

            } else if (i >= (xDim - 1) * yDim && i < (xDim - 1) * yDim + xDim) {            //vertical faces (top boundary)
                meshFile << faceIndex << " " << keywordValues.at("bottom_vel") << " " << keywordValues.at("bottom_u") << " " << keywordValues.at("bottom_v") << std::endl;

            } else if (i % (xDim - 1) == 0 && i <= (xDim - 1) * yDim) {                     //horizontal faces (left boundary)
                meshFile << faceIndex << " " << keywordValues.at("left_vel") << " " << keywordValues.at("left_u") << " " << keywordValues.at("left_v") << std::endl;

            } else if ((i + 1 -(xDim - 1) * yDim) % xDim == 1 && i > (xDim - 1) * yDim) {                          //vertical faces (left boundary)
                meshFile << faceIndex << " " << keywordValues.at("left_vel") << " " << keywordValues.at("left_u") << " " << keywordValues.at("left_v") << std::endl;

            } else if (i % (xDim - 1) == xDim - 2 && i <= (xDim - 1) * yDim) {              //horizontal faces (right boundary)
                meshFile << faceIndex << " " << keywordValues.at("right_vel") << " " << keywordValues.at("right_u") << " " << keywordValues.at("right_v") << std::endl;

            } else if ((i + 1 - (xDim - 1) * yDim) % xDim == 0 && i > (xDim - 1) * yDim) {      //vertical faces (right boundary)
                meshFile << faceIndex << " " << keywordValues.at("right_vel") << " " << keywordValues.at("right_u") << " " << keywordValues.at("right_v") << std::endl;

            } else {
                meshFile << faceIndex << " 0 0 0" << std::endl;
            }
            faceIndex++;
        }

    } else {
        std::cerr << "Unable to create mesh file." << std::endl;
    }
}

int main(int argc, char **argv) {
    std::ifstream meshSettingsFile(argv[1]);

    if (!meshSettingsFile) {
        std::cerr << "Unable to open file." << std::endl;
        return 1;
    }

    // Define a map to store values by keyword
    std::unordered_map<std::string, double> keywordValues = {
        {"x", 0},
        {"y", 0},
        {"dx", 0},
        {"dy", 0},
        {"left_sc", 0},
        {"left_p", 0},
        {"right_sc", 0},
        {"right_p", 0},
        {"top_sc", 0},
        {"top_p", 0},
        {"bottom_sc", 0},
        {"bottom_p", 0},
        {"left_vel", 0},
        {"right_vel", 0},
        {"top_vel", 0},
        {"bottom_vel", 0},
        {"left_u", 0},
        {"right_u", 0},
        {"top_u", 0},
        {"bottom_u", 0},
        {"left_v", 0},
        {"right_v", 0},
        {"top_v", 0},
        {"bottom_v", 0},
        {"pressure_field_type", 0},
        {"pressure", 0},
        {"pressure_left", 0},
        {"pressure_right", 0},
        {"pressure_up", 0},
        {"pressure_down", 0},
        {"scalar_field_type", 0},
        {"scalar", 0},
        {"scalar_left", 0},
        {"scalar_right", 0},
        {"scalar_up", 0},
        {"scalar_down", 0}
    };

    processFile(meshSettingsFile, keywordValues);
    createMesh(keywordValues);
    meshSettingsFile.close();

    return 0;
}
