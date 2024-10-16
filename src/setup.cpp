#include <iostream>
#include <vector>
#include "data.h"
#include <fstream>

bool fSetup(const char* alphaFilePath, Data2D& data){
    std::cout << "Setup" << std::endl;
    std::vector<int> alphaPointIds;

    std::ifstream file(alphaFilePath);
    // put all integers line by line into a new entry in the vector
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << alphaFilePath << std::endl;
        return 1;
    }
    std::string line;
    // Read the file line by line
    while (std::getline(file, line)) {
        std::istringstream iss(line); // Create a string stream for each line
        int alphaPointId;

        // Extract integers from the line and add them to the vector
        while (iss >> alphaPointId) {
            alphaPointIds.push_back(alphaPointId);
        }
        if (alphaPointIds.size() == 0){
            std::cout << "Alpha Point list empty." << std::endl;
        }
    }

    file.close(); // Close the file when done
    

    // face definition
    int k = 0; // face index

    // horizontal faces
    for (int i = 0; i < data.nPoints; i++) {
        if ((i + 1) % data.dimX == 0) continue;
        data.faces[k].id = k;
        data.faces[k].points[0] = &data.points[i];
        data.faces[k].points[1] = &data.points[i + 1];
        k++;
    }
    
    // vertical faces
    for (int i = 0; i < data.nPoints - data.dimX; i++) {
        data.faces[k].id = k;
        data.faces[k].points[0] = &data.points[i];
        data.faces[k].points[1] = &data.points[i + data.dimX];
        k++;
    }

    // cells definition
    int l = 0; // cell index
    for (int i = 0; i < data.nPoints - data.dimX; i++) {
        if ((i + 1) % data.dimX == 0) continue;
        data.cells[l].id = l;
        data.cells[l].points[0] = &data.points[i];
        data.cells[l].points[1] = &data.points[i + 1];
        data.cells[l].points[2] = &data.points[i + data.dimX];
        data.cells[l].points[3] = &data.points[i + data.dimX + 1];
        data.cells[l].alpha_prev = 0;
        l++;
    }
    data.isThereAlpha = false;
    for (int i = 0; i < data.nCells; i++) {
        int pointCounter = 0;
        double p = 0;
        for (int k = 0; k < 4; k++){
            for (int j = 0; j<alphaPointIds.size(); j++) {
                if (data.cells[i].points[k]->id == alphaPointIds[j]) {
                    pointCounter++;
                }
            }
        }
        if (pointCounter == 4){
            data.cells[i].alpha = 1;
            data.isThereAlpha = true;
        }
        if (pointCounter == 2){ 
            data.cells[i].alpha = 0.75;
            data.isThereAlpha = true;
        }
        if (pointCounter == 1){
            data.cells[i].alpha = 0.5;
            data.isThereAlpha = true;
        }
    }
        

    // face-cell connection
    int m=0; // row index
    int n=0; // column index
    for (int i = 0; i < data.nFaces; i++) {
        if (i < data.nhorizontalFaces) {   // horizontal faces
            // check if face is on the bottom boundary
            if (i < data.dimX - 1) {
                data.faces[i].neighCells[0] = nullptr;
                data.faces[i].neighCells[1] = &data.cells[i];
                // std::cout << "Face " << data.faces[i].id << " connected to Cell " << data.faces[i].neighCells[1]->id << std::endl;
            // check if face is in the middle of the grid (horizontal faces)
            } else if (i <= ((data.dimX - 1) * data.dimY) - data.dimX ) {
                data.faces[i].neighCells[0] = &data.cells[i - (data.dimX - 1)];
                data.faces[i].neighCells[1] = &data.cells[i];
                // std::cout << "Face " << data.faces[i].id << " connected to Cell " << data.faces[i].neighCells[0]->id << " and Cell " << data.faces[i].neighCells[1]->id << std::endl;
            // check if face is on the top boundary
            } else if (i < (data.dimX - 1) * data.dimY) {
                data.faces[i].neighCells[0] = &data.cells[i - (data.dimX - 1)];
                data.faces[i].neighCells[1] = nullptr;
                // std::cout << "Face " << data.faces[i].id << " connected to Cell " << data.faces[i].neighCells[0]->id << std::endl;
        }} else { // vertical faces
            // check if face is on the left boundary
            if (i % (data.nhorizontalFaces + data.dimX * m)  == 0) {
                data.faces[i].neighCells[0] = nullptr;
                data.faces[i].neighCells[1] = &data.cells[m * (data.dimX - 1)];
                // std::cout << "Face " << data.faces[i].id << " connected to Cell " << data.faces[i].neighCells[1]->id << std::endl;
            // check if face is on the right boundary
            } else if (i % (data.nhorizontalFaces + data.dimX * (m + 1) - 1) == 0) {
                data.faces[i].neighCells[0] = &data.cells[data.dimX - 2  + m * (data.dimX - 1)];
                data.faces[i].neighCells[1] = nullptr;
                m++; // increment m (row index)
                n = 0; // reset n (column index)
                // std::cout << "Face " << data.faces[i].id << " connected to Cell " << data.faces[i].neighCells[0]->id << std::endl;
            // check if face is in the middle of the grid (vertical faces)
            } else {
                data.faces[i].neighCells[0] = &data.cells[m * (data.dimX - 1) + n];
                data.faces[i].neighCells[1] = &data.cells[m * (data.dimX - 1) + 1 + n];
                n++; // increment n (column index)
                // std::cout << "Face " << data.faces[i].id << " connected to Cell " << data.faces[i].neighCells[0]->id << " and Cell " << data.faces[i].neighCells[1]->id << std::endl;
            }
        }
    }

    // cell-face connection
    int o=0; // face index
    for (int i = 0; i < data.nCells; i++) {
        if (i % (data.dimX - 1) == 0 && i != 0) {
            o++; // increment o (skip vertical face index at the right boundary of the grid)
        }
        data.cells[i].faces[0] = &data.faces[i];
        data.cells[i].faces[1] = &data.faces[i + data.dimX - 1];
        data.cells[i].faces[2] = &data.faces[data.nhorizontalFaces + o];
        data.cells[i].faces[3] = &data.faces[data.nhorizontalFaces + o + 1];
        o++; // increment o (face index)
        // std::cout << "Cell " << data.cells[i].id << " connected to Face " << data.cells[i].faces[0]->id << ", Face " << data.cells[i].faces[1]->id << ", Face " << data.cells[i].faces[2]->id << ", Face " << data.cells[i].faces[3]->id << std::endl;
    }
    
    //cell-cell connection (go through all cells and assign each neighbering cell of the connected faces as neighbouring cell)
    for (int i = 0; i < data.nCells; i++) {
        for (int j = 0; j < 4; j++) {
            if (data.cells[i].faces[j]->neighCells[0] != nullptr) {
                if (data.cells[i].faces[j]->neighCells[0]->id != data.cells[i].id) {
                    data.cells[i].neighCells[j] = data.cells[i].faces[j]->neighCells[0];
                } else {
                    data.cells[i].neighCells[j] = data.cells[i].faces[j]->neighCells[1];
                }
            } else {
                data.cells[i].neighCells[j] = nullptr;
            }
        }

        // std::cout << "Cell " << data.cells[i].id << " connected to Cell ";
        // for (int j = 0; j < 4; j++) {
        //     if (data.cells[i].neighCells[j] != nullptr) {
        //         std::cout << data.cells[i].neighCells[j]->id << ", ";
        //     }
        //     else {
        //         std::cout << "None, ";
        //     }
        // }
        // std::cout << std::endl;

    }
    //point-face connection
    for (int i = 0; i < data.nFaces; i++) {
        Face2D* curFace = &data.faces[i];
        if (i < data.nhorizontalFaces) {
            curFace->points[UP]->faces[SOUTH] = curFace;
            curFace->points[DOWN]->faces[NORTH] = curFace;
        } else {
            curFace->points[LEFT]->faces[EAST] = curFace;
            curFace->points[RIGHT]->faces[WEST] = curFace;
        }
    }
    // face centers, lengths and numFlux intialisation
    for (int i = 0; i < data.nFaces; i++){
        data.faces[i].x = (data.faces[i].points[0]->x + data.faces[i].points[1]->x) / 2;
        data.faces[i].y = (data.faces[i].points[0]->y + data.faces[i].points[1]->y) / 2;
        data.faces[i].dx = abs(data.faces[i].points[0]->x - data.faces[i].points[1]->x);
        data.faces[i].dy = abs(data.faces[i].points[0]->y - data.faces[i].points[1]->y);
    }
    // cell volumes and centers (Cartesian grid only)
    bool isInside = false;
    double p = 0;
    for (int i = 0; i < data.nCells; i++){
        data.cells[i].vol = data.cells[i].faces[0]->dx * data.cells[i].faces[2]->dy;
        data.cells[i].x = data.cells[i].faces[SOUTH]->x;
        data.cells[i].y = data.cells[i].faces[WEST]->y;
        data.cells[i].interfaceVectorLine.origin[0] = 0;
        data.cells[i].interfaceVectorLine.origin[1] = 0;
        data.cells[i].interfaceVectorLine.direction[0] = 0;
        data.cells[i].interfaceVectorLine.direction[1] = 0;
        if (data.cells[data.nCells - 1 - i].alpha == 1){
            if(!isInside){
                isInside = true;
                p += data.rho1*9.81*data.cells[data.nCells - 1 - i].faces[WEST]->dy;
            }
            data.cells[data.nCells - 1 - i].p[INITIAL] = p;
        } else {
            isInside = false;
        }
    }

    return true; 
}
