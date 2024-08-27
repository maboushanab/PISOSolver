
#include "data.h"
#include "interface.h"

std::vector<double> calculateFluxLines(Data2D data, int cellId){
    std::vector<double> fluxLines;
    fluxLines.push_back(0);
    fluxLines.push_back(0);
    Cell2D curCell = data.cells[cellId];

    fluxLines[0] = data.dt * (curCell.faces[WEST]->u[CORRECTED_2] + curCell.faces[EAST]->u[CORRECTED_2]);
    fluxLines[1] = data.dt * (curCell.faces[SOUTH]->v[CORRECTED_2] + curCell.faces[NORTH]->v[CORRECTED_2]);

    return fluxLines;
}

void calculatePolygonArea(double n, double m, Eigen::Vector2d normalVec, Eigen::Vector2d p1, Eigen::Vector2d p2, Eigen::Vector2d p3, Eigen::Vector2d p4, std::vector<Eigen::Vector2d>& polygon, bool& isIntersecting){
    const std::vector<Eigen::Vector2d> cell = {
        p1, p2, p3, p4
    }; 
    std::vector<Eigen::Vector2d> cellEdit = {
        p1, p2, p3, p4
    }; 
    Eigen::Vector2d vertex;
    // create stack starting from p1 and seeing if the line intersects between two points
    std::vector<Eigen::Vector2d> intersectionPoints;
    std::vector<Eigen::Vector2d> polygon1;
    std::vector<Eigen::Vector2d> polygon2;

    if (x(cell[0][1], m, n) > cell[0](0) && x(cell[0][1], m, n) < cell[1](0)) { // top face
        vertex(0) = x(cell[0](1), m, n);
        vertex(1) = cell[0](1);
        intersectionPoints.push_back(vertex);
    }
    if (y(cell[1](0), m, n) > cell[2](1) && y(cell[1](0), m, n) < cell[1](1)) { // right face
        vertex(0) = cell[1](0);
        vertex(1) = y(cell[1](0), m, n);
        intersectionPoints.push_back(vertex);
    }
    if (x(cell[2](1), m, n) > cell[3](0) && x(cell[2](1), m, n) < cell[2](0)) { // bottom face
        vertex(0) = x(cell[2](1), m, n);
        vertex(1) = cell[2](1);
        intersectionPoints.push_back(vertex);
    }
    if (y(cell[3](0), m, n) > cell[2](1) && y(cell[3](0), m, n) < cell[0](1)) { // left face
        vertex(0) = cell[3](0);
        vertex(1) = y(cell[3](0), m, n);
        intersectionPoints.push_back(vertex);
    }

    if (intersectionPoints.size() == 0){
        isIntersecting = false;
        return;
    } else {
        isIntersecting = true;
    }

    if (intersectionPoints[1](0) == cell[1](0)){
        polygon1.push_back(cellEdit[2]);
        cellEdit.erase(cellEdit.begin() + 2);
        if (intersectionPoints[0](0) == cell[0](0)){
            // break;
        } else {
            polygon1.push_back(cellEdit[2]);
            cellEdit.erase(cellEdit.begin() + 2);
            if (intersectionPoints[0](0) == cell[1](0)){
                // break;
            } else {
                polygon1.push_back(cellEdit[0]);
                cellEdit.erase(cellEdit.begin());
            }
        }
    } else if (intersectionPoints[1](1) == cell[3](1)){
        polygon1.push_back(cellEdit[3]);
        cellEdit.erase(cellEdit.begin() + 3);
        if (intersectionPoints[0](0) == cell[0](0)){
            // break;
        } else {
            polygon1.push_back(cellEdit[0]);
            cellEdit.erase(cellEdit.begin());
            if (intersectionPoints[0](1) == cell[0](1)){
                // break;
            } else {
                polygon1.push_back(cellEdit[0]);
                cellEdit.erase(cellEdit.begin());
            }
        }
    } else if (intersectionPoints[1](0) == cell[0](0)){
        polygon1.push_back(cellEdit[0]);
        cellEdit.erase(cellEdit.begin());
        if (intersectionPoints[0](1) == cell[0](1)){
            // break;
        } else {
            polygon1.push_back(cellEdit[0]);
            cellEdit.erase(cellEdit.begin());
            if (intersectionPoints[0](0) == cell[1](0)){
                // break;
            } else {
                polygon1.push_back(cellEdit[0]);
                cellEdit.erase(cellEdit.begin());
            }
        }
    } else if (intersectionPoints[1](1) == cell[0](1)){
        polygon1.push_back(cellEdit[1]);
        cellEdit.erase(cellEdit.begin() + 1);
        if (intersectionPoints[0](0) == cell[1](0)){
            // break;
        } else {
            polygon1.push_back(cellEdit[1]);
            cellEdit.erase(cellEdit.begin() + 1);
            if (intersectionPoints[0](1) == cell[3](1)){
                // break;
            } else {
                polygon1.push_back(cellEdit[1]);
                cellEdit.erase(cellEdit.begin() + 1);
            }
        }
    }
    polygon1.push_back(intersectionPoints[0]);
    polygon1.push_back(intersectionPoints[1]);
        
    for (int i = 0; i < cellEdit.size(); i++){
        polygon2.push_back(cellEdit[i]);
    }
    polygon2.push_back(intersectionPoints[1]);
    polygon2.push_back(intersectionPoints[0]);

    if (intersectionPoints[1] == p1 || intersectionPoints[1] == p2 || intersectionPoints[1] == p3 || intersectionPoints[1] == p4){
        polygon1.erase(polygon1.end() - 1);
        polygon2.erase(polygon2.end() - 2);
    }
    if (intersectionPoints[0] == p1 || intersectionPoints[0] == p2 || intersectionPoints[0] == p3 || intersectionPoints[0] == p4){
        polygon1.erase(polygon1.end() - 2);
        polygon2.erase(polygon2.end() - 1);
    }

    Eigen::Vector2d edgeVec = intersectionPoints[1] - intersectionPoints[0];
    double crossProduct = edgeVec.x()*normalVec.y() - edgeVec.y()*normalVec.x();
    if (crossProduct < 0){
        polygon = polygon2;
    } else {
        polygon = polygon1;
    }
}

double shoeLaceFormula(std::vector<Eigen::Vector2d> polygon){
    double area = 0;
    for (int i = 0; i < polygon.size() - 1; i++){
        area += polygon[i](0)*polygon[i+1](1) - polygon[i+1](0)*polygon[i](1);
    }
    area += polygon[polygon.size() - 1](0)*polygon[0](1) - polygon[0](0)*polygon[polygon.size() - 1](1);
    area = std::abs(area)/2;
    return area;
}

bool isSectionInsideAlpha(Eigen::Vector2d p, Eigen::Vector2d interfaceMidPoint, Eigen::Vector2d normalVec, bool isX){
    double s = 0;
    if (isX == true){
        s = (p(0) - interfaceMidPoint(0))/normalVec(0);
    } else {
        s = (p(1) - interfaceMidPoint(1))/normalVec(1);
    }
    if (s > 0){
        return true;
    } else {
        return false;
    }
}

void advectalphaface(Data2D data, int cellId){
    Cell2D curCell = data.cells[cellId];
    double dx = curCell.faces[WEST]->dy;
    double dy = curCell.faces[SOUTH]->dx;
    std::vector<double> fluxLines = calculateFluxLines(data, cellId);

    // Calculate the intersection points of the flux lines with the cell edges
    // x-Line
    double xLine = 0;
    if (fluxLines[0] < 0){
        // x-Line
        xLine = -(dx/2) + (std::abs(fluxLines[0]));
    } else if (fluxLines[0] > 0){
        // x-Line
        xLine = (dx/2) - (std::abs(fluxLines[0]));
    } else {
        xLine = 0;
    }

    // y-Line
    double yLine = 0;
    if(fluxLines[1] < 0){
        // y-Line
        yLine = -(dy/2) + (std::abs(fluxLines[1]));
    } else if (fluxLines[1] > 0){
        // y-Line
        yLine = (dy/2) - (std::abs(fluxLines[1]));
    } else {
        yLine = 0;
    }


    // std::vector<std::vector<Eigen::Vector2d>> polygons;
    // Determine the alpha sections
    Vector normalVec = {curCell.normalVector[0], curCell.normalVector[1]};
    bool isIntersecting = false;
    double m = curCell.interfaceLine.m;
    double n = curCell.interfaceLine.n;               
    if (xLine == 0){
        // y-line cuts the cell
        if (fluxLines[1] < 0){
            // y-flux is negative => calculate bottom alpha value
            Eigen::Vector2d p1 = {-dx/2, yLine};
            Eigen::Vector2d p2 = {dx/2, yLine};
            Eigen::Vector2d p3 = {dx/2, -dy/2};
            Eigen::Vector2d p4 = {-dx/2, -dy/2};
            std::vector<Eigen::Vector2d> polygon;
            calculatePolygonArea(n, m, normalVec, p1, p2, p3, p4, polygon, isIntersecting);
            if (isIntersecting){
                double alphaArea = shoeLaceFormula(polygon);
                double totalArea = dx*dy;
                curCell.neighCells[SOUTH]->alphaFlux += alphaArea/totalArea;
                curCell.alphaFlux -= alphaArea/totalArea;
            } else {
                bool isInside = isSectionInsideAlpha(p1, curCell.interfaceMidPoint, normalVec, false);
                if (isInside == true){
                    double alphaArea = shoeLaceFormula({p1, p2, p3, p4});
                    double totalArea = dx*dy;
                    curCell.neighCells[SOUTH]->alphaFlux += alphaArea/totalArea;
                    curCell.alphaFlux -= alphaArea/totalArea;
                } else {
                    return;
                }
            }
        } else if (fluxLines[1] > 0){
            // y-flux is positive => calculate top alpha value
            Eigen::Vector2d p1 = {-dx/2, dy/2};
            Eigen::Vector2d p2 = {dx/2, dy/2};
            Eigen::Vector2d p3 = {dx/2, yLine};
            Eigen::Vector2d p4 = {-dx/2, yLine};
            std::vector<Eigen::Vector2d> polygon;
            calculatePolygonArea(n, m, normalVec, p1, p2, p3, p4, polygon, isIntersecting);
            if (isIntersecting){
                double alphaArea = shoeLaceFormula(polygon);
                double totalArea = dx*dy;
                curCell.neighCells[NORTH]->alphaFlux += alphaArea/totalArea;
                curCell.alphaFlux -= alphaArea/totalArea;
            } else {
                bool isInside = isSectionInsideAlpha(p1, curCell.interfaceMidPoint, normalVec, false);
                if (isInside == true){
                    double alphaArea = shoeLaceFormula({p1, p2, p3, p4});
                    double totalArea = dx*dy;
                    curCell.neighCells[NORTH]->alphaFlux += alphaArea/totalArea;
                    curCell.alphaFlux -= alphaArea/totalArea;
                } else {
                    return;
                }
            }
        }       
    } else if (yLine == 0){
        // x-line cuts the cell
        if (fluxLines[0] < 0){
            // x-flux is negative => calculate left alpha value
            Eigen::Vector2d p1 = {-dx/2, dy/2};
            Eigen::Vector2d p2 = {xLine, dy/2};
            Eigen::Vector2d p3 = {xLine, -dy/2};
            Eigen::Vector2d p4 = {-dx/2, -dy/2};
            std::vector<Eigen::Vector2d> polygon;
            calculatePolygonArea(n, m, normalVec, p1, p2, p3, p4, polygon, isIntersecting);
            if (isIntersecting){
                double alphaArea = shoeLaceFormula(polygon);
                double totalArea = dx*dy;
                curCell.neighCells[WEST]->alphaFlux += alphaArea/totalArea;
                curCell.alphaFlux -= alphaArea/totalArea;
            } else {
                bool isInside = isSectionInsideAlpha(p1, curCell.interfaceMidPoint, normalVec, true);
                if (isInside == true){
                    double alphaArea = shoeLaceFormula({p1, p2, p3, p4});
                    double totalArea = dx*dy;
                    curCell.neighCells[WEST]->alphaFlux += alphaArea/totalArea;
                    curCell.alphaFlux -= alphaArea/totalArea;
                } else {
                    return;
                }
            }
        } else if (fluxLines[0] > 0){
            // x-flux is positive => calculate right alpha value
            Eigen::Vector2d p1 = {xLine, dy/2};
            Eigen::Vector2d p2 = {dx/2, dy/2};
            Eigen::Vector2d p3 = {dx/2, -dy/2};
            Eigen::Vector2d p4 = {xLine, -dy/2};
            std::vector<Eigen::Vector2d> polygon;
            calculatePolygonArea(n, m, normalVec, p1, p2, p3, p4, polygon, isIntersecting);
            if (isIntersecting){
                double alphaArea = shoeLaceFormula(polygon);
                double totalArea = dx*dy;
                curCell.neighCells[EAST]->alphaFlux += alphaArea/totalArea;
                curCell.alphaFlux -= alphaArea/totalArea;
            } else {
                bool isInside = isSectionInsideAlpha(p1, curCell.interfaceMidPoint, normalVec, true);
                if (isInside == true){
                    double alphaArea = shoeLaceFormula({p1, p2, p3, p4});
                    double totalArea = dx*dy;
                    curCell.neighCells[EAST]->alphaFlux += alphaArea/totalArea;
                    curCell.alphaFlux -= alphaArea/totalArea;
                } else {
                    return;
                }
            }
        }
    } else if (xLine != 0 && yLine != 0){
        // x-line and y-line intersect
        if (fluxLines[0] < 0 && fluxLines[1] < 0){
            // x-flux and y-flux are negative
            // First quadrant under y-line for south direction
            Eigen::Vector2d p1 = {xLine, yLine};
            Eigen::Vector2d p2 = {dx/2, yLine};
            Eigen::Vector2d p3 = {dx/2, -dy/2};
            Eigen::Vector2d p4 = {xLine, -dy/2};
            std::vector<Eigen::Vector2d> polygonSouth;
            calculatePolygonArea(n, m, normalVec, p1, p2, p3, p4, polygonSouth, isIntersecting);
            if (isIntersecting){
                double alphaArea = shoeLaceFormula(polygonSouth);
                double totalArea = dx*dy;
                curCell.neighCells[SOUTH]->alphaFlux += alphaArea/totalArea;
                curCell.alphaFlux -= alphaArea/totalArea;
            } else {
                bool isInside = isSectionInsideAlpha(p1, curCell.interfaceMidPoint, normalVec, false);
                if (isInside == true){
                    double alphaArea = shoeLaceFormula({p1, p2, p3, p4});
                    double totalArea = dx*dy;
                    curCell.neighCells[SOUTH]->alphaFlux += alphaArea/totalArea;
                    curCell.alphaFlux -= alphaArea/totalArea;
                } else {
                    return;
                }
            }
            // Second quadrant under x-line for west direction
            p1 = {-dx/2, dy/2};
            p2 = {xLine, dy/2};
            p3 = {xLine, yLine};
            p4 = {-dx/2, yLine};
            std::vector<Eigen::Vector2d> polygonWest;
            calculatePolygonArea(n, m, normalVec, p1, p2, p3, p4, polygonWest, isIntersecting);
            if (isIntersecting){
                double alphaArea = shoeLaceFormula(polygonWest);
                double totalArea = dx*dy;
                curCell.neighCells[WEST]->alphaFlux += alphaArea/totalArea;
                curCell.alphaFlux -= alphaArea/totalArea;
            } else {
                bool isInside = isSectionInsideAlpha(p1, curCell.interfaceMidPoint, normalVec, true);
                if (isInside == true){
                    double alphaArea = shoeLaceFormula({p1, p2, p3, p4});
                    double totalArea = dx*dy;
                    curCell.neighCells[WEST]->alphaFlux += alphaArea/totalArea;
                    curCell.alphaFlux -= alphaArea/totalArea;
                } else {
                    return;
                }
            }
            // Third quadrant under y-line for south-west direction
            p1 = {-dx/2, yLine};
            p2 = {xLine, yLine};
            p3 = {xLine, -dy/2};
            p4 = {-dx/2, -dy/2};
            std::vector<Eigen::Vector2d> polygonSouthWest;
            calculatePolygonArea(n, m, normalVec, p1, p2, p3, p4, polygonSouthWest, isIntersecting);
            if (isIntersecting){
                double alphaArea = shoeLaceFormula(polygonSouthWest);
                double totalArea = dx*dy;
                curCell.neighCells[SOUTH]->neighCells[WEST]->alphaFlux += alphaArea/totalArea;
                curCell.alphaFlux -= alphaArea/totalArea;
            } else {
                bool isInside = isSectionInsideAlpha(p1, curCell.interfaceMidPoint, normalVec, false);
                if (isInside == true){
                    double alphaArea = shoeLaceFormula({p1, p2, p3, p4});
                    double totalArea = dx*dy;
                    curCell.neighCells[SOUTH]->neighCells[WEST]->alphaFlux += alphaArea/totalArea;
                    curCell.alphaFlux -= alphaArea/totalArea;
                } else {
                    return;
                }
            }
        } else if (fluxLines[0] > 0 && fluxLines[1] < 0){
            // x-flux is positive and y-flux is negative
            // First quadrant under y-line for south direction
            Eigen::Vector2d p1 = {-dx/2, yLine};
            Eigen::Vector2d p2 = {xLine, yLine};
            Eigen::Vector2d p3 = {xLine, -dy/2};
            Eigen::Vector2d p4 = {-dx/2, -dy/2};
            std::vector<Eigen::Vector2d> polygonSouth;
            calculatePolygonArea(n, m, normalVec, p1, p2, p3, p4, polygonSouth, isIntersecting);
            if (isIntersecting){
                double alphaArea = shoeLaceFormula(polygonSouth);
                double totalArea = dx*dy;
                curCell.neighCells[SOUTH]->alphaFlux += alphaArea/totalArea;
                curCell.alphaFlux -= alphaArea/totalArea;
            } else {
                bool isInside = isSectionInsideAlpha(p1, curCell.interfaceMidPoint, normalVec, false);
                if (isInside == true){
                    double alphaArea = shoeLaceFormula({p1, p2, p3, p4});
                    double totalArea = dx*dy;
                    curCell.neighCells[SOUTH]->alphaFlux += alphaArea/totalArea;
                    curCell.alphaFlux -= alphaArea/totalArea;
                } else {
                    return;
                }
            }
            // Second quadrant under x-line for east direction
            p1 = {xLine, dy/2};
            p2 = {dx/2, dy/2};
            p3 = {dx/2, yLine};
            p4 = {xLine, yLine};
            std::vector<Eigen::Vector2d> polygonEast;
            calculatePolygonArea(n, m, normalVec, p1, p2, p3, p4, polygonEast, isIntersecting);
            if (isIntersecting){
                double alphaArea = shoeLaceFormula(polygonEast);
                double totalArea = dx*dy;
                curCell.neighCells[EAST]->alphaFlux += alphaArea/totalArea;
                curCell.alphaFlux -= alphaArea/totalArea;
            } else {
                bool isInside = isSectionInsideAlpha(p1, curCell.interfaceMidPoint, normalVec, true);
                if (isInside == true){
                    double alphaArea = shoeLaceFormula({p1, p2, p3, p4});
                    double totalArea = dx*dy;
                    curCell.neighCells[EAST]->alphaFlux += alphaArea/totalArea;
                    curCell.alphaFlux -= alphaArea/totalArea;
                } else {
                    return;
                }
            }
            // Third quadrant under y-line for south-east direction
            p1 = {xLine, yLine};
            p2 = {dx/2, yLine};
            p3 = {dx/2, -dy/2};
            p4 = {xLine, -dy/2};
            std::vector<Eigen::Vector2d> polygonSouthEast;
            calculatePolygonArea(n, m, normalVec, p1, p2, p3, p4, polygonSouthEast, isIntersecting);
            if (isIntersecting){
                double alphaArea = shoeLaceFormula(polygonSouthEast);
                double totalArea = dx*dy;
                curCell.neighCells[SOUTH]->neighCells[EAST]->alphaFlux += alphaArea/totalArea;
                curCell.alphaFlux -= alphaArea/totalArea;
            } else {
                bool isInside = isSectionInsideAlpha(p1, curCell.interfaceMidPoint, normalVec, false);
                if (isInside == true){
                    double alphaArea = shoeLaceFormula({p1, p2, p3, p4});
                    double totalArea = dx*dy;
                    curCell.neighCells[SOUTH]->neighCells[EAST]->alphaFlux += alphaArea/totalArea;
                    curCell.alphaFlux -= alphaArea/totalArea;
                } else {
                    return;
                }
            }
        } else if (fluxLines[0] < 0 && fluxLines[1] > 0){
            // x-flux is negative and y-flux is positive
            // First quadrant under x-line for north direction
            Eigen::Vector2d p1 = {xLine, dy/2};
            Eigen::Vector2d p2 = {dx/2, dy/2};
            Eigen::Vector2d p3 = {dx/2, yLine};
            Eigen::Vector2d p4 = {xLine, yLine};
            std::vector<Eigen::Vector2d> polygonNorth;
            calculatePolygonArea(n, m, normalVec, p1, p2, p3, p4, polygonNorth, isIntersecting);
            if (isIntersecting){
                double alphaArea = shoeLaceFormula(polygonNorth);
                double totalArea = dx*dy;
                curCell.neighCells[NORTH]->alphaFlux += alphaArea/totalArea;
                curCell.alphaFlux -= alphaArea/totalArea;
            } else {
                bool isInside = isSectionInsideAlpha(p1, curCell.interfaceMidPoint, normalVec, false);
                if (isInside == true){
                    double alphaArea = shoeLaceFormula({p1, p2, p3, p4});
                    double totalArea = dx*dy;
                    curCell.neighCells[NORTH]->alphaFlux += alphaArea/totalArea;
                    curCell.alphaFlux -= alphaArea/totalArea;
                } else {
                    return;
                }
            } 
            // Second quadrant under y-line for west direction
            p1 = {-dx/2, yLine};
            p2 = {xLine, yLine};
            p3 = {xLine, -dy/2};
            p4 = {-dx/2, -dy/2};
            std::vector<Eigen::Vector2d> polygonWest;
            calculatePolygonArea(n, m, normalVec, p1, p2, p3, p4, polygonWest, isIntersecting);
            if (isIntersecting){
                double alphaArea = shoeLaceFormula(polygonWest);
                double totalArea = dx*dy;
                curCell.neighCells[WEST]->alphaFlux += alphaArea/totalArea;
                curCell.alphaFlux -= alphaArea/totalArea;
            } else {
                bool isInside = isSectionInsideAlpha(p1, curCell.interfaceMidPoint, normalVec, true);
                if (isInside == true){
                    double alphaArea = shoeLaceFormula({p1, p2, p3, p4});
                    double totalArea = dx*dy;
                    curCell.neighCells[WEST]->alphaFlux += alphaArea/totalArea;
                    curCell.alphaFlux -= alphaArea/totalArea;
                } else {
                    return;
                }
            } 
            // Third quadrant under x-line for north-west direction
            p1 = {-dx/2, dy/2};
            p2 = {xLine, dy/2};
            p3 = {xLine, yLine};
            p4 = {-dx/2, yLine};
            std::vector<Eigen::Vector2d> polygonNorthWest;
            calculatePolygonArea(n, m, normalVec, p1, p2, p3, p4, polygonNorthWest, isIntersecting);
            if (isIntersecting){
                double alphaArea = shoeLaceFormula(polygonNorthWest);
                double totalArea = dx*dy;
                curCell.neighCells[NORTH]->neighCells[WEST]->alphaFlux += alphaArea/totalArea;
                curCell.alphaFlux -= alphaArea/totalArea;
            } else {
                bool isInside = isSectionInsideAlpha(p1, curCell.interfaceMidPoint, normalVec, false);
                if (isInside == true){
                    double alphaArea = shoeLaceFormula({p1, p2, p3, p4});
                    double totalArea = dx*dy;
                    curCell.neighCells[NORTH]->neighCells[WEST]->alphaFlux += alphaArea/totalArea;
                    curCell.alphaFlux -= alphaArea/totalArea;
                } else {
                    return;
                }
            }
        } else if (fluxLines[0] > 0 && fluxLines[1] > 0){
            // x-flux and y-flux are positive
            // First quadrant under x-line for north direction
            Eigen::Vector2d p1 = {-dx/2, dy/2};
            Eigen::Vector2d p2 = {xLine, dy/2};
            Eigen::Vector2d p3 = {xLine, yLine};
            Eigen::Vector2d p4 = {-dx/2, yLine};
            std::vector<Eigen::Vector2d> polygonNorth;
            calculatePolygonArea(n, m, normalVec, p1, p2, p3, p4, polygonNorth, isIntersecting);
            if (isIntersecting){
                double alphaArea = shoeLaceFormula(polygonNorth);
                double totalArea = dx*dy;
                curCell.neighCells[NORTH]->alphaFlux += alphaArea/totalArea;
                curCell.alphaFlux -= alphaArea/totalArea;
            } else {
                bool isInside = isSectionInsideAlpha(p1, curCell.interfaceMidPoint, normalVec, false);
                if (isInside == true){
                    double alphaArea = shoeLaceFormula({p1, p2, p3, p4});
                    double totalArea = dx*dy;
                    curCell.neighCells[NORTH]->alphaFlux += alphaArea/totalArea;
                    curCell.alphaFlux -= alphaArea/totalArea;
                } else {
                    return;
                }
            }
            // Second quadrant under y-line for east direction
            p1 = {xLine, yLine};
            p2 = {dx/2, yLine};
            p3 = {dx/2, -dy/2};
            p4 = {xLine, -dy/2};
            std::vector<Eigen::Vector2d> polygonEast;
            calculatePolygonArea(n, m, normalVec, p1, p2, p3, p4, polygonEast, isIntersecting);
            if (isIntersecting){
                double alphaArea = shoeLaceFormula(polygonEast);
                double totalArea = dx*dy;
                curCell.neighCells[EAST]->alphaFlux += alphaArea/totalArea;
                curCell.alphaFlux -= alphaArea/totalArea;
            } else {
                bool isInside = isSectionInsideAlpha(p1, curCell.interfaceMidPoint, normalVec, true);
                if (isInside == true){
                    double alphaArea = shoeLaceFormula({p1, p2, p3, p4});
                    double totalArea = dx*dy;
                    curCell.neighCells[EAST]->alphaFlux += alphaArea/totalArea;
                    curCell.alphaFlux -= alphaArea/totalArea;
                } else {
                    return;
                }
            }
            // Third quadrant under y-line for north-east direction
            p1 = {xLine, dy/2};
            p2 = {dx/2, dy/2};
            p3 = {dx/2, yLine};
            p4 = {xLine, yLine};
            std::vector<Eigen::Vector2d> polygonNorthEast;
            calculatePolygonArea(n, m, normalVec, p1, p2, p3, p4, polygonNorthEast, isIntersecting);
            if (isIntersecting){
                double alphaArea = shoeLaceFormula(polygonNorthEast);
                double totalArea = dx*dy;
                curCell.neighCells[NORTH]->neighCells[EAST]->alphaFlux += alphaArea/totalArea;
                curCell.alphaFlux -= alphaArea/totalArea;
            } else {
                bool isInside = isSectionInsideAlpha(p1, curCell.interfaceMidPoint, normalVec, false);
                if (isInside == true){
                    double alphaArea = shoeLaceFormula({p1, p2, p3, p4});
                    double totalArea = dx*dy;
                    curCell.neighCells[NORTH]->neighCells[EAST]->alphaFlux += alphaArea/totalArea;
                    curCell.alphaFlux -= alphaArea/totalArea;
                } else {
                    return;
                }
            }
        }      
    } 
}
double interfaceFlux(Data2D& data, int cellId){
    Cell2D curCell = data.cells[cellId];
    double dx = curCell.faces[WEST]->dy;
    double dy = curCell.faces[SOUTH]->dx;
    std::vector<double> fluxLines = calculateFluxLines(data, cellId);
    return 0.0;
}

void preformXSweep(Data2D& data){
    std::vector<double> interAlpha;
    for (int i = 0; i < data.nCells; i++){
        Cell2D curCell = data.cells[i];
        double fe = 0;
        double fw = 0;
        if (curCell.bType_sc == DIRICHLET || curCell.bType_sc == NEUMANN){
            interAlpha.push_back(curCell.alpha);
            continue;
        } else {

            if (curCell.neighCells[EAST]->alpha != 0 && curCell.neighCells[EAST]->alpha != 1){
                fe = interfaceFlux(data, curCell.neighCells[EAST]->id);
            } else {
                fe = curCell.neighCells[EAST]->alpha*data.dt*curCell.faces[EAST]->u[CORRECTED_2]*curCell.faces[EAST]->dy;
            }

            if (curCell.neighCells[WEST]->alpha != 0 && curCell.neighCells[WEST]->alpha != 1){
                fw = interfaceFlux(data, curCell.neighCells[WEST]->id);
            } else {
                fw = curCell.neighCells[EAST]->alpha*data.dt*curCell.faces[WEST]->u[CORRECTED_2]*curCell.faces[WEST]->dy;
            }

            double tmp = (curCell.alpha + fw - fe)/(1 - data.dt/curCell.faces[NORTH]->dx*(curCell.faces[WEST]->u[CORRECTED_2] - curCell.faces[EAST]->u[CORRECTED_2]));
            interAlpha.push_back(tmp);
        }
    }
    for (int i = 0; i < data.nCells; i++){
        data.cells[i].alpha = interAlpha[i];
    }
    std::stack<Cell2D> cornerCells;
    for (int i = 0; i < data.nCells; i++){
        Cell2D curCell = data.cells[i];
        if (curCell.bType_sc == NEUMANN) {
            if ((curCell.neighCells[EAST] == nullptr || curCell.neighCells[EAST]->bType_sc == SOLID) && (curCell.neighCells[NORTH] == nullptr || curCell.neighCells[NORTH]->bType_sc == SOLID)) {
                cornerCells.push(curCell);
            } else if ((curCell.neighCells[WEST] == nullptr || curCell.neighCells[WEST]->bType_sc == SOLID) && (curCell.neighCells[NORTH] == nullptr || curCell.neighCells[NORTH]->bType_sc == SOLID)) {
                cornerCells.push(curCell);
            } else if ((curCell.neighCells[EAST] == nullptr || curCell.neighCells[EAST]->bType_sc == SOLID) && (curCell.neighCells[SOUTH] == nullptr || curCell.neighCells[SOUTH]->bType_sc == SOLID)) {
                cornerCells.push(curCell);
            } else if ((curCell.neighCells[WEST] == nullptr || curCell.neighCells[WEST]->bType_sc == SOLID) && (curCell.neighCells[SOUTH] == nullptr || curCell.neighCells[SOUTH]->bType_sc == SOLID)) {
                cornerCells.push(curCell);
            } else if (curCell.neighCells[WEST]->bType_sc == SOLID || curCell.neighCells[WEST] == nullptr){
                curCell.alpha = curCell.neighCells[EAST]->alpha;
            } else if (curCell.neighCells[EAST]->bType_sc == SOLID || curCell.neighCells[EAST] == nullptr){
                curCell.alpha = curCell.neighCells[WEST]->alpha;
            } else if (curCell.neighCells[NORTH]->bType_sc == SOLID || curCell.neighCells[NORTH] == nullptr){
                curCell.alpha = curCell.neighCells[SOUTH]->alpha;
            } else if (curCell.neighCells[SOUTH]->bType_sc == SOLID || curCell.neighCells[SOUTH] == nullptr){
                curCell.alpha = curCell.neighCells[NORTH]->alpha;
            }
        }
    }
    while (!cornerCells.empty()){
        Cell2D curCell = cornerCells.top();
        cornerCells.pop();
        if ((curCell.neighCells[EAST] == nullptr || curCell.neighCells[EAST]->bType_sc == SOLID) && (curCell.neighCells[NORTH] == nullptr || curCell.neighCells[NORTH]->bType_sc == SOLID)) {
            curCell.alpha = 0.5 * (curCell.neighCells[WEST]->alpha + curCell.neighCells[SOUTH]->alpha);
        } else if ((curCell.neighCells[WEST] == nullptr || curCell.neighCells[WEST]->bType_sc == SOLID) && (curCell.neighCells[NORTH] == nullptr || curCell.neighCells[NORTH]->bType_sc == SOLID)) {
            curCell.alpha = 0.5 * (curCell.neighCells[EAST]->alpha + curCell.neighCells[SOUTH]->alpha);
        } else if ((curCell.neighCells[EAST] == nullptr || curCell.neighCells[EAST]->bType_sc == SOLID) && (curCell.neighCells[SOUTH] == nullptr || curCell.neighCells[SOUTH]->bType_sc == SOLID)) {
            curCell.alpha = 0.5 * (curCell.neighCells[WEST]->alpha + curCell.neighCells[NORTH]->alpha);
        } else if ((curCell.neighCells[WEST] == nullptr || curCell.neighCells[WEST]->bType_sc == SOLID) && (curCell.neighCells[SOUTH] == nullptr || curCell.neighCells[SOUTH]->bType_sc == SOLID)) {
            curCell.alpha = 0.5 * (curCell.neighCells[EAST]->alpha + curCell.neighCells[NORTH]->alpha);
        }
    }
}

void preformYSweep(Data2D& data){
    std::vector<double> interAlpha;
    for (int i = 0; i < data.nCells; i++){
        Cell2D curCell = data.cells[i];
        double gn = 0;
        double gs = 0;
        if (curCell.bType_sc == DIRICHLET || curCell.bType_sc == NEUMANN){
            interAlpha.push_back(curCell.alpha);
            continue;
        } else {
                
                if (curCell.neighCells[NORTH]->alpha != 0 && curCell.neighCells[NORTH]->alpha != 1){
                    gn = interfaceFlux(data, curCell.neighCells[NORTH]->id);
                } else {
                    gn = curCell.neighCells[NORTH]->alpha*data.dt*curCell.faces[NORTH]->v[CORRECTED_2]*curCell.faces[NORTH]->dx;
                }
    
                if (curCell.neighCells[SOUTH]->alpha != 0 && curCell.neighCells[SOUTH]->alpha != 1){
                    gs = interfaceFlux(data, curCell.neighCells[SOUTH]->id);
                } else {
                    gs = curCell.neighCells[SOUTH]->alpha*data.dt*curCell.faces[SOUTH]->v[CORRECTED_2]*curCell.faces[SOUTH]->dx;
                }
    
                double tmp = (curCell.alpha + gs - gn + curCell.alpha*data.dt/curCell.faces[WEST]->dy*(curCell.faces[SOUTH]->v[CORRECTED_2] - curCell.faces[NORTH]->v[CORRECTED_2]));
                interAlpha.push_back(tmp);
        }
    }
    for (int i = 0; i < data.nCells; i++){
        data.cells[i].alpha = interAlpha[i];
    }
    std::stack<Cell2D> cornerCells;
    for (int i = 0; i < data.nCells; i++){
        Cell2D curCell = data.cells[i];
        if (curCell.bType_sc == NEUMANN) {
            if ((curCell.neighCells[EAST] == nullptr || curCell.neighCells[EAST]->bType_sc == SOLID) && (curCell.neighCells[NORTH] == nullptr || curCell.neighCells[NORTH]->bType_sc == SOLID)) {
                cornerCells.push(curCell);
            } else if ((curCell.neighCells[WEST] == nullptr || curCell.neighCells[WEST]->bType_sc == SOLID) && (curCell.neighCells[NORTH] == nullptr || curCell.neighCells[NORTH]->bType_sc == SOLID)) {
                cornerCells.push(curCell);
            } else if ((curCell.neighCells[EAST] == nullptr || curCell.neighCells[EAST]->bType_sc == SOLID) && (curCell.neighCells[SOUTH] == nullptr || curCell.neighCells[SOUTH]->bType_sc == SOLID)) {
                cornerCells.push(curCell);
            } else if ((curCell.neighCells[WEST] == nullptr || curCell.neighCells[WEST]->bType_sc == SOLID) && (curCell.neighCells[SOUTH] == nullptr || curCell.neighCells[SOUTH]->bType_sc == SOLID)) {
                cornerCells.push(curCell);
            } else if (curCell.neighCells[WEST]->bType_sc == SOLID || curCell.neighCells[WEST] == nullptr){
                curCell.alpha = curCell.neighCells[EAST]->alpha;
            } else if (curCell.neighCells[EAST]->bType_sc == SOLID || curCell.neighCells[EAST] == nullptr){
                curCell.alpha = curCell.neighCells[WEST]->alpha;
            } else if (curCell.neighCells[NORTH]->bType_sc == SOLID || curCell.neighCells[NORTH] == nullptr){
                curCell.alpha = curCell.neighCells[SOUTH]->alpha;
            } else if (curCell.neighCells[SOUTH]->bType_sc == SOLID || curCell.neighCells[SOUTH] == nullptr){
                curCell.alpha = curCell.neighCells[NORTH]->alpha;
            }
        }
    }
    while (!cornerCells.empty()){
        Cell2D curCell = cornerCells.top();
        cornerCells.pop();
        if ((curCell.neighCells[EAST] == nullptr || curCell.neighCells[EAST]->bType_sc == SOLID) && (curCell.neighCells[NORTH] == nullptr || curCell.neighCells[NORTH]->bType_sc == SOLID)) {
            curCell.alpha = 0.5 * (curCell.neighCells[WEST]->alpha + curCell.neighCells[SOUTH]->alpha);
        } else if ((curCell.neighCells[WEST] == nullptr || curCell.neighCells[WEST]->bType_sc == SOLID) && (curCell.neighCells[NORTH] == nullptr || curCell.neighCells[NORTH]->bType_sc == SOLID)) {
            curCell.alpha = 0.5 * (curCell.neighCells[EAST]->alpha + curCell.neighCells[SOUTH]->alpha);
        } else if ((curCell.neighCells[EAST] == nullptr || curCell.neighCells[EAST]->bType_sc == SOLID) && (curCell.neighCells[SOUTH] == nullptr || curCell.neighCells[SOUTH]->bType_sc == SOLID)) {
            curCell.alpha = 0.5 * (curCell.neighCells[WEST]->alpha + curCell.neighCells[NORTH]->alpha);
        } else if ((curCell.neighCells[WEST] == nullptr || curCell.neighCells[WEST]->bType_sc == SOLID) && (curCell.neighCells[SOUTH] == nullptr || curCell.neighCells[SOUTH]->bType_sc == SOLID)) {
            curCell.alpha = 0.5 * (curCell.neighCells[EAST]->alpha + curCell.neighCells[NORTH]->alpha);
        }
    }
    
}

void handleExcessAlpha(Data2D& data){
    for (int i = 0; i < data.nCells; i++){
        Cell2D curCell = data.cells[i];
        if (curCell.alpha < 0){
            curCell.alpha = 0;
        } else if (curCell.alpha > 1){
            curCell.alpha = 1;
        }
    }
}          

    

void advectAlpha(Data2D& data){
    reconstructInterfaceLines(data);
    preformXSweep(data);
    handleExcessAlpha(data);
    reconstructInterfaceLines(data);
    preformYSweep(data);
    handleExcessAlpha(data);
}