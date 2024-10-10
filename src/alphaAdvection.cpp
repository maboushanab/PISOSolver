
#include "data.h"
#include "interface.h"

double calculateFluxLength(Data2D data, int cellId, std::string direction){
    Cell2D curCell = data.cells[cellId];
    if (direction == "WEST"){
        if (curCell.faces[WEST]->u < 0){
            return -(data.dt * curCell.faces[WEST]->u[CORRECTED_2]); 
        } else {
            return 0;
        }
    } else if (direction == "EAST"){
        if (curCell.faces[EAST]->u > 0){
            return data.dt * curCell.faces[EAST]->u[CORRECTED_2];
        } else {
            return 0;
        }
    } else if (direction == "SOUTH"){
        if (curCell.faces[SOUTH]->v < 0){
            return -(data.dt * curCell.faces[SOUTH]->v[CORRECTED_2]);
        } else {
            return 0;
        }
    } else if (direction == "NORTH"){
        if (curCell.faces[NORTH]->v > 0){
            return data.dt * curCell.faces[NORTH]->v[CORRECTED_2];
        } else {
            return 0;
        }
    } else {
        std::cerr << "Invalid direction" << std::endl;
        return 0;
    }
}

void calculatePolygonArea(Data2D& data, int cellId, Eigen::Vector2d p1, Eigen::Vector2d p2, Eigen::Vector2d p3, Eigen::Vector2d p4, std::vector<Eigen::Vector2d>& polygon, bool& isIntersecting){
    Cell2D *curCell = &data.cells[cellId];
    std::vector<Eigen::Vector2d> cell = {
        p1, p2, p3, p4
    }; 
    Eigen::Vector2d vertex;
    std::vector<Eigen::Vector2d> intersectionPoints;
    std::vector<Eigen::Vector2d> polygon1;
    std::vector<Eigen::Vector2d> polygon2;
    std::vector<bool> intersections;
    // potential intersection at the top face (y = dy/2 and x \in [-dx/2, dx/2])
    double lambdaTop = 0;
    if (curCell->normalVector[0] != 0){ // If the interface is not completely horizontal (to avoid division by zero)
        lambdaTop = -(p1[1] - curCell->interfaceVectorLine.origin[1])/curCell->normalVector[0];
    } else {
        lambdaTop = 1e6; // throws the point outside the interval
    }
    double x_top = curCell->interfaceVectorLine.origin[0] + lambdaTop*curCell->normalVector[1];
    if (x_top >= p1[0] && x_top < p2[0]){
        intersections.push_back(true);
    } else {
        intersections.push_back(false);
    }

    //  potential intersection at the right face (x = dx/2 and y \in [-dy/2, dy/2])
    double lambdaRight = 0;
    if (curCell->normalVector[1] != 0){ // If the interface is not completely vertical (to avoid division by zero)
        lambdaRight = (p2[0] - curCell->interfaceVectorLine.origin[0])/curCell->normalVector[1];
    } else {
        lambdaRight = 1e6; // throws the point outside the interval
    }
    double y_right = curCell->interfaceVectorLine.origin[1] - lambdaRight*curCell->normalVector[0];
    if (y_right > p3[1] && y_right <= p2[1]){
        intersections.push_back(true);
    } else {
        intersections.push_back(false);
    }

    // potential intersection at the bottom face (y = -dy/2 and x \in [-dx/2, dx/2])
    double lambdaBottom = 0;
    if (curCell->normalVector[0] != 0){ // If the interface is not completely horizontal (to avoid division by zero)
        lambdaBottom = -(-p3[1] - curCell->interfaceVectorLine.origin[1])/curCell->normalVector[0];
    } else {
        lambdaBottom = 1e6; // throws the point outside the interval
    }
    double x_bottom = curCell->interfaceVectorLine.origin[0] + lambdaBottom*curCell->normalVector[1];
    if (x_bottom > p1[0] && x_bottom <= p2[0]){
        intersections.push_back(true);
    } else {
        intersections.push_back(false);
    }

    // potential intersection at the left face (x = -dx/2 and y \in [-dy/2, dy/2])
    double lambdaLeft = 0;
    if (curCell->normalVector[1] != 0){ // If the interface is not completely vertical (to avoid division by zero)
        lambdaLeft = (p4[0] - curCell->interfaceVectorLine.origin[0])/curCell->normalVector[1];
    } else {
        lambdaLeft = 1e6; // throws the point outside the interval
    }
    double y_left = curCell->interfaceVectorLine.origin[1] - lambdaLeft*curCell->normalVector[0];
    if (y_left >= p3[1] && y_left < p2[1]){
        intersections.push_back(true);
    } else {
        intersections.push_back(false);
    }
    if (!intersections[0] && !intersections[1] && !intersections[2] && !intersections[3]){
        isIntersecting = false;
        //std::cout << "No intersection" << std::endl;
        return;
    } else {
        isIntersecting = true;  
    }

    int caseId = 0;
    // Combinations of intersections
    if (intersections[0] && intersections[1]){
        caseId = 1;
        intersectionPoints.push_back(Eigen::Vector2d(x_top, p1[1]));
        intersectionPoints.push_back(Eigen::Vector2d(p2[0], y_right));
    } else if (intersections[0] && intersections[2]){
        caseId = 2;
        intersectionPoints.push_back(Eigen::Vector2d(x_top, p1[1]));
        intersectionPoints.push_back(Eigen::Vector2d(x_bottom, p3[1]));
    } else if (intersections[0] && intersections[3]){
        caseId = 3;
        intersectionPoints.push_back(Eigen::Vector2d(x_top, p1[1]));
        intersectionPoints.push_back(Eigen::Vector2d(p1[0], y_left));
    } else if (intersections[1] && intersections[2]){
        caseId = 4;
        intersectionPoints.push_back(Eigen::Vector2d(p2[0], y_right));
        intersectionPoints.push_back(Eigen::Vector2d(x_bottom, p3[1]));
    } else if (intersections[1] && intersections[3]){
        caseId = 5;
        intersectionPoints.push_back(Eigen::Vector2d(p2[0], y_right));
        intersectionPoints.push_back(Eigen::Vector2d(p1[0], y_left));
    } else if (intersections[2] && intersections[3]){
        caseId = 6;
        intersectionPoints.push_back(Eigen::Vector2d(x_bottom, p3[1]));
        intersectionPoints.push_back(Eigen::Vector2d(p1[0], y_left));
    };

    switch (caseId){
        case 1:
            // Clockwise order polygon
            polygon1.push_back(intersectionPoints[0]);
            polygon1.push_back(cell[1]);
            polygon1.push_back(intersectionPoints[1]);
            // Counter-clockwise order polygon
            polygon2.push_back(intersectionPoints[0]);
            polygon2.push_back(cell[0]);
            polygon2.push_back(cell[3]);
            polygon2.push_back(cell[2]);
            polygon2.push_back(intersectionPoints[1]);
            break;
        case 2:
            // Clockwise order polygon
            polygon1.push_back(intersectionPoints[0]);
            polygon1.push_back(cell[1]);
            polygon1.push_back(cell[2]);
            polygon1.push_back(intersectionPoints[1]);
            // Counter-clockwise order polygon
            polygon2.push_back(intersectionPoints[0]);
            polygon2.push_back(cell[0]);
            polygon2.push_back(cell[3]);
            polygon2.push_back(intersectionPoints[1]);
            break;
        case 3:
            // Clockwise order polygon
            polygon1.push_back(intersectionPoints[0]);
            polygon1.push_back(cell[1]);
            polygon1.push_back(cell[2]);
            polygon1.push_back(cell[3]);
            polygon1.push_back(intersectionPoints[1]);
            // Counter-clockwise order polygon
            polygon2.push_back(intersectionPoints[0]);
            polygon2.push_back(cell[0]);
            polygon2.push_back(intersectionPoints[1]);
            break;
        case 4:
            // Clockwise order polygon
            polygon1.push_back(intersectionPoints[0]);
            polygon1.push_back(cell[2]);
            polygon1.push_back(intersectionPoints[1]);

            // Counter-clockwise order polygon
            polygon2.push_back(intersectionPoints[0]);
            polygon2.push_back(cell[1]);
            polygon2.push_back(cell[0]);
            polygon2.push_back(cell[3]);
            polygon2.push_back(intersectionPoints[1]);
            break;
        case 5:
            // Clockwise order polygon
            polygon1.push_back(intersectionPoints[0]);
            polygon1.push_back(cell[2]);
            polygon1.push_back(cell[3]);
            polygon1.push_back(intersectionPoints[1]);

            // Counter-clockwise order polygon
            polygon2.push_back(intersectionPoints[0]);
            polygon2.push_back(cell[1]);
            polygon2.push_back(cell[0]);
            polygon2.push_back(intersectionPoints[1]);
            break;
        case 6:
            // Clockwise order polygon
            polygon1.push_back(intersectionPoints[0]);
            polygon1.push_back(cell[3]);
            polygon1.push_back(intersectionPoints[1]);

            // Counter-clockwise order polygon
            polygon2.push_back(intersectionPoints[0]);
            polygon2.push_back(cell[2]);
            polygon2.push_back(cell[1]);
            polygon2.push_back(cell[0]);
            polygon2.push_back(intersectionPoints[1]);
            break;
    }
    Eigen::Vector2d edgeVec = intersectionPoints[1] - intersectionPoints[0];  //edgevector for going clockwise through polygon 2
    Eigen::Vector2d normalVec = {curCell->normalVector[0], curCell->normalVector[1]};
    double crossProduct = edgeVec.x()*normalVec.y() - edgeVec.y()*normalVec.x();
    if (crossProduct < 0){ // if the cross product is negative, the normal vector is pointing inward into polygon 2 => polygon 1 is the correct polygon
        polygon = polygon1;
    } else {
        polygon = polygon2;
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

bool isSectionInsideAlpha(Data2D data, int cellId, double p, std::string direction){
    double s = 0;
    if (direction == "WEST" || direction == "EAST"){
        s = (p - data.cells[cellId].interfaceVectorLine.origin[0])/data.cells[cellId].normalVector[0];
    } else {
        s = (p - data.cells[cellId].interfaceVectorLine.origin[1])/data.cells[cellId].normalVector[1];
    }
    if (s > 0){
        return false;
    } else {
        return true;
    }
}


double interfaceFlux(Data2D& data, int cellId, std::string direction){
    Cell2D curCell = data.cells[cellId];
    double dx = curCell.faces[WEST]->dy;
    double dy = curCell.faces[SOUTH]->dx;

    double fluxLength = calculateFluxLength(data, cellId, direction);
    std::cout << "fluxLength: " << fluxLength << std::endl;

    if (fluxLength == 0){
        return 0.0;
    }

    double fluxPoint = 0;
    bool isIntersecting = false;
    std::vector<Eigen::Vector2d> polygon;
    Eigen::Vector2d p1 = {0, 0};
    Eigen::Vector2d p2 = {0, 0};
    Eigen::Vector2d p3 = {0, 0};
    Eigen::Vector2d p4 = {0, 0};

    if (direction == "EAST"){
        fluxPoint = dx/2 - fluxLength;
        p1 = {fluxPoint, dy/2};
        p2 = {dx/2, dy/2};
        p3 = {dx/2, -dy/2};
        p4 = {fluxPoint, -dy/2};
    } else if (direction == "WEST"){
        fluxPoint = -dx/2 + fluxLength;
        p1 = {-dx/2, dy/2};
        p2 = {fluxPoint, dy/2};
        p3 = {fluxPoint, -dy/2};
        p4 = {-dx/2, -dy/2};
    } else if (direction == "NORTH"){
        fluxPoint = dy/2 - fluxLength;
        p1 = {-dx/2, dy/2};
        p2 = {dx/2, dy/2};
        p3 = {dx/2, fluxPoint};
        p4 = {-dx/2, fluxPoint};
    } else if (direction == "SOUTH"){
        fluxPoint = -dy/2 + fluxLength;
        p1 = {-dx/2, fluxPoint};
        p2 = {dx/2, fluxPoint};
        p3 = {dx/2, -dy/2};
        p4 = {-dx/2, -dy/2};
    } else {
        std::cerr << "Invalid direction" << std::endl;
        return 0.0;
    }
    // std::cout << "fluxPoint: " << fluxPoint << std::endl;
    // std::cout << "p1: " << p1(0) << " " << p1(1) << std::endl;
    // std::cout << "p2: " << p2(0) << " " << p2(1) << std::endl;
    // std::cout << "p3: " << p3(0) << " " << p3(1) << std::endl;
    // std::cout << "p4: " << p4(0) << " " << p4(1) << std::endl;
    calculatePolygonArea(data, cellId, p1, p2, p3, p4, polygon, isIntersecting);
    // std::cout << "polygon size: " << polygon.size() << std::endl;
    // for (int i = 0; i < polygon.size(); i++){
    //     std::cout << "polygon: " << polygon[i](0) << " " << polygon[i](1) << std::endl;
    // }
    if (isIntersecting){
        double liquidArea = shoeLaceFormula(polygon);
        double totalArea = dx*dy;
        return liquidArea/totalArea;
    } else {
        bool isInside = isSectionInsideAlpha(data, cellId, fluxPoint, direction);
        std::cout << "isInside: " << isInside << std::endl;
        if (isInside == true){
            double liquidArea = shoeLaceFormula({p1, p2, p3, p4});
            double totalArea = dx*dy;
            return liquidArea/totalArea;
        } else {
            return 0.0;
        }
    }
    return 0;
}

void preformXSweep(Data2D& data){
    // std::cout << "Performing X sweep" << std::endl;
    std::vector<double> interAlpha;
    for (int i = 0; i < data.nCells; i++){
        Cell2D *curCell = &data.cells[i];
        double fe = 0;
        double fw = 0;
        if (curCell->bType_p == DIRICHLET || curCell->bType_p == NEUMANN){
            interAlpha.push_back(curCell->alpha);
            continue;
        } else {
            if (curCell->faces[EAST]->u[CORRECTED_2] > 0){
                if (curCell->alpha != 0 && curCell->alpha != 1){
                    fe = interfaceFlux(data, curCell->id, "EAST");
                    std::cout << "fe geometrically: " << fe << std::endl; 
                } else {
                    fe = curCell->alpha*data.dt*curCell->faces[EAST]->u[CORRECTED_2]*curCell->faces[EAST]->dy;
                }
            } else {
                if (curCell->neighCells[EAST] == nullptr){
                    fe = 0;
                } else {
                    if (curCell->neighCells[EAST]->alpha != 0 && curCell->neighCells[EAST]->alpha != 1){
                    fe = interfaceFlux(data, curCell->neighCells[EAST]->id, "WEST");
                    std::cout << "fe geometrically: " << fe << std::endl;
                    } else {
                        fe = curCell->neighCells[EAST]->alpha*data.dt*curCell->faces[EAST]->u[CORRECTED_2]*curCell->faces[EAST]->dy;
                    }
                }   
            }
            if (curCell->faces[WEST]->u[CORRECTED_2] < 0){
                if (curCell->alpha != 0 && curCell->alpha != 1){
                    fw = interfaceFlux(data, curCell->id, "WEST");
                    std::cout << "fw geometrically: " << fw << std::endl;
                } else {
                    fw = curCell->alpha*data.dt*curCell->faces[WEST]->u[CORRECTED_2]*curCell->faces[WEST]->dy;
                }
            } else {
                if (curCell->neighCells[WEST] == nullptr){
                    fw = 0;
                } else {
                    if (curCell->neighCells[WEST]->alpha != 0 && curCell->neighCells[WEST]->alpha != 1){
                    fw = interfaceFlux(data, curCell->neighCells[WEST]->id, "EAST");
                    std::cout << "fw geometrically: " << fw << std::endl;
                    } else {
                        fw = curCell->neighCells[WEST]->alpha*data.dt*curCell->faces[WEST]->u[CORRECTED_2]*curCell->faces[WEST]->dy;
                    }
                } 
            }

            double tmp = (curCell->alpha + fw - fe)/(1 - (data.dt*(curCell->faces[WEST]->u[CORRECTED_2] - curCell->faces[EAST]->u[CORRECTED_2]))/curCell->faces[NORTH]->dx);
            // std::cout << "cellId: " << curCell.id << " alpha: " << tmp << std::endl;
            interAlpha.push_back(tmp);
        }
    }
    for (int i = 0; i < data.nCells; i++){
        data.cells[i].alpha = interAlpha[i];
    }
    std::stack<int> cornerCells;
    for (int i = 0; i < data.nCells; i++){
        Cell2D *curCell = &data.cells[i];
        if (curCell->bType_p == NEUMANN) {
            if ((curCell->neighCells[EAST] == nullptr || curCell->neighCells[EAST]->bType_p == SOLID) && (curCell->neighCells[NORTH] == nullptr || curCell->neighCells[NORTH]->bType_p == SOLID)) {
                cornerCells.push(curCell->id);
            } else if ((curCell->neighCells[WEST] == nullptr || curCell->neighCells[WEST]->bType_p == SOLID) && (curCell->neighCells[NORTH] == nullptr || curCell->neighCells[NORTH]->bType_p == SOLID)) {
                cornerCells.push(curCell->id);
            } else if ((curCell->neighCells[EAST] == nullptr || curCell->neighCells[EAST]->bType_p == SOLID) && (curCell->neighCells[SOUTH] == nullptr || curCell->neighCells[SOUTH]->bType_p == SOLID)) {
                cornerCells.push(curCell->id);
            } else if ((curCell->neighCells[WEST] == nullptr || curCell->neighCells[WEST]->bType_p == SOLID) && (curCell->neighCells[SOUTH] == nullptr || curCell->neighCells[SOUTH]->bType_p == SOLID)) {
                cornerCells.push(curCell->id);
            } else if (curCell->neighCells[WEST] == nullptr ||curCell->neighCells[WEST]->bType_p == SOLID ){
                curCell->alpha = curCell->neighCells[EAST]->alpha;
            } else if (curCell->neighCells[EAST] == nullptr || curCell->neighCells[EAST]->bType_p == SOLID){
                curCell->alpha = curCell->neighCells[WEST]->alpha;
            } else if (curCell->neighCells[NORTH] == nullptr || curCell->neighCells[NORTH]->bType_p == SOLID){
                curCell->alpha = curCell->neighCells[SOUTH]->alpha;
            } else if (curCell->neighCells[SOUTH] == nullptr || curCell->neighCells[SOUTH]->bType_p == SOLID){
                curCell->alpha = curCell->neighCells[NORTH]->alpha;
            }
        }
    }
    while (!cornerCells.empty()){
        int id = cornerCells.top();
        Cell2D *curCell = &data.cells[id];
        cornerCells.pop();
        if ((curCell->neighCells[EAST] == nullptr || curCell->neighCells[EAST]->bType_p == SOLID) && (curCell->neighCells[NORTH] == nullptr || curCell->neighCells[NORTH]->bType_p == SOLID)) {
            curCell->alpha = 0.5 * (curCell->neighCells[WEST]->alpha + curCell->neighCells[SOUTH]->alpha);
        } else if ((curCell->neighCells[WEST] == nullptr || curCell->neighCells[WEST]->bType_p == SOLID) && (curCell->neighCells[NORTH] == nullptr || curCell->neighCells[NORTH]->bType_p == SOLID)) {
            curCell->alpha = 0.5 * (curCell->neighCells[EAST]->alpha + curCell->neighCells[SOUTH]->alpha);
        } else if ((curCell->neighCells[EAST] == nullptr || curCell->neighCells[EAST]->bType_p == SOLID) && (curCell->neighCells[SOUTH] == nullptr || curCell->neighCells[SOUTH]->bType_p == SOLID)) {
            curCell->alpha = 0.5 * (curCell->neighCells[WEST]->alpha + curCell->neighCells[NORTH]->alpha);
        } else if ((curCell->neighCells[WEST] == nullptr || curCell->neighCells[WEST]->bType_p == SOLID) && (curCell->neighCells[SOUTH] == nullptr || curCell->neighCells[SOUTH]->bType_p == SOLID)) {
            curCell->alpha = 0.5 * (curCell->neighCells[EAST]->alpha + curCell->neighCells[NORTH]->alpha);
        }
    }
}

void preformYSweep(Data2D& data){
    // std::cout << "Performing Y sweep" << std::endl;
    std::vector<double> interAlpha;
    for (int i = 0; i < data.nCells; i++){
        Cell2D *curCell = &data.cells[i];
        double gn = 0;
        double gs = 0;
        if (curCell->bType_p == DIRICHLET || curCell->bType_p == NEUMANN){
            interAlpha.push_back(curCell->alpha);
            continue;
        } else {
            if (curCell->faces[NORTH]->v[CORRECTED_2] > 0){
                if (curCell->alpha != 0 && curCell->alpha != 1){
                    gn = interfaceFlux(data, curCell->id, "NORTH");
                } else {
                    gn = curCell->alpha*data.dt*curCell->faces[NORTH]->v[CORRECTED_2]*curCell->faces[NORTH]->dx;
                }
            } else {
                if (curCell->neighCells[NORTH] == nullptr){
                    gn = 0;
                } else {
                    if (curCell->neighCells[NORTH]->alpha != 0 && curCell->neighCells[NORTH]->alpha != 1){
                    gn = interfaceFlux(data, curCell->neighCells[NORTH]->id, "SOUTH");
                    } else {
                        gn = curCell->neighCells[NORTH]->alpha*data.dt*curCell->faces[NORTH]->v[CORRECTED_2]*curCell->faces[NORTH]->dx;
                    }
                }
            }
            if (curCell->faces[SOUTH]->v[CORRECTED_2] < 0){
                if (curCell->alpha != 0 && curCell->alpha != 1){
                    gs = interfaceFlux(data, curCell->id, "SOUTH");
                } else {
                    gs = curCell->alpha*data.dt*curCell->faces[SOUTH]->v[CORRECTED_2]*curCell->faces[SOUTH]->dx;
                }
            } else {
                if (curCell->neighCells[SOUTH] == nullptr){
                    gs = 0;
                } else {
                    if (curCell->neighCells[SOUTH]->alpha != 0 && curCell->neighCells[SOUTH]->alpha != 1){
                        gs = interfaceFlux(data, curCell->neighCells[SOUTH]->id, "NORTH");
                    } else {
                        gs = curCell->neighCells[SOUTH]->alpha*data.dt*curCell->faces[SOUTH]->v[CORRECTED_2]*curCell->faces[SOUTH]->dx;
                    }
                }
            }
    
            double tmp = (curCell->alpha + gs - gn + curCell->alpha*data.dt/curCell->faces[WEST]->dy*(curCell->faces[SOUTH]->v[CORRECTED_2] - curCell->faces[NORTH]->v[CORRECTED_2]));
            // std::cout << "cellId: " << curCell.id << " alpha: " << tmp << std::endl;
            interAlpha.push_back(tmp);
        }
    }
    for (int i = 0; i < data.nCells; i++){
        data.cells[i].alpha = interAlpha[i];
    }
    std::stack<int> cornerCells;
    for (int i = 0; i < data.nCells; i++){
        Cell2D *curCell = &data.cells[i];
        if (curCell->bType_p == NEUMANN) {
            if ((curCell->neighCells[EAST] == nullptr || curCell->neighCells[EAST]->bType_p == SOLID) && (curCell->neighCells[NORTH] == nullptr || curCell->neighCells[NORTH]->bType_p == SOLID)) {
                cornerCells.push(curCell->id);
            } else if ((curCell->neighCells[WEST] == nullptr || curCell->neighCells[WEST]->bType_p == SOLID) && (curCell->neighCells[NORTH] == nullptr || curCell->neighCells[NORTH]->bType_p == SOLID)) {
                cornerCells.push(curCell->id);
            } else if ((curCell->neighCells[EAST] == nullptr || curCell->neighCells[EAST]->bType_p == SOLID) && (curCell->neighCells[SOUTH] == nullptr || curCell->neighCells[SOUTH]->bType_p == SOLID)) {
                cornerCells.push(curCell->id);
            } else if ((curCell->neighCells[WEST] == nullptr || curCell->neighCells[WEST]->bType_p == SOLID) && (curCell->neighCells[SOUTH] == nullptr || curCell->neighCells[SOUTH]->bType_p == SOLID)) {
                cornerCells.push(curCell->id);
            } else if (curCell->neighCells[WEST] == nullptr ||curCell->neighCells[WEST]->bType_p == SOLID ){
                curCell->alpha = curCell->neighCells[EAST]->alpha;
            } else if (curCell->neighCells[EAST] == nullptr || curCell->neighCells[EAST]->bType_p == SOLID){
                curCell->alpha = curCell->neighCells[WEST]->alpha;
            } else if (curCell->neighCells[NORTH] == nullptr || curCell->neighCells[NORTH]->bType_p == SOLID){
                curCell->alpha = curCell->neighCells[SOUTH]->alpha;
            } else if (curCell->neighCells[SOUTH] == nullptr || curCell->neighCells[SOUTH]->bType_p == SOLID){
                curCell->alpha = curCell->neighCells[NORTH]->alpha;
            }
        }
    }
    while (!cornerCells.empty()){
        int id = cornerCells.top();
        Cell2D *curCell = &data.cells[id];
        cornerCells.pop();
        if ((curCell->neighCells[EAST] == nullptr || curCell->neighCells[EAST]->bType_p == SOLID) && (curCell->neighCells[NORTH] == nullptr || curCell->neighCells[NORTH]->bType_p == SOLID)) {
            curCell->alpha = 0.5 * (curCell->neighCells[WEST]->alpha + curCell->neighCells[SOUTH]->alpha);
        } else if ((curCell->neighCells[WEST] == nullptr || curCell->neighCells[WEST]->bType_p == SOLID) && (curCell->neighCells[NORTH] == nullptr || curCell->neighCells[NORTH]->bType_p == SOLID)) {
            curCell->alpha = 0.5 * (curCell->neighCells[EAST]->alpha + curCell->neighCells[SOUTH]->alpha);
        } else if ((curCell->neighCells[EAST] == nullptr || curCell->neighCells[EAST]->bType_p == SOLID) && (curCell->neighCells[SOUTH] == nullptr || curCell->neighCells[SOUTH]->bType_p == SOLID)) {
            curCell->alpha = 0.5 * (curCell->neighCells[WEST]->alpha + curCell->neighCells[NORTH]->alpha);
        } else if ((curCell->neighCells[WEST] == nullptr || curCell->neighCells[WEST]->bType_p == SOLID) && (curCell->neighCells[SOUTH] == nullptr || curCell->neighCells[SOUTH]->bType_p == SOLID)) {
            curCell->alpha = 0.5 * (curCell->neighCells[EAST]->alpha + curCell->neighCells[NORTH]->alpha);
        }
    }
    
}

void handleExcessAlpha(Data2D& data){
    for (int i = 0; i < data.nCells; i++){
        Cell2D *curCell = &data.cells[i];
        if (curCell->alpha < 0){
            curCell->alpha = 0;
            // std::cout << "Negative alpha in cell " << curCell.id << std::endl;
        } else if (curCell->alpha > 1){
            curCell->alpha = 1;
            // std::cout << "Excess alpha in cell " << curCell.id << std::endl;
            // std::cout << "alpha: " << curCell.alpha << std::endl;
        }
    }
}

void assignInitasCorrVel(Data2D& data){
    for (int i = 0; i < data.nCells; i++){
        Cell2D *curCell = &data.cells[i];
        curCell->faces[WEST]->u[CORRECTED_2] = curCell->faces[WEST]->u[0];
        curCell->faces[EAST]->u[CORRECTED_2] = curCell->faces[EAST]->u[0];
        curCell->faces[SOUTH]->v[CORRECTED_2] = curCell->faces[SOUTH]->v[0];
        curCell->faces[NORTH]->v[CORRECTED_2] = curCell->faces[NORTH]->v[0];
    }
}

void resetNormalVectors(Data2D& data){
    for (int i = 0; i < data.nCells; i++){
        data.cells[i].normalVector[0] = 0;
        data.cells[i].normalVector[1] = 0;
    }
}

void advectAlpha(Data2D& data){

    //for testing
    assignInitasCorrVel(data);
    
    reconstructInterfaceLines(data);
    preformXSweep(data);
    resetNormalVectors(data);
    handleExcessAlpha(data);
    reconstructInterfaceLines(data);
    preformYSweep(data);
    handleExcessAlpha(data);
    resetNormalVectors(data);

}