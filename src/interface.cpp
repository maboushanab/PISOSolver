#include "interface.h"
#include "data.h"
#include <iomanip>
#include <limits>

/**
 * Calculates the normal vector and magnitude of the normal vector using a 2x2 stencil for corner cells and a 3x3 stencil for edge cells.
 * The normal vector is calculated using the formula:
 * 
 * @param data The 2D data set.
 * @param cellId The ID of the cell for which to calculate the normal vector.
 * @param pos The position of the cell in the grid.
 */
void calcAlphaNormalVectorBoundary(Data2D& data, int cellId, CellPosition pos){
    Cell2D *curCell = &data.cells[cellId];
    double nx, ny;
    double dx = curCell->faces[NORTH]->dx;
    double dy = curCell->faces[WEST]->dy; 
    Cell2D *a, *b, *c, *d, *e, *f;
    switch(pos) {
        case CellPosition::LeftEdge:
            a = curCell->neighCells[NORTH];
            b = curCell->neighCells[NORTH]->neighCells[EAST];
            c = curCell;
            d = curCell->neighCells[EAST];
            e = curCell->neighCells[SOUTH];
            f = curCell->neighCells[SOUTH]->neighCells[EAST];
            nx = (c->alpha - d->alpha) / dx;
            ny = (e->alpha - a->alpha) / (2 * dy);
            break;
        case CellPosition::RightEdge:
            a = curCell->neighCells[NORTH]->neighCells[WEST];
            b = curCell->neighCells[NORTH];
            c = curCell->neighCells[WEST];
            d = curCell;
            e = curCell->neighCells[SOUTH]->neighCells[WEST];
            f = curCell->neighCells[SOUTH];
            nx = (c->alpha - d->alpha) / dx;
            ny = (f->alpha - b->alpha) / (2 * dy);
            break;
        case CellPosition::TopEdge:
            a = curCell->neighCells[WEST];
            b = curCell;
            c = curCell->neighCells[EAST];
            d = curCell->neighCells[WEST]->neighCells[SOUTH];
            e = curCell->neighCells[SOUTH];
            f = curCell->neighCells[SOUTH]->neighCells[EAST];
            nx = (a->alpha - c->alpha) / (2 * dx);
            ny = (e->alpha - b->alpha) / dy;
            break;
        case CellPosition::BottomEdge:
            a = curCell->neighCells[WEST]->neighCells[NORTH];
            b = curCell->neighCells[NORTH];
            c = curCell->neighCells[EAST]->neighCells[NORTH];
            d = curCell->neighCells[WEST];
            e = curCell;
            f = curCell->neighCells[EAST];
            nx = (d->alpha - f->alpha) / (2 * dx);
            ny = (e->alpha - b->alpha) / dy;
            break;
        case CellPosition::TopLeft:
            a = curCell;
            b = curCell->neighCells[EAST];
            c = curCell->neighCells[SOUTH];
            nx = (a->alpha - b->alpha) / dx;
            ny = (c->alpha - a->alpha) / dy;
            break;
        case CellPosition::TopRight:
            a = curCell->neighCells[WEST];
            b = curCell;
            d = curCell->neighCells[SOUTH];
            nx = (a->alpha - b->alpha) / dx;
            ny = (d->alpha - a->alpha) / dy;
            break;
        case CellPosition::BottomLeft:
            a = curCell->neighCells[NORTH];
            c = curCell;
            d = curCell->neighCells[EAST];
            nx = (c->alpha - d->alpha) / dx;
            ny = (c->alpha - a->alpha) / dy;   
            break;
        case CellPosition::BottomRight:
            b = curCell->neighCells[NORTH];
            c = curCell->neighCells[WEST];
            d = curCell;
            nx = (c->alpha - d->alpha) / dx;
            ny = (d->alpha - b->alpha) / dy;
            break;
    }
    double mag = sqrt(nx*nx + ny*ny);
    curCell->normalVector[0] = nx/mag;
    curCell->normalVector[1] = ny/mag;
}

/**
 * Calculates the normal vector and magnitude of the normal vector using the Mixed Youngs-Centered Method using a 3x3 stencil.
 * The normal vector is calculated using the formula:
 * nx_C = (alpha_E - alpha_W) / (2 * dx)
 * ny_C = (alpha_N - alpha_S) / (2 * dy)
 * nx_Y = (1/(8*dx)) * (alpha_SW + 2*alpha_W + alpha_NW - alpha_SE - 2*alpha_E - alpha_NE)
 * ny_Y = (1/(8*dy)) * (alpha_NW + 2*alpha_N + alpha_NE - alpha_SW - 2*alpha_S - alpha_SE)
 * nx = (nx_C + nx_Y)/2
 * ny = (ny_C + ny_Y)/2
 * 
 * @param data The 2D data set.
 * @param cellId The ID of the cell for which to calculate the normal vector.
 */
void calcAlphaNormalVector3x3(Data2D& data, int cellId){
    Cell2D *curCell = &data.cells[cellId];
    double nx_C = (curCell->neighCells[EAST]->alpha - curCell->neighCells[WEST]->alpha) / (2 * curCell->faces[NORTH]->dx);
    double ny_C = (curCell->neighCells[NORTH]->alpha - curCell->neighCells[SOUTH]->alpha) / (2 *curCell->faces[WEST]->dy);
    double nx_Y = (1/(8*curCell->faces[NORTH]->dx)) 
                    * (curCell->neighCells[SOUTH]->neighCells[EAST]->alpha + 2*curCell->neighCells[EAST]->alpha + curCell->neighCells[NORTH]->neighCells[EAST]->alpha
                    - curCell->neighCells[SOUTH]->neighCells[WEST]->alpha - 2*curCell->neighCells[WEST]->alpha - curCell->neighCells[NORTH]->neighCells[WEST]->alpha);
    double ny_Y = (1/(8*curCell->faces[WEST]->dy)) 
                    * (curCell->neighCells[NORTH]->neighCells[WEST]->alpha + 2*curCell->neighCells[NORTH]->alpha + curCell->neighCells[NORTH]->neighCells[EAST]->alpha
                    - curCell->neighCells[SOUTH]->neighCells[WEST]->alpha - 2*curCell->neighCells[SOUTH]->alpha - curCell->neighCells[SOUTH]->neighCells[EAST]->alpha);
    double nx = (nx_C + nx_Y)/2;
    double ny = (ny_C + ny_Y)/2;
    double mag = sqrt(nx*nx + ny*ny);
    curCell->normalVector[0] = -nx/mag; 
    curCell->normalVector[1] = -ny/mag;
}

/**
 * Initializes the interface line of a given cell in the data structure.
 * The interface line is defined by the equation y = mx, where m is the slope and n is the y-intercept starting at the origin (cell center).
 * and x is the independent variable.
 *
 * @param data The data structure containing the cells.
 * @param cellId The ID of the cell to initialize the interface line for.
 */
void initInterface(Data2D& data, int cellId){
    Cell2D *curCell = &data.cells[cellId];
    curCell->interfaceVectorLine.origin[0] = 0;
    curCell->interfaceVectorLine.origin[1] = 0;
    curCell->interfaceVectorLine.direction[0] = curCell->normalVector[1];
    curCell->interfaceVectorLine.direction[1] = -curCell->normalVector[0];
}


void determinePolygonAsCell(Data2D& data, int cellId, std::vector<Eigen::Vector2d>& polygon){
    Cell2D *curCell = &data.cells[cellId];
    double dx = curCell->faces[NORTH]->dx;
    double dy = curCell->faces[WEST]->dy;
    std::vector<Eigen::Vector2d> cell = {
        Eigen::Vector2d(-dx/2, dy/2),
        Eigen::Vector2d(dx/2, dy/2),
        Eigen::Vector2d(dx/2, -dy/2),
        Eigen::Vector2d(-dx/2, -dy/2)
    }; 
    Eigen::Vector2d vertex;
    std::vector<Eigen::Vector2d> intersectionPoints;
    std::vector<Eigen::Vector2d> polygon1;
    std::vector<Eigen::Vector2d> polygon2;
    std::vector<bool> intersections;
    // potential intersection at the top face (y = dy/2 and x \in [-dx/2, dx/2])
    double lambdaTop = 0;
    if (curCell->normalVector[0] != 0){ // If the interface is not completely horizontal (to avoid division by zero)
        lambdaTop = -(dy/2 - curCell->interfaceVectorLine.origin[1])/curCell->normalVector[0];
    } else {
        lambdaTop = 1e6; // throws the point outside the interval
    }
    double x_top = curCell->interfaceVectorLine.origin[0] + lambdaTop*curCell->normalVector[1];
    if (x_top >= -dx/2 && x_top < dx/2){
        intersections.push_back(true);
    } else {
        intersections.push_back(false);
    }

    //  potential intersection at the right face (x = dx/2 and y \in [-dy/2, dy/2])
    double lambdaRight = 0;
    if (curCell->normalVector[1] != 0){ // If the interface is not completely vertical (to avoid division by zero)
        lambdaRight = (dx/2 - curCell->interfaceVectorLine.origin[0])/curCell->normalVector[1];
    } else {
        lambdaRight = 1e6; // throws the point outside the interval
    }
    double y_right = curCell->interfaceVectorLine.origin[1] - lambdaRight*curCell->normalVector[0];
    if (y_right > -dy/2 && y_right <= dy/2){
        intersections.push_back(true);
    } else {
        intersections.push_back(false);
    }

    // potential intersection at the bottom face (y = -dy/2 and x \in [-dx/2, dx/2])
    double lambdaBottom = 0;
    if (curCell->normalVector[0] != 0){ // If the interface is not completely horizontal (to avoid division by zero)
        lambdaBottom = -(-dy/2 - curCell->interfaceVectorLine.origin[1])/curCell->normalVector[0];
    } else {
        lambdaBottom = 1e6; // throws the point outside the interval
    }
    double x_bottom = curCell->interfaceVectorLine.origin[0] + lambdaBottom*curCell->normalVector[1];
    if (x_bottom > -dx/2 && x_bottom <= dx/2){
        intersections.push_back(true);
    } else {
        intersections.push_back(false);
    }

    // potential intersection at the left face (x = -dx/2 and y \in [-dy/2, dy/2])
    double lambdaLeft = 0;
    if (curCell->normalVector[1] != 0){ // If the interface is not completely vertical (to avoid division by zero)
        lambdaLeft = (-dx/2 - curCell->interfaceVectorLine.origin[0])/curCell->normalVector[1];
    } else {
        lambdaLeft = 1e6; // throws the point outside the interval
    }
    double y_left = curCell->interfaceVectorLine.origin[1] - lambdaLeft*curCell->normalVector[0];
    if (y_left >= -dy/2 && y_left < dy/2){
        intersections.push_back(true);
    } else {
        intersections.push_back(false);
    }
    int caseId = 0;
    // Combinations of intersections
    if (intersections[0] && intersections[1]){
        caseId = 1;
        intersectionPoints.push_back(Eigen::Vector2d(x_top, dy/2));
        intersectionPoints.push_back(Eigen::Vector2d(dx/2, y_right));
    } else if (intersections[0] && intersections[2]){
        caseId = 2;
        intersectionPoints.push_back(Eigen::Vector2d(x_top, dy/2));
        intersectionPoints.push_back(Eigen::Vector2d(x_bottom, -dy/2));
    } else if (intersections[0] && intersections[3]){
        caseId = 3;
        intersectionPoints.push_back(Eigen::Vector2d(x_top, dy/2));
        intersectionPoints.push_back(Eigen::Vector2d(-dx/2, y_left));
    } else if (intersections[1] && intersections[2]){
        caseId = 4;
        intersectionPoints.push_back(Eigen::Vector2d(dx/2, y_right));
        intersectionPoints.push_back(Eigen::Vector2d(x_bottom, -dy/2));
    } else if (intersections[1] && intersections[3]){
        caseId = 5;
        intersectionPoints.push_back(Eigen::Vector2d(dx/2, y_right));
        intersectionPoints.push_back(Eigen::Vector2d(-dx/2, y_left));
    } else if (intersections[2] && intersections[3]){
        caseId = 6;
        intersectionPoints.push_back(Eigen::Vector2d(x_bottom, -dy/2));
        intersectionPoints.push_back(Eigen::Vector2d(-dx/2, y_left));
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

bool checkIntersections(Data2D& data, int cellId){
    Cell2D *curCell = &data.cells[cellId];
    double dx = curCell->faces[NORTH]->dx;
    double dy = curCell->faces[WEST]->dy;
    std::vector<Eigen::Vector2d> cell = {
        Eigen::Vector2d(-dx/2, dy/2),
        Eigen::Vector2d(dx/2, dy/2),
        Eigen::Vector2d(dx/2, -dy/2),
        Eigen::Vector2d(-dx/2, -dy/2)
    }; 
    Eigen::Vector2d vertex;
    std::vector<Eigen::Vector2d> intersectionPoints;
    std::vector<Eigen::Vector2d> polygon1;
    std::vector<Eigen::Vector2d> polygon2;
    std::vector<bool> intersections;
    // potential intersection at the top face (y = dy/2 and x \in [-dx/2, dx/2])
    double lambdaTop = 0;
    if (curCell->normalVector[0] != 0){ // If the interface is not completely horizontal (to avoid division by zero)
        lambdaTop = -(dy/2 - curCell->interfaceVectorLine.origin[1])/curCell->normalVector[0];
    } else {
        lambdaTop = dx/2 + 1; // throws the point outside the interval
    }
    double x_top = curCell->interfaceVectorLine.origin[0] + lambdaTop*curCell->normalVector[1];
    if (x_top >= -dx/2 && x_top < dx/2){
        intersections.push_back(true);
    } else {
        intersections.push_back(false);
    }

    //  potential intersection at the right face (x = dx/2 and y \in [-dy/2, dy/2])
    double lambdaRight = 0;
    if (curCell->normalVector[1] != 0){ // If the interface is not completely vertical (to avoid division by zero)
        lambdaRight = (dx/2 - curCell->interfaceVectorLine.origin[0])/curCell->normalVector[1];
    } else {
        lambdaRight = dy/2 + 1; // throws the point outside the interval
    }
    double y_right = curCell->interfaceVectorLine.origin[1] - lambdaRight*curCell->normalVector[0];
    if (y_right > -dy/2 && y_right <= dy/2){
        intersections.push_back(true);
    } else {
        intersections.push_back(false);
    }

    // potential intersection at the bottom face (y = -dy/2 and x \in [-dx/2, dx/2])
    double lambdaBottom = 0;
    if (curCell->normalVector[0] != 0){ // If the interface is not completely horizontal (to avoid division by zero)
        lambdaBottom = -(dy/2 - curCell->interfaceVectorLine.origin[1])/curCell->normalVector[0];
    } else {
        lambdaBottom = dx/2 + 1; // throws the point outside the interval
    }
    double x_bottom = curCell->interfaceVectorLine.origin[0] + lambdaBottom*curCell->normalVector[1];
    if (x_bottom > -dx/2 && x_bottom <= dx/2){
        intersections.push_back(true);
    } else {
        intersections.push_back(false);
    }

    // potential intersection at the left face (x = -dx/2 and y \in [-dy/2, dy/2])
    double lambdaLeft = 0;
    if (curCell->normalVector[1] != 0){ // If the interface is not completely vertical (to avoid division by zero)
        lambdaLeft = (-dx/2 - curCell->interfaceVectorLine.origin[0])/curCell->normalVector[1];
    } else {
        lambdaLeft = dy/2 + 1; // throws the point outside the interval
    }
    double y_left = curCell->interfaceVectorLine.origin[1] - lambdaLeft*curCell->normalVector[0];
    if (y_left >= -dy/2 && y_left < dy/2){
        intersections.push_back(true);
    } else {
        intersections.push_back(false);
    }
    // if (!intersections[0] && !intersections[1] && !intersections[2] && !intersections[3]){
    //     std::cout << "x_top " << x_top << ", y_right " << y_right << ", x_bottom " << x_bottom << ", y_left " << y_left << std::endl;
    //     std::cout << "origin: [" << curCell->interfaceVectorLine.origin[0] << ", " << curCell->interfaceVectorLine.origin[1] << "]" << std::endl;
    //     std::cout << "normal: [" << curCell->normalVector[0] << ", " << curCell->normalVector[1] << "]" << std::endl;
    // }

    if (intersections[0] && intersections[1]){
        return true;
    } else if (intersections[0] && intersections[2]){
        return true;
    } else if (intersections[0] && intersections[3]){
        return true;
    } else if (intersections[1] && intersections[2]){
        return true;
    } else if (intersections[1] && intersections[3]){
        return true;
    } else if (intersections[2] && intersections[3]){
        return true;
    } else {
        return false;
    }
}

        

/**
 * Calculates the alpha prediction for a given cell in a 2D data set.
 * The alpha prediction is calculated by finding the area of the cell that is above the interface line. (Goldman Area of Planar Polygons)
 * First a stack is created to store the vertices of the cell. The vertices are pushed in a clockwise order starting from the top left corner.
 * Checking between every corner if the line intersects the face. If it does, the intersection point is added to the stack.
 * The area is then calculated by summing the cross products of the vertices using this formula:
 * area = (x1*y2 - x2*y1 + x2*y3 - x3*y2 + ... + xN*y1 - x1*yN) / (dx*dy)
 *
 * @param data The 2D data set.
 * @param cellId The ID of the cell for which to calculate the alpha prediction.
 * @param m The slope of the line.
 * @param n The y-intercept of the line.
 * @return The volume fraction prediction for the specified cell.
 */
double calcAlphaPrediction(Data2D& data, int cellId, bool print){
    Cell2D *curCell = &data.cells[cellId];
    double dx = curCell->faces[NORTH]->dx;
    double dy = curCell->faces[WEST]->dy;
    std::vector<Eigen::Vector2d> polygon;
    determinePolygonAsCell(data, cellId, polygon);
    std::stack<Eigen::Vector2d> vertices;
    // if (print){
    //     std::cout << "Polygon size: " << polygon.size() << " with n: " << n << std::endl;
    // }
    for (int i = 0; i < polygon.size(); i++){
        vertices.push(polygon[i]);
        // if (print){
        //     std::cout << "Polygon vertex: " << polygon[i](0) << " " << polygon[i](1) << std::endl;
        // }
    }
    //Goldman Area of Planar Polygons
    double area = 0;
    if (polygon.size() == 3){
        Eigen::Vector2d vertex1 = polygon[0];
        Eigen::Vector2d vertex2 = polygon[1];
        Eigen::Vector2d vertex3 = polygon[2];
        area = (vertex1(0)*vertex2(1) - vertex2(0)*vertex1(1) + vertex2(0)*vertex3(1) - vertex3(0)*vertex2(1) + vertex3(0)*vertex1(1) - vertex1(0)*vertex3(1))/2;
    } else if (polygon.size() == 4){
        Eigen::Vector2d vertex1 = polygon[0];
        Eigen::Vector2d vertex2 = polygon[1];
        Eigen::Vector2d vertex3 = polygon[2];
        Eigen::Vector2d vertex4 = polygon[3];
        area = (vertex1(0)*vertex2(1) - vertex2(0)*vertex1(1) + vertex2(0)*vertex3(1) - vertex3(0)*vertex2(1) + vertex3(0)*vertex4(1) - vertex4(0)*vertex3(1) + vertex4(0)*vertex1(1) - vertex1(0)*vertex4(1))/2;
    } else if (polygon.size() == 5){
        Eigen::Vector2d vertex1 = polygon[0];
        Eigen::Vector2d vertex2 = polygon[1];
        Eigen::Vector2d vertex3 = polygon[2];
        Eigen::Vector2d vertex4 = polygon[3];
        Eigen::Vector2d vertex5 = polygon[4];
        area = (vertex1(0)*vertex2(1) - vertex2(0)*vertex1(1) + vertex2(0)*vertex3(1) - vertex3(0)*vertex2(1) + vertex3(0)*vertex4(1) - vertex4(0)*vertex3(1) + vertex4(0)*vertex5(1) - vertex5(0)*vertex4(1) + vertex5(0)*vertex1(1) - vertex1(0)*vertex5(1))/2;
    }
    // for (int i = 0; i < polygon.size(); i++){
    //     area += (polygon[i](0)*polygon[(i+1)](1) - polygon[(i+1)](0)*polygon[i](1));
    // }
    // area += (polygon[polygon.size()-1](0)*polygon[0](1) - polygon[0](0)*polygon[polygon.size()-1](1));
        
    return abs(area)/(dx*dy);
}

/**
 * Calculates the sensitivity of the alpha parameter with respect to the input parameters.
 * The sensitivity is calculated using the formula:
 * sensitivity = (alpha(n + deltaN) - alpha(n)) / deltaN
 * 
 * @param data The 2D data set.
 * @param cellId The ID of the cell.
 * @param alpha Volume fraction value.
 * @return The sensitivity of the volume fraction parameter.
 */
double calcAlphaSensitivity(Data2D& data, int cellId, double alpha){
    double deltaN = 1e-6; // Small perturbation
    data.cells[cellId].interfaceVectorLine.origin[0] += deltaN*data.cells[cellId].normalVector[0];
    data.cells[cellId].interfaceVectorLine.origin[1] += deltaN*data.cells[cellId].normalVector[1];
    // std::cout << "Pertrubed basepoint: [" << data.cells[cellId].interfaceVectorLine.origin[0] << ", " <<  data.cells[cellId].interfaceVectorLine.origin[1] << "]" << std::endl;
    bool intersects = checkIntersections(data, cellId);
    if (!intersects){
        return 0;
    }
    double alphaPerturbed = calcAlphaPrediction(data, cellId, false); // Perturb n
    data.cells[cellId].interfaceVectorLine.origin[0] -= deltaN*data.cells[cellId].normalVector[0];
    data.cells[cellId].interfaceVectorLine.origin[1] -= deltaN*data.cells[cellId].normalVector[1];
    // std::cout << (alphaPerturbed - alpha) / deltaN << std::endl;
    return (alphaPerturbed - alpha) / deltaN; // Sensitivity
}

/**
 * Estimates the interface line for a given cell in the data.
 * The interface line is estimated by minimizing the difference between the predicted alpha value and the actual alpha value.
 * The estimation is done using the newton rhapson method.
 * The slope of the line is kept constant while the y-intercept is adjusted to minimize the difference.
 * The sensitivity of the alpha parameter with respect to the y-intercept is calculated and used to adjust the y-intercept.
 * The step size is adjusted based on the magnitude of the update.
 * The process is repeated until the difference between the predicted alpha value and the actual alpha value is minimized.
 * 
 * @param data The 2D data set.
 * @param cellId The ID of the cell for which to estimate the interface line.
 */
void estimateInterfaceLine(Data2D& data, int cellId){
    Cell2D *curCell = &data.cells[cellId];
    double dx = curCell->faces[NORTH]->dx;
    double dy = curCell->faces[WEST]->dy;
    initInterface(data, cellId);
    // std::cout << "origin: [" << curCell->interfaceVectorLine.origin[0] << ", " << curCell->interfaceVectorLine.origin[1] << "]" << std::endl;
    double alphaPrediction = calcAlphaPrediction(data, cellId, false);
    double alpha = curCell->alpha;
    double alphaDiff = alphaPrediction - alpha;
    int i = 0;
    double minUpdate = 1e-6; // Convergence criteria
    double gamma = 0.1; // Step size
    while (abs(alphaDiff) > 1e-6 && i < 10000){
        double sensitivity = calcAlphaSensitivity(data, cellId, alphaPrediction);
        // std::cout << "Sensitivity = " << sensitivity << std::endl;
        double oUpdate = alphaDiff / sensitivity;
        curCell->interfaceVectorLine.origin[0] -= oUpdate * curCell->normalVector[0];
        curCell->interfaceVectorLine.origin[1] -= oUpdate * curCell->normalVector[1];
        
        bool intersects = checkIntersections(data, cellId);
        // Correction if the line goes out of bounds
        while (!intersects){
            // std::cout << "Intersect correction" << std::endl;
            oUpdate = oUpdate/2;
            curCell->interfaceVectorLine.origin[0] += oUpdate * curCell->normalVector[0];
            curCell->interfaceVectorLine.origin[1] += oUpdate * curCell->normalVector[1];
            intersects = checkIntersections(data, cellId);
        }


        // Convergence check based on update magnitude
        if (abs(oUpdate) < minUpdate) {
            std::cout << "Break; oUpdate = " << oUpdate << std::endl; 
            break; // Converged due to small update
        }

        // Adjust gamma if needed
        if (abs(oUpdate) < 1e-5) {
            gamma *= 1.1; // Increase step size if changes are too small
        } else if (abs(oUpdate) > 0.1) {
            gamma *= 0.9; // Decrease step size if changes are too large
        }

        alphaPrediction = calcAlphaPrediction(data, cellId, false);
        alphaDiff = alphaPrediction - alpha;
        // if (i % 100 == 0){
            // std::cout << "Iteration: " << i << " alphaDiff: " << alphaDiff << std::endl;
        // }
        i++;
    }
    if (i == 10000){
        std::cout << "Failed to converge for cell " << cellId << " with difference : " << alphaDiff <<std::endl;
    } else {
        std::cout << "Converged for cell " << cellId << " in " << i << " iterations and alphaDiff: " << alphaDiff << std::endl;
    }
    calcAlphaPrediction(data, cellId, true);
}

/**
 * Estimates the interface lines for the given data.
 * This function calculates the alpha normal vector and estimates the interface line for each cell in the data.
 *
 * @param data The 2D data set.
 */
void reconstructInterfaceLines(Data2D& data){
    for (int i = 0; i < data.nCells; i++){
        Cell2D *curCell = &data.cells[i];
        if (curCell->alpha > 0 && curCell->alpha < 1){
            if (curCell->id == 0){               //Bottom Left
                calcAlphaNormalVectorBoundary(data, i, CellPosition::BottomLeft);
                estimateInterfaceLine(data, i);
                calcAlphaPrediction(data, i, false);
            } else if (curCell->id == data.dimX - 2){   //Bottom Right
                calcAlphaNormalVectorBoundary(data, i, CellPosition::BottomRight);
                estimateInterfaceLine(data, i);
                calcAlphaPrediction(data, i, false);
            } else if (curCell->id == data.nCells - (data.dimX - 1)){ //Top Left
                calcAlphaNormalVectorBoundary(data, i, CellPosition::TopLeft);
                estimateInterfaceLine(data, i);
                calcAlphaPrediction(data, i, false);
            } else if (curCell->id == data.nCells - 1){ //Top Right
                calcAlphaNormalVectorBoundary(data, i, CellPosition::TopRight);
                estimateInterfaceLine(data, i);
                calcAlphaPrediction(data, i, false);
            } else if (curCell->neighCells[NORTH] == nullptr){ //Top Edge
                calcAlphaNormalVectorBoundary(data, i, CellPosition::TopEdge);
                estimateInterfaceLine(data, i);
                calcAlphaPrediction(data, i, false);
            } else if (curCell->neighCells[WEST] == nullptr){ //Left Edge
                calcAlphaNormalVectorBoundary(data, i, CellPosition::LeftEdge);
                estimateInterfaceLine(data, i);
                calcAlphaPrediction(data, i, false);
            } else if (curCell->neighCells[EAST] == nullptr){ //Right Edge
                calcAlphaNormalVectorBoundary(data, i, CellPosition::RightEdge);
                estimateInterfaceLine(data, i);
                calcAlphaPrediction(data, i, false);
            } else if ( curCell->neighCells[SOUTH] == nullptr){ //Bottom Edge
                calcAlphaNormalVectorBoundary(data, i, CellPosition::BottomEdge);
                estimateInterfaceLine(data, i);
                calcAlphaPrediction(data, i, false);
            } else {
                calcAlphaNormalVector3x3(data, i);
                estimateInterfaceLine(data, i);
                calcAlphaPrediction(data, i, false);
            }
            // std::cout << "Cell " << i << " normal vector: " << std::setprecision(std::numeric_limits<double>::max_digits10) << curCell->normalVector[0] << " " << curCell->normalVector[1] << std::setprecision(5) << std::endl;
            // std::cout << "Cell " << i << " alpha: " << curCell->alpha << " alpha prediction: " << calcAlphaPrediction(data, i, false) << std::endl;
        }
    }
}