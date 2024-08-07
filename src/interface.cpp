#include "interface.h"
#include "data.h"

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
    switch(pos) {
        case CellPosition::LeftEdge:
            Cell2D *a = curCell->neighCells[NORTH];
            Cell2D *b = curCell->neighCells[NORTH]->neighCells[EAST];
            Cell2D *c = curCell;
            Cell2D *d = curCell->neighCells[WEST];
            Cell2D *e = curCell->neighCells[SOUTH];
            Cell2D *f = curCell->neighCells[SOUTH]->neighCells[EAST];
            nx = (c->alpha - d->alpha) / dx;
            ny = (e->alpha - a->alpha) / (2 * dy);
            break;
        case CellPosition::RightEdge:
            Cell2D *a = curCell->neighCells[NORTH]->neighCells[WEST];
            Cell2D *b = curCell->neighCells[NORTH];
            Cell2D *c = curCell->neighCells[EAST];
            Cell2D *d = curCell;
            Cell2D *e = curCell->neighCells[SOUTH]->neighCells[WEST];
            Cell2D *f = curCell->neighCells[SOUTH];
            nx = (d->alpha - c->alpha) / dx;
            ny = (f->alpha - b->alpha) / (2 * dy);
            break;
        case CellPosition::TopEdge:
            Cell2D *a = curCell->neighCells[WEST];
            Cell2D *b = curCell;
            Cell2D *c = curCell->neighCells[EAST];
            Cell2D *d = curCell->neighCells[WEST]->neighCells[SOUTH];
            Cell2D *e = curCell->neighCells[SOUTH];
            Cell2D *f = curCell->neighCells[SOUTH]->neighCells[EAST];
            nx = (c->alpha - a->alpha) / (2 * dx);
            ny = (b->alpha - e->alpha) / dy;
            break;
        case CellPosition::BottomEdge:
            Cell2D *a = curCell->neighCells[WEST]->neighCells[NORTH];
            Cell2D *b = curCell->neighCells[NORTH];
            Cell2D *c = curCell->neighCells[EAST]->neighCells[NORTH];
            Cell2D *d = curCell->neighCells[WEST];
            Cell2D *e = curCell;
            Cell2D *f = curCell->neighCells[EAST];
            nx = (f->alpha - d->alpha) / (2 * dx);
            ny = (e->alpha - b->alpha) / dy;
            break;
        case CellPosition::TopLeft:
            Cell2D *a = curCell;
            Cell2D *b = curCell->neighCells[EAST];
            Cell2D *c = curCell->neighCells[SOUTH];
            nx = (b->alpha - a->alpha) / dx;
            ny = (c->alpha - a->alpha) / dy;
            break;
        case CellPosition::TopRight:
            Cell2D *a = curCell->neighCells[WEST];
            Cell2D *b = curCell;
            Cell2D *d = curCell->neighCells[SOUTH];
            nx = (b->alpha - a->alpha) / dx;
            ny = (d->alpha - a->alpha) / dy;
            break;
        case CellPosition::BottomLeft:
            Cell2D *a = curCell->neighCells[NORTH];
            Cell2D *c = curCell;
            Cell2D *d = curCell->neighCells[EAST];
            nx = (d->alpha - c->alpha) / dx;
            ny = (c->alpha - a->alpha) / dy;   
            break;
        case CellPosition::BottomRight:
            Cell2D *b = curCell->neighCells[NORTH];
            Cell2D *c = curCell->neighCells[WEST];
            Cell2D *d = curCell;
            nx = (d->alpha - c->alpha) / dx;
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
    double nx_C = 0;
    double ny_C = 0;
    double nx_Y = 0;
    double ny_Y = 0;
    if (curCell->b_sc == INNERCELL){
        nx_C = (curCell->neighCells[EAST]->alpha - curCell->neighCells[WEST]->alpha) / (2 * curCell->faces[NORTH]->dx);
        ny_C = (curCell->neighCells[NORTH]->alpha - curCell->neighCells[SOUTH]->alpha) / (2 *curCell->faces[WEST]->dy);
        nx_Y = (1/(8*curCell->faces[NORTH]->dx)) 
                    * (curCell->neighCells[SOUTH]->neighCells[EAST]->alpha + 2*curCell->neighCells[EAST]->alpha + curCell->neighCells[NORTH]->neighCells[EAST]->alpha
                    - curCell->neighCells[SOUTH]->neighCells[WEST]->alpha - 2*curCell->neighCells[WEST]->alpha - curCell->neighCells[NORTH]->neighCells[WEST]->alpha);
        ny_Y = (1/(8*curCell->faces[WEST]->dy)) 
                    * (curCell->neighCells[NORTH]->neighCells[WEST]->alpha + 2*curCell->neighCells[NORTH]->alpha + curCell->neighCells[NORTH]->neighCells[EAST]->alpha
                    - curCell->neighCells[SOUTH]->neighCells[WEST]->alpha - 2*curCell->neighCells[SOUTH]->alpha - curCell->neighCells[SOUTH]->neighCells[EAST]->alpha);
    }


    double nx = (nx_C + nx_Y)/2;
    double ny = (ny_C + ny_Y)/2;
    double mag = sqrt(nx*nx + ny*ny);
    curCell->normalVector[0] = nx/mag; 
    curCell->normalVector[1] = ny/mag;
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
    curCell->interfaceLine.n = 0;
    curCell->interfaceLine.m = -(curCell->normalVector[0]/curCell->normalVector[1]); // m = -nx/ny
}

/**
 * Calculates the value of x using the formula x = (y - n)/m.
 *
 * @param y The value of y.
 * @param m The slope of the line.
 * @param n The y-intercept of the line.
 * @return The calculated value of x.
 */
double x(double y, double m, double n){
    return (y - n)/m;
}

/**
 * Calculates the value of y for a given x using the equation y = m*x + n.
 *
 * @param x The input value.
 * @param m The slope of the line.
 * @param n The y-intercept of the line.
 * @return The calculated value of y.
 */
double y(double x, double m, double n){
    return m*x + n;
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
double calcAlphaPrediction(Data2D& data, int cellId, double m, double n){
    Cell2D *curCell = &data.cells[cellId];
    double dx = curCell->faces[NORTH]->dx;
    double dy = curCell->faces[WEST]->dy;
    std::stack<Eigen::Vector2d> vertices;
    Eigen::Vector2d vertex;
    vertex(0) = -dx/2;
    vertex(1) = dy/2;
    vertices.push(vertex);
    if (x(dy/2, m, n) > -dx/2 && x(dy/2, m, n) < dx/2){ // top face 
        vertex(0) = x(dy/2, m, n);
        vertex(1) = dy/2;
        vertices.push(vertex);
    }
    vertex(0) = dx/2;
    vertex(1) = dy/2;
    vertices.push(vertex);
    if (y(dx/2, m, n) > -dy/2 && y(dx/2, m, n) < dy/2){ // right face
        vertex(0) = dx/2;
        vertex(1) = y(dx/2, m, n);
        vertices.push(vertex);
    }
    vertex(0) = dx/2;
    vertex(1) = -dy/2;
    vertices.push(vertex);
    if (x(-dy/2, m, n) > -dx/2 && x(-dy/2, m, n) < dx/2){ // bottom face
        vertex(0) = x(-dy/2, m, n);
        vertex(1) = -dy/2;
        vertices.push(vertex);
    }
    vertex(0) = -dx/2;
    vertex(1) = -dy/2;
    vertices.push(vertex);
    if (y(-dx/2, m, n) > -dy/2 && y(-dx/2, m, n) < dy/2){ // left face
        vertex(0) = -dx/2;
        vertex(1) = y(-dx/2, m, n);
        vertices.push(vertex);
    }
    vertex(0) = -dx/2;
    vertex(1) = dy/2;
    vertices.push(vertex);

    //Goldman Area of Planar Polygons
    double area = 0;
    Eigen::Vector2d vertex1;
    Eigen::Vector2d vertex2;
    while (vertices.size() > 1){
        vertex1 = vertices.top();
        vertices.pop();
        vertex2 = vertices.top();
        area += vertex1(0)*vertex2(1) - vertex2(0)*vertex1(1);
    }
    return area/(dx*dy);
}

/**
 * Calculates the sensitivity of the alpha parameter with respect to the input parameters.
 * The sensitivity is calculated using the formula:
 * sensitivity = (alpha(n + deltaN) - alpha(n)) / deltaN
 * 
 * @param data The 2D data set.
 * @param cellId The ID of the cell.
 * @param m The slope of the line.
 * @param n The y-intercept of the line.
 * @param alpha Volume fraction value.
 * @return The sensitivity of the volume fraction parameter.
 */
double calcAlphaSensitivity(Data2D& data, int cellId, double m, double n, double alpha){
    double deltaN = 1e-6; // Small perturbation
    double alphaPerturbed = calcAlphaPrediction(data, cellId, m, n + deltaN); // Perturb n
    return (alphaPerturbed - alpha) / deltaN; // Sensitivity
}

/**
 * Estimates the interface line for a given cell in the data.
 * The interface line is estimated by minimizing the difference between the predicted alpha value and the actual alpha value.
 * The estimation is done using the steepest descent method.
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
    double m = curCell->interfaceLine.m;
    double n = curCell->interfaceLine.n;

    double alphaPrediction = calcAlphaPrediction(data, cellId, m, n);
    double alpha = curCell->alpha;
    double alphaDiff = alphaPrediction - alpha;
    double maxIter = 0;
    double minUpdate = 1e-8; // Convergence criteria
    double gamma = 0.01; // Step size
    while (abs(alphaDiff) > 1e-6 && maxIter < 100){
        double sensitivity = calcAlphaSensitivity(data, cellId, m, n, alphaPrediction);
        double nUpdate = alphaDiff * sensitivity * gamma;
        n = n - nUpdate;

        // Convergence check based on update magnitude
        if (abs(nUpdate) < minUpdate) {
            break; // Converged due to small update
        }

        // Adjust gamma if needed
        if (abs(nUpdate) < 1e-7) {
            gamma *= 1.1; // Increase step size if changes are too small
        } else if (abs(nUpdate) > 0.1) {
            gamma *= 0.5; // Decrease step size if changes are too large
        }

        alphaPrediction = calcAlphaPrediction(data, cellId, m, n);
        alphaDiff = alphaPrediction - alpha;
        maxIter++;
    }
    curCell->interfaceLine.n = n;
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
            if (curCell->bType_sc == INNERCELL){
                calcAlphaNormalVector3x3(data, i);
                estimateInterfaceLine(data, i);
            } else if (curCell->id == 0){               //Top Left
                calcAlphaNormalVectorBoundary(data, i, CellPosition::TopLeft);
                estimateInterfaceLine(data, i);
            } else if (curCell->id == data.dimX - 1){   //Top Right
                calcAlphaNormalVectorBoundary(data, i, CellPosition::TopRight);
                estimateInterfaceLine(data, i);
            } else if (curCell->id == data.nCells - data.dimX){ //Bottom Left
                calcAlphaNormalVectorBoundary(data, i, CellPosition::BottomLeft);
                estimateInterfaceLine(data, i);
            } else if (curCell->id == data.nCells - 1){ //Bottom Right
                calcAlphaNormalVectorBoundary(data, i, CellPosition::BottomRight);
                estimateInterfaceLine(data, i);
            } else if (curCell->neighCells[NORTH] == nullptr){ //Top Edge
                calcAlphaNormalVectorBoundary(data, i, CellPosition::TopEdge);
                estimateInterfaceLine(data, i);
            } else if (curCell->neighCells[WEST] == nullptr){ //Left Edge
                calcAlphaNormalVectorBoundary(data, i, CellPosition::LeftEdge);
                estimateInterfaceLine(data, i);
            } else if (curCell->neighCells[EAST] == nullptr){ //Right Edge
                calcAlphaNormalVectorBoundary(data, i, CellPosition::RightEdge);
                estimateInterfaceLine(data, i);
            } else if ( curCell->neighCells[SOUTH] == nullptr){ //Bottom Edge
                calcAlphaNormalVectorBoundary(data, i, CellPosition::BottomEdge);
                estimateInterfaceLine(data, i);
            } 
        }
    }
}