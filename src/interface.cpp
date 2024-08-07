#include "interface.h"
#include "data.h"

void setGhostCells(Data2D& data){
    for (int i=0; i<data.nCells; i++){
        Cell2D *curCell = &data.cells[i];
        if (curCell->bType_sc == DIRICHLET){
            if (curCell->neighCells[EAST] == nullptr){
                curCell->alpha_ghost = 2 * curCell->alpha - curCell->neighCells[WEST]->alpha; 
            } else if (curCell->neighCells[WEST] == nullptr){
                curCell->alpha_ghost = 2 * curCell->alpha - curCell->neighCells[EAST]->alpha; 
            } else if (curCell->neighCells[NORTH] == nullptr){
                curCell->alpha_ghost = 2 * curCell->alpha - curCell->neighCells[SOUTH]->alpha; 
            } else if (curCell->neighCells[SOUTH] == nullptr){
                curCell->alpha_ghost = 2 * curCell->alpha - curCell->neighCells[NORTH]->alpha;
            }
        } else if (curCell->bType_sc == NEUMANN){
            if (curCell->neighCells[EAST] == nullptr){
                curCell->alpha_ghost = curCell->neighCells[WEST]->alpha + 2 * curCell->faces[NORTH]->dy * curCell->g_sc;
            } else if (curCell->neighCells[WEST] == nullptr){
                curCell->alpha_ghost = curCell->neighCells[EAST]->alpha - 2 * curCell->faces[NORTH]->dy * curCell->g_sc;
            } else if (curCell->neighCells[NORTH] == nullptr){
                curCell->alpha_ghost = curCell->alpha + 2 * curCell->faces[WEST]->dx * curCell->g_sc;
            } else if (curCell->neighCells[SOUTH] == nullptr){
                curCell->alpha_ghost = curCell->alpha - 2 * curCell->faces[WEST]->dx * curCell->g_sc;
            }
        } 
    }
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
void calcAlphaNormalVector(Data2D& data, int cellId){
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
                    
    } else if (curCell->b_sc == DIRICHLET || curCell->bType_sc == NEUMANN){
        if (curCell->neighCells[EAST] == nullptr){
            nx_C = (curCell->alpha_ghost - curCell->neighCells[WEST]->alpha) / (2 * curCell->faces[NORTH]->dx);
            ny_C = (curCell->neighCells[NORTH]->alpha - curCell->neighCells[SOUTH]->alpha) / (2 *curCell->faces[WEST]->dy); 
            nx_Y = (1/(8*curCell->faces[NORTH]->dx)) 
                    * (curCell->neighCells[SOUTH]->neighCells[WEST]->alpha + 2*curCell->neighCells[WEST]->alpha + curCell->neighCells[NORTH]->neighCells[WEST]->alpha
                    - curCell->neighCells[SOUTH]->alpha_ghost - 2*curCell->alpha_ghost - curCell->neighCells[NORTH]->alpha_ghost);
            ny_Y = (1/(8*curCell->faces[WEST]->dy))
                    * (curCell->neighCells[NORTH]->neighCells[WEST]->alpha + 2*curCell->neighCells[NORTH]->alpha + curCell->neighCells[NORTH]->alpha_ghost
                    - curCell->neighCells[SOUTH]->neighCells[WEST]->alpha - 2*curCell->neighCells[SOUTH]->alpha - curCell->neighCells[SOUTH]->alpha_ghost);

        } else if (curCell->neighCells[WEST] == nullptr){
            nx_C = (curCell->neighCells[EAST]->alpha - curCell->alpha_ghost) / (2 * curCell->faces[NORTH]->dx);
            ny_C = (curCell->neighCells[NORTH]->alpha - curCell->neighCells[SOUTH]->alpha) / (2 *curCell->faces[WEST]->dy); 
            nx_Y = (1/(8*curCell->faces[NORTH]->dx)) 
                    * (curCell->neighCells[SOUTH]->alpha_ghost + 2*curCell->alpha_ghost + curCell->neighCells[NORTH]->alpha_ghost
                    - curCell->neighCells[SOUTH]->neighCells[EAST]->alpha - 2*curCell->neighCells[EAST]->alpha - curCell->neighCells[NORTH]->neighCells[EAST]->alpha);
            ny_Y = (1/(8*curCell->faces[WEST]->dy))
                    * (curCell->neighCells[NORTH]->neighCells[EAST]->alpha + 2*curCell->neighCells[NORTH]->alpha + curCell->neighCells[NORTH]->alpha_ghost
                    - curCell->neighCells[SOUTH]->neighCells[EAST]->alpha - 2*curCell->neighCells[SOUTH]->alpha - curCell->neighCells[SOUTH]->alpha_ghost);

        } else if (curCell->neighCells[NORTH] == nullptr){
            nx_C = (curCell->neighCells[EAST]->alpha - curCell->neighCells[WEST]->alpha) / (2 * curCell->faces[NORTH]->dx);
            ny_C = (curCell->alpha_ghost - curCell->neighCells[SOUTH]->alpha) / (2 *curCell->faces[WEST]->dy); 
            nx_Y = (1/(8*curCell->faces[NORTH]->dx)) 
                    * (curCell->neighCells[SOUTH]->neighCells[WEST]->alpha + 2*curCell->neighCells[WEST]->alpha + curCell->neighCells[WEST]->alpha_ghost
                    - curCell->neighCells[SOUTH]->neighCells[EAST]->alpha - 2*curCell->neighCells[EAST]->alpha - curCell->neighCells[EAST]->alpha_ghost);
            ny_Y = (1/(8*curCell->faces[WEST]->dy))
                    * (curCell->neighCells[WEST]->alpha_ghost + 2*curCell->alpha_ghost + curCell->neighCells[EAST]->alpha_ghost
                    - curCell->neighCells[SOUTH]->neighCells[WEST]->alpha - 2*curCell->neighCells[SOUTH]->alpha - curCell->neighCells[SOUTH]->neighCells[EAST]->alpha);

        } else if (curCell->neighCells[SOUTH] == nullptr){
            nx_C = (curCell->neighCells[EAST]->alpha - curCell->neighCells[WEST]->alpha) / (2 * curCell->faces[NORTH]->dx);
            ny_C = (curCell->neighCells[NORTH]->alpha - curCell->alpha_ghost) / (2 *curCell->faces[WEST]->dy); 
            nx_Y = (1/(8*curCell->faces[NORTH]->dx)) 
                    * (curCell->neighCells[WEST]->alpha_ghost + 2*curCell->neighCells[WEST]->alpha + curCell->neighCells[NORTH]->neighCells[WEST]->alpha
                    - curCell->neighCells[EAST]->alpha_ghost - 2*curCell->neighCells[EAST]->alpha - curCell->neighCells[NORTH]->neighCells[EAST]->alpha);
            ny_Y = (1/(8*curCell->faces[WEST]->dy))
                    * (curCell->neighCells[NORTH]->neighCells[WEST]->alpha + 2*curCell->neighCells[NORTH]->alpha + curCell->neighCells[NORTH]->neighCells[EAST]->alpha
                    - curCell->neighCells[WEST]->alpha_ghost - 2*curCell->alpha_ghost - curCell->neighCells[EAST]->alpha_ghost);
        }
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
            calcAlphaNormalVector(data, i);
            estimateInterfaceLine(data, i);
        }
    }
}