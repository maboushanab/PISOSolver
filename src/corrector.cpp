#include "solve.h"
#include "data.h"
#include "corrector.h"
#include <unsupported/Eigen/IterativeSolvers>

/**
 * Sets the pressure matrix and vector for solving the pressure equation.
 * Ax = b
 * 
 * @param data The data structure containing the cell information.
 * @param pressureMatrix The sparse matrix representing the pressure coefficients. (A)
 * @param pressureVector The vector representing the pressure values. (b)
 * @param nCells The total number of cells.
 * @param step The current time step of the PISO algorithm.
 */

void setPressureMatrix(Data2D& data, SpMat& pressureMatrix, Vector& pressureVector, int nCells, int step) {
    for (int i = 0; i < nCells; i++) {
        Cell2D curCell = data.cells[i];
        if (curCell.bType_p == INNERCELL) {
            pressureVector(i) = -curCell.b;
            pressureMatrix.coeffRef(i, i) = -curCell.a_p;
            if (curCell.neighCells[EAST] != nullptr) {
                pressureMatrix.coeffRef(i, curCell.neighCells[EAST]->id) = curCell.a_e;
            }
            if (curCell.neighCells[WEST] != nullptr) {
                pressureMatrix.coeffRef(i, curCell.neighCells[WEST]->id) = curCell.a_w;
            }
            if (curCell.neighCells[NORTH] != nullptr) {
                pressureMatrix.coeffRef(i, curCell.neighCells[NORTH]->id) = curCell.a_n;
            }
            if (curCell.neighCells[SOUTH] != nullptr) {
                pressureMatrix.coeffRef(i, curCell.neighCells[SOUTH]->id) = curCell.a_s;
            }
        } else if (curCell.bType_p == DIRICHLET || curCell.bType_p == SOLID) {
            pressureMatrix.coeffRef(i, i) = 1.0;
            pressureVector(i) = curCell.p[step - 1];
        } else if (curCell.bType_p == NEUMANN) {
            pressureMatrix.coeffRef(i, i) = 1.0;
            if ((curCell.neighCells[EAST] == nullptr || curCell.neighCells[EAST]->bType_sc == SOLID) && (curCell.neighCells[NORTH] == nullptr || curCell.neighCells[NORTH]->bType_sc == SOLID)) {
                pressureMatrix.coeffRef(i, curCell.neighCells[WEST]->id) = -0.5;
                pressureMatrix.coeffRef(i, curCell.neighCells[SOUTH]->id) = -0.5;
                pressureVector(i) = 0.5 * (curCell.g_p*curCell.faces[EAST]->dx + curCell.g_p*curCell.faces[NORTH]->dy);
            } else if ((curCell.neighCells[WEST] == nullptr || curCell.neighCells[WEST]->bType_sc == SOLID) && (curCell.neighCells[NORTH] == nullptr || curCell.neighCells[NORTH]->bType_sc == SOLID)) {
                pressureMatrix.coeffRef(i, curCell.neighCells[EAST]->id) = -0.5;
                pressureMatrix.coeffRef(i, curCell.neighCells[SOUTH]->id) = -0.5;
                pressureVector(i) = -0.5 * (curCell.g_p*curCell.faces[WEST]->dx + curCell.g_p*curCell.faces[NORTH]->dy);
            } else if ((curCell.neighCells[EAST] == nullptr || curCell.neighCells[EAST]->bType_sc == SOLID) && (curCell.neighCells[SOUTH] == nullptr || curCell.neighCells[SOUTH]->bType_sc == SOLID)) {
                pressureMatrix.coeffRef(i, curCell.neighCells[WEST]->id) = -0.5;
                pressureMatrix.coeffRef(i, curCell.neighCells[NORTH]->id) = -0.5;
                pressureVector(i) = 0.5 * (curCell.g_p*curCell.faces[EAST]->dx - curCell.g_p*curCell.faces[SOUTH]->dy);
            } else if ((curCell.neighCells[WEST] == nullptr || curCell.neighCells[WEST]->bType_sc == SOLID) && (curCell.neighCells[SOUTH] == nullptr || curCell.neighCells[SOUTH]->bType_sc == SOLID)) {
                pressureMatrix.coeffRef(i, curCell.neighCells[EAST]->id) = -0.5;
                pressureMatrix.coeffRef(i, curCell.neighCells[NORTH]->id) = -0.5;
                pressureVector(i) = -0.5 * (curCell.g_p*curCell.faces[WEST]->dx - curCell.g_p*curCell.faces[SOUTH]->dy);
            } else if (curCell.neighCells[EAST] == nullptr || curCell.neighCells[EAST]->bType_sc == SOLID) {
                pressureMatrix.coeffRef(i, curCell.neighCells[WEST]->id) = -1.0;
                pressureVector(i) = curCell.g_p*curCell.faces[EAST]->dx;
            } else if (curCell.neighCells[WEST] == nullptr || curCell.neighCells[WEST]->bType_sc == SOLID) {
                pressureMatrix.coeffRef(i, curCell.neighCells[EAST]->id) = -1.0;
                pressureVector(i) = -curCell.g_p*curCell.faces[WEST]->dx;
            } else if (curCell.neighCells[NORTH] == nullptr || curCell.neighCells[NORTH]->bType_sc == SOLID) {
                pressureMatrix.coeffRef(i, curCell.neighCells[SOUTH]->id) = -1.0;
                pressureVector(i) = curCell.g_p*curCell.faces[NORTH]->dy;
            } else if (curCell.neighCells[SOUTH] == nullptr || curCell.neighCells[SOUTH]->bType_sc == SOLID) {
                pressureMatrix.coeffRef(i, curCell.neighCells[NORTH]->id) = -1.0;
                pressureVector(i) = -curCell.g_p*curCell.faces[SOUTH]->dy;
            }
        }
        
        // // set a central node in the middle of the grid to dirichlet with p = 0
        // if (i == std::round(nCells / 2)) {
        //     for (int j = 0; j < nCells; j++) {
        //         if (j != i) {
        //             pressureMatrix.coeffRef(i, j) = 0.0;
        //         }
        //     }
        //     pressureMatrix.coeffRef(i, i) = 1.0;
        //     pressureVector(i) = 0.0;
        // }

    }
}


/**
 * Computes the pressure coefficients for a given cell.
 * 
 * @param data The data structure containing the simulation data.
 * @param cellId The ID of the cell for which to compute the pressure coefficients.
 * @param step The current time step of the PISO algorithm.
 */
void computePressureCoeff(Data2D& data, int cellId){
    Cell2D *curCell = &data.cells[cellId];	
    double dx = curCell->faces[SOUTH]->dx;
    double dy = curCell->faces[WEST]->dy;

    curCell->a_e = (fRho(data, curCell->neighCells[EAST]->alpha) * dy * dy) / curCell->faces[EAST]->a_p_tilde;
    curCell->a_w = (fRho(data, curCell->neighCells[WEST]->alpha) * dy * dy) / curCell->faces[WEST]->a_p_tilde;
    curCell->a_n = (fRho(data, curCell->neighCells[NORTH]->alpha) * dx * dx) / curCell->faces[NORTH]->a_p_tilde;
    curCell->a_s = (fRho(data, curCell->neighCells[SOUTH]->alpha) * dx * dx) / curCell->faces[SOUTH]->a_p_tilde;
    // curCell->a_e = -dy/(dx*fRho(data, curCell->neighCells[EAST]->alpha));
    // curCell->a_w = -dy/(dx*fRho(data, curCell->neighCells[WEST]->alpha));
    // curCell->a_n = -dx/(dy*fRho(data, curCell->neighCells[NORTH]->alpha));
    // curCell->a_s = -dx/(dy*fRho(data, curCell->neighCells[SOUTH]->alpha));
    curCell->a_p = -curCell->a_e - curCell->a_w - curCell->a_n - curCell->a_s;
    curCell->b = - ((fRho(data, curCell->alpha) + fRho(data, curCell->neighCells[EAST]->alpha)) * 0.5 * curCell->faces[EAST]->u[INTERMEDIATE_1] * dy) + ((fRho(data, curCell->alpha) + fRho(data, curCell->neighCells[WEST]->alpha)) * 0.5 *curCell->faces[WEST]->u[INTERMEDIATE_1] * dy) 
    - ((fRho(data, curCell->alpha) + fRho(data, curCell->neighCells[NORTH]->alpha))* 0.5 *curCell->faces[NORTH]->v[INTERMEDIATE_1] * dx) + ((fRho(data, curCell->alpha) + fRho(data, curCell->neighCells[SOUTH]->alpha))* 0.5 *curCell->faces[SOUTH]->v[INTERMEDIATE_1] * dx);
    // curCell->b = (curCell->faces[EAST]->u[INTERMEDIATE_1] - curCell->faces[WEST]->u[INTERMEDIATE_1]) * dy + (curCell->faces[NORTH]->v[INTERMEDIATE_1] - curCell->faces[SOUTH]->v[INTERMEDIATE_1]) * dx;
    if (data.mode == 0){
        curCell->b += (fRho(data, curCell->alpha_prev) - fRho(data, curCell->alpha)) * dx * dy / data.dt;
    }
}

/**
 * Solves the pressure equation for a given data structure.
 * This function updates the pressure field in the data structure by solving a linear system of equations
 *
 * @param data The data structure containing the grid and cell information.
 * @param step The current step of the PISO algorithm.
 */

void correctPressureEquationBiCGStab(Data2D& data, int step) {
    SpMat pressureCorrMatrix(data.nCells, data.nCells);
    Vector pressureCorrVector(data.nCells);
    pressureCorrMatrix.setZero();
    pressureCorrVector.setZero();

    for (int i = 0; i < data.nCells; i++) {
        Cell2D *curCell = &data.cells[i];
        if (curCell->bType_p == INNERCELL) {
            if (step == INTERMEDIATE_1){
                computePressureCoeff(data, i);
            } else if (step == INTERMEDIATE_2){
                computePressureCoeff2(data, i);
            }
    
        }
    }

    setPressureMatrix(data, pressureCorrMatrix, pressureCorrVector, data.nCells, step);
    // Eigen::MatrixXd denseA = Eigen::MatrixXd(pressureCorrMatrix);
    // Eigen::JacobiSVD<Eigen::MatrixXd> svd(denseA);
    // double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);
    // std::cout << "Condition number of pressure matrix: " << cond << std::endl;

    pressureCorrMatrix.makeCompressed();
    // //calculate condition number of matrix


    try {
        // Solve for pressure correction
        // with preconditioner
        Eigen::BiCGSTAB<SpMat, Eigen::IncompleteLUT<double>> solver;
        // // solve with gmres
        // Eigen::GMRES<SpMat, Eigen::IdentityPreconditioner> solver;
        solver.compute(pressureCorrMatrix);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Decomposition failed for pressureCorrMatrix");
        }
        Vector pressureCorrSolution = solver.solve(pressureCorrVector);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Solving failed for pressureCorrMatrix");
        }
        for (int i = 0; i < data.nCells; i++) {
            data.cells[i].p[step] = pressureCorrSolution(i);
        }
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

void correctPressureEquationSparseLU(Data2D& data, int step) {
    SpMat pressureCorrMatrix(data.nCells, data.nCells);
    Vector pressureCorrVector(data.nCells);
    pressureCorrMatrix.setZero();
    pressureCorrVector.setZero();

    for (int i = 0; i < data.nCells; i++) {
        Cell2D *curCell = &data.cells[i];
        if (curCell->bType_p == INNERCELL) {
            if (step == INTERMEDIATE_1){
                computePressureCoeff(data, i);
            } else if (step == INTERMEDIATE_2){
                computePressureCoeff2(data, i);
            }
    
        }
    }

    setPressureMatrix(data, pressureCorrMatrix, pressureCorrVector, data.nCells, step);
    pressureCorrMatrix.makeCompressed();

    try {
        // Solve for pressure correction
        Eigen::SparseLU<SpMat> solver;
        solver.analyzePattern(pressureCorrMatrix);
        solver.factorize(pressureCorrMatrix);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Factorization failed for pressureCorrMatrix");
        }
        Vector pressureCorrSolution = solver.solve(pressureCorrVector);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Solving failed for pressureCorrMatrix");
        }
        for (int i = 0; i < data.nCells; i++) {
            data.cells[i].p[step] = pressureCorrSolution(i);
        }
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

void correctPressureEquation(Data2D& data, int step) {
    if (data.presSolver == 1) {
        correctPressureEquationSparseLU(data, step);
    } else {
        correctPressureEquationBiCGStab(data, step);
    }
    //fix small domain around central node with cell id = nCells/2
    // recomputePressureSmallDomain(data, step);
    
}

/**
 * Applies the pressure correction and velocity field correction to the given data structure.
 * This function updates the pressure field and velocity field in the data structure based on the correction calculations.
 *
 * @param data The data structure containing the cells and faces.
 * @param step The current step of the PISO algorithm.
 */
void corrector1(Data2D& data) {
    // Update Pressure bzw. cells
    data.alpha_p_relax = 1;
    for (int i = 0; i < data.nCells; i++) {
        Cell2D *curCell = &data.cells[i];
        if (curCell->bType_p == INNERCELL) {
            curCell->p[CORRECTED_1] = curCell->p[INITIAL] + curCell->p[INTERMEDIATE_1] * data.alpha_p_relax;
            // std::cout << "Cell " << i << "; p: " << curCell->p[step + 1] << std::endl;
        }
    }
    // Update Boundary Conditions
    std::stack<Cell2D*> cornerCells; //saves the corner cells to update them after the other cells for average pressure
    for (int i = 0; i < data.nCells; i++){
        Cell2D *curCell = &data.cells[i];
        if (curCell->bType_p == DIRICHLET || curCell->bType_p == SOLID) {
            curCell->p[CORRECTED_1] = curCell->p[INITIAL];
        } else if (curCell->bType_p == NEUMANN){
            if ((curCell->neighCells[EAST] == nullptr || curCell->neighCells[EAST]->bType_sc == SOLID) && (curCell->neighCells[NORTH] == nullptr || curCell->neighCells[NORTH]->bType_sc == SOLID)) {
                cornerCells.push(curCell);
            } else if ((curCell->neighCells[WEST] == nullptr || curCell->neighCells[WEST]->bType_sc == SOLID) && (curCell->neighCells[NORTH] == nullptr || curCell->neighCells[NORTH]->bType_sc == SOLID)) {
                cornerCells.push(curCell);
            } else if ((curCell->neighCells[EAST] == nullptr || curCell->neighCells[EAST]->bType_sc == SOLID) && (curCell->neighCells[SOUTH] == nullptr || curCell->neighCells[SOUTH]->bType_sc == SOLID)) {
                cornerCells.push(curCell);
            } else if ((curCell->neighCells[WEST] == nullptr || curCell->neighCells[WEST]->bType_sc == SOLID) && (curCell->neighCells[SOUTH] == nullptr || curCell->neighCells[SOUTH]->bType_sc == SOLID)) {
                cornerCells.push(curCell);
            } else if (curCell->neighCells[EAST] == nullptr || curCell->neighCells[EAST]->bType_sc == SOLID) {
                curCell->p[CORRECTED_1] = curCell->neighCells[WEST]->p[CORRECTED_1] + curCell->g_p * curCell->faces[EAST]->dx;
            } else if (curCell->neighCells[WEST] == nullptr || curCell->neighCells[WEST]->bType_sc == SOLID) {
                curCell->p[CORRECTED_1] = curCell->neighCells[EAST]->p[CORRECTED_1] - curCell->g_p * curCell->faces[WEST]->dx;
            } else if (curCell->neighCells[NORTH] == nullptr || curCell->neighCells[NORTH]->bType_sc == SOLID) {
                curCell->p[CORRECTED_1] = curCell->neighCells[SOUTH]->p[CORRECTED_1] + curCell->g_p * curCell->faces[NORTH]->dy;
            } else if (curCell->neighCells[SOUTH] == nullptr || curCell->neighCells[SOUTH]->bType_sc == SOLID) {
                curCell->p[CORRECTED_1] = curCell->neighCells[NORTH]->p[CORRECTED_1] - curCell->g_p * curCell->faces[SOUTH]->dy;
            }
        }
    }
    // Update Corner Cells
    while (!cornerCells.empty()){
        Cell2D *curCell = cornerCells.top();
        cornerCells.pop();
        if ((curCell->neighCells[EAST] == nullptr || curCell->neighCells[EAST]->bType_sc == SOLID) && (curCell->neighCells[NORTH] == nullptr || curCell->neighCells[NORTH]->bType_sc == SOLID)) {
            curCell->p[CORRECTED_1] = 0.5 * (curCell->neighCells[WEST]->p[CORRECTED_1] + curCell->neighCells[SOUTH]->p[CORRECTED_1]);
        } else if ((curCell->neighCells[WEST] == nullptr || curCell->neighCells[WEST]->bType_sc == SOLID) && (curCell->neighCells[NORTH] == nullptr || curCell->neighCells[NORTH]->bType_sc == SOLID)) {
            curCell->p[CORRECTED_1] = 0.5 * (curCell->neighCells[EAST]->p[CORRECTED_1] + curCell->neighCells[SOUTH]->p[CORRECTED_1]);
        } else if ((curCell->neighCells[EAST] == nullptr || curCell->neighCells[EAST]->bType_sc == SOLID) && (curCell->neighCells[SOUTH] == nullptr || curCell->neighCells[SOUTH]->bType_sc == SOLID)) {
            curCell->p[CORRECTED_1] = 0.5 * (curCell->neighCells[WEST]->p[CORRECTED_1] + curCell->neighCells[NORTH]->p[CORRECTED_1]);
        } else if ((curCell->neighCells[WEST] == nullptr || curCell->neighCells[WEST]->bType_sc == SOLID) && (curCell->neighCells[SOUTH] == nullptr || curCell->neighCells[SOUTH]->bType_sc == SOLID)) {
            curCell->p[CORRECTED_1] = 0.5 * (curCell->neighCells[EAST]->p[CORRECTED_1] + curCell->neighCells[NORTH]->p[CORRECTED_1]);
        }
    }
        
    // Update Velocities bzw. faces
    for (int i = 0; i < data.nFaces; i++) {
        Face2D *curFace = &data.faces[i];
        if (curFace->bType_u == INNERCELL) {
            if (i < data.nhorizontalFaces) {
                curFace->v[CORRECTED_1] = curFace->dx/curFace->a_p_tilde * (curFace->neighCells[DOWN]->p[INTERMEDIATE_1] - curFace->neighCells[UP]->p[INTERMEDIATE_1]);
                curFace->v[INTERMEDIATE_2] = curFace->v[CORRECTED_1] + curFace->v[INTERMEDIATE_1];
            } else if (i >= data.nhorizontalFaces) {
                curFace->u[CORRECTED_1] = curFace->dy/curFace->a_p_tilde * (curFace->neighCells[LEFT]->p[INTERMEDIATE_1] - curFace->neighCells[RIGHT]->p[INTERMEDIATE_1]);
                curFace->u[INTERMEDIATE_2] = curFace->u[CORRECTED_1] + curFace->u[INTERMEDIATE_1];
            }   
        } 
    }
    // Update Boundary Conditions
    for (int i = 0; i < data.nFaces; i++) {
        Face2D *curFace = &data.faces[i];
        if (curFace->bType_u == DIRICHLET || curFace->bType_u == SOLID){
            if (i < data.nhorizontalFaces) {
                curFace->v[CORRECTED_1] = curFace->v[INITIAL];
                curFace->v[INTERMEDIATE_2] = curFace->v[INITIAL];
            } else if (i >= data.nhorizontalFaces) {
                curFace->u[CORRECTED_1] = curFace->u[INITIAL];
                curFace->u[INTERMEDIATE_2] = curFace->u[INITIAL];
            }
        } else if (curFace->bType_u == NEUMANN) {
            if (i < data.nhorizontalFaces) {
                if(curFace->neighCells[UP] == nullptr || curFace->neighCells[UP]->bType_sc == SOLID)                                                     //TOP BOUNDARY (HORIZONTAL)
                {
                    curFace->v[CORRECTED_1] = curFace->neighCells[DOWN]->faces[SOUTH]->v[CORRECTED_1] + curFace->g_v * curFace->dx;
                    curFace->v[INTERMEDIATE_2] = curFace->neighCells[DOWN]->faces[SOUTH]->v[INTERMEDIATE_2] + curFace->g_v * curFace->dx;
                }
                else if (curFace->neighCells[DOWN] == nullptr || curFace->neighCells[DOWN]->bType_sc == SOLID)                                           //BOTTOM BOUNDARY (HORIZONTAL) 
                {
                    curFace->v[CORRECTED_1] = curFace->neighCells[UP]->faces[NORTH]->v[CORRECTED_1] - curFace->g_v * curFace->dx;
                    curFace->v[INTERMEDIATE_2] = curFace->neighCells[UP]->faces[NORTH]->v[INTERMEDIATE_2] - curFace->g_v * curFace->dx;
                }
            } else {
                if (curFace->neighCells[LEFT] == nullptr || curFace->neighCells[LEFT]->bType_sc == SOLID)                                                 //LEFT BOUNDARY (VERTICAL) 
                {
                    curFace->u[CORRECTED_1] = curFace->neighCells[RIGHT]->faces[EAST]->u[CORRECTED_1] - curFace->g_u * curFace->dy;
                    curFace->u[INTERMEDIATE_2] = curFace->neighCells[RIGHT]->faces[EAST]->u[INTERMEDIATE_2] - curFace->g_u * curFace->dy;
                }
                else if (curFace->neighCells[RIGHT] == nullptr || curFace->neighCells[RIGHT]->bType_sc == SOLID)                                         //RIGHT BOUNDARY (VERTICAL)
                {
                    curFace->u[CORRECTED_1] = curFace->neighCells[LEFT]->faces[WEST]->u[CORRECTED_1] + curFace->g_u * curFace->dy;
                    curFace->u[INTERMEDIATE_2] = curFace->neighCells[LEFT]->faces[WEST]->u[INTERMEDIATE_2] + curFace->g_u * curFace->dy;
                }
            }
        }
    }
}

void computePressureCoeff2(Data2D& data, int cellId){
    Cell2D *curCell = &data.cells[cellId];
    double dx = curCell->faces[SOUTH]->dx;
    double dy = curCell->faces[WEST]->dy;

    double b_e = (0.5 * (fRho(data, curCell->neighCells[EAST]->alpha) + fRho(data, curCell->alpha)) * dy) / curCell->faces[EAST]->a_p_tilde
                * (curCell->neighCells[EAST]->faces[EAST]->u[CORRECTED_1] * curCell->neighCells[EAST]->faces[EAST]->a_p_tilde
                + curCell->faces[WEST]->u[CORRECTED_1] * curCell->faces[WEST]->a_p_tilde
                + curCell->neighCells[NORTH]->faces[EAST]->u[CORRECTED_1] * curCell->neighCells[NORTH]->faces[EAST]->a_p_tilde
                + curCell->neighCells[SOUTH]->faces[EAST]->u[CORRECTED_1] * curCell->neighCells[SOUTH]->faces[EAST]->a_p_tilde); 

    double b_w = (0.5 * (fRho(data, curCell->neighCells[WEST]->alpha) + fRho(data, curCell->alpha)) * dy) / curCell->faces[WEST]->a_p_tilde
                * (curCell->faces[EAST]->u[CORRECTED_1] * curCell->faces[EAST]->a_p_tilde
                + curCell->neighCells[WEST]->faces[WEST]->u[CORRECTED_1] * curCell->neighCells[WEST]->faces[WEST]->a_p_tilde
                + curCell->neighCells[NORTH]->faces[WEST]->u[CORRECTED_1] * curCell->neighCells[NORTH]->faces[WEST]->a_p_tilde
                + curCell->neighCells[SOUTH]->faces[WEST]->u[CORRECTED_1] * curCell->neighCells[SOUTH]->faces[WEST]->a_p_tilde);

    double b_n = (0.5 * (fRho(data, curCell->neighCells[NORTH]->alpha) + fRho(data, curCell->alpha)) * dx) / curCell->faces[NORTH]->a_p_tilde
                * (curCell->neighCells[NORTH]->faces[NORTH]->v[CORRECTED_1] * curCell->neighCells[NORTH]->faces[NORTH]->a_p_tilde
                + curCell->neighCells[WEST]->faces[NORTH]->v[CORRECTED_1] * curCell->neighCells[WEST]->faces[NORTH]->a_p_tilde
                + curCell->neighCells[EAST]->faces[NORTH]->v[CORRECTED_1] * curCell->neighCells[EAST]->faces[NORTH]->a_p_tilde
                + curCell->faces[SOUTH]->v[CORRECTED_1] * curCell->faces[SOUTH]->a_p_tilde);

    double b_s = (0.5 * (fRho(data, curCell->neighCells[SOUTH]->alpha) + fRho(data, curCell->alpha)) * dx) / curCell->faces[SOUTH]->a_p_tilde
                * (curCell->faces[NORTH]->v[CORRECTED_1] * curCell->faces[NORTH]->a_p_tilde
                + curCell->neighCells[WEST]->faces[SOUTH]->v[CORRECTED_1] * curCell->neighCells[WEST]->faces[SOUTH]->a_p_tilde
                + curCell->neighCells[EAST]->faces[SOUTH]->v[CORRECTED_1] * curCell->neighCells[EAST]->faces[SOUTH]->a_p_tilde
                + curCell->faces[SOUTH]->v[CORRECTED_1] * curCell->faces[SOUTH]->a_p_tilde);
    
    curCell->b = -b_e + b_w - b_n + b_s;
    if (data.mode == 0){
        curCell->b += (fRho(data, curCell->alpha_prev) - fRho(data, curCell->alpha)) * dx * dy / data.dt;
    }

    curCell->a_e = (fRho(data, curCell->neighCells[EAST]->alpha) * dy * dy) / curCell->faces[EAST]->a_p_tilde;
    curCell->a_w = (fRho(data, curCell->neighCells[WEST]->alpha) * dy * dy) / curCell->faces[WEST]->a_p_tilde;
    curCell->a_n = (fRho(data, curCell->neighCells[NORTH]->alpha) * dx * dx) / curCell->faces[NORTH]->a_p_tilde;
    curCell->a_s = (fRho(data, curCell->neighCells[SOUTH]->alpha) * dx * dx) / curCell->faces[SOUTH]->a_p_tilde;
    // curCell->a_e = -dy/(dx*fRho(data, curCell->neighCells[EAST]->alpha));
    // curCell->a_w = -dy/(dx*fRho(data, curCell->neighCells[WEST]->alpha));
    // curCell->a_n = -dx/(dy*fRho(data, curCell->neighCells[NORTH]->alpha));
    // curCell->a_s = -dx/(dy*fRho(data, curCell->neighCells[SOUTH]->alpha));
    curCell->a_p = -curCell->a_e - curCell->a_w - curCell->a_n - curCell->a_s;

    // curCell->b = (curCell->faces[EAST]->u[CORRECTED_1] - curCell->faces[WEST]->u[CORRECTED_1]) / dx + (curCell->faces[NORTH]->v[CORRECTED_1] - curCell->faces[SOUTH]->v[CORRECTED_1]) / dy;
}

void corrector2(Data2D& data){
    data.alpha_p_relax = 1;
    for (int i = 0; i < data.nCells; i++) {
        Cell2D *curCell = &data.cells[i];
        if (curCell->bType_p == INNERCELL) {
            curCell->p[CORRECTED_2] = curCell->p[CORRECTED_1] + curCell->p[INTERMEDIATE_2] * data.alpha_p_relax;
        }
    }
    // Update Boundary Conditions 
    std::stack<Cell2D*> cornerCells; //saves the corner cells to update them after the other cells for average pressure
    for (int i = 0; i < data.nCells; i++){
        Cell2D *curCell = &data.cells[i];
        if (curCell->bType_p == DIRICHLET || curCell->bType_p == SOLID) {
            curCell->p[CORRECTED_2] = curCell->p[INITIAL];
        } else if (curCell->bType_p == NEUMANN){
            if ((curCell->neighCells[EAST] == nullptr || curCell->neighCells[EAST]->bType_sc == SOLID) && (curCell->neighCells[NORTH] == nullptr || curCell->neighCells[NORTH]->bType_sc == SOLID)) {
                cornerCells.push(curCell);
            } else if ((curCell->neighCells[WEST] == nullptr || curCell->neighCells[WEST]->bType_sc == SOLID) && (curCell->neighCells[NORTH] == nullptr || curCell->neighCells[NORTH]->bType_sc == SOLID)) {
                cornerCells.push(curCell);
            } else if ((curCell->neighCells[EAST] == nullptr || curCell->neighCells[EAST]->bType_sc == SOLID) && (curCell->neighCells[SOUTH] == nullptr || curCell->neighCells[SOUTH]->bType_sc == SOLID)) {
                cornerCells.push(curCell);
            } else if ((curCell->neighCells[WEST] == nullptr || curCell->neighCells[WEST]->bType_sc == SOLID) && (curCell->neighCells[SOUTH] == nullptr || curCell->neighCells[SOUTH]->bType_sc == SOLID)) {
                cornerCells.push(curCell);
            } else if (curCell->neighCells[EAST] == nullptr || curCell->neighCells[EAST]->bType_sc == SOLID) {
                curCell->p[CORRECTED_2] = curCell->neighCells[WEST]->p[CORRECTED_2] + curCell->g_p * curCell->faces[EAST]->dx;
            } else if (curCell->neighCells[WEST] == nullptr || curCell->neighCells[WEST]->bType_sc == SOLID) {
                curCell->p[CORRECTED_2] = curCell->neighCells[EAST]->p[CORRECTED_2] - curCell->g_p * curCell->faces[WEST]->dx;
            } else if (curCell->neighCells[NORTH] == nullptr || curCell->neighCells[NORTH]->bType_sc == SOLID) {
                curCell->p[CORRECTED_2] = curCell->neighCells[SOUTH]->p[CORRECTED_2] + curCell->g_p * curCell->faces[NORTH]->dy;
            } else if (curCell->neighCells[SOUTH] == nullptr || curCell->neighCells[SOUTH]->bType_sc == SOLID) {
                curCell->p[CORRECTED_2] = curCell->neighCells[NORTH]->p[CORRECTED_2] - curCell->g_p * curCell->faces[SOUTH]->dy;
            }
        }
    }
    // Update Corner Cells
    while (!cornerCells.empty()){
        Cell2D *curCell = cornerCells.top();
        cornerCells.pop();
        if ((curCell->neighCells[EAST] == nullptr || curCell->neighCells[EAST]->bType_sc == SOLID) && (curCell->neighCells[NORTH] == nullptr || curCell->neighCells[NORTH]->bType_sc == SOLID)) {
            curCell->p[CORRECTED_2] = 0.5 * (curCell->neighCells[WEST]->p[CORRECTED_2] + curCell->neighCells[SOUTH]->p[CORRECTED_2]);
        } else if ((curCell->neighCells[WEST] == nullptr || curCell->neighCells[WEST]->bType_sc == SOLID) && (curCell->neighCells[NORTH] == nullptr || curCell->neighCells[NORTH]->bType_sc == SOLID)) {
            curCell->p[CORRECTED_2] = 0.5 * (curCell->neighCells[EAST]->p[CORRECTED_2] + curCell->neighCells[SOUTH]->p[CORRECTED_2]);
        } else if ((curCell->neighCells[EAST] == nullptr || curCell->neighCells[EAST]->bType_sc == SOLID) && (curCell->neighCells[SOUTH] == nullptr || curCell->neighCells[SOUTH]->bType_sc == SOLID)) {
            curCell->p[CORRECTED_2] = 0.5 * (curCell->neighCells[WEST]->p[CORRECTED_2] + curCell->neighCells[NORTH]->p[CORRECTED_2]);
        } else if ((curCell->neighCells[WEST] == nullptr || curCell->neighCells[WEST]->bType_sc == SOLID) && (curCell->neighCells[SOUTH] == nullptr || curCell->neighCells[SOUTH]->bType_sc == SOLID)) {
            curCell->p[CORRECTED_2] = 0.5 * (curCell->neighCells[EAST]->p[CORRECTED_2] + curCell->neighCells[NORTH]->p[CORRECTED_2]);
        }
    }
    // Update Velocities bzw. faces
    for (int i = 0; i < data.nFaces; i++) {
        Face2D *curFace = &data.faces[i];
        if (curFace->bType_u == INNERCELL) {
            if (i < data.nhorizontalFaces) {
                double sum_a_v_nb = curFace->neighCells[NORTH]->neighCells[WEST]->faces[SOUTH]->a_p_tilde * curFace->neighCells[NORTH]->neighCells[WEST]->faces[SOUTH]->v[CORRECTED_1]
                            + curFace->neighCells[NORTH]->neighCells[EAST]->faces[SOUTH]->a_p_tilde * curFace->neighCells[NORTH]->neighCells[EAST]->faces[SOUTH]->v[CORRECTED_1]
                            + curFace->neighCells[NORTH]->faces[NORTH]->a_p_tilde * curFace->neighCells[NORTH]->faces[NORTH]->v[CORRECTED_1]
                            + curFace->neighCells[NORTH]->faces[SOUTH]->a_p_tilde * curFace->neighCells[NORTH]->faces[SOUTH]->v[CORRECTED_1];
                curFace->v[CORRECTED_2] = curFace->v[INTERMEDIATE_2] + ((curFace->neighCells[UP]->p[INTERMEDIATE_2] - curFace->neighCells[DOWN]->p[INTERMEDIATE_2]) * curFace->dx + sum_a_v_nb) / curFace->a_p_tilde;
            } else if (i >= data.nhorizontalFaces) {
                double sum_a_u_nb = curFace->neighCells[LEFT]->faces[WEST]->a_p_tilde * curFace->neighCells[LEFT]->faces[WEST]->u[CORRECTED_1]
                            + curFace->neighCells[RIGHT]->faces[EAST]->a_p_tilde * curFace->neighCells[RIGHT]->faces[EAST]->u[CORRECTED_1]
                            + curFace->neighCells[LEFT]->neighCells[NORTH]->faces[EAST]->a_p_tilde * curFace->neighCells[LEFT]->neighCells[NORTH]->faces[EAST]->u[CORRECTED_1]
                            + curFace->neighCells[RIGHT]->neighCells[SOUTH]->faces[WEST]->a_p_tilde * curFace->neighCells[RIGHT]->neighCells[SOUTH]->faces[WEST]->u[CORRECTED_1];
                curFace->u[CORRECTED_2] = curFace->u[INTERMEDIATE_2] + ((curFace->neighCells[LEFT]->p[INTERMEDIATE_2] - curFace->neighCells[RIGHT]->p[INTERMEDIATE_2]) * curFace->dy + sum_a_u_nb) / curFace->a_p_tilde;
            }
        } else if (curFace->bType_u == DIRICHLET || curFace->bType_u == SOLID){
            if (i < data.nhorizontalFaces) {
                curFace->v[CORRECTED_2] = curFace->v[INITIAL];
            } else if (i >= data.nhorizontalFaces) {
                curFace->u[CORRECTED_2] = curFace->u[INITIAL];
            }
        }
    }
        // Update Boundary Conditions
    for (int i = 0; i < data.nFaces; i++) {
        Face2D *curFace = &data.faces[i];
        if (curFace->bType_u == DIRICHLET || curFace->bType_u == SOLID){
            if (i < data.nhorizontalFaces) {
                // std::cout << "v :" << curFace->v[step] << ", type: " << curFace->bType_u <<std::endl;
                curFace->v[CORRECTED_2] = curFace->v[INITIAL];
            } else if (i >= data.nhorizontalFaces) {
                // std::cout << "u :" << curFace->u[step] << ", type: " << curFace->bType_u <<std::endl;
                curFace->u[CORRECTED_2] = curFace->u[INITIAL];
            }
        } else if (curFace->bType_u == NEUMANN) {
            if (i < data.nhorizontalFaces) {
                if(curFace->neighCells[UP] == nullptr || curFace->neighCells[UP]->bType_sc == SOLID)                                                     //TOP BOUNDARY (HORIZONTAL)
                {
                    curFace->v[CORRECTED_2] = curFace->neighCells[DOWN]->faces[SOUTH]->v[CORRECTED_2] + curFace->g_v * curFace->dx;
                }
                else if (curFace->neighCells[DOWN] == nullptr || curFace->neighCells[DOWN]->bType_sc == SOLID)                                           //BOTTOM BOUNDARY (HORIZONTAL) 
                {
                    curFace->v[CORRECTED_2] = curFace->neighCells[UP]->faces[NORTH]->v[CORRECTED_2] - curFace->g_v * curFace->dx;
                }
            } else {
                if (curFace->neighCells[LEFT] == nullptr || curFace->neighCells[LEFT]->bType_sc == SOLID)                                                 //LEFT BOUNDARY (VERTICAL) 
                {
                    curFace->u[CORRECTED_2] = curFace->neighCells[RIGHT]->faces[EAST]->u[CORRECTED_2] - curFace->g_u * curFace->dy;
                }
                else if (curFace->neighCells[RIGHT] == nullptr || curFace->neighCells[RIGHT]->bType_sc == SOLID)                                         //RIGHT BOUNDARY (VERTICAL)
                {
                    curFace->u[CORRECTED_2] = curFace->neighCells[LEFT]->faces[WEST]->u[CORRECTED_2] + curFace->g_u * curFace->dy;
                }
            }
        }
    }  
}

// Function to recompute the pressure in the small domain around the central node with the boundaries of the domain as dirichlet boundary 
// with p as the value in the previous calculation
void recomputePressureSmallDomain(Data2D& data, int step){
    // Create a vector with the cell ids of the small domain 1/10 (rounded) of the total domain in cells
    int centralCellId = std::round(data.nCells / 2);
    int nCells = data.nCells;
    int dimCellsX = data.dimX - 1;
    int dimCellsY = data.dimY - 1;
    int dimCellsXSmallDomain = std::round(dimCellsX);
    int dimCellsYSmallDomain = std::round(dimCellsY);
    int nCellsSmallDomain = dimCellsXSmallDomain * dimCellsYSmallDomain;
    //int startCellId = centralCellId - std::round(dimCellsXSmallDomain/2) - std::round(dimCellsYSmallDomain/2) * dimCellsX;
    //int endCellId = startCellId + dimCellsXSmallDomain - 1 + (dimCellsYSmallDomain - 1) * dimCellsX;
    int startCellId = 0;
    int endCellId = nCells - 1;
    double xEnd = data.cells[endCellId].points[1]->x;
    std::vector<int> boundaries;
    for (int i = startCellId; i < startCellId + dimCellsXSmallDomain; i++){
        boundaries.push_back(i);
    }
    int i = data.cells[startCellId].neighCells[NORTH]->id;
    while (i < endCellId - dimCellsXSmallDomain + 1){
        boundaries.push_back(i);
        i = data.cells[i].neighCells[NORTH]->id;
    }   
    i = data.cells[startCellId+dimCellsXSmallDomain-1].neighCells[NORTH]->id;
    while (i < endCellId){
        boundaries.push_back(i);
        i = data.cells[i].neighCells[NORTH]->id;
    }
    for (int i = endCellId - dimCellsXSmallDomain + 1; i <= endCellId; i++){
        boundaries.push_back(i);
    }

    // set up a matrix and vector for the pressure correction
    SpMat pressureCorrMatrix(nCellsSmallDomain, nCellsSmallDomain);
    Vector pressureCorrVector(nCellsSmallDomain);
    pressureCorrMatrix.setZero();
    pressureCorrVector.setZero();
    i = startCellId;
    std::vector<int> cellIds;
    bool isBoundary = false;
    int j = 0;
    while (i <= endCellId){
        Cell2D *curCell = &data.cells[i];
        for (int i = 0; i<boundaries.size(); i++){
            if (boundaries[i] == curCell->id){
                pressureCorrMatrix.coeffRef(j, j) = 1;
                pressureCorrVector(j) = curCell->p[step];
                isBoundary = true;
            }
        }
        if (!isBoundary){
            pressureCorrMatrix.coeffRef(j,j) = -curCell->a_p;
            
            pressureCorrMatrix.coeffRef(j, j + 1) = curCell->a_e;
            
            pressureCorrMatrix.coeffRef(j, j - 1) = curCell->a_w;
            
            pressureCorrMatrix.coeffRef(j, j - dimCellsXSmallDomain) = curCell->a_n;
            
            pressureCorrMatrix.coeffRef(j, j + dimCellsXSmallDomain) = curCell->a_s;

            pressureCorrVector(i) = -curCell->b;
        }
        cellIds.push_back(i);
        
        isBoundary = false;
        // i = curCell->neighCells[EAST]->id;
        // if (curCell->neighCells[EAST]->points[0]->x == xEnd){
        //     i = curCell->neighCells[NORTH]->id;
        //     i -= dimCellsXSmallDomain - 1;
        // }
        if (curCell->neighCells[EAST] != nullptr){
            i = curCell->neighCells[EAST]->id;
        } else if (curCell->neighCells[NORTH] != nullptr){
            i = curCell->neighCells[NORTH]->id;
            i -= (dimCellsXSmallDomain - 1);
        } else {
            i = endCellId + 1;
        }
        j++;
    }


    // std::cout << "Pressure Correction Matrix: " << std::endl;
    // std::cout << pressureCorrMatrix << std::endl;

    pressureCorrMatrix.makeCompressed();

    try {
        // Solve for pressure correction
        Eigen::SparseLU<SpMat> solver;
        solver.analyzePattern(pressureCorrMatrix);
        solver.factorize(pressureCorrMatrix);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Factorization failed for pressureCorrMatrix");
        }
        Vector pressureCorrSolution = solver.solve(pressureCorrVector);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Solving failed for pressureCorrMatrix");
        }
        for (int j = 0; j < pressureCorrSolution.size(); j++) {
            int id = cellIds[j];
            data.cells[id].p[step] = pressureCorrSolution(j);
        }
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}