#include "solve.h"
#include "data.h"
#include "corrector.h"

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
            pressureVector(i) = curCell.b;
            pressureMatrix.coeffRef(i, i) = curCell.a_p;
            if (curCell.neighCells[EAST]->bType_sc == INNERCELL) {
                pressureMatrix.coeffRef(i, curCell.neighCells[EAST]->id) = -curCell.a_e;
            }
            if (curCell.neighCells[WEST]->bType_sc == INNERCELL) {
                pressureMatrix.coeffRef(i, curCell.neighCells[WEST]->id) = -curCell.a_w;
            }
            if (curCell.neighCells[NORTH]->bType_sc == INNERCELL) {
                pressureMatrix.coeffRef(i, curCell.neighCells[NORTH]->id) = -curCell.a_n;
            }
            if (curCell.neighCells[SOUTH]->bType_sc == INNERCELL) {
                pressureMatrix.coeffRef(i, curCell.neighCells[SOUTH]->id) = -curCell.a_s;
            }
        } else if (curCell.bType_p == DIRICHLET || curCell.bType_p == SOLID) {
            pressureMatrix.coeffRef(i, i) = 1.0;
            pressureVector(i) = curCell.p[step - 1];
        } else if (curCell.bType_p == NEUMANN) {
            pressureMatrix.coeffRef(i, i) = 1.0;
            if (curCell.neighCells[EAST] == nullptr || curCell.neighCells[EAST]->bType_sc == SOLID) {
                pressureMatrix.coeffRef(i, curCell.neighCells[WEST]->id) = -1.0;
            }
            if (curCell.neighCells[WEST] == nullptr || curCell.neighCells[WEST]->bType_sc == SOLID) {
                pressureMatrix.coeffRef(i, curCell.neighCells[EAST]->id) = -1.0;
            }
            if (curCell.neighCells[NORTH] == nullptr || curCell.neighCells[NORTH]->bType_sc == SOLID) {
                pressureMatrix.coeffRef(i, curCell.neighCells[SOUTH]->id) = -1.0;
            }
            if (curCell.neighCells[SOUTH] == nullptr || curCell.neighCells[SOUTH]->bType_sc == SOLID) {
                pressureMatrix.coeffRef(i, curCell.neighCells[NORTH]->id) = -1.0;
            }
        }
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
    curCell->a_p = curCell->a_e + curCell->a_w + curCell->a_n + curCell->a_s;
    curCell->b = - ((fRho(data, curCell->alpha) + fRho(data, curCell->neighCells[EAST]->alpha))* 0.5 *curCell->faces[EAST]->u[INTERMEDIATE_1] * dy) + ((fRho(data, curCell->alpha) + fRho(data, curCell->neighCells[WEST]->alpha))* 0.5 *curCell->faces[WEST]->u[INTERMEDIATE_1] * dy) 
    - ((fRho(data, curCell->alpha) + fRho(data, curCell->neighCells[NORTH]->alpha))* 0.5 *curCell->faces[NORTH]->v[INTERMEDIATE_1] * dx) + ((fRho(data, curCell->alpha) + fRho(data, curCell->neighCells[SOUTH]->alpha))* 0.5 *curCell->faces[SOUTH]->v[INTERMEDIATE_1] * dx);
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
    pressureCorrMatrix.makeCompressed();

    try {
        // Solve for pressure correction
        Eigen::BiCGSTAB<SpMat> solver;
        solver.setTolerance(1e-6);
        solver.setMaxIterations(10000);
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
    data.alpha_p_relax = 1.0;
    for (int i = 0; i < data.nCells; i++) {
        Cell2D *curCell = &data.cells[i];
        if (curCell->bType_p == INNERCELL) {
            curCell->p[CORRECTED_1] = curCell->p[INITIAL] + curCell->p[INTERMEDIATE_1] * data.alpha_p_relax;
            // std::cout << "Cell " << i << "; p: " << curCell->p[step + 1] << std::endl;
        }
    }
    // Update Boundary Conditions
    for (int i = 0; i < data.nCells; i++){
        Cell2D *curCell = &data.cells[i];
        if (curCell->bType_p == DIRICHLET || curCell->bType_p == SOLID) {
            curCell->p[CORRECTED_1] = curCell->p[INITIAL];
        } else if (curCell->bType_p == NEUMANN){
            if (curCell->neighCells[EAST] == nullptr || curCell->neighCells[EAST]->bType_sc == SOLID) {
                curCell->p[CORRECTED_1] = curCell->neighCells[WEST]->p[CORRECTED_1];
            }
            if (curCell->neighCells[WEST] == nullptr || curCell->neighCells[WEST]->bType_sc == SOLID) {
                curCell->p[CORRECTED_1] = curCell->neighCells[EAST]->p[CORRECTED_1];
            }
            if (curCell->neighCells[NORTH] == nullptr || curCell->neighCells[NORTH]->bType_sc == SOLID) {
                curCell->p[CORRECTED_1] = curCell->neighCells[SOUTH]->p[CORRECTED_1];
            }
            if (curCell->neighCells[SOUTH] == nullptr || curCell->neighCells[SOUTH]->bType_sc == SOLID) {
                curCell->p[CORRECTED_1] = curCell->neighCells[NORTH]->p[CORRECTED_1];
            }
        }
    }
    // Update Velocities bzw. faces
    for (int i = 0; i < data.nFaces; i++) {
        Face2D *curFace = &data.faces[i];
        if (curFace->bType_u == INNERCELL) {
            if (i < data.nhorizontalFaces) {
                curFace->v[CORRECTED_1] = curFace->v[INTERMEDIATE_1] + ((curFace->neighCells[UP]->p[INTERMEDIATE_1] - curFace->neighCells[DOWN]->p[INTERMEDIATE_1]) * curFace->dx) / curFace->a_p_tilde;
                // std::cout << "Face " << i << "; v: " << curFace->v[step + 1] << std::endl;
            } else if (i >= data.nhorizontalFaces) {
                curFace->u[CORRECTED_1] = curFace->u[INTERMEDIATE_1] + ((curFace->neighCells[LEFT]->p[INTERMEDIATE_1] - curFace->neighCells[RIGHT]->p[INTERMEDIATE_1]) * curFace->dy) / curFace->a_p_tilde;
                // std::cout << "Face " << i << "; u: " << curFace->u[step + 1] << std::endl;
            }   
        } 
    }
    // Update Boundary Conditions
    for (int i = 0; i < data.nFaces; i++) {
        Face2D *curFace = &data.faces[i];
        if (curFace->bType_u == DIRICHLET || curFace->bType_u == SOLID){
            if (i < data.nhorizontalFaces) {
                // std::cout << "v :" << curFace->v[step] << ", type: " << curFace->bType_u <<std::endl;
                curFace->v[CORRECTED_1] = curFace->v[INITIAL];
            } else if (i >= data.nhorizontalFaces) {
                // std::cout << "u :" << curFace->u[step] << ", type: " << curFace->bType_u <<std::endl;
                curFace->u[CORRECTED_1] = curFace->u[INITIAL];
            }
        } else if (curFace->bType_u == NEUMANN) {
            if (i < data.nhorizontalFaces) {
                if(curFace->neighCells[UP] == nullptr || curFace->neighCells[UP]->bType_sc == SOLID)                                                     //TOP BOUNDARY (HORIZONTAL)
                {
                    curFace->v[CORRECTED_1] = curFace->neighCells[DOWN]->faces[SOUTH]->v[CORRECTED_1];
                }
                else if (curFace->neighCells[DOWN] == nullptr || curFace->neighCells[DOWN]->bType_sc == SOLID)                                           //BOTTOM BOUNDARY (HORIZONTAL) 
                {
                    curFace->v[CORRECTED_1] = curFace->neighCells[UP]->faces[NORTH]->v[CORRECTED_1];
                }
            } else {
                if (curFace->neighCells[LEFT] == nullptr || curFace->neighCells[LEFT]->bType_sc == SOLID)                                                 //LEFT BOUNDARY (VERTICAL) 
                {
                    curFace->u[CORRECTED_1] = curFace->neighCells[RIGHT]->faces[EAST]->u[CORRECTED_1];
                }
                else if (curFace->neighCells[RIGHT] == nullptr || curFace->neighCells[RIGHT]->bType_sc == SOLID)                                         //RIGHT BOUNDARY (VERTICAL)
                {
                    curFace->u[CORRECTED_1] = curFace->neighCells[LEFT]->faces[WEST]->u[CORRECTED_1];
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
                * ((curCell->neighCells[EAST]->faces[EAST]->u[CORRECTED_1] - curCell->neighCells[EAST]->faces[EAST]->u[INTERMEDIATE_1]) * curCell->neighCells[EAST]->faces[EAST]->a_p_tilde
                + (curCell->faces[WEST]->u[CORRECTED_1] - curCell->faces[WEST]->u[INTERMEDIATE_1]) * curCell->faces[WEST]->a_p_tilde
                + (curCell->neighCells[NORTH]->faces[EAST]->u[CORRECTED_1] - curCell->neighCells[NORTH]->faces[EAST]->u[INTERMEDIATE_1]) * curCell->neighCells[NORTH]->faces[EAST]->a_p_tilde
                + (curCell->neighCells[SOUTH]->faces[EAST]->u[CORRECTED_1] - curCell->neighCells[SOUTH]->faces[EAST]->u[INTERMEDIATE_1]) * curCell->neighCells[SOUTH]->faces[EAST]->a_p_tilde); 

    double b_w = (0.5 * (fRho(data, curCell->neighCells[WEST]->alpha) + fRho(data, curCell->alpha)) * dy) / curCell->faces[WEST]->a_p_tilde
                * ((curCell->faces[EAST]->u[CORRECTED_1] - curCell->faces[EAST]->u[INTERMEDIATE_1]) * curCell->faces[EAST]->a_p_tilde
                + (curCell->neighCells[WEST]->faces[WEST]->u[CORRECTED_1] - curCell->neighCells[WEST]->faces[WEST]->u[INTERMEDIATE_1]) * curCell->neighCells[WEST]->faces[WEST]->a_p_tilde
                + (curCell->neighCells[NORTH]->faces[WEST]->u[CORRECTED_1] - curCell->neighCells[NORTH]->faces[WEST]->u[INTERMEDIATE_1]) * curCell->neighCells[NORTH]->faces[WEST]->a_p_tilde
                + (curCell->neighCells[SOUTH]->faces[WEST]->u[CORRECTED_1] - curCell->neighCells[SOUTH]->faces[WEST]->u[INTERMEDIATE_1]) * curCell->neighCells[SOUTH]->faces[WEST]->a_p_tilde);

    double b_n = (0.5 * (fRho(data, curCell->neighCells[NORTH]->alpha) + fRho(data, curCell->alpha)) * dx) / curCell->faces[NORTH]->a_p_tilde
                * ((curCell->neighCells[NORTH]->faces[NORTH]->v[CORRECTED_1] - curCell->neighCells[NORTH]->faces[NORTH]->v[INTERMEDIATE_1]) * curCell->neighCells[NORTH]->faces[NORTH]->a_p_tilde
                + (curCell->neighCells[WEST]->faces[NORTH]->v[CORRECTED_1] - curCell->neighCells[WEST]->faces[NORTH]->v[INTERMEDIATE_1]) * curCell->neighCells[WEST]->faces[NORTH]->a_p_tilde
                + (curCell->neighCells[EAST]->faces[NORTH]->v[CORRECTED_1] - curCell->neighCells[EAST]->faces[NORTH]->v[INTERMEDIATE_1]) * curCell->neighCells[EAST]->faces[NORTH]->a_p_tilde
                + (curCell->faces[SOUTH]->v[CORRECTED_1] - curCell->faces[SOUTH]->v[INTERMEDIATE_1]) * curCell->faces[SOUTH]->a_p_tilde);

    double b_s = (0.5 * (fRho(data, curCell->neighCells[SOUTH]->alpha) + fRho(data, curCell->alpha)) * dx) / curCell->faces[SOUTH]->a_p_tilde
                * ((curCell->faces[NORTH]->v[CORRECTED_1] - curCell->faces[NORTH]->v[INTERMEDIATE_1]) * curCell->faces[NORTH]->a_p_tilde
                + (curCell->neighCells[WEST]->faces[SOUTH]->v[CORRECTED_1] - curCell->neighCells[WEST]->faces[SOUTH]->v[INTERMEDIATE_1]) * curCell->neighCells[WEST]->faces[SOUTH]->a_p_tilde
                + (curCell->neighCells[EAST]->faces[SOUTH]->v[CORRECTED_1] - curCell->neighCells[EAST]->faces[SOUTH]->v[INTERMEDIATE_1]) * curCell->neighCells[EAST]->faces[SOUTH]->a_p_tilde
                + (curCell->faces[SOUTH]->v[CORRECTED_1] - curCell->faces[SOUTH]->v[INTERMEDIATE_1]) * curCell->faces[SOUTH]->a_p_tilde);
    
    curCell->b = b_e + b_w + b_n + b_s;
    if (data.mode == 0){
        curCell->b += (fRho(data, curCell->alpha_prev) - fRho(data, curCell->alpha)) * dx * dy / data.dt;
    }

    curCell->a_e = (fRho(data, curCell->neighCells[EAST]->alpha) * dy * dy) / curCell->faces[EAST]->a_p_tilde;
    curCell->a_w = (fRho(data, curCell->neighCells[WEST]->alpha) * dy * dy) / curCell->faces[WEST]->a_p_tilde;
    curCell->a_n = (fRho(data, curCell->neighCells[NORTH]->alpha) * dx * dx) / curCell->faces[NORTH]->a_p_tilde;
    curCell->a_s = (fRho(data, curCell->neighCells[SOUTH]->alpha) * dx * dx) / curCell->faces[SOUTH]->a_p_tilde;
}

void corrector2(Data2D& data){
    data.alpha_p_relax = 1.0;
    for (int i = 0; i < data.nCells; i++) {
        Cell2D *curCell = &data.cells[i];
        if (curCell->bType_p == INNERCELL) {
            curCell->p[CORRECTED_2] = curCell->p[CORRECTED_1] + curCell->p[INTERMEDIATE_2] * data.alpha_p_relax;
        }
    }
    for (int i = 0; i < data.nFaces; i++) {
        Face2D *curFace = &data.faces[i];
        if (curFace->bType_u == INNERCELL) {
            if (i < data.nhorizontalFaces) {
                double dv_w = curFace->neighCells[NORTH]->neighCells[WEST]->faces[SOUTH]->v[CORRECTED_1] - curFace->neighCells[NORTH]->neighCells[WEST]->faces[SOUTH]->v[INTERMEDIATE_1];
                double dv_e = curFace->neighCells[NORTH]->neighCells[EAST]->faces[SOUTH]->v[CORRECTED_1] - curFace->neighCells[NORTH]->neighCells[EAST]->faces[SOUTH]->v[INTERMEDIATE_1];
                double dv_n = curFace->neighCells[NORTH]->faces[NORTH]->v[CORRECTED_1] - curFace->neighCells[NORTH]->faces[NORTH]->v[INTERMEDIATE_1];
                double dv_s = curFace->neighCells[NORTH]->faces[SOUTH]->v[CORRECTED_1] - curFace->neighCells[NORTH]->faces[SOUTH]->v[INTERMEDIATE_1];
                double a_dv_nb = curFace->neighCells[NORTH]->neighCells[WEST]->faces[SOUTH]->a_p_tilde * dv_w
                            + curFace->neighCells[NORTH]->neighCells[EAST]->faces[SOUTH]->a_p_tilde * dv_e
                            + curFace->neighCells[NORTH]->faces[NORTH]->a_p_tilde * dv_n
                            + curFace->neighCells[NORTH]->faces[SOUTH]->a_p_tilde * dv_s;
                curFace->v[CORRECTED_2] = curFace->v[CORRECTED_1] + ((curFace->neighCells[UP]->p[INTERMEDIATE_2] - curFace->neighCells[DOWN]->p[INTERMEDIATE_2]) * curFace->dx + a_dv_nb) / curFace->a_p_tilde;
            } else if (i >= data.nhorizontalFaces) {
                double du_w = curFace->neighCells[LEFT]->faces[WEST]->u[CORRECTED_1] - curFace->neighCells[LEFT]->faces[WEST]->u[INTERMEDIATE_1];
                double du_e = curFace->neighCells[RIGHT]->faces[EAST]->u[CORRECTED_1] - curFace->neighCells[RIGHT]->faces[EAST]->u[INTERMEDIATE_1];
                double du_n = curFace->neighCells[LEFT]->neighCells[NORTH]->faces[EAST]->u[CORRECTED_1] - curFace->neighCells[LEFT]->neighCells[NORTH]->faces[EAST]->u[INTERMEDIATE_1];
                double du_s = curFace->neighCells[RIGHT]->neighCells[SOUTH]->faces[WEST]->u[CORRECTED_1] - curFace->neighCells[RIGHT]->neighCells[SOUTH]->faces[WEST]->u[INTERMEDIATE_1];
                double a_du_nb = curFace->neighCells[LEFT]->faces[WEST]->a_p_tilde * du_w
                            + curFace->neighCells[RIGHT]->faces[EAST]->a_p_tilde * du_e
                            + curFace->neighCells[LEFT]->neighCells[NORTH]->faces[EAST]->a_p_tilde * du_n
                            + curFace->neighCells[RIGHT]->neighCells[SOUTH]->faces[WEST]->a_p_tilde * du_s;
                curFace->u[CORRECTED_2] = curFace->u[CORRECTED_1] + ((curFace->neighCells[LEFT]->p[INTERMEDIATE_2] - curFace->neighCells[RIGHT]->p[INTERMEDIATE_2]) * curFace->dy + a_du_nb) / curFace->a_p_tilde;
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
                    curFace->v[CORRECTED_2] = curFace->neighCells[DOWN]->faces[SOUTH]->v[CORRECTED_2];
                }
                else if (curFace->neighCells[DOWN] == nullptr || curFace->neighCells[DOWN]->bType_sc == SOLID)                                           //BOTTOM BOUNDARY (HORIZONTAL) 
                {
                    curFace->v[CORRECTED_2] = curFace->neighCells[UP]->faces[NORTH]->v[CORRECTED_2];
                }
            } else {
                if (curFace->neighCells[LEFT] == nullptr || curFace->neighCells[LEFT]->bType_sc == SOLID)                                                 //LEFT BOUNDARY (VERTICAL) 
                {
                    curFace->u[CORRECTED_2] = curFace->neighCells[RIGHT]->faces[EAST]->u[CORRECTED_2];
                }
                else if (curFace->neighCells[RIGHT] == nullptr || curFace->neighCells[RIGHT]->bType_sc == SOLID)                                         //RIGHT BOUNDARY (VERTICAL)
                {
                    curFace->u[CORRECTED_2] = curFace->neighCells[LEFT]->faces[WEST]->u[CORRECTED_2];
                }
            }
        }
    }  
}