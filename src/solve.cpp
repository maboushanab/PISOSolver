#include "data.h"
#include "solve.h"
#include "output.h"
#include <iostream>

std::string fSolve(Data2D& data) {
    std::string directoryName = createDirectory();
    std::cout << "Solve" << std::endl;
    data.timeStep = 0;
    unstaggerGrid(data, INITIAL);
    fOutputVTKframe(data, directoryName, INITIAL);
    // Time loop
    double time = 0;
    while (time < data.maxTime) {
        std::cout << "Time: " << data.timeStep << std::endl;
        data.timeStep ++;

        // PISO algorithm
        if (data.mode == 0){            // Transient
            iterateTransient(data);
        } else {                        // Steady State
            iterateSteady(data, 1);
        }


        unstaggerGrid(data, CORRECTED_2);
        fOutputVTKframe(data, directoryName, CORRECTED_2);
        resetData(data);
        if (data.mode == 0){
            time += data.dt;
        } else if (data.mode == 1){
            return directoryName;
        }
    }
    return directoryName;
}

void iterateSteady(Data2D& data, int iteration){
    predictVelocityField(data);

    correctPressureEquation(data, INTERMEDIATE_1);
    corrector1(data);

    correctPressureEquation(data, INTERMEDIATE_2);
    corrector2(data);

    calcScalarTransfer(data);
    checkConvergence(data, iteration);
}

void iterateTransient(Data2D& data){
    assignPrevData(data);
    predictVelocityField(data);

    correctPressureEquation(data, INTERMEDIATE_1);
    corrector1(data);

    correctPressureEquation(data, INTERMEDIATE_2);
    corrector2(data);

    calcScalarTransfer(data);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                         SUPPLEMANTRY FUNCTIONS                                                 //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////                      

/**
 * Calculates the value of A based on the Peclet number.
 *
 * @param Pe The Peclet number.
 * @return The calculated value of A.
 */
double A(Data2D data, double Pe){
    double res = 0.0;
    if (data.pecFunc == 0) res = 1.0;                                                       // Upwind
    else if (data.pecFunc == 1) res = std::max(0.0, std::pow(1.0 - 0.1 * std::abs(Pe), 5)); // Potenzgesetz
    else if (data.pecFunc == 2) res = 1 - 0.5 * std::abs(Pe);                               // Central
    else if (data.pecFunc == 3) res = std::abs(Pe)/(std::exp(std::abs(Pe)) - 1);            // Exponential
    else if (data.pecFunc == 4) res = std::max(0.0, 1.0 - 0.5 * std::abs(Pe));              // Hybrid
    return res;
}

/**
 * Calculates the interpolated density value based on the given data and alpha value.
 * 
 * @param data The data structure containing the rho values.
 * @param alpha The phase fraction of the cell
 * @return The weighted density of the cell.
 */
double fRho(Data2D& data, double alpha){
    return (alpha * data.rho1 + (1 - alpha) * data.rho2);
}

/**
 * Calculates the interpolated viscosity value based on the given data and alpha value.
 *
 * @param data The data structure containing the eta values.
 * @param alpha The phase fraction of the cell
 * @return The weighted average of the eta values.
 */
double fEta(Data2D& data, double alpha){
    return (alpha * data.eta1 + (1 - alpha) * data.eta2);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                PISO ALGORITHM                                                  //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                  PREDICTOR                                                     //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





/**
 * Computes the velocity coefficients in the x-direction for a given face to solve the momentum equation.
 * 
 * @param data The data structure containing the simulation data.
 * @param faceId The ID of the face for which to compute the velocity coefficients.
 * @param step The current step of the PISO algorithm.
 */
void computeVelocityCoeff_x(Data2D& data, int faceId) {
    Face2D *curFace = &data.faces[faceId];

    // Velocity Terms (U)
    double u = curFace->u[INITIAL];
    double u_W = curFace->neighCells[LEFT]->faces[WEST]->u[INITIAL];
    double u_E = curFace->neighCells[RIGHT]->faces[EAST]->u[INITIAL];

    // Velocity Terms (V)
    double v_S_L = curFace->neighCells[LEFT]->faces[SOUTH]->v[INITIAL];
    double v_S_R = curFace->neighCells[RIGHT]->faces[SOUTH]->v[INITIAL];
    double v_N_L = curFace->neighCells[LEFT]->faces[NORTH]->v[INITIAL];
    double v_N_R = curFace->neighCells[RIGHT]->faces[NORTH]->v[INITIAL];

    // Cell Properties
    double dx = curFace->neighCells[LEFT]->faces[SOUTH]->dx;
    double dy = curFace->dy;
    double dt = data.dt;

    // Diffusion Terms
    double d_n = (fEta(data, curFace->neighCells[LEFT]->neighCells[NORTH]->alpha) + fEta(data, curFace->neighCells[LEFT]->alpha) + fEta(data, curFace->neighCells[RIGHT]->neighCells[NORTH]->alpha) + fEta(data, curFace->neighCells[RIGHT]->alpha) )/ (4 * dy);
    double d_s = (fEta(data, curFace->neighCells[LEFT]->neighCells[SOUTH]->alpha) + fEta(data, curFace->neighCells[LEFT]->alpha) + fEta(data, curFace->neighCells[RIGHT]->neighCells[SOUTH]->alpha) + fEta(data, curFace->neighCells[RIGHT]->alpha) )/ (4 * dy);
    double d_e = fEta(data, curFace->neighCells[RIGHT]->alpha) / dx;
    double d_w = fEta(data, curFace->neighCells[LEFT]->alpha) / dx;


    // Mass Flux Terms
    double f_e;
    double f_w;

    if (curFace->neighCells[RIGHT]->neighCells[EAST] == nullptr) {
        f_e = 0.5 * (fRho(data, curFace->neighCells[RIGHT]->alpha) * u_E + (fRho(data, curFace->neighCells[LEFT]->alpha) + fRho(data, curFace->neighCells[RIGHT]->alpha)) / 2 * u); //In case face is one away of the right boundary
    } else {
        f_e = 0.5 * ((fRho(data, curFace->neighCells[RIGHT]->alpha) + fRho(data, curFace->neighCells[RIGHT]->neighCells[EAST]->alpha)) / 2 * u_E + (fRho(data, curFace->neighCells[LEFT]->alpha) + fRho(data, curFace->neighCells[RIGHT]->alpha)) / 2 * u); 
    }

    if (curFace->neighCells[LEFT]->neighCells[WEST] == nullptr) {
        f_w = 0.5 * ((fRho(data, curFace->neighCells[LEFT]->alpha) + fRho(data, curFace->neighCells[RIGHT]->alpha)) / 2 * u + fRho(data, curFace->neighCells[LEFT]->alpha) / 2 * u_W); //In case face is one away of the left boundary
    } else {
        f_w = 0.5 * ((fRho(data, curFace->neighCells[LEFT]->alpha) + fRho(data, curFace->neighCells[RIGHT]->alpha)) / 2 * u + (fRho(data, curFace->neighCells[LEFT]->alpha) + fRho(data, curFace->neighCells[LEFT]->neighCells[WEST]->alpha)) / 2 * u_W);
    }

    double g_n = 0.5 * ((fRho(data, curFace->neighCells[RIGHT]->alpha) + fRho(data, curFace->neighCells[RIGHT]->neighCells[NORTH]->alpha)) / 2 * v_N_R + (fRho(data, curFace->neighCells[LEFT]->alpha) + fRho(data, curFace->neighCells[LEFT]->neighCells[NORTH]->alpha)) / 2 * v_N_L);
    double g_s = 0.5 * ((fRho(data, curFace->neighCells[RIGHT]->alpha) + fRho(data, curFace->neighCells[RIGHT]->neighCells[SOUTH]->alpha)) / 2 * v_S_R + (fRho(data, curFace->neighCells[LEFT]->alpha) + fRho(data, curFace->neighCells[LEFT]->neighCells[SOUTH]->alpha)) / 2 * v_S_L);

    // Peclet Number
    double Pe_e = f_e / d_e;
    double Pe_w = f_w / d_w;
    double Pe_n = g_n / d_n;
    double Pe_s = g_s / d_s;

    // Compute the coefficients
    curFace->a_e = d_e * dy * A(data, Pe_e) + std::max(-f_e, 0.0);
    curFace->a_w = d_w * dy * A(data, Pe_w) + std::max(f_w, 0.0);
    curFace->a_n = d_n * dx * A(data, Pe_n) + std::max(-g_n, 0.0);
    curFace->a_s = d_s * dx * A(data, Pe_s) + std::max(g_s, 0.0);
    curFace->a_p_v = 0.5 * (fRho(data, curFace->neighCells[LEFT]->alpha) + fRho(data, curFace->neighCells[RIGHT]->alpha)) * dx * dy / dt;
    curFace->a_p_tilde = curFace->a_e + curFace->a_w + curFace->a_n + curFace->a_s;
    curFace->b = (curFace->neighCells[LEFT]->p[INITIAL] - curFace->neighCells[RIGHT]->p[INITIAL]) * dy;
    if (data.mode == 0){
        curFace->b += curFace->a_p_v * curFace->u_prev;
    }
}

/**
 * Computes the velocity coefficients in the y-direction for a given face to solve the momentum equation.
 * 
 * @param data The data structure containing the simulation data.
 * @param faceId The ID of the face for which to compute the velocity coefficients.
 * @param step The current step of the PISO algorithm.
 */
void computeVelocityCoeff_y(Data2D& data, int faceId) {
    Face2D *curFace = &data.faces[faceId];

    // Velocity Terms (U)
    double u_E_U = curFace->neighCells[UP]->faces[EAST]->u[INITIAL];
    double u_E_D = curFace->neighCells[DOWN]->faces[EAST]->u[INITIAL];
    double u_W_U = curFace->neighCells[UP]->faces[WEST]->u[INITIAL];
    double u_W_D = curFace->neighCells[DOWN]->faces[WEST]->u[INITIAL];

    // Velocity Terms (V)
    double v = curFace->v[INITIAL];
    double v_S = curFace->neighCells[DOWN]->faces[SOUTH]->v[INITIAL];
    double v_N = curFace->neighCells[UP]->faces[NORTH]->v[INITIAL];

    // Cell Properties
    double dx = curFace->dx;
    double dy = curFace->neighCells[UP]->faces[WEST]->dy;
    double dt = data.dt;

    // Diffusion Terms
    double d_n = fEta(data, curFace->neighCells[UP]->alpha) / dy;
    double d_s = fEta(data, curFace->neighCells[DOWN]->alpha) / dy;
    double d_e = (fEta(data, curFace->neighCells[UP]->neighCells[EAST]->alpha) + fEta(data, curFace->neighCells[UP]->alpha) + fEta(data, curFace->neighCells[DOWN]->neighCells[EAST]->alpha) + fEta(data, curFace->neighCells[DOWN]->alpha) ) / (4 * dx);
    double d_w = (fEta(data, curFace->neighCells[UP]->neighCells[WEST]->alpha) + fEta(data, curFace->neighCells[UP]->alpha) + fEta(data, curFace->neighCells[DOWN]->neighCells[WEST]->alpha) + fEta(data, curFace->neighCells[DOWN]->alpha) ) / (4 * dx);

    // Mass Flux Terms
    double f_e = 0.5 * ((fRho(data, curFace->neighCells[UP]->alpha) + fRho(data, curFace->neighCells[UP]->neighCells[EAST]->alpha)) / 2 * u_E_U + (fRho(data, curFace->neighCells[DOWN]->alpha) + fRho(data, curFace->neighCells[DOWN]->neighCells[EAST]->alpha)) / 2 * u_E_D);
    double f_w = 0.5 * ((fRho(data, curFace->neighCells[UP]->alpha) + fRho(data, curFace->neighCells[UP]->neighCells[WEST]->alpha)) / 2 * u_W_U + (fRho(data, curFace->neighCells[DOWN]->alpha) + fRho(data, curFace->neighCells[DOWN]->neighCells[WEST]->alpha)) / 2 * u_W_D);

    double g_n;
    double g_s;

    if (curFace->neighCells[UP]->neighCells[NORTH] == nullptr) {
        g_n = 0.5 * (fRho(data, curFace->neighCells[UP]->alpha) / 2 * v_N + (fRho(data, curFace->neighCells[DOWN]->alpha) + fRho(data, curFace->neighCells[UP]->alpha)) / 2 * v); //In case face is one away of the top boundary
    } else {
        g_n = 0.5 * ((fRho(data, curFace->neighCells[UP]->alpha) + fRho(data, curFace->neighCells[UP]->neighCells[NORTH]->alpha)) / 2 * v_N + (fRho(data, curFace->neighCells[DOWN]->alpha) + fRho(data, curFace->neighCells[UP]->alpha)) / 2 * v);
    }

    if (curFace->neighCells[DOWN]->neighCells[SOUTH] == nullptr) {
        g_s = 0.5 * (fRho(data, curFace->neighCells[DOWN]->alpha) * v + (fRho(data, curFace->neighCells[UP]->alpha) + fRho(data, curFace->neighCells[DOWN]->alpha)) / 2 * v_S); //In case face is one away of the bottom boundary
    } else {
        g_s = 0.5 * ((fRho(data, curFace->neighCells[DOWN]->alpha) + fRho(data, curFace->neighCells[DOWN]->neighCells[SOUTH]->alpha)) / 2 * v + (fRho(data, curFace->neighCells[UP]->alpha) + fRho(data, curFace->neighCells[DOWN]->alpha)) / 2 * v_S);
    }

    // Peclet Number
    double Pe_e = f_e / d_e;
    double Pe_w = f_w / d_w;
    double Pe_n = g_n / d_n;
    double Pe_s = g_s / d_s;

    // Compute the coefficients
    curFace->a_e = d_e * dy * A(data, Pe_e) + std::max(-f_e, 0.0);
    curFace->a_w = d_w * dy * A(data, Pe_w) + std::max(f_w, 0.0);
    curFace->a_n = d_n * dx * A(data, Pe_n) + std::max(-g_n, 0.0);
    curFace->a_s = d_s * dx * A(data, Pe_s) + std::max(g_s, 0.0);
    curFace->a_p_v = 0.5 * (fRho(data, curFace->neighCells[UP]->alpha) + fRho(data, curFace->neighCells[DOWN]->alpha)) * dx * dy / dt;
    curFace->a_p_tilde = curFace->a_e + curFace->a_w + curFace->a_n + curFace->a_s;
    curFace->b = (curFace->neighCells[DOWN]->p[INITIAL] - curFace->neighCells[UP]->p[INITIAL]) * dx;
    if (data.mode == 0){
        curFace->b += curFace->a_p_v * curFace->v_prev;
    }
}

/**
 * Sets the momentum equation matrix and vector for a given data set.
 * Ax = b
 *
 * @param data The data structure containing the simulation data.
 * @param momMatrix The momentum equation matrix to be set. (A)
 * @param momVector The momentum equation vector to be set. (b)
 * @param step The current step of the PISO algorithm.
 */
void setMomentumEquationMatrix(Data2D& data, SpMat& momMatrix, Vector& momVector){
    for (int i = 0; i < data.nFaces; i++) {
        Face2D* curFace = &data.faces[i];
        if (curFace->bType_u == INNERCELL) {
            momMatrix.insert(i, i) = -curFace->a_p_tilde;
            momVector(i) = curFace->b;
            if (i > 0) {
                momMatrix.insert(i, i - 1) = curFace->a_w;
            }
            if (i < data.nFaces - 1) {
                momMatrix.insert(i, i + 1) = curFace->a_e;
            }
            if (i < data.nhorizontalFaces){
                if (curFace->neighCells[UP] != nullptr) {
                    int j = curFace->neighCells[UP]->faces[NORTH]->id;
                    momMatrix.insert(i, j) = curFace->a_n;
                }
                if (curFace->neighCells[DOWN] != nullptr) {
                    int j = curFace->neighCells[DOWN]->faces[SOUTH]->id;
                    momMatrix.insert(i, j) = curFace->a_s;
                }
            } else {
                if (curFace->neighCells[LEFT]->neighCells[NORTH] != nullptr) {
                    int j = curFace->neighCells[LEFT]->neighCells[NORTH]->faces[EAST]->id;
                    momMatrix.insert(i, j) = curFace->a_n;
                }
                if (curFace->neighCells[RIGHT]->neighCells[SOUTH] != nullptr) {
                    int j = curFace->neighCells[RIGHT]->neighCells[SOUTH]->faces[WEST]->id;
                    momMatrix.insert(i, j) = curFace->a_s;
                }
            }
        } else if (curFace->bType_u == DIRICHLET || curFace->bType_u == SOLID) {
            momMatrix.insert(i, i) = 1.0;
            if (i < data.nhorizontalFaces) {
                momVector(i) = curFace->v[INITIAL];
            } else {
                momVector(i) = curFace->u[INITIAL];
            }
        } else if (curFace->bType_u == NEUMANN) {
            momMatrix.insert(i, i) = -1.0;
            if (i < data.nhorizontalFaces) {
                if (curFace->neighCells[UP] != nullptr) {
                    int j = curFace->neighCells[UP]->faces[NORTH]->id;
                    momMatrix.insert(i, j) = 1.0;
                }
                if (curFace->neighCells[DOWN] != nullptr) {
                    int j = curFace->neighCells[DOWN]->faces[SOUTH]->id;
                    momMatrix.insert(i, j) = 1.0;
                }
            } else {
                if (curFace->neighCells[LEFT] != nullptr) {
                    int j = curFace->neighCells[LEFT]->faces[WEST]->id;
                    momMatrix.insert(i, j) = 1.0;
                }
                if (curFace->neighCells[RIGHT] != nullptr) {
                    int j = curFace->neighCells[RIGHT]->faces[EAST]->id;
                    momMatrix.insert(i, j) = 1.0;
                }
            }
        }

    }
}

/**
 * Predicts the velocity field for the given data.
 *
 * This function iterates over the faces in the data and computes the velocity coefficients
 * based on the face type. It then sets up the momentum equation matrices and vectors for
 * the u and v components of the velocity field. The matrices are solved using the BiCGSTAB
 * solver, and the resulting solutions are stored in the data structure.
 *
 * @param data The data structure containing the face and velocity information.
 */
///////////////////////////////////////////////
//  IMPLICIT GAUSS SEIDEL SOLVER METHOD      //
///////////////////////////////////////////////
void gaussSeidelSolver(const SpMat& A, Vector& x, const Vector& b, int maxIterations, double tolerance) {
    int n = A.rows();
    Vector x_old = x;
    double sum;
    
    for (int k = 0; k < maxIterations; ++k) {
        for (int i = 0; i < n; ++i) {
            sum = 0.0;
            for (Eigen::SparseMatrix<double>::InnerIterator it(A, i); it; ++it) {
                if (it.index() != i) {
                    sum += it.value() * x[it.index()];
                }
            }
            x[i] = (b[i] - sum) / A.coeff(i, i);
        }
        if ((x - x_old).norm() < tolerance) {
            break;
        }
        x_old = x;
    }
}
void predictVelocityFieldGaussSeidl(Data2D& data) {
    for (int i = 0; i < data.nFaces; i++) {
        Face2D *curFace = &data.faces[i];
        if (curFace->bType_u == INNERCELL) {                        
            if (i < data.nhorizontalFaces) {
                computeVelocityCoeff_y(data, i);
            } else if (i >= data.nhorizontalFaces) {
                computeVelocityCoeff_x(data, i);
            }
        }
    }
    
    // Initialize matrices and vectors
    SpMat momMatrix(data.nFaces, data.nFaces);
    Vector momVector(data.nFaces);
    momVector.setZero();

    // Set up the momentum equation matrices
    setMomentumEquationMatrix(data, momMatrix, momVector);

    // Compress matrices
    momMatrix.makeCompressed();
    
    // if (data.timeStep == 30) {
    //     std::cout << "uMatrix: " << uMatrix << std::endl;
    //     std::cout << "vMatrix: " << vMatrix << std::endl;
    // }

    // Initialize solution vectors
    Vector solution = Vector::Zero(data.nFaces);

    // Gauss-Seidel solver parameters
    int maxIterations = 1000;
    double tolerance = 1e-6;
    
    // Solve for u and v velocities using Gauss-Seidel method
    gaussSeidelSolver(momMatrix, solution, momVector, maxIterations, tolerance);

        // Update the faces with the solution
        for (int i = 0; i < data.nFaces; i++) {
            Face2D *curFace = &data.faces[i];
            if (i < data.nhorizontalFaces) {
                curFace->v[INTERMEDIATE_1] = solution(i);
            } else {
                curFace->u[INTERMEDIATE_1] = solution(i);
            }
        }
}

///////////////////////////////////////////////
//         IMPLICIT EIGEN SOLVER METHOD      //
///////////////////////////////////////////////
void predictVelocityFieldBiCGStab(Data2D& data) {
    for (int i = 0; i < data.nFaces; i++) {
        Face2D *curFace = &data.faces[i];
        if (curFace->bType_u == INNERCELL) {                        
            if (i < data.nhorizontalFaces) {
                computeVelocityCoeff_y(data, i);
            } else if (i >= data.nhorizontalFaces) {
                computeVelocityCoeff_x(data, i);
            }
        }
    }

    // Initialize matrices and vectors
    SpMat momMatrix(data.nFaces, data.nFaces);
    Vector momVector(data.nFaces);
    momVector.setZero();

    // Set up the momentum equation matrices
    setMomentumEquationMatrix(data, momMatrix, momVector);

    // Compress matrices
    momMatrix.makeCompressed();

    try {
        // Solve for uMatrix
        Eigen::BiCGSTAB<SpMat> solver;
        solver.setTolerance(1e-6);
        solver.setMaxIterations(10000);
        solver.compute(momMatrix);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Decomposition failed for uMatrix");
        }
        Vector solution = solver.solve(momVector);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Solving failed for uMatrix");
        }

        // Update the faces with the solution
        for (int i = 0; i < data.nFaces; i++) {
            Face2D *curFace = &data.faces[i];
            if (i < data.nhorizontalFaces) {
                curFace->v[INTERMEDIATE_1] = solution(i);
            } else {
                curFace->u[INTERMEDIATE_1] = solution(i);
            }
        }
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

void predictVelocityFieldSparseLU(Data2D& data) {
    for (int i = 0; i < data.nFaces; i++) {
        Face2D *curFace = &data.faces[i];
        if (curFace->bType_u == INNERCELL) {                        
            if (i < data.nhorizontalFaces) {
                computeVelocityCoeff_y(data, i);
            } else if (i >= data.nhorizontalFaces) {
                computeVelocityCoeff_x(data, i);
            }
        }
    }

    // Initialize matrices and vectors
    SpMat momMatrix(data.nFaces, data.nFaces);
    Vector momVector(data.nFaces);
    momVector.setZero();

    // Set up the momentum equation matrices
    setMomentumEquationMatrix(data, momMatrix, momVector);

    // Compress matrices
    momMatrix.makeCompressed();

    //if (data.timeStep == 1) {
    //      std::cout << "momMatrix: " << momMatrix << std::endl;
    //}

    try {
        // Solve for momMatrix
        Eigen::SparseLU<SpMat> solver;
        solver.analyzePattern(momMatrix);
        solver.factorize(momMatrix);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Factorization failed for momMatrix");
        }
        Vector solution = solver.solve(momVector);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Solving failed for momMatrix");
        }

        // Update the faces with the solution
        for (int i = 0; i < data.nFaces; i++) {
            Face2D *curFace = &data.faces[i];
            if (i < data.nhorizontalFaces) {
                curFace->v[INTERMEDIATE_1] = solution(i);
            } else {
                curFace->u[INTERMEDIATE_1] = solution(i);
            }
        }
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

///////////////////////////////////////////////
//           EXPLICIT METHOD                 //
///////////////////////////////////////////////
void predictVelocityFieldExplicit(Data2D& data) {
    for (int i = 0; i < data.nFaces; i++) {
        Face2D *curFace = &data.faces[i];
        if (curFace->bType_u == INNERCELL) {                        
            if (i < data.nhorizontalFaces) {
                computeVelocityCoeff_y(data, i);
            }
            else if (i >= data.nhorizontalFaces) {
                computeVelocityCoeff_x(data, i);
            }
        }
    }

    for(int i = 0; i < data.nFaces; i++){
        Face2D *curFace = &data.faces[i];
        if (curFace->bType_u == INNERCELL) {
            if (i < data.nhorizontalFaces) {
                curFace->v[INTERMEDIATE_1] = (1/curFace->a_p_tilde) * (curFace->b 
                    + curFace->a_w * curFace->neighCells[UP]->neighCells[WEST]->faces[SOUTH]->v[INITIAL] 
                    + curFace->a_e * curFace->neighCells[UP]->neighCells[EAST]->faces[SOUTH]->v[INITIAL] 
                    + curFace->a_n * curFace->neighCells[UP]->faces[NORTH]->v[INITIAL] 
                    + curFace->a_s * curFace->neighCells[DOWN]->faces[SOUTH]->v[INITIAL]); ;
            } else {
                curFace->u[INTERMEDIATE_1] = (1/curFace->a_p_tilde) * (curFace->b 
                    + curFace->a_w * curFace->neighCells[LEFT]->faces[WEST]->u[INITIAL] 
                    + curFace->a_e * curFace->neighCells[RIGHT]->faces[EAST]->u[INITIAL] 
                    + curFace->a_n * curFace->neighCells[LEFT]->neighCells[NORTH]->faces[EAST]->u[INITIAL] 
                    + curFace->a_s * curFace->neighCells[RIGHT]->neighCells[SOUTH]->faces[WEST]->u[INITIAL]);
            }
        } else if (curFace->bType_u == DIRICHLET || curFace->bType_u == SOLID) {
            if (i < data.nhorizontalFaces) {
                curFace->v[INTERMEDIATE_1] = curFace->v[INITIAL];
            } else {
                curFace->u[INTERMEDIATE_1] = curFace->u[INITIAL];
            }
        } 
    }
    for(int i = 0; i < data.nFaces; i++){
    Face2D *curFace = &data.faces[i];
        if (curFace->bType_u == NEUMANN) {
            if (i < data.nhorizontalFaces) {
                if (curFace->neighCells[DOWN] == nullptr || curFace->neighCells[DOWN]->bType_sc == SOLID) {
                    curFace->v[INTERMEDIATE_1] = curFace->neighCells[UP]->faces[NORTH]->v[INTERMEDIATE_1];
                } else if (curFace->neighCells[UP] == nullptr || curFace->neighCells[UP]->bType_sc == SOLID) {
                    curFace->v[INTERMEDIATE_1] = curFace->neighCells[DOWN]->faces[SOUTH]->v[INTERMEDIATE_1];
                }
            } else if (i >= data.nhorizontalFaces) {
                if (curFace->neighCells[LEFT] == nullptr || curFace->neighCells[LEFT]->bType_sc == SOLID) {
                    curFace->u[INTERMEDIATE_1] = curFace->neighCells[RIGHT]->faces[EAST]->u[INTERMEDIATE_1];
                } else if (curFace->neighCells[RIGHT] == nullptr || curFace->neighCells[RIGHT]->bType_sc == SOLID) {
                    curFace->u[INTERMEDIATE_1] = curFace->neighCells[LEFT]->faces[WEST]->u[INTERMEDIATE_1];
                }
            }
        }
    }
}

void predictVelocityField(Data2D& data){
    if (data.velSolver == 3){
        predictVelocityFieldExplicit(data);
    } else if (data.velSolver == 0){
        predictVelocityFieldBiCGStab(data);
    } else if (data.velSolver == 1){
        predictVelocityFieldSparseLU(data);
    } else if (data.velSolver == 2){
        predictVelocityFieldGaussSeidl(data);
    }
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                             CORRECTOR STEP 1                                                   //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

void computeScalarCoeff(Data2D& data, int cellId){
    Cell2D *curCell = &data.cells[cellId];
    double dx = curCell->faces[SOUTH]->dx;
    double dy = curCell->faces[WEST]->dy;
    double dt = data.dt;
    double d_e = fEta(data, curCell->neighCells[EAST]->alpha) / dx;
    double d_w = fEta(data, curCell->neighCells[WEST]->alpha) / dx;
    double d_n = fEta(data, curCell->neighCells[NORTH]->alpha) / dy;
    double d_s = fEta(data, curCell->neighCells[SOUTH]->alpha) / dy;
    double f_e = 0.5 * (fRho(data, curCell->neighCells[EAST]->alpha) + fRho(data, curCell->alpha)) * curCell->faces[EAST]->u[CORRECTED_2];
    double f_w = 0.5 * (fRho(data, curCell->neighCells[WEST]->alpha) + fRho(data, curCell->alpha)) * curCell->faces[WEST]->u[CORRECTED_2];
    double g_n = 0.5 * (fRho(data, curCell->neighCells[NORTH]->alpha) + fRho(data, curCell->alpha)) * curCell->faces[NORTH]->v[CORRECTED_2];
    double g_s = 0.5 * (fRho(data, curCell->neighCells[SOUTH]->alpha) + fRho(data, curCell->alpha)) * curCell->faces[SOUTH]->v[CORRECTED_2];
    curCell->a_e_sc = d_e * dy * A(data, f_e / d_e) + std::max(-f_e*dy, 0.0);
    curCell->a_w_sc = d_w * dy * A(data, f_w / d_w) + std::max(f_w*dy, 0.0);
    curCell->a_n_sc = d_n * dx * A(data, g_n / d_n) + std::max(-g_n*dx, 0.0);
    curCell->a_s_sc = d_s * dx * A(data, g_s / d_s) + std::max(g_s*dx, 0.0);
    curCell->b_sc = curCell->a_p_v_sc*curCell->sc;
    curCell->a_p_v_sc = fRho(data, curCell->neighCells[NORTH]->alpha) * dx * dy / dt;
    curCell->a_p_sc = curCell->a_e_sc + curCell->a_w_sc + curCell->a_n_sc + curCell->a_s_sc + curCell->a_p_v_sc;
}

void calcScalarTransfer(Data2D& data){
    for (int i = 0; i < data.nCells; i++) {
        Cell2D *curCell = &data.cells[i];
        if (curCell->bType_sc == INNERCELL) {
            computeScalarCoeff(data, i);
        }
    }
    for (int i = 0; i < data.nCells; i++) {
        Cell2D *curCell = &data.cells[i];
        if (curCell->bType_sc == INNERCELL) {
            curCell->sc = (1/curCell->a_p_sc) * (curCell->b_sc + curCell->a_w_sc * curCell->neighCells[WEST]->sc + curCell->a_e_sc * curCell->neighCells[EAST]->sc + curCell->a_n_sc * curCell->neighCells[NORTH]->sc + curCell->a_s_sc * curCell->neighCells[SOUTH]->sc);
        } else if (curCell->bType_sc == DIRICHLET || curCell->bType_sc == SOLID) {
            curCell->sc = curCell->sc;
        } 
    }
    for (int i = 0; i < data.nCells; i++) {
        Cell2D *curCell = &data.cells[i];
        if (curCell->bType_sc == NEUMANN) {
            if (curCell->neighCells[EAST] == nullptr || curCell->neighCells[EAST]->bType_sc == SOLID) {
                curCell->sc = curCell->neighCells[WEST]->sc;
            }
            if (curCell->neighCells[WEST] == nullptr || curCell->neighCells[WEST]->bType_sc == SOLID) {
                curCell->sc = curCell->neighCells[EAST]->sc;
            }
            if (curCell->neighCells[NORTH] == nullptr || curCell->neighCells[NORTH]->bType_sc == SOLID) {
                curCell->sc = curCell->neighCells[SOUTH]->sc;
            }
            if (curCell->neighCells[SOUTH] == nullptr || curCell->neighCells[SOUTH]->bType_sc == SOLID) {
                curCell->sc = curCell->neighCells[NORTH]->sc;
            }
        }
    }
}

double continutyResidual(Data2D& data, int cellId){
    Cell2D *curCell = &data.cells[cellId];
    double rho_e;
    double rho_w;
    double rho_n;
    double rho_s;

    if (curCell->neighCells[EAST] == nullptr) {
        rho_e = fRho(data, curCell->alpha);
    } else {
        rho_e = 0.5 * (fRho(data, curCell->neighCells[EAST]->alpha) + fRho(data, curCell->alpha));
    }
    if (curCell->neighCells[WEST] == nullptr) {
        rho_w = fRho(data, curCell->alpha);
    } else {
        rho_w = 0.5 * (fRho(data, curCell->neighCells[WEST]->alpha) + fRho(data, curCell->alpha));
    }

    if (curCell->neighCells[NORTH] == nullptr) {
        rho_n = fRho(data, curCell->alpha);
    } else {
        rho_n = 0.5 * (fRho(data, curCell->neighCells[NORTH]->alpha) + fRho(data, curCell->alpha));
    }
    if (curCell->neighCells[SOUTH] == nullptr) {
        rho_s = fRho(data, curCell->alpha);
    } else {
        rho_s = 0.5 * (fRho(data, curCell->neighCells[SOUTH]->alpha) + fRho(data, curCell->alpha));
    }
 
    double flux_e = rho_e * curCell->faces[EAST]->u[CORRECTED_2] * curCell->faces[EAST]->dy;
    double flux_w = rho_w * curCell->faces[WEST]->u[CORRECTED_2] * curCell->faces[WEST]->dy;
    double flux_n = rho_n * curCell->faces[NORTH]->v[CORRECTED_2] * curCell->faces[NORTH]->dx;
    double flux_s = rho_s * curCell->faces[SOUTH]->v[CORRECTED_2] * curCell->faces[SOUTH]->dx;
    double imbalance = (flux_e - flux_w) + (flux_n - flux_s);
    return imbalance;
}

void checkConvergence(Data2D& data, int iteration){
    double res = 0.0;
    for (int i=0; i < data.nCells; i++){
        res += continutyResidual(data, i);
    }
    double rmsRes = sqrt(abs(res)/data.nCells);
    data.continuityResiduals.push(rmsRes);
    if (rmsRes > 1e-6 && iteration < data.maxIteration){
        data.stackOfContinuityResiduals.push(data.continuityResiduals);
        iteration++;
        resetData(data);
        iterateSteady(data, iteration);
    } else {
        data.stackOfContinuityResiduals.push(data.continuityResiduals);
        while(!data.continuityResiduals.empty()){
            data.continuityResiduals.pop();
        }
        if (iteration == data.maxIteration){
            std::cout << "Maximum number of iterations reached." << std::endl;
        } else {
            std::cout << "Converged after " << iteration << " iterations." << std::endl;
        }
    }
}

void unstaggerGrid(Data2D& data, int step) {
    for (int i = 0; i < data.nCells; i++) {
        Cell2D *curCell = &data.cells[i];
        curCell->u[step] = 0.5 * (curCell->faces[EAST]->u[step]+ curCell->faces[WEST]->u[step]);
        curCell->v[step] = 0.5 * (curCell->faces[NORTH]->v[step] + curCell->faces[SOUTH]->v[step]);
    }
}

void assignPrevData(Data2D& data){
    for (int i = 0; i < data.nFaces; i++) {
        Face2D *curFace = &data.faces[i];
        if (i < data.nhorizontalFaces) {
            curFace->v_prev = curFace->v[CORRECTED_2];
        } else {
            curFace->u_prev = curFace->u[CORRECTED_2];
        }
    }
    for (int i = 0; i < data.nCells; i++) {
        Cell2D *curCell = &data.cells[i];
        curCell->alpha_prev = curCell->p[CORRECTED_2];
    }
}

void resetData(Data2D& data) {
    for (int i = 0; i < data.nCells; i++) {
        Cell2D *curCell = &data.cells[i];
        if (curCell->bType_p == INNERCELL) {
            curCell->p[INITIAL] = curCell->p[CORRECTED_2];
        } else if (curCell->bType_p == DIRICHLET || curCell->bType_p == SOLID) {
            continue;
        } else if (curCell->bType_p == NEUMANN) {
            curCell->p[INITIAL] = curCell->p[CORRECTED_2];
        }

    }
    for (int i = 0; i < data.nFaces; i++) {
        Face2D *curFace = &data.faces[i];
        if (curFace->bType_u == INNERCELL) {
            curFace->u[INITIAL] = curFace->u[CORRECTED_2];
            curFace->v[INITIAL] = curFace->v[CORRECTED_2];
        } else if (curFace->bType_u == DIRICHLET || curFace->bType_u == SOLID) {
            continue;
        } else if (curFace->bType_u == NEUMANN) {
            curFace->u[INITIAL] = curFace->u[CORRECTED_2];
            curFace->v[INITIAL] = curFace->v[CORRECTED_2];
        }
    }
}