#include "solve.h"
#include "data.h"
#include "velPredict.h"

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
        f_e = fRho(data, curFace->neighCells[RIGHT]->alpha) * u_E * dy; //In case face is one away of the right boundary
    } else {
        f_e = (fRho(data, curFace->neighCells[RIGHT]->alpha) + fRho(data, curFace->neighCells[RIGHT]->neighCells[EAST]->alpha)) / 2 * u_E * dy; 
    }

    if (curFace->neighCells[LEFT]->neighCells[WEST] == nullptr) {
        f_w = fRho(data, curFace->neighCells[LEFT]->alpha) * u_W * dy; //In case face is one away of the left boundary
    } else {
        f_w = (fRho(data, curFace->neighCells[LEFT]->alpha) + fRho(data, curFace->neighCells[LEFT]->neighCells[WEST]->alpha)) / 2 * u_W * dy;
    }

    double g_n = 0.5 * dx * ((fRho(data, curFace->neighCells[RIGHT]->alpha) + fRho(data, curFace->neighCells[RIGHT]->neighCells[NORTH]->alpha)) / 2 * v_N_R + (fRho(data, curFace->neighCells[LEFT]->alpha) + fRho(data, curFace->neighCells[LEFT]->neighCells[NORTH]->alpha)) / 2 * v_N_L);
    double g_s = 0.5 * dx * ((fRho(data, curFace->neighCells[RIGHT]->alpha) + fRho(data, curFace->neighCells[RIGHT]->neighCells[SOUTH]->alpha)) / 2 * v_S_R + (fRho(data, curFace->neighCells[LEFT]->alpha) + fRho(data, curFace->neighCells[LEFT]->neighCells[SOUTH]->alpha)) / 2 * v_S_L);

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
    double f_e = 0.5 * dy * ((fRho(data, curFace->neighCells[UP]->alpha) + fRho(data, curFace->neighCells[UP]->neighCells[EAST]->alpha)) / 2 * u_E_U + (fRho(data, curFace->neighCells[DOWN]->alpha) + fRho(data, curFace->neighCells[DOWN]->neighCells[EAST]->alpha)) / 2 * u_E_D);
    double f_w = 0.5 * dy * ((fRho(data, curFace->neighCells[UP]->alpha) + fRho(data, curFace->neighCells[UP]->neighCells[WEST]->alpha)) / 2 * u_W_U + (fRho(data, curFace->neighCells[DOWN]->alpha) + fRho(data, curFace->neighCells[DOWN]->neighCells[WEST]->alpha)) / 2 * u_W_D);

    double g_n;
    double g_s;

    if (curFace->neighCells[UP]->neighCells[NORTH] == nullptr) {
        g_n = fRho(data, curFace->neighCells[UP]->alpha) * dx * v_N; //In case face is one away of the top boundary
    } else {
        g_n = (fRho(data, curFace->neighCells[UP]->alpha) + fRho(data, curFace->neighCells[UP]->neighCells[NORTH]->alpha)) / 2 * v_N * dx;
    }

    if (curFace->neighCells[DOWN]->neighCells[SOUTH] == nullptr) {
        g_s = fRho(data, curFace->neighCells[DOWN]->alpha) * v_S * dx; //In case face is one away of the bottom boundary
    } else {
        g_s = (fRho(data, curFace->neighCells[DOWN]->alpha) + fRho(data, curFace->neighCells[DOWN]->neighCells[SOUTH]->alpha)) / 2 * v_S * dx;
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
    // curFace->b += 9.81 * (fRho(data, curFace->neighCells[UP]->alpha) * dx * curFace->neighCells[UP]->faces[WEST]->dy + fRho(data, curFace->neighCells[DOWN]->alpha) * dx * curFace->neighCells[DOWN]->faces[WEST]->dy) / 2;
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
            momVector(i) = -curFace->b;
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
            momMatrix.insert(i, i) = 1.0;
            if (i < data.nhorizontalFaces) {
                if (curFace->neighCells[UP] != nullptr) {               // Bottom Boundary
                    int j = curFace->neighCells[UP]->faces[NORTH]->id;
                    momMatrix.insert(i, j) = -1.0;
                    momVector(i) = -curFace->g_v*curFace->dy;
                }
                if (curFace->neighCells[DOWN] != nullptr) {             // Top Boundary
                    int j = curFace->neighCells[DOWN]->faces[SOUTH]->id;
                    momMatrix.insert(i, j) = -1.0;
                    momVector(i) = curFace->g_v*curFace->dy;
                }
            } else {
                if (curFace->neighCells[LEFT] != nullptr) {             // Right Boundary
                    int j = curFace->neighCells[LEFT]->faces[WEST]->id;
                    momMatrix.insert(i, j) = -1.0;
                    momVector(i) = curFace->g_u*curFace->dx;
                }
                if (curFace->neighCells[RIGHT] != nullptr) {            // Left Boundary
                    int j = curFace->neighCells[RIGHT]->faces[EAST]->id;
                    momMatrix.insert(i, j) = -1.0;
                    momVector(i) = -curFace->g_u*curFace->dx;
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
            throw std::runtime_error("Decomposition failed for momMatrix");
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