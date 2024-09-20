#include "solve.h"
#include "data.h"
#include "velPredict.h"
#include <unsupported/Eigen/IterativeSolvers>

/**
 * Computes the velocity coefficients in the x-direction for a given face to solve the momentum equation.
 * 
 * @param data The data structure containing the simulation data.
 * @param faceId The ID of the face for which to compute the velocity coefficients.
 * @param step The current step of the PISO algorithm.
 */
void computeVelocityCoeff_x(Data2D& data, int faceId) {
    Face2D *curFace = &data.faces[faceId];

    double rho_W = 0;
    double rho_E = 0;
    double rho_N = 0;
    double rho_S = 0;
    double eta_W = 0;
    double eta_E = 0;
    double eta_N = 0;
    double eta_S = 0;

    // Velocity Terms (U)
    double u_W = (curFace->neighCells[LEFT]->faces[WEST]->u[INITIAL] + curFace->u[INITIAL]) / 2;
    double u_E = (curFace->neighCells[RIGHT]->faces[EAST]->u[INITIAL] + curFace->u[INITIAL]) / 2;

    // Velocity Terms (V)
    double v_S = (curFace->neighCells[LEFT]->faces[SOUTH]->v[INITIAL] + curFace->neighCells[RIGHT]->faces[SOUTH]->v[INITIAL]) / 2;
    double v_N = (curFace->neighCells[LEFT]->faces[NORTH]->v[INITIAL] + curFace->neighCells[RIGHT]->faces[NORTH]->v[INITIAL]) / 2;

    // Density terms
    rho_W = fRho(data, curFace->neighCells[LEFT]->alpha);
    rho_E = fRho(data, curFace->neighCells[RIGHT]->alpha);
    if (curFace->neighCells[LEFT]->neighCells[NORTH] == nullptr && curFace->neighCells[RIGHT]->neighCells[NORTH] == nullptr) {
        rho_N = (fRho(data, curFace->neighCells[LEFT]->alpha) + fRho(data, curFace->neighCells[RIGHT]->alpha)) / 2;
    } else {
        rho_N = (fRho(data, curFace->neighCells[LEFT]->alpha) + fRho(data, curFace->neighCells[RIGHT]->alpha) + fRho(data, curFace->neighCells[RIGHT]->neighCells[NORTH]->alpha) + fRho(data, curFace->neighCells[LEFT]->neighCells[NORTH]->alpha)) / 4;
    }
    if (curFace->neighCells[LEFT]->neighCells[SOUTH] == nullptr && curFace->neighCells[RIGHT]->neighCells[SOUTH] == nullptr) {
        rho_S = (fRho(data, curFace->neighCells[LEFT]->alpha) + fRho(data, curFace->neighCells[RIGHT]->alpha)) / 2;
    } else {
        rho_S = (fRho(data, curFace->neighCells[LEFT]->alpha) + fRho(data, curFace->neighCells[RIGHT]->alpha) + fRho(data, curFace->neighCells[RIGHT]->neighCells[SOUTH]->alpha) + fRho(data, curFace->neighCells[LEFT]->neighCells[SOUTH]->alpha)) / 4;
    }


    // Viscosity Terms
    eta_W = fEta(data, curFace->neighCells[LEFT]->alpha);
    eta_E = fEta(data, curFace->neighCells[RIGHT]->alpha);
    if (curFace->neighCells[LEFT]->neighCells[NORTH] == nullptr && curFace->neighCells[RIGHT]->neighCells[NORTH] == nullptr) {
        eta_N = (fEta(data, curFace->neighCells[LEFT]->alpha) + fEta(data, curFace->neighCells[RIGHT]->alpha)) / 2;
    } else {
        eta_N = (fEta(data, curFace->neighCells[LEFT]->alpha) + fEta(data, curFace->neighCells[RIGHT]->alpha) + fEta(data, curFace->neighCells[RIGHT]->neighCells[NORTH]->alpha) + fEta(data, curFace->neighCells[LEFT]->neighCells[NORTH]->alpha)) / 4;
    }
    if (curFace->neighCells[LEFT]->neighCells[SOUTH] == nullptr && curFace->neighCells[RIGHT]->neighCells[SOUTH] == nullptr) {
        eta_S = (fEta(data, curFace->neighCells[LEFT]->alpha) + fEta(data, curFace->neighCells[RIGHT]->alpha)) / 2;
    } else {
        eta_S = (fEta(data, curFace->neighCells[LEFT]->alpha) + fEta(data, curFace->neighCells[RIGHT]->alpha) + fEta(data, curFace->neighCells[RIGHT]->neighCells[SOUTH]->alpha) + fEta(data, curFace->neighCells[LEFT]->neighCells[SOUTH]->alpha)) / 4;
    }


    // Cell Properties
    double dx = curFace->neighCells[LEFT]->faces[SOUTH]->dx;
    double dy = curFace->dy;
    double dt = data.dt;

    double Pe_e = (rho_E * u_E * dx) / eta_E;
    double Pe_w = (rho_W * u_W * dx) / eta_W;
    double Pe_n = (rho_N * v_N * dy) / eta_N;
    double Pe_s = (rho_S * v_S * dy) / eta_S;

    double f_e = rho_E * u_E * dy;
    double f_w = rho_W * u_W * dy;
    double f_n = rho_N * v_N * dx;
    double f_s = rho_S * v_S * dx;

    // Compute the coefficients
    curFace->a_e = -std::max(-f_e, 0.0) - A(data, Pe_e) * eta_E * dy / dx;
    curFace->a_w = -std::max(f_w, 0.0) - A(data, Pe_w) * eta_W * dy / dx;
    curFace->a_n = -std::max(-f_n, 0.0) - A(data, Pe_n) * eta_N * dx / dy;
    curFace->a_s = -std::max(f_s, 0.0) - A(data, Pe_s) * eta_S * dx / dy;

    curFace->a_p_v = 0.5 * (fRho(data, curFace->neighCells[LEFT]->alpha) + fRho(data, curFace->neighCells[RIGHT]->alpha)) * dx * dy / dt;

    curFace->a_p_tilde = -(curFace->a_e + curFace->a_w + curFace->a_n + curFace->a_s);

    curFace->b = (curFace->neighCells[LEFT]->p[INITIAL] - curFace->neighCells[RIGHT]->p[INITIAL]) * dy;

    if (data.mode == 0) {
        curFace->a_p_tilde += curFace->a_p_v;
        curFace->b += curFace->a_p_v * curFace->u_prev;
    }

    if (curFace->neighCells[LEFT]->faces[SOUTH]->bType_u == DIRICHLET) {
        curFace->a_p_tilde += curFace->a_s;
        curFace->a_p_tilde += (3 * eta_S * dx) / dy;
        curFace->a_s = 0;
        curFace->a_n -= (eta_S * dx) / (3 * dy);
        curFace->b += (8 * eta_S * dx * data.u_bottom) / (3 * dy);
        // curFace->a_e = 0;
        // curFace->a_w = 0;
        // curFace->a_s = 0;
        // curFace->a_n = 0;
        // curFace->a_p_tilde = 1;
        // curFace->b = data.u_bottom;
    }

    if (curFace->neighCells[RIGHT]->faces[NORTH]->bType_u == DIRICHLET) {
        curFace->a_p_tilde += curFace->a_n;
        curFace->a_p_tilde += (3 * eta_N * dx) / dy;
        curFace->a_n = 0;
        curFace->a_s -= (eta_N * dx) / (3 * dy);
        curFace->b += (8 * eta_N * dx * data.u_top) / (3 * dy);
    //     curFace->a_e = 0;
    //     curFace->a_w = 0;
    //     curFace->a_s = 0;
    //     curFace->a_n = 0;
    //     curFace->a_p_tilde = 1;
    //     curFace->b = data.u_top;
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

    double rho_W = 0;
    double rho_E = 0;
    double rho_N = 0;
    double rho_S = 0;

    double eta_W = 0;
    double eta_E = 0;
    double eta_N = 0;
    double eta_S = 0;

    // Velocity Terms (U)
    double u_W = (curFace->neighCells[DOWN]->faces[WEST]->u[INITIAL] + curFace->neighCells[UP]->faces[WEST]->u[INITIAL]) / 2;
    double u_E = (curFace->neighCells[DOWN]->faces[EAST]->u[INITIAL] + curFace->neighCells[UP]->faces[EAST]->u[INITIAL]) / 2;

    // Velocity Terms (V)
    double v_S = (curFace->v[INITIAL] + curFace->neighCells[DOWN]->faces[SOUTH]->v[INITIAL]) / 2;
    double v_N = (curFace->v[INITIAL] + curFace->neighCells[UP]->faces[NORTH]->v[INITIAL]) / 2;

    // Density terms
    if (curFace->neighCells[UP]->neighCells[WEST] == nullptr && curFace->neighCells[DOWN]->neighCells[WEST] == nullptr) {
        rho_W = (fRho(data, curFace->neighCells[DOWN]->alpha) + fRho(data, curFace->neighCells[UP]->alpha)) / 2;
    } else {
        rho_W = (fRho(data, curFace->neighCells[DOWN]->alpha) + fRho(data, curFace->neighCells[UP]->alpha) + fRho(data, curFace->neighCells[UP]->neighCells[WEST]->alpha) + fRho(data, curFace->neighCells[DOWN]->neighCells[WEST]->alpha)) / 4;
    }
    
    if (curFace->neighCells[UP]->neighCells[EAST] == nullptr && curFace->neighCells[DOWN]->neighCells[EAST] == nullptr) {
        rho_E = (fRho(data, curFace->neighCells[DOWN]->alpha) + fRho(data, curFace->neighCells[UP]->alpha)) / 2;
    } else {
        rho_E = (fRho(data, curFace->neighCells[DOWN]->alpha) + fRho(data, curFace->neighCells[UP]->alpha) + fRho(data, curFace->neighCells[UP]->neighCells[EAST]->alpha) + fRho(data, curFace->neighCells[DOWN]->neighCells[EAST]->alpha)) / 4;
    }
    
    rho_N = fRho(data, curFace->neighCells[UP]->alpha);
    rho_S = fRho(data, curFace->neighCells[DOWN]->alpha);

    // Viscosity Terms
    if (curFace->neighCells[UP]->neighCells[WEST] == nullptr && curFace->neighCells[DOWN]->neighCells[WEST] == nullptr) {
        eta_W = (fEta(data, curFace->neighCells[DOWN]->alpha) + fEta(data, curFace->neighCells[UP]->alpha)) / 2;
    } else {
        eta_W = (fEta(data, curFace->neighCells[DOWN]->alpha) + fEta(data, curFace->neighCells[UP]->alpha) + fEta(data, curFace->neighCells[UP]->neighCells[WEST]->alpha) + fEta(data, curFace->neighCells[DOWN]->neighCells[WEST]->alpha)) / 4;
    }

    if (curFace->neighCells[UP]->neighCells[EAST] == nullptr && curFace->neighCells[DOWN]->neighCells[EAST] == nullptr) {
        eta_E = (fEta(data, curFace->neighCells[DOWN]->alpha) + fEta(data, curFace->neighCells[UP]->alpha)) / 2;
    } else {
        eta_E = (fEta(data, curFace->neighCells[DOWN]->alpha) + fEta(data, curFace->neighCells[UP]->alpha) + fEta(data, curFace->neighCells[UP]->neighCells[EAST]->alpha) + fEta(data, curFace->neighCells[DOWN]->neighCells[EAST]->alpha)) / 4;
    }
    
    eta_N = fEta(data, curFace->neighCells[UP]->alpha);
    eta_S = fEta(data, curFace->neighCells[DOWN]->alpha);

    // Cell Properties
    double dx = curFace->dx;
    double dy = curFace->neighCells[UP]->faces[WEST]->dy;
    double dt = data.dt;

    double Pe_e = (rho_E * u_E * dx) / eta_E;
    double Pe_w = (rho_W * u_W * dx) / eta_W;
    double Pe_n = (rho_N * v_N * dy) / eta_N;
    double Pe_s = (rho_S * v_S * dy) / eta_S;

    double f_e = rho_E * u_E * dy;
    double f_w = rho_W * u_W * dy;
    double f_n = rho_N * v_N * dx;
    double f_s = rho_S * v_S * dx;

    curFace->a_e = -std::max(-f_e, 0.0) - A(data, Pe_e) * eta_E * dy / dx;
    curFace->a_w = -std::max(f_w, 0.0) - A(data, Pe_w) * eta_W * dy / dx;
    curFace->a_n = -std::max(-f_n, 0.0) - A(data, Pe_n) * eta_N * dx / dy;
    curFace->a_s = -std::max(f_s, 0.0) - A(data, Pe_s) * eta_S * dx / dy;

    curFace->a_p_v = 0.5 * (fRho(data, curFace->neighCells[LEFT]->alpha) + fRho(data, curFace->neighCells[RIGHT]->alpha)) * dx * dy / dt;
    curFace->a_p_tilde =  -(curFace->a_e + curFace->a_w + curFace->a_n + curFace->a_s);
    
    // Compute source term with pressure gradients
    curFace->b = (curFace->neighCells[DOWN]->p[INITIAL] - curFace->neighCells[UP]->p[INITIAL]) * dx;
    curFace->b -= 9.81 * (fRho(data, curFace->neighCells[UP]->alpha) * dx * curFace->neighCells[UP]->faces[WEST]->dy + fRho(data, curFace->neighCells[DOWN]->alpha) * dx * curFace->neighCells[DOWN]->faces[WEST]->dy) / 2;

    if (data.mode == 0) {
        curFace->a_p_tilde += curFace->a_p_v;
        curFace->b += curFace->a_p_v * curFace->v_prev;
    }

    // Handle Dirichlet boundary conditions
    if (curFace->neighCells[UP]->faces[EAST]->bType_u == DIRICHLET) {
        curFace->a_p_tilde += curFace->a_e;
        curFace->a_p_tilde += (3 * eta_E * dy) / dx;
        curFace->a_e = 0;
        curFace->a_w -= (eta_E * dy) / (3 * dx);
        curFace->b += (8 * eta_E * dy * data.v_right) / (3 * dx);
    }
    if (curFace->neighCells[UP]->faces[WEST]->bType_u == DIRICHLET) {
        curFace->a_p_tilde += curFace->a_w;
        curFace->a_p_tilde += (3 * eta_W * dy) / dx;
        curFace->a_w = 0;
        curFace->a_e -= (eta_W * dy) / (3 * dx);
        curFace->b += (8 * eta_W * dy * data.v_left) / (3 * dx);
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
void setMomentumEquationYMatrix(Data2D& data, SpMat& momMatrix, Vector& momVector){
    for (int i = 0; i < data.nhorizontalFaces; i++) {
        Face2D* curFace = &data.faces[i];
        if (curFace->bType_u == INNERCELL) {
            momMatrix.insert(i, i) = curFace->a_p_tilde/data.yInertiaDamper;
            momVector(i) = curFace->b;
            if (i > 0) {
                momMatrix.insert(i, i - 1) = curFace->a_w;
            }
            if (i < data.nhorizontalFaces - 1) {
                momMatrix.insert(i, i + 1) = curFace->a_e;
            }
            if (curFace->neighCells[UP] != nullptr) {
                int j = curFace->neighCells[UP]->faces[NORTH]->id;
                momMatrix.insert(i, j) = curFace->a_n;
            }
            if (curFace->neighCells[DOWN] != nullptr) {
                int j = curFace->neighCells[DOWN]->faces[SOUTH]->id;
                momMatrix.insert(i, j) = curFace->a_s;
            }
        } else if (curFace->bType_u == DIRICHLET || curFace->bType_u == SOLID) {
            momMatrix.insert(i, i) = 1.0;
            momVector(i) = curFace->v[INITIAL];
        } else if (curFace->bType_u == NEUMANN) {
            momMatrix.insert(i, i) = 1.0;
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
        }
    }
}

void setMomentumEquationXMatrix(Data2D& data, SpMat& momMatrix, Vector& momVector){
for (int k = data.nhorizontalFaces; k < data.nFaces; k++) {
        int i = k - data.nhorizontalFaces;
        Face2D* curFace = &data.faces[k];
        if (curFace->bType_u == INNERCELL) {
            momMatrix.insert(i, i) = curFace->a_p_tilde/data.xInertiaDamper;
            momVector(i) = curFace->b;
            if (i > 0) {
                momMatrix.insert(i, i - 1) = curFace->a_w;
            }
            if (i < data.nFaces - 1) {
                momMatrix.insert(i, i + 1) = curFace->a_e;
            }
            if (curFace->neighCells[LEFT]->neighCells[NORTH] != nullptr) {
                int j = curFace->neighCells[LEFT]->neighCells[NORTH]->faces[EAST]->id - data.nhorizontalFaces;
                momMatrix.insert(i, j) = curFace->a_n;
            }
            if (curFace->neighCells[RIGHT]->neighCells[SOUTH] != nullptr) {
                int j = curFace->neighCells[RIGHT]->neighCells[SOUTH]->faces[WEST]->id- data.nhorizontalFaces;
                momMatrix.insert(i, j) = curFace->a_s;
            }
        } else if (curFace->bType_u == DIRICHLET || curFace->bType_u == SOLID) {
            momMatrix.insert(i, i) = 1.0;
            momVector(i) = curFace->u[INITIAL];
        } else if (curFace->bType_u == NEUMANN) {
            momMatrix.insert(i, i) = 1.0;
            if (curFace->neighCells[LEFT] != nullptr) {             // Right Boundary
                int j = curFace->neighCells[LEFT]->faces[WEST]->id - data.nhorizontalFaces;
                momMatrix.insert(i, j) = -1.0;
                momVector(i) = curFace->g_u*curFace->dx;
            }
            if (curFace->neighCells[RIGHT] != nullptr) {            // Left Boundary
                int j = curFace->neighCells[RIGHT]->faces[EAST]->id - data.nhorizontalFaces;
                momMatrix.insert(i, j) = -1.0;
                momVector(i) = -curFace->g_u*curFace->dx;
            
            }
        }
    }
}

void setMomentumEquationMatrix(Data2D& data, SpMat& momMatrix, Vector& momVector){
    for (int i = 0; i < data.nFaces; i++) {
        Face2D* curFace = &data.faces[i];
        if (curFace->bType_u == INNERCELL) {
            momMatrix.insert(i, i) = curFace->a_p_tilde;
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
void predictXVelocityFieldBiCGStab(Data2D& data) {
    std::cout << std::endl;
    std::cout << "Velocity X Prediction using BiCGStab" << std::endl;
    for (int i = data.nhorizontalFaces; i < data.nFaces; i++) {
        Face2D *curFace = &data.faces[i];
        if (curFace->bType_u == INNERCELL) {                        
            computeVelocityCoeff_x(data, i);
        }
    }

    // Initialize matrices and vectors
    SpMat momMatrix(data.nFaces - data.nhorizontalFaces, data.nFaces - data.nhorizontalFaces);
    Vector momVector(data.nFaces - data.nhorizontalFaces);
    momVector.setZero();
    Vector momGuess = Vector::Zero(data.nFaces - data.nhorizontalFaces);
    for (int i = data.nhorizontalFaces; i < data.nFaces; i++) {
        momGuess(i - data.nhorizontalFaces) = data.faces[i].u[INITIAL];
    }

    // Set up the momentum equation matrices
    setMomentumEquationXMatrix(data, momMatrix, momVector);

    
    //calculate initial residual
    Vector residual = momVector - momMatrix * momGuess;
    //L2-norm of the residual
    double residual_norm = residual.norm();
    std::cout << "Initial residual: " << residual_norm << std::endl;
    data.momentumXResidual = residual_norm;

    // Compress matrices
    momMatrix.makeCompressed();
    // std::cout << "momMatrix: " << momMatrix << std::endl;
    // std::cout << "momVector: " << momVector << std::endl;

    try {
        // Solve for uMatrix
        Eigen::BiCGSTAB<SpMat, Eigen::DiagonalPreconditioner<double>> solver;
        //GMRES
        // Eigen::GMRES<SpMat, Eigen::DiagonalPreconditioner<double>> solver;
        solver.setMaxIterations(10000); // Increase iteration count if needed
        solver.setTolerance(1e-6);     // Adjust the tolerance to balance accuracy and speed
        solver.compute(momMatrix);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Decomposition failed for momMatrix X");
        }
        Vector solution = solver.solveWithGuess(momVector, momGuess);
        if (solver.info() == Eigen::NoConvergence) {
            throw std::runtime_error("Solving failed for momMatrix X: No convergence");
        } else if (solver.info() == Eigen::NumericalIssue) {
            throw std::runtime_error("Solving failed for momMatrix X: Numerical issues");
        } else if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Solving failed for momMatrix X");
        }
        std::cout << "Iterations: " << solver.iterations();
        std::cout << " || Estimated error: " << solver.error() << std::endl;


        if (solver.info() == Eigen::Success) {
            for (int i = data.nhorizontalFaces; i < data.nFaces; i++) {
                int j = i - data.nhorizontalFaces;
                Face2D *curFace = &data.faces[i];
                if(curFace->bType_u != DIRICHLET){
                    curFace->u[INTERMEDIATE_1] = solution(j);
                }
            }
        } else {
            std::cerr << "Solving failed for momMatrix X" << std::endl;
        }
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

void predictYVelocityFieldBiCGStab(Data2D& data) {
    std::cout << std::endl;
    std::cout << "Velocity Y Prediction using BiCGStab" << std::endl;
    for (int i = 0; i < data.nhorizontalFaces; i++) {
        Face2D *curFace = &data.faces[i];
        if (curFace->bType_u == INNERCELL) {                        
            computeVelocityCoeff_y(data, i);
        }
    }

    // Initialize matrices and vectors
    SpMat momMatrix(data.nhorizontalFaces, data.nhorizontalFaces);
    Vector momVector(data.nhorizontalFaces);
    momVector.setZero();
    Vector momGuess = Vector::Zero(data.nhorizontalFaces);
    for (int i = 0; i < data.nhorizontalFaces; i++) {
        momGuess(i) = data.faces[i].v[INITIAL];
    }

    // Set up the momentum equation matrices
    setMomentumEquationYMatrix(data, momMatrix, momVector);

    //calculate initial residual
    Vector residual = momVector - momMatrix * momGuess;
    //L2-norm of the residual
    double residual_norm = residual.norm();
    std::cout << "Initial residual: " << residual_norm << std::endl;
    data.momentumYResidual = residual_norm;

    // Compress matrices
    momMatrix.makeCompressed();
    // std::cout << "momMatrix: " << momMatrix << std::endl;
    // std::cout << "momVector: " << momVector << std::endl;

    try {
        // Solve for uMatrix
        Eigen::BiCGSTAB<SpMat, Eigen::DiagonalPreconditioner<double>> solver;
        //GMRES
        // Eigen::GMRES<SpMat, Eigen::DiagonalPreconditioner<double>> solver;
        solver.setMaxIterations(10000); // Increase iteration count if needed
        solver.setTolerance(1e-6);     // Adjust the tolerance to balance accuracy and speed
        solver.compute(momMatrix);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Decomposition failed for momMatrix Y");
        }
        Vector solution = solver.solveWithGuess(momVector, momGuess);
        if (solver.info() == Eigen::NoConvergence) {
            throw std::runtime_error("Solving failed for momMatrix Y: No convergence");
        } else if (solver.info() == Eigen::NumericalIssue) {
            throw std::runtime_error("Solving failed for momMatrix Y: Numerical issues");
        } else if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Solving failed for momMatrix Y");
        }
        std::cout << "Iterations: " << solver.iterations();
        std::cout << " || Estimated error: " << solver.error() << std::endl;

        
        // Update the faces with the solution
        if (solver.info() == Eigen::Success) {
            for (int i = 0; i < data.nhorizontalFaces; i++) {
                Face2D *curFace = &data.faces[i];
                if(curFace->bType_u != DIRICHLET){
                    curFace->v[INTERMEDIATE_1] = solution(i);
                }
            }
        } else {
            std::cerr << "Solving failed for momMatrix Y" << std::endl;
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
                if (curFace->neighCells[DOWN] == nullptr || curFace->neighCells[DOWN]->bType_p == SOLID) {
                    curFace->v[INTERMEDIATE_1] = curFace->neighCells[UP]->faces[NORTH]->v[INTERMEDIATE_1];
                } else if (curFace->neighCells[UP] == nullptr || curFace->neighCells[UP]->bType_p == SOLID) {
                    curFace->v[INTERMEDIATE_1] = curFace->neighCells[DOWN]->faces[SOUTH]->v[INTERMEDIATE_1];
                }
            } else if (i >= data.nhorizontalFaces) {
                if (curFace->neighCells[LEFT] == nullptr || curFace->neighCells[LEFT]->bType_p == SOLID) {
                    curFace->u[INTERMEDIATE_1] = curFace->neighCells[RIGHT]->faces[EAST]->u[INTERMEDIATE_1];
                } else if (curFace->neighCells[RIGHT] == nullptr || curFace->neighCells[RIGHT]->bType_p == SOLID) {
                    curFace->u[INTERMEDIATE_1] = curFace->neighCells[LEFT]->faces[WEST]->u[INTERMEDIATE_1];
                }
            }
        }
    }
}

void predictXVelocityField(Data2D& data){
    if (data.velSolver == 3){
        predictVelocityFieldExplicit(data);
    } else if (data.velSolver == 0){
        predictXVelocityFieldBiCGStab(data);
    } else if (data.velSolver == 1){
        predictVelocityFieldSparseLU(data);
    } else if (data.velSolver == 2){
        predictVelocityFieldGaussSeidl(data);
    }
}
void predictYVelocityField(Data2D& data){
    if (data.velSolver == 3){
        predictVelocityFieldExplicit(data);
    } else if (data.velSolver == 0){
        predictYVelocityFieldBiCGStab(data);
    } else if (data.velSolver == 1){
        predictVelocityFieldSparseLU(data);
    } else if (data.velSolver == 2){
        predictVelocityFieldGaussSeidl(data);
    }
}