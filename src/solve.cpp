#include "data.h"
#include "solve.h"
#include "output.h"
#include "velPredict.h"
#include "corrector.h"
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
    if (data.fixedPressure)
    {
        assignVelocities(data, INTERMEDIATE_1);
    } else {
        correctPressureEquation(data, INTERMEDIATE_1);
        corrector1(data);

        correctPressureEquation(data, INTERMEDIATE_2);
        corrector2(data);
    }

    //calcScalarTransfer(data);
    checkConvergence(data, iteration);
}

void iterateTransient(Data2D& data){
    assignPrevData(data);
    predictVelocityField(data);

    correctPressureEquation(data, INTERMEDIATE_1);
    corrector1(data);

    correctPressureEquation(data, INTERMEDIATE_2);
    corrector2(data);

    //calcScalarTransfer(data);
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
    curCell->b_sc = curCell->a_p_v_sc*curCell->alpha;
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
            curCell->alpha = (1/curCell->a_p_sc) * (curCell->b_sc + curCell->a_w_sc * curCell->neighCells[WEST]->sc + curCell->a_e_sc * curCell->neighCells[EAST]->sc + curCell->a_n_sc * curCell->neighCells[NORTH]->sc + curCell->a_s_sc * curCell->neighCells[SOUTH]->sc);
        } else if (curCell->bType_sc == DIRICHLET || curCell->bType_sc == SOLID) {
            curCell->alpha = curCell->alpha;
        } 
    }
    for (int i = 0; i < data.nCells; i++) {
        Cell2D *curCell = &data.cells[i];
        if (curCell->bType_sc == NEUMANN) {
            if (curCell->neighCells[EAST] == nullptr || curCell->neighCells[EAST]->bType_sc == SOLID) {
                curCell->alpha = curCell->neighCells[WEST]->alpha + curCell->g_sc * curCell->faces[WEST]->dx;
            }
            if (curCell->neighCells[WEST] == nullptr || curCell->neighCells[WEST]->bType_sc == SOLID) {
                curCell->sc = curCell->neighCells[EAST]->alpha - curCell->g_sc * curCell->faces[EAST]->dx;
            }
            if (curCell->neighCells[NORTH] == nullptr || curCell->neighCells[NORTH]->bType_sc == SOLID) {
                curCell->sc = curCell->neighCells[SOUTH]->alpha + curCell->g_sc * curCell->faces[SOUTH]->dy;
            }
            if (curCell->neighCells[SOUTH] == nullptr || curCell->neighCells[SOUTH]->bType_sc == SOLID) {
                curCell->sc = curCell->neighCells[NORTH]->alpha - curCell->g_sc * curCell->faces[NORTH]->dy;
            }
        }
    }
}

double continutyResidual(Data2D& data, int cellId, int step){
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
 
    double flux_e = rho_e * curCell->faces[EAST]->u[step] * curCell->faces[EAST]->dy;
    double flux_w = rho_w * curCell->faces[WEST]->u[step] * curCell->faces[WEST]->dy;
    double flux_n = rho_n * curCell->faces[NORTH]->v[step] * curCell->faces[NORTH]->dx;
    double flux_s = rho_s * curCell->faces[SOUTH]->v[step] * curCell->faces[SOUTH]->dx;
    double imbalance = (flux_e - flux_w) + (flux_n - flux_s);
    return imbalance;
}

void checkConvergence(Data2D& data, int iteration){
    double res = 0.0;
    for (int i=0; i < data.nCells; i++){
        res += continutyResidual(data, i, CORRECTED_1);
    }
    double rmsRes = sqrt(abs(res)/data.nCells);
    data.continuityResiduals.push(rmsRes);
    if (rmsRes > 1e-8 && iteration < data.maxIteration){
    //if (iteration < data.maxIteration){
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

void assignVelocities(Data2D& data, int step){
    for (int i = 0; i < data.nFaces; i++) {
        data.faces[i].u[CORRECTED_2] = data.faces[i].u[step];
        data.faces[i].v[CORRECTED_2] = data.faces[i].v[step];
    }
    for (int i = 0; i < data.nCells; i++) {
        data.cells[i].p[CORRECTED_2] = data.cells[i].p[INITIAL];
    }
}