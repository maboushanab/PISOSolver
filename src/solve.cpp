#include "data.h"
#include "solve.h"
#include "output.h"
#include "velPredict.h"
#include "corrector.h"
#include "alphaAdvection.h"
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
            iterateTransient(data, 1, true);
        } else {                        // Steady State
            iterateSteady(data, 1);
        }


        unstaggerGrid(data, CORRECTED_2);
        fOutputVTKframe(data, directoryName, CORRECTED_2);
        calcCFL(data);
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
    predictXVelocityField(data);
    predictYVelocityField(data);
    if (data.fixedPressure)
    {
        assignVelocities(data, INTERMEDIATE_1);
    } else {
        correctPressureEquation(data, INTERMEDIATE_1);
        corrector1(data);

        // calcNeighbourSums(data, CORRECTED_1);
        // correctPressureEquation(data, INTERMEDIATE_2);
        corrector2(data);
    }

    checkConvergence(data, iteration);
    // advectAlpha(data);
}

void iterateTransient(Data2D& data, int iteration, bool iterateTime){
    if (iterateTime == true){
        assignPrevData(data);
    }
    predictXVelocityField(data);
    predictYVelocityField(data);
    if (data.fixedPressure)
    {
        assignVelocities(data, INTERMEDIATE_1);
    } else {
        correctPressureEquation(data, INTERMEDIATE_1);
        corrector1(data);

        calcNeighbourSums(data);
        correctPressureEquation(data, INTERMEDIATE_2);
        corrector2(data);
    }
        checkConvergence(data, iteration);
        if (data.isThereAlpha = true){
            advectAlpha(data);
        }
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
    else if (data.pecFunc == 3){
        if (std::abs(Pe) < 1e-6) res = 1.0;                                                  // Exponential
        else res = std::abs(Pe)/(std::exp(std::abs(Pe)) - 1);
    }
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
    double pressureRest = 0;
    double momentumXRest = 0;
    double momentumYRest = 0;
    for (int i = 0; i < data.nCells; i++) {
        Cell2D *curCell = &data.cells[i];
        pressureRest += std::pow((curCell->p[CORRECTED_2] - curCell->p[INITIAL]), 2.0);
    }
    for (int i = 0; i < data.nFaces; i++) {
        Face2D *curFace = &data.faces[i];
        momentumXRest += std::pow((curFace->u[CORRECTED_2] - curFace->u[INITIAL]), 2.0);
        momentumYRest += std::pow((curFace->v[CORRECTED_2] - curFace->v[INITIAL]), 2.0);
    }
    pressureRest = std::sqrt(pressureRest);
    momentumXRest = std::sqrt(momentumXRest);
    momentumYRest = std::sqrt(momentumYRest);
    std::cout << std::endl << "Pressure iteration residual: " << pressureRest << std::endl;
    std::cout << "Momentum X iteration residual: " << momentumXRest << std::endl;
    std::cout << "Momentum Y iteration residual: " << momentumYRest << std::endl;
    std::cout << std::endl << "======================================" << std::endl;
    data.pIterationRes.push_back(pressureRest);
    data.uIterationRes.push_back(momentumXRest);
    data.vIterationRes.push_back(momentumYRest);
        
    if ((std::abs(data.continuityResidual) > 1e-6 || std::abs(data.momentumXResidual) > 1e-4 || std::abs(data.momentumYResidual) > 1e-4 || pressureRest > 1e-7 || momentumXRest > 1e-7 || momentumYRest > 1e-7) &&  iteration < data.maxIteration){
        data.continuityResiduals.push_back(std::abs(data.continuityResidual));
        data.momentumXResiduals.push_back(std::abs(data.momentumXResidual));
        data.momentumYResiduals.push_back(std::abs(data.momentumYResidual));
        iteration++;
        resetData(data);
        if (data.mode == 0){
            iterateTransient(data, iteration, false);
        } else {
            iterateSteady(data, iteration);
        }
    } else {
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
        curCell->u[step] = 0.5 * (curCell->faces[EAST]->u[step] + curCell->faces[WEST]->u[step]);
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
        curCell->alpha_prev = curCell->alpha;
    }
}

void resetData(Data2D& data) {
    for (int i = 0; i < data.nCells; i++) {
        Cell2D *curCell = &data.cells[i];
        if (curCell->bType_p == INNERCELL || curCell->bType_p == NEUMANN) {
            curCell->p[INITIAL] = curCell->p[CORRECTED_2];
        } else if (curCell->bType_p == DIRICHLET || curCell->bType_p == SOLID) {
            continue;
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

void calcCFL(Data2D& data){
    double maxCFL = 0;
    double maxU = 0;
    double maxDx = 0;
    for (int i = 0; i < data.nFaces; i++) {
        Face2D *curFace = &data.faces[i];
        if(i < data.nhorizontalFaces){
            double v = curFace->v[CORRECTED_2];
            double dy = 0;
            if (curFace->neighCells[UP] != nullptr){
                dy = curFace->neighCells[UP]->faces[WEST]->dy;
            } else {
                dy = curFace->neighCells[DOWN]->faces[WEST]->dy;
            }
            if (v > maxU){
                maxU = v;
                maxDx = dy;
            } 
        } else {
            double u = curFace->u[CORRECTED_2];
            double dx = 0;
            if (curFace->neighCells[LEFT] != nullptr){
                dx = curFace->neighCells[LEFT]->faces[NORTH]->dx;
            } else {
                dx = curFace->neighCells[RIGHT]->faces[NORTH]->dx;
            }
            if (u > maxU){
                maxDx = dx;
                maxU = u;
            }
        }
    }
    maxCFL = std::abs(maxU) * data.dt / maxDx;
    std::cout << "CFL: " << maxCFL << std::endl;
    if (maxCFL > 0.1){
        data.dt = 0.1 * maxDx / std::abs(maxU);
    }
}