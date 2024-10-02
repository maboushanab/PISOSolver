#ifndef SOLVE_H
#define SOLVE_H

#include "data.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>


typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> Triplet;
typedef Eigen::VectorXd Vector;

std::string fSolve(Data2D& data);
void iterateSteady(Data2D& data, int iteration);
void iterateTransient(Data2D& data, int iteration, bool iterateTime);

// Suplementary functions
double A(Data2D data, double Pe);
double fRho(Data2D& data, double alpha);
double fEta(Data2D& data, double alpha);

void calcCFL(Data2D& data);

// Convergence
double continutyResidual(Data2D& data, int cellId, int step);
void checkConvergence(Data2D& data, int iteration);
// Reset
void resetData(Data2D& data);
void assignPrevData(Data2D& data);

// Output
void unstaggerGrid(Data2D& data, int step);
void assignVelocities(Data2D& data, int step);

#endif