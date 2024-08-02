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
void iterateTransient(Data2D& data);

// Suplementary functions
double A(Data2D data, double Pe);
double fRho(Data2D& data, double alpha);
double fEta(Data2D& data, double alpha);

// PISO algorithm
// Predictor
void setMomentumEquationMatrix(Data2D& data, SpMat& momMatrix, Vector& momVector);
void computeVelocityCoeff_x(Data2D& data, int faceId);
void computeVelocityCoeff_y(Data2D& data, int faceId);
void predictVelocityFieldExplicit(Data2D& data);
void predictVelocityFieldSparseLU(Data2D& data);
void predictVelocityFieldBiCGStab(Data2D& data);
void predictVelocityFieldGaussSeidl(Data2D& data);
void predictVelocityField(Data2D& data);
void setPressureMatrix(Data2D& data, SpMat& pressureMatrix, Vector& pressureVector, int nCells, int step);

// Corrector 1
void correctPressureEquationBiCGStab(Data2D& data, int step);
void correctPressureEquationSparseLU(Data2D& data, int step);
void computePressureCoeff(Data2D& data, int cellId);
void correctPressureEquation(Data2D& data, int step);
void corrector1(Data2D& data);

void computePressureCoeff2(Data2D& data, int cellId);
void corrector2(Data2D& data);

void computeScalarCoeff(Data2D& data, int cellId);
void calcScalarTransfer(Data2D& data);

// Convergence
double continutyResidual(Data2D& data, int cellId);
void checkConvergence(Data2D& data, int iteration);
// Reset
void resetData(Data2D& data);
void assignPrevData(Data2D& data);

// Output
void unstaggerGrid(Data2D& data, int step);

#endif