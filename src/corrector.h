#ifndef CORRECTOR_H
#define CORRECTOR_H

#include "data.h"

// Corrector 1
void setPressureMatrix(Data2D& data, SpMat& pressureMatrix, Vector& pressureVector, int nCells, int step);
void correctPressureEquationBiCGStab(Data2D& data, int step);
void correctPressureEquationSparseLU(Data2D& data, int step);
void computePressureCoeff(Data2D& data, int cellId);
void correctPressureEquation(Data2D& data, int step);
void corrector1(Data2D& data);
void computePressureCoeff2(Data2D& data, int cellId);
void corrector2(Data2D& data);

void calcNeighbourSums(Data2D& data, int step);

#endif