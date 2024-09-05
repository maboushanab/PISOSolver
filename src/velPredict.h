#ifndef VELPREDICT_H
#define VELPREDICT_H

#include "data.h"

void setMomentumEquationMatrix(Data2D& data, SpMat& momMatrix, Vector& momVector);
void setMomentumEquationXMatrix(Data2D& data, SpMat& momMatrix, Vector& momVector);
void setMomentumEquationYMatrix(Data2D& data, SpMat& momMatrix, Vector& momVector);
void computeVelocityCoeff_x(Data2D& data, int faceId);
void computeVelocityCoeff_y(Data2D& data, int faceId);
void predictVelocityFieldExplicit(Data2D& data);
void predictVelocityFieldSparseLU(Data2D& data);
void predictVelocityFieldBiCGStab(Data2D& data);
void predictVelocityFieldGaussSeidl(Data2D& data);
void predictVelocityField(Data2D& data);
void predictXVelocityField(Data2D& data);
void predictYVelocityField(Data2D& data);
void predictXVelocityFieldBiCGStab(Data2D& data);
void predictYVelocityFieldBiCGStab(Data2D& data);

#endif