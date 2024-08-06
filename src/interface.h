#ifndef INTERFACE_H
#define INTERFACE_H

#include "data.h"

void calcAlphaNormalVector(Data2D& data, int cellId);
void initInterface(Data2D& data, int cellId);
double x(double y, double m, double n);
double y(double x, double m, double n);
double calcAlphaPrediction(Data2D& data, int cellId, double m, double n);
double calcAlphaSensitivity(Data2D& data, int cellId, double m, double n, double alpha);
void estimateInterfaceLine(Data2D& data, int cellId);
void reconstructInterfaceLines(Data2D& data);


#endif