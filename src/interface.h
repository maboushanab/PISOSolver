#ifndef INTERFACE_H
#define INTERFACE_H

#include "data.h"

void calcAlphaNormalVector3x3(Data2D& data, int cellId);
void initInterface(Data2D& data, int cellId);
double calcAlphaPrediction(Data2D& data, int cellId, bool print);
double calcAlphaSensitivity(Data2D& data, int cellId, double alpha);
void estimateInterfaceLine(Data2D& data, int cellId);
void reconstructInterfaceLines(Data2D& data);


#endif