#ifndef OUTPUT_H
#define OUTPUT_H

#include "data.h"


bool fOutputVTKframe(Data2D& data, std::string finalDirectoryName, int step);
std::string createDirectory();
void jsonOutput(Data2D& data);
void residualPlot(Data2D& data);

#endif