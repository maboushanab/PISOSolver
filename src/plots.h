#ifndef PLOT_H
#define PLOT_H

#include "data.h"

void initResidualPlot(Data2D& data, std::string finalDirectoryName);
void outerLoopResidualPlot(Data2D& data, std::string finalDirectoryName);

#endif