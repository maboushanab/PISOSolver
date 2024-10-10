#ifndef ADVECTION_H
#define ADVECTION_H

#include "data.h"

void advectAlpha(Data2D& data);
void handleExcessAlpha(Data2D& data);
void preformYSweep(Data2D& data);
double calculateFluxLength(Data2D data, int cellId, std::string direction);
void calculatePolygonArea(Data2D& data, int cellId, Eigen::Vector2d p1, Eigen::Vector2d p2, Eigen::Vector2d p3, Eigen::Vector2d p4, std::vector<Eigen::Vector2d>& polygon, bool& isIntersecting);
double shoeLaceFormula(std::vector<Eigen::Vector2d> polygon);
bool isSectionInsideAlpha(Data2D data, int cellId, double p, std::string direction);
double interfaceFlux(Data2D& data, int cellId, std::string direction);
void preformXSweep(Data2D& data);
void resetNormalVectors(Data2D& data);


#endif