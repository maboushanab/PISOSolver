#ifndef ADVECTION_H
#define ADVECTION_H

#include "data.h"

void advectAlpha(Data2D& data);
void handleExcessAlpha(Data2D& data);
void preformYSweep(Data2D& data);
double calculateFluxLength(Data2D data, int cellId, std::string direction);
void calculatePolygonArea(double n, double m, Eigen::Vector2d normalVec, Eigen::Vector2d p1, Eigen::Vector2d p2, Eigen::Vector2d p3, Eigen::Vector2d p4, std::vector<Eigen::Vector2d>& polygon, bool& isIntersecting, bool xCoordinates);
void determinePolygonAsCell(Data2D& data, int cellId, std::vector<Eigen::Vector2d>& polygon, double m, double n);
double shoeLaceFormula(std::vector<Eigen::Vector2d> polygon);
bool isSectionInsideAlpha(double p, Eigen::Vector2d interfaceMidPoint, Eigen::Vector2d normalVec, std::string direction);
double interfaceFlux(Data2D& data, int cellId, std::string direction);
void preformXSweep(Data2D& data);


#endif