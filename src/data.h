#ifndef DATA_H
#define DATA_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <stack>


struct Data2D;
struct Point2D;
struct Cell2D;
struct Face2D;

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> Triplet;
typedef Eigen::VectorXd Vector;

//Epsilon
#define EPS 1e-12;

//Boundary Types

#define INNERCELL 0
#define DIRICHLET 1
#define NEUMANN 2
#define SOLID 3

//Cell Navigation

#define NORTH 0
#define SOUTH 1
#define WEST 2
#define EAST 3

//Face Navigation

#define UP 0
#define DOWN 1

#define LEFT 0
#define RIGHT 1

//PISO Algorithm Steps

#define INITIAL 0
#define INTERMEDIATE_1 1
#define CORRECTED_1 2
#define INTERMEDIATE_2 3
#define CORRECTED_2 4


struct Point2D {
    int id;
    Face2D* faces[4];       //point faces
   
    double x;
    double y; 
};

struct Face2D {
    int id;
    
    Point2D* points[2];     //face points
    Cell2D* neighCells[2];  //neighbour cells
   
    double x;               //face center x-coordinate
    double y;               //face center y-coordinate
    double dx;              //face x-length
    double dy;              //face y-length

    // boundary settings
    int    bType_u;           // velocity boundary type

    // numerical coefficients
    double a_p_tilde;
    double a_e;
    double a_w;
    double a_n;
    double a_s;
    double a_p_v;         
    double b;
    
    // physical settings
    double u[5];                 // velocity
    double v[5];                 // velocity
    double u_prev;               // u(t-1)
    double v_prev;               // v(t-1)
};

struct Cell2D {
    int id;
    Face2D* faces[4];       //cell faces
    Point2D* points[4];     //cell points
    Cell2D* neighCells[4];  //neighbour cells


    double x;              //cell center x-coordinate
    double y;              //cell center y-coordinate

    double alpha;           //phase fraction
    double alpha_prev;      //alpha(t-1)

    // numeric quantities
    int     bType_sc;         // scalar boundary type
    int     bType_p;   	       // pressure boundary type

    // physical quantities
    double	 vol;           // cell volume
    double	 sc;             // scalar at cell center
    double 	 p[5];		    // pressure at cell center


    // velocity only for display (ParaView)
    double u[5];
    double v[5];

    // numerical coefficients for pressure
    double a_p;
    double a_e;
    double a_w;
    double a_n;
    double a_s;
    double b;
    
    // numerical coefficients for scalar
    double a_p_sc;
    double a_p_v_sc;
    double a_e_sc;
    double a_w_sc;
    double a_n_sc;
    double a_s_sc;
    double b_sc;


};

struct Data2D {
    Point2D* points;
    Face2D* faces;
    Cell2D* cells;


    int dimX;                                                   //grid dimensions x
    int dimY;                                                   //grid dimensions y
    int nPoints;                                                //number of points
    int nFaces;                                                 //number of faces
    int nCells;                                                 //number of cells
    int nhorizontalFaces;                                       //number of horizontal faces bzw. v velocities
    int nverticalFaces;                                         //number of vertical faces bzw. u velocities

    double maxTime;                                             //maximum time
    double dt;                                                  //time step
    int maxIteration;
    int timeStep;                                               //number of time steps 
    int pecFunc;                                                //peclet function
    int velSolver;                                              //velocity solver type
    int presSolver;                                             //pressure solver type
    int mode;                                                   //mode

    double alpha_p_relax;                                       //pressure relaxation factor
    double alpha_sc_relax;                                      //scalar relaxation factor
    double alpha_u_relax;                                       //velocity relaxation factor
    double alpha_v_relax;                                       //velocity relaxation factor

    double rho1;                                                //density phase 1
    double rho2;                                                //density phase 2
    double eta1;                                                //viscosity phase 1
    double eta2;                                                //viscosity phase 2

    std::stack <double> continuityResiduals;                    //continuity residuals
    std::stack <std::stack<double>> stackOfContinuityResiduals; //stack of continuity residuals
};

#endif