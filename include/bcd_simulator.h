#ifndef __BCD_SIMULATOR__
#define __BCD_SIMULATOR__

//Reference
//http://graphics.berkeley.edu/papers/Liu-FSM-2013-11/Liu-FSM-2013-11.pdf
//https://github.com/alecjacobson/computer-graphics-mass-spring-systems

#include "simulator.h"
#include <Eigen/SparseCholesky>

//Block Coordinate Descent Simulator
class BCDSimulator: public Simulator{
    private:
        Eigen::SparseMatrix<double> A; //signed incidence matrix
        Eigen::SparseMatrix<double> C; //pinned vertices
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

        template <typename Derived> void d(const Eigen::MatrixBase<Derived>& q, Eigen::VectorXd& d);
        
        double omega = 100.0; //strength of penalty for pinned vertices
    public:
        BCDSimulator(const MeshObject& obj, double k, double dt);
        void step();
};

#endif //__BCD_SIMULATOR___