#include "liti_simulator.h"

using Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXi, Eigen::MatrixXi, Eigen::SparseMatrix, Eigen::SparseVector, Eigen::Vector3d;

LITISimulator::LITISimulator(const MeshObject& obj, double dt):Simulator(obj, dt){

}

void LITISimulator::forward_one_step() {  
    VectorXd dVdq;
    SparseMatrix<double> d2Vdq2;
    this->dVdq(dVdq, q_curr);
    this->d2Vdq2(d2Vdq2, q_curr);

    SparseMatrix<double> A = P * (M - (dt*dt) * (-d2Vdq2)) * P.transpose();
    Eigen::VectorXd b =  P * (M * q_dot_curr + dt * (-dVdq));

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);

    VectorXd q_dot_next = P.transpose() * solver.solve(b);
    VectorXd q_next = q_prev + dt * q_dot_next;

    q_prev.swap(q_curr);
    q_curr.swap(q_next);

    q_dot_prev.swap(q_dot_curr);
    q_dot_curr.swap(q_dot_next);
}
