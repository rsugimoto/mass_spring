#include "rk4_simulator.h"
#include <unsupported/Eigen/KroneckerProduct>

using Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXi, Eigen::MatrixXi, Eigen::SparseMatrix, Eigen::SparseVector, Eigen::Vector3d;

RK4Simulator::RK4Simulator(const MeshObject& obj, double dt):Simulator(obj, dt){
    SparseMatrix<double> I3(3,3);
    I3.setIdentity();
    M_inverse = Eigen::kroneckerProduct(I3, static_cast<SparseMatrix<double>>(obj.M.asDiagonal().inverse()));
}

void RK4Simulator::forward_one_step() {  
    VectorXd dVdq;
    
    this->dVdq(dVdq, q_curr);
    VectorXd dq_dot_1 = P.transpose() * P * dt*M_inverse*(-dVdq);
    VectorXd dq_1 = dt*q_dot_curr;

    this->dVdq(dVdq, q_curr+dq_1/2.0);
    VectorXd dq_dot_2 = P.transpose() * P * dt*M_inverse*(-dVdq);
    VectorXd dq_2 = dt*(q_dot_curr+dq_dot_1/2.0);

    this->dVdq(dVdq, q_curr+dq_2/2.0);
    VectorXd dq_dot_3 = P.transpose() * P * dt*M_inverse*(-dVdq);
    VectorXd dq_3 = dt*(q_dot_curr+dq_dot_2/2.0);

    this->dVdq(dVdq, q_curr+dq_3);
    VectorXd dq_dot_4 = P.transpose() * P * dt*M_inverse*(-dVdq);
    VectorXd dq_4 = dt*(q_dot_curr+dq_dot_3);

    VectorXd q_dot_next = q_dot_prev + (dq_dot_1 + 2.0*dq_dot_2 + 2.0*dq_dot_3 + dq_dot_4)/6.0;
    VectorXd q_next = q_prev + (dq_1 + 2.0*dq_2 + 2.0*dq_3 + dq_4)/6.0;
    
    q_prev.swap(q_curr);
    q_curr.swap(q_next);

    q_dot_prev.swap(q_dot_curr);
    q_dot_curr.swap(q_dot_next);
}
