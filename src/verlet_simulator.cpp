#include "verlet_simulator.h"
#include <unsupported/Eigen/KroneckerProduct>

using Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXi, Eigen::MatrixXi, Eigen::SparseMatrix, Eigen::SparseVector, Eigen::Vector3d;

VerletSimulator::VerletSimulator(const MeshObject& obj, double k, double dt):Simulator(obj, k, dt){
    VectorXd dVdq;
    this->dVdq(dVdq, q_curr);
    q_curr = q_prev + P.transpose() * P *  (dt*dt*0.5*M_inverse*(-dVdq));

    q_dot_curr = (q_curr - q_prev)/dt;
}

 void VerletSimulator::set_initial_config(){
    Simulator::set_initial_config();
    VectorXd dVdq;
    this->dVdq(dVdq, q_curr);
    q_curr = q_prev + P.transpose() * P *  (dt*dt*0.5*M_inverse*(-dVdq));
    q_dot_curr = (q_curr - q_prev)/dt;
 }

void VerletSimulator::step() {  
    VectorXd dVdq;
    this->dVdq(dVdq, q_curr);
    VectorXd q_next = q_curr + P.transpose() * P * (q_curr - q_prev + dt*dt*M_inverse*(-dVdq));
    // VectorXd q_dot_next = (q_next - q_curr)/dt;
    VectorXd q_dot_next = (q_next - q_prev)/(2.0*dt);

    q_prev.swap(q_curr);
    q_curr.swap(q_next);

    q_dot_prev.swap(q_dot_curr);
    q_dot_curr.swap(q_dot_next);
}
