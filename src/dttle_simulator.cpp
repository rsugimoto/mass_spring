#include "dttle_simulator.h"

using Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXi, Eigen::MatrixXi, Eigen::SparseMatrix, Eigen::SparseVector, Eigen::Vector3d;

DTTLESimulator::DTTLESimulator(const MeshObject& obj, double a, double b, double dt):
    Simulator(obj, 0.0, dt), a(a), b(b)
{

}

void DTTLESimulator::step() {
    const MatrixXi& S = obj.S;
    const VectorXd& M = obj.M;
    auto q = Eigen::Map<const MatrixXd>(q_curr.data(), q_curr.rows()/3, 3);

    MatrixXd update = MatrixXd::Zero(q.rows(), 3);
    for(auto s=0; s<S.rows(); ++s) {
        int i0 = S(s, 0); //index of mass 0
        int i1 = S(s, 1); //index of mass 1
        auto q0 = q.row(i0); //position of mass 0
        auto q1 = q.row(i1); //position of mass 1
        double l0 = L(s); //rest length of the spring
        double l = (q0-q1).norm(); //current length of the spring
        double exp_b_l0 = std::exp(b*l0);
        double ds_0_squared = dt*dt*(a*b*exp_b_l0)/M(i0);
        double ds_1_squared = dt*dt*(a*b*exp_b_l0)/M(i1);
        double exp_b_l = std::exp(b*l);
        
        double update_0 = std::log((1.0 + ds_0_squared/exp_b_l)/(1.0 + ds_0_squared/exp_b_l0))/b;
        double update_1 = std::log((1.0 + ds_1_squared/exp_b_l)/(1.0 + ds_1_squared/exp_b_l0))/b;

        auto n = (q0-q1).normalized();
        update.row(i0) += update_0 * n;
        update.row(i1) -= update_1 * n;
    }

    VectorXd ext_force_term;
    this->f_ext(ext_force_term);
    ext_force_term = dt*dt*M_inverse*ext_force_term;

    VectorXd q_next = q_curr + P.transpose() * P * (q_curr - q_prev + Eigen::Map<const VectorXd>(update.data(), update.rows()*3) + ext_force_term);
    VectorXd q_dot_next = (q_next - q_prev)/(2.0*dt);

    q_prev.swap(q_curr);
    q_curr.swap(q_next);

    q_dot_prev.swap(q_dot_curr);
    q_dot_curr.swap(q_dot_next);
}
