#include "bcd_simulator.h"

#include <iostream>

using Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXi, Eigen::MatrixXi, Eigen::SparseMatrix, Eigen::SparseVector, Eigen::Vector3d;

BCDSimulator::BCDSimulator(const MeshObject& obj, double dt):Simulator(obj, dt){
    int num_masses = q_curr.rows()/3;

    const MatrixXi& S = obj.S;
    std::vector<Eigen::Triplet<double>> A_entries;
    for(auto s=0; s<S.rows(); ++s) {
        auto i0 = S(s, 0); //index of mass 0
        auto i1 = S(s, 1); //index of mass 1
        for(int i=0; i<3; ++i){
            A_entries.push_back(Eigen::Triplet<double>(S.rows()*i+s, num_masses*i+i0, 1.0));
            A_entries.push_back(Eigen::Triplet<double>(S.rows()*i+s, num_masses*i+i1, -1.0));
        }
    }
    A = SparseMatrix<double>(S.rows()*3, q_curr.rows());
    A.setFromTriplets(A_entries.begin(), A_entries.end());

    const MatrixXi& _A = obj.A;
    std::vector<Eigen::Triplet<double>> C_entries;
    for(auto a=0; a<_A.rows(); ++a) {
        for(int i=0; i<3; ++i){
            C_entries.push_back(Eigen::Triplet<double>(_A.rows()*i+a, num_masses*i+_A(a), 1.0));
        }
    }
    C = SparseMatrix<double>(S.rows()*3, q_curr.rows());
    C.setFromTriplets(C_entries.begin(), C_entries.end());

    auto  Q = omega*C.transpose()*C  + obj.k*A.transpose()*A + (1.0/(dt*dt))*M;
    solver.compute(Q);
}

template <typename Derived>
void BCDSimulator::d(const Eigen::MatrixBase<Derived>& q, Eigen::VectorXd& d){
    Eigen::VectorXd Aq = A*q;
    auto q_diff = Eigen::Map<MatrixXd>((Aq).data(), Aq.rows()/3, 3);
    d = Eigen::VectorXd(Aq.rows());
    auto d_mat = Eigen::Map<MatrixXd>(d.data(), Aq.rows()/3, 3);
    d_mat = L.asDiagonal() * (q_diff.rowwise().normalized());
}

void BCDSimulator::f_ext(Eigen::VectorXd& f_ext) {
    double g = 9.8;
    f_ext = VectorXd::Zero(q_curr.rows());
    f_ext.middleRows(q_curr.rows()/3, q_curr.rows()/3) -= obj.M*g;
}

void BCDSimulator::forward_one_step() {
    auto q_rest = Eigen::Map<const VectorXd>(obj.V.data(), obj.V.rows()*3);
    VectorXd q_next = q_curr;
    VectorXd f_ext;
    this->f_ext(f_ext);
    VectorXd y = ((1.0/(dt*dt))*M*(2.0*q_curr-q_prev) + f_ext);
    VectorXd d;
    for(auto i=0; i<50; ++i){
        this->d(q_next, d);
        auto b = (obj.k*A.transpose()*d) + omega*C.transpose()*C*q_rest + y;
        q_next = solver.solve(b);
    }

    VectorXd q_dot_next = (q_next - q_prev)/(2.0*dt);

    q_prev.swap(q_curr);
    q_curr.swap(q_next);

    q_dot_prev.swap(q_dot_curr);
    q_dot_curr.swap(q_dot_next);
}
