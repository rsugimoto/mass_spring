#include "ltiti_simulator.h"

#include <unsupported/Eigen/KroneckerProduct>
#include <igl/edge_lengths.h>
#include <igl/colon.h>
#include <igl/setdiff.h>
#include <cmath>
#include <Eigen/SparseCholesky>

using Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXi, Eigen::MatrixXi, Eigen::SparseMatrix, Eigen::SparseVector, Eigen::Vector3d;
using Vector6d = Eigen::Matrix<double, 6, 1>;
using Matrix6d = Eigen::Matrix<double, 6, 6>;

LTITISimulator::LTITISimulator(const MeshObject& obj, double dt):obj(obj){
    this->dt = dt;
   
    q_curr = Eigen::Map<const VectorXd>(obj.V.data(), obj.V.rows()*3);
    q_prev = q_curr; //deep copy
    q_dot_curr = VectorXd::Zero(obj.V.rows()*3);
    q_dot_prev = VectorXd::Zero(obj.V.rows()*3);

    SparseMatrix<double> I3(3,3);
    I3.setIdentity();
    this->M = Eigen::kroneckerProduct(I3, (SparseMatrix<double>)obj.M.asDiagonal());

    VectorXi V_all, notA, diff_idx;
    igl::colon<int>(0, obj.V.rows()-1, V_all);
    igl::setdiff(V_all, obj.A, notA, diff_idx);
    std::vector<Eigen::Triplet<double>> P_entries;

    for(int i=0; i<notA.rows(); ++i) {
        P_entries.push_back(Eigen::Triplet<double>(i, notA(i),1.0));
        P_entries.push_back(Eigen::Triplet<double>(notA.rows() + i, obj.V.rows() + notA(i),1.0));
        P_entries.push_back(Eigen::Triplet<double>(2*notA.rows() + i, 2*obj.V.rows() + notA(i),1.0));
    }
    P = SparseMatrix<double>(notA.rows()*3, obj.V.rows()*3);
    P.setFromTriplets(P_entries.begin(), P_entries.end());

    igl::edge_lengths(obj.V, obj.S, L);
}

void LTITISimulator::set_initial_config() {
    q_curr = Eigen::Map<const VectorXd>(obj.V.data(), obj.V.rows()*3);
    q_prev = q_curr; //deep copy
    q_dot_curr = VectorXd::Zero(obj.V.rows()*3);
    q_dot_prev = VectorXd::Zero(obj.V.rows()*3);
}

void LTITISimulator::d2Vdq2(SparseMatrix<double>& d2Vdq2) {
    auto q = Eigen::Map<const MatrixXd>(q_curr.data(), q_curr.rows()/3, 3);
    const MatrixXi& S = obj.S;
    const VectorXd& K = obj.K;

    std::vector<Eigen::Triplet<double>> entries;

    //For each spring
    for(auto s=0; s<S.rows(); ++s) {
        auto q0 = q.row(S(s, 0)); //position of mass 0
        auto q1 = q.row(S(s, 1)); //position of mass 1
        double l = L(s); //rest length of the spring
        int k = K(s); //spring constant

        auto q_diff = q0 - q1;
        Vector3d q_squared = q_diff.colwise().squaredNorm();
        double q_length_squared = q_squared.sum();
        double one_over_q_length_power_3 = std::pow(q_length_squared, -1.5);
        double diag_common = k*(1.-l*q_length_squared*one_over_q_length_power_3);

        Vector6d q_diff_extended;
        q_diff_extended.topRows(3) = q0 - q1;
        q_diff_extended.bottomRows(3) = q1 - q0;

        Matrix6d _d2Vdq2 = (k*l*one_over_q_length_power_3)*q_diff_extended*(q_diff_extended.transpose());
        _d2Vdq2.diagonal().array() += diag_common;
        _d2Vdq2.topRightCorner<3, 3>().diagonal().array() -= diag_common;
        _d2Vdq2.bottomLeftCorner<3, 3>().diagonal().array() -= diag_common;

        for(int p0=0; p0<2; ++p0){
            for(int p1=0; p1<2; ++p1){
                for(int i=0; i<3; ++i){
                    for(int j=0; j<3; ++j) {
                        double val = _d2Vdq2(p0*3+i, p1*3+j);
                        entries.push_back(Eigen::Triplet<double>(q.rows()*i+S(s, p0), q.rows()*j+S(s, p1), val));
                    }
                }
            }
        }
    }

    d2Vdq2 = SparseMatrix<double>(q_curr.rows(), q_curr.rows());
    d2Vdq2.setFromTriplets(entries.begin(), entries.end());
}

void LTITISimulator::dVdq(VectorXd& dVdq) {
    auto q = Eigen::Map<const MatrixXd>(q_curr.data(), q_curr.rows()/3, 3);
    const MatrixXi& S = obj.S;
    const VectorXd& K = obj.K;

    dVdq = VectorXd::Zero(q_curr.rows());

    //For each spring
    for(auto s=0; s<S.rows(); ++s) {
        auto q0 = q.row(S(s, 0)); //position of mass 0
        auto q1 = q.row(S(s, 1)); //position of mass 1
        double l = L(s); //rest length of the spring
        int k = K(s); //spring constant

        double q_length = (q0 - q1).norm();
        double common_factor = k*(1.-l/q_length);

        auto dVdq0 = common_factor*(q0-q1);
        auto dVdq1 = common_factor*(q1-q0);

        for(int i=0; i<3; ++i) {
            dVdq(q.rows()*i +  S(s, 0)) += dVdq0(i);
            dVdq(q.rows()*i +  S(s, 1)) += dVdq1(i);
        }
    }

    //Gravity force
    double g = 9.8;
    dVdq.middleRows(q.rows(), q.rows()) += obj.M*g;
}

void LTITISimulator::forward_one_step() {  
    VectorXd dVdq;
    SparseMatrix<double> d2Vdq2;
    this->dVdq(dVdq);
    this->d2Vdq2(d2Vdq2);

    SparseMatrix<double> A = P * (M - (dt*dt) * (-d2Vdq2)) * (P.transpose());
    Eigen::VectorXd b =  P * (M * q_dot_curr + dt * (-dVdq));

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);

    q_dot_prev.swap(q_dot_curr);
    q_dot_curr = P.transpose() * solver.solve(b);

    q_prev.swap(q_curr);
    q_curr = q_prev + dt * q_dot_curr;
}

Eigen::Map<const MatrixXd> LTITISimulator::V() {
    return Eigen::Map<const MatrixXd>(q_curr.data(), q_curr.rows()/3, 3);
}

