#include "cloth.h"
#include <igl/triangle/triangulate.h>
#include <igl/adjacency_list.h>
#include <igl/unique.h>
#include <vector>
#include <array>
#include <cmath>

using Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXi, Eigen::MatrixXi;

Cloth::Cloth(double k, double size) {
  this->k = k;
  this->size = size;
  init();
}

void Cloth::set_mesh() {
  auto _V =  (MatrixXd(4, 2) <<
    -size/2, -size/2,
    size/2, -size/2,
    size/2, size/2,
    -size/2, size/2
  ).finished();

  auto _E =  (MatrixXi(4, 2) <<
    0, 1, 1, 2, 2, 3, 3, 0
  ).finished();

  auto _H = MatrixXd(0,2);

  MatrixXd V_;
  igl::triangle::triangulate(_V,_E,_H,"a0.0005q",V_,F);

  V = MatrixXd::Zero(V_.rows(), 3);
  V.leftCols(2) = V_;

  for(auto i=0; i<V.rows(); ++i) {
    auto a = size/2 - V(i, 1);
    V(i, 2) += a*std::sqrt(3)/2.0;
    V(i, 1) += a/2.0;
  }
}

void Cloth::set_anchor_points() {
  A = VectorXi(2);

  Eigen::RowVector3d p0; p0 <<-size/2, size/2, 0;
  Eigen::RowVector3d p1; p1 << size/2, size/2, 0;

  (V.rowwise() - p0).rowwise().squaredNorm().minCoeff(&A(0));
  (V.rowwise() - p1).rowwise().squaredNorm().minCoeff(&A(1));
}

void Cloth::set_springs() {
  std::vector<std::vector<int>> A;
  igl::adjacency_list(F, A);

  std::vector<std::array<int, 2>> springs;

  for(auto i = 0; i<V.rows(); ++i){
    for (int p1: A[i]){
      if(i<p1)springs.push_back(std::array<int, 2>{i, p1});
      for (int p2: A[p1]){
        if(i<p2)springs.push_back(std::array<int, 2>{i, p2});
      }
    }
  }

  S = MatrixXi(springs.size(),2);
  K = VectorXd(springs.size());
  for(auto i=0; i<springs.size(); ++i){
    S(i, 0) = springs[i][0];
    S(i, 1) = springs[i][1];
    K(i) = k;
  }
}