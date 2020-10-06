#include "tetrahedral_mesh_object.h"
#include <igl/readOFF.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/massmatrix.h>
#include <igl/edges.h>
#include <igl/sort.h>
#include <vector>
#include <array>
#include <cmath>
#include <iostream>

using Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXi, Eigen::MatrixXi;

TetrahedralMeshObject::TetrahedralMeshObject(char const* file, double k, double size) {
  this->file = file;
  this->k = k;
  this->size = size;
  init();
}

#include <igl/isdiag.h>
void TetrahedralMeshObject::set_massmatrix() {
    Eigen::SparseMatrix<double> _M;
    igl::massmatrix(V, T, igl::MASSMATRIX_TYPE_BARYCENTRIC, _M);
    M = _M.diagonal();
}

void TetrahedralMeshObject::set_mesh() {
  MatrixXd _V; 
  MatrixXi _F;

  igl::readOFF(this->file, _V, _F);
  _V *= size;

  igl::copyleft::tetgen::tetrahedralize(_V, _F, "pq1.414a0.1", V, T, F);

  std::cout<<"x: "<<V.col(0).minCoeff()<<" "<<V.col(0).maxCoeff()<<std::endl
           <<"y: "<<V.col(1).minCoeff()<<" "<<V.col(1).maxCoeff()<<std::endl
           <<"z: "<<V.col(2).minCoeff()<<" "<<V.col(2).maxCoeff()<<std::endl;
}

void TetrahedralMeshObject::set_anchor_points() {
  A = VectorXi((int)(V.rows()*0.1));

  VectorXd sorted;
  VectorXi sorted_indices;
  igl::sort(V.col(1), 1, false, sorted, sorted_indices);

  for(auto i=0; i<A.rows(); ++i) A(i) = sorted_indices(i);
}

void TetrahedralMeshObject::set_springs() {
  igl::edges(F, S);
  K = VectorXd::Ones(S.rows())*k;
}