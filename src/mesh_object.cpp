#include "mesh_object.h"
#include <igl/massmatrix.h>

void MeshObject::set_mass() {
    Eigen::SparseMatrix<double> _M;
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, _M);
    M = _M.diagonal();
}

void MeshObject::init() {
    set_mesh();
    set_anchor_points();
    set_springs();
    set_mass();
}

MeshObject::~MeshObject(){}