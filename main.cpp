#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/barycenter.h>
#include <imgui/imgui.h>
#include "cloth.h"
#include "tetrahedral_mesh_object.h"
#include <iostream>
#include <cmath>
#include <vector>

using Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXi, Eigen::MatrixXi, Eigen::SparseMatrix, Eigen::SparseVector;

void forward_one_step(MatrixXd& V_curr, MatrixXd& V_prev, const MatrixXd& V_init, const VectorXi& A, const VectorXd& M, const MatrixXi& S, const VectorXd& K, const VectorXd& RL, double t=0, double dt=0.001, double damping_coeff=1.0) {
  MatrixXd forces = MatrixXd::Zero(V_curr.rows(), 3);

  for (auto s=0; s<S.rows(); ++s) {
    auto v_i = S(s, 0), v_j = S(s, 1);
    auto x_i = V_curr.row(v_i), x_j = V_curr.row(v_j);
    auto x = x_i - x_j;
    auto x_norm = x.norm();
    auto f = - K(s) * (x_norm-RL(s)) * x/x_norm;
    forces.row(v_i) += f;
    forces.row(v_j) -= f;
  }

  double sin=std::sin(t)+std::cos(3*t);
  for (auto i=0; i<V_curr.rows(); ++i) {
    double m = M(i);
    forces(i, 1) -= 9.8*m; //gravitational force
    forces(i, 2) += m*sin; //wind
    forces.row(i) -=  damping_coeff * (V_curr.row(i) - V_prev.row(i)); //damping force
  }

  MatrixXd V_next = (M.asDiagonal().inverse())*forces*dt*dt + 2*V_curr - V_prev;
  V_prev.swap(V_curr);
  V_curr.swap(V_next);

  //Fix anchor points
  for (auto i=0; i<A.rows(); ++i) V_curr.row(A(i)) = V_init.row(A(i));
}

int main(int argc, char *argv[]) {
  double dt=0.001, damping_coeff=1.0, k=10.0, scale=1.0;
  char const* data_src = "../data/bunny.off";
  if(argc>=2) data_src = argv[1];
  if(argc>=3) dt = std::stod(argv[2]);
  if(argc>=5) k = std::stod(argv[3]);
  if(argc>=4) damping_coeff = std::stod(argv[4]);
  if(argc>=5) scale = std::stod(argv[5]);

  std::cout<<"src: \""<<data_src<<"\""<<std::endl<<"dt: "<<dt<<std::endl<<"spring constant: "<<k<<std::endl
    <<"damping coeff: "<<damping_coeff<<std::endl<<"scale: "<<scale<<std::endl<<std::endl;

  // auto obj = Cloth(k);
  auto obj = TetrahedralMeshObject(data_src, k, scale);
  MatrixXd V_curr = obj.V;
  MatrixXd V_prev = obj.V;

  double t = 0.0;
  auto update = [&](igl::opengl::glfw::Viewer & viewer)->bool {
    if(viewer.core().is_animating) {
      for (auto i=0; i<100; i++) {
        t += dt;
        forward_one_step(V_curr, V_prev, obj.V, obj.A, obj.M, obj.S, obj.K, obj.RL, t, dt, damping_coeff);
      }
      viewer.data().set_vertices(V_curr);
      viewer.data().compute_normals();
    }
    return false;
  };

  igl::opengl::glfw::Viewer viewer;
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  viewer.plugins.push_back(&menu);
  viewer.core().is_animating = true;
  viewer.core().animation_max_fps = 30.0;
  viewer.data().set_mesh(V_curr, obj.F);
  viewer.data().show_lines = true;
  viewer.callback_pre_draw = update;
  viewer.callback_key_down = [&](igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)->bool {
    if (key >= '1' && key <= '9') {
      double t = double((key - '1')+1) / 9.0;

      MatrixXd B;
      igl::barycenter(V_curr, obj.T, B);

      VectorXd v = B.col(2).array() - B.col(2).minCoeff();
      v /= v.maxCoeff();

      std::vector<int> s;
      for (unsigned i=0; i<v.size();++i)
        if (v(i) < t) s.push_back(i);

      MatrixXi _F(s.size()*4,3);
      for (unsigned i=0; i<s.size();++i) {
        int v0 = obj.T(s[i], 0), v1 = obj.T(s[i], 1), v2 = obj.T(s[i], 2), v3 = obj.T(s[i], 3);
        _F.row(i*4+0) << v0, v1, v3;
        _F.row(i*4+1) << v0, v2, v1;
        _F.row(i*4+2) << v3, v2, v0;
        _F.row(i*4+3) << v1, v2, v3;
      }
      viewer.data().clear();
      viewer.data().set_mesh(V_curr, _F);
      viewer.data().set_face_based(true);
    }
    return false;
  };
  viewer.callback_key_down(viewer, '9', 0);
  viewer.launch();
}
