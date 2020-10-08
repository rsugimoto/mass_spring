#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/barycenter.h>
#include <imgui/imgui.h>
#include <iostream>
#include <vector>
#include <cstring>
#include <thread>

#include "cloth.h"
#include "tetrahedral_mesh_object.h"
#include "ltiti_simulator.h"

int main(int argc, char *argv[]) {
  double dt=0.01, k=20.0, scale=1.0, decimate=0;
  char const* data_src = "cloth";
  if(argc>=2) data_src = argv[1];
  if(argc>=3) dt = std::stod(argv[2]);
  if(argc>=4) k = std::stod(argv[3]);
  if(argc>=5) scale = std::stod(argv[4]);
  if(argc>=6) decimate = std::stoi(argv[5]);

  std::cout<<"src: \""<<data_src<<"\""<<std::endl<<"dt: "<<dt<<std::endl<<"spring constant: "<<k<<std::endl
    <<"scale: "<<scale<<std::endl<<"decimate: "<<decimate<<std::endl<<std::endl;

  const MeshObject& obj = std::strcmp(data_src,"cloth")==0?
    static_cast<const MeshObject&>(Cloth(k, scale, decimate))
    : static_cast<const MeshObject&>(TetrahedralMeshObject(data_src, k, scale, decimate));
  auto simulator = LTITISimulator(obj, dt);
  std::thread([&](){while(true)simulator.forward_one_step();}).detach();

  igl::opengl::glfw::Viewer viewer;
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  viewer.plugins.push_back(&menu);
  viewer.core().is_animating = true;
  viewer.core().animation_max_fps = 30.0;
  viewer.data().set_mesh(simulator.V(), obj.F);
  viewer.data().show_lines = true;
  viewer.callback_pre_draw =  [&](igl::opengl::glfw::Viewer & viewer)->bool {
    if(viewer.core().is_animating) {
      viewer.data().set_vertices(simulator.V());
      viewer.data().compute_normals();
    }
    return false;
  };

  if (const TetrahedralMeshObject* obj_p = dynamic_cast<const TetrahedralMeshObject*>(&obj)) {
    viewer.callback_key_down = [&](igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)->bool {
      if (key >= '1' && key <= '9') {
        double t = double((key - '1')+1) / 9.0;

        Eigen::MatrixXd B;
        igl::barycenter(simulator.V(), obj_p->T, B);

        Eigen::VectorXd v = B.col(2).array() - B.col(2).minCoeff();
        v /= v.maxCoeff();

        std::vector<int> s;
        for (unsigned i=0; i<v.size();++i)
          if (v(i) <= t) s.push_back(i);

        Eigen::MatrixXi _F(s.size()*4,3);
        for (unsigned i=0; i<s.size();++i) {
          int v0 = obj_p->T(s[i], 0), v1 = obj_p->T(s[i], 1), v2 = obj_p->T(s[i], 2), v3 = obj_p->T(s[i], 3);
          _F.row(i*4+0) << v0, v1, v3;
          _F.row(i*4+1) << v0, v2, v1;
          _F.row(i*4+2) << v3, v2, v0;
          _F.row(i*4+3) << v1, v2, v3;
        }
        viewer.data().clear();
        viewer.data().set_mesh(simulator.V(), _F);
        viewer.data().set_face_based(true);
      }else if(key==' '){
        simulator.set_initial_config();
      }
      return false;
    };
    viewer.callback_key_down(viewer, '9', 0);
  }
  viewer.launch();
}
