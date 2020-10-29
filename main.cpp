#include <chrono>
#include <cstring>
#include <iostream>
#include <memory>
#include <string>
#include <thread>
#include <vector>

#include <cxxopts.hpp>
#include <igl/barycenter.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <imgui/imgui.h>

#include "cloth.h"
#include "tetrahedral_mesh_object.h"

#include "bcd_simulator.h"
#include "dttle_simulator.h"
#include "liti_simulator.h"
#include "rk4_simulator.h"
#include "verlet_simulator.h"

auto parse_args(int argc, char *argv[]){
  using namespace cxxopts;
  Options options{"Mass spring system simulator"};
  options.add_options()
    ("m,simulator", "type of simulator. BCD, DTTLE, LITI, RK4, or Verlet", value<std::string>()->default_value("BCD"))
    ("s,src", "input mesh in off format, or \"cloth\" for default cloth mesh", value<std::string>()->default_value("cloth"))
    ("c,scale", "scale of geometry", value<double>()->default_value("1.0"))
    ("t,dt", "time step of simulation", value<double>()->default_value("0.01"))
    ("k", "spring constant (not used for DTTLE simulator)", value<double>()->default_value("20.0"))
    ("a", "constant a for DTTLE simulator", value<double>()->default_value("-5.0"))
    ("b", "constant b for DTTLE simulator", value<double>()->default_value("-0.1"))
    ("d,decimate", "when a non-zero integer is given, decimate the mesh to the given number of faces", value<int>()->default_value("0"))
    ("h,help", "Print Usage");
  auto result = options.parse(argc, argv);
  if(result.count("help")) {std::cout<<options.help()<<std::endl; exit(0);}
  
  std::cout<<"simulator: "<<result["simulator"].as<std::string>()
  <<std::endl<<"src: \""<<result["src"].as<std::string>()<<"\""
  <<std::endl<<"scale: "<<result["scale"].as<double>()
  <<std::endl<<"dt: "<<result["dt"].as<double>()
  <<std::endl<<"k: "<<result["k"].as<double>()
  <<std::endl<<"a: "<<result["a"].as<double>()
  <<std::endl<<"b: "<<result["b"].as<double>()
  <<std::endl<<"decimate: "<<result["decimate"].as<int>()
  <<std::endl<<std::endl;
  return result;
}

int main(int argc, char *argv[]) {
  auto options = parse_args(argc, argv);

  const MeshObject& obj = options["src"].as<std::string>() == "cloth"?
    static_cast<const MeshObject&>(Cloth(options["scale"].as<double>(), options["decimate"].as<int>()))
    : static_cast<const MeshObject&>(TetrahedralMeshObject(options["src"].as<std::string>().c_str(), options["scale"].as<double>(), options["decimate"].as<int>()));

  std::unique_ptr<Simulator> simulator;
  
  std::string sim = options["simulator"].as<std::string>();
  if (sim == "BCD") simulator.reset(new BCDSimulator(obj, options["k"].as<double>(), options["dt"].as<double>()));
  else if (sim == "DTTLE") simulator.reset(new DTTLESimulator(obj, options["a"].as<double>(), options["b"].as<double>(), options["dt"].as<double>()));
  else if (sim == "LITI") simulator.reset(new LITISimulator(obj, options["k"].as<double>(), options["dt"].as<double>()));
  else if (sim == "RK4") simulator.reset(new RK4Simulator(obj, options["k"].as<double>(), options["dt"].as<double>())); 
  else if (sim == "Verlet") simulator.reset(new VerletSimulator(obj, options["k"].as<double>(), options["dt"].as<double>()));
  else {
    std::cout<<"No simulator option\""<<options["simulator"].as<std::string>()<<"\". Using default BCD simulator."<<std::endl;
    simulator.reset(new BCDSimulator(obj, options["k"].as<double>(), options["dt"].as<double>()));
  }

  std::mutex mtx;
  volatile bool terminate_thread = false;
  auto simulation_thread = std::thread([&](){
    auto time = std::chrono::system_clock::now();
    while(!terminate_thread){
      auto new_time = std::chrono::system_clock::now();
      if (static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(new_time - time).count() / 1000.0) < options["dt"].as<double>()) continue;
      else time = new_time;

      std::lock_guard lock(mtx);
      simulator->step();
    }
  });

  igl::opengl::glfw::Viewer viewer;
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  viewer.plugins.push_back(&menu);
  viewer.core().is_animating = true;
  viewer.core().animation_max_fps = 30.0;
  viewer.data().set_mesh(simulator->V(), obj.F);
  viewer.data().show_lines = true;
  viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer & viewer)->bool {
    if(viewer.core().is_animating) {
      viewer.data().set_vertices(simulator->V());
      viewer.data().compute_normals();
    }
    return false;
  };

  viewer.callback_key_down = [&](igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)->bool {
    if (key >= '1' && key <= '9') {
      //Visualization of internal tetrahedra
      if (const TetrahedralMeshObject* obj_p = dynamic_cast<const TetrahedralMeshObject*>(&obj)){
        double t = double((key - '1')+1) / 9.0;

        Eigen::MatrixXd B;
        igl::barycenter(simulator->V(), obj_p->T, B);

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
        viewer.data().set_mesh(simulator->V(), _F);
        viewer.data().set_face_based(true);
      }
    }else if(key==' '){
      std::lock_guard lock(mtx);
      simulator->set_initial_config();
    }
    return false;
  };
  viewer.callback_key_down(viewer, '9', 0);
  viewer.launch();

  terminate_thread = true;
  simulation_thread.join();
}
