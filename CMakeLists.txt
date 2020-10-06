cmake_minimum_required(VERSION 3.1)
project(main)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"         ON)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"       ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use IMGUI"          ON)

find_package(LIBIGL REQUIRED QUIET)

# Add your project files
file(GLOB SRCFILES *.cpp src/*cpp)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
add_executable(${PROJECT_NAME} ${SRCFILES})
target_link_libraries(${PROJECT_NAME} igl::core igl::opengl_glfw igl::tetgen igl::triangle igl::opengl_glfw_imgui)
