# Mass Spring System

Mass spring system using libigl.

## Dependencies

The only dependencies are stl, eigen, [libigl](http://libigl.github.io/libigl/) and
the dependencies of the `igl::opengl::glfw::Viewer`.

The cmake build system will attempt to find libigl according to environment variables (e.g., `LIBIGL`) and searching in common desitinations (e.g., `/usr/local/libigl/`). If you haven't installed libigl before, we recommend you to clone a copy of libigl right here:

    cd libigl-example-project/
    git clone https://github.com/libigl/libigl.git

## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

## Run

From within the `build` directory just issue:

    ./main
