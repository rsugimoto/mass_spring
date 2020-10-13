# Mass Spring System

Mass spring system with libigl.
The following integrators are inpremented.
- linearly-implicit time integrator
- Runge-Kutta fourth-order integrator
- Verlet integrator


## Dependencies
- libigl
    
    Run the following command in the project folder.
    ```
    git clone https://github.com/libigl/libigl.git
    ```

## Compile
    
    mkdir build
    cd build
    cmake ..
    make

## Run

From within the `build` directory just issue:

    ./main

Read main.cpp for command line options.