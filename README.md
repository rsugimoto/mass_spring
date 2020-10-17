# Mass Spring System

Mass spring system with libigl.
The following methods are inpremented.
- linearly-implicit time integrator
- Runge-Kutta fourth-order integrator
- Verlet integrator
- :new: the method presented in ["Fast Simulation of Mass-Spring Systems" \[Tiantian Liu et al. 2013\]](http://graphics.berkeley.edu/papers/Liu-FSM-2013-11/Liu-FSM-2013-11.pdf)


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

## References
- https://github.com/dilevin/CSC417-a2-mass-spring-3d
- https://github.com/alecjacobson/computer-graphics-mass-spring-systems
- ["Fast Simulation of Mass-Spring Systems" \[Tiantian Liu et al. 2013\]](http://graphics.berkeley.edu/papers/Liu-FSM-2013-11/Liu-FSM-2013-11.pdf)