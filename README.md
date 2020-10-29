# Mass Spring System

Mass spring system with libigl.
The following methods are implemented.
- linearly-implicit time integrator
- Runge-Kutta fourth-order integrator
- Verlet integrator
- Method presented in ["Fast Simulation of Mass-Spring Systems" \[Tiantian Liu et al. 2013\]](http://graphics.berkeley.edu/papers/Liu-FSM-2013-11/Liu-FSM-2013-11.pdf)
- :new: Method with Discrete-time Toda Lattice Equation. (Note this method uses a different spring model)

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

Run with `-h` option for help.

## References
- https://github.com/dilevin/CSC417-a2-mass-spring-3d
- https://github.com/alecjacobson/computer-graphics-mass-spring-systems
- [Tiantian Liu, Adam W. Bargteil, James F. O'Brien, and Ladislav Kavan. 2013. Fast simulation of mass-spring systems. *ACM Trans. Graph.* 32, 6, Article 214 (November 2013), 7 pages. DOI:https://doi.org/10.1145/2508363.2508406](http://graphics.berkeley.edu/papers/Liu-FSM-2013-11/Liu-FSM-2013-11.pdf)
- Masafumi Yamamoto. 2018. *Energy-preserving Integrator using Toda Potential.* Master's thesis. The University of Tokyo, Tokyo, Japan.