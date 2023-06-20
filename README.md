# Projet: Grands systemés lineares

  >_Developed by Patrik Rác (21214903)_

---

## Overview

The Implementation covers the exercises specified in the project of the course _Grands systemés lineares_.

Specifically, it implements a `templated LU solver` that works with any type (relevant here for values of type `Cplx`).  
It implements a `dense least square` solver based on the LU-solver.  
Lastly, it implements a `restarted GMRES` solver for complex-valued types.

The program features a test file and the main file that plots the convergence of the given problems. The program required the _matrix-archive_ folder to be in the same directory.

## Compilation

The project comes with a `CMakeLists.txt` that allows for the creation of the build system.
The steps needed to compile:

1. `mkdir build`
2. `cd build`
3. `cmake ..`
4. `make`

## Execution

The compilation created an executable `test` and an executable `proj`.  
The `test` executable runs a few test methods to ensure the method produces correct results.  
The `proj `executable runs the convergence benchmarks using the matrices in the folder _matrix-archive_ that has to be in the same folder as the _src_ and _include_ folders.  
Both can be run using

- `./test`
- `./proj`

The execution of the `proj` creates multiple output files named 

- `convergence_<i>_<restart>.txt`
  - _i_ is the matrix and right-hand side from the matrix archive
  - _restart_ is the restart value used in the benchmark

Each line of these files contains the data

- <_iteration_> <_Quadratic norm of the relative residual_>

By default, the benchmarks are executed with the _restart_ values of 20, 40, and 80.

Additionally, I provide a `plot.gnu` script that allows for the creation of plots in pdf format of the default benchmarks. For this simply exit the build folder and run

- `gnuplot plot.gnu`

Which will create the files `plot_<i>.pdf`. Pre-computed plots can be also found in the plot-archive folder.
