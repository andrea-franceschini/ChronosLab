# ChronosLab

ChronosLab is a MATLAB package implementing **Classical AMG**. It leverages on
C/C++ source codes called from MATLAB through a MEX interface to speedup the
more time consuming kernels.

## Description ##

Multigrid methods are considered scalable or optimal because they can
efficiently solve a linear system with *N* unknowns using *O(N)* computational
effort. This optimality is achieved through the use of two complementary
processes: *smoothing* and *coarse-grid correction*. Smoothing involves the
application of a smoother, such as Gauss-Seidel, which is a simple iterative
method. Coarse-grid correction involves transferring information to a coarser
grid through restriction, solving the system of equations on the coarse grid,
and then transferring the solution back to the finer grid through
interpolation. Generally speaking, smoothing reduces high-frequency errors,
while coarse-grid correction eliminates low-frequency errors.

The key ingredients of a classical AMG are:

* **smoothing**: a simple iterative method to dump high-frequency errors;
* **test space computation**: an approximation of the near null kernel, i.e.,
  the set of vectors associated to low-frequency errors;
* **coarsening**: the process of identifying which fine nodes become coarse;
* **prolongation**: the operator transferring information from a coarse to a
  fine grid.

## Package content ##

The package is organized as follows. In the **source** folder, there are:

* **amg**: it contains the main recursive kernel;
* **coarsen**: with all the possible coarsening strategies;
* **prolong**: it allows for the computation of the interpolation operator;
* **smoother**: with several effective smoothers;
* **tspace**: it collects all the functionalities to approximate the near null
  space;
* **utils**: with some generic functions;
* **driver**: a simple driver.

Moreover, there is the *compileAll.m* file, useful to create all the MEX
objects. It has to be run just once.

In the **example** folder, two test cases are provided. They are:

* **mech**: a simple structural problem;
* **lapl**: a simple flow problem.

This folder contains also *initChronosLab.m*, a file to setup the path to run
ChronosLab. It has to be called each time MATLAB is opened.

## Authors ##
The main contributors to this package are:

* **Carlo Janna** (carlo.janna@unipd.it and c.janna@m3eweb.it), associate
  professor at [University of Padova](https://www.unipd.it/en) and president
and chief scientific officer of [M3E S.r.l.](https://www.m3eweb.it) (via
Giambellino 7, 35129 Padova, Italy);
* **Andrea Franceschini** (andrea.franceschini@unipd.it), assistant professor
  at [University of Padova](https://www.unipd.it/en).

## Parameters ##
The code uses the following parameters:

* problem ...

## Terms of usage ##
**ChronosLab** is an open source software under the MIT license and can be
freely downloaded and modified, provided the changes are distributed under the
same license. See the details in [LICENSE](LICENSE).

If you use **ChronosLab** in your research, please cite it. A few overview
references are provided here for your convenience:

* C. Janna, A. Franceschini, J. B. Schroder, L. Olson (2023). [Parallel
  Energy-Minimization Prolongation for Algebraic
Multigrid](https://doi.org/10.1137/22M1513794). SIAM Journal on Scientific
Computing 45 (5), A2561-A2584.
* G. Isotton, M. Frigo, N. Spiezia, C. Janna (2021). [Chronos: a general
  purpose classical AMG solver for high performance
computing](https://doi.org/10.1137/21M1398586). SIAM Journal on Scientific
Computing 43 (5), C335-C357.
* A. Franceschini, V. A. P. Magri, G. Mazzucco, N. Spiezia, C. Janna (2019). [A
  robust adaptive algebraic multigrid linear solver for structural
mechanics](https://doi.org/10.1016/j.cma.2019.04.034). Computer Methods in
Applied Mechanics and Engineering 352, 389-416.
* V. A. Paludetto Magri, A. Franceschini, C. Janna (2019). [A novel algebraic
  multigrid approach based on adaptive smoothing and prolongation for
ill-conditioned systems](https://doi.org/10.1137/17M1161178). SIAM Journal on
Scientific Computing 41 (1), A190-A219.

Finally, please note that this software is provided "as is", without any
guarantees or warranties of any kind. The user assumes all risk for any damage
or data loss resulting from the use of this software.

## Bugs ##
Please, report the bugs either opening an issue on Github or emailing the
details to one of the authors.
