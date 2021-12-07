# ParallelFluidSolver-3D
Numerically solve the Boussinesq equations in parallel on a three dimensional domain.

The domain is a rectangle that is decomposed into three dimensional columns that extend the length of the vertical domain using cartesian topologies in MPI.  The elliptic pressure equation is solved with a parallel conjugate gradient method.

To compile the code:

make fluid_par
 
To run the code:
 
mpirun -n (# of processors) fluid_par
  
