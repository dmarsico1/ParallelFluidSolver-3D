# ParallelFluidSolver-3D
Numerically solve the Boussinesq equations in parallel on a three dimensional domain.

The domain is a rectangle that is decomposed into three dimensional columns that extend the length of the vertical domain.

To compile the code:
  make fluid_par
 
 To run the code:
  ./fluid_par -n (# of processors)
  
