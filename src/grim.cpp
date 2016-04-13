#include "grim.hpp"
#include "params.hpp"
#include "boundary/boundary.hpp"

/* Returns memory bandwidth in GB/sec */
double memoryBandwidth(const double numReads,
                       const double numWrites,
                       const double numEvals,
                       const double timeElapsed
                      )
{
  switch (params::dim)
  {
    case 1:
    return   (double)(params::N1 + 2*params::numGhost) 
           * 8. * (numReads + numWrites) * numEvals /timeElapsed/1e9;

    case 2:
    return   (double)(params::N1 + 2*params::numGhost) 
           * (double)(params::N2 + 2*params::numGhost)
           * 8. * (numReads + numWrites) * numEvals /timeElapsed/1e9;

    case 3:
    return   (double)(params::N1 + 2*params::numGhost) 
           * (double)(params::N2 + 2*params::numGhost)
           * (double)(params::N3 + 2*params::numGhost)
           * 8. * (numReads + numWrites) * numEvals /timeElapsed/1e9;

  }
}

void bandwidthTest()
{
//  grid prim(params::N1, params::N2, params::N3,
//            params::numGhost, params::dim, 3,
//            DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
//            DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
//            DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED
//           );
//
//  prim.vars[0] = af::randu(prim.vars[0].dims(0),
//                           prim.vars[0].dims(1),
//                           prim.vars[0].dims(2),
//                           f64
//                          );
//  prim.vars[1] = af::randu(prim.vars[1].dims(0),
//                           prim.vars[1].dims(1),
//                           prim.vars[1].dims(2),
//                           f64
//                          );
//  prim.vars[2] = af::randu(prim.vars[2].dims(0),
//                           prim.vars[2].dims(1),
//                           prim.vars[2].dims(2),
//                           f64
//                          );
//  prim.vars[2] = prim.vars[1] + prim.vars[0];
//  prim.vars[2].eval();
//  af::sync();
//
//  int numEvals = 10000;
//  af::timer::start();
//  for (int n=0; n < numEvals; n++)
//  {
//    prim.vars[2] = prim.vars[1] + prim.vars[0];
//    prim.vars[2].eval();
//  }
//  af::sync();
//  double timeElapsed = af::timer::stop();
//  PetscPrintf(PETSC_COMM_WORLD, 
//              "Summation kernel: num evals = %d, time taken = %g secs, memory bandwidth = %g GB/sec\n",
//              numEvals, timeElapsed, 
//              memoryBandwidth(2, 1, numEvals, timeElapsed)
//             );
}

int main(int argc, char **argv)
{ 
  PetscInitialize(&argc, &argv, NULL, help);
  int rank, size;
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  PetscPrintf(PETSC_COMM_WORLD, 
              "#### System info ####\n"
             );
  //if (rank==0)
  {
    //af::setDevice(1);
    af::info();
  };

  /* Local scope so that destructors of all classes are called before
   * PetscFinalize() */
  {
    timeStepper ts(params::N1, params::N2, params::N3,
                   params::dim, vars::dof, params::numGhost,
                   params::Time, params::dt,
                   params::boundaryLeft,  params::boundaryRight,
                   params::boundaryTop,   params::boundaryBottom,
                   params::boundaryFront, params::boundaryBack,
                   params::metric, params::blackHoleSpin, params::hSlope,
                   params::X1Start, params::X1End,
                   params::X2Start, params::X2End,
                   params::X3Start, params::X3End
                  );

    int numReads, numWrites;
    ts.timeStep(numReads, numWrites);

    af::sync();

    PetscPrintf(PETSC_COMM_WORLD, "\nKernel compilation complete\n");

    af::timer::start();
    int n=0;
    while(ts.time<params::finalTime)
    {
      n++;
      PetscPrintf(PETSC_COMM_WORLD, "\n|----Time step %d----|  t = %e\n", n,ts.time);
      ts.timeStep(numReads, numWrites);
    }
    double timeElapsed = af::timer::stop();
    PetscPrintf(PETSC_COMM_WORLD, "Time taken for %d time steps = %g\n",
               n, timeElapsed
		);

  }

  PetscFinalize();  
  return(0);
}
