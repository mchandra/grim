#include "timestepper.hpp"

void timeStepper::solve(grid &primGuess)
{
  /* Get the domain of the bulk */
  af::seq domainX1 = *residual->domainX1;
  af::seq domainX2 = *residual->domainX2;
  af::seq domainX3 = *residual->domainX3;

  int world_rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(PETSC_COMM_WORLD, &world_size);

  for (int nonLinearIter=0;
       nonLinearIter < params::maxNonLinearIter; nonLinearIter++
      )
  {
    /* True residual, with explicit terms (not needed for Jacobian) */
    int numReadsResidual, numWritesResidual;
    computeResidual(primGuess, *residual, true, 
                    numReadsResidual, numWritesResidual
                   );
    for (int var=0; var < residual->numVars; var++)
    {
      /* Need residualSoA to compute norms */
      residualSoA(span, span, span, var) = residual->vars[var];

      /* Initialize primGuessPlusEps. Needed to numerically assemble the
       * Jacobian */
      primGuessPlusEps->vars[var] = primGuess.vars[var];
    }

    /* Sum along last dim:vars to get L2 norm */
    array l2Norm  = 
      af::sum(af::pow(residualSoA(domainX1, domainX2, domainX3), 2.), 3);
    l2Norm.eval();
    array notConverged      = l2Norm > params::nonlinearsolve_atol;
    array conditionIndices  = where(notConverged > 0);

    /* Communicate residual */
    double localresnorm = 
      af::norm(af::flat(residualSoA(domainX1, domainX2, domainX3)),AF_NORM_VECTOR_1);
    double globalresnorm = localresnorm;
    int localNonConverged = conditionIndices.elements();
    int globalNonConverged = localNonConverged;
    if (world_rank == 0)
      {
	double temp;
	int Nel;
	for(int i=1;i<world_size;i++)
	  {
	    MPI_Recv(&temp, 1, MPI_DOUBLE, i, i, PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
	    MPI_Recv(&Nel, 1, MPI_INT, i, i+world_size, PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
	    globalresnorm+=temp;
	    globalNonConverged+=Nel;
	  }
      }
    else
      {
	MPI_Send(&localresnorm, 1, MPI_DOUBLE, 0, world_rank, PETSC_COMM_WORLD);
	MPI_Send(&localNonConverged, 1, MPI_INT, 0, world_rank+world_size, PETSC_COMM_WORLD);
      }
    MPI_Barrier(PETSC_COMM_WORLD);
    MPI_Bcast(&globalresnorm,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
    MPI_Barrier(PETSC_COMM_WORLD);
    MPI_Bcast(&globalNonConverged,1,MPI_INT,0,PETSC_COMM_WORLD);
    MPI_Barrier(PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, " ||Residual|| = %g; %i pts haven't converged\n", 
                globalresnorm,globalNonConverged
               );
    if (globalNonConverged == 0)
    {
      break;
    }


    /* Residual without explicit terms, for faster Jacobian assembly */
    computeResidual(primGuess, *residual, false,
                    numReadsResidual, numWritesResidual
                   );

    /* Assemble the Jacobian in Struct of Arrays format where the physics
     * operations are all vectorized */
    for (int row=0; row < numVars; row++)
    {
      /* Recommended value of Jacobian differencing parameter to achieve fp64
       * machine precision */
      double epsilon = params::JacobianAssembleEpsilon;

      array smallPrim = af::abs(primGuess.vars[row])<.5*epsilon;

      primGuessPlusEps->vars[row]  = 
          (1. + epsilon)*primGuess.vars[row]*(1.-smallPrim)
	      + smallPrim*epsilon; 

      computeResidual(*primGuessPlusEps, *residualPlusEps, false,
                      numReadsResidual, numWritesResidual
                     );

      for (int column=0; column < numVars; column++)
      {
        jacobianSoA(span, span, span, column + numVars*row)
          = (  residualPlusEps->vars[column] 
             - residual->vars[column]
            )
            /(primGuessPlusEps->vars[row]-primGuess.vars[row]);
      }
      /* reset */
      primGuessPlusEps->vars[row]  = primGuess.vars[row]; 
    }
    /* Jacobian assembly complete */

    /* Solve the linear system Jacobian * deltaPrim = -residual for the
     * correction deltaPrim */

    array jacobianAoS = af::reorder(jacobianSoA, 3, 0, 1, 2);

    /* RHS of Ax = b in Array of Structs format */
    array bAoS = -af::reorder(residualSoA, 3, 0, 1, 2);


    /* Now solve Ax = b using direct inversion, where
     * A = Jacobian
     * x = deltaPrim
     * b = -residual 
     *
     * Currently inverting locally by looping over individual zones. Need to
     * call the batch function magma_dgesv_batched() from the MAGMA library
     * for optimal use on NVIDIA cards */
    batchLinearSolve(jacobianAoS, bAoS, deltaPrimAoS);

    /* Done with the solve. Now rearrange from AoS -> SoA */
    array deltaPrimSoA = af::reorder(deltaPrimAoS, 1, 2, 3, 0);

    /* Quadratic backtracking :
     We minimize f(u+stepLength*du) = 0.5*sqr(residual[u+stepLength*du]).
     We use
       f0 = f(u)
       fPrime0 = df/d(stepLength)(u) = -2*f0
       [because stepLength=1 is the solution of the linearized problem...]
       f1 = f(u+stepLength0*du)
     to reconstruct
       f = (f1-f0-fPrime0*stepLength0)*(stepLength/stepLength0)^2 + fPrime0*stepLength + f0
     which has a minimum at the new value of stepLength,
       stepLength = -fPrime0*stepLength0^2 / (f1-f0-fPrime0*stepLength0)/2
     */
    array f0      = 0.5 * l2Norm;
    array fPrime0 = -2.*f0;
    
    /* Start with a full step */
    stepLength = 1.;
    int lineSearchIter=0;
    for (;lineSearchIter < params::maxLineSearchIters; lineSearchIter++
        )
    {
      /* 1) First take current step stepLength */
      for (int var=0; var<numVars; var++)
      {
        primGuessLineSearchTrial->vars[var] =  
          primGuess.vars[var] + stepLength*deltaPrimSoA(span, span, span, var);
      } 

      /* ...and then compute the norm */
      computeResidual(*primGuessLineSearchTrial, *residual, true,
                      numReadsResidual, numWritesResidual
                     );
      for (int var=0; var<numVars; var++)
      {
        residualSoA(span, span, span, var) = residual->vars[var];
      }
      l2Norm = af::sum(af::pow(residualSoA(domainX1, domainX2, domainX3), 2.), 3);
      array f1 = 0.5 * l2Norm;

      /* We have 3 pieces of information:
       * a) f(0)
       * b) f'(0) 
       * c) f(stepLength) 
       */
    
      const double alpha    = 1e-4;
      const double EPS      = params::linesearchfloor;
      array stepLengthNoGhost = stepLength(domainX1, domainX2, domainX3);
      array condition = f1 > (f0*(1. - alpha*stepLengthNoGhost) +EPS);
      array denom     =   (f1-f0-fPrime0*stepLengthNoGhost) * condition 
                        + (1.-condition);
      array nextStepLengthNoGhost =
        -fPrime0*stepLengthNoGhost*stepLengthNoGhost/denom/2.;
      stepLength(domainX1, domainX2, domainX3)
        = stepLengthNoGhost*(1. - condition) + condition*nextStepLengthNoGhost;
      
      array conditionIndices = where(condition > 0);
      if (conditionIndices.elements() == 0)
      {
        break;
      }
    }

    /* stepLength has now been set */
    for (int var=0; var<numVars; var++)
    {
      primGuess.vars[var] = 
        primGuess.vars[var] + stepLength*deltaPrimSoA(span, span, span, var);
    }
  }

  
  /* Diagnostics */
  /*{
    for (int var=0; var<vars::dof; var++)
      {
	array res = af::abs(residual->vars[var]((domainX1, domainX2, domainX3)));
	array badPoints = where(res>1.);
	if(badPoints.elements()>0)
	  {
	    printf("Found bad residual for variable %i on proc %i\n",var,world_rank);
	    af_print(badPoints);
	    array flatR = af::flat(XCoords->vars[0](domainX1, domainX2, domainX3));
	    af_print(af::exp(flatR(badPoints)),12);
	    array flatRho = af::flat(primGuess.vars[vars::RHO](domainX1, domainX2, domainX3));
	    array flatU = af::flat(primGuess.vars[vars::U](domainX1, domainX2, domainX3));
	    array flatDP = af::flat(primGuess.vars[vars::DP](domainX1, domainX2, domainX3));
	    af_print(flatRho(badPoints),12);
	    af_print(flatU(badPoints),12);
	    af_print(flatDP(badPoints),12);
	    array flatRes = af::flat(res);
	    af_print(flatRes(badPoints),12);
	    exit(1);
	  }
      }
      }*/

}

void timeStepper::batchLinearSolve(const array &A, const array &b, array &x)
{
  A.host(AHostPtr);
  b.host(bHostPtr);
  x.host(xHostPtr);

  int numVars = residual->numVars;
  int N1Total = residual->N1Total;
  int N2Total = residual->N2Total;
  int N3Total = residual->N3Total;

  #pragma omp parallel for
  for (int k=0; k<N3Total; k++)
  {
    for (int j=0; j<N2Total; j++)
    {
      for (int i=0; i<N1Total; i++)
      {
        double ALocal[numVars*numVars];
        double bLocal[numVars];
        int pivot[numVars];

        const int spatialIndex = 
          i +  N1Total*(j + (N2Total*k) );

        /* Assemble ALocal */
        for (int row=0; row < numVars; row++)
        {
          for (int column=0; column < numVars; column++)
          {
            const int indexALocal = column + (numVars*row);
            const int indexAHost  = 
              column + numVars*(row + (numVars*spatialIndex) );

            ALocal[indexALocal] = AHostPtr[indexAHost];
          }
        }

        /* Assemble bLocal */
        for (int column=0; column < numVars; column++)
        {
          const int indexbLocal = column;
          const int indexbHost  = column + numVars*spatialIndex;

          bLocal[indexbLocal] = bHostPtr[indexbHost];
        }

        LAPACKE_dgesv(LAPACK_COL_MAJOR, numVars, 1, ALocal, numVars, 
                      pivot, bLocal, numVars
                     );

        /* Copy solution to xHost */
        for (int column=0; column < numVars; column++)
        {
          const int indexbLocal = column;
          const int indexbHost  = column + numVars*spatialIndex;

          xHostPtr[indexbHost] = bLocal[indexbLocal];
        }

      }
    }
  }

  /* Copy solution to x on device */
  x = array(numVars, N1Total, N2Total, N3Total, xHostPtr);
}
