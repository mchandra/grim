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

  array l2Norm,notConverged;
  int localNonConverged,globalNonConverged;
  double localresnorm,globalresnorm;
  notConverged = GZmask;

  for (int nonLinearIter=0;
       nonLinearIter < params::maxNonLinearIter; nonLinearIter++
      )
  {
    /* True residual, with explicit terms (not needed for Jacobian) */
    int numReadsResidual, numWritesResidual;

    array IdxToComputeInit = where(notConverged > 0);
    IdxToComputeInit.eval();
    computeResidual(primGuess, *residual, true, 
		    numReadsResidual, numWritesResidual,
		    IdxToComputeInit
		    );

    for (int var=0; var < vars::dof; var++)
    {
      /* Need residualSoA to compute norms */
      residualSoA(span, span, span, var) = residual->vars[var]*GZmask;

      /* Initialize primGuessPlusEps. Needed to numerically assemble the
       * Jacobian */
      primGuessPlusEps->vars[var]        = primGuess.vars[var];
    }

    /* Sum along last dim:vars to get L2 norm */
    l2Norm  = 
      af::sum(af::pow(residualSoA, 2.), 3);
    l2Norm.eval();
    notConverged      = l2Norm > params::nonlinearsolve_atol;
    array conditionIndices  = where(notConverged > 0);
    localNonConverged = conditionIndices.elements();

    /* Communicate residual */
    localresnorm = 
      af::norm(af::flat(residualSoA),
               AF_NORM_VECTOR_1
              );
    globalresnorm = localresnorm;
    globalNonConverged = localNonConverged;
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

    /*if(nonLinearIter==2 && localNonConverged>100)
      {
	array l2Norm_t  =
	  af::sum(af::pow(residualSoA, 2.), 3);
	array l2Norm2 = l2Norm_t*0.;
	l2Norm2(domainX1, domainX2, domainX3) = l2Norm_t(domainX1, domainX2, domainX3);
	l2Norm2.eval();
	array notConverged2      = l2Norm2 > 100.;
	array conditionIndices2  = where(notConverged2 > 0);

	PetscPrintf(PETSC_COMM_WORLD, "Points which have not converged:\n");
	array xCoords[3];
	geomCenter->getxCoords(xCoords);
	af_print(xCoords[0](conditionIndices2),12);
	af_print(xCoords[1](conditionIndices2),12);
	af_print(xCoords[2](conditionIndices2),12);
	for(int var=0;var<vars::dof;var++)
	  af_print(prim->vars[var](conditionIndices2),12);
	for(int var=0;var<vars::dof;var++)
	  af_print(residual->vars[var](conditionIndices2),12);
	af_print(elem->bSqr(conditionIndices2),12);
	af_print(elem->gammaLorentzFactor(conditionIndices2),12);
	af_print(elem->tau(conditionIndices2),12);
	af_print(elem->chi_emhd(conditionIndices2),12);
	af_print(elem->nu_emhd(conditionIndices2),12);
	}*/

    // Figure out which points should be solved for.
    // On the first step, solve even if the residual is small
    array IdxToComputeMain = where(notConverged > 0);
    if(nonLinearIter==0)
      IdxToComputeMain = IdxToComputeInit;
    IdxToComputeMain.eval();
    /* Residual without explicit terms, for faster Jacobian assembly */
    computeResidual(primGuess, *residual, false,
                    numReadsResidual, numWritesResidual,
		    IdxToComputeMain
		    );

    /* Assemble the Jacobian in Struct of Arrays format where the physics
     * operations are all vectorized */
    for (int row=0; row < vars::dof; row++)
    {
      /* Recommended value of Jacobian differencing parameter to achieve fp64
       * machine precision */
      double epsilon = params::JacobianAssembleEpsilon;

      array smallPrim = af::abs(primGuess.vars[row](IdxToComputeMain))<.5*epsilon;

      primGuessPlusEps->vars[row](IdxToComputeMain)  = 
	(1. + epsilon)*primGuess.vars[row](IdxToComputeMain)*(1.-smallPrim)
	+ smallPrim*epsilon; 

      computeResidual(*primGuessPlusEps, *residualPlusEps, false,
                      numReadsResidual, numWritesResidual,
		      IdxToComputeMain
                     );

      for (int column=0; column < vars::dof; column++)
      {
	array temp = jacobianSoA(span, span, span, column + vars::dof*row);
	temp(IdxToComputeMain) = (  residualPlusEps->vars[column](IdxToComputeMain)
				    -residual->vars[column](IdxToComputeMain)
				    )
	  /(primGuessPlusEps->vars[row](IdxToComputeMain)-primGuess.vars[row](IdxToComputeMain));
	  
        jacobianSoA(span, span, span, column + vars::dof*row)
          = temp;
      }
      /* reset */
      primGuessPlusEps->vars[row](IdxToComputeMain)  = primGuess.vars[row](IdxToComputeMain); 
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
    batchLinearSolve(jacobianAoS, bAoS, deltaPrimAoS,notConverged);

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
    array f0      = 0.5 * l2Norm(IdxToComputeMain);
    array fPrime0 = -2.*f0;
    
    /* Start with a full step */
    stepLength = 1.;
    int lineSearchIter=0;
    for (;lineSearchIter < params::maxLineSearchIters; lineSearchIter++
        )
    {
      /* 1) First take current step stepLength */
      for (int var=0; var<vars::dof; var++)
      {
	array dP = deltaPrimSoA(span, span, span, var);
        primGuessLineSearchTrial->vars[var](IdxToComputeMain) =  
          primGuess.vars[var](IdxToComputeMain) + stepLength(IdxToComputeMain)*dP(IdxToComputeMain);
      } 

      /* ...and then compute the norm */
      computeResidual(*primGuessLineSearchTrial, *residual, true,
                      numReadsResidual, numWritesResidual,
		      IdxToComputeMain
                     );
      for (int var=0; var<vars::dof; var++)
      {
        residualSoA(span, span, span, var) = residual->vars[var]*GZmask;
      }
      l2Norm = af::sum(af::pow(residualSoA, 2.), 3);
      array f1 = 0.5 * l2Norm(IdxToComputeMain);

      /* We have 3 pieces of information:
       * a) f(0)
       * b) f'(0) 
       * c) f(stepLength) 
       */
    
      const double alpha    = 1e-4;
      const double EPS      = params::linesearchfloor;
      array condition = f1 > (f0*(1. - alpha*stepLength(IdxToComputeMain)) +EPS);
      array denom     =   (f1-f0-fPrime0*stepLength(IdxToComputeMain)) * condition 
	+ (1.-condition);
      array nextStepLength =
        -fPrime0*stepLength(IdxToComputeMain)*stepLength(IdxToComputeMain)/denom/2.;
      stepLength(IdxToComputeMain) = stepLength(IdxToComputeMain)*(1. - condition) + condition*nextStepLength;
      
      array conditionIndices = where(condition > 0);
      if (conditionIndices.elements() == 0)
      {
        break;
      }
    }

    /* stepLength has now been set */
    for (int var=0; var<vars::dof; var++)
    {
      array dP = deltaPrimSoA(span, span, span, var);
      primGuess.vars[var](IdxToComputeMain) = 
        primGuess.vars[var](IdxToComputeMain) + stepLength(IdxToComputeMain)*dP(IdxToComputeMain);
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

void timeStepper::batchLinearSolve(const array &A, const array &b, array &x, const array UsePt)
{
  array Idx = where(UsePt > 0);
  int numVars = residual->numVars;
  int numPoints = Idx.elements();

  array UsePtVec = b*0;
  for(int v=0;v<numVars;v++)
    {
      array temp = UsePtVec(v,span,span,span);
      temp(Idx) = 1.;
      UsePtVec(v,span,span,span)=temp;
    }
  array IdxPtVec = where(UsePtVec>0);
  IdxPtVec.eval();
  array UsePtMat = A*0.;
  for(int v=0;v<numVars*numVars;v++)
    {
      array temp = UsePtMat(v,span,span,span);
      temp(Idx) = 1.;
      UsePtMat(v,span,span,span) = temp;
    }
  array IdxPtMat = where(UsePtMat>0);
  IdxPtMat.eval();

  array AIdx = A(IdxPtMat);
  array bIdx = b(IdxPtVec);
  array xIdx = x(IdxPtVec);
  AIdx.eval();
  bIdx.eval();
  xIdx.eval();
  AIdx.host(AHostPtr); 
  bIdx.host(bHostPtr);
  xIdx.host(xHostPtr);
  
  #pragma omp parallel for
  for (int k=0; k<numPoints; k++)
    {
      double ALocal[numVars*numVars];
      double bLocal[numVars];
      int pivot[numVars];

      const int spatialIndex = k;
	
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
  
  /* Copy solution to x on device */
  xIdx = array(numVars*numPoints, xHostPtr);
  x(IdxPtVec) = xIdx;
  x.eval();
}
