#include "../problem.hpp"

namespace problemParams
{
  std::string rhoInputFile      = "atmosphere_soln_rho.txt";
  std::string uInputFile        = "atmosphere_soln_u.txt"  ;
  std::string phiInputFile      = "atmosphere_soln_phi.txt";
  std::string rCoordsInputFile  = "atmosphere_soln_rCoords.txt";
  array rhoLeft, rhoRight;
  array uLeft, uRight;
  array u1Left, u1Right;
};

void fluidElement::setFluidElementParameters(const geometry &geom)
{
  array xCoords[3],XCoords[3];
  geom.getXCoords(XCoords);
  geom.XCoordsToxCoords(XCoords,xCoords);

  array r = xCoords[directions::X1];

  chi_emhd = rho * af::pow(r, 0.5);
  tau      = 1.2 * chi_emhd/ (soundSpeed * soundSpeed);
}

void timeStepper::initialConditions(int &numReads, int &numWrites)
{
  PetscPrintf(PETSC_COMM_WORLD, "Setting up atmosphere...");

  double rho1D[N1+2*numGhost];
  double u1D[N1+2*numGhost];
  double phi1D[N1+2*numGhost];
  double rCoords1D[N1+2*numGhost];

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank); /* get current proc id */

  if (rank==0)
  {
    FILE *rhoFile;
    FILE *uFile;
    FILE *phiFile;
    FILE *rCoordsFile;
  
    char *rhoLine     = NULL;
    char *uLine       = NULL;
    char *phiLine     = NULL; 
    char *rCoordsLine = NULL;

    size_t rhoLen     = 0; ssize_t rhoRead;
    size_t uLen       = 0; ssize_t uRead;
    size_t phiLen     = 0; ssize_t phiRead;
    size_t rCoordsLen = 0; ssize_t rCoordsRead;

    rhoFile     = fopen(problemParams::rhoInputFile.c_str(), "r");
    uFile       = fopen(problemParams::uInputFile.c_str()  , "r");
    phiFile     = fopen(problemParams::phiInputFile.c_str(), "r");
    rCoordsFile = fopen(problemParams::rCoordsInputFile.c_str(), "r");

    if (   rhoFile      == NULL 
        || uFile        == NULL
        || phiFile      == NULL
        || rCoordsFile  == NULL
       )
    {
      PetscPrintf(PETSC_COMM_WORLD, "Input data files not found!\n");
      exit(1);
    }

    for (int i=-numGhost; i<N1+numGhost; i++)
    {
      rhoRead     = getline(&rhoLine    , &rhoLen     , rhoFile);
      uRead       = getline(&uLine      , &uLen       , uFile);
      phiRead     = getline(&phiLine    , &phiLen     , phiFile);
      rCoordsRead = getline(&rCoordsLine, &rCoordsLen , rCoordsFile);

      if (   rhoRead      == -1
          || uRead        == -1
          || phiRead      == -1
          || rCoordsRead  == -1
         )
      {
        PetscPrintf(PETSC_COMM_WORLD,
        "Found the input data files but error in reading them! Check number of grid zones in python script\n");
        exit(1);
      }

      rho1D[i+numGhost]     = atof(rhoLine);
      u1D[i+numGhost]       = atof(uLine);
      phi1D[i+numGhost]     = atof(phiLine);
      rCoords1D[i+numGhost] = atof(rCoordsLine);
    }

    free(rhoLine);
    free(uLine);
    free(phiLine);
    free(rCoordsLine);

    fclose(rhoFile);
    fclose(uFile);
    fclose(phiFile);
    fclose(rCoordsFile);
  }

  PetscPrintf(PETSC_COMM_WORLD, "file read complete\n");

  /* Broadcast the data from rank 0 proc to all other procs */
  MPI_Bcast(&rho1D[0]     , N1+2*numGhost, MPIU_SCALAR, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&u1D[0]       , N1+2*numGhost, MPIU_SCALAR, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&phi1D[0]     , N1+2*numGhost, MPIU_SCALAR, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&rCoords1D[0] , N1+2*numGhost, MPIU_SCALAR, 0, PETSC_COMM_WORLD);
  /* Broadcast complete */
  PetscPrintf(PETSC_COMM_WORLD, "Broadcast complete\n");

  array &rho = primOld->vars[vars::RHO];
  array &u   = primOld->vars[vars::U];
  array &u1  = primOld->vars[vars::U1];
  array &u2  = primOld->vars[vars::U2];
  array &u3  = primOld->vars[vars::U3];
  array &B1  = primOld->vars[vars::B1];
  array &B2  = primOld->vars[vars::B2];
  array &B3  = primOld->vars[vars::B3];

  array xCoords[3], XCoords[3];
  geomCenter->getXCoords(XCoords);
  geomCenter->XCoordsToxCoords(XCoords,xCoords);
  array r = xCoords[directions::X1];

  for (int i=0; i<primOld->N1Total; i++)
  {
    if (fabs(r(i, 0, 0).scalar<double>() - rCoords1D[i]) > 1e-10)
    {
      PetscPrintf(PETSC_COMM_WORLD,
                 "Mismatch in rCoords! Check r coords in python script. r = %.18f, rCoords = %.18f, i = %d\n",
                 r(i, 0, 0).scalar<double>(), rCoords1D[i], i
                 );
      exit(1);
    }
    rho(i, span, span) = rho1D[i];
    u(i, span, span)   = u1D[i];

    if (params::conduction)
    {
      array &q = primOld->vars[vars::Q];
      q(i, span, span) = phi1D[i];
    }
  }

  array uCon0 = 1./af::sqrt(-geomCenter->gCov[0][0]);
  array uCon1 = 0.*u1;
  array uCon2 = 0.*u2;
  array uCon3 = 0.*u3;

  /* Formula to output vUpr in Modified KerrSchild from uUpr in Boyer-Lindquist 
   * from Ben Ryan */
  array a = geomCenter->gCov[1][1];
  array b = geomCenter->gCon[0][1];
  array c = geomCenter->gCon[0][0];
  array vUprMKS = 
    (  c*uCon1/r 
     - sqrt(-a*b*b*b*b
            -a*b*b*c*uCon1*uCon1/(r*r) - b*b*c
           )
    )/(a*b*b + c);

  /* Set values in both bulk and boundaries. Later set values in bulk where we
   * need to add random fluctuations */
  u1 = vUprMKS;

  PetscRandom randNumGen;
  PetscRandomCreate(PETSC_COMM_SELF, &randNumGen);
  PetscRandomSetType(randNumGen, PETSCRAND48);
  int numGhostX1 = primOld->numGhostX1;
  int numGhostX2 = primOld->numGhostX2;
  int numGhostX3 = primOld->numGhostX3;
  int N1Local    = primOld->N1Local;
  int N2Local    = primOld->N2Local;
  int N3Local    = primOld->N3Local;

  for (int k=numGhostX3; k < N3Local+numGhostX3; k++)
  {
    for (int j=numGhostX2; j < N2Local+numGhostX2; j++)
    {
      for (int i=numGhostX1; i < N1Local+numGhostX1; i++)
      {
        double randNum = 0.5;
        PetscRandomGetValue(randNumGen, &randNum);

        u1(i,j,k) = 
          vUprMKS(i,j,k)
        * (1. + params::InitialPerturbationAmplitude*(randNum-0.5));
      }
    }
  }
  PetscRandomDestroy(&randNumGen);

  //B2BottomOld->vars[0] = 1e-5/geomBottom->g;
  B1LeftOld->vars[0] = 1e-5/geomLeft->g;

  af::seq leftBoundary(0, numGhost);
  af::seq rightBoundary(primOld->N1Local + numGhost, 
                        primOld->N1Local + 2*numGhost-1
                       );
  if (primOld->iLocalStart == 0)
  {
    problemParams::rhoLeft = rho(leftBoundary, span, span);
    problemParams::uLeft   = u(leftBoundary,   span, span);
    problemParams::u1Left  = u1(leftBoundary,  span, span);
  }

  if (primOld->iLocalEnd == N1)
  {
    problemParams::rhoRight = rho(rightBoundary, span, span); 
    problemParams::uRight   = u(rightBoundary,   span, span); 
    problemParams::u1Right  = u1(rightBoundary,  span, span); 
  }
  
  PetscPrintf(PETSC_COMM_WORLD, "done\n");
  /* Done with setting the initial conditions */

  af::seq domainX1 = *primOld->domainX1;
  af::seq domainX2 = *primOld->domainX2;
  af::seq domainX3 = *primOld->domainX3;
  
  computeDivB(*B1LeftOld, *B2BottomOld, *B3BackOld, numReads, numWrites);
  double divBNorm = 
    af::norm(af::flat(divB->vars[0](domainX1, domainX2, domainX3)), AF_NORM_VECTOR_1);
  printf("divB = %g\n", divBNorm);
}



void timeStepper::setProblemSpecificBCs(grid &prim,
                                        grid &B1Left,
                                        grid &B2Bottom,
                                        grid &B3Back,
                                        int &numReads, 
                                        int &numWrites
                                       )
{
  af::seq leftBoundary(0, numGhost);
  af::seq rightBoundary(primOld->N1Local + numGhost, 
                        primOld->N1Local + 2*numGhost-1
                       );

  af::seq bottomBoundary(0, numGhost);
  af::seq topBoundary(primOld->N2Local + numGhost,
                      primOld->N2Local + 2*numGhost-1
                     );

  if (primOld->iLocalStart == 0)
  {
    prim.vars[vars::RHO](leftBoundary, span, span) = problemParams::rhoLeft;
    prim.vars[vars::U](leftBoundary,   span, span) = problemParams::uLeft; 
    prim.vars[vars::U1](leftBoundary,  span, span) = problemParams::u1Left;
    B1Left.vars[0](leftBoundary, span, span)   =1e-5/geomLeft->g(leftBoundary, span, span) ;
    B2Bottom.vars[0](leftBoundary, span, span) = 0e-5/geomBottom->g(leftBoundary, span, span);
    B3Back.vars[0](leftBoundary, span, span)   = 
      0e-5/geomCenter->g(leftBoundary, span, span);
    
    if (params::conduction)
    {
      prim.vars[vars::Q](leftBoundary, span, span) = 0.;
    }
  }

  if (primOld->iLocalEnd == N1)
  {
    prim.vars[vars::RHO](rightBoundary, span, span) =  problemParams::rhoRight;
    prim.vars[vars::U](rightBoundary,   span, span) =  problemParams::uRight;
    prim.vars[vars::U1](rightBoundary,  span, span) =  problemParams::u1Right;
    B1Left.vars[0](rightBoundary, span, span)   = 1e-5/geomLeft->g(rightBoundary, span, span);
    B2Bottom.vars[0](rightBoundary, span, span) = 0e-5/geomBottom->g(rightBoundary, span, span);
    B3Back.vars[0](rightBoundary, span, span)   =
      0e-5/geomCenter->g(rightBoundary, span, span);

    if (params::conduction)
    {
      prim.vars[vars::Q](rightBoundary, span, span) = 0.;
    }
  }

  if (primOld->jLocalStart == 0)
  {
    B1Left.vars[0](span, bottomBoundary, span) = 0.;
    B2Bottom.vars[0](span, bottomBoundary, span) =
      0.*1e-5/geomBottom->g(span, bottomBoundary, span); ;
    B3Back.vars[0](span, bottomBoundary, span) =
      0e-5/geomCenter->g(span, bottomBoundary, span);

    if (params::conduction)
    {
      prim.vars[vars::Q](rightBoundary, span, span) = 0.;
    }
  }

  if (primOld->jLocalEnd == N2)
  {
    B1Left.vars[0](span, topBoundary, span) = 0.;
    B2Bottom.vars[0](span, topBoundary, span) =
      0.*1e-5/geomBottom->g(span, topBoundary, span);
    B3Back.vars[0](span, topBoundary, span) =
      0e-5/geomCenter->g(span, topBoundary, span);

    if (params::conduction)
    {
      prim.vars[vars::Q](rightBoundary, span, span) = 0.;
    }
  }

}

void timeStepper::halfStepDiagnostics(int &numReads,int &numWrites)
{
}

void timeStepper::fullStepDiagnostics(int &numReads,int &numWrites)
{

  int world_rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(PETSC_COMM_WORLD, &world_size);

  af::seq domainX1 = *primOld->domainX1;
  af::seq domainX2 = *primOld->domainX2;
  af::seq domainX3 = *primOld->domainX3;
  
  // Time step control
  array minSpeedTemp,maxSpeedTemp;
  array minSpeed,maxSpeed;
  elemOld->computeMinMaxCharSpeeds(*geomCenter,directions::X1,minSpeedTemp,maxSpeedTemp,numReads,numWrites);
  minSpeedTemp = minSpeedTemp/XCoords->dX1;
  maxSpeedTemp = maxSpeedTemp/XCoords->dX1;
  minSpeed=minSpeedTemp;
  maxSpeed=maxSpeedTemp;
  if(params::dim>1)
    {
      elemOld->computeMinMaxCharSpeeds(*geomCenter,directions::X2,minSpeedTemp,maxSpeedTemp,numReads,numWrites);
      minSpeedTemp = minSpeedTemp/XCoords->dX2;
      maxSpeedTemp = maxSpeedTemp/XCoords->dX2;
      minSpeed=af::min(minSpeed,minSpeedTemp);
      maxSpeed=af::max(maxSpeed,maxSpeedTemp);
    }
  if(params::dim>2)
    {
      elemOld->computeMinMaxCharSpeeds(*geomCenter,directions::X3,minSpeedTemp,maxSpeedTemp,numReads,numWrites);
      minSpeedTemp = minSpeedTemp/XCoords->dX3;
      maxSpeedTemp = maxSpeedTemp/XCoords->dX3;
      minSpeed=af::min(minSpeed,minSpeedTemp);
      maxSpeed=af::max(maxSpeed,maxSpeedTemp);
    }
  maxSpeed = af::max(maxSpeed,af::abs(minSpeed));
  maxSpeed.eval();
  array maxInvDt_af = af::max(af::max(af::max(maxSpeed,2),1),0);
  double maxInvDt = maxInvDt_af.host<double>()[0];
  /* Use MPI to find minimum over all processors */
  if (world_rank == 0) 
    {
      double temp; 
      for(int i=1;i<world_size;i++)
	{
	  MPI_Recv(&temp, 1, MPI_DOUBLE, i, i, PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
	  if( maxInvDt < temp)
	    maxInvDt = temp;
	}
    }
  else
    {
      MPI_Send(&maxInvDt, 1, MPI_DOUBLE, 0, world_rank, PETSC_COMM_WORLD);
    }
  MPI_Barrier(PETSC_COMM_WORLD);
  MPI_Bcast(&maxInvDt,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
  MPI_Barrier(PETSC_COMM_WORLD);
  dt = params::CourantFactor/maxInvDt;
  PetscPrintf(PETSC_COMM_WORLD,"New dt = %e\n",dt);

   bool WriteData = (floor(time/params::WriteDataEveryDt) != floor((time-dt)/params::WriteDataEveryDt));
  if(WriteData)
  {
    const int WriteIdx = floor(time/params::WriteDataEveryDt);
    if(WriteIdx==1)
	  {
	    PetscPrintf(PETSC_COMM_WORLD, "Printing gCov\n");
	    geomCenter->setgCovGrid();
	    geomCenter->gCovGrid->dump("gCov","gCovCenter.h5");
	    geomCenter->setgConGrid();
	    geomCenter->gConGrid->dump("gCon","gConCenter.h5");
	    geomCenter->setgGrid();
	    geomCenter->gGrid->dump("sqrtDetg","sqrtDetgCenter.h5");
	    geomCenter->setxCoordsGrid();
	    geomCenter->xCoordsGrid->dump("xCoords","xCoordsCenter.h5");

	    geomLeft->setgCovGrid();
	    geomLeft->gCovGrid->dump("gCov","gCovLeft.h5");
	    geomLeft->setgConGrid();
	    geomLeft->gConGrid->dump("gCon","gConLeft.h5");
	    geomLeft->setgGrid();
	    geomLeft->gGrid->dump("sqrtDetg","sqrtDetgLeft.h5");
	    geomLeft->setxCoordsGrid();
	    geomLeft->xCoordsGrid->dump("xCoords","xCoordsLeft.h5");


	    geomBottom->setgCovGrid();
	    geomBottom->gCovGrid->dump("gCov","gCovBottom.h5");
	    geomBottom->setgConGrid();
	    geomBottom->gConGrid->dump("gCon","gConBottom.h5");
	    geomBottom->setgGrid();
	    geomBottom->gGrid->dump("sqrtDetg","sqrtDetgBottom.h5");
	    geomBottom->setxCoordsGrid();
	    geomBottom->xCoordsGrid->dump("xCoords","xCoordsBottom.h5");
	  }
      
    std::string filename   = "primVarsT";
    std::string B1filename = "B1Left";
    std::string B2filename = "B2Bottom";
    std::string B3filename = "B3Back";
    std::string s_idx = std::to_string(WriteIdx);
      
    for(int i=0;i<6-s_idx.size();i++)
    {
	    filename=filename+"0";
	    B1filename=B1filename+"0";
	    B2filename=B2filename+"0";
	    B3filename=B3filename+"0";
    }
      
    filename=filename+s_idx;
    filename=filename+".h5";
    primOld->dump("primitives",filename);

    B1filename=B1filename+s_idx;
    B1filename=B1filename+".h5";
    B1LeftOld->dump("B1", B1filename);

    B2filename=B2filename+s_idx;
    B2filename=B2filename+".h5";
    B2BottomOld->dump("B2", B2filename);

    B3filename=B3filename+s_idx;
    B3filename=B3filename+".h5";
    B3BackOld->dump("B3", B3filename);


  }

  computeDivB(*B1LeftOld, *B2BottomOld, *B3BackOld, numReads, numWrites);
  double divBNorm = 
    af::norm(af::flat(divB->vars[0](domainX1, domainX2, domainX3)), AF_NORM_VECTOR_1);
  printf("divB = %g\n", divBNorm);

}
