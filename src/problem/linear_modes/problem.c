#include "../problem.h"

#if (CONDUCTION)
  REAL kappaProblem;
  REAL tauProblem; 

  void setConductionParameters(const struct geometry geom[ARRAY_ARGS 1],
                               struct fluidElement elem[ARRAY_ARGS 1])
  {
    elem->kappa = kappaProblem;
    elem->tau   = tauProblem;
  }
#endif

#if (VISCOSITY)
  REAL etaProblem;
  REAL tauVisProblem; 

  void setViscosityParameters(const struct geometry geom[ARRAY_ARGS 1],
                               struct fluidElement elem[ARRAY_ARGS 1])
  {
    elem->eta     = etaProblem;
    elem->tauVis  = tauVisProblem;
  }
#endif

void initialConditions(struct timeStepper ts[ARRAY_ARGS 1])
{
  ARRAY(primOldGlobal);
  DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, 
                     ts->primPetscVecOld, &primOldGlobal);

  REAL randNum;
  PetscRandom randNumGen;
  PetscRandomCreate(PETSC_COMM_SELF, &randNumGen);
  PetscRandomSetType(randNumGen, PETSCRAND48);


  LOOP_OVER_TILES(ts->X1Size, ts->X2Size)
  {
    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);

      REAL XCoords[NDIM];
      getXCoords(&zone, CENTER, XCoords);


      #if (MODE==HYDRO_SOUND_MODE_1D) /* Eigenvalue = 3.09362659024*I */
        REAL primVars0[DOF];
        REAL complex deltaPrimVars[DOF];

        primVars0[RHO] = 1.;
        primVars0[UU]  = 2.;
        primVars0[U1]  = 0.;
        primVars0[U2]  = 0.;
        primVars0[U3]  = 0.;
        primVars0[B1]  = 0.;
        primVars0[B2]  = 0.;
        primVars0[B3]  = 0.;

        deltaPrimVars[RHO] = 0.345991032308;
        deltaPrimVars[UU]  = 0.922642752822;
        deltaPrimVars[U1]  = -0.170354208129;
        deltaPrimVars[U2]  = 0.;
        deltaPrimVars[U3]  = 0.;
        deltaPrimVars[B1]  = 0.;
        deltaPrimVars[B2]  = 0.;
        deltaPrimVars[B3]  = 0.;

        REAL k1 = 2*M_PI;
        REAL k2 = 0.;

        REAL complex mode = cexp(I*(k1*XCoords[1] + k2*XCoords[2]) );

        for (int var=0; var<DOF; var++)
        {
          INDEX_PETSC(primOldGlobal, &zone, var) =  
            primVars0[var] + AMPLITUDE*creal(deltaPrimVars[var]*mode);
        }

      #elif (MODE==ENTROPY_WAVE_1D) /* Eigenvalue = 0 */
        REAL primVars0[DOF];
        REAL complex deltaPrimVars[DOF];

        primVars0[RHO] = 1.;
        primVars0[UU]  = 2.;
        primVars0[U1]  = 0.;
        primVars0[U2]  = 0.;
        primVars0[U3]  = 0.;
        primVars0[B1]  = 0.;
        primVars0[B2]  = 0.;
        primVars0[B3]  = 0.;

        deltaPrimVars[RHO] = 1.;
        deltaPrimVars[UU]  = 0.;
        deltaPrimVars[U1]  = 0.;
        deltaPrimVars[U2]  = 0.;
        deltaPrimVars[U3]  = 0.;
        deltaPrimVars[B1]  = 0.;
        deltaPrimVars[B2]  = 0.;
        deltaPrimVars[B3]  = 0.;

        REAL k1 = 2*M_PI;
        REAL k2 = 0.;

        REAL complex mode = cexp(I*(k1*XCoords[1] + k2*XCoords[2]) );

        for (int var=0; var<DOF; var++)
        {
          INDEX_PETSC(primOldGlobal, &zone, var) =  
            primVars0[var] + AMPLITUDE*creal(deltaPrimVars[var]*mode);
        }

      #elif (MODE==CONDUCTION_STABLE_1D) 
        /* Eigenvalue = -0.498597173331 - 0.857832357798*I
         * kappa = 0.1
         * tau = 1.01818181818182 */
        kappaProblem = .1;
        tauProblem   = 1.01818181818182;

        REAL primVars0[DOF];
        REAL complex deltaPrimVars[DOF];

        primVars0[RHO] = 1.;
        primVars0[UU]  = 2.;
        primVars0[U1]  = 0.;
        primVars0[U2]  = 0.;
        primVars0[U3]  = 0.;
        primVars0[B1]  = 0.0000000001;
        primVars0[B2]  = 0.;
        primVars0[B3]  = 0.;
        primVars0[PHI] = 0.;
        
        deltaPrimVars[RHO] = 0.911376933183;
        deltaPrimVars[UU]  = 0.030751595371 - 0.0635975709194*I;
        deltaPrimVars[U1]  = 0.124428706971 - 0.0723215917578*I;
        deltaPrimVars[U2]  = 0.;
        deltaPrimVars[U3]  = 0.;
        deltaPrimVars[B1]  = 0.;
        deltaPrimVars[B2]  = 0.;
        deltaPrimVars[B3]  = 0.;
        deltaPrimVars[PHI] = -0.332658158111 + 0.181734443922*I;

        REAL k1 = 2*M_PI;
        REAL k2 = 0.;

        REAL complex mode = cexp(I*(k1*XCoords[1] + k2*XCoords[2]) );

        for (int var=0; var<DOF; var++)
        {
          INDEX_PETSC(primOldGlobal, &zone, var) =  
            primVars0[var] + AMPLITUDE*creal(deltaPrimVars[var]*mode);
        }
      #elif (MODE==CONDUCTION_STABLE_2D) 
        /* Eigenvalue = -0.498689052044 - 1.23434614343*I
         * kappa = 0.1
         * tau = 1.01818181818182 */
        kappaProblem = .1;
        tauProblem   = 1.01818181818182;

        REAL primVars0[DOF];
        REAL complex deltaPrimVars[DOF];

        primVars0[RHO] = 1.;
        primVars0[UU]  = 2.;
        primVars0[U1]  = 0.;
        primVars0[U2]  = 0.;
        primVars0[U3]  = 0.;
        primVars0[B1]  = 0.01;
        primVars0[B2]  = 0.;
        primVars0[B3]  = 0.;
        primVars0[PHI] = 0.;

        deltaPrimVars[RHO] = 0.912960868047;
        deltaPrimVars[UU]  = 0.0441633411305 - 0.0470501442451*I;
        deltaPrimVars[U1]  = 0.068161459988 - 0.0280266780212*I;
        deltaPrimVars[U2]  = 0.111191793414 - 0.0444339558103*I;
        deltaPrimVars[U3]  = 0.;
        deltaPrimVars[B1]  = -0.00130516937615 + 6.41592876069e-05*I;
        deltaPrimVars[B2]  = 0.00130516937615 - 6.41592876069e-05*I;
        deltaPrimVars[B3]  = 0.;
        deltaPrimVars[PHI] = -0.352802085664 + 0.134521891027*I;

        REAL k1 = 2*M_PI;
        REAL k2 = 2*M_PI;

        REAL complex mode = cexp(I*(k1*XCoords[1] + k2*XCoords[2]) );

        for (int var=0; var<DOF; var++)
        {
          INDEX_PETSC(primOldGlobal, &zone, var) =  
            primVars0[var] + AMPLITUDE*creal(deltaPrimVars[var]*mode);
        }
      #elif (MODE==VISCOSITY_2D) 
        /* Eigenvalue = -0.488186438695 + 0.516141145988*I
         * eta = 0.1
         * tau = 1.01818181818182 */
        etaProblem    = .1;
        tauVisProblem = 1.01818181818182;

        REAL primVars0[DOF];
        REAL complex deltaPrimVars[DOF];

        primVars0[RHO] = 1.;
        primVars0[UU]  = 2.;
        primVars0[U1]  = 0.;
        primVars0[U2]  = 0.;
        primVars0[U3]  = 0.;
        primVars0[B1]  = 0.01;
        primVars0[B2]  = 0.;
        primVars0[B3]  = 0.;
        primVars0[PSI] = 0.;

        deltaPrimVars[RHO] = 0.34457248379 - 8.32667268469e-17*I;
        deltaPrimVars[UU]  = 0.918859956773;
        deltaPrimVars[U1]  = -0.180848290246 - 0.00323785674963*I;
        deltaPrimVars[U2]  = 0.;
        deltaPrimVars[U3]  = 0.;
        deltaPrimVars[B1]  = 0.;
        deltaPrimVars[B2]  = 0.;
        deltaPrimVars[B3]  = 0.;
        deltaPrimVars[PSI] = 0.0624512526648 + 0.0186932205343*I;

        REAL k1 = 2*M_PI;
        REAL k2 = 0.;

        REAL complex mode = cexp(I*(k1*XCoords[1] + k2*XCoords[2]) );

        for (int var=0; var<DOF; var++)
        {
          INDEX_PETSC(primOldGlobal, &zone, var) =  
            primVars0[var] + AMPLITUDE*creal(deltaPrimVars[var]*mode);
        }
      #elif (MODE==VISCOSITY_1D)
        etaProblem    = .1;
        tauVisProblem = 1.01818181818182;

        REAL primVars0[DOF];
        REAL complex deltaPrimVars[DOF];

        primVars0[RHO] = 1.;
        primVars0[UU]  = 2.;
        primVars0[U1]  = 0.;
        primVars0[U2]  = 0.;
        primVars0[U3]  = 0.;
        primVars0[B1]  = 0.01;
        primVars0[B2]  = 0.02;
        primVars0[B3]  = 0.;
        primVars0[PSI] = 0.;

        deltaPrimVars[RHO] = 0.34457248379 - 8.32667268469e-17*I;
        deltaPrimVars[UU]  = 0.918859956773;
        deltaPrimVars[U1]  = -0.180848290246 - 0.00323785674963*I;
        deltaPrimVars[U2]  = 0.;
        deltaPrimVars[U3]  = 0.;
        deltaPrimVars[B1]  = 0.;
        deltaPrimVars[B2]  = 0.;
        deltaPrimVars[B3]  = 0.;
        deltaPrimVars[PSI] = 0.0624512526648 + 0.0186932205343*I;

        REAL k1 = 2*M_PI;

        REAL complex mode = cexp(I*k1*XCoords[1]);

        for (int var=0; var<DOF; var++)
        {
          INDEX_PETSC(primOldGlobal, &zone, var) =  
            primVars0[var] + AMPLITUDE*creal(deltaPrimVars[var]*mode);
        }
      #elif (MODE==ALFVEN_2D) 
        /* Eigenvalue = -0.488186438695 + 0.516141145988*I
         * eta = 0.1
         * tau = 1.01818181818182 */
        
        REAL primVars0[DOF];
        REAL complex deltaPrimVars[DOF];

        primVars0[RHO] = 1.;
        primVars0[UU]  = 2.;
        primVars0[U1]  = 0.;
        primVars0[U2]  = 0.;
        primVars0[U3]  = 0.;
        primVars0[B1]  = 0.01;
        primVars0[B2]  = 0.;
        primVars0[B3]  = 0.;
	
        deltaPrimVars[RHO] = 0.;
        deltaPrimVars[UU]  = 0.;
        deltaPrimVars[U1]  = 0.;
        deltaPrimVars[U2]  = 1.;
        deltaPrimVars[U3]  = 0.;
        deltaPrimVars[B1]  = 0.;
        deltaPrimVars[B2]  = 1.00005;
        deltaPrimVars[B3]  = 0.;

	etaProblem    = .1;
        tauVisProblem = 1.01818181818182;
	primVars0[PSI]  = 0.;
	deltaPrimVars[PSI]  = 0.;

        REAL k1 = 2*M_PI;
        REAL k2 = 0.;

        REAL complex mode = cexp(I*(k1*XCoords[1] + k2*XCoords[2]) );

        for (int var=0; var<DOF; var++)
        {
          INDEX_PETSC(primOldGlobal, &zone, var) =  
            primVars0[var] + AMPLITUDE*creal(deltaPrimVars[var]*mode);
        }
      #elif (MODE==FIREHOSE) 
        REAL primVars0[DOF];
        REAL complex deltaPrimVars[DOF];

        primVars0[RHO] = 1.;
        primVars0[UU]  = 0.01;
        primVars0[U1]  = 0.;
        primVars0[U2]  = 0.;
        primVars0[U3]  = 0.;
        primVars0[B1]  = 0.01;
        primVars0[B2]  = 0.;
        primVars0[B3]  = 0.;
	primVars0[PSI]  = 0.;
	
        deltaPrimVars[RHO] = 0.991060736449;
        deltaPrimVars[UU]  = 0.0132141431527 - 4.60785923306e-19*I;
        deltaPrimVars[U1]  = -0.131266091803 - 0.000591490468813*I;
        deltaPrimVars[U2]  = 1.e-4;
        deltaPrimVars[U3]  = 0.;
        deltaPrimVars[B1]  = 0.;
        deltaPrimVars[B2]  = -1.e-4;
        deltaPrimVars[B3]  = 0.;
	deltaPrimVars[PSI]  = 0.0198194260853 + 0.000238162630751*I;

	etaProblem    = 1.;
        tauVisProblem = 100.;

        REAL k1 = 2*M_PI;
        REAL k2 = 0.;

        REAL complex mode = cexp(I*(k1*XCoords[1] + k2*XCoords[2]) );

        for (int var=0; var<DOF; var++)
        {
	  if(var!=U2 && var!=B2)
	    INDEX_PETSC(primOldGlobal, &zone, var) =  
	      primVars0[var] + AMPLITUDE*creal(deltaPrimVars[var]*mode);
        }
	PetscRandomGetValue(randNumGen, &randNum);
	INDEX_PETSC(primOldGlobal, &zone, U2) = 1.e-8*randNum;
	PetscRandomGetValue(randNumGen, &randNum);
        INDEX_PETSC(primOldGlobal, &zone, B2) =1.e-8*randNum;
      #endif

    }
  }
  
  DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, 
                         ts->primPetscVecOld, &primOldGlobal);
}


void applyFloor(const int iTile, const int jTile,
                const int X1Start, const int X2Start,
                const int X1Size, const int X2Size,
                REAL primTile[ARRAY_ARGS TILE_SIZE])
{
  LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
  {
    struct gridZone zone;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start,
                X1Size, X2Size,
                &zone);

    if (primTile[INDEX_TILE(&zone, RHO)] < RHO_FLOOR)
    {
        primTile[INDEX_TILE(&zone, RHO)] = RHO_FLOOR;
    }

    if (primTile[INDEX_TILE(&zone, UU)] < UU_FLOOR)
    {
        primTile[INDEX_TILE(&zone, UU)] = UU_FLOOR;
    }
  }
}

void applyFloorInsideNonLinearSolver(const int iTile, const int jTile,
                                     const int X1Start, const int X2Start,
                                     const int X1Size, const int X2Size,
                                     REAL primTile[ARRAY_ARGS TILE_SIZE])
{

}

void applyAdditionalProblemSpecificBCs
(
  const int iTile, const int jTile,
  const int X1Start, const int X2Start,
  const int X1Size, const int X2Size,
  const struct problemData problemSpecificData[ARRAY_ARGS 1],
  REAL primTile[ARRAY_ARGS TILE_SIZE]
)
{

}

void applyProblemSpecificFluxFilter
(
  const int iTile, const int jTile,
  const int X1Start, const int X2Start,
  const int X1Size, const int X2Size,
  const struct problemData problemSpecificData[ARRAY_ARGS 1],
  REAL fluxX1Tile[ARRAY_ARGS TILE_SIZE],
  REAL fluxX2Tile[ARRAY_ARGS TILE_SIZE]
)
{

}

void halfStepDiagnostics(struct timeStepper ts[ARRAY_ARGS 1])
{

}

void fullStepDiagnostics(struct timeStepper ts[ARRAY_ARGS 1])
{

}

