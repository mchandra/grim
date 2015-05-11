#include "timestepper.h"

PetscErrorCode computeResidual(SNES snes, 
                               Vec primPetscVec,
                               Vec residualPetscVec,
                               void *ptr)
{
  struct timeStepper *ts = (struct timeStepper*)ptr;

  int X1Start, X2Start, X3Start;
  int X1Size, X2Size, X3Size;

  DMDAGetCorners(ts->dmdaWithGhostZones,
                 &X1Start, &X2Start, &X3Start,
                 &X1Size, &X2Size, &X3Size);

  ARRAY(primGlobal);
  ARRAY(primOldGlobal);
  ARRAY(primHalfStepGlobal);
  ARRAY(divFluxOldGlobal);
  ARRAY(sourceTermsOldGlobal);
  ARRAY(conservedVarsOldGlobal);
  ARRAY(connectionGlobal);
  ARRAY(dtGlobal);
  ARRAY(residualGlobal);

  DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, primPetscVec, 
                     &primGlobal);
  DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVecOld,
                     &primOldGlobal); 
  DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVecHalfStep,
                     &primHalfStepGlobal); 
  DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, residualPetscVec,
                     &residualGlobal);

  DMDAVecGetArrayDOF(ts->dmdaFluxes, ts->divFluxPetscVecOld,
                     &divFluxOldGlobal);
  DMDAVecGetArrayDOF(ts->dmdaFluxes, ts->sourceTermsPetscVecOld,
                     &sourceTermsOldGlobal);           
  DMDAVecGetArrayDOF(ts->dmdaFluxes, ts->conservedVarsPetscVecOld,
                     &conservedVarsOldGlobal); 

  DMDAVecGetArrayDOF(ts->connectionDMDA, ts->connectionPetscVec,
                     &connectionGlobal);

  DMDAVecGetArrayDOF(ts->dmdaDt, ts->dtPetscVec, &dtGlobal);

  #if (CONDUCTION)
    ARRAY(gradTGlobal);
    ARRAY(graduConGlobal);
    ARRAY(graduConHigherOrderTerm1Global);
    ARRAY(graduConHigherOrderTerm2Global);

    DMDAVecGetArrayDOF(ts->gradTDM, ts->gradTPetscVec, &gradTGlobal);
    DMDAVecGetArrayDOF(ts->graduConDM, ts->graduConPetscVec, 
                       &graduConGlobal);
    DMDAVecGetArrayDOF(ts->graduConHigherOrderTermsDM, 
                       ts->graduConHigherOrderTerm1PetscVec,
                       &graduConHigherOrderTerm1Global);
    DMDAVecGetArrayDOF(ts->graduConHigherOrderTermsDM, 
                       ts->graduConHigherOrderTerm2PetscVec,
                       &graduConHigherOrderTerm2Global);
  #endif

  #if (REAPER)
    ARRAY(momentsGlobal);
    DMDAVecGetArrayDOF(ts->momentsDM, ts->momentsVec, &momentsGlobal);
  #endif

  if (ts->computeOldSourceTermsAndOldDivOfFluxes)
  {
    Vec primPetscVecOldLocal, primPetscVecHalfStepLocal;
    DMGetLocalVector(ts->dmdaWithGhostZones, &primPetscVecOldLocal);
    DMGetLocalVector(ts->dmdaWithGhostZones, &primPetscVecHalfStepLocal);

    /* Exchange ghost zone data. */
    if (ts->computeDivOfFluxAtTimeN)
    {
      /* Compute Div(flux) at t=n */
      DMGlobalToLocalBegin(ts->dmdaWithGhostZones, 
                           ts->primPetscVecOld,
                           INSERT_VALUES,
                           primPetscVecOldLocal);
      DMGlobalToLocalEnd(ts->dmdaWithGhostZones,
                         ts->primPetscVecOld,
                         INSERT_VALUES,
                         primPetscVecOldLocal);
    }
    else if (ts->computeDivOfFluxAtTimeNPlusHalf)
    {
      /* Compute Div(flux) at t=n+1/2 */
      DMGlobalToLocalBegin(ts->dmdaWithGhostZones, 
                           ts->primPetscVecHalfStep,
                           INSERT_VALUES,
                           primPetscVecHalfStepLocal);
      DMGlobalToLocalEnd(ts->dmdaWithGhostZones,
                         ts->primPetscVecHalfStep,
                         INSERT_VALUES,
                         primPetscVecHalfStepLocal);
    }

    ARRAY(primOldLocal);
    ARRAY(primHalfStepLocal);

    DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, primPetscVecOldLocal,
                       &primOldLocal);
    DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, primPetscVecHalfStepLocal,
                       &primHalfStepLocal);


    /* Loop through tiles. We use tiles to maximize cache usage.*/
    #if (USE_OPENMP)
      #pragma omp parallel for
    #endif
    LOOP_OVER_TILES(X1Size, X2Size)
    {
      REAL primTile[TILE_SIZE];
      REAL fluxX1Tile[TILE_SIZE], fluxX2Tile[TILE_SIZE];

      /* Load data from the global memory on RAM onto a tile small enough
        * to reside on the cache */
      if (ts->computeDivOfFluxAtTimeN)
      {
        LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
        {
          struct gridZone zone;
          setGridZone(iTile, jTile,
                      iInTile, jInTile,
                      X1Start, X2Start, 
                      X1Size, X2Size, 
                      &zone);
          for (int var=0; var<DOF; var++)
          {
            primTile[INDEX_TILE(&zone, var)] =
            INDEX_PETSC(primOldLocal, &zone, var);
          }
        }
      } 
      else if (ts->computeDivOfFluxAtTimeNPlusHalf)
      {
        LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
        {
          struct gridZone zone;
          setGridZone(iTile, jTile,
                      iInTile, jInTile,
                      X1Start, X2Start, 
                      X1Size, X2Size, 
                      &zone);
          for (int var=0; var<DOF; var++)
          {
            primTile[INDEX_TILE(&zone, var)] =
            INDEX_PETSC(primHalfStepLocal, &zone, var);
          }
        }
      }

      /* Sync point */
      
      applyTileBoundaryConditions(iTile, jTile,
                                  X1Start, X2Start,
                                  X1Size, X2Size,
                                  primTile);

      applyAdditionalProblemSpecificBCs(iTile, jTile,
                                        X1Start, X2Start,
                                        X1Size, X2Size,
                                        ts->problemSpecificData,
                                        primTile);

      /* Sync point */
  
      /* Work on the tiles.*/
      computeFluxesOverTile(primTile,
                            iTile, jTile,
                            X1Start, X2Start,
                            X1Size, X2Size,
                            fluxX1Tile, fluxX2Tile,
                            dtGlobal);

      applyProblemSpecificFluxFilter(iTile, jTile,
                                     X1Start, X2Start,
                                     X1Size, X2Size,
                                     ts->problemSpecificData,
                                     fluxX1Tile, fluxX2Tile);

      LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
      {
        struct gridZone zone;
        setGridZone(iTile, jTile,
                    iInTile, jInTile,
                    X1Start, X2Start, 
                    X1Size, X2Size, 
                    &zone);

        REAL XCoords[NDIM], sourceTerms[NUM_FLUXES], conservedVars[NUM_FLUXES];

        getXCoords(&zone, CENTER, XCoords);
        struct geometry geom; setGeometry(XCoords, &geom);

        /* Now we need to compute conservedVarsOld using data from
         * primOldLocal. */
        struct fluidElement elem;
        setFluidElement(&INDEX_PETSC(primOldLocal, &zone, 0), &geom, &elem);
        computeFluxes(&elem, &geom, 0, conservedVars);

        #if (REAPER)
          for (int mu=0; mu<NDIM; mu++)
          {
            INDEX_PETSC(momentsGlobal, &zone, N_UP(mu)) = 
              elem.moments[N_UP(mu)];

            for (int nu=0; nu<NDIM; nu++)
            {
              INDEX_PETSC(momentsGlobal, &zone, T_UP_DOWN(mu, nu)) =
                elem.moments[T_UP_DOWN(mu, nu)];

              #if (REAPER_MOMENTS==15)
              	for (int lambda=0; lambda<NDIM; lambda++)
              	{
                  INDEX_PETSC(momentsGlobal, &zone, M_UP_UP_UP(mu, nu, lambda)) =
                    elem.moments[M_UP_UP_UP(mu, nu, lambda)];
                }
              #endif
            }
          }
        #endif

        for (int var=0; var<NUM_FLUXES; var++)
        {
          INDEX_PETSC(conservedVarsOldGlobal, &zone, var) = 
            conservedVars[var];

          INDEX_PETSC(divFluxOldGlobal, &zone, var) = 
            (  fluxX1Tile[INDEX_TILE_PLUS_ONE_X1(&zone, var)]
             - fluxX1Tile[INDEX_TILE(&zone, var)]
            )/zone.dX1
          #if (COMPUTE_DIM==2)
          + 
            (  fluxX2Tile[INDEX_TILE_PLUS_ONE_X2(&zone, var)]
             - fluxX2Tile[INDEX_TILE(&zone, var)]
            )/zone.dX2
          #endif
            ;
        }

        if (ts->computeSourceTermsAtTimeN)
        {
          computeSourceTerms(&elem, &geom,
                             &INDEX_PETSC(connectionGlobal, &zone, 0),
                             sourceTerms);

          for (int var=0; var<NUM_FLUXES; var++)
          {
            INDEX_PETSC(sourceTermsOldGlobal, &zone, var) = 
              sourceTerms[var];
          }
        }
        else if (ts->computeSourceTermsAtTimeNPlusHalf)
        {
          setFluidElement(&INDEX_PETSC(primHalfStepLocal, &zone, 0),
                          &geom, &elem);

          computeSourceTerms(&elem, &geom,
                             &INDEX_PETSC(connectionGlobal, &zone, 0),
                             sourceTerms);

          for (int var=0; var<NUM_FLUXES; var++)
          {
            INDEX_PETSC(sourceTermsOldGlobal, &zone, var) = 
              sourceTerms[var];
          }
        }

      }

      #if (CONDUCTION)
        addConductionSourceTermsToResidual
          (primTile,
           primGlobal, primHalfStepGlobal, primOldGlobal,
           connectionGlobal, 
           gradTGlobal, graduConGlobal, 
           graduConHigherOrderTerm1Global,
           graduConHigherOrderTerm2Global,
           ts->dt,
           ts->computeOldSourceTermsAndOldDivOfFluxes,
           ts->computeDivOfFluxAtTimeN,
           ts->computeDivOfFluxAtTimeNPlusHalf,
           iTile, jTile, X1Start, X2Start, X1Size, X2Size,
           residualGlobal
          );
      #endif

    } /* End of LOOP_OVER_TILES */

    DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, primPetscVecOldLocal,
                           &primOldLocal);
    DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, primPetscVecHalfStepLocal,
                           &primHalfStepLocal);

    DMRestoreLocalVector(ts->dmdaWithGhostZones, &primPetscVecOldLocal);
    DMRestoreLocalVector(ts->dmdaWithGhostZones, &primPetscVecHalfStepLocal);

    /* All old sources and divFluxes have now been computed */
    ts->computeOldSourceTermsAndOldDivOfFluxes = 0;
  }

  /* The following computation requires no communication*/

  #if (TIME_STEPPING==IMPLICIT)
    Vec primPetscVecLocal;
    DMGetLocalVector(ts->dmdaWithGhostZones, &primPetscVecLocal);

    /* Exchange ghost zone data. */
    DMGlobalToLocalBegin(ts->dmdaWithGhostZones, 
                         primPetscVec,
                         INSERT_VALUES,
                         primPetscVecLocal);
    DMGlobalToLocalEnd(ts->dmdaWithGhostZones,
                       primPetscVec,
                       INSERT_VALUES,
                       primPetscVecLocal);

    ARRAY(primLocal);

    DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, primPetscVecLocal, &primLocal);
  #endif

  #if (USE_OPENMP)
    #pragma omp parallel for
  #endif
  LOOP_OVER_TILES(X1Size, X2Size)
  {
    REAL primTile[TILE_SIZE];
    REAL fluxX1Tile[TILE_SIZE], fluxX2Tile[TILE_SIZE];

    #if (TIME_STEPPING==IMPLICIT)
      LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
      {
        struct gridZone zone;
        setGridZone(iTile, jTile,
                    iInTile, jInTile,
                    X1Start, X2Start, 
                    X1Size, X2Size, 
                    &zone);
        for (int var=0; var<DOF; var++)
        {
          primTile[INDEX_TILE(&zone, var)] =
          INDEX_PETSC(primLocal, &zone, var);
        }
      }
      /* Sync point */
    
      /* Apply boundary conditions on each tile */
      applyTileBoundaryConditions(iTile, jTile,
                                  X1Start, X2Start,
                                  X1Size, X2Size,
                                  primTile);

      applyAdditionalProblemSpecificBCs(iTile, jTile,
                                        X1Start, X2Start,
                                        X1Size, X2Size,
                                        ts->problemSpecificData,
                                        primTile);

      computeFluxesOverTile(primTile, 
                            iTile, jTile,
                            X1Start, X2Start,
                            X1Size, X2Size,
                            fluxX1Tile, fluxX2Tile,
                            dtGlobal);

      applyProblemSpecificFluxFilter(iTile, jTile,
                                     X1Start, X2Start,
                                     X1Size, X2Size,
                                     ts->problemSpecificData,
                                     fluxX1Tile, fluxX2Tile);
    #endif

    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  X1Start, X2Start, 
                  X1Size, X2Size, 
                  &zone);

      REAL XCoords[NDIM];

      getXCoords(&zone, CENTER, XCoords);
      struct geometry geom; setGeometry(XCoords, &geom);
      struct fluidElement elem;
      setFluidElement(&INDEX_PETSC(primGlobal, &zone, 0), &geom, &elem);

      REAL conservedVars[NUM_FLUXES];
      computeFluxes(&elem, &geom, 0, conservedVars);

      #if (TIME_STEPPING==IMEX || TIME_STEPPING==IMPLICIT)
        REAL sourceTerms[NUM_FLUXES];
        computeSourceTerms(&elem, &geom,
                           &INDEX_PETSC(connectionGlobal, &zone, 0),
                           sourceTerms);
      #endif

      REAL g = sqrt(-geom.gDet);
      REAL norm = g;

      REAL residual[NUM_FLUXES]; 
      for (int var=0; var<NUM_FLUXES; var++)
      {
        #if (TIME_STEPPING==EXPLICIT)

          residual[var] =
          ( (  conservedVars[var]
             - INDEX_PETSC(conservedVarsOldGlobal, &zone, var)
            )/ts->dt
            + INDEX_PETSC(divFluxOldGlobal, &zone, var)
            - INDEX_PETSC(sourceTermsOldGlobal, &zone, var)
          )/norm;

        #elif (TIME_STEPPING==IMEX)

          residual[var] =
          ( (  conservedVars[var]
             - INDEX_PETSC(conservedVarsOldGlobal, &zone, var)
            )/ts->dt
            + INDEX_PETSC(divFluxOldGlobal, &zone, var)
            - 0.5*(  INDEX_PETSC(sourceTermsOldGlobal, &zone, var)
                   + sourceTerms[var]
                  )
          )/norm;

        #elif (TIME_STEPPING==IMPLICIT)

          residual[var] =
          ( (  conservedVars[var]
             - INDEX_PETSC(conservedVarsOldGlobal, &zone, var)
            )/ts->dt
            + 0.5*(  INDEX_PETSC(divFluxOldGlobal, &zone, var)
                   + 
                    (  fluxX1Tile[INDEX_TILE_PLUS_ONE_X1(&zone, var)]
                     - fluxX1Tile[INDEX_TILE(&zone, var)]
                    )/zone.dX1
                  #if (COMPUTE_DIM==2)
                   + 
                    (  fluxX2Tile[INDEX_TILE_PLUS_ONE_X2(&zone, var)]
                     - fluxX2Tile[INDEX_TILE(&zone, var)]
                    )/zone.dX2
                  #endif
                  )
            - 0.5*(  INDEX_PETSC(sourceTermsOldGlobal, &zone, var)
                   + sourceTerms[var]
                  )
          )/norm;

        #endif
      }
      #if (REAPER_MOMENTS==15 && GYROAVERAGING)
        /* NUM_FLUXES > DOF. Need to contract so that final number of residuals
         * == DOF */
        REAL divMUpUp[NDIM][NDIM];
        divMUpUp[0][0] = 0.;
        divMUpUp[0][1] = residual[B01_FLUX];
        divMUpUp[0][2] = residual[B02_FLUX];
        divMUpUp[0][3] = residual[B03_FLUX];
        divMUpUp[1][0] = residual[B01_FLUX];
        divMUpUp[1][1] = residual[B11_FLUX];
        divMUpUp[1][2] = residual[B12_FLUX];
        divMUpUp[1][3] = residual[B13_FLUX];
        divMUpUp[2][0] = residual[B02_FLUX];
        divMUpUp[2][1] = residual[B12_FLUX];
        divMUpUp[2][2] = residual[B22_FLUX];
        divMUpUp[2][3] = residual[B23_FLUX];
        divMUpUp[3][0] = residual[B03_FLUX];
        divMUpUp[3][1] = residual[B13_FLUX];
        divMUpUp[3][2] = residual[B23_FLUX];
        divMUpUp[3][3] = residual[B33_FLUX];

        residual[B01_FLUX] = 0.;
        residual[B02_FLUX] = 0.;
        residual[B11_FLUX] = 0.;
        residual[B12_FLUX] = 0.;
        residual[B22_FLUX] = 0.;
        for (int mu=0; mu<NDIM; mu++)
        {
          for (int nu=0; nu<NDIM; nu++)
          {
            residual[B01_FLUX] +=
              elem.eDownNoHatUpHat[mu][0]
            * elem.eDownNoHatUpHat[nu][1]
            * divMUpUp[mu][nu];

            residual[B02_FLUX] +=
              elem.eDownNoHatUpHat[mu][0]
            * elem.eDownNoHatUpHat[nu][2]
            * divMUpUp[mu][nu];

            residual[B11_FLUX] +=
              elem.eDownNoHatUpHat[mu][1]
            * elem.eDownNoHatUpHat[nu][1]
            * divMUpUp[mu][nu];

            residual[B12_FLUX] +=
              elem.eDownNoHatUpHat[mu][1]
            * elem.eDownNoHatUpHat[nu][2]
            * divMUpUp[mu][nu];

            residual[B22_FLUX] +=
              elem.eDownNoHatUpHat[mu][2]
            * elem.eDownNoHatUpHat[nu][2]
            * divMUpUp[mu][nu];
          }
        }

        INDEX_PETSC(residualGlobal, &zone, ALPHA) =
          residual[ALPHA_FLUX];

        INDEX_PETSC(residualGlobal, &zone, A0) = 
          residual[A0_FLUX];

        INDEX_PETSC(residualGlobal, &zone, U1) = 
          residual[U1_FLUX];

        INDEX_PETSC(residualGlobal, &zone, U2) = 
          residual[U2_FLUX];

        INDEX_PETSC(residualGlobal, &zone, U3) = 
          residual[U3_FLUX];

        INDEX_PETSC(residualGlobal, &zone, B01) = 
          residual[B01_FLUX];

        INDEX_PETSC(residualGlobal, &zone, B02) = 
          residual[B02_FLUX];

        INDEX_PETSC(residualGlobal, &zone, B11) = 
          residual[B11_FLUX];

        INDEX_PETSC(residualGlobal, &zone, B12) = 
          residual[B12_FLUX];

        INDEX_PETSC(residualGlobal, &zone, B22) = 
          residual[B22_FLUX];

        INDEX_PETSC(residualGlobal, &zone, B1) = 
          residual[B1_FLUX];

        INDEX_PETSC(residualGlobal, &zone, B2) = 
          residual[B2_FLUX];

        INDEX_PETSC(residualGlobal, &zone, B3) = 
          residual[B3_FLUX];

      #else
        /* No GYROAVERAGING. NUM_FLUXES == DOF. Generic case. Easy peasy */
        for (int var=0; var<DOF; var++)
        {
          INDEX_PETSC(residualGlobal, &zone, var) = residual[var];
        }
        INDEX_PETSC(residualGlobal, &zone, F0_ALPHA) = 
          elem.collisionIntegrals[0];
        INDEX_PETSC(residualGlobal, &zone, F0_A0) = 
          elem.collisionIntegrals[1];
        INDEX_PETSC(residualGlobal, &zone, F0_A1) = 
          elem.collisionIntegrals[2];
        INDEX_PETSC(residualGlobal, &zone, F0_A2) = 
          elem.collisionIntegrals[3];
        INDEX_PETSC(residualGlobal, &zone, F0_A3) = 
          elem.collisionIntegrals[4];
      #endif /* No GYROAVERAGING */
    }

    #if (CONDUCTION)
      addConductionSourceTermsToResidual
        (primTile,
         primGlobal, primHalfStepGlobal, primOldGlobal,
         connectionGlobal,
         gradTGlobal, graduConGlobal, 
         graduConHigherOrderTerm1Global,
         graduConHigherOrderTerm2Global,
         ts->dt,
         ts->computeOldSourceTermsAndOldDivOfFluxes,
         ts->computeDivOfFluxAtTimeN,
         ts->computeDivOfFluxAtTimeNPlusHalf,
         iTile, jTile, X1Start, X2Start, X1Size, X2Size,
         residualGlobal
        );
    #endif


  } /* End of LOOP_OVER_TILES */


  #if (TIME_STEPPING==IMPLICIT)
    DMRestoreLocalVector(ts->dmdaWithGhostZones, &primPetscVecLocal);

    DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, primPetscVecLocal,
                           &primLocal);
  #endif

  DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, primPetscVec,
                         &primGlobal);
  DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVecOld,
                         &primOldGlobal); 
  DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, ts->primPetscVecHalfStep,
                         &primHalfStepGlobal);
  DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, residualPetscVec,
                         &residualGlobal);

  DMDAVecRestoreArrayDOF(ts->dmdaFluxes, ts->divFluxPetscVecOld,
                         &divFluxOldGlobal);
  DMDAVecRestoreArrayDOF(ts->dmdaFluxes, ts->sourceTermsPetscVecOld,
                         &sourceTermsOldGlobal); 
  DMDAVecRestoreArrayDOF(ts->dmdaFluxes, ts->conservedVarsPetscVecOld,
                         &conservedVarsOldGlobal); 
  
  DMDAVecRestoreArrayDOF(ts->connectionDMDA, ts->connectionPetscVec,
                         &connectionGlobal);
  DMDAVecRestoreArrayDOF(ts->dmdaDt, ts->dtPetscVec, &dtGlobal);

  #if (CONDUCTION)
    DMDAVecRestoreArrayDOF(ts->gradTDM, ts->gradTPetscVec, &gradTGlobal);
    DMDAVecRestoreArrayDOF(ts->graduConDM, ts->graduConPetscVec, 
                           &graduConGlobal);
    DMDAVecRestoreArrayDOF(ts->graduConHigherOrderTermsDM, 
                           ts->graduConHigherOrderTerm1PetscVec,
                           &graduConHigherOrderTerm1Global);
    DMDAVecRestoreArrayDOF(ts->graduConHigherOrderTermsDM, 
                           ts->graduConHigherOrderTerm2PetscVec,
                           &graduConHigherOrderTerm2Global);
  #endif

  #if (REAPER)
    DMDAVecRestoreArrayDOF(ts->momentsDM, ts->momentsVec, &momentsGlobal);
  #endif

  return(0);
}
