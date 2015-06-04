#include "../problem.h"
#include "torus.h"

/* Calculate the constant angular momentum per unit inertial mass (l = u_phi *
 * u^t) for a given black hole spin and a radius of the accretion disk.  Eqn 3.8
 * of Fishbone and Moncrief, 1976 */
REAL lFishboneMoncrief(REAL a, REAL r, REAL theta)
{
  return sqrt(M/pow(r, 3.)) \
        *(  pow(r, 4.) + r*r*a*a - 2.*M*r*a*a \
          - a*sqrt(M*r)*(r*r - a*a) \
         )/ \
         (r*r - 3*M*r + 2.*a*sqrt(M*r));
}

REAL lnOfhTerm1(REAL a,
                REAL r, REAL theta, 
                REAL l)
{
  REAL Delta = computeDelta(a, r, theta);
  REAL Sigma = computeSigma(a, r, theta);
  REAL A     = computeA(a, r, theta);

  return 0.5*log( (1. + sqrt(1. + (4.*l*l*Sigma*Sigma*Delta)/ \
                                  (A*sin(theta)*A*sin(theta))
                            )
                  ) / (Sigma*Delta/A)
                );
}

REAL lnOfhTerm2(REAL a,
                REAL r, REAL theta, 
                REAL l)
{
  REAL Delta = computeDelta(a, r, theta);
  REAL Sigma = computeSigma(a, r, theta);
  REAL A     = computeA(a, r, theta);

  return -0.5*sqrt(1. + (4.*l*l*Sigma*Sigma*Delta) /
                        (A*A*sin(theta)*sin(theta))
                  );

}

REAL lnOfhTerm3(REAL a,
                REAL r, REAL theta, 
                REAL l)
{
  REAL A     = computeA(a, r, theta);

  return -2*a*M*r*l/A;
}

REAL computeDelta(REAL a, REAL r, REAL theta)
{
  return r*r - 2*M*r + a*a;
}

REAL computeSigma(REAL a, REAL r, REAL theta)
{
  return r*r + a*a*cos(theta)*cos(theta);
}

REAL computeA(REAL a, REAL r, REAL theta)
{
  REAL Delta = computeDelta(a, r, theta);

  return pow(r*r + a*a, 2.) - Delta*a*a*sin(theta)*sin(theta);
}


REAL computeLnOfh(REAL a, REAL r, REAL theta)
{
  REAL l = lFishboneMoncrief(a, R_PRESSURE_MAX, M_PI/2.);

  REAL term1 = lnOfhTerm1(a, r, theta, l);
  REAL term2 = lnOfhTerm2(a, r, theta, l);
  REAL term3 = lnOfhTerm3(a, r, theta, l);

  REAL term1InnerEdge = lnOfhTerm1(a, R_INNER_EDGE, M_PI/2., l);
  REAL term2InnerEdge = lnOfhTerm2(a, R_INNER_EDGE, M_PI/2., l);
  REAL term3InnerEdge = lnOfhTerm3(a, R_INNER_EDGE, M_PI/2., l);

  return  term1 + term2 + term3 \
        - term1InnerEdge - term2InnerEdge - term3InnerEdge;

}

void initialConditions(struct timeStepper ts[ARRAY_ARGS 1])
{
  PetscPrintf(PETSC_COMM_WORLD, "Initializing torus...");

  REAL randNum;
  PetscRandom randNumGen;
  PetscRandomCreate(PETSC_COMM_SELF, &randNumGen);
  PetscRandomSetType(randNumGen, PETSCRAND48);

  ARRAY(primOldGlobal);
  DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, 
                     ts->primPetscVecOld, &primOldGlobal);

  LOOP_OVER_TILES(ts->X1Size, ts->X2Size)
  {
    REAL primTile[TILE_SIZE];

    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);

      REAL XCoords[NDIM], xCoords[NDIM];
      getXCoords(&zone, CENTER, XCoords);
      XTox(XCoords, xCoords);

      /* The Fishbone-Moncrief solution is a stable hydrodynamic disk with a
       * constant angular momentum per unit inertial mass (l = u_phi * u^t). The
       * value of l is chosen by first choosing the radius at which the pressure
       * is maximum in the disk. Then the solution is characterized by two
       * free parameters:
       * 1) R_PRESSURE_MAX : The location of the pressure maximum of the disk.
       * 2) R_INNER_EDGE : The location of the inner edge of the disk.
       *
       * In addition, one has the choice of the adiabat C in the equation of
       * state: P = C rho^Gamma
       */

      REAL r = xCoords[1], theta = xCoords[2];

      REAL lnOfh;
      if (r >= R_INNER_EDGE)
      {
        lnOfh = computeLnOfh(BH_SPIN, r, theta);
      } else
      {
        lnOfh = 1.;
      }

      /* Region outside the torus */
      if (lnOfh < 0. || r < R_INNER_EDGE)
      {
        primTile[INDEX_TILE(&zone, RHO)] = RHO_FLOOR_MIN;
        primTile[INDEX_TILE(&zone, UU)] = UU_FLOOR_MIN;
        primTile[INDEX_TILE(&zone, U1)] = 0.;
        primTile[INDEX_TILE(&zone, U2)] = 0.;
        primTile[INDEX_TILE(&zone, U3)] = 0.;
      } else
      {
        REAL h = exp(lnOfh);

        /* Solve for rho using the definition of h = (rho + u + P)/rho where rho
         * here is the rest mass energy density and P = C * rho^Gamma */

        REAL rho = pow((h-1)*(ADIABATIC_INDEX-1.)/(ADIABAT*ADIABATIC_INDEX), 
                       1./(ADIABATIC_INDEX-1.)
                      );
        REAL u =  ADIABAT * pow(rho, ADIABATIC_INDEX)
                / (ADIABATIC_INDEX-1.);

        PetscRandomGetValue(randNumGen, &randNum);
        u = u*(1. + PERTURBATIONS_AMPLITUDE*(randNum-0.5));

        /* Fishbone-Moncrief u_phi is given in the Boyer-Lindquist coordinates.
         * Need to transform to (modified) Kerr-Schild */
        REAL A = computeA(BH_SPIN, r, theta);
        REAL Sigma = computeSigma(BH_SPIN, r, theta);
        REAL Delta = computeDelta(BH_SPIN, r, theta);
        REAL l = lFishboneMoncrief(BH_SPIN, R_PRESSURE_MAX, M_PI/2.);

			  REAL expOfMinus2Chi = Sigma*Sigma*Delta/(A*A*sin(theta)*sin(theta)) ;
        REAL uCovPhiBL = sqrt((-1. + sqrt(1. + 4*l*l*expOfMinus2Chi)
                               )/2.
                             );
        REAL uConPhiBL =   2.*BH_SPIN*r*sqrt(1. + uCovPhiBL*uCovPhiBL)
                          /sqrt(A*Sigma*Delta)
                        + sqrt(Sigma/A)*uCovPhiBL/sin(theta);

        REAL uConBL[NDIM];
        uConBL[0] = 0.;
        uConBL[1] = 0.;
        uConBL[2] = 0.;
        uConBL[3] = uConPhiBL;

        REAL gCovBL[NDIM][NDIM], gConBL[NDIM][NDIM];
        REAL transformBLToMKS[NDIM][NDIM];
      
        for (int alpha=0; alpha<NDIM; alpha++)
        {
          for (int beta=0; beta<NDIM; beta++)
          {
            gCovBL[alpha][beta] = 0.;
            gConBL[alpha][beta] = 0.;
            transformBLToMKS[alpha][beta] = 0.;
          }
        }

        REAL mu = 1 + BH_SPIN*BH_SPIN*cos(theta)*cos(theta)/(r*r);

        gCovBL[0][0] = -(1. - 2./(r*mu));
        gCovBL[0][3] = -2.*BH_SPIN*sin(theta)*sin(theta)/(r*mu);
        gCovBL[3][0] = gCovBL[0][3];
        gCovBL[1][1] = mu*r*r/Delta;
        gCovBL[2][2] = r*r*mu;
        gCovBL[3][3] = r*r*sin(theta)*sin(theta)*\
                       (1. + BH_SPIN*BH_SPIN/(r*r) +
                        2.*BH_SPIN*BH_SPIN*sin(theta)*sin(theta)/\
                        (r*r*r*mu)
                       );

        gConBL[0][0] = -1. -2.*(1 + BH_SPIN*BH_SPIN/(r*r))/(Delta*mu/r);
        gConBL[0][3] = -2.*BH_SPIN/(r*Delta*mu);
        gConBL[3][0] = gConBL[0][3];
        gConBL[1][1] = Delta/(r*r*mu);
        gConBL[2][2] = 1./(r*r*mu);
        gConBL[3][3] = (1. - 2./(r*mu))/(sin(theta)*sin(theta)*Delta);

        transformBLToMKS[0][0] = 1.;
        transformBLToMKS[1][1] = 1.;
        transformBLToMKS[2][2] = 1.;
        transformBLToMKS[3][3] = 1.;
        transformBLToMKS[0][1] = 2.*r/Delta;
        transformBLToMKS[3][1] = BH_SPIN/Delta; 

        /* Need to get uConBL[0] using u^mu u_mu = -1 */
        REAL AA = gCovBL[0][0];
        REAL BB = 2.*(gCovBL[0][1]*uConBL[1] +
                      gCovBL[0][2]*uConBL[2] +
                      gCovBL[0][3]*uConBL[3]
                     );
        REAL CC = 1. + gCovBL[1][1]*uConBL[1]*uConBL[1] +
                       gCovBL[2][2]*uConBL[2]*uConBL[2] +
                       gCovBL[3][3]*uConBL[3]*uConBL[3] +
                   2.*(gCovBL[1][2]*uConBL[1]*uConBL[2] +
                       gCovBL[1][3]*uConBL[1]*uConBL[3] +
                       gCovBL[2][3]*uConBL[2]*uConBL[3]);
  
        REAL discriminent = BB*BB - 4.*AA*CC;
        uConBL[0] = -(BB + sqrt(discriminent))/(2.*AA);

        REAL uConKS[NDIM];

        for (int alpha=0; alpha<NDIM; alpha++)
        {
          uConKS[alpha] = 0.;

          for (int beta=0; beta<NDIM; beta++)
          {
            uConKS[alpha] += transformBLToMKS[alpha][beta]*uConBL[beta];
          }
        }

        /* Finally get the four-velocity in the X coordinates, which is modified
         * Kerr-Schild */
        REAL uConMKS[NDIM];
        REAL rFactor = r;
        REAL hFactor = M_PI + (1. - H_SLOPE)*M_PI*cos(2.*M_PI*XCoords[2]);
        uConMKS[0] = uConKS[0];
        uConMKS[1] = uConKS[1]/rFactor;
        uConMKS[2] = uConKS[2]/hFactor;
        uConMKS[3] = uConKS[3];

        primTile[INDEX_TILE(&zone, RHO)] = rho;
        primTile[INDEX_TILE(&zone, UU)] = u;

        struct geometry geom; 
        setGeometry(XCoords, &geom);
        primTile[INDEX_TILE(&zone, U1)] =   
          uConMKS[1] + pow(geom.alpha, 2.)*geom.gCon[0][1]*uConMKS[0];
        primTile[INDEX_TILE(&zone, U2)] =   
          uConMKS[2] + pow(geom.alpha, 2.)*geom.gCon[0][2]*uConMKS[0];
        primTile[INDEX_TILE(&zone, U3)] =   
          uConMKS[3] + pow(geom.alpha, 2.)*geom.gCon[0][3]*uConMKS[0];


      }

      primTile[INDEX_TILE(&zone, B1)] = 0.;
      primTile[INDEX_TILE(&zone, B2)] = 0.;
      primTile[INDEX_TILE(&zone, B3)] = 0.;
     
      #if (CONDUCTION)
        primTile[INDEX_TILE(&zone, PHI)] = 0.;
      #endif
    }

    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);
      
      for (int var=0; var<DOF; var++)
      {
        INDEX_PETSC(primOldGlobal, &zone, var) =
          primTile[INDEX_TILE(&zone, var)];
      }
    }

  }

  REAL rhoMax;
  VecStrideMax(ts->primPetscVecOld, RHO, NULL, &rhoMax);

  LOOP_OVER_TILES(ts->X1Size, ts->X2Size)
  {
    REAL primTile[TILE_SIZE];

    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);

      for (int var=0; var<DOF; var++)
      {
        primTile[INDEX_TILE(&zone, var)] = 
          INDEX_PETSC(primOldGlobal, &zone, var);
      }

      primTile[INDEX_TILE(&zone, RHO)] = 
        primTile[INDEX_TILE(&zone, RHO)]/rhoMax;
      primTile[INDEX_TILE(&zone, UU)] =
        primTile[INDEX_TILE(&zone, UU)]/rhoMax;
    }

    applyFloor(iTile, jTile, 
               ts->X1Start, ts->X2Start,
               ts->X1Size, ts->X2Size,
               primTile);

    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);

      for (int var=0; var<DOF; var++)
      {
        INDEX_PETSC(primOldGlobal, &zone, var) =
          primTile[INDEX_TILE(&zone, var)];
      }
    }

  }

  Vec primPetscVecOldLocal, bSqrPetscVecGlobal;
  DMGetLocalVector(ts->dmdaWithGhostZones, &primPetscVecOldLocal);
  DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &bSqrPetscVecGlobal);

  ARRAY(primOldLocal);
  ARRAY(bSqrGlobal);
  DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, primPetscVecOldLocal,
                     &primOldLocal);
  DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones, bSqrPetscVecGlobal,
                     &bSqrGlobal);

  DMGlobalToLocalBegin(ts->dmdaWithGhostZones, 
                       ts->primPetscVecOld,
                       INSERT_VALUES,
                       primPetscVecOldLocal);
  DMGlobalToLocalEnd(ts->dmdaWithGhostZones,
                     ts->primPetscVecOld,
                     INSERT_VALUES,
                     primPetscVecOldLocal);

  /* Now set the magnetic field using the magnetic vector potential */
  LOOP_OVER_TILES(ts->X1Size, ts->X2Size)
  {
    REAL primTile[TILE_SIZE], AVectorTile[TILE_SIZE];

    LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);

      for (int var=0; var<DOF; var++)
      {
        primTile[INDEX_TILE(&zone, var)] = 
          INDEX_PETSC(primOldLocal, &zone, var);
      }
    }

    LOOP_INSIDE_TILE(-1, TILE_SIZE_X1+1, -1, TILE_SIZE_X2+1)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);

      REAL rhoAvg = 
        0.25*(  primTile[INDEX_TILE_MANUAL(zone.iInTile,   zone.jInTile,   RHO)]
              + primTile[INDEX_TILE_MANUAL(zone.iInTile-1, zone.jInTile,   RHO)]
              + primTile[INDEX_TILE_MANUAL(zone.iInTile,   zone.jInTile-1, RHO)]
              + primTile[INDEX_TILE_MANUAL(zone.iInTile-1, zone.jInTile-1, RHO)]
             );

      REAL AVec = rhoAvg - 0.2;

      AVectorTile[INDEX_TILE(&zone, 0)] = 0.;

      if (AVec > 0.)
      {
        AVectorTile[INDEX_TILE(&zone, 0)] = AVec;
      }
      
      /* Two loop initial magnetic field */
      REAL XCoords[NDIM], xCoords[NDIM];
      getXCoords(&zone, CENTER, XCoords);
      XTox(XCoords, xCoords);

      REAL theta = xCoords[2];
      AVectorTile[INDEX_TILE(&zone, 0)] = 
        cos(theta) * AVectorTile[INDEX_TILE(&zone, 0)];
    }

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
      struct geometry geom; setGeometry(XCoords, &geom);
      REAL g = sqrt(-geom.gDet);

      INDEX_PETSC(primOldGlobal, &zone, B1) = 
        -(  AVectorTile[INDEX_TILE_MANUAL(zone.iInTile,   zone.jInTile,   0)]
          - AVectorTile[INDEX_TILE_MANUAL(zone.iInTile,   zone.jInTile+1, 0)]
          + AVectorTile[INDEX_TILE_MANUAL(zone.iInTile+1, zone.jInTile,   0)]
          - AVectorTile[INDEX_TILE_MANUAL(zone.iInTile+1, zone.jInTile+1, 0)]
         )/(2.*zone.dX2*g);

      INDEX_PETSC(primOldGlobal, &zone, B2) = 
         (  AVectorTile[INDEX_TILE_MANUAL(zone.iInTile,   zone.jInTile,   0)]
          + AVectorTile[INDEX_TILE_MANUAL(zone.iInTile,   zone.jInTile+1, 0)]
          - AVectorTile[INDEX_TILE_MANUAL(zone.iInTile+1, zone.jInTile,   0)]
          - AVectorTile[INDEX_TILE_MANUAL(zone.iInTile+1, zone.jInTile+1, 0)]
         )/(2.*zone.dX1*g);

      INDEX_PETSC(primOldGlobal, &zone, B3) = 0.; 

      struct fluidElement elem;
      setFluidElement(&INDEX_PETSC(primOldGlobal, &zone, 0), &geom, &elem);

      REAL bCov[NDIM];
      conToCov(elem.bCon, &geom, bCov);
      REAL bSqr = covDotCon(bCov, elem.bCon);

      INDEX_PETSC(bSqrGlobal, &zone, 0) = bSqr;
    } 

  }

  DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones, bSqrPetscVecGlobal,
                         &bSqrGlobal);
  DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, primPetscVecOldLocal,
                         &primOldLocal);
  DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, 
                         ts->primPetscVecOld, &primOldGlobal);

  DMRestoreLocalVector(ts->dmdaWithGhostZones, &primPetscVecOldLocal);

  REAL bSqrMax, uMax;
  VecStrideMax(bSqrPetscVecGlobal, 0, NULL, &bSqrMax);
  VecStrideMax(ts->primPetscVecOld, UU, NULL, &uMax);

  REAL betaActual = (ADIABATIC_INDEX-1.)*uMax/(0.5*bSqrMax);
  REAL norm = sqrt(betaActual/PLASMA_BETA);

  VecStrideScale(ts->primPetscVecOld, B1, norm);
  VecStrideScale(ts->primPetscVecOld, B2, norm);
  VecStrideScale(ts->primPetscVecOld, B3, norm);

  VecDestroy(&bSqrPetscVecGlobal);
  PetscRandomDestroy(&randNumGen);

  PetscPrintf(PETSC_COMM_WORLD, "done\n");
  /* Done with setting the initial conditions */

  PetscPrintf(PETSC_COMM_WORLD, "Dumping initial fluxes and source terms...");

  /* Compute and output the fluxes for the initial data */
  Vec fluxX1PetscVec, fluxX2PetscVec, sourcesPetscVec;
  DMCreateGlobalVector(ts->dmdaWithGhostZones, &fluxX1PetscVec);
  DMCreateGlobalVector(ts->dmdaWithGhostZones, &fluxX2PetscVec);
  DMCreateGlobalVector(ts->dmdaWithoutGhostZones, &sourcesPetscVec);

  ARRAY(fluxX1Global);
  ARRAY(fluxX2Global);
  ARRAY(sourcesGlobal);
  ARRAY(connectionGlobal);
  ARRAY(dtGlobal);

  DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, 
                     fluxX1PetscVec, &fluxX1Global);
  DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, 
                     fluxX2PetscVec, &fluxX2Global);
  DMDAVecGetArrayDOF(ts->dmdaWithoutGhostZones,
                     sourcesPetscVec, &sourcesGlobal);
  DMDAVecGetArrayDOF(ts->connectionDMDA, ts->connectionPetscVec,
                     &connectionGlobal);
  DMDAVecGetArrayDOF(ts->dmdaDt, ts->dtPetscVec, &dtGlobal);

  Vec primPetscVecLocal;
  DMGetLocalVector(ts->dmdaWithGhostZones, &primPetscVecLocal);

  DMGlobalToLocalBegin(ts->dmdaWithGhostZones, 
                       ts->primPetscVecOld,
                       INSERT_VALUES,
                       primPetscVecLocal);
  DMGlobalToLocalEnd(ts->dmdaWithGhostZones,
                     ts->primPetscVecOld,
                     INSERT_VALUES,
                     primPetscVecLocal);

  ARRAY(primLocal);
  DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, primPetscVecLocal,
                     &primLocal);

  LOOP_OVER_TILES(ts->X1Size, ts->X2Size)
  {
    REAL primTile[TILE_SIZE];
    REAL fluxX1Tile[TILE_SIZE], fluxX2Tile[TILE_SIZE];

    LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start, 
                  ts->X1Size, ts->X2Size, 
                  &zone);
      for (int var=0; var<DOF; var++)
      {
        primTile[INDEX_TILE(&zone, var)] =
        INDEX_PETSC(primLocal, &zone, var);
      }
    }

    applyTileBoundaryConditions(iTile, jTile,
                                ts->X1Start, ts->X2Start,
                                ts->X1Size, ts->X2Size,
                                primTile);

    computeFluxesOverTile(primTile,
                          iTile, jTile,
                          ts->X1Start, ts->X2Start,
                          ts->X1Size, ts->X2Size,
                          fluxX1Tile, fluxX2Tile,
                          dtGlobal);

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
      struct geometry geom; setGeometry(XCoords, &geom);
      struct fluidElement elem;
      setFluidElement(&primTile[INDEX_TILE(&zone, 0)], &geom, &elem);

      REAL sourceTerms[DOF];
      computeSourceTerms(&elem, &geom,
                         &INDEX_PETSC(connectionGlobal, &zone, 0),
                         sourceTerms);

      for (int var=0; var<DOF; var++)
      {
        INDEX_PETSC(fluxX1Global, &zone, var) =
          fluxX1Tile[INDEX_TILE(&zone, var)];

        INDEX_PETSC(fluxX2Global, &zone, var) =
          fluxX2Tile[INDEX_TILE(&zone, var)];

        INDEX_PETSC(sourcesGlobal, &zone, var) =
          sourceTerms[var];
      }
    }
  }

  DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones,
                         primPetscVecLocal, &primLocal);
  DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, 
                         fluxX1PetscVec, &fluxX1Global);
  DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, 
                         fluxX2PetscVec, &fluxX2Global);
  DMDAVecRestoreArrayDOF(ts->dmdaWithoutGhostZones,
                         sourcesPetscVec, &sourcesGlobal);
  DMDAVecRestoreArrayDOF(ts->connectionDMDA, ts->connectionPetscVec,
                         &connectionGlobal);
  DMDAVecRestoreArrayDOF(ts->dmdaDt, ts->dtPetscVec, &dtGlobal);

  char fluxesX1FileName[50], fluxesX2FileName[50], sourcesFileName[50];
  sprintf(fluxesX1FileName, "%s.h5", "fluxX1Initial");
  sprintf(fluxesX2FileName, "%s.h5", "fluxX2Initial");
  sprintf(sourcesFileName, "%s.h5", "sourcesInitial");

  PetscViewer viewer;
  PetscViewerHDF5Open(PETSC_COMM_WORLD, fluxesX1FileName,
                      FILE_MODE_WRITE, &viewer);
  PetscObjectSetName((PetscObject) fluxX1PetscVec, "fluxX1");
  VecView(fluxX1PetscVec, viewer);
  PetscViewerDestroy(&viewer);

  PetscViewerHDF5Open(PETSC_COMM_WORLD, fluxesX2FileName,
                      FILE_MODE_WRITE, &viewer);
  PetscObjectSetName((PetscObject) fluxX2PetscVec, "fluxX2");
  VecView(fluxX2PetscVec, viewer);
  PetscViewerDestroy(&viewer);

  PetscViewerHDF5Open(PETSC_COMM_WORLD, sourcesFileName,
                      FILE_MODE_WRITE, &viewer);
  PetscObjectSetName((PetscObject) sourcesPetscVec, "sources");
  VecView(sourcesPetscVec, viewer);
  PetscViewerDestroy(&viewer);

  VecDestroy(&fluxX1PetscVec);
  VecDestroy(&fluxX2PetscVec); 
  VecDestroy(&sourcesPetscVec); 

  PetscPrintf(PETSC_COMM_WORLD, "done\n");
}

void applyFloorInsideNonLinearSolver(const int iTile, const int jTile,
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

    if (primTile[INDEX_TILE(&zone, RHO)] < RHO_FLOOR_MIN)
    {
      primTile[INDEX_TILE(&zone, RHO)] = RHO_FLOOR_MIN;
    }

    if (primTile[INDEX_TILE(&zone, UU)] < UU_FLOOR_MIN)
    {
      primTile[INDEX_TILE(&zone, UU)] = UU_FLOOR_MIN;
    }
  }
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

    REAL XCoords[NDIM], xCoords[NDIM];
    getXCoords(&zone, CENTER, XCoords);
    XTox(XCoords, xCoords);

    REAL r = xCoords[1];

    REAL rhoFloor = RHO_FLOOR*pow(r, RHO_FLOOR_FALLOFF);
    REAL uFloor = UU_FLOOR*pow(r, UU_FLOOR_FALLOFF);

    if (rhoFloor < RHO_FLOOR_MIN)
    {
      rhoFloor = RHO_FLOOR_MIN;
    }

    if (uFloor < UU_FLOOR_MIN)
    {
      uFloor = UU_FLOOR_MIN;
    }

    REAL rho = primTile[INDEX_TILE(&zone, RHO)];
    REAL u   = primTile[INDEX_TILE(&zone, UU)];

    if (rho < rhoFloor)
    {
      primTile[INDEX_TILE(&zone, RHO)] = rhoFloor;
    }

    if (u < uFloor)
    {
      primTile[INDEX_TILE(&zone, UU)] = uFloor;
    }

    #if (CONDUCTION)
      REAL P   = (ADIABATIC_INDEX-1.)*u;
      REAL cs  = sqrt((ADIABATIC_INDEX-1.)*P/(rho + u));
      REAL phi = primTile[INDEX_TILE(&zone, PHI)];

      REAL phiCeil = 0.9 * rho * pow(cs, 3.);

      if (abs(phi) > phiCeil)
      {
        primTile[INDEX_TILE(&zone, PHI)] = copysign(phiCeil, phi);
      }
    #endif

    struct geometry geom;
    setGeometry(XCoords, &geom);

    struct fluidElement elem;
    setFluidElement(&primTile[INDEX_TILE(&zone, 0)], &geom, &elem);

    if (elem.gamma > GAMMA_MAX)
    {
      REAL factor = sqrt( (GAMMA_MAX*GAMMA_MAX-1.)
                         /(elem.gamma*elem.gamma-1.)
                        );

      primTile[INDEX_TILE(&zone, U1)] *= factor;
      primTile[INDEX_TILE(&zone, U2)] *= factor;
      primTile[INDEX_TILE(&zone, U3)] *= factor;
    }
  }

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

  LOOP_INSIDE_TILE(-NG, TILE_SIZE_X1+NG, -NG, TILE_SIZE_X2+NG)
  {
    struct gridZone zone;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start,
                X1Size, X2Size,
                &zone);

    /* Inflow check at the inner and outer radial boundary */
    if (zone.i < 0 || zone.i > N1-1)
    {
      inflowCheck(&zone, primTile);
    }

    if (zone.j < 0 || zone.j > N2-1)
    {
      primTile[INDEX_TILE(&zone, U2)] *= -1;
      primTile[INDEX_TILE(&zone, B2)] *= -1;
    }
  }

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

  LOOP_INSIDE_TILE(-2, TILE_SIZE_X1+2, -2, TILE_SIZE_X2+2)
  {
    struct gridZone zone;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start,
                X1Size, X2Size,
                &zone);

    /* Set flux on the polar boundaries to zero 
     *
     *  index = 0      N2-1  N2
     *  |              |     |
     *  v              v     v
     *  |  o  |        |  o  | */

    if (zone.j == 0 || zone.j == N2) 
    {
      for (int var=0; var<DOF; var++)
      {
        fluxX2Tile[INDEX_TILE(&zone, var)] = 0.;
      }
    }

    if (zone.j == 0)
    {
      fluxX1Tile[INDEX_TILE_MINUS_ONE_X2(&zone, B2)] = 
        -fluxX1Tile[INDEX_TILE(&zone, B2)];
    }

    if (zone.j == N2)
    {
      fluxX1Tile[INDEX_TILE(&zone, B2)] = 
        -fluxX1Tile[INDEX_TILE_MINUS_ONE_X2(&zone, B2)];
    }

    /* Make sure there is no inflow from the radial boundaries */
    if (zone.i==0) /* Inner radial boundary */
    {
      if (fluxX1Tile[INDEX_TILE(&zone, RHO)] > 0.)
      {
        fluxX1Tile[INDEX_TILE(&zone, RHO)] = 0.;
      }
    }

    /* Outer radial boundary. Note that the index is N1 and not N1-1 */
    if (zone.i==N1)     
    {
      if (fluxX1Tile[INDEX_TILE(&zone, RHO)] < 0.)
      {
        fluxX1Tile[INDEX_TILE(&zone, RHO)] = 0.;
      }
    }

  }

}

void halfStepDiagnostics(struct timeStepper ts[ARRAY_ARGS 1])
{
  ARRAY(primHalfStepGlobal);
  DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, 
                     ts->primPetscVecHalfStep, &primHalfStepGlobal);

  #if (USE_OPENMP)
    #pragma omp parallel for
  #endif
  LOOP_OVER_TILES(ts->X1Size, ts->X2Size)
  {
    REAL primTile[TILE_SIZE];

    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);

      for (int var=0; var<DOF; var++)
      {
        primTile[INDEX_TILE(&zone, var)] =
        INDEX_PETSC(primHalfStepGlobal, &zone, var);
      }
    }

    applyFloor(iTile, jTile,
               ts->X1Start, ts->X2Start,
               ts->X1Size, ts->X2Size,
               primTile);

    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);

      for (int var=0; var<DOF; var++)
      {
        INDEX_PETSC(primHalfStepGlobal, &zone, var) =
        primTile[INDEX_TILE(&zone, var)];
      }
    }

  }

  DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, 
                         ts->primPetscVecHalfStep, &primHalfStepGlobal);
}

void fullStepDiagnostics(struct timeStepper ts[ARRAY_ARGS 1])
{
  /* Stricty speaking, need to use ts->dmdaWithoutGhostZones for EXPLICIT, IMEX
   * and ts->dmdaWithGhostZones for IMPLICIT. I'm not going to bother coding
   * that right now but this is a note to make sure you don't access the ghost
   * zones when using EXPLICT and IMEX */
  ARRAY(primGlobal);
  DMDAVecGetArrayDOF(ts->dmdaWithGhostZones, 
                     ts->primPetscVec, &primGlobal);

  #if (USE_OPENMP)
    #pragma omp parallel for
  #endif
  LOOP_OVER_TILES(ts->X1Size, ts->X2Size)
  {
    REAL primTile[TILE_SIZE];

    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);

      for (int var=0; var<DOF; var++)
      {
        primTile[INDEX_TILE(&zone, var)] =
        INDEX_PETSC(primGlobal, &zone, var);
      }
    }

    applyFloor(iTile, jTile,
               ts->X1Start, ts->X2Start,
               ts->X1Size, ts->X2Size,
               primTile);

    LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  ts->X1Start, ts->X2Start,
                  ts->X1Size, ts->X2Size,
                  &zone);

      for (int var=0; var<DOF; var++)
      {
        INDEX_PETSC(primGlobal, &zone, var) =
        primTile[INDEX_TILE(&zone, var)];
      }
    }

  }

  DMDAVecRestoreArrayDOF(ts->dmdaWithGhostZones, 
                         ts->primPetscVec, &primGlobal);
}

void inflowCheck(const struct gridZone zone[ARRAY_ARGS 1],
                 REAL primTile[ARRAY_ARGS TILE_SIZE])
{
  REAL XCoords[NDIM];
  getXCoords(zone, CENTER, XCoords);
  struct geometry geom;
  setGeometry(XCoords, &geom);

  struct fluidElement elem;
  setFluidElement(&primTile[INDEX_TILE(zone, 0)], &geom, &elem);

  int inflowInnerBoundary = (zone->i < 0    && elem.uCon[1] > 0.);
  int inflowOuterBoundary = (zone->i > N1-1 && elem.uCon[1] < 0.);

  if (inflowInnerBoundary || inflowOuterBoundary)
  {
    /* Remove gamma from the primitives. Reason? */
    primTile[INDEX_TILE(zone, U1)] /= elem.gamma;
    primTile[INDEX_TILE(zone, U2)] /= elem.gamma;
    primTile[INDEX_TILE(zone, U3)] /= elem.gamma;

    /* Set new radial velocity to be zero */
    primTile[INDEX_TILE(zone, U1)] = geom.gCon[0][1]*geom.alpha;

    /* Now calculate the new gamma corresponding to the new primitives */
    REAL vSqr = 0.;
    for (int a=1; a<NDIM; a++)
    {
      for (int b=1; b<NDIM; b++)
      {
        vSqr +=  primTile[INDEX_TILE(zone, UU+a)]
               * primTile[INDEX_TILE(zone, UU+b)]
               * geom.gCov[a][b];
      }
    }

    if (fabs(vSqr) < 1e-13) vSqr = 1e-13;
    if (vSqr >= 1.) vSqr = 1. - 1./(GAMMA_MAX*GAMMA_MAX);

    REAL gamma = 1./sqrt(1. - vSqr);

    /* Put this new gamma into the primitives */
    primTile[INDEX_TILE(zone, U1)] *= gamma;
    primTile[INDEX_TILE(zone, U2)] *= gamma;
    primTile[INDEX_TILE(zone, U3)] *= gamma;

  }

}

#if (CONDUCTION)
void setConductionParameters(const struct geometry geom[ARRAY_ARGS 1],
                             struct fluidElement elem[ARRAY_ARGS 1])
{
  REAL xCoords[NDIM];
  XTox(geom->XCoords, xCoords);

  REAL r = xCoords[1];

  REAL Rho = elem->primVars[RHO];
  if(Rho<RHO_FLOOR_MIN)
     Rho=RHO_FLOOR_MIN;
  REAL U = elem->primVars[UU];
  if(U<UU_FLOOR_MIN)
     U = UU_FLOOR_MIN;

  REAL P   = (ADIABATIC_INDEX-1.)*U;
  REAL cs  = sqrt(  (ADIABATIC_INDEX-1.)*P
                  / (Rho + U)
                 );

  REAL phiCeil = Rho * pow(cs, 3.);

  REAL tauDynamical = pow(r, 3./2.);
  REAL lambda       = 0.1;
  //REAL y            = elem->primVars[PHI]/phiCeil;
  //REAL fermiDirac   = 1./(exp((y-1.)/lambda) + 1.) + 1e-10;

  REAL y = fabs(elem->primVars[PHI])/fabs(phiCeil);
  y = (y-1)/lambda;
  REAL fermiDirac = exp(-y)/(exp(-y) + 1.)+1.e-10;
    
  REAL tau    = tauDynamical*fermiDirac;
  elem->kappa = 10.*cs*cs*tau*elem->primVars[RHO];
  elem->tau   = tau + 0.01; 

//  elem->kappa = 0.2 * pow(r, 0.5) * fabs(elem->primVars[RHO]);
//  elem->tau   = 1.2 * fabs(elem->kappa/elem->primVars[RHO]/ T) + 0.1;  
}
#endif
