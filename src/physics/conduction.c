#include "physics.h"

void addConductionSourceTermsToResidual
(
  const REAL primTile[ARRAY_ARGS TILE_SIZE],
  ARRAY(primGlobal), ARRAY(primHalfStepGlobal), ARRAY(primOldGlobal),
  ARRAY(connectionGlobal),
  ARRAY(gradTGlobal), ARRAY(graduConGlobal), 
  ARRAY(graduConHigherOrderTerm1Global),
  ARRAY(graduConHigherOrderTerm2Global),
  REAL dt,
  int computeOldSourceTermsAndOldDivOfFluxes,
  int computeDivOfFluxAtTimeN,
  int computeDivOfFluxAtTimeNPlusHalf,
  const int iTile, const int jTile,
  const int X1Start, const int X2Start,
  const int X1Size, const int X2Size,
  ARRAY(residualGlobal)
)
{
#if (CONDUCTION)
  LOOP_INSIDE_TILE(0, TILE_SIZE_X1, 0, TILE_SIZE_X2)
  {
    struct gridZone zoneCenter;
    setGridZone(iTile, jTile,
                iInTile, jInTile,
                X1Start, X2Start, 
                X1Size, X2Size, 
                &zoneCenter);
    
    /* Values needed to compute spatial derivative.  Depending on where the
     * spatial derivatives are to be computed (t=n or t=n+1/2), the values of
     * primTile is set (already done, outside this routine) */

    if (computeOldSourceTermsAndOldDivOfFluxes)
    {
      REAL gradT[COMPUTE_DIM];
      REAL graduCon[COMPUTE_DIM*NDIM];
      REAL graduConHigherOrderTerm1[COMPUTE_DIM];
      REAL graduConHigherOrderTerm2[COMPUTE_DIM];
      computeConductionSpatialGradientTerms
        (primTile, 
         iTile, jTile, 
         iInTile, jInTile,
         X1Start, X2Start,
         X1Size, X2Size, 
         gradT, graduCon, 
         graduConHigherOrderTerm1,
         graduConHigherOrderTerm2
        );
      
      INDEX_PETSC(gradTGlobal, &zoneCenter, 0) = gradT[0];

      for (int mu=0; mu<NDIM; mu++)
      {
        INDEX_PETSC(graduConGlobal, &zoneCenter, mu) = graduCon[mu];
      }

      INDEX_PETSC(graduConHigherOrderTerm1Global, &zoneCenter, 0) = 
        graduConHigherOrderTerm1[0];

      INDEX_PETSC(graduConHigherOrderTerm2Global, &zoneCenter, 0) = 
        graduConHigherOrderTerm2[0];

      #if (COMPUTE_DIM==2)
        INDEX_PETSC(gradTGlobal, &zoneCenter, 1) = gradT[1];
        
        for (int mu=0; mu<NDIM; mu++)
        {
          INDEX_PETSC(graduConGlobal, &zoneCenter, mu+NDIM) = graduCon[mu+NDIM];
        }

        INDEX_PETSC(graduConHigherOrderTerm1Global, &zoneCenter, 1) = 
          graduConHigherOrderTerm1[1];

        INDEX_PETSC(graduConHigherOrderTerm2Global, &zoneCenter, 1) = 
          graduConHigherOrderTerm2[1];
      #endif
    } 
    else
    {
      #if (TIME_STEPPING==IMPLICIT)
        REAL gradT[COMPUTE_DIM];
        REAL graduCon[COMPUTE_DIM*NDIM];
        REAL graduConHigherOrderTerm1[COMPUTE_DIM];
        REAL graduConHigherOrderTerm2[COMPUTE_DIM];
        computeConductionSpatialGradientTerms
          (primTile, 
           iTile, jTile, 
           iInTile, jInTile,
           X1Start, X2Start,
           X1Size, X2Size, 
           gradT, graduCon,
           graduConHigherOrderTerm1,
           graduConHigherOrderTerm2
          );
      #endif

      /* Setup all the data structures needed 
       * Note: EXPLICIT and IMEX time stepping algorithms have a half step
       * involved whereas the IMPLICIT time step is a single step algorithm 
       * 
       * Legend: elem       --> fluid element at the (n+1)th time step 
       *
       *         elemCenter --> fluid element at the nth time step when in the
       *                        first  half of EXPLICIT/IMEX time step and
       *                        (n+1/2) time step when in the second half the
       *                        the EXPLICIT/IMEX time step. SHOULD NOT BE USED
       *                        WITH IMPLICIT time stepping since it is not set.
       *
       *         elemOld    --> fluid element at the nth time step.
       *
       *         zoneCenter --> zone at (iInTile, jInTile).
       *
       *         geomCenter --> geometry at the (iInTile, jInTile) zone.
       */
      struct fluidElement elem, elemCenter, elemOld;
      struct geometry geomCenter;

      REAL XCoords[NDIM];
      getXCoords(&zoneCenter, CENTER, XCoords);
      setGeometry(XCoords, &geomCenter);
      setFluidElement(&INDEX_PETSC(primGlobal, &zoneCenter, 0), &geomCenter,
                      &elem);

      setFluidElement(&INDEX_PETSC(primOldGlobal, &zoneCenter, 0), &geomCenter,
                      &elemOld);

      if (computeDivOfFluxAtTimeN)
      {
        setFluidElement(&INDEX_PETSC(primOldGlobal, &zoneCenter, 0), 
                        &geomCenter, &elemCenter);
      } 
      else if (computeDivOfFluxAtTimeNPlusHalf)
      {
        setFluidElement(&INDEX_PETSC(primHalfStepGlobal, &zoneCenter, 0),
                        &geomCenter, &elemCenter);
      }

      /* Values needed to compute time derivative */
      REAL T    =   (ADIABATIC_INDEX-1.)
                  * INDEX_PETSC(primGlobal, &zoneCenter, UU)
                  / INDEX_PETSC(primGlobal, &zoneCenter, RHO);

      REAL TOld =   (ADIABATIC_INDEX-1.)
                  * INDEX_PETSC(primOldGlobal, &zoneCenter, UU)
                  / INDEX_PETSC(primOldGlobal, &zoneCenter, RHO);

      #if (TIME_STEPPING == EXPLICIT || TIME_STEPPING == IMEX)

        REAL TCenter =   (ADIABATIC_INDEX-1.)
                       * elemCenter.primVars[UU]
                       / elemCenter.primVars[RHO];

      #elif (TIME_STEPPING == IMPLICIT)

        REAL TCenter = (T + TOld)/2.;

      #endif

      REAL dT[NDIM];

      dT[0] = (T - TOld)/dt;

      #if (TIME_STEPPING==EXPLICIT || TIME_STEPPING==IMEX)
        dT[1] = INDEX_PETSC(gradTGlobal, &zoneCenter, 0);

        dT[2] = 0.; dT[3] = 0.;
    
        #if (COMPUTE_DIM==2)
          dT[2] = INDEX_PETSC(gradTGlobal, &zoneCenter, 1);
        #endif

      #elif (TIME_STEPPING == IMPLICIT)
        dT[1] = 0.5*(  INDEX_PETSC(gradTGlobal, &zoneCenter, 0)
                     + gradT[0]
                    );

        dT[2] = 0.; dT[3] = 0.;
    
        #if (COMPUTE_DIM==2)
          dT[2] = 0.5*(  INDEX_PETSC(gradTGlobal, &zoneCenter, 1)
                       + gradT[1]
                      );
        #endif

      #endif /* TIME_STEPPING options for computing dT */

      /* Higher order term 1 
       * phi*(g*u^\mu)_{,\mu} */

      REAL dHigherOrderTerm1[NDIM];
      
      #if (TIME_STEPPING == EXPLICIT || TIME_STEPPING == IMEX)

        dHigherOrderTerm1[0] = 
          elemCenter.primVars[PHI] 
        * sqrt(-geomCenter.gDet) 
        * (elem.uCon[0] - elemOld.uCon[0])/dt;
      
        dHigherOrderTerm1[1] = 
          elemCenter.primVars[PHI]
        * INDEX_PETSC(graduConHigherOrderTerm1Global, &zoneCenter, 0);

        dHigherOrderTerm1[2] = 0.; dHigherOrderTerm1[3] = 0.;

        #if (COMPUTE_DIM==2)
          dHigherOrderTerm1[2] = 
            elemCenter.primVars[PHI]
          * INDEX_PETSC(graduConHigherOrderTerm1Global, &zoneCenter, 1);
        #endif

      #elif (TIME_STEPPING == IMPLICIT)

        dHigherOrderTerm1[0] = 
          0.5 * (elem.primVars[PHI] + elemOld.primVars[PHI])
              * sqrt(-geomCenter.gDet) 
              * (elem.uCon[0] - elemOld.uCon[0])/dt;
      
        dHigherOrderTerm1[1] = 
          0.5*(  (  elemOld.primVars[PHI]
                  * INDEX_PETSC(graduConHigherOrderTerm1Global, &zoneCenter, 0)
                 )
               + 
                 (  elem.primVars[PHI]
                  * graduConHigherOrderTerm1[0]
                 )
              );

        dHigherOrderTerm1[2] = 0.; dHigherOrderTerm1[3] = 0.;

        #if (COMPUTE_DIM==2)
          dHigherOrderTerm1[2] = 
          0.5*(  (  elemOld.primVars[PHI]
                  * INDEX_PETSC(graduConHigherOrderTerm1Global, &zoneCenter, 1)
                 )
               + 
                 (  elem.primVars[PHI]
                  * graduConHigherOrderTerm1[1]
                 )
              );
        #endif

      #endif /* TIME_STEPPING options for higherOrderTerm1 */

      /* Higher order term 2 
       * (phi*T)/(2*beta) * (beta*g*u^\mu/T)_{,\mu}
       * where beta = tau/kappa */

      REAL dHigherOrderTerm2[NDIM];
      
      REAL beta       = elem.tau/(elem.kappa*T);
      REAL betaOld    = elemOld.tau/(elem.kappa*TOld);
      #if (TIME_STEPPING == EXPLICIT || TIME_STEPPING == IMEX)

        REAL betaCenter = elemCenter.tau/(elemCenter.kappa*TCenter);

      #elif (TIME_STEPPING == IMPLICIT)

        REAL betaCenter = 0.5*(  (elem.tau/(elem.kappa*T) )
                               + (elemOld.tau/(elemOld.kappa*TOld) )
                              );  
      #endif

      #if (TIME_STEPPING == EXPLICIT || TIME_STEPPING == IMEX)

        dHigherOrderTerm2[0] = 
          (elemCenter.primVars[PHI] * TCenter)
        / (2.*betaCenter) 
        * sqrt(-geomCenter.gDet) 
        * ((beta*elem.uCon[0]/T) - (betaOld*elemOld.uCon[0]/TOld))/dt;
      
        dHigherOrderTerm2[1] = 
          (elemCenter.primVars[PHI] * TCenter)
        / (2.*betaCenter) 
        * INDEX_PETSC(graduConHigherOrderTerm2Global, &zoneCenter, 0);

        dHigherOrderTerm2[2] = 0.; dHigherOrderTerm2[3] = 0.;

        #if (COMPUTE_DIM==2)
          dHigherOrderTerm2[2] = 
            (elemCenter.primVars[PHI] * TCenter)
          / (2.*betaCenter) 
          * INDEX_PETSC(graduConHigherOrderTerm2Global, &zoneCenter, 1);
        #endif

      #elif (TIME_STEPPING == IMPLICIT)

        dHigherOrderTerm2[0] = 
          0.5*(  (elem.primVars[PHI] * T/(2.*beta) )
               + (elemOld.primVars[PHI] * TOld/(2.*betaOld) )
              )
        * sqrt(-geomCenter.gDet) 
        * ((beta*elem.uCon[0]/T) - (betaOld*elemOld.uCon[0]/TOld))/dt;


        dHigherOrderTerm2[1] = 
          0.5*(  (  elemOld.primVars[PHI] * TOld/(2.*betaOld) 
                  * INDEX_PETSC(graduConHigherOrderTerm2Global, &zoneCenter, 0)
                 )
               + 
                 (  elem.primVars[PHI] * T/(2.*beta)
                  * graduConHigherOrderTerm2[0]
                 )
              );

        dHigherOrderTerm2[2] = 0.; dHigherOrderTerm2[3] = 0.;

        #if (COMPUTE_DIM==2)
          dHigherOrderTerm2[2] = 
            0.5*(  (  elemOld.primVars[PHI] * TOld/(2.*betaOld) 
                    * INDEX_PETSC(graduConHigherOrderTerm2Global, &zoneCenter, 1)
                   )
                + 
                   (  elem.primVars[PHI] * T/(2.*beta)
                    * graduConHigherOrderTerm2[1]
                   )
                );
        #endif

      #endif

      REAL higherOrderTerm1 = 0.;
      REAL higherOrderTerm2 = 0.;
      for (int mu=0; mu<NDIM; mu++)
      {
        higherOrderTerm1 += dHigherOrderTerm1[mu];
        higherOrderTerm2 += dHigherOrderTerm2[mu];
      }

      REAL aCon[NDIM];

      for (int mu=0; mu<NDIM; mu++)
      {
        #if (TIME_STEPPING == EXPLICIT || TIME_STEPPING == IMEX)

          aCon[mu] =   elemCenter.uCon[0]
                     * (elem.uCon[mu] - elemOld.uCon[mu])/dt;
      
          aCon[mu] +=   elemCenter.uCon[1]
                      * INDEX_PETSC(graduConGlobal, &zoneCenter, mu);

          #if (COMPUTE_DIM==2)
            aCon[mu] +=   elemCenter.uCon[2]
                        * INDEX_PETSC(graduConGlobal, &zoneCenter, mu+NDIM);
          #endif

        #elif (TIME_STEPPING == IMPLICIT)

          aCon[mu] =  0.5 * (elem.uCon[0] + elemOld.uCon[0])
                          * (elem.uCon[mu] - elemOld.uCon[mu])/dt;
          
          aCon[mu] += 
            0.5*(  (elemOld.uCon[1]*INDEX_PETSC(graduConGlobal, &zoneCenter, mu) )
                 + (elem.uCon[1]*graduCon[mu])
                );

          #if (COMPUTE_DIM==2)
            aCon[mu] += 
              0.5*(  (elemOld.uCon[2]*INDEX_PETSC(graduConGlobal, &zoneCenter, mu+NDIM) )
                   + (elem.uCon[2]*graduCon[mu+NDIM])
                  );
          #endif

        #endif

        /* Add the connection terms now */
        #if (TIME_STEPPING == IMEX || TIME_STEPPING == IMPLICIT)
          for (int alpha=0; alpha<NDIM; alpha++)
          {
            for (int beta=0; beta<NDIM; beta++)
            {
              aCon[mu] += 
              INDEX_PETSC(connectionGlobal, &zoneCenter, GAMMA_UP_DOWN_DOWN(mu, alpha, beta))
              * (   (elem.uCon[alpha]    * elem.uCon[beta])
                  + (elemOld.uCon[alpha] * elemOld.uCon[beta])
                )/2.;
            }
          }
        #elif (TIME_STEPPING == EXPLICIT)
          for (int alpha=0; alpha<NDIM; alpha++)
          {
            for (int beta=0; beta<NDIM; beta++)
            {
              aCon[mu] +=  
              INDEX_PETSC(connectionGlobal, &zoneCenter, GAMMA_UP_DOWN_DOWN(mu, alpha, beta))
              * elemCenter.uCon[alpha] * elemCenter.uCon[beta];
            }
          }
        #endif
      }

      REAL aCov[NDIM];
      for (int mu=0; mu<NDIM; mu++)
      {
        aCov[mu] = 0.;
        for (int nu=0; nu<NDIM; nu++)
        {
          aCov[mu] += geomCenter.gCov[mu][nu]*aCon[nu];
        }
      }

      REAL qConEckart[NDIM];
      for (int mu=0; mu<NDIM; mu++)
      {
        qConEckart[mu] = 0.;
        for (int nu=0; nu<NDIM; nu++)
        {
          qConEckart[mu] +=
            -elemCenter.kappa
            * (  elemCenter.uCon[mu]*elemCenter.uCon[nu] 
               + geomCenter.gCon[mu][nu]
              )
            * (dT[nu] + TCenter*aCov[nu]);
        }
      }
      
      REAL bDotq = 0., bSqr, bCov[NDIM];
      conToCov(elemCenter.bCon, &geomCenter, bCov);
      bSqr = getbSqr(&elemCenter, &geomCenter);

      for (int mu=0; mu<NDIM; mu++)
      {
        bDotq += bCov[mu]*qConEckart[mu]/sqrt(bSqr);
      }
    
      
//      REAL P   = (ADIABATIC_INDEX-1.)*elemCenter.primVars[UU];
//      REAL cs  =  sqrt((ADIABATIC_INDEX-1.)*P
//                 /(elemCenter.primVars[RHO] + elemCenter.primVars[UU]));
//      REAL phiCeil = 0.1 * elemCenter.primVars[RHO] * pow(cs, 3.);
//
//      if (abs(bDotq) > phiCeil)
//      {
//        bDotq = copysign(phiCeil, bDotq);
//      }
//

      REAL g = sqrt(-geomCenter.gDet);
      REAL norm = g;
      #if (TIME_STEPPING == EXPLICIT)

//        INDEX_PETSC(residualGlobal, &zoneCenter, PHI) += 
//          (- higherOrderTerm1 + 0.*higherOrderTerm2
//           + g*( elemCenter.primVars[PHI]
//                 - bDotq
//               )/elemCenter.tau
//          )/norm;
        
        INDEX_PETSC(residualGlobal, &zoneCenter, PHI) *= elemCenter.tau;

        INDEX_PETSC(residualGlobal, &zoneCenter, PHI) +=
              (    (- higherOrderTerm1 + 0.*higherOrderTerm2)*elemCenter.tau
               + g*( elemCenter.primVars[PHI] - bDotq
                   )
              )/norm;
      #elif (TIME_STEPPING == IMEX || TIME_STEPPING == IMPLICIT)
  
//        INDEX_PETSC(residualGlobal, &zoneCenter, PHI) += 
//          (- higherOrderTerm1 + higherOrderTerm2
//           + g*( 0.5*(elem.primVars[PHI] + elemOld.primVars[PHI])
//                - bDotq
//               )/elemCenter.tau
//          )/norm;

        INDEX_PETSC(residualGlobal, &zoneCenter, PHI) *= elemCenter.tau;

        INDEX_PETSC(residualGlobal, &zoneCenter, PHI) +=
              (    (- higherOrderTerm1 + higherOrderTerm2)*elemCenter.tau
               + g*( 0.5*(elem.primVars[PHI] + elemOld.primVars[PHI])
                    - bDotq
                   )
              )/norm;

      #endif

    }
  }
#endif /* CONDUCTION */
}

void computeConductionSpatialGradientTerms
(
  const REAL primTile[ARRAY_ARGS TILE_SIZE],
  const int iTile, const int jTile,
  const int iInTile, const int jInTile,
  const int X1Start, const int X2Start,
  const int X1Size, const int X2Size,
  REAL gradT[COMPUTE_DIM],
  REAL graduCon[COMPUTE_DIM*NDIM],
  REAL graduConHigherOrderTerm1[COMPUTE_DIM],
  REAL graduConHigherOrderTerm2[COMPUTE_DIM]
)
{
#if (CONDUCTION)
  REAL XCoords[NDIM];
  struct gridZone zone, zoneCenter;
  struct geometry geom;
  struct fluidElement elem;
  REAL uConLeft[NDIM], uConCenter[NDIM], uConRight[NDIM];
  REAL gLeft, gCenter, gRight;
  REAL betaLeft, betaCenter, betaRight; /* Israel-Stewart's Beta = tau/(kappa*T) */

  setGridZone(iTile, jTile,
              iInTile, jInTile,
              X1Start, X2Start, 
              X1Size, X2Size, 
              &zoneCenter);

  REAL TCenter  =   (ADIABATIC_INDEX-1.)
                  * primTile[INDEX_TILE(&zoneCenter, UU)]
                  / primTile[INDEX_TILE(&zoneCenter, RHO)];

  REAL TLeftX1  =   (ADIABATIC_INDEX-1.)
                  * primTile[INDEX_TILE_MINUS_ONE_X1(&zoneCenter, UU)]
                  / primTile[INDEX_TILE_MINUS_ONE_X1(&zoneCenter, RHO)];

  REAL TRightX1 =  (ADIABATIC_INDEX-1.)
                  * primTile[INDEX_TILE_PLUS_ONE_X1(&zoneCenter, UU)]
                  / primTile[INDEX_TILE_PLUS_ONE_X1(&zoneCenter, RHO)];

  /* gradT[0] = dT_dX1 */
  gradT[0] = 
    slopeLimitedDerivative(TLeftX1, TCenter, TRightX1)/zoneCenter.dX1;

  #if (COMPUTE_DIM==2)
    REAL TLeftX2  =   (ADIABATIC_INDEX-1.)
                    * primTile[INDEX_TILE_MINUS_ONE_X2(&zoneCenter, UU)]
                    / primTile[INDEX_TILE_MINUS_ONE_X2(&zoneCenter, RHO)];

    REAL TRightX2 =   (ADIABATIC_INDEX-1.)
                    * primTile[INDEX_TILE_PLUS_ONE_X2(&zoneCenter, UU)]
                    / primTile[INDEX_TILE_PLUS_ONE_X2(&zoneCenter, RHO)];

    /* gradT[1] = dT_dX2 */
    gradT[1] = 
      slopeLimitedDerivative(TLeftX2, TCenter, TRightX2)/zoneCenter.dX2;
  #endif

  /* uConCenter */
  getXCoords(&zoneCenter, CENTER, XCoords);
  setGeometry(XCoords, &geom); gCenter = sqrt(-geom.gDet);
  setFluidElement(&primTile[INDEX_TILE(&zoneCenter, 0)], &geom, &elem);
  betaCenter = elem.tau/(elem.kappa*TCenter);
  for (int mu=0; mu<NDIM; mu++)
  {
    uConCenter[mu] = elem.uCon[mu];
  }

  /* uConLeft */
  setGridZone(iTile, jTile,
              iInTile-1, jInTile,
              X1Start, X2Start, 
              X1Size, X2Size, 
              &zone);
  getXCoords(&zone, CENTER, XCoords);
  setGeometry(XCoords, &geom); gLeft = sqrt(-geom.gDet);
  setFluidElement(&primTile[INDEX_TILE(&zone, 0)], &geom, &elem);
  betaLeft = elem.tau/(elem.kappa*TLeftX1);
  for (int mu=0; mu<NDIM; mu++)
  {
    uConLeft[mu] = elem.uCon[mu];
  }

  /* uConRight */
  setGridZone(iTile, jTile,
              iInTile+1, jInTile,
              X1Start, X2Start, 
              X1Size, X2Size, 
              &zone);
  getXCoords(&zone, CENTER, XCoords);
  setGeometry(XCoords, &geom); gRight = sqrt(-geom.gDet);
  setFluidElement(&primTile[INDEX_TILE(&zone, 0)], &geom, &elem);
  betaRight = elem.tau/(elem.kappa*TRightX1);
  for (int mu=0; mu<NDIM; mu++)
  {
    uConRight[mu] = elem.uCon[mu];
  }

  for (int mu=0; mu<NDIM; mu++)
  {
    graduCon[mu] = 
      slopeLimitedDerivative(uConLeft[mu],
                             uConCenter[mu],
                             uConRight[mu])/zoneCenter.dX1;
  }

  graduConHigherOrderTerm1[0] = 
    slopeLimitedDerivative(gLeft*uConLeft[1],
                           gCenter*uConCenter[1],
                           gRight*uConRight[1])/zoneCenter.dX1;

  graduConHigherOrderTerm2[0] = 
    slopeLimitedDerivative
      (betaLeft*gLeft*uConLeft[1]/TLeftX1,
       betaCenter*gCenter*uConCenter[1]/TCenter,
       betaRight*gRight*uConRight[1]/TRightX1
      )/zoneCenter.dX1;

  #if (COMPUTE_DIM==2)
    /* uConLeft */
    setGridZone(iTile, jTile,
                iInTile, jInTile-1,
                X1Start, X2Start, 
                X1Size, X2Size, 
                &zone);
    getXCoords(&zone, CENTER, XCoords);
    setGeometry(XCoords, &geom); gLeft = sqrt(-geom.gDet);
    setFluidElement(&primTile[INDEX_TILE(&zone, 0)], &geom, &elem);
    betaLeft = elem.tau/(elem.kappa*TLeftX2);
    for (int mu=0; mu<NDIM; mu++)
    {
      uConLeft[mu] = elem.uCon[mu];
    }

    /* uConRight */
    setGridZone(iTile, jTile,
                iInTile, jInTile+1,
                X1Start, X2Start, 
                X1Size, X2Size, 
                &zone);
    getXCoords(&zone, CENTER, XCoords);
    setGeometry(XCoords, &geom); gRight = sqrt(-geom.gDet);
    setFluidElement(&primTile[INDEX_TILE(&zone, 0)], &geom, &elem);
    betaRight = elem.tau/(elem.kappa*TRightX2);
    for (int mu=0; mu<NDIM; mu++)
    {
      uConRight[mu] = elem.uCon[mu];
    }

    for (int mu=0; mu<NDIM; mu++)
    {
      graduCon[mu+NDIM] = 
        slopeLimitedDerivative(uConLeft[mu],
                               uConCenter[mu],
                               uConRight[mu])/zoneCenter.dX2;
    }

    graduConHigherOrderTerm1[1] = 
      slopeLimitedDerivative(gLeft*uConLeft[2],
                             gCenter*uConCenter[2],
                             gRight*uConRight[2])/zoneCenter.dX2;

    graduConHigherOrderTerm2[1] = 
      slopeLimitedDerivative
        (betaLeft*gLeft*uConLeft[2]/TLeftX2,
         betaCenter*gCenter*uConCenter[2]/TCenter,
         betaRight*gRight*uConRight[2]/TRightX2
        )/zoneCenter.dX2;

  #endif

#endif /* CONDUCTION */
}
