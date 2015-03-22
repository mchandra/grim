#include "physics.h"

#if (REAPER)
  /* Public functions needed by the REAPER scheme in order to set initial 
   * conditions */
  REAL getAlpha(REAL rho, REAL temperature)
  {
    REAL besselK2;
    gsl_sf_result result;
    gsl_sf_bessel_Kn_e(2, 1./theta, &result);
    besselK2 = result.val;
    return (rho/(4.*M_PI*temperature*K2));
  }
  
  REAL getA0(REAL temperature)
  {
    return (1./temperature);
  }
#endif

void setGamma(const struct geometry geom[ARRAY_ARGS 1],
              struct fluidElement elem[ARRAY_ARGS 1])
{
  elem->gamma = 
    sqrt(1 + geom->gCov[1][1]*elem->primVars[U1]*elem->primVars[U1]
           + geom->gCov[2][2]*elem->primVars[U2]*elem->primVars[U2] 
           + geom->gCov[3][3]*elem->primVars[U3]*elem->primVars[U3]

         + 2*(  geom->gCov[1][2]*elem->primVars[U1]*elem->primVars[U2]
              + geom->gCov[1][3]*elem->primVars[U1]*elem->primVars[U3] 
              + geom->gCov[2][3]*elem->primVars[U2]*elem->primVars[U3]
             )
        );
}

void setUCon(const struct geometry geom[ARRAY_ARGS 1],
             struct fluidElement elem[ARRAY_ARGS 1])
{
  elem->uCon[0] = elem->gamma/geom->alpha;

  for (int i=1; i<NDIM; i++)
  {
    elem->uCon[i] =  elem->primVars[UU+i] 
                   - elem->gamma*geom->gCon[0][i]*geom->alpha;
  }
}

void setBCon(const struct geometry geom[ARRAY_ARGS 1],
             struct fluidElement elem[ARRAY_ARGS 1])
{
  REAL uCov[NDIM];

  conToCov(elem->uCon, geom, uCov);

  elem->bCon[0] =   elem->primVars[B1]*uCov[1]
                  + elem->primVars[B2]*uCov[2] 
                  + elem->primVars[B3]*uCov[3];
  
  for (int i=1; i<NDIM; i++)
  {
    elem->bCon[i] = (  elem->primVars[U3+i] 
                     + elem->bCon[0]*elem->uCon[i]
                    )/elem->uCon[0];
  }
}

void setFluidElement(const REAL primVars[ARRAY_ARGS DOF],
                     const struct geometry geom[ARRAY_ARGS 1],
                     struct fluidElement elem[ARRAY_ARGS 1])
{
  for (int var=0; var<DOF; var++)
  {
    elem->primVars[var] = primVars[var];
  }

  /* Need to be set in exactly the following order because of the dependecy
     structure */
  setGamma(geom, elem);
  setUCon(geom, elem);
  setBCon(geom, elem);
  #if (CONDUCTION)
    setConductionParameters(geom, elem);
  #endif
  computeMoments(geom, elem);
}

void computeMoments(const struct geometry geom[ARRAY_ARGS 1],
                    struct fluidElement elem[ARRAY_ARGS 1])
{
#if (REAPER && REAPER_MOMENTS==5)
  REAL theta = elem->primVars[PRESSURE]/elem->primVars[RHO];
  if (theta<1e-5) theta = 1e-5;

  /* When the temperature is very high, the particle velocity v approaches c
   * and the distribution function (as a function of the new scaled variables
   * t) is strongly peaked at 1. This can make the quadrature inaccurate.
   * Therefore we introduce a scale factor and change the transformation
   * depending on whether the dimensionless temp is < 1 (non-relativistic) or
   * > 1 (highly relativistic) 
   */
  REAL scaleFactor;
  if (theta < 1)
  {
    scaleFactor = 1.;
  }
  else
  {
    scaleFactor = theta;
  }

  /* Accurately calculate Bessel K2 (expensive) OUTSIDE of the integration! */
  REAL besselK2;
  gsl_sf_result result;
  gsl_sf_bessel_Kn_e(2, 1./theta, &result);
  besselK2 = result.val;

  fixedQuadIntegration5Moments(elem, geom, theta, besselK2, scaleFactor,
                               moments);

#else
  REAL pressure = (ADIABATIC_INDEX - 1.)*elem->primVars[UU];
  REAL bCov[NDIM], bSqr, uCov[NDIM];
  
  bSqr = getbSqr(elem, geom);

  conToCov(elem->uCon, geom, uCov);
  conToCov(elem->bCon, geom, bCov);

  for (int mu=0; mu<NDIM; mu++)
  {
    elem->moments[N_UP(mu)] = elem->primVars[RHO]*elem->uCon[mu];

    for (int nu=0; nu<NDIM; nu++)
    {
      elem->moments[T_UP_DOWN(mu,nu)] =   
                          (  elem->primVars[RHO] + elem->primVars[UU]
                           + pressure + bSqr
                          )*elem->uCon[mu]*uCov[nu]

                        + (pressure + 0.5*bSqr)*DELTA(mu, nu)

                        - elem->bCon[mu]*bCov[nu]
      #if (CONDUCTION) 
        + elem->primVars[PHI]/sqrt(bSqr)
        * (elem->uCon[mu]*bCov[nu] + elem->bCon[mu]*uCov[nu])
      #endif                  
                        ;

    }
  }
#endif

}

void computeFluxes(const struct fluidElement elem[ARRAY_ARGS 1],
                   const struct geometry geom[ARRAY_ARGS 1],
                   const int dir,
                   REAL fluxes[ARRAY_ARGS DOF])
{
  REAL g = sqrt(-geom->gDet);

  fluxes[RHO] = g*elem->moments[N_UP(dir)];

  fluxes[UU] = g*elem->moments[T_UP_DOWN(dir, 0)] + fluxes[RHO];
  fluxes[U1] = g*elem->moments[T_UP_DOWN(dir, 1)];
  fluxes[U2] = g*elem->moments[T_UP_DOWN(dir, 2)];
  fluxes[U3] = g*elem->moments[T_UP_DOWN(dir, 3)];

  fluxes[B1] = g*(elem->bCon[1]*elem->uCon[dir] - elem->bCon[dir]*elem->uCon[1]);
  fluxes[B2] = g*(elem->bCon[2]*elem->uCon[dir] - elem->bCon[dir]*elem->uCon[2]);
  fluxes[B3] = g*(elem->bCon[3]*elem->uCon[dir] - elem->bCon[dir]*elem->uCon[3]);

  #if (CONDUCTION)
    fluxes[PHI] = g*(elem->uCon[dir]*elem->primVars[PHI]);
  #endif
}

/*  Returns sqrt(-gDet) * T^kappa^lamda * Gamma_lamda_nu_kappa for each nu.
 *  (Eqn (4) of HARM paper)
 * 
 */

void computeSourceTerms(const struct fluidElement elem[ARRAY_ARGS 1],
                        const struct geometry geom[ARRAY_ARGS 1],
                        const REAL christoffel[ARRAY_ARGS 64],
                        REAL sourceTerms[ARRAY_ARGS DOF])
{
  for (int var=0; var<DOF; var++)
  {
    sourceTerms[var] = 0.;
  }

  #if (METRIC==KERRSCHILD)
    REAL g = sqrt(-geom->gDet);

    for (int nu=0; nu<NDIM; nu++)
    {
      for (int kappa=0; kappa<NDIM; kappa++)
      {
        for (int lamda=0; lamda<NDIM; lamda++)
        {
          sourceTerms[UU+nu] +=   
                g * elem->moments[T_UP_DOWN(kappa, lamda)]
                  * christoffel[GAMMA_UP_DOWN_DOWN(lamda, kappa, nu)];
        }
      }
    }
  #endif
}

REAL getbSqr(const struct fluidElement elem[ARRAY_ARGS 1],
             const struct geometry geom[ARRAY_ARGS 1])
{
  REAL bCov[NDIM], bSqr;
  conToCov(elem->bCon, geom, bCov);
  bSqr = covDotCon(bCov, elem->bCon);
  if (bSqr < 1e-25)
  {
    bSqr = 1e-25;
  }

  return bSqr;
}
