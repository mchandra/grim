#include "physics.h"

#if (REAPER)
  /* Public functions needed by the REAPER scheme in order to set initial 
   * conditions */
  REAL getAlpha(REAL rho, REAL temperature)
  {
    REAL besselK2;
    gsl_sf_result result;
    gsl_sf_bessel_Kn_e(2, 1./temperature, &result);
    besselK2 = result.val;
    return (rho/(4.*M_PI*temperature*besselK2));
  }
  
  REAL getA0(REAL temperature)
  {
    return (-1./temperature);
  }

  REAL getTemperature(const struct fluidElement elem[ARRAY_ARGS 1])
  {
    return (-1./elem->primVars[A0]);
  }

  REAL getDensity(const struct fluidElement elem[ARRAY_ARGS 1])
  {
    double temperature = getTemperature(elem);
    REAL besselK2;
    gsl_sf_result result;
    gsl_sf_bessel_Kn_e(2, 1./temperature, &result);
    besselK2 = result.val;
    return (elem->primVars[ALPHA]*4.*M_PI*temperature*besselK2);
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

  elem->uCon[1] =   elem->primVars[U1] 
                  - elem->gamma*geom->gCon[0][1]*geom->alpha;
  elem->uCon[2] =   elem->primVars[U2] 
                  - elem->gamma*geom->gCon[0][2]*geom->alpha;
  elem->uCon[3] =   elem->primVars[U3] 
                  - elem->gamma*geom->gCon[0][3]*geom->alpha;
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
    elem->bCon[i] = (  elem->primVars[B1+i-1] 
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
  #if (REAPER)
    setTetrad(geom, elem);
  #endif
  computeMoments(geom, elem);
}

void computeMoments(const struct geometry geom[ARRAY_ARGS 1],
                    struct fluidElement elem[ARRAY_ARGS 1])
{
#if (REAPER)
  REAL temperature = getTemperature(elem);

  /* When the temperature is very high, the particle velocity v approaches c
   * and the distribution function (as a function of the new scaled variables
   * t) is strongly peaked at 1. This can make the quadrature inaccurate.
   * Therefore we introduce a scale factor and change the transformation
   * depending on whether the dimensionless temp is < 1 (non-relativistic) or
   * > 1 (highly relativistic) 
   */
  REAL scaleFactor;
  if (temperature < 1)
  {
    scaleFactor = 1.;
  }
  else
  {
    scaleFactor = temperature;
  }

  REAL moments[NUM_ALL_COMPONENTS];
  fixedQuadIntegration(elem, geom, scaleFactor, moments);

  /* Now lower .._UP to .._DOWN (except for the number flux vector)*/
  for (int mu=0; mu<NDIM; mu++)
  {
    elem->moments[N_UP(mu)] = moments[N_UP(mu)];

    for (int nu=0; nu<NDIM; nu++)
    {
      elem->moments[T_UP_DOWN(mu, nu)] = 0.;

      for (int alpha=0; alpha<NDIM; alpha++)
      {
        elem->moments[T_UP_DOWN(mu, nu)] +=
	      moments[T_UP_UP(mu, alpha)]*geom->gCov[alpha][nu];
      }
      #if (REAPER_MOMENTS==15)
      	for (int lambda=0; lambda<NDIM; lambda++)
      	{
          elem->moments[M_UP_UP_DOWN(mu, nu, lambda)] = 0.;
      
          for (int alpha=0; alpha<NDIM; alpha++)
          {
            elem->moments[M_UP_UP_DOWN(mu, nu, lambda)] +=
                  moments[M_UP_UP_UP(mu, nu, alpha)]*geom->gCov[alpha][lambda];
          }
        }
      #endif
    }
  }
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
    }
  }
#endif

}

void computeFluxes(const struct fluidElement elem[ARRAY_ARGS 1],
                   const struct geometry geom[ARRAY_ARGS 1],
                   const int dir,
                   REAL fluxes[ARRAY_ARGS NUM_EQNS])
{
  REAL g = sqrt(-geom->gDet);
  #if (REAPER)
    fluxes[ALPHA] = g*elem->moments[N_UP(dir)];

    fluxes[A0] = g*elem->moments[T_UP_DOWN(dir, 0)] + fluxes[ALPHA];
    fluxes[U1] = g*elem->moments[T_UP_DOWN(dir, 1)];
    fluxes[U2] = g*elem->moments[T_UP_DOWN(dir, 2)];
    fluxes[U3] = g*elem->moments[T_UP_DOWN(dir, 3)];

    #if (REAPER_MOMENTS == 15)
      fluxes[B00] = g*elem->moments[M_UP_UP_DOWN(dir, 0, 0)];
      fluxes[B01] = g*elem->moments[M_UP_UP_DOWN(dir, 0, 1)]; 
      fluxes[B02] = g*elem->moments[M_UP_UP_DOWN(dir, 0, 2)]; 
      fluxes[B03] = g*elem->moments[M_UP_UP_DOWN(dir, 0, 3)]; 
      fluxes[B11] = g*elem->moments[M_UP_UP_DOWN(dir, 1, 1)]; 
      fluxes[B12] = g*elem->moments[M_UP_UP_DOWN(dir, 1, 2)]; 
      fluxes[B13] = g*elem->moments[M_UP_UP_DOWN(dir, 1, 3)]; 
      fluxes[B22] = g*elem->moments[M_UP_UP_DOWN(dir, 2, 2)]; 
      fluxes[B23] = g*elem->moments[M_UP_UP_DOWN(dir, 2, 3)]; 
      fluxes[B33] = g*elem->moments[M_UP_UP_DOWN(dir, 3, 3)]; 
    #endif

    fluxes[B1] = g*(elem->bCon[1]*elem->uCon[dir] - elem->bCon[dir]*elem->uCon[1]);
    fluxes[B2] = g*(elem->bCon[2]*elem->uCon[dir] - elem->bCon[dir]*elem->uCon[2]);
    fluxes[B3] = g*(elem->bCon[3]*elem->uCon[dir] - elem->bCon[dir]*elem->uCon[3]);
    
  #else

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

  #endif
}

/*  Returns sqrt(-gDet) * T^kappa^lamda * Gamma_lamda_nu_kappa for each nu.
 *  (Eqn (4) of HARM paper)
 * 
 */

void computeSourceTerms(const struct fluidElement elem[ARRAY_ARGS 1],
                        const struct geometry geom[ARRAY_ARGS 1],
                        const REAL christoffel[ARRAY_ARGS 64],
                        REAL sourceTerms[ARRAY_ARGS NUM_EQNS])
{
  /* Source terms not coded in for the REAPER scheme yet */

  for (int var=0; var<NUM_EQNS; var++)
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
