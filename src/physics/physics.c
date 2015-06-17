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
  fixedQuadIntegration(elem, geom, scaleFactor, moments, 
                       elem->collisionIntegrals);

  /* Now lower .._UP to .._DOWN for the stress-tensor to exactly conserve
   * angular momentum and energy in the Kerr metric. */
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
          elem->moments[M_UP_UP_UP(mu, nu, lambda)] = 
                moments[M_UP_UP_UP(mu, nu, lambda)];
        }
      #endif

      #if (REAPER_MOMENTS==35)
      	for (int lambda=0; lambda<NDIM; lambda++)
      	{
          elem->moments[M_UP_UP_UP(mu, nu, lambda)] = 
                moments[M_UP_UP_UP(mu, nu, lambda)];
      	  
          for (int eta=0; eta<NDIM; eta++)
      	  {
            elem->moments[R_UP_UP_UP_UP(mu, nu, lambda, eta)] = 
                  moments[R_UP_UP_UP_UP(mu, nu, lambda, eta)];
          }
        }
      #endif

    }
  }
#else
  /* Using canonical GRIM scheme */

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
#endif /* GRIM or REAPER? */

}

void computeFluxes(const struct fluidElement elem[ARRAY_ARGS 1],
                   const struct geometry geom[ARRAY_ARGS 1],
                   const int dir,
                   REAL fluxes[ARRAY_ARGS NUM_FLUXES])
{
  REAL g = sqrt(-geom->gDet);
  #if (REAPER)
    fluxes[ALPHA_FLUX] = g*elem->moments[N_UP(dir)];

    fluxes[A0_FLUX] = g*elem->moments[T_UP_DOWN(dir, 0)] + fluxes[ALPHA_FLUX];
    fluxes[U1_FLUX] = g*elem->moments[T_UP_DOWN(dir, 1)];
    fluxes[U2_FLUX] = g*elem->moments[T_UP_DOWN(dir, 2)];
    fluxes[U3_FLUX] = g*elem->moments[T_UP_DOWN(dir, 3)];

    #if (REAPER_MOMENTS == 15)
      fluxes[B01_FLUX] = g*elem->moments[M_UP_UP_UP(dir, 0, 1)]; 
      fluxes[B02_FLUX] = g*elem->moments[M_UP_UP_UP(dir, 0, 2)]; 
      fluxes[B03_FLUX] = g*elem->moments[M_UP_UP_UP(dir, 0, 3)]; 
      fluxes[B11_FLUX] = g*elem->moments[M_UP_UP_UP(dir, 1, 1)]; 
      fluxes[B12_FLUX] = g*elem->moments[M_UP_UP_UP(dir, 1, 2)]; 
      fluxes[B13_FLUX] = g*elem->moments[M_UP_UP_UP(dir, 1, 3)]; 
      fluxes[B22_FLUX] = g*elem->moments[M_UP_UP_UP(dir, 2, 2)]; 
      fluxes[B23_FLUX] = g*elem->moments[M_UP_UP_UP(dir, 2, 3)]; 
      fluxes[B33_FLUX] = g*elem->moments[M_UP_UP_UP(dir, 3, 3)]; 

      fluxes[F0_ALPHA_FLUX] = 0.;
      fluxes[F0_A0_FLUX]    = 0.;
      fluxes[F0_A1_FLUX]    = 0.;
      fluxes[F0_A2_FLUX]    = 0.;
      fluxes[F0_A3_FLUX]    = 0.;

    #elif (REAPER_MOMENTS == 35)

      fluxes[B01_FLUX] = g*elem->moments[M_UP_UP_UP(dir, 0, 1)]; 
      fluxes[B02_FLUX] = g*elem->moments[M_UP_UP_UP(dir, 0, 2)]; 
      fluxes[B03_FLUX] = g*elem->moments[M_UP_UP_UP(dir, 0, 3)]; 
      fluxes[B11_FLUX] = g*elem->moments[M_UP_UP_UP(dir, 1, 1)]; 
      fluxes[B12_FLUX] = g*elem->moments[M_UP_UP_UP(dir, 1, 2)]; 
      fluxes[B13_FLUX] = g*elem->moments[M_UP_UP_UP(dir, 1, 3)]; 
      fluxes[B22_FLUX] = g*elem->moments[M_UP_UP_UP(dir, 2, 2)]; 
      fluxes[B23_FLUX] = g*elem->moments[M_UP_UP_UP(dir, 2, 3)]; 
      fluxes[B33_FLUX] = g*elem->moments[M_UP_UP_UP(dir, 3, 3)]; 

      fluxes[C011_FLUX] = g*elem->moments[R_UP_UP_UP_UP(dir,0,1,1)];
      fluxes[C012_FLUX] = g*elem->moments[R_UP_UP_UP_UP(dir,0,1,2)];
      fluxes[C013_FLUX] = g*elem->moments[R_UP_UP_UP_UP(dir,0,1,3)];
      fluxes[C022_FLUX] = g*elem->moments[R_UP_UP_UP_UP(dir,0,2,2)];
      fluxes[C023_FLUX] = g*elem->moments[R_UP_UP_UP_UP(dir,0,2,3)];
      fluxes[C033_FLUX] = g*elem->moments[R_UP_UP_UP_UP(dir,0,3,3)];
      fluxes[C111_FLUX] = g*elem->moments[R_UP_UP_UP_UP(dir,1,1,1)];
      fluxes[C112_FLUX] = g*elem->moments[R_UP_UP_UP_UP(dir,1,1,2)];
      fluxes[C113_FLUX] = g*elem->moments[R_UP_UP_UP_UP(dir,1,1,3)];
      fluxes[C122_FLUX] = g*elem->moments[R_UP_UP_UP_UP(dir,1,2,2)];
      fluxes[C123_FLUX] = g*elem->moments[R_UP_UP_UP_UP(dir,1,2,3)];
      fluxes[C133_FLUX] = g*elem->moments[R_UP_UP_UP_UP(dir,1,3,3)];
      fluxes[C222_FLUX] = g*elem->moments[R_UP_UP_UP_UP(dir,2,2,2)];
      fluxes[C223_FLUX] = g*elem->moments[R_UP_UP_UP_UP(dir,2,2,3)];
      fluxes[C233_FLUX] = g*elem->moments[R_UP_UP_UP_UP(dir,2,3,3)];
      fluxes[C333_FLUX] = g*elem->moments[R_UP_UP_UP_UP(dir,3,3,3)];

      fluxes[F0_ALPHA_FLUX] = 0.;
      fluxes[F0_A0_FLUX]    = 0.;
      fluxes[F0_A1_FLUX]    = 0.;
      fluxes[F0_A2_FLUX]    = 0.;
      fluxes[F0_A3_FLUX]    = 0.;

    #endif

    fluxes[B1_FLUX] = g*(  elem->bCon[1]*elem->uCon[dir] 
                         - elem->bCon[dir]*elem->uCon[1]
                        );
    fluxes[B2_FLUX] = g*(  elem->bCon[2]*elem->uCon[dir] 
                         - elem->bCon[dir]*elem->uCon[2]
                        );
    fluxes[B3_FLUX] = g*(  elem->bCon[3]*elem->uCon[dir] 
                         - elem->bCon[dir]*elem->uCon[3]
                        );
  #else

    /* Using canonical GRIM scheme */
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
                        REAL sourceTerms[ARRAY_ARGS NUM_FLUXES])
{
  /* Source terms not coded in for the REAPER scheme yet */

  for (int var=0; var<NUM_FLUXES; var++)
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

  #if (REAPER && REAPER_MOMENTS==15)

    sourceTerms[F0_ALPHA_FLUX] = elem->collisionIntegrals[0];
    sourceTerms[F0_A0_FLUX]    = elem->collisionIntegrals[1];
    sourceTerms[F0_A1_FLUX]    = elem->collisionIntegrals[2];
    sourceTerms[F0_A2_FLUX]    = elem->collisionIntegrals[3];
    sourceTerms[F0_A3_FLUX]    = elem->collisionIntegrals[4];

    sourceTerms[B01_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_2(0,1)];
    sourceTerms[B02_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_2(0,2)];
    sourceTerms[B03_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_2(0,3)];
    sourceTerms[B11_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_2(1,1)];
    sourceTerms[B12_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_2(1,2)];
    sourceTerms[B13_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_2(1,3)];
    sourceTerms[B22_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_2(2,2)];
    sourceTerms[B23_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_2(2,3)];
    sourceTerms[B33_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_2(3,3)];

  #elif (REAPER && REAPER_MOMENTS==35)

    sourceTerms[F0_ALPHA_FLUX] = elem->collisionIntegrals[0];
    sourceTerms[F0_A0_FLUX]    = elem->collisionIntegrals[1];
    sourceTerms[F0_A1_FLUX]    = elem->collisionIntegrals[2];
    sourceTerms[F0_A2_FLUX]    = elem->collisionIntegrals[3];
    sourceTerms[F0_A3_FLUX]    = elem->collisionIntegrals[4];

    sourceTerms[B01_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_2(0,1)];
    sourceTerms[B02_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_2(0,2)];
    sourceTerms[B03_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_2(0,3)];
    sourceTerms[B11_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_2(1,1)];
    sourceTerms[B12_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_2(1,2)];
    sourceTerms[B13_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_2(1,3)];
    sourceTerms[B22_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_2(2,2)];
    sourceTerms[B23_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_2(2,3)];
    sourceTerms[B33_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_2(3,3)];

    sourceTerms[C011_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_3(0,1,1)];
    sourceTerms[C012_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_3(0,1,2)];
    sourceTerms[C013_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_3(0,1,3)];
    sourceTerms[C022_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_3(0,2,2)];
    sourceTerms[C023_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_3(0,2,3)];
    sourceTerms[C033_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_3(0,3,3)];
    sourceTerms[C111_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_3(1,1,1)];
    sourceTerms[C112_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_3(1,1,2)];
    sourceTerms[C113_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_3(1,1,3)];
    sourceTerms[C122_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_3(1,2,2)];
    sourceTerms[C123_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_3(1,2,3)];
    sourceTerms[C133_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_3(1,3,3)];
    sourceTerms[C222_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_3(2,2,2)];
    sourceTerms[C223_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_3(2,2,3)];
    sourceTerms[C233_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_3(2,3,3)];
    sourceTerms[C333_FLUX] = elem->collisionIntegrals[COLLISION_INTEGRAL_3(3,3,3)];

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
