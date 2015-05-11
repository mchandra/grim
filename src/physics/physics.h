#ifndef GRIM_PHYSICS_H_
#define GRIM_PHYSICS_H_

#include "../inputs.h"

/* Primitive variable mnemonics */
#if (REAPER)
  #if (REAPER_MOMENTS==5)

    #define ALPHA           (0)
    #define A0              (1) /* A0 = -1/temperature */
    #define U1              (2)
    #define U2              (3)
    #define U3              (4)
    #define B1              (5)
    #define B2              (6)
    #define B3              (7)
    #define DOF             (8)
    #define NUM_FLUXES      (DOF) /* NUM_FLUXES same as DOF when not using GYROAVERAGING */

  #elif (REAPER_MOMENTS==15)

    #define ALPHA_FLUX      (0)
    #define A0_FLUX         (1)
    #define U1_FLUX         (2)
    #define U2_FLUX         (3)
    #define U3_FLUX         (4)
    #define B01_FLUX        (5)
    #define B02_FLUX        (6)
    #define B03_FLUX        (7)
    #define B11_FLUX        (8)
    #define B12_FLUX        (9)
    #define B13_FLUX        (10)
    #define B22_FLUX        (11)
    #define B23_FLUX        (12)
    #define B33_FLUX        (13)
    #define B1_FLUX         (14)
    #define B2_FLUX         (15)
    #define B3_FLUX         (16)
    #define F0_ALPHA_FLUX   (17)
    #define F0_A0_FLUX      (18)
    #define F0_A1_FLUX      (19)
    #define F0_A2_FLUX      (20)
    #define F0_A3_FLUX      (21)
    #define NUM_FLUXES      (22)

    #if (GYROAVERAGING)
      /* 1 --> Parallel
         2 --> Perpendicular

         When GYROAVERAGING, number of primitive variables is lesser than the
         number of fluxes. So we have seperate mnemonics for each. First the
         primitive variables and then the fluxes.
      */
      #define ALPHA           (0)
      #define A0              (1)
      #define U1              (2)
      #define U2              (3)
      #define U3              (4)
      #define B01             (5)
      #define B02             (6)
      #define B11             (7)
      #define B12             (8)
      #define B22             (9)
      #define B1              (10)
      #define B2              (11)
      #define B3              (12)
      #define DOF             (13)

    #else
      /* Not GYROAVERAGING. Number of primitive variables == Number of fluxes */
      #define ALPHA           (0)
      #define A0              (1)
      #define U1              (2)
      #define U2              (3)
      #define U3              (4)
      #define B01             (5)
      #define B02             (6)
      #define B03             (7)
      #define B11             (8)
      #define B12             (9)
      #define B13             (10)
      #define B22             (11)
      #define B23             (12)
      #define B33             (13)
      #define B1              (14)
      #define B2              (15)
      #define B3              (16)
      #define F0_ALPHA        (17)
      #define F0_A0           (18)
      #define F0_A1           (19)
      #define F0_A2           (20)
      #define F0_A3           (21)
      #define DOF             (22)

    #endif /* GYROAVERAGING? */
  
  #endif /* REAPER 15 moments */
  
#else

  /* Not using REAPER scheme */

  #define RHO             (0)
  #define UU              (1)
  #define U1              (2)
  #define U2              (3)
  #define U3              (4)
  #define B1              (5)
  #define B2              (6)
  #define B3              (7)
  #if (CONDUCTION)
    #define PHI           (8)
    #define DOF           (9)
  #else
    #define DOF           (8)
  #endif
  #define NUM_FLUXES      (DOF)

#endif

/* Contains all the variables needed for physics. Independent variables are only
 * primVars. The rest are auxiliary variables stored for convenience.
 */
#if (REAPER_MOMENTS==5 || (REAPER==OFF) )
  /*   4  components of the number flux vector
   * + 16 components of the stress tensor.
   * Both grim and reaper 5 moment scheme have the same number of components */
  #define NUM_ALL_COMPONENTS  (20)
  #define NUM_COLLISION_PARAMETERS    (0)
#elif (REAPER_MOMENTS==15)
  /*   4  components of the number flux vector
   * + 16 components of the stress tensor
   * + 64 components of the rank 3 moment of f */
  #define NUM_ALL_COMPONENTS          (84) 
  #define NUM_COLLISION_PARAMETERS    (5)
  #define NUM_ALL_COLLISION_INTEGRALS (21) /* 5 parameters + 16 */
#endif

#define DELTA(mu, nu) (mu==nu ? 1 : 0)

/* Indices for elem.moments[] */
#define N_UP(mu) (mu)
#define T_UP_DOWN(mu, nu) (nu + NDIM*(mu) + NDIM)
#define T_UP_UP(mu, nu)   (nu + NDIM*(mu) + NDIM)
#define M_UP_UP_DOWN(mu, nu, lambda) (lambda + NDIM*(nu) + NDIM*NDIM*(mu) + NDIM*NDIM + NDIM)
#define M_UP_UP_UP(mu, nu, lambda)   (lambda + NDIM*(nu) + NDIM*NDIM*(mu) + NDIM*NDIM + NDIM)

/* Indices for elem.collisionIntegrals */
#define COLLISION_INTEGRAL(mu, nu) (nu + NDIM*(mu) + NUM_COLLISION_PARAMETERS)

/* Indices for the Christoffel symbols */
#define GAMMA_UP_DOWN_DOWN(eta,mu,nu) (eta+NDIM*(mu+NDIM*(nu) ) )

/* Macros for the REAPER scheme */
/* Quad integration constants and mnemonics*/
#define MINUS_INFINITY    (-1)
#define ZERO              (0)
#define PLUS_INFINITY     (1)

/* Integration is over pDown0, pDown1, pDown2 */
#define NDIM_INTEGRATION  (3)
#define NUM_QUAD          (21)

/* Macros used in integrands.c 
 * FULL functions are used when integrating from -inf to inf 
 * HALF functions are used when integrating from -inf to 0 or 0 to inf
 * We only use FULL functions */
#define JACOBIAN_FULL(t) (  (1. + ( (t)*(t) ))/\
                            (\
                              (1. - ( (t)*(t) ) )\
                             *(1. - ( (t)*(t) ) ) \
                            )\
                         )

#define JACOBIAN_PLUS_HALF(t) (1./(\
                                    (1. - (t) )*(1. - (t) )\
                                  )\
                              )

#define JACOBIAN_MINUS_HALF(t) (1./(\
                                     (1. + (t) )*(1. + (t) )\
                                   )\
                               )

#define P_DOWN_FROM_T_FULL(t) ( (t)/(1. - (t)*(t) ) )

#define P_DOWN_FROM_T_PLUS_HALF(t) ( (t)/(1. - (t) ) )

#define P_DOWN_FROM_T_MINUS_HALF(t) ( (t)/(1. + (t) ) )

struct fluidElement
{
  REAL gamma;
  REAL primVars[DOF];
  REAL uCon[NDIM];
  REAL bCon[NDIM];
  REAL moments[NUM_ALL_COMPONENTS]; 
  /* TODO: We have 4 independent components of N^\mu and 10
  independent components of T^\mu^\nu, so do something later to exploit the
  symmetry of T^\mu^\nu */

  #if (CONDUCTION)
    REAL kappa, tau;
  #endif

  #if (REAPER)
    /* Orthonormal tetrad in the coordinate basis 
     *
     * $\overheadarrow{e}_\hat{\mu}$ =   eDownHatUpNoHat[mu][nu]
     *                                 * $\overheadarrow{e}_{\nu}$
     */
    REAL eDownHatUpNoHat[NDIM][NDIM];
    
    /* Coordinate basis in the orthonormal tetrad
     *
     * $\overheadarrow{e}_{\mu}$ =   eDownNoHatUpHat[mu][nu]
     *                             * $\overheadarrow{e}_\hat{nu}$
     */
    REAL eDownNoHatUpHat[NDIM][NDIM]; 

    #if (REAPER_MOMENTS==15)
      REAL collisionIntegrals[NUM_ALL_COLLISION_INTEGRALS];
    #endif
  #endif
};

#include "../geometry/geometry.h"
#include "../reconstruct/reconstruct.h"
#include "../timestepper/timestepper.h"
#include "../problem/problem.h"
#include <gsl/gsl_sf_bessel.h>

/* Public functions: */
void setFluidElement(const REAL primVars[ARRAY_ARGS DOF],
                     const struct geometry geom[ARRAY_ARGS 1],
                     struct fluidElement elem[ARRAY_ARGS 1]);

void computeMoments(const struct geometry geom[ARRAY_ARGS 1],
                    struct fluidElement elem[ARRAY_ARGS 1]);

void computeFluxes(const struct fluidElement elem[ARRAY_ARGS 1],
                   const struct geometry geom[ARRAY_ARGS 1],
                   const int dir,
                   REAL fluxes[ARRAY_ARGS NUM_FLUXES]);

void computeSourceTerms(const struct fluidElement elem[ARRAY_ARGS 1],
                        const struct geometry geom[ARRAY_ARGS 1],
                        const REAL christoffel[ARRAY_ARGS 64],
                        REAL sourceTerms[ARRAY_ARGS NUM_FLUXES]);

#if (CONDUCTION)
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
);

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
);
#endif

#if (REAPER)
  /* Public functions needed by the REAPER scheme in order to set initial 
   * conditions */
  REAL getAlpha(REAL rho, REAL temperature);
  
  REAL getA0(REAL temperature);

  REAL getTemperature(const struct fluidElement elem[ARRAY_ARGS 1]);

  REAL getDensity(const struct fluidElement elem[ARRAY_ARGS 1]);
#endif

/* Internal functions */
void setGamma(const struct geometry geom[ARRAY_ARGS 1],
              struct fluidElement elem[ARRAY_ARGS 1]);

void setUCon(const struct geometry geom[ARRAY_ARGS 1],
             struct fluidElement elem[ARRAY_ARGS 1]);

void setBCon(const struct geometry geom[ARRAY_ARGS 1],
             struct fluidElement elem[ARRAY_ARGS 1]);

REAL getbSqr(const struct fluidElement elem[ARRAY_ARGS 1],
             const struct geometry geom[ARRAY_ARGS 1]);

#if (REAPER)
/* Internal functions used by the reaper scheme */
void fixedQuadIntegration(const struct fluidElement *elem,
                          const struct geometry *geom,
                          REAL scaleFactor,
                          REAL *moments,
                          REAL *collisionIntegrals);

void computefAndPUpHat
(
  REAL pDownHat[NDIM],
  const struct geometry* geom,
  const struct fluidElement* elem,
  REAL pUpHat[NDIM],
  REAL *f,
  REAL *collisionOperator
);

/* Internal functions used in tetrad construction, also needed for the repaer
 * scheme */
void setTetrad(const struct geometry *geom,
               struct fluidElement *elem);

void makeTetrad(const REAL eDown0Hat[NDIM],
                const REAL eDown1Hat[NDIM],
                const struct geometry *geom,
                REAL eDownMuHatUpNu[NDIM][NDIM],
                REAL eDownMuUpNuHat[NDIM][NDIM]);

void normalize(const struct geometry *geom, 
               REAL vecCon[NDIM]);

void orthogonalize(const REAL inputVecCon[NDIM],
                   const struct geometry *geom,
                   REAL trialVecCon[NDIM]);
#endif

#endif /* GRIM_PHYSICS_H_ */
