#ifndef GRIM_RECONSTRUCT_H_
#define GRIM_RECONSTRUCT_H_

#include "../inputs.h"
#include "../physics/physics.h"

#define X1 (1)
#define X2 (2)

/* RECONSTRUCTION options */
#define MONOTONIZED_CENTRAL (0)
#define MIN_MOD             (1)
#define MP5                 (2)

/* The following macros taken from PLUTO. Needed for MP5 */
#define MINMOD(a,b)  ((a)*(b) > 0.0 ? (fabs(a) < fabs(b) ? (a):(b)):0.0)
#define MY_MAX(a,b)  ( (a) >= (b) ? (a) : (b) ) 
#define MY_MIN(a,b)  ( (a) <= (b) ? (a) : (b) ) 
#define MEDIAN(a,b,c) (a + MINMOD(b-a,c-a))


/* Only Monotonized central slope limiter implemented so far. */
void slopeLim(const REAL left[ARRAY_ARGS NUM_FLUXES],
              const REAL mid[ARRAY_ARGS NUM_FLUXES],
              const REAL right[ARRAY_ARGS NUM_FLUXES],
              REAL ans[ARRAY_ARGS NUM_FLUXES]);

REAL slopeLimitedDerivative(const REAL left, const REAL mid, const REAL right);

void reconstruct(const REAL primTile[ARRAY_ARGS TILE_SIZE],
                 const int dir,
                 const int iTile, const int jTile,
                 const int X1Start, const int X2Start,
                 const int X1Size, const int X2Size,
                 REAL primVarsLeft[ARRAY_ARGS TILE_SIZE],
                 REAL primVarsRight[ARRAY_ARGS TILE_SIZE]);

/* MP5 reconstruction from Suresh and Huynh, 1997. Code taken from PLUTO */
REAL MP5_Reconstruct(const REAL Fjm2,
                     const REAL Fjm1,
                     const REAL Fj,
                     const REAL Fjp1,
                     const REAL Fjp2);

#endif /* GRIM_RECONSTRUCT_H_ */
