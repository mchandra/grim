#ifndef GRIM_RIEMANNSOLVER_H_
#define GRIM_RIEMANNSOLVER_H_

#include "../inputs.h"
#include "../physics/physics.h"
#include "../geometry/geometry.h"

REAL riemannSolver(const REAL fluxLeft[ARRAY_ARGS NUM_FLUXES],
                   const REAL fluxRight[ARRAY_ARGS NUM_FLUXES],
                   const REAL conservedVarsLeft[ARRAY_ARGS NUM_FLUXES],
                   const REAL conservedVarsRight[ARRAY_ARGS NUM_FLUXES],
                   const REAL primVarsLeft[ARRAY_ARGS NUM_FLUXES],
                   const REAL primVarsRight[ARRAY_ARGS NUM_FLUXES],
                   const struct geometry geom[ARRAY_ARGS NUM_FLUXES],
                   const int dir, REAL fluxes[ARRAY_ARGS NUM_FLUXES]);

void waveSpeeds(const struct fluidElement elem[ARRAY_ARGS 1],
                const struct geometry geom[ARRAY_ARGS 1],
                const int dir,
                REAL cMin[ARRAY_ARGS 1], REAL cMax[ARRAY_ARGS 1]);

#endif /* GRIM_TIMESTEPPER_H_ */

