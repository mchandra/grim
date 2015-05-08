#ifndef GRIM_TIMESTEPPER_H_
#define GRIM_TIMESTEPPER_H_

#include <petsc.h>
#include <petscviewerhdf5.h>
#include <unistd.h>
#include "../inputs.h"

#define EXPLICIT (0)
#define IMPLICIT (1)
#define IMEX     (2)

struct problemData;

struct timeStepper
{
  REAL t, dt, tDump;
  int timeStepCounter;
  int dumpCounter;

  SNES snes;
  DM dmdaWithGhostZones;
  DM dmdaWithoutGhostZones;
  DM dmdaFluxes;
  DM connectionDMDA;
  DM dmdaDt;

  Vec primPetscVec;
  Vec residualPetscVec;
  Vec primPetscVecOld;
  Vec primPetscVecHalfStep;
  Vec sourceTermsPetscVecOld;
  Vec divFluxPetscVecOld;
  Vec conservedVarsPetscVecOld;
  Vec connectionPetscVec;
  Vec dtPetscVec;

  #if (CONDUCTION)
    DM  gradTDM;
    DM  graduConDM;
    DM  graduConHigherOrderTermsDM;
    Vec gradTPetscVec;
    Vec graduConPetscVec;
    Vec graduConHigherOrderTerm1PetscVec;
    Vec graduConHigherOrderTerm2PetscVec;
  #endif

  #if (REAPER)
    DM momentsDM;
    Vec momentsVec;
  #endif

  struct problemData *problemSpecificData;

  int computeOldSourceTermsAndOldDivOfFluxes;
  int computeDivOfFluxAtTimeN;
  int computeDivOfFluxAtTimeNPlusHalf;
  int computeSourceTermsAtTimeN;
  int computeSourceTermsAtTimeNPlusHalf;

  int X1Start, X1Size;
  int X2Start, X2Size;
  int X3Start, X3Size;
};

/* The problem library requires timeStepper struct, so put this include after
 * the definition of the struct */
#include "../geometry/geometry.h"
#include "../gridzone/gridzone.h"
#include "../boundary/boundary.h"
#include "../reconstruct/reconstruct.h"
#include "../riemannsolver/riemannsolver.h"
#include "../physics/physics.h"
#include "../problem/problem.h"

/* User functions */

void timeStepperInit(struct timeStepper ts[ARRAY_ARGS 1]);

void timeStep(struct timeStepper ts[ARRAY_ARGS 1]);

void timeStepperDestroy(struct timeStepper ts[ARRAY_ARGS 1]);

/* Internal functions */

PetscErrorCode computeResidual(SNES snes, 
                               Vec primPetscVec,
                               Vec residalPetscVec,
                               void *ptr);

void computeFluxesOverTile(const REAL primTile[ARRAY_ARGS TILE_SIZE],
                           const int iTile, const int jTile,
                           const int X1Start, const int X2Start,
                           const int X1Size, const int X2Size,
                           REAL fluxX1Tile[ARRAY_ARGS TILE_SIZE],
                           REAL fluxX2Tile[ARRAY_ARGS TILE_SIZE],
                           ARRAY(dtGlobal));

void setChristoffelSymbols(struct timeStepper ts[ARRAY_ARGS 1]);

void diagnostics(struct timeStepper ts[ARRAY_ARGS 1]);

#if (CONDUCTION)
  void initConductionDataStructures(struct timeStepper ts[ARRAY_ARGS 1]);
  void destroyConductionDataStructures(struct timeStepper ts[ARRAY_ARGS 1]);
#endif

#if (REAPER)
  void initReaperDataStructures(struct timeStepper ts[ARRAY_ARGS 1]);
  void destroyReaperDataStructures(struct timeStepper ts[ARRAY_ARGS 1]);
#endif

#endif /* GRIM_TIMESTEPPER_H_ */
