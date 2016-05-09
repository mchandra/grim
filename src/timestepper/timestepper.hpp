#ifndef GRIM_TIMESTEPPER_H_
#define GRIM_TIMESTEPPER_H_

#include "../params.hpp"
#include "../grid/grid.hpp"
#include "../physics/physics.hpp"
#include "../geometry/geometry.hpp"
#include "../boundary/boundary.hpp"
#include "mkl.h"

namespace timeStepperSwitches
{
  enum
  {
    HALF_STEP, FULL_STEP
  };
};

class timeStepper
{
  grid *primGuessLineSearchTrial;
  grid *primGuessPlusEps;
  grid *residual;
  grid *residualPlusEps;

  array residualSoA;
  array jacobianSoA;
  array deltaPrimAoS;
  array stepLength;

  double *AHostPtr, *bHostPtr, *xHostPtr;

  void solve(grid &primGuess);
  void computeResidual(const grid &prim, grid &residual,
		                   const bool ComputeExplicitTerms,
                       int &numReads,
                       int &numWrites
                      );
  void batchLinearSolve(const array &A, const array &b, array &x);

  geometry *geomEdge;
  public:
    double dt, time;
    int N1, N2, N3, numGhost;
    int numVars;

    int timeStepCounter;
    int dumpCounter;

    int boundaryLeft, boundaryRight;
    int boundaryTop,  boundaryBottom;
    int boundaryFront, boundaryBack;

    coordinatesGrid *XCoords;
    grid *prim, *primHalfStep, *primOld;
    grid *cons, *consOld;
    grid *sourcesExplicit;
    grid *sourcesImplicit;
    grid *sourcesImplicitOld;
    grid *sourcesTimeDer;
    grid *primLeft, *primRight;
    grid *fluidFluxesX1,    *fluidFluxesX2,    *fluidFluxesX3;
    grid *magneticFluxesX1, *magneticFluxesX2, *magneticFluxesX3;

    grid *divFluxes;
    grid *divB;
    /* Primary magnetic field variables, indicated by their respective
     * locations */
    grid *B1Left,         *B2Bottom,          *B3Back;
    grid *B1LeftOld,      *B2BottomOld,       *B3BackOld;
    grid *B1LeftHalfStep, *B2BottomHalfStep,  *B3BackHalfStep;
    
    /* Cell centered magnetic fields, computed using interpolation from 
     * face-centered fields */
    grid *B1Center, *B2Center, *B3Center;

    /* Reconstructed face-centered magnetic fields, computed using
     * reconstruction of the cell-centered magnetic fields */
    grid *B1Bottom, *B1Top,   *B1Back,   *B1Front;
    grid *B2Left,   *B2Right, *B2Back,   *B2Front;
    grid *B3Left,   *B3Right, *B3Bottom, *B3Top;

    /* Primary electric field variables, indicated by their locations */
    grid *E1BottomBack, *E2LeftBack, *E3LeftBottom;

    /* Auxiliary grids needed for inputs into the Riemann solver */
    grid *magneticFieldsLeft, *magneticFieldsRight;
    grid *magneticFieldsCenter;

    geometry *geomLeft,   *geomRight;
    geometry *geomBottom, *geomTop;
    geometry *geomCenter;
    array gLeftBottomEdge,  gLeftTopEdge;
    array gRightBottomEdge, gRightTopEdge;
    array gBottomBackEdge,  gBottomFrontEdge;
    array gTopBackEdge,     gTopFrontEdge;
    array gLeftBackEdge,    gLeftFrontEdge;
    array gRightBackEdge,   gRightFrontEdge;

    fluidElement *elem, *elemOld, *elemHalfStep;

    riemannSolver *riemann;

    void computeCellCenteredMagneticFields(const grid &B1Left,
                                           const grid &B2Bottom,
                                           const grid &B3Back,          
                                           int &numReads,
                                           int &numWrites
                                          );
    void computeEdgeElectricFields(int &numReads,
                                   int &numWrites
                                  );

    void computeDivOfFluxes(const grid &primFlux,
                            const grid &B1Left,
                            const grid &B2Bottom,
                            const grid &B3Back,
                            int &numReads,
                            int &numWrites
                           );

    int currentStep;

    timeStepper(const int N1, 
                const int N2,
                const int N3,
                const int dim,
                const int numVars, 
                const int numGhost,
                const double time,
                const double dt,
                const int boundaryLeft,  const int boundaryRight,
                const int boundaryTop,   const int boundaryBottom,
                const int boundaryFront, const int boundaryBack,
                const int metric,
                const double blackHoleSpin,
                const double hSlope,
                const double X1Start, const double X1End,
                const double X2Start, const double X2End,
                const double X3Start, const double X3End
               );
    ~timeStepper();

    void timeStep(int &numReads, int &numWrites);
    void timeStepBFields(const double dt,
                         int &numReads,
                         int &numWrites
                        );

    void fluxCT(int &numReads, int &numWrites);
    void computeEMF(int &numReadsEMF, int &numWritesEMF);
    void computeDivB(const grid &B1Left,
                     const grid &B2Bottom,
                     const grid &B3Back,
                     int &numReads,
                     int &numWrites
                    );

    /* Function definitions in the problem folder */
    void initialConditions(int &numReads, int &numWrites);
    void halfStepDiagnostics(int &numReads, int &numWrites);
    void fullStepDiagnostics(int &numReads, int &numWrites);
    void setProblemSpecificBCs(grid &prim,
                               grid &B1Left,
                               grid &B2Bottom,
                               grid &B3Back,
                               int &numReads, 
                               int &numWrites
                              );
};

#endif /* GRIM_TIMESTEPPER_H_ */
