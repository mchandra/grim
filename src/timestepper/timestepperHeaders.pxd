from gridHeaders cimport grid, coordinatesGrid
from geometryHeaders cimport geometry
from physicsHeaders cimport fluidElement

cdef extern from "timestepper.hpp":
  cdef enum:
    STAGE_HALF_STEP "timeStepperSwitches::HALF_STEP"
    STAGE_FULL_STEP "timeStepperSwitches::FULL_STEP"

  cdef cppclass timeStepper:
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
               )
    coordinatesGrid *XCoords
    grid *prim
    grid *primOld
    grid *primHalfStep
    grid *fluidFluxesX1
    grid *fluidFluxesX2
    grid *fluidFluxesX3
    grid *magneticFluxesX1
    grid *magneticFluxesX2
    grid *magneticFluxesX3
    grid *divFluxes
    grid *divB
    grid *B1Left
    grid *B2Bottom
    grid *B3Back
    grid *B1LeftHalfStep
    grid *B2BottomHalfStep
    grid *B3BackHalfStep
    grid *B1LeftOld
    grid *B2BottomOld
    grid *B3BackOld
    grid *B1Center
    grid *B2Center
    grid *B3Center
    grid *E1BottomBack
    grid *E2LeftBack
    grid *E3LeftBottom
    grid *cons
    grid *consOld
    grid *sourcesExplicit
    grid *sourcesImplicit
    grid *sourcesImplicitOld
    grid *sourcesTimeDer

    fluidElement *elem
    fluidElement *elemOld
    fluidElement *elemHalfStep

    geometry *geomCenter
    geometry *geomLeft
    geometry *geomRight
    geometry *geomBottom
    geometry *geomTop

    void timeStep(int &numReads, int &numWrites)

    void computeEdgeElectricFields(int &numReads, int &numWrites)
    void computeDivOfFluxes(const grid &primFlux,
                            const grid &B1Left,
                            const grid &B2Bottom,
                            const grid &B3Back,
                            int &numReads,
                            int &numWrites
                           );
    void computeDivB(const grid &B1Left,
                     const grid &B2Bottom,
                     const grid &B3Back,
                     int &numReads,
                     int &numWrites
                    )

