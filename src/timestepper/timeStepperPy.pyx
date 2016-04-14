"""
timeStepperPy
======

Python interface to the timeStepper class.
"""

import numpy as np
from gridPy cimport gridPy, coordinatesGridPy
from geometryPy cimport geometryPy
from physicsPy cimport fluidElementPy
from timestepperHeaders cimport timeStepper
from timestepperHeaders cimport STAGE_HALF_STEP
from timestepperHeaders cimport STAGE_FULL_STEP

# Time stepping stage macros
HALF_STEP = STAGE_HALF_STEP
FULL_STEP = STAGE_FULL_STEP

cdef class timeStepperPy(object):

  def __cinit__(self, const int N1,
                      const int N2,
                      const int N3,
                      const int dim,
                      const int numVars,
                      const int numGhost,
                      const double time,
                      const double dt,
                      const int boundaryLeft, const int boundaryRight,
                      const int boundaryTop,  const int boundaryBottom,
                      const int boundaryFront, const int boundaryBack,
                      const int metric,
                      const double blackHoleSpin,
                      const double hSlope,
                      const double X1Start, const double X1End,
                      const double X2Start, const double X2End,
                      const double X3Start, const double X3End
               ):
    self.timeStepperPtr = \
        new timeStepper(N1, N2, N3, dim, numVars, numGhost,
                        time, dt,
                        boundaryLeft, boundaryRight,
                        boundaryTop,  boundaryBottom,
                        boundaryFront, boundaryBack,
                        metric, blackHoleSpin, hSlope,
                        X1Start, X1End,
                        X2Start, X2End,
                        X3Start, X3End
                       )

    self.XCoords = coordinatesGridPy()
    self.XCoords.setGridPtr(self.timeStepperPtr.XCoords)

    self.prim         = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.prim)
    self.primOld      = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.primOld)
    self.primHalfStep = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.primHalfStep)

    self.fluidFluxesX1 = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.fluidFluxesX1)
    self.fluidFluxesX2 = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.fluidFluxesX2)
    self.fluidFluxesX3 = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.fluidFluxesX3)

    self.magneticFluxesX1 = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.magneticFluxesX1)
    self.magneticFluxesX2 = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.magneticFluxesX2)
    self.magneticFluxesX3 = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.magneticFluxesX3)

    self.divFluxes = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.divFluxes)

    self.B1Left   = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.B1Left)
    self.B2Bottom = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.B2Bottom)
    self.B3Back   = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.B3Back)

    self.B1LeftHalfStep   = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.B1LeftHalfStep)
    self.B2BottomHalfStep = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.B2BottomHalfStep)
    self.B3BackHalfStep   = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.B3BackHalfStep)

    self.B1LeftOld   = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.B1LeftOld)
    self.B2BottomOld = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.B2BottomOld)
    self.B3BackOld   = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.B3BackOld)

    self.E1BottomBack = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.E1BottomBack)
    self.E2LeftBack = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.E2LeftBack)
    self.E3LeftBottom = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.E3LeftBottom)

    self.sourcesExplicit = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.sourcesExplicit)

    self.divB = \
        gridPy.createGridPyFromGridPtr(self.timeStepperPtr.divB)

    self.geomCenter = \
        geometryPy.createGeometryPyFromGeometryPtr(self.timeStepperPtr.geomCenter)
    self.geomLeft = \
        geometryPy.createGeometryPyFromGeometryPtr(self.timeStepperPtr.geomLeft)
    self.geomRight = \
        geometryPy.createGeometryPyFromGeometryPtr(self.timeStepperPtr.geomRight)
    self.geomBottom = \
        geometryPy.createGeometryPyFromGeometryPtr(self.timeStepperPtr.geomBottom)
    self.geomTop = \
        geometryPy.createGeometryPyFromGeometryPtr(self.timeStepperPtr.geomTop)

    self.elem = \
        fluidElementPy.createFluidElementPyFromElemPtr(self.timeStepperPtr.elem)


  def __dealloc__(self):
    del self.timeStepperPtr

  property XCoords:
    def __get__(self):
     return self.XCoords

  property elem:
    def __get__(self):
     return self.elem

  property geomCenter:
    def __get__(self):
     return self.geomCenter

  property geomRight:
    def __get__(self):
     return self.geomRight

  property geomLeft:
    def __get__(self):
     return self.geomLeft

  property geomBottom:
    def __get__(self):
     return self.geomBottom

  property geomTop:
    def __get__(self):
     return self.geomTop

  property prim:
    def __get__(self):
     return self.prim

  property primOld:
    def __get__(self):
     return self.primOld

  property primHalfStep:
    def __get__(self):
     return self.primHalfStep

  property fluidFluxesX1:
    def __get__(self):
     return self.fluidFluxesX1

  property fluidFluxesX2:
    def __get__(self):
     return self.fluidFluxesX2

  property fluidFluxesX3:
    def __get__(self):
     return self.fluidFluxesX3

  property divFluxes:
    def __get__(self):
     return self.divFluxes

  property B1Left:
    def __get__(self):
     return self.B1Left

  property B2Bottom:
    def __get__(self):
     return self.B2Bottom

  property B3Back:
    def __get__(self):
     return self.B3Back

  property B1LeftHalfStep:
    def __get__(self):
     return self.B1LeftHalfStep

  property B2BottomHalfStep:
    def __get__(self):
     return self.B2BottomHalfStep

  property B3BackHalfStep:
    def __get__(self):
     return self.B3BackHalfStep

  property B1LeftOld:
    def __get__(self):
     return self.B1LeftOld

  property B2BottomOld:
    def __get__(self):
     return self.B2BottomOld

  property B3BackOld:
    def __get__(self):
     return self.B3BackOld

  property E1BottomBack:
    def __get__(self):
     return self.E1BottomBack

  property E2LeftBack:
    def __get__(self):
     return self.E2LeftBack

  property E3LeftBottom:
    def __get__(self):
     return self.E3LeftBottom

  property sourcesExplicit:
    def __get__(self):
     return self.sourcesExplicit

  property divB:
    def __get__(self):
     return self.divB

  def computeDivB(self, gridPy B1Left, gridPy B2Bottom, gridPy B3Back):
     cdef int numReads  = 0
     cdef int numWrites = 0
     self.timeStepperPtr.computeDivB(B1Left.getGridPtr()[0],
                                     B2Bottom.getGridPtr()[0],
                                     B3Back.getGridPtr()[0],
                                     numReads, numWrites
                                    )
     return numReads, numWrites

  def computeEdgeElectricFields(self):
     cdef int numReads  = 0
     cdef int numWrites = 0
     self.timeStepperPtr.computeEdgeElectricFields(numReads, numWrites)
     return numReads, numWrites

  def timeStep(self):
      cdef int numReads  = 0
      cdef int numWrites = 0
      self.timeStepperPtr.timeStep(numReads, numWrites)
      return numReads, numWrites

  def computeDivOfFluxes(self, gridPy prim, 
                               gridPy B1Left,
                               gridPy B2Bottom,
                               gridPy B3Back
                        ):
    cdef int numReads  = 0
    cdef int numWrites = 0
    self.timeStepperPtr.computeDivOfFluxes(prim.getGridPtr()[0],
                                           B1Left.getGridPtr()[0],
                                           B2Bottom.getGridPtr()[0],
                                           B3Back.getGridPtr()[0],
                                           numReads, numWrites
                                          )
    return numReads, numWrites
