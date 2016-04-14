import numpy as np
cimport numpy as np
from timestepperHeaders cimport timeStepper
from gridPy cimport gridPy, coordinatesGridPy
from geometryPy cimport geometryPy
from physicsPy cimport fluidElementPy

cdef class timeStepperPy(object):
  cdef timeStepper *timeStepperPtr
  cdef coordinatesGridPy XCoords
  cdef geometryPy geomCenter
  cdef geometryPy geomLeft, geomRight
  cdef geometryPy geomTop, geomBottom
  cdef geometryPy geomFront, geomBack
  cdef gridPy prim, primOld, primHalfStep
  cdef gridPy B1Left, B2Bottom, B3Back
  cdef gridPy B1LeftHalfStep, B2BottomHalfStep, B3BackHalfStep
  cdef gridPy B1LeftOld, B2BottomOld, B3BackOld
  cdef gridPy fluidFluxesX1,    fluidFluxesX2, fluidFluxesX3
  cdef gridPy magneticFluxesX1, magneticFluxesX2, magneticFluxesX3
  cdef gridPy divFluxes
  cdef gridPy E1BottomBack, E2LeftBack, E3LeftBottom
  cdef gridPy sourcesExplicit
  cdef gridPy divB
  cdef fluidElementPy elem, elemOld, elemHalfStep
