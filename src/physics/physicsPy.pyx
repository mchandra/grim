"""
physicsPy
======

Python interface to the fluidElement and riemannSolver class.
"""

import numpy as np
from gridPy cimport gridPy, coordinatesGridPy
from geometryPy cimport geometryPy
from physicsHeaders cimport fluidElement

cdef class fluidElementPy(object):

  def __cinit__(self, gridPy prim = gridPy(),
                      gridPy magneticFields = gridPy(),
                      geometryPy geom = geometryPy()
               ):
    cdef int numReads  = 0
    cdef int numWrites = 0
    if (prim.usingExternalPtr):
      self.usingExternalPtr = 1
      self.elemPtr = NULL
      return
      
    self.usingExternalPtr = 0
    self.elemPtr = new fluidElement(prim.getGridPtr()[0],
                                    magneticFields.getGridPtr()[0],
                                    geom.getGeometryPtr()[0],
                                    numReads, numWrites
                                   )

  def computeFluidFluxes(self, geometryPy geom,
                         const int direction,
                         gridPy flux
                        ):
    cdef int numReads  = 0
    cdef int numWrites = 0
    self.elemPtr.computeFluidFluxes(geom.getGeometryPtr()[0],
                                    direction,
                                    flux.getGridPtr()[0],
                                    numReads, numWrites
                                   )
    return numReads, numWrites

  def __dealloc__(self):
    if (self.usingExternalPtr):
      return
    del self.elemPtr

  cdef setElemPtr(self, fluidElement *elemPtr):
    self.elemPtr = elemPtr

  cdef fluidElement* getElemPtr(self, fluidElement *elemPtr):
    return self.elemPtr

  @staticmethod
  cdef createFluidElementPyFromElemPtr(fluidElement *elemPtr):
    cdef fluidElementPy fluidElementPyObject = fluidElementPy()
    fluidElementPyObject.usingExternalPtr = 1
    fluidElementPyObject.setElemPtr(elemPtr)
    return fluidElementPyObject
