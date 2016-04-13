from gridHeaders cimport grid, coordinatesGrid
from geometryHeaders cimport geometry

cdef extern from "physics.hpp":
  cdef cppclass fluidElement:
    fluidElement(const grid &prim,
                 const grid &magneticFields,
                 const geometry &geom,
                 int &numReads,
                 int &numWrites
                )
    void computeFluidFluxes(const geometry &geom, 
                            const int direction,
                            grid &flux,
                            int &numReads,
                            int &numWrites
                           )
