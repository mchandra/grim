import mpi4py, petsc4py
from petsc4py import PETSc
import numpy as np
import pylab as pl
import pytest
import gridPy
import geometryPy
import boundaryPy
import timeStepperPy

petsc4py.init()
petscComm  = petsc4py.PETSc.COMM_WORLD
comm = petscComm.tompi4py()
rank = comm.Get_rank()
numProcs = comm.Get_size()
PETSc.Sys.Print("Using %d procs" % numProcs)

N1  = int(pytest.config.getoption('N1'))
N2  = int(pytest.config.getoption('N2'))
N3  = int(pytest.config.getoption('N3'))
dim = int(pytest.config.getoption('dim'))

# Geometry parameters
blackHoleSpin = float(pytest.config.getoption('blackHoleSpin'))
hSlope        = float(pytest.config.getoption('hSlope'))
numGhost = 3

Rin = 0.98*(1.+np.sqrt(1.-blackHoleSpin*blackHoleSpin));
Rout = 40.

#X1Start = np.log(Rin); X1End = np.log(Rout)
#X2Start = 1e-8; X2End = 1.-1e-8
#X3Start = 0.; X3End = 2.*np.pi
#boundaryLeft   = boundaryPy.OUTFLOW
#boundaryRight  = boundaryPy.OUTFLOW
#boundaryTop    = boundaryPy.OUTFLOW
#boundaryBottom = boundaryPy.OUTFLOW
#boundaryFront  = boundaryPy.PERIODIC
#boundaryBack   = boundaryPy.PERIODIC

X1Start = 0.; X1End = 1.
X2Start = 0.; X2End = 1.
X3Start = 0.; X3End = 1.

boundaryLeft   = boundaryPy.PERIODIC
boundaryRight  = boundaryPy.PERIODIC
boundaryTop    = boundaryPy.PERIODIC
boundaryBottom = boundaryPy.PERIODIC
boundaryFront  = boundaryPy.PERIODIC
boundaryBack   = boundaryPy.PERIODIC


time = 0.
dt   = 0.0005
numVars = 5
#metric = geometryPy.MODIFIED_KERR_SCHILD
metric = geometryPy.MINKOWSKI
ts = timeStepperPy.timeStepperPy(N1, N2, N3,
                                 dim, numVars, numGhost,
                                 time, dt,
                                 boundaryLeft, boundaryRight,
                                 boundaryTop,  boundaryBottom,
                                 boundaryFront, boundaryBack,
                                 metric, blackHoleSpin, hSlope,
                                 X1Start, X1End,
                                 X2Start, X2End,
                                 X3Start, X3End
                                )
#ts.computeDivB(ts.primOld)
#print "divB.shape = ", ts.divB.shape
#print "Div B = ", np.max(np.abs(ts.divB.getVars()[0, 0, numGhost:N2+numGhost,
#  numGhost:N1+numGhost]))
#pl.contourf(np.log10(np.abs(ts.divB.getVars()[0, 0,numGhost:N2+numGhost,
#                                              numGhost:N1+numGhost])), 100)
#pl.colorbar()
#pl.savefig("divB_initial_conditions.png")
#pl.clf()
#ts.timeStep()
#ts.computeDivB(ts.primOld)
#pl.contourf(np.log10(np.abs(ts.divB.getVars()[0, 0,numGhost:N2+numGhost,
#                                              numGhost:N1+numGhost])), 100)
#pl.colorbar()
#pl.savefig("divB_timestepped.png")

#for n in xrange(10):
#    ts.timeStep()
#ts.computeDivOfFluxes(ts.primHalfStep)
#pl.contourf(ts.emfX3.getVars()[0, 0,numGhost:N2+numGhost,
#                                               numGhost:N1+numGhost], 100)
#pl.colorbar()
#pl.savefig("EMF_initial_conditions.png")
#np.savetxt("emf.txt", ts.divFluxes.getVars()[5, 0, numGhost:N2+numGhost,
#                                             numGhost:N1+numGhost
#                                      ]
#          )


for n in xrange(4000):
  print "Time step = ", n
  ts.timeStep()
  if (n%10 == 0):
    ts.computeDivB(ts.B1LeftHalfStep, ts.B2BottomHalfStep, ts.B3BackHalfStep)
    print "Div B primHalfStep = ", \
    np.max(np.abs(ts.divB.getVars()[0, 0, numGhost:N2+numGhost,
                                              numGhost:N1+numGhost]
                     )
              )
    ts.computeDivB(ts.B1LeftOld, ts.B2BottomOld, ts.B3BackOld)
    print "Div B primOld = ", \
        np.max(np.abs(ts.divB.getVars()[0, 0, numGhost:N2+numGhost, numGhost:N1+numGhost]))
  
    pl.contourf(np.log10(np.abs(ts.divB.getVars()[0, 0,numGhost:N2+numGhost,
                                            numGhost:N1+numGhost])), 100)
    pl.colorbar()
    pl.title("Time = " + str(n*dt))
    pl.savefig("divB_" + str(n) + ".png")
    pl.clf()

    pl.contourf(ts.primOld.getVars()[0, 0, :, :], 100)
    pl.title("Time = " + str(n*dt))
    pl.savefig("rho_" + str(n) + ".png")
    pl.clf()
