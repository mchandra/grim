#include "grid.hpp"

grid::grid(const int N1,
           const int N2,
           const int N3,
           const int dim, 
           const int numVars,
           const int numGhost,
           const int periodicBoundariesX1,
           const int periodicBoundariesX2,
           const int periodicBoundariesX3
          )
{
  this->numVars  = numVars;
  this->numGhost = numGhost;
  this->N1 = N1;
  this->N2 = N2;
  this->N3 = N3;
  this->dim = dim;
  this->periodicBoundariesX1 = periodicBoundariesX1;
  this->periodicBoundariesX2 = periodicBoundariesX2;
  this->periodicBoundariesX3 = periodicBoundariesX3;

  hasHostPtrBeenAllocated = 0;

  /* Implementations for MIRROR, OUTFLOW in boundary.cpp and DIRICHLET in
   * problem.cpp */
  DMBoundaryLeft  = DM_BOUNDARY_GHOSTED; DMBoundaryRight  = DM_BOUNDARY_GHOSTED;
  DMBoundaryTop   = DM_BOUNDARY_GHOSTED; DMBoundaryBottom = DM_BOUNDARY_GHOSTED;
  DMBoundaryFront = DM_BOUNDARY_GHOSTED; DMBoundaryBack   = DM_BOUNDARY_GHOSTED;

  if (periodicBoundariesX1)
  {
    DMBoundaryLeft  = DM_BOUNDARY_PERIODIC;
    DMBoundaryRight = DM_BOUNDARY_PERIODIC;
  }

  if (periodicBoundariesX2)
  {
    DMBoundaryTop    = DM_BOUNDARY_PERIODIC;
    DMBoundaryBottom = DM_BOUNDARY_PERIODIC;
  }

  if (periodicBoundariesX3)
  {
    DMBoundaryFront = DM_BOUNDARY_PERIODIC;
    DMBoundaryBack  = DM_BOUNDARY_PERIODIC;
  }

  switch (dim)
  {
    case 1:
      numGhostX1 = numGhost;
      numGhostX2 = 0;
      numGhostX3 = 0;

      domainX1 = new af::seq(numGhost, af::end - numGhost);
      domainX2 = new af::seq(span);
      domainX3 = new af::seq(span);

      DMDACreate1d(PETSC_COMM_WORLD, DMBoundaryLeft, 
                   N1, numVars, numGhostX1, NULL,
                   &dm
                  );

      break;
  
    case 2:
      numGhostX1 = numGhost;
      numGhostX2 = numGhost;
      numGhostX3 = 0;

      domainX1 = new af::seq(numGhost, af::end - numGhost);
      domainX2 = new af::seq(numGhost, af::end - numGhost);
      domainX3 = new af::seq(span);

      DMDACreate2d(PETSC_COMM_WORLD, 
                   DMBoundaryLeft, DMBoundaryBottom,
                   DMDA_STENCIL_BOX,
                   N1, N2,
                   PETSC_DECIDE, PETSC_DECIDE,
                   numVars, numGhostX1,
                   PETSC_NULL, PETSC_NULL,
                   &dm
                  );

      break;

    case 3:
      numGhostX1 = numGhost;
      numGhostX2 = numGhost;
      numGhostX3 = numGhost;

      domainX1 = new af::seq(numGhost, af::end - numGhost);
      domainX2 = new af::seq(numGhost, af::end - numGhost);
      domainX3 = new af::seq(numGhost, af::end - numGhost);

      DMDACreate3d(PETSC_COMM_WORLD, 
                   DMBoundaryLeft, DMBoundaryBottom, DMBoundaryBack,
                   DMDA_STENCIL_BOX,
                   N1, N2, N3,
                   PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                   numVars, numGhostX1,
                   PETSC_NULL, PETSC_NULL, PETSC_NULL,
                   &dm
                  );

      break;
  }

  DMDAGetCorners
  (dm, &iLocalStart, &jLocalStart, &kLocalStart,
       &N1Local,     &N2Local,     &N3Local
  );

  N1Total = N1Local + 2*numGhostX1;
  N2Total = N2Local + 2*numGhostX2;
  N3Total = N3Local + 2*numGhostX3;

  DMCreateGlobalVector(dm, &globalVec);
  DMCreateLocalVector(dm, &localVec);
  VecSet(globalVec, 0.);
  VecSet(localVec,  0.);

  iLocalEnd = iLocalStart + N1Local;
  jLocalEnd = jLocalStart + N2Local;
  kLocalEnd = kLocalStart + N3Local;

  vars = new array[numVars];

  /* Initialize vars to localVec */
  copyLocalVecToVars();
}

coordinatesGrid::coordinatesGrid
  (const int N1, 
   const int N2,
   const int N3,
   const int dim,
   const int numGhost,
   const double X1Start, const double X1End,
   const double X2Start, const double X2End,
   const double X3Start, const double X3End
  ) : grid(N1, N2, N3, dim, 3, numGhost, false, false, false)
{
  this->X1Start = X1Start; this->X1End = X1End;
  this->X2Start = X2Start; this->X2End = X2End;
  this->X3Start = X3Start; this->X3End = X3End;

  dX1 = (X1End - X1Start)/N1;
  dX2 = (X2End - X2Start)/N2;
  dX3 = (X3End - X3Start)/N3;
}

void coordinatesGrid::setXCoords(const int location)
{
  this->location = location;

  array indicesX1
    = af::range(N1Total, /* number of total zones in X1 */
                N2Total, /* number of total zones in X2 */
                N3Total, /* number of total zones in X3 */
                1,                /* number of variables */
                directions::X1 ,  /* Vary indices in X1 direction */
                f64               /* Double precision */
               ) - numGhostX1 + iLocalStart; /* Offsets for MPI */

  array indicesX2
    = af::range(N1Total, /* number of total zones in X1 */
                N2Total, /* number of total zones in X2 */
                N3Total, /* number of total zones in X3 */
                1,                /* number of variables */
                directions::X2 ,  /* Vary indices in X1 direction */
                f64               /* Double precision */
               ) - numGhostX2 + jLocalStart; /* Offsets for MPI */

  array indicesX3
    = af::range(N1Total, /* number of total zones in X1 */
                N2Total, /* number of total zones in X2 */
                N3Total, /* number of total zones in X3 */
                1,                /* number of variables */
                directions::X3 ,  /* Vary indices in X1 direction */
                f64               /* Double precision */
               ) - numGhostX3 + kLocalStart; /* Offsets for MPI */


  switch (location)
  {
    case locations::CENTER:
      vars[directions::X1] = X1Start + (indicesX1 + 0.5)*dX1;
      vars[directions::X2] = X2Start + (indicesX2 + 0.5)*dX2;
      vars[directions::X3] = X3Start + (indicesX3 + 0.5)*dX3;

      break;

    case locations::LEFT:
      vars[directions::X1] = X1Start + (indicesX1      )*dX1;
      vars[directions::X2] = X2Start + (indicesX2 + 0.5)*dX2;
      vars[directions::X3] = X3Start + (indicesX3 + 0.5)*dX3;

      break;

    case locations::RIGHT:
      vars[directions::X1] = X1Start + (indicesX1 +   1)*dX1;
      vars[directions::X2] = X2Start + (indicesX2 + 0.5)*dX2;
      vars[directions::X3] = X3Start + (indicesX3 + 0.5)*dX3;

      break;

    case locations::BOTTOM:
      vars[directions::X1] = X1Start + (indicesX1 + 0.5)*dX1;
      vars[directions::X2] = X2Start + (indicesX2      )*dX2;
      vars[directions::X3] = X3Start + (indicesX3 + 0.5)*dX3;

      break;

    case locations::TOP:
      vars[directions::X1] = X1Start + (indicesX1 + 0.5)*dX1;
      vars[directions::X2] = X2Start + (indicesX2 +   1)*dX2;
      vars[directions::X3] = X3Start + (indicesX3 + 0.5)*dX3;

      break;

    case locations::FRONT:
      vars[directions::X1] = X1Start + (indicesX1 + 0.5)*dX1;
      vars[directions::X2] = X2Start + (indicesX2 + 0.5)*dX2;
      vars[directions::X3] = X3Start + (indicesX3 +   1)*dX3;

      break;

    case locations::BACK:
      vars[directions::X1] = X1Start + (indicesX1 + 0.5)*dX1;
      vars[directions::X2] = X2Start + (indicesX2 + 0.5)*dX2;
      vars[directions::X3] = X3Start + (indicesX3      )*dX3;

      break;
  }
}

void grid::copyLocalVecToVars()
{
  array varsCopiedFromLocalVec;

  double *pointerToLocalVec;
  VecGetArray(localVec, &pointerToLocalVec);

  varsCopiedFromLocalVec = 
    array(numVars,
          N1Total, N2Total, N3Total, 
          pointerToLocalVec
         );

  VecRestoreArray(localVec, &pointerToLocalVec);

  varsSoA = af::reorder(varsCopiedFromLocalVec, 1, 2, 3, 0);
  for (int var=0; var<numVars; var++)
  {
    vars[var] = varsSoA(span, span, span, var);
  }
}

void grid::copyVarsToHostPtr()
{
  for (int var=0; var < numVars; var++)
  {
            /*i,    j,    k*/
    varsSoA(span, span, span, var) = vars[var];
  }

  /* TODO: Find out if I need to delete and then copy again or can I just copy
   * directly */
  /* If already called once, free the existing memory */
  if (hasHostPtrBeenAllocated)
  {
    delete hostPtr;
  }

  /* Allocate space for hostPtr and then copy data from device to host */
  hostPtr =  varsSoA.host<double>();
  hasHostPtrBeenAllocated = 1;
}

void grid::copyHostPtrToVars(const double *hostPtr)
{
  varsSoA = array(N1Total, N2Total, N3Total, numVars,
                  hostPtr
                 );

  for (int var=0; var < numVars; var++)
  {
    vars[var] = varsSoA(span, span, span, var);
  }
}

void grid::communicate()
{
  copyVarsToGlobalVec();

  DMGlobalToLocalBegin(dm, globalVec, INSERT_VALUES, localVec);
  DMGlobalToLocalEnd(dm, globalVec, INSERT_VALUES, localVec);

  /* Now copy data from localVec to vars */
  copyLocalVecToVars();
}

void grid::copyVarsToGlobalVec()
{
  /* Get data into Array of Structs format as needed by Petsc */
  for (int var=0; var < numVars; var++)
  {
    varsSoA(span, span, span, var) = vars[var];
  }

  array varsAoS = af::reorder(varsSoA, 3, 0, 1, 2);

  /* Copy part of varsAoS contained in [0, N1) x [0, N2) x [0, N3) to global vec */
  double *pointerToLocalVec;
  double *pointerToGlobalVec;
  VecGetArray(localVec, &pointerToLocalVec);
  VecGetArray(globalVec, &pointerToGlobalVec);

  varsAoS.host(pointerToLocalVec);

  for (int k=0; k<N3Local; k++)
  {
    for (int j=0; j<N2Local; j++)
    {
      for (int i=0; i<N1Local; i++)
      {
        for (int var=0; var<numVars; var++)
        {
	        /* Note on indices: Petsc uses non-ghosted vectors for the global vector,
	        * but ghosted indices for the local vector! */
          const int globalindex = var + numVars*(i + N1Local*(j + N2Local*(k) ) );
          const int localindex  = var + numVars*(i + numGhostX1 
                                                   + N1Total*(j + numGhostX2
                                                              + N2Total*(k + numGhostX3)
                                                             ) 
                                                );

          pointerToGlobalVec[globalindex] = pointerToLocalVec[localindex];
        }
      }
    }
  }

  VecRestoreArray(localVec, &pointerToLocalVec);
  VecRestoreArray(globalVec, &pointerToGlobalVec);
}

coordinatesGrid::~coordinatesGrid()
{

}

void grid::dump(const std::string varsName, const std::string fileName)
{
  copyVarsToGlobalVec();

  PetscViewer viewer;
  PetscViewerHDF5Open(PETSC_COMM_WORLD,
                      fileName.c_str(), FILE_MODE_WRITE, &viewer
                     );

  /* Output the variables */
  PetscObjectSetName((PetscObject) globalVec, varsName.c_str());
  VecView(globalVec, viewer);

  PetscViewerDestroy(&viewer);
}

grid::~grid()
{
  if (hasHostPtrBeenAllocated)
  {
    delete hostPtr;
  }
  delete [] vars;
  VecDestroy(&globalVec);
  VecDestroy(&localVec);

  DMDestroy(&dm);

  delete domainX1, domainX2, domainX3;
}
