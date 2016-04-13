#include "timestepper.hpp"

void timeStepper::computeDivOfFluxes(const grid &primFlux,
                                     int &numReads,
                                     int &numWrites
                                    )
{
  int numReadsReconstruction, numWritesReconstruction;
  int numReadsRiemann, numWritesRiemann;
  int numReadsCT, numWritesCT;

  /* 1) Compute cell-centered B-fields */
  computeCellCenteredMagneticFields(numReads, numWrites);

  /*---- directions::X1 ----*/
  /* 2) Reconstuct fluid variables from cell-averages to i-1/2+eps and
   *    i+1/2-eps */
  reconstruction::reconstruct(primFlux, directions::X1,
                              *primLeft, *primRight,
                              numReadsReconstruction,
                              numWritesReconstruction
                             );
  numReads  += numReadsReconstruction;
  numWrites += numWritesReconstruction;

  /* 3) Reconstruction B2 from cell centers to i-1/2+eps and i+1/2-eps */
  reconstruction::reconstruct(*B2Center, directions::X1,
                              *B2Left, *B2Right,
                              numReadsReconstruction,
                              numWritesReconstruction
                             );
  numReads  += numReadsReconstruction;
  numWrites += numWritesReconstruction;

  /* 4) Reconstruction B3 from cell centers to i-1/2+eps and i+1/2-eps */
  reconstruction::reconstruct(*B3Center, directions::X1,
                              *B3Left, *B3Right,
                              numReadsReconstruction,
                              numWritesReconstruction
                             );
  numReads  += numReadsReconstruction;
  numWrites += numWritesReconstruction;

  /* 5) Set magneticFieldsLeft and magneticFieldsRight with the appropriate
   * magnetic fields. Note that we don't need to reconstruct B1 along X1
   * since it lives on the left faces : just shift to get the right face
   * values */
  magneticFieldsLeft->vars[vars::B1] = B1Left->vars[0];
  magneticFieldsLeft->vars[vars::B2] = B2Left->vars[0];
  magneticFieldsLeft->vars[vars::B3] = B3Left->vars[0];

  magneticFieldsRight->vars[vars::B1] = 
    shift(B1Left->vars[0], -1, 0, 0);
  magneticFieldsRight->vars[vars::B2] = B2Right->vars[0];
  magneticFieldsRight->vars[vars::B3] = B3Right->vars[0];

  /* 6) Finally, solve the Riemann problem. Inputs are values of prim, geom
   * and the magnetic fields at i-1/2+eps and i+1/2-eps */
  riemann->solve(*primLeft,           *primRight,
                 *magneticFieldsLeft, *magneticFieldsRight,
                 *geomLeft,           *geomRight,
                 directions::X1, 
                 *fluidFluxesX1, *magneticFluxesX1,
                 numReadsRiemann, numWritesRiemann
                );
  numReads  += numReadsRiemann;
  numWrites += numWritesRiemann;

  if (primFlux.dim >= 2)
  {
    /*---- directions::X2 ----*/
    /* 7) Reconstuct fluid variables from cell-averages to j-1/2+eps and
    *    j+1/2-eps */
    reconstruction::reconstruct(primFlux, directions::X2,
                                *primLeft, *primRight,
                                numReadsReconstruction,
                                numWritesReconstruction
                               );
    numReads  += numReadsReconstruction;
    numWrites += numWritesReconstruction;

    /* 8) Reconstruction B1 from cell centers to j-1/2+eps and j+1/2-eps */
    reconstruction::reconstruct(*B1Center, directions::X2,
                                *B1Bottom, *B1Top,
                                numReadsReconstruction,
                                numWritesReconstruction
                               );
    numReads  += numReadsReconstruction;
    numWrites += numWritesReconstruction;

    /* 9) Reconstruction B3 from cell centers to j-1/2+eps and j+1/2-eps */
    reconstruction::reconstruct(*B3Center, directions::X2,
                                *B3Bottom, *B3Top,
                                numReadsReconstruction,
                                numWritesReconstruction
                               );
    numReads  += numReadsReconstruction;
    numWrites += numWritesReconstruction;

    /* 10) Set magneticFieldsLeft and magneticFieldsRight with the appropriate
     * magnetic fields. Note that we don't need to reconstruct B2 along X2
     * since it lives on the bottom faces : just shift to get the top face
     * values */
    magneticFieldsLeft->vars[vars::B1] = B1Bottom->vars[0];
    magneticFieldsLeft->vars[vars::B2] = B2Bottom->vars[0];
    magneticFieldsLeft->vars[vars::B3] = B3Bottom->vars[0];

    magneticFieldsRight->vars[vars::B1] = B1Top->vars[0];
    magneticFieldsRight->vars[vars::B2] =
      af::shift(B2Bottom->vars[0], 0, -1, 0);
    magneticFieldsRight->vars[vars::B3] = B3Top->vars[0];

    /* 11) Finally, solve the Riemann problem. Inputs are values of prim, geom
     * and the magnetic fields at j-1/2+eps and j+1/2-eps */
    riemann->solve(*primLeft,           *primRight,
                   *magneticFieldsLeft, *magneticFieldsRight,
                   *geomBottom,         *geomTop,
                   directions::X2,
                   *fluidFluxesX2, *magneticFluxesX2,
                   numReadsRiemann, numWritesRiemann
                  );
    numReads  += numReadsRiemann;
    numWrites += numWritesRiemann;

    if (primFlux.dim == 3)
    {
      /*---- directions::X3 ----*/
      /* 12) Reconstuct fluid variables from cell-averages to k-1/2+eps and
       *     k+1/2-eps */
      reconstruction::reconstruct(primFlux, directions::X3,
                                  *primLeft, *primRight,
                                  numReadsReconstruction,
                                  numWritesReconstruction
                                 );
      numReads  += numReadsReconstruction;
      numWrites += numWritesReconstruction;

      /* 8) Reconstruction B1 from cell centers to k-1/2+eps and k+1/2-eps */
      reconstruction::reconstruct(*B1Center, directions::X2,
                                  *B1Back, *B1Front,
                                  numReadsReconstruction,
                                  numWritesReconstruction
                                 );
      numReads  += numReadsReconstruction;
      numWrites += numWritesReconstruction;

      /* 9) Reconstruction B2 from cell centers to k-1/2+eps and k+1/2-eps */
      reconstruction::reconstruct(*B2Center, directions::X2,
                                  *B2Back, *B2Front,
                                  numReadsReconstruction,
                                  numWritesReconstruction
                                 );
      numReads  += numReadsReconstruction;
      numWrites += numWritesReconstruction;

      /* 10) Set magneticFieldsLeft and magneticFieldsRight with the appropriate
       * magnetic fields. Note that we don't need to reconstruct B3 along X3
       * since it lives on the back faces : just shift to get the front face
       * values */
      magneticFieldsLeft->vars[vars::B1] = B1Back->vars[0];
      magneticFieldsLeft->vars[vars::B2] = B2Back->vars[0];
      magneticFieldsLeft->vars[vars::B3] = B3Back->vars[0];

      magneticFieldsRight->vars[vars::B1] = B1Front->vars[0];
      magneticFieldsRight->vars[vars::B2] = B2Front->vars[0];
      magneticFieldsRight->vars[vars::B3] =
        shift(B3Back->vars[0], 0, 0, -1);

      /* 11) Finally, solve the Riemann problem. Inputs are values of prim, geom
       * and the magnetic fields at k-1/2+eps and k+1/2-eps */
      riemann->solve(*primLeft,           *primRight,
                     *magneticFieldsLeft, *magneticFieldsRight,
                     *geomCenter,         *geomCenter,
                     directions::X3,
                     *fluidFluxesX3, *magneticFluxesX3,
                     numReadsRiemann, numWritesRiemann
                    );
      numReads  += numReadsRiemann;
      numWrites += numWritesRiemann;
    }

    /* 12) Constrained transport : Compute electric fields at the corners
     */
    computeEdgeElectricFields(numReadsCT, numWritesCT);
    numReads  += numReadsCT;
    numWrites += numWritesCT;
  }

  for (int var=0; var < primFlux.numVars; var++)
  {
    double filter1D[] = {1, -1, 0}; /* Forward difference */
 
    array filterX1 = array(3, 1, 1, 1, filter1D)/(XCoords->dX1);

    array dFluxX1_dX1 = convolve(fluidFluxesX1->vars[var], filterX1);
    divFluxes->vars[var] = dFluxX1_dX1;

    if (primFlux.dim >= 2)
    {
      array filterX2 = array(1, 3, 1, 1, filter1D)/(XCoords->dX2);
      array dFluxX2_dX2 = convolve(fluidFluxesX2->vars[var], filterX2);
      divFluxes->vars[var] += dFluxX2_dX2;

      if (primFlux.dim == 3)
      {
        array filterX3 = array(1, 1, 3, 1, filter1D)/(XCoords->dX3);
        array dFluxX3_dX3 = convolve(fluidFluxesX3->vars[var], filterX3);
        divFluxes->vars[var] += dFluxX3_dX3;
      }
    }

    divFluxes->vars[var].eval();
  }

}
