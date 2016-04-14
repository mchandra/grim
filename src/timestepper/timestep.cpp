#include "timestepper.hpp"

void timeStepper::timeStep(int &numReads, int &numWrites)
{
  PetscPrintf(PETSC_COMM_WORLD, "Time = %f, dt = %f\n", time, dt);
  /* First take a half step */
  PetscPrintf(PETSC_COMM_WORLD, "  ---Half step--- \n");

  currentStep = timeStepperSwitches::HALF_STEP;
  /* 1) Apply boundary conditions on primOld, B1LeftOld, B2BottomOld and
   * B3BackOld */
  boundaries::applyBoundaryConditions(boundaryLeft, boundaryRight,
                                      boundaryTop,  boundaryBottom,
                                      boundaryFront, boundaryBack,
                                      *primOld
                                     );
  boundaries::applyBoundaryConditions(boundaryLeft, boundaryRight,
                                      boundaryTop,  boundaryBottom,
                                      boundaryFront, boundaryBack,
                                      *B1LeftOld
                                     );
  boundaries::applyBoundaryConditions(boundaryLeft, boundaryRight,
                                      boundaryTop,  boundaryBottom,
                                      boundaryFront, boundaryBack,
                                      *B2BottomOld
                                     );
  boundaries::applyBoundaryConditions(boundaryLeft, boundaryRight,
                                      boundaryTop,  boundaryBottom,
                                      boundaryFront, boundaryBack,
                                      *B3BackOld
                                     );
  setProblemSpecificBCs(numReads,numWrites);

  /* 2) Compute B1Center, B2Center and B3Center using B1LeftOld, B2BottomOld and
   * B3BackOld */
  computeCellCenteredMagneticFields(*B1LeftOld, *B2BottomOld, *B3BackOld,
                                    numReads, numWrites
                                   );
  magneticFieldsCenter->vars[vars::B1] = B1Center->vars[0];
  magneticFieldsCenter->vars[vars::B2] = B2Center->vars[0];
  magneticFieldsCenter->vars[vars::B3] = B3Center->vars[0];

  /* 3) Set elemOld using primOld and B1Center, B2Center, B3Center which have
   * been computed using B1LeftOld, B2BottomOld and B3BackOld */
  int numReadsElemSet, numWritesElemSet;
  int numReadsComputeFluxes, numWritesComputeFluxes;
  elemOld->set(*primOld, *magneticFieldsCenter, *geomCenter,
               numReadsElemSet, numWritesElemSet
              );
  elemOld->computeFluidFluxes(*geomCenter, 0, *consOld, 
                              numReadsComputeFluxes, numWritesComputeFluxes
                             );
  numReads  = numReadsElemSet  + numReadsComputeFluxes;
  numWrites = numWritesElemSet + numWritesComputeFluxes; 

  /* 4) Compute EMHD gradients gradT and graduCov using elemOld */
  if (params::viscosity || params::conduction)
  {
    int numReadsEMHDGradients, numWritesEMHDGradients;
    double dX[3];
    dX[0] = XCoords->dX1;
    dX[1] = XCoords->dX2;
    dX[2] = XCoords->dX3;
    elemOld->computeEMHDGradients(*geomCenter, dX,
                                  numReadsEMHDGradients,
                                  numWritesEMHDGradients
                                 );
    numReads  += numReadsEMHDGradients;
    numWrites += numWritesEMHDGradients;
  }
  
  /* 5) Compute div.F and the edge electric fields using primOld, B1LeftOld,
   * B2BottomOld and B3BackOld */
  int numReadsDivFluxes, numWritesDivFluxes;
  computeDivOfFluxes(*primOld, 
                     *B1LeftOld, *B2BottomOld, *B3BackOld,
                     numReadsDivFluxes, numWritesDivFluxes
                    );
  numReads  += numReadsDivFluxes;
  numWrites += numWritesDivFluxes;

  /* 6) Solve the induction equation for B1Left, B2Bottom and B3Back at n+1/2 */
  int numReadsBFields = 0, numWritesBFields = 0;
  timeStepBFields(0.5*dt, numReadsBFields, numWritesBFields);
  numReads  += numReadsBFields;
  numWrites += numWritesBFields;

  /* 7) Compute new B1Center, B2Center, B3Center, using the B1Left, B2Bottom,
   * B3Back at n+1/2 */
  computeCellCenteredMagneticFields(*B1Left, *B2Bottom, *B3Back,
                                    numReads, numWrites
                                   );
  magneticFieldsCenter->vars[vars::B1] = B1Center->vars[0];
  magneticFieldsCenter->vars[vars::B2] = B2Center->vars[0];
  magneticFieldsCenter->vars[vars::B3] = B3Center->vars[0];

  /* 8) Solve the fluid conservation equations. Note that in the residual
   * computation, elemOld uses magneticFieldsCenter computed using B1LeftOld,
   * B2BottomOld and B3BackOld, and elem uses magneticFieldsCenter computed
   * using B1Left, B2Bottom and B3Back */
  /* Set a guess for prim */
  for (int var=0; var < prim->numVars; var++)
  {
    prim->vars[var] = primOld->vars[var];
    numReads  += 1;
    numWrites += 1;
  }
  /* Solve dU/dt + div.F - S = 0 to get prim at n+1/2 */
  solve(*prim);

  /* 9) Copy solution to primHalfStep, B1LeftHalfStep, B2BottomHalfStep and
   * B3BackHalfStep */
  for (int var=0; var < prim->numVars; var++)
  {
    primHalfStep->vars[var] = prim->vars[var];
    numReads  += 1;
    numWrites += 1;
  }
  B1LeftHalfStep->vars[0]   = B1Left->vars[0];
  B2BottomHalfStep->vars[0] = B2Bottom->vars[0];
  B3BackHalfStep->vars[0]   = B3Back->vars[0];

  /* 10) Communicate the ghost zone data */
  primHalfStep->communicate();
  B1LeftHalfStep->communicate();
  B2BottomHalfStep->communicate();
  B3BackHalfStep->communicate();

  /* 11) Half Step Diagnostics */
  halfStepDiagnostics(numReads,numWrites);
  /* Half step complete */

  /* Now take the full step */
  PetscPrintf(PETSC_COMM_WORLD, "  ---Full step--- \n");

  currentStep = timeStepperSwitches::FULL_STEP;
  /* 12) Apply boundary conditions on primHalfStep, B1LeftHalfStep,
   * B2BottomHalfStep and B3BackHalfStep */
  boundaries::applyBoundaryConditions(boundaryLeft, boundaryRight,
                                      boundaryTop,  boundaryBottom,
                                      boundaryFront, boundaryBack,
                                      *primHalfStep
                                     );
  boundaries::applyBoundaryConditions(boundaryLeft, boundaryRight,
                                      boundaryTop,  boundaryBottom,
                                      boundaryFront, boundaryBack,
                                      *B1LeftHalfStep
                                     );
  boundaries::applyBoundaryConditions(boundaryLeft, boundaryRight,
                                      boundaryTop,  boundaryBottom,
                                      boundaryFront, boundaryBack,
                                      *B2BottomHalfStep
                                     );
  boundaries::applyBoundaryConditions(boundaryLeft, boundaryRight,
                                      boundaryTop,  boundaryBottom,
                                      boundaryFront, boundaryBack,
                                      *B3BackHalfStep
                                     );
  setProblemSpecificBCs(numReads,numWrites);

  /* 13) Compute B1Center, B2Center and B3Center using B1LeftHalfStep, 
   * B2BottomHalfStep and B3BackHalfStep */
  computeCellCenteredMagneticFields(*B1LeftHalfStep,
                                    *B2BottomHalfStep,
                                    *B3BackHalfStep,
                                    numReads, numWrites
                                   );
  magneticFieldsCenter->vars[vars::B1] = B1Center->vars[0];
  magneticFieldsCenter->vars[vars::B2] = B2Center->vars[0];
  magneticFieldsCenter->vars[vars::B3] = B3Center->vars[0];

  /* 14) Set elemHalfStep using primHalfStep and B1Center, B2Center, B3Center
   * which have been computed using B1LeftHalfStep, B2BottomHalfStep and 
   * B3BackHalfStep */
  elemHalfStep->set(*primHalfStep, *magneticFieldsCenter,
                    *geomCenter,
                    numReadsElemSet, numWritesElemSet
                   );
  numReads  += numReadsElemSet;
  numWrites += numWritesElemSet; 

  /* 15) Compute EMHD gradients gradT and graduCov using elemOld */
  if (params::viscosity || params::conduction)
  {
    int numReadsEMHDGradients, numWritesEMHDGradients;
    double dX[3];
    dX[0] = XCoords->dX1;
    dX[1] = XCoords->dX2;
    dX[2] = XCoords->dX3;
    elemHalfStep->computeEMHDGradients(*geomCenter, dX,
                                       numReadsEMHDGradients,
                                       numWritesEMHDGradients
                                      );
    numReads  += numReadsEMHDGradients;
    numWrites += numWritesEMHDGradients;
  }

  /* 16) Compute div.F and the edge electric fields using primHalfStep,
   * B1LeftHalfStep, B2BottomHalfStep and B3BackHalfStep */
  computeDivOfFluxes(*primHalfStep, 
                     *B1LeftHalfStep, *B2BottomHalfStep, *B3BackHalfStep,
                     numReadsDivFluxes, numWritesDivFluxes
                    );
  numReads  += numReadsDivFluxes;
  numWrites += numWritesDivFluxes;

  /* 17) Solve the induction equation for B1Left, B2Bottom and B3Back at n+1 */
  timeStepBFields(dt, numReadsBFields, numWritesBFields);
  numReads  += numReadsBFields;
  numWrites += numWritesBFields;

  /* 18) Compute new B1Center, B2Center, B3Center, using the B1Left, B2Bottom,
   * B3Back at n+1 */
  computeCellCenteredMagneticFields(*B1Left, *B2Bottom, *B3Back,
                                    numReads, numWrites
                                   );
  magneticFieldsCenter->vars[vars::B1] = B1Center->vars[0];
  magneticFieldsCenter->vars[vars::B2] = B2Center->vars[0];
  magneticFieldsCenter->vars[vars::B3] = B3Center->vars[0];

  /* 19) Solve dU/dt + div.F - S = 0 to get prim at n+1. NOTE: prim already has
   * primHalfStep as a guess */
  solve(*prim);

  /* 20) Copy solution to primOld, B1LeftOld, B2BottomOld, and B3BackOld */
  for (int var=0; var < prim->numVars; var++)
  {
    primOld->vars[var] = prim->vars[var];
    numReads  += 1;
    numWrites += 1;
  }
  B1LeftOld->vars[0]   = B1Left->vars[0];
  B2BottomOld->vars[0] = B2Bottom->vars[0];
  B3BackOld->vars[0]   = B3Back->vars[0];

  /* 21) Communicate the ghost zone data */
  primOld->communicate();
  B1LeftOld->communicate();
  B2BottomOld->communicate();
  B3BackOld->communicate();

  /* 22) Compute full step diagnostics */
  time += dt;
  fullStepDiagnostics(numReads,numWrites);

  /* done */
}

void timeStepper::timeStepBFields(const double dt,
                                  int &numReads,
                                  int &numWrites
                                 )
{
  double dX1 = XCoords->dX1;
  double dX2 = XCoords->dX2;
  double dX3 = XCoords->dX3;

  array E3LeftBottomShiftedUp  = shift(E3LeftBottom->vars[0], 0, -1,  0);
  array E2LeftBackShiftedFront = shift(E2LeftBack->vars[0],   0,  0, -1);
  B1Left->vars[0] = 
    B1LeftOld->vars[0] - dt*(  E2LeftBack->vars[0]    * dX2
                             + E3LeftBottomShiftedUp  * dX3
                             - E2LeftBackShiftedFront * dX2
                             - E3LeftBottom->vars[0]  * dX3
                            ) / (geomLeft->g * dX2 * dX3);
  B1Left->vars[0].eval();

  array E1BottomBackShiftedFront = shift(E1BottomBack->vars[0],  0,  0, -1);
  array E3LeftBottomShiftedRight = shift(E3LeftBottom->vars[0], -1,  0,  0);
  B2Bottom->vars[0] =
    B2BottomOld->vars[0] - dt*(  E1BottomBackShiftedFront * dX1
                               - E3LeftBottomShiftedRight * dX3
                               - E1BottomBack->vars[0]    * dX1
                               + E3LeftBottom->vars[0]    * dX3
                              ) / (geomBottom->g * dX1 * dX3);
  B2Bottom->vars[0].eval();

  array E2LeftBackShiftedRight = shift(E2LeftBack->vars[0],   -1,  0,  0);
  array E1BottomBackShiftedUp  = shift(E1BottomBack->vars[0],  0,  0, -1);
  B3Back->vars[0] =
    B3BackOld->vars[0] - dt*(  E2LeftBackShiftedRight * dX2
                             - E1BottomBackShiftedUp  * dX1
                             - E2LeftBack->vars[0]    * dX2
                             + E1BottomBack->vars[0]  * dX1
                            ) / (geomCenter->g * dX1 * dX2);
}
