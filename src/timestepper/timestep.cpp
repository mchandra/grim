#include "timestepper.hpp"

void timeStepper::timeStep(int &numReads, int &numWrites)
{
  PetscPrintf(PETSC_COMM_WORLD, "Time = %f, dt = %f\n", time, dt);
  /* First take a half step */
  PetscPrintf(PETSC_COMM_WORLD, "  ---Half step--- \n");

  currentStep = timeStepperSwitches::HALF_STEP;
  /* Apply boundary conditions on primOld */
  boundaries::applyBoundaryConditions(boundaryLeft, boundaryRight,
                                      boundaryTop,  boundaryBottom,
                                      boundaryFront, boundaryBack,
                                      *primOld
                                     );
  boundaries::applyBoundaryConditions(boundaryLeft, boundaryRight,
                                      boundaryTop,  boundaryBottom,
                                      boundaryFront, boundaryBack,
                                      *B1Left
                                     );
  boundaries::applyBoundaryConditions(boundaryLeft, boundaryRight,
                                      boundaryTop,  boundaryBottom,
                                      boundaryFront, boundaryBack,
                                      *B2Bottom
                                     );
  boundaries::applyBoundaryConditions(boundaryLeft, boundaryRight,
                                      boundaryTop,  boundaryBottom,
                                      boundaryFront, boundaryBack,
                                      *B3Back
                                     );
  setProblemSpecificBCs(numReads,numWrites);

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
  
  int numReadsDivFluxes, numWritesDivFluxes;
  computeDivOfFluxes(*primOld, numReadsDivFluxes, numWritesDivFluxes);
  numReads  += numReadsDivFluxes;
  numWrites += numWritesDivFluxes;

  /* Set a guess for prim */
  for (int var=0; var < prim->numVars; var++)
  {
    prim->vars[var] = primOld->vars[var];
    numReads  += 1;
    numWrites += 1;
  }

  int numReadsBFields = 0, numWritesBFields = 0;
  timeStepBFields(0.5*dt, numReadsBFields, numWritesBFields);
  numReads  += numReadsBFields;
  numWrites += numWritesBFields;

  /* Solve dU/dt + div.F - S = 0 to get prim at n+1/2 */
  solve(*prim);

  /* Copy solution to primHalfStepGhosted. WARNING: Right now
   * primHalfStep->vars[var] points to prim->vars[var]. Might need to do a deep
   * copy. */
  for (int var=0; var < prim->numVars; var++)
  {
    primHalfStep->vars[var] = prim->vars[var];
    numReads  += 1;
    numWrites += 1;
  }
  primHalfStep->communicate();
  halfStepDiagnostics(numReads,numWrites);
  /* Half step complete */

  /* Now take the full step */
  PetscPrintf(PETSC_COMM_WORLD, "  ---Full step--- \n");

  currentStep = timeStepperSwitches::FULL_STEP;
  /* apply boundary conditions on primHalfStep */
  boundaries::applyBoundaryConditions(boundaryLeft, boundaryRight,
                                      boundaryTop,  boundaryBottom,
                                      boundaryFront, boundaryBack,
                                      *primHalfStep
                                     );
  boundaries::applyBoundaryConditions(boundaryLeft, boundaryRight,
                                      boundaryTop,  boundaryBottom,
                                      boundaryFront, boundaryBack,
                                      *B1Left
                                     );
  boundaries::applyBoundaryConditions(boundaryLeft, boundaryRight,
                                      boundaryTop,  boundaryBottom,
                                      boundaryFront, boundaryBack,
                                      *B2Bottom
                                     );
  boundaries::applyBoundaryConditions(boundaryLeft, boundaryRight,
                                      boundaryTop,  boundaryBottom,
                                      boundaryFront, boundaryBack,
                                      *B3Back
                                     );
  setProblemSpecificBCs(numReads,numWrites);

  elemHalfStep->set(*primHalfStep, *magneticFieldsCenter,
                    *geomCenter,
                    numReadsElemSet, numWritesElemSet
                   );
  numReads  += numReadsElemSet;
  numWrites += numWritesElemSet; 

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
  computeDivOfFluxes(*primHalfStep, numReadsDivFluxes, numWritesDivFluxes);
  numReads  += numReadsDivFluxes;
  numWrites += numWritesDivFluxes;

  timeStepBFields(dt, numReadsBFields, numWritesBFields);
  numReads  += numReadsBFields;
  numWrites += numWritesBFields;

  /* Solve dU/dt + div.F - S = 0 to get prim at n+1/2. NOTE: prim already has
   * primHalfStep as a guess */
  solve(*prim);

  /* Copy solution to primOldGhosted */
  for (int var=0; var < prim->numVars; var++)
  {
    primOld->vars[var] = prim->vars[var];
    numReads  += 1;
    numWrites += 1;
  }
  /* Compute diagnostics */
  primOld->communicate();
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
    B1Left->vars[0] - dt*(  E2LeftBack->vars[0]    * dX2
                          + E3LeftBottomShiftedUp  * dX3
                          - E2LeftBackShiftedFront * dX2
                          - E3LeftBottom->vars[0]  * dX3
                         ) / (geomLeft->g * dX2 * dX3);
  B1Left->vars[0].eval();

  array E1BottomBackShiftedFront = shift(E1BottomBack->vars[0],  0,  0, -1);
  array E3LeftBottomShiftedRight = shift(E3LeftBottom->vars[0], -1,  0,  0);
  B2Bottom->vars[0] =
    B2Bottom->vars[0] - dt*(  E1BottomBackShiftedFront * dX1
                            - E3LeftBottomShiftedRight * dX3
                            - E1BottomBack->vars[0]    * dX1
                            + E3LeftBottom->vars[0]    * dX3
                           ) / (geomBottom->g * dX1 * dX3);
  B2Bottom->vars[0].eval();

  array E2LeftBackShiftedRight = shift(E2LeftBack->vars[0],   -1,  0,  0);
  array E1BottomBackShiftedUp  = shift(E1BottomBack->vars[0],  0,  0, -1);
  B3Back->vars[0] =
    B3Back->vars[0] - dt*(  E2LeftBackShiftedRight * dX2
                          - E1BottomBackShiftedUp  * dX1
                          - E2LeftBack->vars[0]    * dX2
                          + E1BottomBack->vars[0]  * dX1
                         ) / (geomCenter->g * dX1 * dX2);
}
