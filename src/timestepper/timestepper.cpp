#include "timestepper.hpp"

timeStepper::timeStepper(const int N1, 
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
{
  this->time = time;
  this->dt = dt;
  this->N1 = N1;
  this->N2 = N2;
  this->N3 = N3;
  this->numVars = numVars;

  this->boundaryLeft   = boundaryLeft;
  this->boundaryRight  = boundaryRight;
  this->boundaryTop    = boundaryTop;
  this->boundaryBottom = boundaryBottom;
  this->boundaryFront  = boundaryFront;
  this->boundaryBack   = boundaryBack;

  int periodicBoundariesX1;
  int periodicBoundariesX2;
  int periodicBoundariesX3;

  if (   boundaryLeft  == boundaries::PERIODIC
      || boundaryRight == boundaries::PERIODIC
     )
  {
    periodicBoundariesX1 = 1;
  }

  if (   boundaryTop    == boundaries::PERIODIC
      || boundaryBottom == boundaries::PERIODIC
     )
  {
    periodicBoundariesX2 = 1;
  }

  if (   boundaryFront == boundaries::PERIODIC
      || boundaryBack  == boundaries::PERIODIC
     )
  {
    periodicBoundariesX3 = 1;
  }

  prim         = new grid(N1, N2, N3,
                          dim, numVars, numGhost,
                          periodicBoundariesX1,
                          periodicBoundariesX2,
                          periodicBoundariesX3
                         );

  primHalfStep = new grid(N1, N2, N3,
                          dim, numVars, numGhost,
                          periodicBoundariesX1,
                          periodicBoundariesX2,
                          periodicBoundariesX3
                         );

  primOld      = new grid(N1, N2, N3,
                          dim, numVars, numGhost,
                          periodicBoundariesX1,
                          periodicBoundariesX2,
                          periodicBoundariesX3
                         );

  cons         = new grid(N1, N2, N3,
                          dim, numVars, numGhost,
                          periodicBoundariesX1,
                          periodicBoundariesX2,
                          periodicBoundariesX3
                         );

  consOld      = new grid(N1, N2, N3,
                          dim, numVars, numGhost,
                          periodicBoundariesX1,
                          periodicBoundariesX2,
                          periodicBoundariesX3
                         );

  sourcesExplicit    = new grid(N1, N2, N3,
                                dim, numVars, numGhost,
                                periodicBoundariesX1,
                                periodicBoundariesX2,
                                periodicBoundariesX3
                               );

  sourcesImplicit    = new grid(N1, N2, N3,
                                dim, numVars, numGhost,
                                periodicBoundariesX1,
                                periodicBoundariesX2,
                                periodicBoundariesX3
                               );

  sourcesImplicitOld = new grid(N1, N2, N3,
                                dim, numVars, numGhost,
                                periodicBoundariesX1,
                                periodicBoundariesX2,
                                periodicBoundariesX3
                               );

  sourcesTimeDer     = new grid(N1, N2, N3,
                                dim, numVars, numGhost,
                                periodicBoundariesX1,
                                periodicBoundariesX2,
                                periodicBoundariesX3
                               );

  primLeft  = new grid(N1, N2, N3,
                       dim, numVars, numGhost,
                       periodicBoundariesX1,
                       periodicBoundariesX2,
                       periodicBoundariesX3
                      );

  primRight = new grid(N1, N2, N3,
                       dim, numVars, numGhost,
                       periodicBoundariesX1,
                       periodicBoundariesX2,
                       periodicBoundariesX3
                      );

  fluidFluxesX1  = new grid(N1, N2, N3,
                            dim, numVars, numGhost,
                            periodicBoundariesX1,
                            periodicBoundariesX2,
                            periodicBoundariesX3
                           );

  fluidFluxesX2  = new grid(N1, N2, N3,
                            dim, numVars, numGhost,
                            periodicBoundariesX1,
                            periodicBoundariesX2,
                            periodicBoundariesX3
                           );

  fluidFluxesX3  = new grid(N1, N2, N3,
                            dim, numVars, numGhost,
                            periodicBoundariesX1,
                            periodicBoundariesX2,
                            periodicBoundariesX3
                           );

  magneticFluxesX1  = new grid(N1, N2, N3,
                               dim, 3, numGhost,
                               periodicBoundariesX1,
                               periodicBoundariesX2,
                               periodicBoundariesX3
                              );

  magneticFluxesX2  = new grid(N1, N2, N3,
                               dim, 3, numGhost,
                               periodicBoundariesX1,
                               periodicBoundariesX2,
                               periodicBoundariesX3
                              );

  magneticFluxesX3  = new grid(N1, N2, N3,
                               dim, 3, numGhost,
                               periodicBoundariesX1,
                               periodicBoundariesX2,
                               periodicBoundariesX3
                              );

  divFluxes = new grid(N1, N2, N3,
                       dim, numVars, numGhost,
                       periodicBoundariesX1,
                       periodicBoundariesX2,
                       periodicBoundariesX3
                      );

  divB = new grid(N1, N2, N3,
                  dim, 1, numGhost,
                  periodicBoundariesX1,
                  periodicBoundariesX2,
                  periodicBoundariesX3
                 );

  /* Primary magnetic field variables, indicated by their respective
   * locations */
  B1Left = new grid(N1, N2, N3,
                    dim, 1, numGhost,
                    periodicBoundariesX1,
                    periodicBoundariesX2,
                    periodicBoundariesX3
                   );

  B2Bottom = new grid(N1, N2, N3,
                      dim, 1, numGhost,
                      periodicBoundariesX1,
                      periodicBoundariesX2,
                      periodicBoundariesX3
                     );

  B3Back = new grid(N1, N2, N3,
                    dim, 1, numGhost,
                    periodicBoundariesX1,
                    periodicBoundariesX2,
                    periodicBoundariesX3
                   );

  /* Cell centered magnetic fields, computed using interpolation from 
   * face-centered fields */
  B1Center = new grid(N1, N2, N3,
                      dim, 1, numGhost,
                      periodicBoundariesX1,
                      periodicBoundariesX2,
                      periodicBoundariesX3
                     );

  B2Center = new grid(N1, N2, N3,
                      dim, 1, numGhost,
                      periodicBoundariesX1,
                      periodicBoundariesX2,
                      periodicBoundariesX3
                     );

  B3Center = new grid(N1, N2, N3,
                      dim, 1, numGhost,
                      periodicBoundariesX1,
                      periodicBoundariesX2,
                      periodicBoundariesX3
                     );

  /* Reconstructed face-centered magnetic fields, computed using
   * reconstruction of the cell-centered magnetic fields */
  B1Bottom = new grid(N1, N2, N3,
                      dim, 1, numGhost,
                      periodicBoundariesX1,
                      periodicBoundariesX2,
                      periodicBoundariesX3
                     );

  B1Top = new grid(N1, N2, N3,
                   dim, 1, numGhost,
                   periodicBoundariesX1,
                   periodicBoundariesX2,
                   periodicBoundariesX3
                  );

  B1Back = new grid(N1, N2, N3,
                    dim, 1, numGhost,
                    periodicBoundariesX1,
                    periodicBoundariesX2,
                    periodicBoundariesX3
                   );

  B1Front = new grid(N1, N2, N3,
                     dim, 1, numGhost,
                     periodicBoundariesX1,
                     periodicBoundariesX2,
                     periodicBoundariesX3
                    );

  B2Left = new grid(N1, N2, N3,
                    dim, 1, numGhost,
                    periodicBoundariesX1,
                    periodicBoundariesX2,
                    periodicBoundariesX3
                   );

  B2Right = new grid(N1, N2, N3,
                     dim, 1, numGhost,
                     periodicBoundariesX1,
                     periodicBoundariesX2,
                     periodicBoundariesX3
                    );

  B2Back = new grid(N1, N2, N3,
                    dim, 1, numGhost,
                    periodicBoundariesX1,
                    periodicBoundariesX2,
                    periodicBoundariesX3
                   );

  B2Front = new grid(N1, N2, N3,
                     dim, 1, numGhost,
                     periodicBoundariesX1,
                     periodicBoundariesX2,
                     periodicBoundariesX3
                    );

  B3Left = new grid(N1, N2, N3,
                    dim, 1, numGhost,
                    periodicBoundariesX1,
                    periodicBoundariesX2,
                    periodicBoundariesX3
                   );

  B3Right = new grid(N1, N2, N3,
                     dim, 1, numGhost,
                     periodicBoundariesX1,
                     periodicBoundariesX2,
                     periodicBoundariesX3
                    );

  B3Bottom = new grid(N1, N2, N3,
                      dim, 1, numGhost,
                      periodicBoundariesX1,
                      periodicBoundariesX2,
                      periodicBoundariesX3
                     );

  B3Top = new grid(N1, N2, N3,
                   dim, 1, numGhost,
                   periodicBoundariesX1,
                   periodicBoundariesX2,
                   periodicBoundariesX3
                  );

  /* Primary electric field variables, indicated by their locations */
  E1BottomBack  = new grid(N1, N2, N3,
                           dim, 1, numGhost,
                           periodicBoundariesX1,
                           periodicBoundariesX2,
                           periodicBoundariesX3
                          );

  E2LeftBack  = new grid(N1, N2, N3,
                         dim, 1, numGhost,
                         periodicBoundariesX1,
                         periodicBoundariesX2,
                         periodicBoundariesX3
                        );

  E3LeftBottom  = new grid(N1, N2, N3,
                           dim, 1, numGhost,
                           periodicBoundariesX1,
                           periodicBoundariesX2,
                           periodicBoundariesX3
                          );

  /* Auxiliary grids needed for inputs into the Riemann solver */
  magneticFieldsLeft  = new grid(N1, N2, N3,
                                 dim, 3, numGhost,
                                 periodicBoundariesX1,
                                 periodicBoundariesX2,
                                 periodicBoundariesX3
                                );

  magneticFieldsRight  = new grid(N1, N2, N3,
                                  dim, 3, numGhost,
                                  periodicBoundariesX1,
                                  periodicBoundariesX2,
                                  periodicBoundariesX3
                                 );

  magneticFieldsCenter  = new grid(N1, N2, N3,
                                   dim, 3, numGhost,
                                   periodicBoundariesX1,
                                   periodicBoundariesX2,
                                   periodicBoundariesX3
                                  );

  XCoords = new coordinatesGrid(N1, N2, N3,
                                dim, numGhost,
                                X1Start, X1End,
                                X2Start, X2End,
                                X3Start, X3End
                               );

  XCoords->setXCoords(locations::LEFT);
  geomLeft    = new geometry(metric,
                             blackHoleSpin,
                             hSlope, 
                             *XCoords
                            );

  XCoords->setXCoords(locations::RIGHT);
  geomRight   = new geometry(metric,
                             blackHoleSpin,
                             hSlope,
                             *XCoords
                            );

  XCoords->setXCoords(locations::BOTTOM);
  geomBottom  = new geometry(metric,
                             blackHoleSpin,
                             hSlope,
                             *XCoords
                            );

  XCoords->setXCoords(locations::TOP);
  geomTop     = new geometry(metric,
                             blackHoleSpin,
                             hSlope, 
                             *XCoords
                            );

  XCoords->setXCoords(locations::CENTER);
  geomCenter  = new geometry(metric,
                             blackHoleSpin,
                             hSlope,
                             *XCoords
                            );
  geomCenter->computeConnectionCoeffs();
  /* XCoords set to locations::CENTER */

  int numReads, numWrites;
  elem          = new fluidElement(*prim, *magneticFieldsCenter,
                                   *geomCenter,
                                   numReads, numWrites
                                  ); /* n+1   */
  elemOld       = new fluidElement(*primOld, *magneticFieldsCenter,
                                   *geomCenter,
                                   numReads, numWrites
                                  ); /* n     */
  elemHalfStep  = new fluidElement(*primHalfStep, *magneticFieldsCenter,
                                   *geomCenter,
                                   numReads, numWrites
                                  ); /* n+1/2 */

  riemann = new riemannSolver(*prim, *magneticFieldsCenter, *geomCenter);

  /* Data structures needed for the nonlinear solver */
  residual        = new grid(N1, N2, N3,
                             dim, numVars, numGhost,
                             periodicBoundariesX1,
                             periodicBoundariesX2,
                             periodicBoundariesX3
                            );

  residualPlusEps  = new grid(N1, N2, N3,
                              dim, numVars, numGhost,
                              periodicBoundariesX1,
                              periodicBoundariesX2,
                              periodicBoundariesX3
                             );

  primGuessPlusEps = new grid(N1, N2, N3,
                              dim, numVars, numGhost,
                              periodicBoundariesX1,
                              periodicBoundariesX2,
                              periodicBoundariesX3
                             );

  primGuessLineSearchTrial = new grid(N1, N2, N3,
                                      dim, numVars, numGhost,
                                      periodicBoundariesX1,
                                      periodicBoundariesX2,
                                      periodicBoundariesX3
                                     );

  /* The grid data structure arranges data in Struct of Arrays format. Need to
   * rearrange to Array of Structs format in order to solve the linear system Ax
   * = b */
  array zero = 0.*residual->vars[0];

  /* Arrays are copy-on-write. Hence, can use residual->varsSoA to initialize */
  residualSoA  = residual->varsSoA;

  /* Jacobian \partial residual/ \prim in Struct of Arrays format */
  jacobianSoA  = af::constant(1., residual->vars[0].dims(0),
                                  residual->vars[0].dims(1),
                                  residual->vars[0].dims(2),
                                  numVars * numVars,
                                  f64
                             );

  /* Correction dP_k in P_{k+1} = P_k + lambda*dP_k in Array of Structs format */
  deltaPrimAoS = af::constant(0., 
                              numVars,
                              residual->vars[0].dims(0),
                              residual->vars[0].dims(1),
                              residual->vars[0].dims(2),
                              f64
                             );

  /* Steplength lambda in P_{k+1} = P_k + lambda*dP_k */ 
  stepLength  = zero;
                     
  int N1Total = residual->N1Total;
  int N2Total = residual->N2Total;
  int N3Total = residual->N3Total;

  AHostPtr = new double [numVars*numVars*N1Total*N2Total*N3Total];
  bHostPtr = new double [numVars*N1Total*N2Total*N3Total];
  xHostPtr = new double [numVars*N1Total*N2Total*N3Total];

  int a,b;
  initialConditions(a,b);
  XCoords->setXCoords(locations::CENTER);
  /* Reset XCoords to CENTER in case it has been changed in initalConditions()
   */
}

timeStepper::~timeStepper()
{
  delete prim, primHalfStep, primOld;
  delete cons, consOld;
  delete sourcesExplicit, sourcesImplicit, sourcesImplicitOld, sourcesTimeDer;
  delete primLeft, primRight;
  delete fluidFluxesX1,    fluidFluxesX2,    fluidFluxesX3;
  delete magneticFluxesX1, magneticFluxesX2, magneticFluxesX3;
  delete divFluxes;
  delete divB;

  delete B1Left, B2Bottom, B3Back;
  delete B1Center, B2Center, B3Center;
  delete B1Bottom, B1Top,   B1Back,   B1Front;
  delete B2Left,   B2Right, B2Back,   B2Front;
  delete B3Left,   B3Right, B3Bottom, B3Top;
  delete magneticFieldsLeft, magneticFieldsRight, magneticFieldsCenter;
  delete E1BottomBack, E2LeftBack, E3LeftBottom;

  delete elem, elemOld, elemHalfStep;
  delete riemann;
  delete geomLeft, geomRight, geomBottom, geomTop, geomCenter;

  delete primGuessLineSearchTrial;
  delete primGuessPlusEps;
  delete residual;
  delete residualPlusEps;

  delete[] AHostPtr;
  delete[] bHostPtr;
  delete[] xHostPtr;
}
