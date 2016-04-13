#include "timestepper.hpp"

void timeStepper::computeCellCenteredMagneticFields
  (int &numReads, 
   int &numWrites
  )
{
  B1Center->vars[0] = 0.5*(  B1Left->vars[0] 
                           + shift(B1Left->vars[0], -1, 0, 0)
                         );
  B2Center->vars[0] = 0.5*(  B2Bottom->vars[0] 
                           + shift(B2Bottom->vars[0], 0, -1, 0)
                          );
  B3Center->vars[0] = 0.5*(  B3Back->vars[0] 
                           + shift(B3Back->vars[0], 0, 0, -1)
                          );
}

void timeStepper::computeEdgeElectricFields(int &numReads,
                                            int &numWrites
                                           )
{
  if (magneticFluxesX1->dim >= 2)
  {
    array E3Left   = -magneticFluxesX1->vars[vars::B2];
    array E3Bottom =  magneticFluxesX2->vars[vars::B1];
    array E3LeftShiftedDown   = shift(E3Left,   0, 1, 0);
    array E3BottomShiftedLeft = shift(E3Bottom, 1, 0, 0);

    E3LeftBottom->vars[0] =
      0.25*(  E3Left   + E3LeftShiftedDown
            + E3Bottom + E3BottomShiftedLeft
           );
    E3LeftBottom->vars[0].eval();


    if (magneticFluxesX1->dim == 3)
    {
      array E1Back   =  magneticFluxesX3->vars[vars::B2];
      array E1Bottom = -magneticFluxesX2->vars[vars::B3];
      array E1BackShiftedDown   = shift(E1Back,   0, 1, 0);
      array E1BottomShiftedBack = shift(E1Bottom, 0, 0, 1);

      E1BottomBack->vars[0] =
        0.25*(  E1Back   + E1BackShiftedDown
              + E1Bottom + E1BottomShiftedBack
             );
      E1BottomBack->vars[0].eval();

      array E2Left =  magneticFluxesX1->vars[vars::B3];
      array E2Back = -magneticFluxesX3->vars[vars::B1];
      array E2LeftShiftedBack = shift(E2Left, 0, 0, 1);
      array E2BackShiftedLeft = shift(E2Back, 1, 0, 0);

      E2LeftBack->vars[0] =
        0.25*(  E2Left + E2LeftShiftedBack
              + E2Back + E2BackShiftedLeft
             );
      E2LeftBack->vars[0].eval();
    }
  }
}

void timeStepper::computeDivB(const grid &prim,
                              int &numReads,
                              int &numWrites
                             )
{
  array B1Right = shift(B1Left->vars[0],   -1,  0,  0);
  array B2Top   = shift(B2Bottom->vars[0],  0, -1,  0);
  array B3Front = shift(B3Back->vars[0],    0,  0, -1);

  double dX1 = XCoords->dX1;
  double dX2 = XCoords->dX2;
  double dX3 = XCoords->dX3;

  if (prim.dim >= 2)
  {
    divB->vars[0] =
      ( (geomRight->g * B1Right) - (geomLeft->g   * B1Left->vars[0]  ) )/dX1
    + ( (geomTop->g   * B2Top  ) - (geomBottom->g * B2Bottom->vars[0]) )/dX2;

    if (prim.dim == 3)
    {
      divB->vars[0] += geomCenter->g*(B3Front - B3Back->vars[0])/dX3;
    }
  }
}
