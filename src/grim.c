#include "grim.h"

int main(int argc, char **argv)
{ 
  PetscInitialize(&argc, &argv, NULL, help);
  
  struct timeStepper ts;
  timeStepperInit(&ts);

  while (ts.t + ts.dt < FINAL_TIME)
  {
    timeStep(&ts);
  }
  
  /* One final step */
  if (ts.t < FINAL_TIME)
  {
    ts.dt = FINAL_TIME - ts.t;
    timeStep(&ts);
  }

  timeStepperDestroy(&ts);

//  REAL rho = 1.;
//  REAL temperature = 100.;
//  REAL primVars[DOF];
//  primVars[ALPHA] = getAlpha(rho, temperature);
//  primVars[A0]    = getA0(temperature);
//  primVars[U1]    = -0.170*1e-4;
//  primVars[U2]    = 0.;
//  primVars[U3]    = 0.;
//  primVars[B1]    = 0.;
//  primVars[B2]    = 0.;
//  primVars[B3]    = 0.;
//
//  struct gridZone zone;
//  setGridZone(0, 0,
//              0, 0,
//              0, 0,
//              32, 32,
//              &zone);
//
//  REAL XCoords[NDIM];
//  getXCoords(&zone, CENTER, XCoords);
//  struct fluidElement elem;
//  struct geometry geom;
//
//  setGeometry(XCoords, &geom);
//  setFluidElement(primVars, &geom, &elem);
//
//  REAL momentsUpUp[NUM_ALL_COMPONENTS], moments[NUM_ALL_COMPONENTS];
//  fixedQuadIntegration5Moments(&elem, &geom, temperature, momentsUpUp);
//
//  for (int mu=0; mu<NDIM; mu++)
//  {
//    for (int nu=0; nu<NDIM; nu++)
//    {
//      moments[T_UP_DOWN(mu, nu)] = 0.;
//      for (int alpha=0; alpha<NDIM; alpha++)
//      {
//        moments[T_UP_DOWN(mu, nu)] +=
//	        momentsUpUp[T_UP_UP(mu, alpha)]*geom.gCov[alpha][nu];
//      }
//    }
//  }
//
//  PetscPrintf(PETSC_COMM_WORLD, "T^0_0 = %f\n", moments[T_UP_DOWN(0, 0)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^0_1 = %f\n", moments[T_UP_DOWN(0, 1)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^0_2 = %f\n", moments[T_UP_DOWN(0, 2)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^0_3 = %f\n", moments[T_UP_DOWN(0, 3)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^1_0 = %f\n", moments[T_UP_DOWN(1, 0)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^1_1 = %f\n", moments[T_UP_DOWN(1, 1)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^1_2 = %f\n", moments[T_UP_DOWN(1, 2)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^1_3 = %f\n", moments[T_UP_DOWN(1, 3)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^2_0 = %f\n", moments[T_UP_DOWN(2, 0)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^2_1 = %f\n", moments[T_UP_DOWN(2, 1)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^2_2 = %f\n", moments[T_UP_DOWN(2, 2)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^2_3 = %f\n", moments[T_UP_DOWN(2, 3)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^3_0 = %f\n", moments[T_UP_DOWN(3, 0)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^3_1 = %f\n", moments[T_UP_DOWN(3, 1)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^3_2 = %f\n", moments[T_UP_DOWN(3, 2)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^3_3 = %f\n", moments[T_UP_DOWN(3, 3)]);
//
//
//  REAL pressure = rho * temperature;
//  REAL uu = pressure/(ADIABATIC_INDEX - 1.);
//  REAL bCov[NDIM], bSqr, uCov[NDIM];
//
//  bSqr = getbSqr(&elem, &geom);
//
//  conToCov(elem.uCon, &geom, uCov);
//  conToCov(elem.bCon, &geom, bCov);
//
//  REAL TUpDownIdeal[20];
//
//  for (int mu=0; mu<NDIM; mu++)
//  {
//    for (int nu=0; nu<NDIM; nu++)
//    {
//      TUpDownIdeal[T_UP_DOWN(mu,nu)] =   
//                          (  rho + uu + pressure + bSqr
//                          )*elem.uCon[mu]*uCov[nu]
//
//                        + (pressure + 0.5*bSqr)*DELTA(mu, nu)
//
//                        - elem.bCon[mu]*bCov[nu];      
//    }
//  }
//
//  PetscPrintf(PETSC_COMM_WORLD, "\n");
//  PetscPrintf(PETSC_COMM_WORLD, "-------------------------\n");
//  PetscPrintf(PETSC_COMM_WORLD, "\n");
//
//  PetscPrintf(PETSC_COMM_WORLD, "T^0_0 = %f\n", TUpDownIdeal[T_UP_DOWN(0, 0)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^0_1 = %f\n", TUpDownIdeal[T_UP_DOWN(0, 1)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^0_2 = %f\n", TUpDownIdeal[T_UP_DOWN(0, 2)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^0_3 = %f\n", TUpDownIdeal[T_UP_DOWN(0, 3)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^1_0 = %f\n", TUpDownIdeal[T_UP_DOWN(1, 0)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^1_1 = %f\n", TUpDownIdeal[T_UP_DOWN(1, 1)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^1_2 = %f\n", TUpDownIdeal[T_UP_DOWN(1, 2)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^1_3 = %f\n", TUpDownIdeal[T_UP_DOWN(1, 3)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^2_0 = %f\n", TUpDownIdeal[T_UP_DOWN(2, 0)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^2_1 = %f\n", TUpDownIdeal[T_UP_DOWN(2, 1)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^2_2 = %f\n", TUpDownIdeal[T_UP_DOWN(2, 2)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^2_3 = %f\n", TUpDownIdeal[T_UP_DOWN(2, 3)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^3_0 = %f\n", TUpDownIdeal[T_UP_DOWN(3, 0)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^3_1 = %f\n", TUpDownIdeal[T_UP_DOWN(3, 1)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^3_2 = %f\n", TUpDownIdeal[T_UP_DOWN(3, 2)]);
//  PetscPrintf(PETSC_COMM_WORLD, "T^3_3 = %f\n", TUpDownIdeal[T_UP_DOWN(3, 3)]);

  PetscFinalize();  
  return(0);
}
