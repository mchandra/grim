#include "grim.h"

int main(int argc, char **argv)
{ 
  PetscInitialize(&argc, &argv, NULL, help);
  
  struct timeStepper ts;
  timeStepperInit(&ts);
  
  timeStepperDestroy(&ts);

  PetscFinalize();  
  return(0);
}