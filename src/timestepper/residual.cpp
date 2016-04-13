#include "timestepper.hpp"

void timeStepper::computeResidual(const grid &primGuess,
                                  grid &residualGuess, 
                                  const bool computeExplicitTerms,
                                  int &numReads,
                                  int &numWrites
                                 )
{
  numReads = 0; numWrites = 0;
  int numReadsElemSet, numWritesElemSet;
  int numReadsComputeFluxes, numWritesComputeFluxes;
  elem->set(primGuess, *magneticFieldsCenter, *geomCenter,
            numReadsElemSet, numWritesElemSet
           );
  elem->computeFluidFluxes(*geomCenter, 0, *cons,
                           numReadsComputeFluxes,
                           numWritesComputeFluxes
                          );
  numReads  += numReadsElemSet  + numReadsComputeFluxes;
  numWrites += numWritesElemSet + numWritesComputeFluxes;

  if (currentStep == timeStepperSwitches::HALF_STEP)
  {
    if (computeExplicitTerms)
    {
      int numReadsExplicitSouces, numWritesExplicitSouces;
      elemOld->computeExplicitSources(*geomCenter, *sourcesExplicit,
                                      numReadsExplicitSouces, 
                                      numWritesExplicitSouces
                                     );
      numReads  += numReadsExplicitSouces;
      numWrites += numWritesExplicitSouces;
    }
    else
    {
      for (int var=0; var < numVars; var++)
      {
	      sourcesExplicit->vars[var] = 0.;
      }
    }

    int numReadsImplicitSources, numReadsTimeDerivSources;
    int numWritesImplicitSources, numWritesTimeDerivSources;
    elemOld->computeImplicitSources(*geomCenter, *sourcesImplicitOld,
				                            elemOld->tau,
                                    numReadsImplicitSources,
                                    numWritesImplicitSources
                                   );
    elem->computeImplicitSources(*geomCenter, *sourcesImplicit,
				                         elemOld->tau,
                                 numReadsImplicitSources,
                                 numWritesImplicitSources
                                );
    elemOld->computeTimeDerivSources(*geomCenter,
				                             *elemOld, *elem,
		                         		     dt/2,
                        				     *sourcesTimeDer,
                                     numReadsTimeDerivSources,
                                     numWritesTimeDerivSources
                                    );
    numReads  += 2*numReadsImplicitSources  + numReadsTimeDerivSources;
    numWrites += 2*numWritesImplicitSources + numWritesTimeDerivSources; 

    for (int var=0; var<residualGuess.numVars; var++)
    {
      residualGuess.vars[var] = 
        (cons->vars[var] - consOld->vars[var])/(dt/2.)
  	  + divFluxes->vars[var]
      + sourcesExplicit->vars[var]
  	  + 0.5*(sourcesImplicitOld->vars[var] + sourcesImplicit->vars[var])
	    + sourcesTimeDer->vars[var];
    }
    /* Reads:
     * -----
     *  cons[var], consOld[var], divFluxes[var]       : 3*numVars
     *  sourcesExplicit[var], sourcesTimeDer[var]     : 2*numVars
     *  sourcesImplicitOld[var], sourcesImplicit[var] : 2*numVars
     *
     * Writes:
     * ------
     * residualGuess[var] : numVars */
    numReads  += 7*residualGuess.numVars;

    /* Normalization of the residualGuess */
    //for (int var=0; var<vars::dof; var++)
    //  residualGuess.vars[var] = residualGuess.vars[var]/geomCenter->g;
    if (params::conduction)
    {
      if(params::highOrderTermsConduction)
      {
	      residualGuess.vars[vars::Q] *=
           elemOld->temperature 
	       * af::sqrt(elemOld->rho*elemOld->chi_emhd*elemOld->tau);

        numReads += 4;
      }
      else
      {
	      residualGuess.vars[vars::Q] *= elemOld->tau;
        numReads += 1;
      }
    }

    if (params::viscosity)
    {
      if(params::highOrderTermsViscosity)
      {
	      residualGuess.vars[vars::DP] *=
          af::sqrt(   elemOld->rho*elemOld->nu_emhd
		                * elemOld->temperature*elemOld->tau
		              );
        numReads += 4;
      }
      else
      {
	      residualGuess.vars[vars::DP] *= elemOld->tau;
        numReads += 1;
      }
    }

  } /* End of timeStepperSwitches::HALF_STEP */

  else if (currentStep == timeStepperSwitches::FULL_STEP)
  {
    if (computeExplicitTerms)
    {
      int numReadsExplicitSouces, numWritesExplicitSouces;
      elemHalfStep->computeExplicitSources(*geomCenter, *sourcesExplicit,
                                           numReadsExplicitSouces,
                                           numWritesExplicitSouces
                                          );
      numReads  += numReadsExplicitSouces;
      numWrites += numWritesExplicitSouces;
    }
    else
    {
      for (int var=0; var < numVars; var++)
      {
	      sourcesExplicit->vars[var]=0.;
      }
    }

    int numReadsImplicitSources, numReadsTimeDerivSources;
    int numWritesImplicitSources, numWritesTimeDerivSources;
    elemOld->computeImplicitSources(*geomCenter, *sourcesImplicitOld,
				                            elemHalfStep->tau,
                                    numReadsImplicitSources,
                                    numWritesImplicitSources
                                   );
    elem->computeImplicitSources(*geomCenter, *sourcesImplicit,
				                         elemHalfStep->tau,
                                 numReadsImplicitSources,
                                 numWritesImplicitSources
                                );
    elemHalfStep->computeTimeDerivSources(*geomCenter,
	                                			  *elemOld, *elem,
					                                dt,
					                                *sourcesTimeDer,
                                          numReadsTimeDerivSources,
                                          numWritesTimeDerivSources
                                         );
    numReads  += 2*numReadsImplicitSources  + numReadsTimeDerivSources;
    numWrites += 2*numWritesImplicitSources + numWritesTimeDerivSources; 

    for (int var=0; var<numVars; var++)
    {
      residualGuess.vars[var] = 
        (cons->vars[var] - consOld->vars[var])/dt
    	+ divFluxes->vars[var]
	    + sourcesExplicit->vars[var]
	    + 0.5*(sourcesImplicitOld->vars[var] + sourcesImplicit->vars[var])
	    + sourcesTimeDer->vars[var];
    }
    /* Reads:
     * -----
     *  cons[var], consOld[var], divFluxes[var]       : 3*numVars
     *  sourcesExplicit[var], sourcesTimeDer[var]     : 2*numVars
     *  sourcesImplicitOld[var], sourcesImplicit[var] : 2*numVars
     *
     * Writes:
     * ------
     * residualGuess[var] : numVars */
    numReads  += 7*residualGuess.numVars;

    /* Normalization of the residualGuess */
    if (params::conduction)
    {
	    if(params::highOrderTermsConduction)
      {
	      residualGuess.vars[vars::Q] *= 
          elemHalfStep->temperature
        * af::sqrt(elemHalfStep->rho*elemHalfStep->chi_emhd*elemHalfStep->tau);

        numReads += 4;
      }
    	else
      {
	      residualGuess.vars[vars::Q] *= elemHalfStep->tau;
        numReads += 1;
      }
    }

    if (params::viscosity)
    {
      if (params::highOrderTermsViscosity)
      {
	      residualGuess.vars[vars::DP] *= 
          af::sqrt(   elemHalfStep->rho*elemHalfStep->nu_emhd
                    * elemHalfStep->temperature*elemHalfStep->tau
                  );

        numReads += 4;
      }
	    else
      {
	      residualGuess.vars[vars::DP] *= elemHalfStep->tau;
        numReads += 1;
      }
    }

  } /* End of timeStepperSwitches::FULL_STEP */

  //Zero the residual in global ghost zones
  // TODO: prefill a mask with zero in the ghost zones and one in the bulk and
  // and then do residual = residual*mask;
  int N1Total = residualGuess.N1Total;
  int N2Total = residualGuess.N2Total;
  int N3Total = residualGuess.N3Total;

  int N1Local = residualGuess.N1Local;
  int N2Local = residualGuess.N2Local;
  int N3Local = residualGuess.N3Local;

  int numGhostX1 = residualGuess.numGhostX1;
  int numGhostX2 = residualGuess.numGhostX2;
  int numGhostX3 = residualGuess.numGhostX3;

  for (int var=0; var<residualGuess.numVars; var++) 
  {
    /* Left boundary */
    for(int i=0; i<numGhostX1; i++)
    {
	    residualGuess.vars[var](i,span,span)=0.;
    }

    /* Right boundary */
    for(int i=N1Local+numGhostX1; i<N1Local + 2*numGhostX1; i++)
    {
	    residualGuess.vars[var](i,span,span)=0.;
    }
    
    if(params::dim==1)
	    continue;
      
    /* Bottom boundary */
    for(int j=0; j < numGhostX2; j++)
    {
	    residualGuess.vars[var](span,j,span)=0.;
    }
      
    /* Top boundary */
    for(int j=N2Local+numGhostX2; j < N2Local+2*numGhostX2; j++)
    {
	    residualGuess.vars[var](span,j,span)=0.;
    }
      
    if(params::dim==2)
	    continue;
      
    /* Back boundary */
    for(int k=0; k<numGhostX3; k++)
    {
	    residualGuess.vars[var](span,span,k)=0.;
    }
      
    /* Front boundary */
    for(int k=N3Local+numGhostX3; k<N3Local + 2*numGhostX3; k++)
    {
	    residualGuess.vars[var](span,span,k)=0.;
    }
  }

  for (int var=0; var < residualGuess.numVars; var++)
  {
    residualGuess.vars[var].eval();
  }
  numWrites += residualGuess.numVars;
}
