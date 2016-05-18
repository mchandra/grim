#include "timestepper.hpp"

void timeStepper::computeResidual(const grid &primGuess,
                                  grid &residualGuess, 
                                  const bool computeExplicitTerms,
                                  int &numReads,
                                  int &numWrites,
				  const array Idx
                                 )
{
  numReads = 0; numWrites = 0;
  int numReadsElemSet, numWritesElemSet;
  int numReadsComputeFluxes, numWritesComputeFluxes;
  elem->set(primGuess, *geomCenter, numReadsElemSet, numWritesElemSet,Idx);
  elem->computeFluxes(*geomCenter, 0, *cons, 
                      numReadsComputeFluxes, numWritesComputeFluxes,
		      Idx
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
                                      numWritesExplicitSouces,
				      Idx
                                     );
      numReads  += numReadsExplicitSouces;
      numWrites += numWritesExplicitSouces;
    }
    else
    {
      for (int var=0; var<vars::dof; var++)
      {
	sourcesExplicit->vars[var](Idx) = 0.;
      }
    }

    int numReadsImplicitSources, numReadsTimeDerivSources;
    int numWritesImplicitSources, numWritesTimeDerivSources;
    elemOld->computeImplicitSources(*geomCenter, *sourcesImplicitOld,
				                            elemOld->tau,
                                    numReadsImplicitSources,
                                    numWritesImplicitSources,
				    Idx
                                   );
    elem->computeImplicitSources(*geomCenter, *sourcesImplicit,
				                         elemOld->tau,
                                 numReadsImplicitSources,
                                 numWritesImplicitSources,
				 Idx
                                );
    elemOld->computeTimeDerivSources(*geomCenter,
				     *elemOld, *elem,
				     dt/2,
				     *sourcesTimeDer,
                                     numReadsTimeDerivSources,
                                     numWritesTimeDerivSources,
				     Idx
                                    );
    numReads  += 2*numReadsImplicitSources  + numReadsTimeDerivSources;
    numWrites += 2*numWritesImplicitSources + numWritesTimeDerivSources; 

    for (int var=0; var<residualGuess.numVars; var++)
    {
      residualGuess.vars[var](Idx) = 
        (cons->vars[var](Idx) - consOld->vars[var](Idx))/(dt/2.)
	+ divFluxes->vars[var](Idx)
	+ sourcesExplicit->vars[var](Idx)
	+ 0.5*(sourcesImplicitOld->vars[var](Idx) + sourcesImplicit->vars[var](Idx))
	+ sourcesTimeDer->vars[var](Idx);
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
	residualGuess.vars[vars::Q](Idx) *=
	  elemOld->temperature(Idx) 
	  * af::sqrt(elemOld->rho(Idx)*elemOld->chi_emhd(Idx)*elemOld->tau(Idx));

        numReads += 4;
      }
      else
      {
	residualGuess.vars[vars::Q](Idx) *= elemOld->tau(Idx);
        numReads += 1;
      }
    }

    if (params::viscosity)
    {
      if(params::highOrderTermsViscosity)
      {
	residualGuess.vars[vars::DP](Idx) *=
          af::sqrt(   elemOld->rho(Idx)*elemOld->nu_emhd(Idx)
		      * elemOld->temperature(Idx)*elemOld->tau(Idx)
		      );
        numReads += 4;
      }
      else
      {
	residualGuess.vars[vars::DP](Idx) *= elemOld->tau(Idx);
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
                                           numWritesExplicitSouces,
					   Idx
                                          );
      numReads  += numReadsExplicitSouces;
      numWrites += numWritesExplicitSouces;
    }
    else
    {
      for (int var=0; var<vars::dof; var++)
      {
	sourcesExplicit->vars[var](Idx)=0.;
      }
    }

    int numReadsImplicitSources, numReadsTimeDerivSources;
    int numWritesImplicitSources, numWritesTimeDerivSources;
    elemOld->computeImplicitSources(*geomCenter, *sourcesImplicitOld,
				    elemHalfStep->tau,
                                    numReadsImplicitSources,
                                    numWritesImplicitSources,
				    Idx
                                   );
    elem->computeImplicitSources(*geomCenter, *sourcesImplicit,
				 elemHalfStep->tau,
                                 numReadsImplicitSources,
                                 numWritesImplicitSources,
				 Idx
                                );
    elemHalfStep->computeTimeDerivSources(*geomCenter,
					  *elemOld, *elem,
					  dt,
					  *sourcesTimeDer,
                                          numReadsTimeDerivSources,
                                          numWritesTimeDerivSources,
					  Idx
                                         );
    numReads  += 2*numReadsImplicitSources  + numReadsTimeDerivSources;
    numWrites += 2*numWritesImplicitSources + numWritesTimeDerivSources; 

    for (int var=0; var<vars::dof; var++)
    {
      residualGuess.vars[var](Idx) = 
        (cons->vars[var](Idx) - consOld->vars[var](Idx))/dt
    	+ divFluxes->vars[var](Idx)
	+ sourcesExplicit->vars[var](Idx)
	+ 0.5*(sourcesImplicitOld->vars[var](Idx) + sourcesImplicit->vars[var](Idx))
	+ sourcesTimeDer->vars[var](Idx);
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
	residualGuess.vars[vars::Q](Idx) *= 
          elemHalfStep->temperature(Idx)
	  * af::sqrt(elemHalfStep->rho(Idx)*elemHalfStep->chi_emhd(Idx)*elemHalfStep->tau(Idx));

        numReads += 4;
      }
    	else
      {
	residualGuess.vars[vars::Q](Idx) *= elemHalfStep->tau(Idx);
        numReads += 1;
      }
    }

    if (params::viscosity)
    {
      if (params::highOrderTermsViscosity)
      {
	residualGuess.vars[vars::DP](Idx) *= 
          af::sqrt(   elemHalfStep->rho(Idx)*elemHalfStep->nu_emhd(Idx)
		      * elemHalfStep->temperature(Idx)*elemHalfStep->tau(Idx)
                  );

        numReads += 4;
      }
	    else
      {
	residualGuess.vars[vars::DP](Idx) *= elemHalfStep->tau(Idx);
        numReads += 1;
      }
    }

  } /* End of timeStepperSwitches::FULL_STEP */

  //Zero the residual in global ghost zones
  for (int var=0; var < residualGuess.numVars; var++)
  {
    residualGuess.vars[var] *= GZmask;
    residualGuess.vars[var].eval();
  }
  numWrites += residualGuess.numVars;
}
