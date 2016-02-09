#include "../problem.hpp"

namespace vars
{
  int Q = 8;
  int DP = 9;
  int dof = 10;
};

namespace params
{
  int N1 = 40;
  int N2 = 40;
  int N3 = 1;

  int dim = 2;
  int numGhost = 3;

  int timeStepper = timeStepping::EXPLICIT;
  double dt = .5/N1;
  double Time = 0.;
  double finalTime = 0.5;
  int metric = metrics::MINKOWSKI;

  double X1Start = 0., X1End = 1.;
  double X2Start = 0., X2End = 1.;
  double X3Start = 0., X3End = 1.;

  int boundaryLeft   = boundaries::PERIODIC;
  int boundaryRight  = boundaries::PERIODIC;

  int boundaryTop    = boundaries::PERIODIC;
  int boundaryBottom = boundaries::PERIODIC;

  int boundaryFront  = boundaries::PERIODIC;
  int boundaryBack   = boundaries::PERIODIC;

  double rhoFloorInFluidElement         = 1e-20;
  double uFloorInFluidElement           = 1e-20;
  double bSqrFloorInFluidElement        = 1e-20;
  double temperatureFloorInFluidElement = 1e-20;

  int conduction = 1;
  int viscosity  = 1;
  int highOrderTermsConduction = 0;
  int highOrderTermsViscosity = 0;

  double adiabaticIndex = 4./3;
  double Aw = 1.e-5;
  double k1 = 2.*M_PI;
  double k2 = 4.*M_PI;
  double Gamma = - 0.5533585207638141;
  double Omega = - 3.6262571286888425;

  double slopeLimTheta = 2;
  int reconstruction = reconstructionOptions::WENO5;
  int riemannSolver  = riemannSolvers::HLL;

  int maxNonLinearIter = 3;
  int maxLineSearchIters = 3;

  //Parameters controlling accuracy of nonlinear solver
  double nonlinearsolve_atol = 1.e-10;
  double JacobianAssembleEpsilon = 4.e-8;
  double linesearchfloor = 1.e-24;

  //Unused params - do we need to define them?
  double hSlope = 0.3;
  double blackHoleSpin = 0.9375;

  
};

void fluidElement::setFluidElementParameters(const geometry &geom)
{
  tau = one;
  chi_emhd = soundSpeed*soundSpeed*tau;
  nu_emhd  = soundSpeed*soundSpeed*tau;
}

void timeStepper::initialConditions(int &numReads,int &numWrites)
{
  array xCoords[3];
  for(int d=0;d<3;d++)
    xCoords[d]=XCoords->vars[d];
  geomCenter->XCoordsToxCoords(XCoords->vars,xCoords);

  array cphi = af::cos(  params::k1*xCoords[directions::X1]
                  		 + params::k2*xCoords[directions::X2]
                      );

  array sphi = af::sin(  params::k1*xCoords[directions::X1]
                  		 + params::k2*xCoords[directions::X2]
                      );

  /* Initial conditions */

  //Alfven Wave
  /*primOld->vars[vars::RHO] = 1.;
  primOld->vars[vars::U]   = 2.;
  primOld->vars[vars::U1]  = 0.;
  primOld->vars[vars::U2]  = params::Aw*0.462905090215*cphi;
  primOld->vars[vars::U3]  = 0.;
  primOld->vars[vars::B1]  = 0.01;
  primOld->vars[vars::B2]  = params::Aw*0.886407850514*cphi;
  primOld->vars[vars::B3]  = 0.;*/

  //Sound wave
  /*primOld->vars[vars::RHO] = 1.+params::Aw*0.345991032308*cphi;
  primOld->vars[vars::U]   = 2.+params::Aw*0.922642752822*cphi;
  primOld->vars[vars::U1]  = 0.-params::Aw*0.170354208129*cphi; 
  primOld->vars[vars::U2]  = 0.; 
  primOld->vars[vars::U3]  = 0.; 
  primOld->vars[vars::B1]  = 0.01; 
  primOld->vars[vars::B2]  = 0.;
  primOld->vars[vars::B3]  = 0.;*/

  //Full EMHD mode (from grim2D)
  primOld->vars[vars::RHO] = 1.;
  primOld->vars[vars::U]   = 2.;
  primOld->vars[vars::U1]  = 0.;
  primOld->vars[vars::U2]  = 0.; 
  primOld->vars[vars::U3]  = 0.; 
  primOld->vars[vars::B1]  = 0.1; 
  primOld->vars[vars::B2]  = 0.3;
  primOld->vars[vars::B3]  = 0.;
  primOld->vars[vars::Q]   = 0.;
  primOld->vars[vars::DP]  = 0.;

  primOld->vars[vars::RHO] += params::Aw*cphi*(-0.518522524082246)
                    +params::Aw*sphi*0.1792647678001878;

  primOld->vars[vars::U]   += params::Aw*cphi*0.5516170736393813;

  primOld->vars[vars::U1]  += params::Aw*cphi*0.008463122479547856
                    +params::Aw*sphi*(-0.011862022608466367);

  primOld->vars[vars::U2]  += params::Aw*cphi*(-0.16175466371870734)
                    +params::Aw*sphi*(0.034828080823603294);

  primOld->vars[vars::B1]  += params::Aw*cphi*(-0.05973794979640743)
                    +params::Aw*sphi*0.03351707506150924;

  primOld->vars[vars::B2]  += params::Aw*cphi*0.02986897489820372
                    -params::Aw*sphi*0.016758537530754618;

  primOld->vars[vars::Q]   += params::Aw*cphi*0.5233486841539436
                    -params::Aw*sphi*0.04767672501939603;

  primOld->vars[vars::DP]  += params::Aw*cphi*0.2909106062057657
                    -params::Aw*sphi*0.02159452055336572;

  for (int var=0; var < numVars; var++)
  {
    primOld->vars[var].eval();
  }
  af::sync();
}

void timeStepper::halfStepDiagnostics(int &numReads,int &numWrites)
{

}

void timeStepper::fullStepDiagnostics(int &numReads,int &numWrites)
{
  /* Compute the errors for the different modes */

  array xCoords[3];
  for(int d=0;d<3;d++)
    xCoords[d]=XCoords->vars[d];
  geomCenter->XCoordsToxCoords(XCoords->vars,xCoords);

  //Alfven wave
  /*cphi = af::cos(2*M_PI*ts.geom->xCoords[locations::CENTER][1]
    +0.0328124176673*params::Time);
    array u2an = params::Aw*0.462905090215*cphi;
    double error = af::norm(af::flat((ts.primOld->vars[vars::U2]
    - u2an)));*/
  
  //MHD Sound wave
  /*cphi = af::cos(2*M_PI*ts.geom->xCoords[locations::CENTER][1]
    +3.09362659024*params::Time);
    array rhoan = 1.+params::Aw*0.345991032308*cphi;
    double error = af::norm(af::flat((ts.primOld->vars[vars::RHO]
    - rhoan)));*/
  
  //EMHD Sound wave
  array cphi = af::cos(  params::k1*xCoords[directions::X1]
			 + params::k2*xCoords[directions::X2]
			 + params::Omega*time
			 );
  
  array sphi = af::sin(  params::k1*xCoords[directions::X1]
			 + params::k2*xCoords[directions::X2]
			 + params::Omega*time
			 );
  
  array rhoAnalytic = 1. + (  params::Aw*cphi*(-0.518522524082246)
			      + params::Aw*sphi*0.1792647678001878
			      )*exp(params::Gamma*time);
  
  array dPAnalytic = 0. + ( params::Aw*cphi*0.2909106062057657
			    -params::Aw*sphi*0.02159452055336572
			    )*exp(params::Gamma*time);
  
  double errorRho = af::norm(af::flat(primOld->vars[vars::RHO] - rhoAnalytic));
  
  double errordP  = af::norm(af::flat(primOld->vars[vars::DP]  - dPAnalytic));
  
  //af_print(primOld->vars[vars::RHO] - rhoAnalytic,12);
  
  errorRho = errorRho/N1/N2/N3;
  errordP  = errordP/N1/N2/N3;
  PetscPrintf(PETSC_COMM_WORLD, 
              "Time = %e; dt = %e; Error in rho = %e; Error in dP = %e\n",
              time, dt, errorRho, errordP
              );
}

void timeStepper::setProblemSpecificBCs(int &numReads,int &numWrites)
{

};
