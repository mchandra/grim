#include "../../params.hpp"
#include <arrayfire.h>

namespace params
{
  int N1 = 128;
  int N2 = 128;
  int N3 = 1;
  int dim = 2;
  int numGhost = 3;

  int timeStepper = timeStepping::EXPLICIT;
  double dt = .002;
  double CourantFactor = 0.9;
  double Time = 0.;
  double finalTime = 50.;
  int metric = metrics::MINKOWSKI;
  double hSlope = 0.3;

  int ObserveEveryNSteps = 10;
  int StepNumber = 0;

  double adiabaticIndex = 5./3;
  double blackHoleSpin = 0.9375;
  double InnerEdgeRadius = 6.;
  double PressureMaxRadius = 12.;
  double MinPlasmaBeta  = 15.;
  double MagneticLoops = 1;
  double Adiabat = 0.001;

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
  
  //Atmosphere parameters
  // Floors are Ampl*pow(radius,power)
  double RhoFloorAmpl = 1.e-3;
  double UFloorAmpl = 1.e-5;
  double RhoFloorSlope = -1.5;
  double UFloorSlope = -2.5;
  // Floors for magnetically dominated regions
  double BsqrOverRhoMax = 10.;
  double BsqrOverUMax = 500.;

  int conduction = 0;
  int viscosity  = 0;
  int highOrderTermsConduction = 1.;
  int highOrderTermsViscosity = 1.;
  double ConductionAlpha = 1.;
  double ViscosityAlpha = 1.;
  double ConductionClosureFactor = 1.;
  double ViscosityClosureFactor = 1.;

  double slopeLimTheta = 2.;
  int reconstruction = reconstructionOptions::WENO5;
  int riemannSolver  = riemannSolvers::LOCAL_LAX_FRIEDRICH;

  int maxNonLinearIter = 3;
  int maxLineSearchIters = 3;

  //Parameters controlling accuracy of nonlinear solver
  double nonlinearsolve_atol = 1.e-12;
  double JacobianAssembleEpsilon = 4.e-8;
  double linesearchfloor = 1.e-24;
  
  double InitialPerturbationAmplitude = 4e-2;
  double ObserveEveryDt = 1.;
  double WriteDataEveryDt = 100.;
};

namespace vars
{
  int dof = 5+params::conduction+params::viscosity;
  int Q = 5;
  int DP = 6;
};

