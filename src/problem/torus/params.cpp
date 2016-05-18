#include "torus.hpp"

namespace params
{
  // 4 GPUs on SAVIO
  int numDevices = 4;

  // Grid size options
  int N1 = 140;
  int N2 = 128;
  int N3 = 16;
  int dim = 3;
  int numGhost = 3;

  // (Re)start options
  double Time = 0.0;
  double InitialDt = 0.01;
  int restart = 0;
  std::string restartFile = "restartFile.h5";

  // Observation / checkpointing intervals
  double ObserveEveryDt = .1;
  double WriteDataEveryDt = 10.;

  // Timestepper opts
  int timeStepper = timeStepping::EXPLICIT;
  double CourantFactor = 0.4;
  double finalTime = 2000.;
  int metric = metrics::MODIFIED_KERR_SCHILD;

  // Grid size options
  double Rin = 0.98*(1.+sqrt(1.-blackHoleSpin*blackHoleSpin));
  double Rout = 63.;

  // Initial conditions
  double adiabaticIndex = 5./3;
  double blackHoleSpin = 0.9375;
  double InnerEdgeRadius = 6.;
  double PressureMaxRadius = 12.;
  double MinPlasmaBeta  = 15.;
  double MagneticLoops = 1;
  double Adiabat = 0.001;
  double InitialPerturbationAmplitude = 4.e-2;

  // EMHD model
  int conduction = 1;
  int viscosity  = 1;
  int highOrderTermsConduction = 1.;
  int highOrderTermsViscosity = 1.;
  // Phi and Psi from Chandra et al. 2015
  // Note that when using Phi=Psi, approximate
  // limits for the characteristic speeds to remain
  // subluminal are 0.29 (Gamma=5/3) and 1.3 (Gamma=4/3).
  double ConductionAlpha = 0.25;
  double ViscosityAlpha = 0.25;
  double ConductionClosureFactor = 1.;
  double ViscosityClosureFactor = 1.;

  //Atmosphere parameters
  double MaxLorentzFactor = 10.;
  // Floors are Ampl*pow(radius,power)
  double RhoFloorAmpl = 1.e-3;
  double UFloorAmpl = 1.e-5;
  double RhoFloorSlope = -1.5;
  double UFloorSlope = -2.5;
  // Floors for magnetically dominated regions
  double BsqrOverRhoMax = 10.;
  double BsqrOverUMax = 500.;

  // Automatic grid boundaries - do not change
  double X1Start = log(Rin), X1End = log(Rout);
  double X2Start = 0.+1.e-8, X2End = 1.-1.e-8;
  double X3Start = 0., X3End = M_PI*1.;

  // Grid parameters
  double hSlope = 0.3;
  int DerefineThetaHorizon = 1;
  int DoCylindrify = 1;
  double X1cyl = log(7.*Rin);
  double X2cyl = 3./N2;

  // Boundary Conditions
  // Radial
  int boundaryLeft   = boundaries::OUTFLOW;
  int boundaryRight  = boundaries::OUTFLOW;
  // Theta
  int boundaryTop    = boundaries::OUTFLOW;
  int boundaryBottom = boundaries::OUTFLOW;
  // Phi
  int boundaryFront  = boundaries::PERIODIC;
  int boundaryBack   = boundaries::PERIODIC;

  // Basic thresholds in fluid element
  double rhoFloorInFluidElement         = 1e-20;
  double uFloorInFluidElement           = 1e-20;
  double bSqrFloorInFluidElement        = 1e-20;
  double temperatureFloorInFluidElement = 1e-20;
  
  // Reconstruction options
  double slopeLimTheta = 2;
  int reconstruction = reconstructionOptions::MINMOD;
  int riemannSolver  = riemannSolvers::LOCAL_LAX_FRIEDRICH;


  //Parameters controlling accuracy of nonlinear solver
  int maxNonLinearIter = 3;
  int maxLineSearchIters = 3;
  double nonlinearsolve_atol = 1.e-6;
  double JacobianAssembleEpsilon = 4.e-8;
  double linesearchfloor = 1.e-24;
  
};

namespace vars
{
  int dof = 8+params::conduction+params::viscosity;
  int Q = 8;
  int DP = 8+params::conduction;
};

