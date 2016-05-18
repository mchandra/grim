#ifndef GRIM_PARAMS_H_
#define GRIM_PARAMS_H_
#include <string>

const int NDIM = 4;
const int LOCATIONS = 7;

namespace vars
{
  enum
  {
    RHO, U, U1, U2, U3, B1, B2, B3
  };
  extern int Q, DP;
  extern int dof;
};

namespace locations
{
  enum
  {
    CENTER, LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK
  };
};

namespace directions
{
  enum
  {
    X1, X2, X3
  };
};

namespace boundaries
{
  enum
  {
    PERIODIC, OUTFLOW, MIRROR, DIRICHLET
  };
};

namespace metrics
{
  enum
  {
    MINKOWSKI, MODIFIED_KERR_SCHILD
  };
};

namespace timeStepping
{
  enum
  {
    EXPLICIT, IMEX, IMPLICIT
  };
};

namespace reconstructionOptions
{
  enum
  {
    MINMOD, WENO5
  };
};

namespace riemannSolvers
{
  enum
  {
    HLL, LOCAL_LAX_FRIEDRICH
  };
};

namespace params
{
  extern int numDevices;

  extern int N1;
  extern int N2;
  extern int N3;
  extern int dim;
  extern int numGhost;

  extern int timeStepper;
  extern double InitialDt;
  extern double Time;
  extern double finalTime;
  extern int metric;
  extern double hSlope;
  extern double blackHoleSpin;

  extern int restart;
  extern std::string restartFile;

  extern double X1Start, X1End;
  extern double X2Start, X2End;
  extern double X3Start, X3End;

  extern int DerefineThetaHorizon;
  extern int DoCylindrify;
  extern double X1cyl, X2cyl;

  extern int boundaryLeft;
  extern int boundaryRight;

  extern int boundaryTop;
  extern int boundaryBottom;

  extern int boundaryFront;
  extern int boundaryBack;

  extern double rhoFloorInFluidElement;
  extern double uFloorInFluidElement;
  extern double bSqrFloorInFluidElement;
  extern double temperatureFloorInFluidElement;

  extern int conduction;
  extern int viscosity;
  extern int highOrderTermsConduction;
  extern int highOrderTermsViscosity;
  extern double adiabaticIndex;
  extern double ConductionAlpha;
  extern double ViscosityAlpha;

  extern double slopeLimTheta;
  extern int reconstruction;
  extern int riemannSolver;

  extern int maxNonLinearIter;
  extern int maxLineSearchIters;

  extern double nonlinearsolve_atol;
  extern double JacobianAssembleEpsilon;
  extern double linesearchfloor;

  //Atmosphere parameters
  extern double MaxLorentzFactor;
  // Floors are Ampl*pow(radius,power)
  extern double RhoFloorAmpl;
  extern double UFloorAmpl;
  extern double RhoFloorSlope;
  extern double UFloorSlope;
  // Floors for magnetically dominated regions
  extern double BsqrOverRhoMax;
  extern double BsqrOverUMax;

  extern double ConductionClosureFactor;
  extern double ViscosityClosureFactor;

  extern int ObserveEveryNSteps;
  extern int StepNumber;

  extern double InnerEdgeRadius ;
  extern double PressureMaxRadius;
  extern double MinPlasmaBeta;
  extern double MagneticLoops;
  extern double Adiabat;
  extern double CourantFactor;
  extern double InitialPerturbationAmplitude;
  extern double ObserveEveryDt;
  extern double WriteDataEveryDt;

  /* Linear modes parameters */
  extern double Aw;
  extern double k1;
  extern double k2;
  extern double Gamma;
  extern double Omega;
};

#endif /* GRIM_PARAMS_H_ */
