#include "torus.hpp"

void fluidElement::setFluidElementParameters(const geometry &geom, const array Idx)
{
  array xCoords[3];
  geom.getxCoords(xCoords);
  array Radius = xCoords[0](Idx);
  array DynamicalTimescale = af::pow(Radius,1.5);
  DynamicalTimescale.eval();
  tau(Idx) = DynamicalTimescale;

  if(params::conduction)
    {
      array Qmax = params::ConductionClosureFactor*rho(Idx)*pow(soundSpeed(Idx),3.);
      double lambda = 0.01;
      array yCon = af::abs(q(Idx))/Qmax;
      yCon = af::exp(-(yCon-1.)/lambda);
      yCon.eval();
      array fdCon = yCon/(yCon+1.)+1.e-5;
      tau(Idx)=af::min(tau(Idx),DynamicalTimescale*fdCon);
    }
  if(params::viscosity)
    {
      array dPmod = af::max(pressure(Idx)-2./3.*deltaP(Idx),params::bSqrFloorInFluidElement)/af::max(pressure(Idx)+1./3.*deltaP(Idx),params::bSqrFloorInFluidElement);
      array dPmaxPlus = af::min(params::ViscosityClosureFactor*bSqr(Idx)*0.5*dPmod,1.49*pressure(Idx)/1.07);
      array dPmaxMinus = af::max(-params::ViscosityClosureFactor*bSqr(Idx),-2.99*pressure(Idx)/1.07);

      array condition = deltaP(Idx)>0.;
      array dPmax = condition*dPmaxPlus + (1.-condition)*dPmaxMinus;

      double lambda = 0.01;
      array yVis = af::abs(deltaP(Idx))/(af::abs(dPmax)+params::bSqrFloorInFluidElement);
      yVis = af::exp(-(yVis-1.)/lambda);
      yVis.eval();
      array fdVis = yVis/(yVis+1.)+1.e-5;
      tau(Idx)=af::min(tau(Idx),DynamicalTimescale*fdVis);
    }
  tau(Idx).eval();
  chi_emhd(Idx) = params::ConductionAlpha*soundSpeed(Idx)*soundSpeed(Idx)*tau(Idx);
  nu_emhd(Idx)  = params::ViscosityAlpha*soundSpeed(Idx)*soundSpeed(Idx)*tau(Idx);
  chi_emhd(Idx).eval();
  nu_emhd(Idx).eval();
  af::sync();
}


/* Calculate the constant angular momentum per unit inertial mass (l = u_phi *
 * u^t) for a given black hole spin and a radius of the accretion disk.  Eqn 3.8
 * of Fishbone and Moncrief, 1976 */
double lFishboneMoncrief(double a, double r, double theta)
{
  double M = 1.;
  return sqrt(M/pow(r, 3.)) \
        *(  pow(r, 4.) + r*r*a*a - 2.*M*r*a*a \
          - a*sqrt(M*r)*(r*r - a*a) \
         )/ \
         (r*r - 3*M*r + 2.*a*sqrt(M*r));
}

double lnOfhTerm1(double a,
                double r, double theta, 
                double l)
{
  double Delta = computeDelta(a, r, theta);
  double Sigma = computeSigma(a, r, theta);
  double A     = computeA(a, r, theta);

  return 0.5*log( (1. + sqrt(1. + (4.*l*l*Sigma*Sigma*Delta)/ \
                                  (A*sin(theta)*A*sin(theta))
                            )
                  ) / (Sigma*Delta/A)
                );
}

double lnOfhTerm2(double a,
                double r, double theta, 
                double l)
{
  double Delta = computeDelta(a, r, theta);
  double Sigma = computeSigma(a, r, theta);
  double A     = computeA(a, r, theta);

  return -0.5*sqrt(1. + (4.*l*l*Sigma*Sigma*Delta) /
                        (A*A*sin(theta)*sin(theta))
                  );

}

double lnOfhTerm3(double a,
                double r, double theta, 
                double l)
{
  double A     = computeA(a, r, theta);
  double M = 1.;
  return -2*a*M*r*l/A;
}

double computeDelta(double a, double r, double theta)
{
  double M = 1.;
  return r*r - 2*M*r + a*a;
}

double computeSigma(double a, double r, double theta)
{
  return r*r + a*a*cos(theta)*cos(theta);
}

double computeA(double a, double r, double theta)
{
  double Delta = computeDelta(a, r, theta);

  return pow(r*r + a*a, 2.) - Delta*a*a*sin(theta)*sin(theta);
}


double computeLnOfh(double a, double r, double theta)
{
  double l = lFishboneMoncrief(a, params::PressureMaxRadius, M_PI/2.);

  double term1 = lnOfhTerm1(a, r, theta, l);
  double term2 = lnOfhTerm2(a, r, theta, l);
  double term3 = lnOfhTerm3(a, r, theta, l);

  double term1InnerEdge = lnOfhTerm1(a, params::InnerEdgeRadius, M_PI/2., l);
  double term2InnerEdge = lnOfhTerm2(a, params::InnerEdgeRadius, M_PI/2., l);
  double term3InnerEdge = lnOfhTerm3(a, params::InnerEdgeRadius, M_PI/2., l);

  return  term1 + term2 + term3 \
        - term1InnerEdge - term2InnerEdge - term3InnerEdge;

}

void timeStepper::initialConditions(int &numReads,int &numWrites)
{
  int world_rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(PETSC_COMM_WORLD, &world_size);

  // Random number generator from PeTSC
  double randNum;
  PetscRandom randNumGen;
  PetscRandomCreate(PETSC_COMM_WORLD, &randNumGen);
  PetscRandomSetType(randNumGen, PETSCRAND48);


  //Let's ignore ArrayFire here - it's only the initial
  //conditions...
  array xCoords[3];
  geomCenter->getxCoords(xCoords);

  const int N1g = primOld->N1Total;
  const int N2g = primOld->N2Total;
  const int N3g = primOld->N3Total;

  array& Rho = primOld->vars[vars::RHO];
  array& U = primOld->vars[vars::U];
  array& U1 = primOld->vars[vars::U1];
  array& U2 = primOld->vars[vars::U2];
  array& U3 = primOld->vars[vars::U3];
  array& B1 = primOld->vars[vars::B1];
  array& B2 = primOld->vars[vars::B2];
  array& B3 = primOld->vars[vars::B3];


  PetscPrintf(PETSC_COMM_WORLD, "Running on %i procs\n", world_size);
  for(int proc=0;proc<world_size;proc++)
    {
      if(world_rank==proc)
	{
	  printf("Local size on proc %i : %i x %i x %i\n",proc,N1g,N2g,N3g);
	  //af_print(xCoords[directions::X1],5);
	}
      MPI_Barrier(PETSC_COMM_WORLD);
    }
  
  double aBH = params::blackHoleSpin;

  for(int k=0;k<N3g;k++)
    for(int j=0;j<N2g;j++)
      for(int i=0;i<N1g;i++)
	{
	  const int p = i+j*N1g+k*N2g;
	  const double& r = xCoords[directions::X1](i,j,k).scalar<double>();
	  const double& theta = xCoords[directions::X2](i,j,k).scalar<double>();
	  const double& phi = xCoords[directions::X3](i,j,k).scalar<double>();
	  const double& X2 = XCoords->vars[directions::X2](i,j,k).scalar<double>();

	  const double& lapse = geomCenter->alpha(i,j,k).scalar<double>();
	  const double& beta1 = geomCenter->gCon[0][1](i,j,k).scalar<double>();
	  const double& beta2 = geomCenter->gCon[0][2](i,j,k).scalar<double>();
	  const double& beta3 = geomCenter->gCon[0][3](i,j,k).scalar<double>();


	  double lnOfh = 1.;
	  if(r>=params::InnerEdgeRadius)
	    lnOfh = computeLnOfh(aBH,r,theta);
	  
	  /* Region outside the torus */
	  if(lnOfh<0. || r<params::InnerEdgeRadius)
	    {
	      Rho(i,j,k)=params::rhoFloorInFluidElement;
	      U(i,j,k)=params::uFloorInFluidElement;
	      U1(i,j,k)=0.;
	      U2(i,j,k)=0.;
	      U3(i,j,k)=0.;
	    }
	  else
	    {
	      double h = exp(lnOfh);
	      double Gamma = params::adiabaticIndex;
	      double Kappa = params::Adiabat;

	      /* Solve for rho using the definition of h = (rho + u + P)/rho where rho
	       * here is the rest mass energy density and P = C * rho^Gamma */
	      Rho(i,j,k) = pow((h-1)*(Gamma-1.)/(Kappa*Gamma), 
			       1./(Gamma-1.));
	      PetscRandomGetValue(randNumGen, &randNum);
	      U(i,j,k) =  Kappa * pow(Rho(i,j,k), Gamma)/(Gamma-1.)
		*(1. + params::InitialPerturbationAmplitude*(randNum-0.5));
	  
	      /* TODO: Should add random noise here */
	  
	      
	      /* Fishbone-Moncrief u_phi is given in the Boyer-Lindquist coordinates.
	       * Need to transform to (modified) Kerr-Schild */
	      double A = computeA(aBH, r, theta);
	      double Sigma = computeSigma(aBH, r, theta);
	      double Delta = computeDelta(aBH, r, theta);
	      double l = lFishboneMoncrief(aBH, params::PressureMaxRadius, M_PI/2.);
	      double expOfMinus2Chi = Sigma*Sigma*Delta/(A*A*sin(theta)*sin(theta)) ;
	      double uCovPhiBL = sqrt((-1. + sqrt(1. + 4*l*l*expOfMinus2Chi))/2.);
	      double uConPhiBL =   2.*aBH*r*sqrt(1. + uCovPhiBL*uCovPhiBL)
		/sqrt(A*Sigma*Delta)+ sqrt(Sigma/A)*uCovPhiBL/sin(theta);

	      double uConBL[NDIM];
	      uConBL[0] = 0.;
	      uConBL[1] = 0.;
	      uConBL[2] = 0.;
	      uConBL[3] = uConPhiBL;
	      
	      double gCovBL[NDIM][NDIM], gConBL[NDIM][NDIM];
	      double transformBLToMKS[NDIM][NDIM];
	      
	      for (int alpha=0; alpha<NDIM; alpha++)
		{
		  for (int beta=0; beta<NDIM; beta++)
		    {
		      gCovBL[alpha][beta] = 0.;
		      gConBL[alpha][beta] = 0.;
		      transformBLToMKS[alpha][beta] = 0.;
		    }
		}
	      
	      double mu = 1 + aBH*aBH*cos(theta)*cos(theta)/(r*r);
	      
	      gCovBL[0][0] = -(1. - 2./(r*mu));
	      gCovBL[0][3] = -2.*aBH*sin(theta)*sin(theta)/(r*mu);
	      gCovBL[3][0] = gCovBL[0][3];
	      gCovBL[1][1] = mu*r*r/Delta;
	      gCovBL[2][2] = r*r*mu;
	      gCovBL[3][3] = r*r*sin(theta)*sin(theta)*
		(1. + aBH*aBH/(r*r) + 2.*aBH*aBH*sin(theta)*sin(theta)/(r*r*r*mu));
	      
	      gConBL[0][0] = -1. -2.*(1 + aBH*aBH/(r*r))/(Delta*mu/r);
	      gConBL[0][3] = -2.*aBH/(r*Delta*mu);
	      gConBL[3][0] = gConBL[0][3];
	      gConBL[1][1] = Delta/(r*r*mu);
	      gConBL[2][2] = 1./(r*r*mu);
	      gConBL[3][3] = (1. - 2./(r*mu))/(sin(theta)*sin(theta)*Delta);
	      
	      transformBLToMKS[0][0] = 1.;
	      transformBLToMKS[1][1] = 1.;
	      transformBLToMKS[2][2] = 1.;
	      transformBLToMKS[3][3] = 1.;
	      transformBLToMKS[0][1] = 2.*r/Delta;
	      transformBLToMKS[3][1] = aBH/Delta; 
	      
	      /* Need to get uConBL[0] using u^mu u_mu = -1 */
	      double AA = gCovBL[0][0];
	      double BB = 2.*(gCovBL[0][1]*uConBL[1] +
			      gCovBL[0][2]*uConBL[2] +
			      gCovBL[0][3]*uConBL[3]
			      );
	      double CC = 1. + gCovBL[1][1]*uConBL[1]*uConBL[1] +
		gCovBL[2][2]*uConBL[2]*uConBL[2] +
		gCovBL[3][3]*uConBL[3]*uConBL[3] +
		2.*(gCovBL[1][2]*uConBL[1]*uConBL[2] +
		    gCovBL[1][3]*uConBL[1]*uConBL[3] +
		    gCovBL[2][3]*uConBL[2]*uConBL[3]);
	      
	      double discriminent = BB*BB - 4.*AA*CC;
	      uConBL[0] = -(BB + sqrt(discriminent))/(2.*AA);
	      
	      double uConKS[NDIM];
	      
	      for (int alpha=0; alpha<NDIM; alpha++)
		{
		  uConKS[alpha] = 0.;
		  
		  for (int beta=0; beta<NDIM; beta++)
		    {
		      uConKS[alpha] += transformBLToMKS[alpha][beta]*uConBL[beta];
		    }
		}
	      
	      /* Finally get the four-velocity in the X coordinates, which is modified
	       * Kerr-Schild */
	      double uConMKS[NDIM];
	      double rFactor = r;
	      double hFactor = M_PI + (1. - params::hSlope)*M_PI*cos(2.*M_PI*X2);
	      uConMKS[0] = uConKS[0];
	      uConMKS[1] = uConKS[1]/rFactor;
	      uConMKS[2] = uConKS[2]/hFactor;
	      uConMKS[3] = uConKS[3];
	      
	      U1(i,j,k) = uConMKS[1] + pow(lapse, 2.)*beta1*uConMKS[0];
	      U2(i,j,k) = uConMKS[2] + pow(lapse, 2.)*beta2*uConMKS[0];
	      U3(i,j,k) = uConMKS[3] + pow(lapse, 2.)*beta3*uConMKS[0];
	    }
	  B1(i,j,k) = 0.;
	  B2(i,j,k) = 0.;
	  B3(i,j,k) = 0.;
	}

  array rhoMax_af = af::max(af::max(af::max(Rho,2),1),0);
  double rhoMax = rhoMax_af.host<double>()[0];

  /* Communicate rhoMax to all processors */
  if (world_rank == 0) 
    {
      double temp; 
      for(int i=1;i<world_size;i++)
	{
	  MPI_Recv(&temp, 1, MPI_DOUBLE, i, i, PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
	  if(rhoMax < temp)
	    rhoMax = temp;
	}
      }
    else
      {
        MPI_Send(&rhoMax, 1, MPI_DOUBLE, 0, world_rank, PETSC_COMM_WORLD);
      }
  MPI_Barrier(PETSC_COMM_WORLD);
  MPI_Bcast(&rhoMax,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
  MPI_Barrier(PETSC_COMM_WORLD);

  PetscPrintf(PETSC_COMM_WORLD,"rhoMax = %e\n",rhoMax);
  Rho=Rho/rhoMax;
  U=U/rhoMax;

  Rho.eval();
  U.eval();
  U1.eval();
  U2.eval();
  U3.eval();
  B1.eval();
  B2.eval();
  B3.eval();

  /* TODO : apply floors */

  /* Set magnetic field. This is MUCH easier to do using ArrayFire... */
  // Set vector potential
  const array& Rho_af = primOld->vars[vars::RHO];
  array rhoAvg = 
    (af::shift(Rho_af,1,0,0)+af::shift(Rho_af,-1,0,0)
     +af::shift(Rho_af,0,1,0)+af::shift(Rho_af,0,-1,0));
  if(params::dim>2)
    {
      rhoAvg = (rhoAvg
		+af::shift(Rho_af,0,0,1)+af::shift(Rho_af,0,0,-1))/6.;
    }
  else
    {
      rhoAvg = rhoAvg/4.;
    }
  array zero = rhoAvg*0.;
  array Avec = af::max(rhoAvg-0.2,zero)*
    af::cos(xCoords[directions::X2]*
	    (params::MagneticLoops-1));
  Avec.eval();
 
  // Compute magnetic field 
  const array& g = geomCenter->g;
  const double dX1 = XCoords->dX1;
  const double dX2 = XCoords->dX2;

  primOld->vars[vars::B1] = 
    (af::shift(Avec,0,-1,0)-af::shift(Avec,0,0,0)
     +af::shift(Avec,-1,-1,0)-af::shift(Avec,-1,0,0))/
    (2.*dX2*g);
  primOld->vars[vars::B2] = 
    (af::shift(Avec,0,0,0)-af::shift(Avec,-1,0,0)
     +af::shift(Avec,0,-1,0)-af::shift(Avec,-1,-1,0))/
    (2.*dX1*g);
  primOld->vars[vars::B1].eval();
  primOld->vars[vars::B2].eval();

  // Set fields to zero in ghost zones
  // (Communication below only sets fields to correct value in inner GZ,
  //  not on inner / outer boundary)
  for(int i=0;i<params::numGhost;i++)
    {
      primOld->vars[vars::B1](i,span,span)=0.;
      primOld->vars[vars::B2](i,span,span)=0.;
      primOld->vars[vars::B1](N1g-1-i,span,span)=0.;
      primOld->vars[vars::B2](N1g-1-i,span,span)=0.;
      if(params::dim>1)
	{
	  primOld->vars[vars::B1](span,i,span)=0.;
	  primOld->vars[vars::B2](span,i,span)=0.;
	  primOld->vars[vars::B1](span,N2g-1-i,span)=0.;
          primOld->vars[vars::B2](span,N2g-1-i,span)=0.;
	}
      if(params::dim>2)
	{
	  primOld->vars[vars::B1](span,span,i)=0.;
	  primOld->vars[vars::B2](span,span,i)=0.;
	  primOld->vars[vars::B1](span,span,N3g-1-i)=0.;
          primOld->vars[vars::B2](span,span,N3g-1-i)=0.;
	}
   }
  


  // We have used ghost zones in the previous steps -> need to communicate
  // before computing global quantities
  primOld->communicate();

  // Need to set fluid element to get b^2...
  {
    elemOld->set(*primOld,*geomCenter,numReads,numWrites);
    const array& bSqr = elemOld->bSqr;
    const array& Pgas = elemOld->pressure;
    array PlasmaBeta = 2.*(Pgas+1.e-13)/(bSqr+1.e-18);
    array BetaMin_af = af::min(af::min(af::min(PlasmaBeta,2),1),0);
    double BFactor = BetaMin_af.host<double>()[0];
    BFactor = sqrt(BFactor/params::MinPlasmaBeta);

    /* Use MPI to find minimum over all processors */
    if (world_rank == 0) 
      {
	double temp; 
	for(int i=1;i<world_size;i++)
	  {
	    MPI_Recv(&temp, 1, MPI_DOUBLE, i, i, PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
	    if(BFactor > temp)
	      BFactor = temp;
	  }
      }
    else
      {
        MPI_Send(&BFactor, 1, MPI_DOUBLE, 0, world_rank, PETSC_COMM_WORLD);
      }
    MPI_Barrier(PETSC_COMM_WORLD);
    MPI_Bcast(&BFactor,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
    MPI_Barrier(PETSC_COMM_WORLD);
    
    PetscPrintf(PETSC_COMM_WORLD,"Bfactor = %e\n",BFactor);

    primOld->vars[vars::B1] *= BFactor;
    primOld->vars[vars::B2] *= BFactor;
    primOld->vars[vars::B3] *= BFactor;
    primOld->vars[vars::B1].eval();
    primOld->vars[vars::B2].eval();
    primOld->vars[vars::B3].eval();
  }

  if(params::conduction)
    {
      primOld->vars[vars::Q] = zero;
      primOld->vars[vars::Q].eval();
    }
  if(params::viscosity)
    {
      primOld->vars[vars::DP] = zero;
      primOld->vars[vars::DP].eval();
    }

  applyFloor(primOld,elemOld,geomCenter,numReads,numWrites);

  for (int var=0; var<vars::dof; var++) 
    {
        primOld->vars[var].eval();
    }

  af::sync();
  fullStepDiagnostics(numReads,numWrites);
}

void applyFloor(grid* prim, fluidElement* elem, geometry* geom, int &numReads,int &numWrites)
{
  array xCoords[3];
  geom->getxCoords(xCoords);
  array Radius = xCoords[0];

  array minRho = af::pow(Radius,params::RhoFloorSlope)*params::RhoFloorAmpl;
  array minU = af::pow(Radius,params::UFloorSlope)*params::UFloorAmpl;

  // Save pre-floor values for later
  array rho_prefloor = prim->vars[vars::RHO]*1.;
  array u_prefloor   = prim->vars[vars::U]*1.;
  rho_prefloor.eval();
  u_prefloor.eval();

  // Apply floors when needed
  array condition = prim->vars[vars::RHO]<minRho;
  array usefloor  = condition;
  prim->vars[vars::RHO] = condition*minRho+(1.-condition)*prim->vars[vars::RHO];

  condition = prim->vars[vars::U]<minU;
  usefloor = af::max(condition,usefloor);
  prim->vars[vars::U] = condition*minU+(1.-condition)*prim->vars[vars::U];

  elem->set(*prim,*geom,numReads,numWrites);
  const array& bSqr = elem->bSqr;
  condition = bSqr>params::BsqrOverRhoMax*prim->vars[vars::RHO];
  usefloor = af::max(condition,usefloor);
  prim->vars[vars::RHO] = prim->vars[vars::RHO]*(1.-condition)+condition*bSqr/params::BsqrOverRhoMax;
  condition = bSqr>params::BsqrOverUMax*prim->vars[vars::U];
  usefloor = af::max(condition,usefloor);
  prim->vars[vars::U] = prim->vars[vars::U]*(1.-condition)+condition*bSqr/params::BsqrOverUMax;
  
  prim->vars[vars::RHO].eval();
  prim->vars[vars::U].eval();
  usefloor.eval();

  array zero = usefloor*0.;
  array trans = af::max(af::min(zero+1.,(bSqr-0.1*prim->vars[vars::RHO])/prim->vars[vars::RHO]),zero)
    *(usefloor>zero);
  trans.eval();
  
  // Sasha's floor method
  array betapar = -elem->bCon[0]/bSqr/elem->uCon[0];
  array betasqr = betapar*betapar*bSqr;
  double betasqrmax = 1. - 1./params::MaxLorentzFactor/params::MaxLorentzFactor;
  betasqr = af::min(betasqr,betasqrmax);
  array gamma = 1./af::sqrt(1.-betasqr);
  gamma.eval();
  array ucondr[NDIM];
  for(int m=0;m<NDIM;m++)
    {
      ucondr[m]=gamma*(elem->uCon[m]+betapar*elem->bCon[m]);
      ucondr[m].eval();
    }
  // Need b-field in the inertial frame
  array Bcon[NDIM];
  Bcon[0]= zero;
  Bcon[1] = prim->vars[vars::B1];
  Bcon[2] = prim->vars[vars::B2];
  Bcon[3] = prim->vars[vars::B3];
  array Bcov[NDIM];
  for(int m=0;m<NDIM;m++)
    {
      Bcov[m]=zero;
      for(int n=0;n<NDIM;n++)
	Bcov[m]+=geom->gCov[n][m]*Bcon[n];
      Bcov[m].eval();
    }
  array udotB = zero;
  array bSqr_fl  = zero;
  for(int n=0;n<NDIM;n++)
    {
      udotB += Bcov[n]*elem->uCon[n];
      bSqr_fl += Bcov[n]*Bcon[n];
    }
  udotB.eval();
  array bNorm = af::sqrt(bSqr_fl);
  double bMin = sqrt(params::bSqrFloorInFluidElement);
  bNorm = af::max(bNorm,bMin);
  bNorm.eval();

  //New velocity
  array wold = rho_prefloor+u_prefloor*params::adiabaticIndex;
  array QdotB = udotB*wold*elem->uCon[0];
  array wnew = prim->vars[vars::RHO]+prim->vars[vars::U]*params::adiabaticIndex;
  array x = 2.*QdotB/(bNorm*wnew*ucondr[0]);
  x.eval();
  array vpar = x/( ucondr[0]*(1.+af::sqrt(1.+x*x)) );
  vpar.eval();
  array one_over_ucondr_t = 1./ucondr[0];
  //new contravariant 3-velocity, v^i
  array vcon[NDIM];
  vcon[0] = zero + 1.;
  for (int m=1; m<NDIM; m++) {
    //parallel (to B) plus perpendicular (to B) velocities
    vcon[m] = vpar*Bcon[m]/bNorm + ucondr[m]*one_over_ucondr_t;
    vcon[m].eval();
  }
  array vsqr=zero;
  for(int m=0;m<NDIM;m++)
    for(int n=0;n<NDIM;n++)
      vsqr+=geom->gCov[m][n]*vcon[m]*vcon[n];
  vsqr.eval();
  array condition_v = (vsqr>=0.|| vsqr<1./geom->gCon[0][0]);
  vsqr=(1.-condition_v)*vsqr+condition_v/geom->gCon[0][0];
  vsqr.eval();
  array ut = af::sqrt(-1./vsqr);
  ut.eval();
  array utcon[NDIM];
  utcon[0]=zero;
  for(int m=1;m<NDIM;m++) {
    utcon[m] = ut*(vcon[m]-geom->gCon[0][m]/geom->gCon[0][0]);
    utcon[m].eval();
  }

  prim->vars[vars::U1] = prim->vars[vars::U1]*(1.-trans)+trans*utcon[1];
  prim->vars[vars::U2] = prim->vars[vars::U2]*(1.-trans)+trans*utcon[2];
  prim->vars[vars::U3] = prim->vars[vars::U3]*(1.-trans)+trans*utcon[3];
  prim->vars[vars::U1].eval();
  prim->vars[vars::U2].eval();
  prim->vars[vars::U3].eval();
  

  // Reset element to get Lorentz factor
  elem->set(*prim, *geom, numReads,numWrites);
  const array& lorentzFactor = elem->gammaLorentzFactor;
  array lorentzFactorSqr = lorentzFactor*lorentzFactor;

  // Now, we impose the maximum lorentz factor...
  condition = (lorentzFactorSqr>params::MaxLorentzFactor*params::MaxLorentzFactor);
  array conditionIndices = where(condition > 0);
  if(conditionIndices.elements()>0)
    {
      array MultFac = af::sqrt(af::max((lorentzFactorSqr-1.)/(params::MaxLorentzFactor*params::MaxLorentzFactor-1.),1.));
      MultFac = 1./MultFac;
      MultFac.eval();
      prim->vars[vars::U1]*=MultFac;
      prim->vars[vars::U2]*=MultFac;
      prim->vars[vars::U3]*=MultFac;
      prim->vars[vars::U1].eval();
      prim->vars[vars::U2].eval();
      prim->vars[vars::U3].eval();
    }  

  elem->set(*prim, *geom, numReads,numWrites);
  
  if(params::conduction)
    {
      const array& rho = prim->vars[vars::RHO];
      const array& cs = elem->soundSpeed;
      array Qmax = 1.07*params::ConductionClosureFactor*rho*pow(cs,3.);
      array LimFac = af::max(af::abs(elem->q)/Qmax,1.);
      prim->vars[vars::Q]=prim->vars[vars::Q]/LimFac;
      prim->vars[vars::Q].eval();
    }
  if(params::viscosity)
    {
      const array& pressure = elem->pressure;
      const array& deltaP = elem->deltaP;
      array dPmod = af::max(pressure-2./3.*deltaP,0.01*params::bSqrFloorInFluidElement)/af::max(pressure+1./3.*deltaP,params::bSqrFloorInFluidElement);
      array dPmaxPlus = af::min(1.07*params::ViscosityClosureFactor*bSqr*0.5*dPmod,1.49*pressure);
      array dPmaxMinus = af::max(-1.07*params::ViscosityClosureFactor*bSqr,-2.99*pressure);

      condition = deltaP>0.;
      prim->vars[vars::DP] = prim->vars[vars::DP]*
	(condition/af::max(deltaP/dPmaxPlus,1.)+(1.-condition)/af::max(deltaP/dPmaxMinus,1.));
      prim->vars[vars::DP].eval();
    }
  if(params::conduction || params::viscosity) elem->set(*prim, *geom, numReads,numWrites);

  af::sync();
}

void timeStepper::halfStepDiagnostics(int &numReads,int &numWrites)
{
  applyFloor(primHalfStep,elemHalfStep,geomCenter,numReads,numWrites);
}

void timeStepper::fullStepDiagnostics(int &numReads,int &numWrites)
{
  applyFloor(primOld,elemOld,geomCenter,numReads,numWrites);

  int world_rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(PETSC_COMM_WORLD, &world_size);

  af::seq domainX1 = *primOld->domainX1;
  af::seq domainX2 = *primOld->domainX2;
  af::seq domainX3 = *primOld->domainX3;
  
  // Time step control
  array minSpeedTemp,maxSpeedTemp;
  array minSpeed,maxSpeed;
  elemOld->computeMinMaxCharSpeeds(*geomCenter,directions::X1,minSpeedTemp,maxSpeedTemp,numReads,numWrites);
  minSpeedTemp = minSpeedTemp/XCoords->dX1;
  maxSpeedTemp = maxSpeedTemp/XCoords->dX1;
  minSpeed=minSpeedTemp;
  maxSpeed=maxSpeedTemp;
  if(params::dim>1)
    {
      elemOld->computeMinMaxCharSpeeds(*geomCenter,directions::X2,minSpeedTemp,maxSpeedTemp,numReads,numWrites);
      minSpeedTemp = minSpeedTemp/XCoords->dX2;
      maxSpeedTemp = maxSpeedTemp/XCoords->dX2;
      minSpeed=af::min(minSpeed,minSpeedTemp);
      maxSpeed=af::max(maxSpeed,maxSpeedTemp);
    }
  if(params::dim>2)
    {
      elemOld->computeMinMaxCharSpeeds(*geomCenter,directions::X3,minSpeedTemp,maxSpeedTemp,numReads,numWrites);
      minSpeedTemp = minSpeedTemp/XCoords->dX3;
      maxSpeedTemp = maxSpeedTemp/XCoords->dX3;
      minSpeed=af::min(minSpeed,minSpeedTemp);
      maxSpeed=af::max(maxSpeed,maxSpeedTemp);
    }
  maxSpeed = af::max(maxSpeed,af::abs(minSpeed));
  maxSpeed.eval();
  array maxInvDt_af = af::max(af::max(af::max(maxSpeed,2),1),0);
  double maxInvDt = maxInvDt_af.host<double>()[0];
  /* Use MPI to find minimum over all processors */
  if (world_rank == 0) 
    {
      double temp; 
      for(int i=1;i<world_size;i++)
	{
	  MPI_Recv(&temp, 1, MPI_DOUBLE, i, i, PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
	  if( maxInvDt < temp)
	    maxInvDt = temp;
	}
    }
  else
    {
      MPI_Send(&maxInvDt, 1, MPI_DOUBLE, 0, world_rank, PETSC_COMM_WORLD);
    }
  MPI_Barrier(PETSC_COMM_WORLD);
  MPI_Bcast(&maxInvDt,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
  MPI_Barrier(PETSC_COMM_WORLD);
  dt = params::CourantFactor/maxInvDt;
  PetscPrintf(PETSC_COMM_WORLD,"New dt = %e\n",dt);

  // On-the-fly observers
  bool ObserveData = (floor(time/params::ObserveEveryDt) != floor((time-dt)/params::ObserveEveryDt));
  bool WriteData   = (floor(time/params::WriteDataEveryDt) != floor((time-dt)/params::WriteDataEveryDt));
  if(ObserveData)
    {
      //Find maximum density
      array rhoMax_af = af::max(af::max(af::max(primOld->vars[vars::RHO](domainX1,domainX2,domainX3),2),1),0);
      double rhoMax = rhoMax_af.host<double>()[0];
      /* Communicate rhoMax to all processors */
      if (world_rank == 0) 
	{
	  double temp; 
	  for(int i=1;i<world_size;i++)
	    {
	      MPI_Recv(&temp, 1, MPI_DOUBLE, i, i, PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
	      if(rhoMax < temp)
		rhoMax = temp;
	    }
	}
      else
	{
	  MPI_Send(&rhoMax, 1, MPI_DOUBLE, 0, world_rank, PETSC_COMM_WORLD);
	}
      MPI_Barrier(PETSC_COMM_WORLD);
      MPI_Bcast(&rhoMax,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
      MPI_Barrier(PETSC_COMM_WORLD);
      
      // Find minimum beta
      const array& bSqr = elemOld->bSqr;
      const array& Pgas = elemOld->pressure;
      array PlasmaBeta = 2.*(Pgas+1.e-13)/(bSqr+1.e-18);
      array BetaMin_af = af::min(af::min(af::min(PlasmaBeta(domainX1,domainX2,domainX3),2),1),0);
      double betaMin = BetaMin_af.host<double>()[0];
      /* Use MPI to find minimum over all processors */
      if (world_rank == 0) 
	{
	  double temp; 
	  for(int i=1;i<world_size;i++)
	    {
	      MPI_Recv(&temp, 1, MPI_DOUBLE, i, i, PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
	      if( betaMin > temp)
		betaMin = temp;
	  }
	}
      else
	{
	  MPI_Send(&betaMin, 1, MPI_DOUBLE, 0, world_rank, PETSC_COMM_WORLD);
	}
      MPI_Barrier(PETSC_COMM_WORLD);
      MPI_Bcast(&betaMin,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
      MPI_Barrier(PETSC_COMM_WORLD);
      
      /* Get the domain of the bulk and volume element*/
      elemOld->computeFluxes(*geomCenter,0, *consOld,numReads,numWrites);
      double volElem = XCoords->dX1;
      if(params::dim>1)
	volElem*=XCoords->dX2;
      if(params::dim>2)
	volElem*=XCoords->dX3;

      // Integrate baryon mass
      array MassIntegrand = consOld->vars[vars::RHO]*volElem;
      array BaryonMass_af = af::sum(af::flat(MassIntegrand(domainX1, domainX2, domainX3)),0);
      double BaryonMass = BaryonMass_af.host<double>()[0];
      /* Communicate to all processors */
      if (world_rank == 0) 
	{
	  double temp; 
	  for(int i=1;i<world_size;i++)
	    {
	      MPI_Recv(&temp, 1, MPI_DOUBLE, i, i, PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
	      BaryonMass += temp;
	    }
	}
      else
	{
	  MPI_Send(&BaryonMass, 1, MPI_DOUBLE, 0, world_rank, PETSC_COMM_WORLD);
	}
      MPI_Barrier(PETSC_COMM_WORLD);
      MPI_Bcast(&BaryonMass,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
      MPI_Barrier(PETSC_COMM_WORLD);
      // Integrate magnetic energy
      array EMagIntegrand = volElem*elemOld->bSqr*0.5*geomCenter->g;
      array EMag_af = af::sum(af::flat(EMagIntegrand(domainX1, domainX2, domainX3)),0);
      double EMag = EMag_af.host<double>()[0];
      /* Communicate to all processors */
      if (world_rank == 0) 
	{
	  double temp; 
	  for(int i=1;i<world_size;i++)
	    {
	      MPI_Recv(&temp, 1, MPI_DOUBLE, i, i, PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
	      EMag += temp;
	    }
	}
      else
	{
	  MPI_Send(&EMag, 1, MPI_DOUBLE, 0, world_rank, PETSC_COMM_WORLD);
	}
      MPI_Barrier(PETSC_COMM_WORLD);
      MPI_Bcast(&EMag,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
      MPI_Barrier(PETSC_COMM_WORLD);

      PetscPrintf(PETSC_COMM_WORLD,"Global quantities at t = %e\n",time);
      PetscPrintf(PETSC_COMM_WORLD,"rhoMax = %e; betaMin = %e;\n",rhoMax,betaMin);
      PetscPrintf(PETSC_COMM_WORLD,"Baryon Mass = %e; Magnetic Energy = %e\n",BaryonMass,EMag);
    }
  if(WriteData)
    {
      const int WriteIdx = floor(time/params::WriteDataEveryDt);
      if(WriteIdx==0)
	{
	  PetscPrintf(PETSC_COMM_WORLD, "Printing gCov\n");
	  geomCenter->setgCovGrid();
	  geomCenter->gCovGrid->dump("gCov","gCov.h5");
	  geomCenter->setgConGrid();
	  geomCenter->gConGrid->dump("gCon","gCon.h5");
	  geomCenter->setgGrid();
	  geomCenter->gGrid->dump("sqrtDetg","sqrtDetg.h5");
	  geomCenter->setxCoordsGrid();
	  geomCenter->xCoordsGrid->dump("xCoords","xCoords.h5");
	}
      std::string filename = "primVarsT";
      std::string s_idx = std::to_string(WriteIdx);
      for(int i=0;i<6-s_idx.size();i++)
	filename=filename+"0";
      filename=filename+s_idx;
      filename=filename+".h5";
      primOld->dump("primitives",filename);
    }
}


void inflowCheck(grid& primBC,fluidElement& elemBC,
		 const geometry& geom, int &numReads,int &numWrites)
{
  const int numGhost = params::numGhost;
  if(primBC.iLocalEnd == primBC.N1)
    {
      af::seq domainX1RightBoundary(primBC.N1Local+numGhost,
				    primBC.N1Local+2*numGhost-1
				    );
      elemBC.set(primBC,geom,numReads,numWrites);
      
      // Prefactor the lorentz factor
      primBC.vars[vars::U1](domainX1RightBoundary,span,span)
	= primBC.vars[vars::U1](domainX1RightBoundary,span,span)
	/elemBC.gammaLorentzFactor(domainX1RightBoundary,span,span);
      primBC.vars[vars::U2](domainX1RightBoundary,span,span)
	= primBC.vars[vars::U2](domainX1RightBoundary,span,span)
	/elemBC.gammaLorentzFactor(domainX1RightBoundary,span,span);
      primBC.vars[vars::U3](domainX1RightBoundary,span,span)
	= primBC.vars[vars::U3](domainX1RightBoundary,span,span)
	/elemBC.gammaLorentzFactor(domainX1RightBoundary,span,span);
      // Reset radial velocity if it is too small (i.e. incoming)
      primBC.vars[vars::U1].eval();
      primBC.vars[vars::U1](domainX1RightBoundary,span,span)
	= af::max(geom.gCon[0][1]*geom.alpha,
		  primBC.vars[vars::U1])
	(domainX1RightBoundary,span,span);
      primBC.vars[vars::U1].eval();
      primBC.vars[vars::U2].eval();
      primBC.vars[vars::U3].eval();
      // Recompute lorentz factor
      array vSqr = primBC.vars[vars::U1](domainX1RightBoundary,span,span)*0.;
      for(int i=0;i<3;i++)
	for(int j=0;j<3;j++)
	  vSqr += geom.gCov[i+1][j+1](domainX1RightBoundary,span,span)*
	    primBC.vars[vars::U1+i](domainX1RightBoundary,span,span)*
	    primBC.vars[vars::U1+j](domainX1RightBoundary,span,span);
      double vSqrMax = 1.-1./params::MaxLorentzFactor/params::MaxLorentzFactor;
      double vSqrMin = 1.e-13;
      vSqr = af::max(af::min(vSqr,vSqrMax),vSqrMin);
      array newLorentzFactor = 1./sqrt(1.-vSqr);
      newLorentzFactor.eval();
      for(int i=0;i<3;i++)
	{
	  primBC.vars[vars::U1+i](domainX1RightBoundary,span,span)=
	    primBC.vars[vars::U1+i](domainX1RightBoundary,span,span)
	    *newLorentzFactor;
	  primBC.vars[vars::U1+i].eval();
	}
      elemBC.set(primBC,geom,numReads,numWrites);
    }
  if(primBC.iLocalStart == 0)
    {
      af::seq domainX1LeftBoundary(0,
                                    numGhost-1
                                    );
      elemBC.set(primBC,geom,numReads,numWrites);

      // Prefactor the lorentz factor
      primBC.vars[vars::U1](domainX1LeftBoundary,span,span)
        = primBC.vars[vars::U1](domainX1LeftBoundary,span,span)
        /elemBC.gammaLorentzFactor(domainX1LeftBoundary,span,span);
      primBC.vars[vars::U2](domainX1LeftBoundary,span,span)
        = primBC.vars[vars::U2](domainX1LeftBoundary,span,span)
        /elemBC.gammaLorentzFactor(domainX1LeftBoundary,span,span);
      primBC.vars[vars::U3](domainX1LeftBoundary,span,span)
        = primBC.vars[vars::U3](domainX1LeftBoundary,span,span)
        /elemBC.gammaLorentzFactor(domainX1LeftBoundary,span,span);
      // Reset radial velocity if it is too small (i.e. incoming)
      primBC.vars[vars::U1].eval();
      primBC.vars[vars::U1](domainX1LeftBoundary,span,span)
        = af::min(geom.gCon[0][1]*geom.alpha,
                  primBC.vars[vars::U1])
        (domainX1LeftBoundary,span,span);
      primBC.vars[vars::U1].eval();
      primBC.vars[vars::U2].eval();
      primBC.vars[vars::U3].eval();
      // Recompute lorentz factor
      array vSqr = primBC.vars[vars::U1](domainX1LeftBoundary,span,span)*0.;
      for(int i=0;i<3;i++)
	for(int j=0;j<3;j++)
	  vSqr += geom.gCov[i+1][j+1](domainX1LeftBoundary,span,span)*
	    primBC.vars[vars::U1+i](domainX1LeftBoundary,span,span)*
	    primBC.vars[vars::U1+j](domainX1LeftBoundary,span,span);
      double vSqrMax = 1.-1./params::MaxLorentzFactor/params::MaxLorentzFactor;
      double vSqrMin = 1.e-13;
      vSqr = af::max(af::min(vSqr,vSqrMax),vSqrMin);
      array newLorentzFactor = 1./sqrt(1.-vSqr);
      newLorentzFactor.eval();
      for(int i=0;i<3;i++)
	{
	  primBC.vars[vars::U1+i](domainX1LeftBoundary,span,span)
	    =primBC.vars[vars::U1+i](domainX1LeftBoundary,span,span)
	    *newLorentzFactor;
	  primBC.vars[vars::U1+i].eval();
	}
      elemBC.set(primBC,geom,numReads,numWrites);
    }
}

void fixPoles(grid& primBC, int &numReads,int &numWrites)
{
  const int numGhost = params::numGhost;
  if(primBC.jLocalStart == 0)
    {
      int idx0 = numGhost;
      int idx1 = numGhost+1;
      int idx2 = numGhost+2;
      double ic0 = 0.2;
      double ic1 = 0.6;
      primBC.vars[vars::RHO](span,idx0,span)=
	primBC.vars[vars::RHO](span,idx2,span);
      primBC.vars[vars::U](span,idx0,span)=
	primBC.vars[vars::U](span,idx2,span);
      primBC.vars[vars::U1](span,idx0,span)=
	primBC.vars[vars::U1](span,idx2,span);
      primBC.vars[vars::U2](span,idx0,span)=
        primBC.vars[vars::U2](span,idx2,span)*ic0;
      primBC.vars[vars::U3](span,idx0,span)=
        primBC.vars[vars::U3](span,idx2,span);
      if(params::conduction)
	{
	  primBC.vars[vars::Q](span,idx0,span)=
	    primBC.vars[vars::Q](span,idx2,span)*ic0;
	}
      if(params::viscosity)
	{
	  primBC.vars[vars::DP](span,idx0,span)=
	    primBC.vars[vars::DP](span,idx2,span)*ic0;
	}

      primBC.vars[vars::RHO](span,idx1,span)=
	primBC.vars[vars::RHO](span,idx2,span);
      primBC.vars[vars::U](span,idx1,span)=
	primBC.vars[vars::U](span,idx2,span);
      primBC.vars[vars::U1](span,idx1,span)=
	primBC.vars[vars::U1](span,idx2,span);
      primBC.vars[vars::U2](span,idx1,span)=
        primBC.vars[vars::U2](span,idx2,span)*ic1;
      primBC.vars[vars::U3](span,idx1,span)=
        primBC.vars[vars::U3](span,idx2,span);
      if(params::conduction)
	{
	  primBC.vars[vars::Q](span,idx1,span)=
	    primBC.vars[vars::Q](span,idx2,span)*ic1;
	}
      if(params::viscosity)
        {
	  primBC.vars[vars::DP](span,idx1,span)=
	    primBC.vars[vars::DP](span,idx2,span)*ic1;
	}

      af::seq domainX2LeftBoundary(0,
                                    numGhost-1
				   );
      primBC.vars[vars::U2](span,domainX2LeftBoundary,span)
	= primBC.vars[vars::U2](span,domainX2LeftBoundary,span)*(-1.0);
      primBC.vars[vars::B2](span,domainX2LeftBoundary,span)
	= primBC.vars[vars::B2](span,domainX2LeftBoundary,span)*(-1.0);

      for(int var=0;var<vars::dof;var++)
	primBC.vars[var].eval();
    }
  if(primBC.jLocalEnd == primBC.N2)
    {
      int idx0 = primBC.N2Local+numGhost-1;
      int idx1 = idx0-1;
      int idx2 = idx0-2;
      double ic0 = 0.2;
      double ic1 = 0.6;
      primBC.vars[vars::RHO](span,idx0,span)=
	primBC.vars[vars::RHO](span,idx2,span);
      primBC.vars[vars::U](span,idx0,span)=
	primBC.vars[vars::U](span,idx2,span);
      primBC.vars[vars::U1](span,idx0,span)=
	primBC.vars[vars::U1](span,idx2,span);
      primBC.vars[vars::U2](span,idx0,span)=
        primBC.vars[vars::U2](span,idx2,span)*ic0;
      primBC.vars[vars::U3](span,idx0,span)=
        primBC.vars[vars::U3](span,idx2,span);
      if(params::conduction)
        {
          primBC.vars[vars::Q](span,idx0,span)=
            primBC.vars[vars::Q](span,idx2,span)*ic0;
        }
      if(params::viscosity)
        {
          primBC.vars[vars::DP](span,idx0,span)=
            primBC.vars[vars::DP](span,idx2,span)*ic0;
        }

      primBC.vars[vars::RHO](span,idx1,span)=
	primBC.vars[vars::RHO](span,idx2,span);
      primBC.vars[vars::U](span,idx1,span)=
	primBC.vars[vars::U](span,idx2,span);
      primBC.vars[vars::U1](span,idx1,span)=
	primBC.vars[vars::U1](span,idx2,span);
      primBC.vars[vars::U2](span,idx1,span)=
        primBC.vars[vars::U2](span,idx2,span)*ic1;
      primBC.vars[vars::U3](span,idx1,span)=
        primBC.vars[vars::U3](span,idx2,span);
      if(params::conduction)
        {
          primBC.vars[vars::Q](span,idx1,span)=
            primBC.vars[vars::Q](span,idx2,span)*ic1;
        }
      if(params::viscosity)
        {
          primBC.vars[vars::DP](span,idx1,span)=
            primBC.vars[vars::DP](span,idx2,span)*ic1;
        }

      af::seq domainX2RightBoundary(primBC.N2Local+numGhost,
                                    primBC.N2Local+2*numGhost-1
				   );
      primBC.vars[vars::U2](span,domainX2RightBoundary,span)
	= primBC.vars[vars::U2](span,domainX2RightBoundary,span)*(-1.0);
      primBC.vars[vars::B2](span,domainX2RightBoundary,span)
	= primBC.vars[vars::B2](span,domainX2RightBoundary,span)*(-1.0);

      for(int var=0;var<vars::dof;var++)
	primBC.vars[var].eval();
    }
}

void timeStepper::setProblemSpecificBCs(int &numReads,int &numWrites)
{
  // 1) Choose which primitive variables are corrected.
  grid* primBC;
  fluidElement* elemBC;
  if(currentStep == timeStepperSwitches::HALF_STEP)
    {
      primBC = primOld;
      elemBC = elemOld;
    }
  else
    {
      primBC = primHalfStep;
      elemBC = elemHalfStep;
    }
  // 2) Check that there is no inflow at the radial boundaries
  inflowCheck(*primBC,*elemBC,*geomCenter,
  	      numReads,numWrites);

  // 3) 'Fix' the polar regions by correcting the firs
  // two active zones
  fixPoles(*primBC,numReads,numWrites);
  af::sync();
};

void timeStepper::applyProblemSpecificFluxFilter(int &numReads,int &numWrites)
{
  const int numGhost = params::numGhost;
  // Prevents matter from flowing into the computational domain
  if(primOld->iLocalStart == 0)
    {
      int idx = numGhost;
      fluxesX1->vars[vars::RHO](idx,span,span)=
	af::min(fluxesX1->vars[vars::RHO](idx,span,span),0.);
      fluxesX1->vars[vars::RHO].eval();
    }
  if(primOld->iLocalEnd == primOld->N1)
    {
      int idx = primOld->N1Local+numGhost;
      fluxesX1->vars[vars::RHO](idx,span,span)=
        af::max(fluxesX1->vars[vars::RHO](idx,span,span),0.);
      fluxesX1->vars[vars::RHO].eval();
    }
  
  // Set fluxes to 0 on the polar axis
  if(primOld->jLocalStart == 0)
    {
      int idx = numGhost;
      for(int var=0;var<vars::dof;var++)
	{
	  fluxesX2->vars[var](span,idx,span)=0.;
	  fluxesX2->vars[var].eval();
	}
      fluxesX1->vars[vars::B2](span,idx-1,span)=fluxesX1->vars[vars::B2](span,idx,span)*(-1.0);
      fluxesX1->vars[vars::B2].eval();
    }
  if(primOld->jLocalEnd == primOld->N2)
    {
      int idx = primOld->N2Local+numGhost;
      for(int var=0;var<vars::dof;var++)
	{
	  fluxesX2->vars[var](span,idx,span)=0.;
	  fluxesX2->vars[var].eval();
	}
      fluxesX1->vars[vars::B2](span,idx,span)=fluxesX1->vars[vars::B2](span,idx-1,span)*(-1.0);
      fluxesX1->vars[vars::B2].eval();
    }
}
