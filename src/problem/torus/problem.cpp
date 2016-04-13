#include "torus.hpp"

void fluidElement::setFluidElementParameters(const geometry &geom)
{
  array xCoords[3],XCoords[3];
  geom.getXCoords(XCoords);
  geom.XCoordsToxCoords(XCoords,xCoords);
  array Radius = xCoords[0];
  array DynamicalTimescale = af::pow(Radius,1.5);
  DynamicalTimescale.eval();
  tau = DynamicalTimescale;

  if(params::conduction)
    {
      array Qmax = params::ConductionClosureFactor*rho*pow(soundSpeed,3.);
      double lambda = 0.01;
      array yCon = q/Qmax;
      yCon = af::exp(-(yCon-1.)/lambda);
      yCon.eval();
      array fdCon = yCon/(yCon+1.)+1.e-5;
      tau=af::min(tau,DynamicalTimescale*fdCon);
    }
  if(params::viscosity)
    {
      array dPmod = af::max(pressure-2./3.*deltaP,params::bSqrFloorInFluidElement)/af::max(pressure+1./3.*deltaP,params::bSqrFloorInFluidElement);
      array dPmaxPlus = af::min(params::ViscosityClosureFactor*bSqr*0.5*dPmod,1.49*pressure/1.07);
      array dPmaxMinus = af::max(-params::ViscosityClosureFactor*bSqr,-2.99*pressure/1.07);

      array condition = deltaP>0.;
      array dPmax = condition*dPmaxPlus + (1.-condition)*dPmaxMinus;

      double lambda = 0.01;
      array yVis = af::abs(deltaP)/(af::abs(dPmax)+params::bSqrFloorInFluidElement);
      yVis = af::exp(-(yVis-1.)/lambda);
      yVis.eval();
      array fdVis = yVis/(yVis+1.)+1.e-5;
      tau=af::min(tau,DynamicalTimescale*fdVis);
    }
  tau.eval();
  chi_emhd = params::ConductionAlpha*soundSpeed*soundSpeed*tau;
  nu_emhd  = params::ViscosityAlpha*soundSpeed*soundSpeed*tau;
  chi_emhd.eval();
  nu_emhd.eval();
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
  geomCenter->XCoordsToxCoords(XCoords->vars, xCoords);

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
  /*for(int proc=0;proc<world_size;proc++)
    {
      if(world_rank==proc)
	{
	  printf("Local size on proc %i : %i x %i x %i\n",proc,N1g,N2g,N3g);
	  af_print(xCoords[directions::X1],5);
	}
      MPI_Barrier(PETSC_COMM_WORLD);
      }*/
  
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
    +af::shift(Rho_af,0,1,0)+af::shift(Rho_af,0,-1,0)
     +af::shift(Rho_af,0,0,1)+af::shift(Rho_af,0,0,-1))
    /6.;
  array zero = rhoAvg*0.;
  array Avec = af::max(rhoAvg-0.2,zero)*
    af::cos(xCoords[directions::X2]*
	    (params::MagneticLoops-1));
  Avec.eval();
 
  // Compute magnetic field 
  const array& g = geomCenter->g;
  const double dX1 = XCoords->dX1;
  const double dX2 = XCoords->dX2;
//  primOld->vars[vars::B1] = 
//    (af::shift(Avec,0,-1,0)-af::shift(Avec,0,0,0)
//     +af::shift(Avec,-1,-1,0)-af::shift(Avec,-1,0,0))/
//    (2.*dX2*g);
//  primOld->vars[vars::B2] = 
//    (af::shift(Avec,0,0,0)-af::shift(Avec,-1,0,0)
//     +af::shift(Avec,0,-1,0)-af::shift(Avec,-1,-1,0))/
//    (2.*dX1*g);
//  primOld->vars[vars::B1].eval();
//  primOld->vars[vars::B2].eval();

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
    array PlasmaBeta = (Pgas+1.e-13)/(bSqr+1.e-18);
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

  applyFloor(primOld,elemOld,geomCenter,XCoords,numReads,numWrites);

  for (int var=0; var<vars::dof; var++) 
    {
        primOld->vars[var].eval();
    }

  af::sync();
}

void applyFloor(grid* prim, fluidElement* elem, geometry* geom, grid* XCoords, int &numReads,int &numWrites)
{
  array xCoords[3];
  geom->XCoordsToxCoords(XCoords->vars,xCoords);
  array Radius = xCoords[0];

  array minRho = af::pow(Radius,params::RhoFloorSlope)*params::RhoFloorAmpl;
  array minU = af::pow(Radius,params::UFloorSlope)*params::UFloorAmpl;

  array condition = prim->vars[vars::RHO]<minRho;
  prim->vars[vars::RHO] = condition*minRho+(1.-condition)*prim->vars[vars::RHO];

  condition = prim->vars[vars::U]<minU;
  prim->vars[vars::U] = condition*minU+(1.-condition)*prim->vars[vars::U];

  elem->set(*prim,*geom,numReads,numWrites);

  const array& bSqr = elem->bSqr;
  condition = bSqr>params::BsqrOverRhoMax*prim->vars[vars::RHO];
  prim->vars[vars::RHO] = prim->vars[vars::RHO]*(1.-condition)+condition*bSqr/params::BsqrOverRhoMax;
  condition = bSqr>params::BsqrOverUMax*prim->vars[vars::U];
  prim->vars[vars::U] = prim->vars[vars::U]*(1.-condition)+condition*bSqr/params::BsqrOverUMax;
  
  prim->vars[vars::RHO].eval();
  prim->vars[vars::U].eval();

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

      array condition = deltaP>0.;
      prim->vars[vars::DP] = prim->vars[vars::DP]*
	(condition/af::max(deltaP/dPmaxPlus,1.)+(1.-condition)/af::max(deltaP/dPmaxMinus,1.));
      prim->vars[vars::DP].eval();
    }

  elem->set(*prim, *geom, numReads,numWrites);
  af::sync();
}

void timeStepper::halfStepDiagnostics(int &numReads,int &numWrites)
{
  applyFloor(primHalfStep,elemHalfStep,geomCenter,XCoords,numReads,numWrites);
}

void timeStepper::fullStepDiagnostics(int &numReads,int &numWrites)
{
  applyFloor(primOld,elemOld,geomCenter,XCoords,numReads,numWrites);

  int world_rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(PETSC_COMM_WORLD, &world_size);
  
  // Time step control
  array minSpeedTemp,maxSpeedTemp;
  array minSpeed,maxSpeed;
  elemOld->computeMinMaxCharSpeeds(*geomCenter,1,minSpeedTemp,maxSpeedTemp,numReads,numWrites);
  array dX = XCoords->dX1*af::sqrt(geomCenter->gCov[1][1]);
  minSpeedTemp = minSpeedTemp/dX;
  maxSpeedTemp = maxSpeedTemp/dX;
  minSpeed=minSpeedTemp;
  maxSpeed=maxSpeedTemp;
  if(params::dim>1)
    {
      elemOld->computeMinMaxCharSpeeds(*geomCenter,2,minSpeedTemp,maxSpeedTemp,numReads,numWrites);
      dX = XCoords->dX2*af::sqrt(geomCenter->gCov[2][2]);
      minSpeedTemp = minSpeedTemp/dX;
      maxSpeedTemp = maxSpeedTemp/dX;
      minSpeed=af::min(minSpeed,minSpeedTemp);
      maxSpeed=af::max(maxSpeed,maxSpeedTemp);
    }
  if(params::dim>2)
    {
      elemOld->computeMinMaxCharSpeeds(*geomCenter,3,minSpeedTemp,maxSpeedTemp,numReads,numWrites);
      dX = XCoords->dX3*af::sqrt(geomCenter->gCov[3][3]);
      minSpeedTemp = minSpeedTemp/dX;
      maxSpeedTemp = maxSpeedTemp/dX;
      minSpeed=af::min(minSpeed,minSpeedTemp);
      maxSpeed=af::max(maxSpeed,maxSpeedTemp);
    }
  maxSpeed = af::max(maxSpeed,af::abs(minSpeed));
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
      array rhoMax_af = af::max(af::max(af::max(primOld->vars[vars::RHO],2),1),0);
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
      array PlasmaBeta = (Pgas+1.e-13)/(bSqr+1.e-18);
      array BetaMin_af = af::min(af::min(af::min(PlasmaBeta,2),1),0);
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
      af::seq domainX1 = *primOld->domainX1;
      af::seq domainX2 = *primOld->domainX2;
      af::seq domainX3 = *primOld->domainX3;
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
      if(WriteIdx==1)
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
      std::string filename = "primVarsT"+std::to_string(WriteIdx)+".h5";
      primOld->dump("primitives",filename);
    }
}

void timeStepper::setProblemSpecificBCs(int &numReads,int &numWrites)
{
};
