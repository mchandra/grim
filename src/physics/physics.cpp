#include "physics.hpp"

fluidElement::fluidElement(const grid &prim,
                           const geometry &geom,
                           int &numReads,
                           int &numWrites
                          )
{
  /* Use this to set various fluid parameters. For ex: tau = 0.1*one etc..*/
  one = af::constant(1, 
                     prim.vars[0].dims(directions::X1),
                     prim.vars[0].dims(directions::X2),
                     prim.vars[0].dims(directions::X3),
		                 f64
           		      );

  array zero = 0.*one;
  rho = zero;
  u = zero;
  u1 = zero; u2=zero; u3=zero;
  B1=zero; B2=zero; B3=zero;
  pressure = zero; temperature = zero;
  gammaLorentzFactor = zero;
  bSqr = zero; soundSpeed = zero;
  bNorm = zero;
  
  for(int d=0;d<NDIM;d++)
    {
      uCon[d]=zero; uCov[d]=zero;
      bCon[d]=zero; bCov[d]=zero;
      NUp[d]=zero;
      for(int dd=0;dd<NDIM;dd++)
	{
	  TUpDown[d][dd]=zero;
	}
    }
  
  /* Allocate memory for gradients used in EMHD */
  if (params::conduction || params::viscosity)
  {
    divuCov = zero;
    
    for(int mu=0;mu<NDIM;mu++)
    {
      gradT[mu] = zero;
      dtuCov[mu] = zero;
      for(int nu=0;nu<NDIM;nu++)
	{
	  graduCov[nu][mu] = zero;
	}
    }
    
    deltaP0 = zero;
    q0 = zero;
    tau = zero;
    chi_emhd = zero;
    nu_emhd = zero;
    qTilde=zero; deltaPTilde=zero;
    q=zero; deltaP=zero;
  }

  set(prim, geom, numReads, numWrites);
}


void fluidElement::set(const grid &prim,
		       const geometry &geom,
		       int &numReads,
		       int &numWrites)
{
  array AllIdx = where(one>0.);
  AllIdx.eval();
  set(prim,geom,numReads,numWrites,AllIdx);
}

void fluidElement::set(const grid &prim,
                       const geometry &geom,
                       int &numReads,
                       int &numWrites,
		       const array Idx
                      )
{
  rho(Idx) = af::max(prim.vars[vars::RHO](Idx),params::rhoFloorInFluidElement);
  u(Idx)   = af::max(prim.vars[vars::U  ](Idx),params::uFloorInFluidElement);
  u1(Idx)  = prim.vars[vars::U1 ](Idx);
  u2(Idx)  = prim.vars[vars::U2 ](Idx);
  u3(Idx)  = prim.vars[vars::U3 ](Idx);
  B1(Idx)  = prim.vars[vars::B1 ](Idx);
  B2(Idx)  = prim.vars[vars::B2 ](Idx);
  B3(Idx)  = prim.vars[vars::B3 ](Idx);

  pressure(Idx)    = (params::adiabaticIndex - 1.)*u(Idx);
  temperature(Idx) = af::max(pressure(Idx)/rho(Idx),params::temperatureFloorInFluidElement);

  soundSpeed(Idx)  = af::sqrt( params::adiabaticIndex*pressure(Idx)
			       /(rho(Idx)+params::adiabaticIndex*u(Idx))
                        );
  
  gammaLorentzFactor(Idx) =
    af::sqrt(1 + geom.gCov[1][1](Idx) * u1(Idx) * u1(Idx)
	     + geom.gCov[2][2](Idx) * u2(Idx) * u2(Idx)
	     + geom.gCov[3][3](Idx) * u3(Idx) * u3(Idx)

             + 2*(  geom.gCov[1][2](Idx) * u1(Idx) * u2(Idx)
		    + geom.gCov[1][3](Idx) * u1(Idx) * u3(Idx)
		    + geom.gCov[2][3](Idx) * u2(Idx) * u3(Idx)
                 )
            );
  gammaLorentzFactor(Idx).eval(); 
  /* Reads:
   * -----
   * gCov[1][1], gCov[2][2], gCov[3][3], gCov[1][2], gCov[1][3], gCov[2][3]: 6
   * u1, u2, u3 : 3
   *
   * Writes:
   * ------
   * gammaLorentzFactor : 1 */

  uCon[0](Idx) = gammaLorentzFactor(Idx)/geom.alpha(Idx);
  uCon[0](Idx).eval(); 
  /* Reads:
   * -----
   * gammaLorentzFactor, alpha : 2
   *
   * Writes:
   * ------
   * uCon[0] : 1 */


  uCon[1](Idx) = u1(Idx) - gammaLorentzFactor(Idx)*geom.gCon[0][1](Idx)*geom.alpha(Idx);
  uCon[1](Idx).eval(); 
  /* Reads:
   * -----
   * u1, gammaLorentzFactor, gCon[0][1], alpha : 4
   *
   * Writes:
   * ------
   * uCon[1] : 1 */

  uCon[2](Idx) = u2(Idx) - gammaLorentzFactor(Idx)*geom.gCon[0][2](Idx)*geom.alpha(Idx);
  uCon[2](Idx).eval(); 
  /* Reads:
   * -----
   * u2, gammaLorentzFactor, gCon[0][2], alpha : 4
   *
   * Writes:
   * ------
   * uCon[2] : 1 */

  uCon[3](Idx) = u3(Idx) - gammaLorentzFactor(Idx)*geom.gCon[0][3](Idx)*geom.alpha(Idx);
  uCon[3](Idx).eval();
  /* Reads:
   * -----
   * u3, gammaLorentzFactor, gCon[0][3], alpha : 4
   *
   * Writes:
   * ------
   * uCon[3] : 1 */

  for (int mu=0; mu < NDIM; mu++)
  {
    uCov[mu](Idx) =  0.;
    for(int nu=0;nu<NDIM;nu++)
      {
	uCov[mu](Idx)+=geom.gCov[mu][nu](Idx)*uCon[nu](Idx);
      }
    uCov[mu](Idx).eval(); 
  } 
  /* Reads:
   * -----
   * gCov[mu][0], gCov[mu][1], gCov[mu][2], gCov[mu][3]: 16
   * uCon[0], uCon[1], uCon[2], uCon[3]: 4 x 4 = 16
   *
   * Writes:
   * ------
   * uCov[mu] : 4 */

  bCon[0](Idx) =  B1(Idx)*uCov[1](Idx) + B2(Idx)*uCov[2](Idx) + B3(Idx)*uCov[3](Idx);
  bCon[0](Idx).eval();
  /* Reads:
   * -----
   * B1, B2, B3, uCov[1], uCov[2], uCov[3] : 6
   *
   * Writes:
   * ------
   * bCon[0] : 1 */

  bCon[1](Idx) = (B1(Idx) + bCon[0](Idx) * uCon[1](Idx))/uCon[0](Idx);
  bCon[1](Idx).eval();
  /* Reads:
   * -----
   * B1, bCon[0], uCon[1], uCon[0]: 4
   *
   * Writes:
   * ------
   * bCon[1] : 1 */

  bCon[2](Idx) = (B2(Idx) + bCon[0](Idx) * uCon[2](Idx))/uCon[0](Idx);
  bCon[2](Idx).eval();
  /* Reads:
   * -----
   * B2, bCon[0], uCon[2], uCon[0]: 4
   *
   * Writes:
   * ------
   * bCon[2] : 1 */

  bCon[3](Idx) = (B3(Idx) + bCon[0](Idx) * uCon[3](Idx))/uCon[0](Idx);
  bCon[3](Idx).eval();
  /* Reads:
   * -----
   * B3, bCon[0], uCon[3], uCon[0]: 4
   *
   * Writes:
   * ------
   * bCon[3] : 1 */

  for (int mu=0; mu < NDIM; mu++)
  {
    bCov[mu](Idx)=0.;
    for(int nu=0;nu<NDIM;nu++)
      {
	bCov[mu](Idx)+=geom.gCov[mu][nu](Idx)*bCon[nu](Idx);
      }
    bCov[mu](Idx).eval();
  }
  /* Reads:
   * -----
   * gCov[mu][0], gCov[mu][1], gCov[mu][2], gCov[mu][3]: 16
   * bCon[0], bCon[1], bCon[2], bCon[3]: 4 x 4 = 16
   *
   * Writes:
   * ------
   * bCov[mu] : 4 */

  bSqr(Idx) = 0.;
  for (int mu=0; mu < NDIM; mu++)
    {
      bSqr(Idx)+=bCon[mu](Idx)*bCov[mu](Idx);
    }
  bSqr(Idx) += params::bSqrFloorInFluidElement;
  bSqr.eval();
  /* Reads:
   * -----
   * bCon[0], bCon[1], bCon[2], bCon[3]: 4
   * bCov[0], bCov[1], bCov[2], bCov[3]: 4
   *
   * Writes:
   * ------
   * bSqr : 1 */

  bNorm(Idx) = af::sqrt(bSqr(Idx));

  // Note: this needs to be before setFluidElementParameters
  // because the closure relation uses q, deltaP!
  if (params::conduction==1)
  {
    qTilde(Idx) = prim.vars[vars::Q](Idx);

    if (params::highOrderTermsConduction==1)
    {
      q(Idx) = qTilde(Idx) * temperature(Idx) * af::sqrt(rho(Idx)*params::ConductionAlpha*soundSpeed(Idx)*soundSpeed(Idx));
      q(Idx).eval();
    }
    else
    {
      q(Idx) = qTilde(Idx);
    }
  }

  if (params::viscosity==1)
  {
    deltaPTilde(Idx) = prim.vars[vars::DP](Idx);

    if (params::highOrderTermsViscosity == 1)
    {
      deltaP(Idx) = deltaPTilde(Idx) * af::sqrt(temperature(Idx) * rho(Idx) * params::ViscosityAlpha*soundSpeed(Idx)*soundSpeed(Idx));
      deltaP(Idx).eval();
    }
    else
    {
      deltaP(Idx) = deltaPTilde(Idx);
    }
  }

  // Note: this uses q, deltaP, bSqr!
  setFluidElementParameters(geom,Idx);

  for (int mu=0; mu < NDIM; mu++)
  {
    NUp[mu](Idx) = rho(Idx) * uCon[mu](Idx);

    for (int nu=0; nu < NDIM; nu++)
    {
      TUpDown[mu][nu](Idx) =   (rho(Idx) + u(Idx) + pressure(Idx) + bSqr(Idx))*uCon[mu](Idx)*uCov[nu](Idx)
	+ (pressure(Idx) + 0.5*bSqr(Idx))*DELTA(mu, nu)
	- bCon[mu](Idx) * bCov[nu](Idx);

      
      if (params::conduction==1)
      {
        TUpDown[mu][nu](Idx) += q(Idx)/bNorm(Idx) * (uCon[mu](Idx)*bCov[nu](Idx) + bCon[mu](Idx)*uCov[nu](Idx));
      }

      if (params::viscosity==1)
      {
        TUpDown[mu][nu](Idx) += (deltaP(Idx)*(-1.0))       
	  * (  bCon[mu](Idx) * bCov[nu](Idx)/bSqr(Idx)
	       - (1./3.)*(DELTA(mu, nu) + uCon[mu](Idx)*uCov[nu](Idx))
                             );
      }
      
      TUpDown[mu][nu](Idx).eval();
      /* Reads:
       * -----
       * rho, u, bSqr, q, deltaP: 5 x 16 = 80
       * uCon[mu], uCov[nu], bCon[mu], bCov[nu]: 4 x 16 = 64
       *
       * Writes:
       * ------
       * TUpDown[mu][nu] : 16 */
    }
    NUp[mu](Idx).eval();
    /* Reads:
     * -----
     * rho : 1 x 4 = 4
     * uCon[mu] : 4
     *
     * Writes:
     * ------
     * NUp[mu] : 4 */
  }
  /* Total reads : 265
   * Total writes: 38 */ 

  numReads  = 265;
  numWrites = 38;
  
  if (params::highOrderTermsConduction)
  {
    numReads  += 1;
    numWrites += 1;
  }

  if (params::highOrderTermsViscosity)
  {
    numReads  += 1;
    numWrites += 1;
  }
}

void fluidElement::computeFluxes(const geometry &geom, 
                                 const int dir,
                                 grid &flux,
                                 int &numReads,
                                 int &numWrites
                                )
{
  array AllIdx = where(one>0.);
  AllIdx.eval();
  computeFluxes(geom,dir,flux,numReads,numWrites,AllIdx);
}

void fluidElement::computeFluxes(const geometry &geom, 
                                 const int dir,
                                 grid &flux,
                                 int &numReads,
                                 int &numWrites,
				 const array Idx
                                )
{
  array g = geom.g;

  flux.vars[vars::RHO](Idx) = g(Idx)*NUp[dir](Idx);
  flux.vars[vars::RHO](Idx).eval();
  /* Reads:
   * -----
   * g, NUp[dir] : 2
   *
   * Writes:
   * ------
   * flux[vars::RHO] : 1 */

  flux.vars[vars::U](Idx)   = g(Idx)*TUpDown[dir][0](Idx) + flux.vars[vars::RHO](Idx);
  flux.vars[vars::U](Idx).eval();
  /* Reads:
   * -----
   * g, TUpDown[dir][0], flux[vars::RHO] : 3
   *
   * Writes:
   * ------
   * flux[vars::U] : 1 */

  flux.vars[vars::U1](Idx)  = g(Idx)*TUpDown[dir][1](Idx);
  flux.vars[vars::U1](Idx).eval();
  /* Reads:
   * -----
   * g, TUpDown[dir][1] : 2
   *
   * Writes:
   * ------
   * flux[vars::U1] : 1 */

  flux.vars[vars::U2](Idx)  = g(Idx)*TUpDown[dir][2](Idx);
  flux.vars[vars::U2](Idx).eval();
  /* Reads:
   * -----
   * g, TUpDown[dir][2] : 2
   *
   * Writes:
   * ------
   * flux[vars::U2] : 1 */

  flux.vars[vars::U3](Idx)  = g(Idx)*TUpDown[dir][3](Idx);
  flux.vars[vars::U3](Idx).eval();
  /* Reads:
   * -----
   * g, TUpDown[dir][3] : 2
   *
   * Writes:
   * ------
   * flux[vars::U3] : 1 */

  flux.vars[vars::B1](Idx)  = g(Idx)*(bCon[1](Idx)*uCon[dir](Idx) - bCon[dir](Idx)*uCon[1](Idx));
  flux.vars[vars::B1](Idx).eval();
  /* Reads:
   * -----
   * g, bCon[1], bCon[dir], uCon[1], uCon[dir] : 5
   *
   * Writes:
   * ------
   * flux[vars::B1] : 1 */

  flux.vars[vars::B2](Idx)  = g(Idx)*(bCon[2](Idx)*uCon[dir](Idx) - bCon[dir](Idx)*uCon[2](Idx));
  flux.vars[vars::B2](Idx).eval();
  /* Reads:
   * -----
   * g, bCon[2], bCon[dir], uCon[2], uCon[dir] : 5
   *
   * Writes:
   * ------
   * flux[vars::B2] : 1 */

  flux.vars[vars::B3](Idx) = g(Idx)*(bCon[3](Idx)*uCon[dir](Idx) - bCon[dir](Idx)*uCon[3](Idx));
  flux.vars[vars::B3](Idx).eval();
  /* Reads:
   * -----
   * g, bCon[3], bCon[dir], uCon[3], uCon[dir] : 5
   *
   * Writes:
   * ------
   * flux[vars::B3] : 1 */

  if (params::conduction)
  {
    flux.vars[vars::Q](Idx) = g(Idx)*(uCon[dir](Idx) * qTilde(Idx));
    flux.vars[vars::Q](Idx).eval();
    /* Reads:
     * -----
     * g, uCon[dir], qTilde : 3
     *
     * Writes:
     * ------
     * flux[vars::Q] : 1 */
  }

  if (params::viscosity)
  {
    flux.vars[vars::DP](Idx) = g(Idx)*(uCon[dir](Idx) * deltaPTilde(Idx));
    flux.vars[vars::DP](Idx).eval();
    /* Reads:
     * -----
     * g, uCon[dir], deltaPTilde : 3
     *
     * Writes:
     * ------
     * flux[vars::DP] : 1 */
  }
  /* Total reads : 32
   * Total writes: 10 */ 

  numReads  = 32;
  numWrites = 10;
}

void fluidElement::computeTimeDerivSources(const geometry &geom,
                                           const fluidElement &elemOld,
                                           const fluidElement &elemNew,
                                           const double dt,
                                           grid &sources,
                                           int &numReads,
                                           int &numWrites,
					   const array Idx
                                          )
{
  for (int var=0; var<vars::dof; var++)
  {
    sources.vars[var](Idx) = 0.;
  }
  numReads = 0;
  numWrites = 0;

  if (params::conduction || params::viscosity)
  {
    // Non-ideal pieces. Note that we only compute
    // the terms proportional to time derivatives. 

    // First, compute part of the source terms
    // shared by conduction and viscosity, i.e. 
    // u_{\mu;\nu} and u^{\mu}_{;\mu}
    // We only include the time derivative here!
    for(int mu=0; mu<NDIM; mu++)
    {
      dtuCov[mu](Idx) = (elemNew.uCov[mu](Idx)-elemOld.uCov[mu](Idx))/dt;
    }
    // Compute divergence. We could compute it from the derivatives of uCon,
    //    but we already have u_{\mu;\nu}, so let's make use of it
    // Naturally, this is not truly the divergence. Only the part proportional
    // to dt_u
    divuCov(Idx) = 0.;
    for(int mu=0; mu<NDIM; mu++)
    {
      divuCov(Idx) += geom.gCon[0][mu](Idx)*dtuCov[mu](Idx);
    }
      
    // -------------------------------------
    // Now, look at viscosity-specific terms
    if(params::viscosity)
    {
    	// Compute target deltaP (time deriv part)
      deltaP0(Idx) = divuCov(Idx)*rho(Idx)*nu_emhd(Idx)*(-1.0);	    
    
      for(int mu=0; mu<NDIM; mu++)
      {
        deltaP0(Idx) += 3. * rho(Idx)*nu_emhd(Idx)*bCon[0](Idx)*bCon[mu](Idx)
	  / bSqr(Idx)*dtuCov[mu](Idx);
      }

    	if (params::highOrderTermsViscosity == 1)
      {
	deltaP0(Idx) *= af::sqrt(tau(Idx)/rho(Idx)/nu_emhd(Idx)/temperature(Idx));
      }
	
	    //Note on sign: we put the sources on the LHS when
    	//computing the residual!
    	sources.vars[vars::DP](Idx) = geom.g(Idx)*(deltaP0(Idx))/tau(Idx)*(-1.0);

    	if (params::highOrderTermsViscosity == 1)
      {
	sources.vars[vars::DP](Idx) -= 0.5*geom.g(Idx)*divuCov(Idx)*deltaPTilde(Idx);
      }
	sources.vars[vars::DP](Idx).eval();
      /* Reads:
       * -----
       *  dtuCov[mu](elemNew.uCov[mu] : 4, elemOld.uCov[mu] :4) : 8
       *  divuCov(geom.gCon[0][mu] : 4, dtuCov[mu] : 0)         : 4
       *  deltaP0(divuCov : 0 (already accounted),
       *          rho : 1,
       *          nu_emhd(rho:0, u:1) : 1
       *          bCon[mu] : 4,
       *          bSqr : 1,
       *          dtuCov[mu] : 0 (already accounted),
       *          temperature(rho:0, u:0) : 0,
       *         )                                              : 7
       *  geom.g                                                : 1
       *  tau                                                   : 1
       *  deltaPTilde                                           : 1
       *
       * Writes:
       * ------
       * sources[vars::DP] : 1 */

      numReads = 21;
      if (params::highOrderTermsViscosity)
      {
        numReads += 1;
      }
      numWrites = 1;

    } /* End of viscosity specific terms */
    
    // -------------------------------------
    // Finally, look at conduction-specific terms (time deriv terms)
    if(params::conduction)
    {
      q0(Idx) = rho(Idx)*chi_emhd(Idx)*bCon[0](Idx)*(-1.0)
	/ bNorm(Idx) *(elemNew.temperature(Idx) - elemOld.temperature(Idx))/dt;
    
      for(int nu=0;nu<NDIM;nu++)
      {
	q0(Idx) -= rho(Idx)*chi_emhd(Idx)*temperature(Idx)*
	  bCon[nu](Idx)/bNorm(Idx)*uCon[0](Idx)*dtuCov[nu](Idx);
      }
	
    	if (params::highOrderTermsConduction == 1)
      {
	q0(Idx) *= af::sqrt(tau(Idx)/rho(Idx)/chi_emhd(Idx))/temperature(Idx);
      }
	    
	    //Note on sign: we put the sources on the LHS when
    	//computing the residual!
    	sources.vars[vars::Q](Idx) = geom.g(Idx)*(q0(Idx))/tau(Idx)*(-1.0);
	
    	if (params::highOrderTermsConduction == 1)
      {
	sources.vars[vars::Q](Idx) -= 0.5*geom.g(Idx)*divuCov(Idx)*qTilde(Idx);
      }
	sources.vars[vars::Q](Idx).eval();
      /* Reads:
       * -----
       *  q0(rho                      : 1,
       *     chi_emhd                 : 1
       *     temperature(rho:0, u:0)  : 0,
       *     bCon[mu]                 : 4,
       *     bNorm                    : 1,
       *     dtuCov[mu]               : 8 ,
       *     uCon0                    : 1
       *    )                                           : 16
       *  geom.g                                        : 1
       *  tau                                           : 1
       *  qTilde                                        : 1
       *  divuCov(geom.gCon[0][mu] : 4, dtuCov[mu] : 0) : 4
       *
       * Writes:
       * ------
       * sources[vars::Q] : 1 */

      numReads += 18;
      if (params::highOrderTermsConduction)
      {
        numReads += 5;
      }
      numWrites += 1;

    } /* End of conduction */

  } /* End of EMHD: viscosity || conduction */
}


void fluidElement::computeImplicitSources(const geometry &geom,
					  grid &sources,
					  array &tauDamp,
                                          int &numReads,
                                          int &numWrites,
					  const array Idx
                                         )
{
  for (int var=0; var<vars::dof; var++)
  {
    sources.vars[var](Idx) = 0.;
  }

  numReads = 0;
  numWrites = 0;

  // Non-ideal pieces. Note that we only compute
  // the terms treated implicitly. 
  // Look at viscosity-specific terms
  if(params::viscosity)
  {
    //Note on sign: we put the sources on the LHS when
    //computing the residual!
    sources.vars[vars::DP](Idx) = geom.g(Idx)*(deltaPTilde(Idx))/tauDamp(Idx);
    sources.vars[vars::DP](Idx).eval();
    /* Reads:
     * -----
     *  geom.g          : 1
     *  deltaPTilde     : 1
     *  tau             : 1
     *
     * Writes:
     * ------
     * sources[vars::DP] : 1 */
    numReads  += 3;
    numWrites += 1;

  } /* End of viscosity specific terms */
  
  // -------------------------------------
  // Finally, look at conduction-specific terms (implicit terms)
  if(params::conduction)
  {
    //Note on sign: we put the sources on the LHS when
    //computing the residual!
    sources.vars[vars::Q](Idx) = geom.g(Idx)*(qTilde(Idx))/tauDamp(Idx);
    sources.vars[vars::Q](Idx).eval();
    /* Reads:
     * -----
     *  geom.g          : 1
     *  qTilde          : 1
     *  tau             : 1
     *
     * Writes:
     * ------
     * sources[vars::DP] : 1 */
    numReads  += 3;
    numWrites += 1;

  } /* End of conduction */
  
}

void fluidElement::computeExplicitSources(const geometry &geom,
					  grid &sources,
                                          int &numReads,
                                          int &numWrites,
					  const array Idx
					  )
{
  for (int var=0; var<vars::dof; var++)
  {
    sources.vars[var](Idx) = 0.;
  }

  //Note on sign: residual computation places
  // the source terms on the LHS of the equation!
  // All ideal MHD terms are treated explicitly.
  numReads = 0; numWrites = 0;
  if (params::metric != metrics::MINKOWSKI)
  {
    for (int nu=0; nu<NDIM; nu++)
    {
      for (int kappa=0; kappa<NDIM; kappa++)
      {
        for (int lamda=0; lamda<NDIM; lamda++)
        {
          sources.vars[vars::U + nu](Idx) -=
            geom.g(Idx)
	    * TUpDown[kappa][lamda](Idx)
	    * geom.gammaUpDownDown[lamda][kappa][nu](Idx);
        }
      }
      sources.vars[vars::U + nu](Idx).eval();
      /* Reads:
       * -----
       *  geom.g                                         : 1
       *  TUpDown[kappa 0-3][lambda 0-3]                 : 16
       *  geom.gammaUpDownDown[lamda 0-3][kappa 0-3][nu] : 16
       *
       * Writes:
       * ------
       * sources[vars::U + nu] : 1 */

      numReads  += 33;
      numWrites += 1;
    }
  }

  if (params::conduction || params::viscosity)
  {
    // Non-ideal pieces (explicit parts)
    // First, compute part of the source terms
    // shared by conduction and viscosity, i.e. 
    // u_{\mu;\nu} and u^{\mu}_{;\mu}
    // Note that the derivatives are all precomputed 
    // in computeEMHDGradients

    // Compute divergence. We could compute it from the derivatives of uCon,
    //    but we already have u_{\mu;\nu}, so let's make use of it
    // Note that this does NOT include terms proportional to dt_u
    divuCov(Idx) = 0.;
    for(int mu=0;mu<NDIM;mu++)
    {  
      for(int nu=0;nu<NDIM;nu++)
      {
	divuCov(Idx) += geom.gCon[mu][nu](Idx)*graduCov[mu][nu](Idx);
      }
    }
      
    // -------------------------------------
    // Now, look at viscosity-specific terms
    if(params::viscosity)
    {
	    // Compute target deltaP (explicit part)
      deltaP0(Idx) = divuCov(Idx)*rho(Idx)*nu_emhd(Idx)*(-1.0);
	
	    for(int mu=0;mu<NDIM;mu++)
      {
	      for(int nu=0;nu<NDIM;nu++)
        {
	  deltaP0(Idx) += 3. * rho(Idx) * nu_emhd(Idx) 
	    * bCon[mu](Idx) * bCon[nu](Idx) / bSqr(Idx)
	    * graduCov[mu][nu](Idx);
         }
      }

      array deltaP0Tilde = deltaP0;
      if (params::highOrderTermsViscosity == 1)
	{
	  deltaP0Tilde(Idx) = deltaP0(Idx) * af::sqrt(tau(Idx)/rho(Idx)/nu_emhd(Idx)/temperature(Idx));
	}
	
	    //Note on sign: we put the sources on the LHS when
    	//computing the residual!
    	//The damping term proportional to deltaPTilde is in the implicit sector.
      sources.vars[vars::DP](Idx) = geom.g(Idx)*(deltaP0Tilde(Idx))/tau(Idx)*(-1.0);
    
      if (params::highOrderTermsViscosity == 1)
	    {
	      sources.vars[vars::DP](Idx) -= 0.5*geom.g(Idx)*divuCov(Idx)*deltaPTilde(Idx);
	    }
      sources.vars[vars::DP](Idx).eval();

      /* Reads:
       * -----
       *  divuCov(geom.gCon[mu][nu] : 16, graduCov[mu][nu] : 16): 32
       *  deltaP0(divuCov : 0 (already accounted),
       *          rho : 1,
       *          nu_emhd(rho:0, u:1) : 1
       *          bCon[mu] : 4,
       *          bCon[nu] : 0,(already accounted)
       *          bSqr : 1,
       *          graduCov[mu][nu] : 0 (already accounted),
       *          temperature(rho:0, u:1) : 1,
       *         )                                              : 7
       *  geom.g                                                : 1
       *  tau                                                   : 1
       *  deltaPTilde                                           : 1
       *
       * Writes:
       * ------
       * sources[vars::DP] : 1 */
      numReads += 41;
      if (params::highOrderTermsViscosity)
      {
        numReads += 1;
      }
      numWrites += 1;

    } /* End of viscosity specific terms */
    
    // -------------------------------------
    // Finally, look at conduction-specific terms (explicit part)
    if(params::conduction)
    {
      q0(Idx) = 0.;
    	//q0 is not exactly targetQ, as the time derivative parts
    	// are in the implicit sector
    	for(int mu=0;mu<NDIM;mu++)
      {
	q0(Idx) -= rho(Idx)*chi_emhd(Idx)*bCon[mu](Idx)/bNorm(Idx)*gradT[mu](Idx);		
	  
        for(int nu=0;nu<NDIM;nu++)
        {
	  q0(Idx) -=  rho(Idx)*chi_emhd(Idx)*temperature(Idx)*bCon[nu](Idx)/bNorm(Idx)
	    * uCon[mu](Idx)*graduCov[mu][nu](Idx);
        }
      }
    
      array q0Tilde = q0;
      if (params::highOrderTermsConduction == 1)
      {
	q0Tilde(Idx) = q0(Idx) * af::sqrt(  tau(Idx)/rho(Idx)/chi_emhd(Idx))/temperature(Idx);
      }

      //Note on sign: we put the sources.vars on the LHS when
    	//computing the residual!
    	// The damping term proportional to qTilde is in the implicit sector. 
      sources.vars[vars::Q](Idx) = geom.g(Idx)*(q0Tilde(Idx))/tau(Idx)*(-1.0);

    	if (params::highOrderTermsConduction == 1)
      {
	sources.vars[vars::Q](Idx) -= 0.5*geom.g(Idx)*divuCov(Idx)*qTilde(Idx);
      }
	sources.vars[vars::Q](Idx).eval();
      /* Reads:
       * -----
       *  q0(rho                      : 1,
       *     chi_emhd                 : 1
       *     temperature(rho:0, u:0)  : 0,
       *     bCon[mu]                 : 4,
       *     bNorm                    : 1,
       *     gradT[mu]                : 4 ,
       *     uCon[mu]                 : 4 ,
       *     graduCov[mu][nu]         : 16
       *    )                                 : 31
       *  geom.g                              : 1
       *  tau                                 : 1
       *  qTilde                              : 1
       *  divuCov(gCon[mu][nu]     : 16,
       *          graduCov[mu][nu] : 0
       *        )                             : 16
       *
       * Writes:
       * ------
       * sources[vars::Q] : 1 */
      numReads += 33;
      if (params::highOrderTermsConduction)
      {
        numReads += 17;
      }
      numWrites += 1;

    } /* End of conduction */
  } /* End of EMHD: viscosity || conduction */

}

void fluidElement::computeEMHDGradients(const geometry &geom,
                                        const double dX[3],
                                        int &numReads,
                                        int &numWrites
                                       )
{
  double dX1 = dX[directions::X1];
  double dX2 = dX[directions::X2];
  double dX3 = dX[directions::X3];
  
  numReads = 0, numWrites = 0;
  for(int mu=0;mu<NDIM;mu++)
  {
    //Time derivative needs to be reset for reach residual computation,
    // so not computed here.
    graduCov[0][mu] = 0.;
  
    int numReadsTmp, numWritesTmp;
    graduCov[1][mu] = reconstruction::slope(directions::X1,dX1,uCov[mu],
                                            numReadsTmp, numWritesTmp
                                           );
    numReads  += numReadsTmp;
    numWrites += numWritesTmp;

    graduCov[2][mu] = 0.;
    if(params::dim>1)
    {
      graduCov[2][mu] = reconstruction::slope(directions::X2,dX2,uCov[mu],
                                              numReadsTmp, numWritesTmp
                                             );
      numReads  += numReadsTmp;
      numWrites += numWritesTmp;
    }
  
    graduCov[3][mu] = 0.;
    if(params::dim>2)
    {
      graduCov[3][mu] = reconstruction::slope(directions::X3,dX3,uCov[mu],
                                              numReadsTmp, numWritesTmp
                                             );
      numReads  += numReadsTmp;
      numWrites += numWritesTmp;
    }

    for(int nu=0;nu<NDIM;nu++)
    {  
	    for(int lambda=0;lambda<NDIM;lambda++)
	    {
	      graduCov[nu][mu] -= geom.gammaUpDownDown[lambda][nu][mu]*uCov[lambda];
	    }

      graduCov[nu][mu].eval();
      /* Reads:
       * -----
       *  geom.gammaUpDownDown[lamda 0-3][mu][nu] : 4
       *  uCon[lambda 0-3]                        : 4
       *
       * Writes:
       * ------
       * graduCov[nu][mu] : 1 */
      numReads  += 8;
      numWrites += 1;
    }
  }
  
  if(params::conduction)
  {
    /* Time derivative not computed here */
    gradT[0] = 0.;

    int numReadsTmp, numWritesTmp;
    gradT[1] = reconstruction::slope(directions::X1,dX1,temperature,
                                     numReadsTmp, numWritesTmp
                                    );
    numReads  += numReadsTmp;
    numWrites += numWritesTmp;

    gradT[2] = 0.;
    if(params::dim>1)
    {
      gradT[2] = reconstruction::slope(directions::X2,dX2,temperature,
                                       numReadsTmp, numWritesTmp
                                      );
      numReads  += numReadsTmp;
      numWrites += numWritesTmp;
    }  

    gradT[3] = 0.;
    if(params::dim>2)
    {
      gradT[3] = reconstruction::slope(directions::X3,dX3,temperature,
                                       numReadsTmp, numWritesTmp
                                      );
      numReads  += numReadsTmp;
      numWrites += numWritesTmp;
    }
  } /* End of conduction specific terms */

}
