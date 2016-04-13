#include "physics.hpp"

fluidElement::fluidElement(const grid &prim,
                           const grid &magneticFields,
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
  gammaLorentzFactor = zero;

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
    q0      = zero;
  }

  set(prim, magneticFields, geom, numReads, numWrites);
}

void fluidElement::set(const grid &prim,
                       const grid &magneticFields,
                       const geometry &geom,
                       int &numReads,
                       int &numWrites
                      )
{
  rho = af::max(prim.vars[vars::RHO], params::rhoFloorInFluidElement);
  u   = af::max(prim.vars[vars::U  ], params::uFloorInFluidElement);
  u1  = prim.vars[vars::U1];
  u2  = prim.vars[vars::U2];
  u3  = prim.vars[vars::U3];
  B1  = magneticFields.vars[vars::B1];
  B2  = magneticFields.vars[vars::B2];
  B3  = magneticFields.vars[vars::B3];

  pressure    = (params::adiabaticIndex - 1.)*u;
  temperature = af::max(pressure/rho,params::temperatureFloorInFluidElement);

  soundSpeed  = af::sqrt( params::adiabaticIndex*pressure
                         /(rho+params::adiabaticIndex*u)
                        );
  
  gammaLorentzFactor =
    af::sqrt(1 + geom.gCov[1][1] * u1 * u1
               + geom.gCov[2][2] * u2 * u2
               + geom.gCov[3][3] * u3 * u3

             + 2*(  geom.gCov[1][2] * u1 * u2
                  + geom.gCov[1][3] * u1 * u3
                  + geom.gCov[2][3] * u2 * u3
                 )
            );
  gammaLorentzFactor.eval(); 
  /* Reads:
   * -----
   * gCov[1][1], gCov[2][2], gCov[3][3], gCov[1][2], gCov[1][3], gCov[2][3]: 6
   * u1, u2, u3 : 3
   *
   * Writes:
   * ------
   * gammaLorentzFactor : 1 */

  uCon[0] = gammaLorentzFactor/geom.alpha;
  uCon[0].eval(); 
  /* Reads:
   * -----
   * gammaLorentzFactor, alpha : 2
   *
   * Writes:
   * ------
   * uCon[0] : 1 */


  uCon[1] = u1 - gammaLorentzFactor*geom.gCon[0][1]*geom.alpha;
  uCon[1].eval(); 
  /* Reads:
   * -----
   * u1, gammaLorentzFactor, gCon[0][1], alpha : 4
   *
   * Writes:
   * ------
   * uCon[1] : 1 */

  uCon[2] = u2 - gammaLorentzFactor*geom.gCon[0][2]*geom.alpha;
  uCon[2].eval(); 
  /* Reads:
   * -----
   * u2, gammaLorentzFactor, gCon[0][2], alpha : 4
   *
   * Writes:
   * ------
   * uCon[2] : 1 */

  uCon[3] = u3 - gammaLorentzFactor*geom.gCon[0][3]*geom.alpha;
  uCon[3].eval();
  /* Reads:
   * -----
   * u3, gammaLorentzFactor, gCon[0][3], alpha : 4
   *
   * Writes:
   * ------
   * uCon[3] : 1 */

  for (int mu=0; mu < NDIM; mu++)
  {
    uCov[mu] =  geom.gCov[mu][0] * uCon[0]
              + geom.gCov[mu][1] * uCon[1]
              + geom.gCov[mu][2] * uCon[2]
              + geom.gCov[mu][3] * uCon[3];
    uCov[mu].eval(); 
  } 
  /* Reads:
   * -----
   * gCov[mu][0], gCov[mu][1], gCov[mu][2], gCov[mu][3]: 16
   * uCon[0], uCon[1], uCon[2], uCon[3]: 4 x 4 = 16
   *
   * Writes:
   * ------
   * uCov[mu] : 4 */

  bCon[0] =  B1*uCov[1] + B2*uCov[2] + B3*uCov[3];
  bCon[0].eval();
  /* Reads:
   * -----
   * B1, B2, B3, uCov[1], uCov[2], uCov[3] : 6
   *
   * Writes:
   * ------
   * bCon[0] : 1 */

  bCon[1] = (B1 + bCon[0] * uCon[1])/uCon[0];
  bCon[1].eval();
  /* Reads:
   * -----
   * B1, bCon[0], uCon[1], uCon[0]: 4
   *
   * Writes:
   * ------
   * bCon[1] : 1 */

  bCon[2] = (B2 + bCon[0] * uCon[2])/uCon[0];
  bCon[2].eval();
  /* Reads:
   * -----
   * B2, bCon[0], uCon[2], uCon[0]: 4
   *
   * Writes:
   * ------
   * bCon[2] : 1 */

  bCon[3] = (B3 + bCon[0] * uCon[3])/uCon[0];
  bCon[3].eval();
  /* Reads:
   * -----
   * B3, bCon[0], uCon[3], uCon[0]: 4
   *
   * Writes:
   * ------
   * bCon[3] : 1 */

  for (int mu=0; mu < NDIM; mu++)
  {
    bCov[mu] =  geom.gCov[mu][0] * bCon[0]
              + geom.gCov[mu][1] * bCon[1]
              + geom.gCov[mu][2] * bCon[2]
              + geom.gCov[mu][3] * bCon[3];
    bCov[mu].eval();
  }
  /* Reads:
   * -----
   * gCov[mu][0], gCov[mu][1], gCov[mu][2], gCov[mu][3]: 16
   * bCon[0], bCon[1], bCon[2], bCon[3]: 4 x 4 = 16
   *
   * Writes:
   * ------
   * bCov[mu] : 4 */

  bSqr =  bCon[0]*bCov[0] + bCon[1]*bCov[1]
        + bCon[2]*bCov[2] + bCon[3]*bCov[3] + params::bSqrFloorInFluidElement;
  bSqr.eval();
  /* Reads:
   * -----
   * bCon[0], bCon[1], bCon[2], bCon[3]: 4
   * bCov[0], bCov[1], bCov[2], bCov[3]: 4
   *
   * Writes:
   * ------
   * bSqr : 1 */

  bNorm = af::sqrt(bSqr);

  // Note: this needs to be before setFluidElementParameters
  // because the closure relation uses q, deltaP!
  if (params::conduction==1)
  {
    qTilde = prim.vars[vars::Q];

    if (params::highOrderTermsConduction==1)
    {
      q = qTilde * temperature * af::sqrt(rho*params::ConductionAlpha*soundSpeed*soundSpeed);
      q.eval();
    }
    else
    {
      q = qTilde;
    }
  }

  if (params::viscosity==1)
  {
    deltaPTilde = prim.vars[vars::DP];

    if (params::highOrderTermsViscosity == 1)
    {
      deltaP = deltaPTilde * af::sqrt(temperature * rho * params::ViscosityAlpha*soundSpeed*soundSpeed);
      deltaP.eval();
    }
    else
    {
      deltaP = deltaPTilde;
    }
  }

  // Note: this uses q, deltaP, bSqr!
  setFluidElementParameters(geom);

  for (int mu=0; mu < NDIM; mu++)
  {
    NUp[mu] = rho * uCon[mu];

    for (int nu=0; nu < NDIM; nu++)
    {
      TUpDown[mu][nu] =   (rho + u + pressure + bSqr)*uCon[mu]*uCov[nu]
                        + (pressure + 0.5*bSqr)*DELTA(mu, nu)
                        - bCon[mu] * bCov[nu];

      if (params::conduction==1)
      {
        TUpDown[mu][nu] += q/bNorm * (uCon[mu]*bCov[nu] + bCon[mu]*uCov[nu]);
      }

      if (params::viscosity==1)
      {
        TUpDown[mu][nu] += (- deltaP)       
                           * (  bCon[mu] * bCov[nu]/bSqr
                              - (1./3.)*(DELTA(mu, nu) + uCon[mu]*uCov[nu])
                             );
      }
      TUpDown[mu][nu].eval();
      /* Reads:
       * -----
       * rho, u, bSqr, q, deltaP: 5 x 16 = 80
       * uCon[mu], uCov[nu], bCon[mu], bCov[nu]: 4 x 16 = 64
       *
       * Writes:
       * ------
       * TUpDown[mu][nu] : 16 */
    }
    NUp[mu].eval();
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

void fluidElement::computeMagneticFluxes(const geometry &geom, 
                                         const int dir,
                                         grid &flux,
                                         int &numReads,
                                         int &numWrites
                                        )
{
  array g = geom.g;

  flux.vars[vars::B1]  = g*(bCon[1]*uCon[dir] - bCon[dir]*uCon[1]);
  flux.vars[vars::B1].eval();
  /* Reads:
   * -----
   * g, bCon[1], bCon[dir], uCon[1], uCon[dir] : 5
   *
   * Writes:
   * ------
   * flux[vars::B1] : 1 */
  numReads  = 5;
  numWrites = 1;

  flux.vars[vars::B2]  = g*(bCon[2]*uCon[dir] - bCon[dir]*uCon[2]);
  flux.vars[vars::B2].eval();
  /* Reads:
   * -----
   * g, bCon[2], bCon[dir], uCon[2], uCon[dir] : 5
   *
   * Writes:
   * ------
   * flux[vars::B2] : 1 */
  numReads  += 5;
  numWrites += 1;

  flux.vars[vars::B3]  = g*(bCon[3]*uCon[dir] - bCon[dir]*uCon[3]);
  flux.vars[vars::B3].eval();
  /* Reads:
   * -----
   * g, bCon[3], bCon[dir], uCon[3], uCon[dir] : 5
   *
   * Writes:
   * ------
   * flux[vars::B3] : 1 */
  numReads  += 5;
  numWrites += 1;
}

void fluidElement::computeFluidFluxes(const geometry &geom, 
                                      const int dir,
                                      grid &flux,
                                      int &numReads,
                                      int &numWrites
                                      )
{
  array g = geom.g;

  flux.vars[vars::RHO] = g*NUp[dir];
  flux.vars[vars::RHO].eval();
  /* Reads:
   * -----
   * g, NUp[dir] : 2
   *
   * Writes:
   * ------
   * flux[vars::RHO] : 1 */

  flux.vars[vars::U]   = g*TUpDown[dir][0] + flux.vars[vars::RHO];
  flux.vars[vars::U].eval();
  /* Reads:
   * -----
   * g, TUpDown[dir][0], flux[vars::RHO] : 3
   *
   * Writes:
   * ------
   * flux[vars::U] : 1 */

  flux.vars[vars::U1]  = g*TUpDown[dir][1];
  flux.vars[vars::U1].eval();
  /* Reads:
   * -----
   * g, TUpDown[dir][1] : 2
   *
   * Writes:
   * ------
   * flux[vars::U1] : 1 */

  flux.vars[vars::U2]  = g*TUpDown[dir][2];
  flux.vars[vars::U2].eval();
  /* Reads:
   * -----
   * g, TUpDown[dir][2] : 2
   *
   * Writes:
   * ------
   * flux[vars::U2] : 1 */

  flux.vars[vars::U3]  = g*TUpDown[dir][3];
  flux.vars[vars::U3].eval();
  /* Reads:
   * -----
   * g, TUpDown[dir][3] : 2
   *
   * Writes:
   * ------
   * flux[vars::U3] : 1 */

  if (params::conduction)
  {
    flux.vars[vars::Q] = g*(uCon[dir] * qTilde);
    flux.vars[vars::Q].eval();
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
    flux.vars[vars::DP] = g*(uCon[dir] * deltaPTilde);
    flux.vars[vars::DP].eval();
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
                                           int &numWrites
                                          )
{
  for (int var=0; var<vars::dof; var++)
  {
    sources.vars[var] = 0.;
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
      dtuCov[mu] = (elemNew.uCov[mu]-elemOld.uCov[mu])/dt;
    }
    // Compute divergence. We could compute it from the derivatives of uCon,
    //    but we already have u_{\mu;\nu}, so let's make use of it
    // Naturally, this is not truly the divergence. Only the part proportional
    // to dt_u
    divuCov = 0.;
    for(int mu=0; mu<NDIM; mu++)
    {
      divuCov += geom.gCon[0][mu]*dtuCov[mu];
    }
      
    // -------------------------------------
    // Now, look at viscosity-specific terms
    if(params::viscosity)
    {
    	// Compute target deltaP (time deriv part)
    	deltaP0 = -divuCov*rho*nu_emhd;	    
    
      for(int mu=0; mu<NDIM; mu++)
      {
        deltaP0 += 3. * rho*nu_emhd*bCon[0]*bCon[mu]
	                / bSqr*dtuCov[mu];
      }

    	if (params::highOrderTermsViscosity == 1)
      {
	      deltaP0 *= af::sqrt(tau/rho/nu_emhd/temperature);
      }
	
	    //Note on sign: we put the sources on the LHS when
    	//computing the residual!
    	sources.vars[vars::DP] = -geom.g*(deltaP0)/tau;

    	if (params::highOrderTermsViscosity == 1)
      {
	      sources.vars[vars::DP] -= 0.5*geom.g*divuCov*deltaPTilde;
      }
      sources.vars[vars::DP].eval();
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
    	q0 =  - rho*chi_emhd*bCon[0]
	          / bNorm *(elemNew.temperature - elemOld.temperature)/dt;
    
      for(int nu=0;nu<NDIM;nu++)
      {
	      q0 -= rho*chi_emhd*temperature*
	            bCon[nu]/bNorm*uCon[0]*dtuCov[nu];
      }
	
    	if (params::highOrderTermsConduction == 1)
      {
	      q0 *= af::sqrt(tau/rho/chi_emhd)/temperature;
      }
	    
	    //Note on sign: we put the sources on the LHS when
    	//computing the residual!
    	sources.vars[vars::Q] = -geom.g*(q0)/tau;
	
    	if (params::highOrderTermsConduction == 1)
      {
	      sources.vars[vars::Q] -= 0.5*geom.g*divuCov*qTilde;
      }
      sources.vars[vars::Q].eval();
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
                                          int &numWrites
                                         )
{
  for (int var=0; var<vars::dof; var++)
  {
    sources.vars[var] = 0.;
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
    sources.vars[vars::DP] = geom.g*(deltaPTilde)/tauDamp;
    sources.vars[vars::DP].eval();
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
    sources.vars[vars::Q] = geom.g*(qTilde)/tauDamp;
    sources.vars[vars::Q].eval();
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
                                          int &numWrites
                             					   )
{
  for (int var=0; var<vars::dof; var++)
  {
    sources.vars[var] = 0.;
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
          sources.vars[vars::U + nu] -=
            geom.g
          * TUpDown[kappa][lamda]
          * geom.gammaUpDownDown[lamda][kappa][nu];
        }
      }
      sources.vars[vars::U + nu].eval();
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
    divuCov = 0.;
    for(int mu=0;mu<NDIM;mu++)
    {  
      for(int nu=0;nu<NDIM;nu++)
      {
	      divuCov += geom.gCon[mu][nu]*graduCov[mu][nu];
      }
    }
      
    // -------------------------------------
    // Now, look at viscosity-specific terms
    if(params::viscosity)
    {
	    // Compute target deltaP (explicit part)
    	deltaP0 = -divuCov*rho*nu_emhd;
	
	    for(int mu=0;mu<NDIM;mu++)
      {
	      for(int nu=0;nu<NDIM;nu++)
        {
    	    deltaP0 += 3. * rho * nu_emhd 
                        * bCon[mu] * bCon[nu] / bSqr
	                      * graduCov[mu][nu];
         }
      }

      array deltaP0Tilde;
    	if (params::highOrderTermsViscosity == 1)
      {
	      deltaP0Tilde = deltaP0 * af::sqrt(tau/rho/nu_emhd/temperature);
      }
      else
      {
        deltaP0Tilde = deltaP0;
      }
	
	    //Note on sign: we put the sources on the LHS when
    	//computing the residual!
    	//The damping term proportional to deltaPTilde is in the implicit sector.
    	sources.vars[vars::DP] = -geom.g*(deltaP0Tilde)/tau;
    
      if (params::highOrderTermsViscosity == 1)
	    {
	      sources.vars[vars::DP] -= 0.5*geom.g*divuCov*deltaPTilde;
	    }
      sources.vars[vars::DP].eval();

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
	    q0 = 0.;
    	//q0 is not exactly targetQ, as the time derivative parts
    	// are in the implicit sector
    	for(int mu=0;mu<NDIM;mu++)
      {
	      q0 -= rho*chi_emhd*bCon[mu]/bNorm*gradT[mu];		
	  
        for(int nu=0;nu<NDIM;nu++)
        {
	        q0 -=  rho*chi_emhd*temperature*bCon[nu]/bNorm
               * uCon[mu]*graduCov[mu][nu];
        }
	    }
    
      array q0Tilde;
      if (params::highOrderTermsConduction == 1)
      {
	      q0Tilde = q0 * af::sqrt(  tau/rho/chi_emhd)/temperature;
      }
      else
      {
        q0Tilde = q0;
      }

      //Note on sign: we put the sources.vars on the LHS when
    	//computing the residual!
    	// The damping term proportional to qTilde is in the implicit sector. 
      sources.vars[vars::Q] = -geom.g*(q0Tilde)/tau;

    	if (params::highOrderTermsConduction == 1)
      {
    	  sources.vars[vars::Q] -= 0.5*geom.g*divuCov*qTilde;
      }
      sources.vars[vars::Q].eval();
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

    if(params::dim>1)
    {
      graduCov[2][mu] = reconstruction::slope(directions::X2,dX2,uCov[mu],
                                              numReadsTmp, numWritesTmp
                                             );
      numReads  += numReadsTmp;
      numWrites += numWritesTmp;
    }
  
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

    if(params::dim>1)
    {
      gradT[2] = reconstruction::slope(directions::X2,dX2,temperature,
                                       numReadsTmp, numWritesTmp
                                      );
      numReads  += numReadsTmp;
      numWrites += numWritesTmp;
    }  

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
