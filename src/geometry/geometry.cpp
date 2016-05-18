#include "geometry.hpp"
#include "CoordinateChangeFunctionsArray.hpp"

geometry::geometry(const int metric,
                   const double blackHoleSpin,
                   const double hSlope,
                   const coordinatesGrid &XCoordsGrid
                  )
{
  PetscPrintf(PETSC_COMM_WORLD, "Setup metric\n");   
  GAMMA_EPS = 1.e-5;
  N1        = XCoordsGrid.N1;
  N2        = XCoordsGrid.N2;
  N3        = XCoordsGrid.N3;
  dim       = XCoordsGrid.dim;
  numGhost  = XCoordsGrid.numGhost;

  this->metric = metric;
  this->blackHoleSpin = blackHoleSpin;
  this->hSlope = hSlope;

  XCoords[directions::X1] = XCoordsGrid.vars[directions::X1];
  XCoords[directions::X2] = XCoordsGrid.vars[directions::X2];
  XCoords[directions::X3] = XCoordsGrid.vars[directions::X3];

  XCoordsToxCoords(XCoords,xCoords);

  /* Allocate space */
  zero       = 0.*XCoords[0];
  g          = zero;
  array gDet = zero;
  for (int mu=0; mu<NDIM; mu++)
  {
    for (int nu=0; nu<NDIM; nu++)
    {
      gCov[mu][nu] = zero;
      gCon[mu][nu] = zero;
    }
  }
  
  setgCovInXCoords(XCoords, gCov);
  setgDetAndgConFromgCov(gCov, gDet, gCon);

  g = af::sqrt(-gDet);
  alpha = 1./af::sqrt(-gCon[0][0]);
  g.eval();
  alpha.eval();

  af::sync();
  PetscPrintf(PETSC_COMM_WORLD, "Done initializing geometry!\n");   
}

void geometry::computeConnectionCoeffs()
{
  PetscPrintf(PETSC_COMM_WORLD, "Compute connections\n");
  array gammaDownDownDown[NDIM][NDIM][NDIM];
  array gCovPlus[NDIM-1][NDIM][NDIM];
  array gCovMinus[NDIM-1][NDIM][NDIM];
  {
    array XCoordsPlus[NDIM-1];
    array XCoordsMinus[NDIM-1];
    for (int d = 0; d<NDIM-1; d++)
      {
	for(int dd=0;dd<NDIM-1;dd++)
	  {
	    XCoordsPlus[dd] = this->XCoords[dd];
	    XCoordsMinus[dd] = this->XCoords[dd];
	  }
	XCoordsPlus[d]+=GAMMA_EPS;
	XCoordsMinus[d]-=GAMMA_EPS;
	setgCovInXCoords(XCoordsPlus,gCovPlus[d]);
	setgCovInXCoords(XCoordsMinus,gCovMinus[d]);
      }
  }

  for (int mu=0; mu<NDIM; mu++)
  {
    for (int nu=0; nu<NDIM; nu++)
    {
      for (int lamda = 0; lamda<NDIM; lamda++)
      {
      	gammaDownDownDown[mu][nu][lamda] = zero;
        computeGammaDownDownDown(mu,nu,lamda,
                                 gammaDownDownDown[mu][nu][lamda],
				 gCovPlus,gCovMinus
                                );
      }
    }
  }

  for (int mu=0; mu<NDIM; mu++)
  {
    for (int nu=0; nu<NDIM; nu++)
    {
      for (int lamda = 0; lamda<NDIM; lamda++)
      {
        gammaUpDownDown[mu][nu][lamda] = zero;
      
        for(int eta=0; eta<NDIM; eta++)
        {
          gammaUpDownDown[mu][nu][lamda] += 
             gCon[mu][eta]
           * gammaDownDownDown[eta][nu][lamda];
        }

        gammaUpDownDown[mu][nu][lamda].eval();
      }
    }
  }

  af::sync();
  PetscPrintf(PETSC_COMM_WORLD, "Done with connections\n");
}

void geometry::setgDetAndgConFromgCov(const array gCov[NDIM][NDIM],
                                      array &gDet,
                                      array gCon[NDIM][NDIM]
                                     )
{
  gDet = 
        gCov[0][0]*gCov[1][1]*gCov[2][2]*gCov[3][3] 
      + gCov[0][0]*gCov[1][2]*gCov[2][3]*gCov[3][1]
      + gCov[0][0]*gCov[1][3]*gCov[2][1]*gCov[3][2] 
      + gCov[0][1]*gCov[1][0]*gCov[2][3]*gCov[3][2] 
      + gCov[0][1]*gCov[1][2]*gCov[2][0]*gCov[3][3] 
      + gCov[0][1]*gCov[1][3]*gCov[2][2]*gCov[3][0] 
      + gCov[0][2]*gCov[1][0]*gCov[2][1]*gCov[3][3] 
      + gCov[0][2]*gCov[1][1]*gCov[2][3]*gCov[3][0] 
      + gCov[0][2]*gCov[1][3]*gCov[2][0]*gCov[3][1] 
      + gCov[0][3]*gCov[1][0]*gCov[2][2]*gCov[3][1]
      + gCov[0][3]*gCov[1][1]*gCov[2][0]*gCov[3][2] 
      + gCov[0][3]*gCov[1][2]*gCov[2][1]*gCov[3][0]
      - gCov[0][0]*gCov[1][1]*gCov[2][3]*gCov[3][2]
      - gCov[0][0]*gCov[1][2]*gCov[2][1]*gCov[3][3]
      - gCov[0][0]*gCov[1][3]*gCov[2][2]*gCov[3][1]
      - gCov[0][1]*gCov[1][0]*gCov[2][2]*gCov[3][3]
      - gCov[0][1]*gCov[1][2]*gCov[2][3]*gCov[3][0]
      - gCov[0][1]*gCov[1][3]*gCov[2][0]*gCov[3][2]
      - gCov[0][2]*gCov[1][0]*gCov[2][3]*gCov[3][1]
      - gCov[0][2]*gCov[1][1]*gCov[2][0]*gCov[3][3]
      - gCov[0][2]*gCov[1][3]*gCov[2][1]*gCov[3][0]
      - gCov[0][3]*gCov[1][0]*gCov[2][1]*gCov[3][2]
      - gCov[0][3]*gCov[1][1]*gCov[2][2]*gCov[3][0]
      - gCov[0][3]*gCov[1][2]*gCov[2][0]*gCov[3][1];

  gDet.eval();

  gCon[0][0] = 
      (  gCov[1][1]*gCov[2][2]*gCov[3][3]
       + gCov[1][2]*gCov[2][3]*gCov[3][1]
       + gCov[1][3]*gCov[2][1]*gCov[3][2]
       - gCov[1][1]*gCov[2][3]*gCov[3][2]
       - gCov[1][2]*gCov[2][1]*gCov[3][3]
       - gCov[1][3]*gCov[2][2]*gCov[3][1])/gDet;

  gCon[0][1] = 
      (  gCov[0][1]*gCov[2][3]*gCov[3][2]
       + gCov[0][2]*gCov[2][1]*gCov[3][3]
       + gCov[0][3]*gCov[2][2]*gCov[3][1]
       - gCov[0][1]*gCov[2][2]*gCov[3][3]
       - gCov[0][2]*gCov[2][3]*gCov[3][1] 
       - gCov[0][3]*gCov[2][1]*gCov[3][2])/gDet;

  gCon[0][2] = 
      (  gCov[0][1]*gCov[1][2]*gCov[3][3]
       + gCov[0][2]*gCov[1][3]*gCov[3][1]
       + gCov[0][3]*gCov[1][1]*gCov[3][2]
       - gCov[0][1]*gCov[1][3]*gCov[3][2]
       - gCov[0][2]*gCov[1][1]*gCov[3][3]
       - gCov[0][3]*gCov[1][2]*gCov[3][1])/gDet;

  gCon[0][3] = 
      (  gCov[0][1]*gCov[1][3]*gCov[2][2]
       + gCov[0][2]*gCov[1][1]*gCov[2][3]
       + gCov[0][3]*gCov[1][2]*gCov[2][1]
       - gCov[0][1]*gCov[1][2]*gCov[2][3]
       - gCov[0][2]*gCov[1][3]*gCov[2][1]
       - gCov[0][3]*gCov[1][1]*gCov[2][2])/gDet;

  gCon[1][0] = gCon[0][1];
  
  gCon[1][1] = 
      (  gCov[0][0]*gCov[2][2]*gCov[3][3]
       + gCov[0][2]*gCov[2][3]*gCov[3][0]
       + gCov[0][3]*gCov[2][0]*gCov[3][2]
       - gCov[0][0]*gCov[2][3]*gCov[3][2]
       - gCov[0][2]*gCov[2][0]*gCov[3][3]
       - gCov[0][3]*gCov[2][2]*gCov[3][0])/gDet;

  gCon[1][2] = 
      (  gCov[0][0]*gCov[1][3]*gCov[3][2]
       + gCov[0][2]*gCov[1][0]*gCov[3][3]
       + gCov[0][3]*gCov[1][2]*gCov[3][0]
       - gCov[0][0]*gCov[1][2]*gCov[3][3]
       - gCov[0][2]*gCov[1][3]*gCov[3][0]
       - gCov[0][3]*gCov[1][0]*gCov[3][2])/gDet;

  gCon[1][3] = 
      (  gCov[0][0]*gCov[1][2]*gCov[2][3]
       + gCov[0][2]*gCov[1][3]*gCov[2][0]
       + gCov[0][3]*gCov[1][0]*gCov[2][2]
       - gCov[0][0]*gCov[1][3]*gCov[2][2]
       - gCov[0][2]*gCov[1][0]*gCov[2][3]
       - gCov[0][3]*gCov[1][2]*gCov[2][0])/gDet;

  gCon[2][0] = gCon[0][2];
  gCon[2][1] = gCon[1][2];

  gCon[2][2] =
      (  gCov[0][0]*gCov[1][1]*gCov[3][3]
       + gCov[0][1]*gCov[1][3]*gCov[3][0]
       + gCov[0][3]*gCov[1][0]*gCov[3][1]
       - gCov[0][0]*gCov[1][3]*gCov[3][1]
       - gCov[0][1]*gCov[1][0]*gCov[3][3]
       - gCov[0][3]*gCov[1][1]*gCov[3][0])/gDet;

  gCon[2][3] =
      (  gCov[0][0]*gCov[1][3]*gCov[2][1]
       + gCov[0][1]*gCov[1][0]*gCov[2][3]
       + gCov[0][3]*gCov[1][1]*gCov[2][0]
       - gCov[0][0]*gCov[1][1]*gCov[2][3]
       - gCov[0][1]*gCov[1][3]*gCov[2][0]
       - gCov[0][3]*gCov[1][0]*gCov[2][1])/gDet;

  gCon[3][0] = gCon[0][3];
  gCon[3][1] = gCon[1][3];
  gCon[3][2] = gCon[2][3];

  gCon[3][3] =
      (  gCov[0][0]*gCov[1][1]*gCov[2][2]  
       + gCov[0][1]*gCov[1][2]*gCov[2][0]  
       + gCov[0][2]*gCov[1][0]*gCov[2][1]  
       - gCov[0][0]*gCov[1][2]*gCov[2][1]  
       - gCov[0][1]*gCov[1][0]*gCov[2][2]  
       - gCov[0][2]*gCov[1][1]*gCov[2][0])/gDet;

  for(int d=0;d<4;d++)
    for(int dd=0;dd<4;dd++)
      gCon[d][dd].eval();

}

void geometry::setgCovInXCoords(const array XCoords[3], 
                                array gCov[NDIM][NDIM]
                               )
{
  switch (metric)
  {
    case metrics::MINKOWSKI:
      
      gCov[0][0] = -1.;
      gCov[1][1] = 1.;
      gCov[2][2] = 1.;
      gCov[3][3] = 1.;

      break;

    case metrics::MODIFIED_KERR_SCHILD:
      /* x^mu = {t, r, theta, phi}, X^mu = {t, X1, X2, phi} */

      array xCoords[3];

      XCoordsToxCoords(XCoords, xCoords);

      /* Easier to read with (r, theta) than (x[1], x[2]) */
      array r     = xCoords[directions::X1];
      array theta = xCoords[directions::X2];

      /* r = exp(X1) => dr/dX = exp(X1) = r */
      array dr_dX1 = r;
      array dr_dX2 = zero;
      /* theta = pi*X2 + 0.5*(1 - H_SLOPE)*sin(2*pi*X2) 
	  -         => dtheta/dX2 = pi + pi*(1 - H_SLOPE)*cos(2*pi*X2) */
      array dtheta_dX2 =  M_PI 
	+ M_PI*(1 - hSlope)
	*af::cos(2*M_PI*XCoords[directions::X2]);
      array dtheta_dX1 = zero;
      
      if(params::DoCylindrify)
	{
	  array XCoordsPlusEps[3];
	  array xCoordsPlusEps[3];
	  array XCoordsMinusEps[3];
	  array xCoordsMinusEps[3];
	  
	  for(int d=0;d<3;d++)
	    {
	      XCoordsPlusEps[d] = XCoords[d];
	      xCoordsPlusEps[d] = xCoords[d];
	      XCoordsMinusEps[d] = XCoords[d];
	      xCoordsMinusEps[d] = xCoords[d];
	    }
	  XCoordsPlusEps[directions::X1]+=GAMMA_EPS;
	  XCoordsMinusEps[directions::X1]-=GAMMA_EPS;
	  XCoordsToxCoords(XCoordsPlusEps,xCoordsPlusEps);      
	  XCoordsToxCoords(XCoordsMinusEps,xCoordsMinusEps);      
	  dr_dX1 = (xCoordsPlusEps[directions::X1]-xCoordsMinusEps[directions::X1])/(2.*GAMMA_EPS);
	  dtheta_dX1 = (xCoordsPlusEps[directions::X2]-xCoordsMinusEps[directions::X2])/(2.*GAMMA_EPS);

	  for(int d=0;d<3;d++)
	    {
	      XCoordsPlusEps[d] = XCoords[d];
	      xCoordsPlusEps[d] = xCoords[d];
	      XCoordsMinusEps[d] = XCoords[d];
              xCoordsMinusEps[d] = xCoords[d];
	    }
	  XCoordsPlusEps[directions::X2]+=GAMMA_EPS;
	  XCoordsMinusEps[directions::X2]-=GAMMA_EPS;
	  XCoordsToxCoords(XCoordsPlusEps,xCoordsPlusEps);      
	  XCoordsToxCoords(XCoordsMinusEps,xCoordsMinusEps);
	  dr_dX2 = (xCoordsPlusEps[directions::X1]-xCoordsMinusEps[directions::X1])/(2.*GAMMA_EPS);
	  dtheta_dX2 = (xCoordsPlusEps[directions::X2]-xCoordsMinusEps[directions::X2])/(2.*GAMMA_EPS);
	}

      dr_dX1.eval();
      dtheta_dX1.eval();
      dr_dX2.eval();
      dtheta_dX2.eval();

      array sigma =  r*r + af::pow(blackHoleSpin * af::cos(theta), 2.);
      sigma.eval();

      /* -(1 - 2*r/sigma) dt^2 */
      gCov[0][0] = -(1. - 2.*r/sigma);     

      /* (4*r/sigma * dr/dX1) dt dX1 */
      gCov[0][1] = (2.*r/sigma) * dr_dX1; 

      /* (4*r/sigma * dr/dX2) dt dX2 */
      gCov[0][2] = (2.*r/sigma) * dr_dX2;

      /* -(4*a*r*sin(theta)^2/sigma) dt dphi */
      gCov[0][3] = -(2.*blackHoleSpin*r*af::pow(af::sin(theta), 2.)/sigma);

      /* (4*r/sigma * dr/dX1) dX1 dt */
      gCov[1][0] = gCov[0][1];

      /* ( (1 + 2*r/sigma)*dr/dX1*dr/dX1 + sigma*dtheta_dX1*dtheta_dX1) dX1 dX1 */
      gCov[1][1] = (1. + 2*r/sigma) * dr_dX1 * dr_dX1
	+ (sigma) * dtheta_dX1 * dtheta_dX1;
  
      /* 2*((1 + 2*r/sigma)*dr/dX1*dr/dX2 + sigma * dtheta_dX1 * dtheta_dX2) dX1 dX2 */
      gCov[1][2] = (1. + 2*r/sigma) * dr_dX1 * dr_dX2
	+ (sigma) * dtheta_dX1 * dtheta_dX2;

      /* -(2*a*(1 + 2.*r/sigma)*sin(theta)^2*dr/dX1) dX1 dphi */
      gCov[1][3] =
        -blackHoleSpin*(1. + 2.*r/sigma)*af::pow(af::sin(theta), 2.)*dr_dX1;

      /* (0) dX2 dt */
      gCov[2][0] = gCov[0][2];

      /* (0) dX2 dX1 */
      gCov[2][1] = gCov[1][2];

      /* (sigma*dtheta/dX2*dtheta/dX2 + (1 + 2*r/sigma)*dr/dX2*dr/dX2) dX2 dX2 */
      gCov[2][2] = sigma*dtheta_dX2*dtheta_dX2
	+ (1. + 2*r/sigma) * dr_dX2 * dr_dX2;

      /* -(2*a*(1 + 2.*r/sigma)*sin(theta)^2*dr/dX2) dX2 dphi */
      gCov[2][3] = 
	-blackHoleSpin*(1. + 2.*r/sigma)*af::pow(af::sin(theta), 2.)*dr_dX2;

      /* -(4*a*r*sin(theta)^2/sigma) dphi dt */
      gCov[3][0] = gCov[0][3];

      /* -(2*a*(1 + 2.*r/sigma)*sin(theta)^2*dr/dX1) dX1 dphi */
      gCov[3][1] = gCov[1][3];

      /* (0) dphi dX2 */
      gCov[3][2] = gCov[2][3];

      /* (sin(theta)^2*(sigma + a^2*(1. + 2*r/sigma)*sin(theta)^2) dphi dphi */
      gCov[3][3] = (  af::pow(af::sin(theta), 2.)
                    * (sigma +   af::pow(blackHoleSpin*af::sin(theta), 2.)
                               * (1. + 2*r/sigma)
                      )
                   );

      for(int d=0;d<4;d++)
	for(int dd=0;dd<4;dd++)
	  gCov[d][dd].eval();
      
      break;
  }


}


void geometry::XCoordsToxCoords(const array XCoords[3], 
                                array xCoords[3]
                               ) const
{
  switch (metric)
  {
    case metrics::MINKOWSKI:
      
      xCoords[directions::X1] = XCoords[directions::X1];
      xCoords[directions::X2] = XCoords[directions::X2];
      xCoords[directions::X3] = XCoords[directions::X3];

      break;

    case metrics::MODIFIED_KERR_SCHILD:
      
      xCoords[directions::X1] = GammieRadius(XCoords[directions::X1]);
      xCoords[directions::X2] = ThetaNoCyl(XCoords[directions::X1],XCoords[directions::X2]); 
      xCoords[directions::X3] = XCoords[directions::X3];

      xCoords[directions::X1].eval();
      xCoords[directions::X2].eval();
      xCoords[directions::X3].eval();
      
      // 'Cylindrify' the coordinates, following a slightly
      // adapted version of Sasha Tchekhovskoy's HARM code.
      // This modifies the map Theta(X1,X2) so that X2=cst surfaces
      // close to the pole axis and the black hole are wide
      // cylinders instead of cones.
      if(params::DoCylindrify)
	{

	  array Xgrid[3];
	  array xGam[3];
	  for(int d=0;d<3;d++)
	    {
	      Xgrid[d]=XCoords[d];
	      xGam[d]=xCoords[d];
	    }
	  // Transform lower hemisphere and ghost zones
	  // to upper hemisphere, 0<theta<pi/2
	  array IsMirrored = Xgrid[directions::X2]>0.5;
	  Xgrid[directions::X2] = Xgrid[directions::X2]*(1-IsMirrored)
	    + IsMirrored*(1.-Xgrid[directions::X2]);
	  xGam[directions::X2] = xGam[directions::X2]*(1-IsMirrored)
	    + IsMirrored* (M_PI-xGam[directions::X2]);
	  array IsGhost = Xgrid[directions::X2]<0.;
	  Xgrid[directions::X2] = Xgrid[directions::X2]*(1-IsGhost)
            + IsGhost*(-(1.0)*Xgrid[directions::X2]);
	  xGam[directions::X2] = xGam[directions::X2]*(1-IsGhost)
            + IsGhost* (-(1.0)*xGam[directions::X2]);

	  // We cylindrify the region with 
	  // X1 < params::X1cyl
	  // X2 < params::X2cyl
	  array X0[3];
	  for(int d=0;d<3;d++)
	    X0[d] = Xgrid[d]*1.;
	  X0[directions::X1] = params::X1cyl;
	  X0[directions::X2] = params::X2cyl;
	  X0[directions::X3]=0.;
	  array R_cyl = GammieRadius(X0[directions::X1]);
	  array Theta_cyl = ThetaNoCyl(X0[directions::X1],X0[directions::X2]);
	  array X1_tr = af::log(0.5*(R_cyl+exp(params::X1Start)));
	  
	  array Xtr[3];
	  array xcyl[3];
	  for(int d=0;d<3;d++)
	    {
	      Xtr[d]=Xgrid[d];
	      xcyl[d]=xGam[d];
	    }
	  Xtr[directions::X1]=X1_tr;
	  
	  // Sasha's new theta
	  array f1 = func1(X0,Xgrid);
	  array f2 = func2(X0,Xgrid);
	  array dftr = func2(X0,Xtr)-func1(X0,Xtr);
	  array sinth = maxs(xGam[directions::X1]*f1, xGam[directions::X1]*f2, GammieRadius(Xtr[directions::X1])*af::abs(dftr)+1.e-20 ) / xGam[directions::X1]; 
	  array th = af::asin(sinth);
	  th.eval();

	  // Move ghost zones and lower hemisphere points
	  // back to their correct theta
	  th = th*(1-IsGhost*2);
	  xcyl[directions::X2] = IsMirrored*(M_PI-th)
	    +(1-IsMirrored)*th;
	  for(int d=0;d<3;d++)
	    {
	      xCoords[d]=xcyl[d];
	      xCoords[d].eval();
	    }
	}
      
      break;
  }
}

void geometry::computeGammaDownDownDown(const int eta,
                                        const int mu,
                                        const int nu,
                                        array& out,
					const array gCovPlus[NDIM-1][NDIM][NDIM],
					const array gCovMinus[NDIM-1][NDIM][NDIM]
                                       )
{
  out = zero;

  /* First, take care of d(g_eta_mu)/dX^nu. To compute this numerically we
   * first do XEpsilon^alpha[nu] = X^alpha[nu] + EPS and then compute the
   * metric corresponding to XEpsilon^alpha coordinates and add this to *ans
   * Now we do XEpsilon^alpha[nu] = X^alpha[nu] - EPS and then compute the
   * metric corresponding to XEpsilon^alpha coordinates and then subtract it
   * from *ans. And this entire thing is divided by 2.*EPS, corresponding to a
   * centered difference around X^alpha. By doing the computations seperately
   * for +EPS and -EPS, we cut storage requirements by half.*/

  if(nu>0)
    {
      /* Handle +EPS first */
      out += 0.5*gCovPlus[nu-1][eta][mu]/(2.*GAMMA_EPS); 
      /* Now do -EPS */
      out -= 0.5*gCovMinus[nu-1][eta][mu]/(2.*GAMMA_EPS);
      /* End of d(g_eta_mu)/dX^nu */
    }

  /* Now, d(g_eta_nu)/dX^mu */
  if(mu>0)
    {
      /* +EPS first */
      out += 0.5*gCovPlus[mu-1][eta][nu]/(2.*GAMMA_EPS);
      /* Now do -EPS */
      out -= 0.5*gCovMinus[mu-1][eta][nu]/(2.*GAMMA_EPS);
      /* End of d(g_eta_nu)/dX^mu */
    }

  /* Finally, d(g_mu_nu)/dX^eta */
  if(eta>0)
    {
      /* +EPS first */
      out -= 0.5*gCovPlus[eta-1][mu][nu]/(2.*GAMMA_EPS);
      /* Now do -EPS */
      out += 0.5*gCovMinus[eta-1][mu][nu]/(2.*GAMMA_EPS);
      /* End of d(g_mu_nu)/dX^eta */
    }

  out.eval();
}

geometry::~geometry()
{
}

void geometry::setgCovGrid()
{
  gCovGrid = new grid(N1, N2, N3, dim, 16, numGhost,
                      false, false, false
                     );

  for (int mu=0; mu<NDIM; mu++)
  {
    for (int nu=0; nu<NDIM; nu++)
    {
      gCovGrid->vars[nu + NDIM*mu] = gCov[mu][nu];
    }
  }
}

void geometry::setgConGrid()
{
  gConGrid = new grid(N1, N2, N3, dim, 16, numGhost,
                      false, false, false
                     );

  for (int mu=0; mu<NDIM; mu++)
  {
    for (int nu=0; nu<NDIM; nu++)
    {
      gConGrid->vars[nu + NDIM*mu] = gCon[mu][nu];
    }
  }
}

void geometry::setgGrid()
{
  gGrid = new grid(N1, N2, N3, dim, 1, numGhost,
                   false, false, false
                  );

  gGrid->vars[0] = g;
}

void geometry::setalphaGrid()
{
  alphaGrid = new grid(N1, N2, N3, dim, 1, numGhost,
                       false, false, false
                      );

  alphaGrid->vars[0] = alpha;
}

void geometry::setgammaUpDownDownGrid()
{
  gammaUpDownDownGrid = new grid(N1, N2, N3, dim, 64, numGhost,
                                 false, false, false
                                );

  for (int mu=0; mu<NDIM; mu++)
  {
    for (int nu=0; nu<NDIM; nu++)
    {
      for (int lamda = 0; lamda<NDIM; lamda++)
      {
        gammaUpDownDownGrid->vars[lamda + NDIM*(nu + NDIM*(mu))]
          = gammaUpDownDown[mu][nu][lamda];
      }
    }
  }
}

void geometry::setxCoordsGrid()
{
  xCoordsGrid = new grid(N1, N2, N3, dim, 3, numGhost,
                         false, false, false
                        );
  array xCoords[3];
  XCoordsToxCoords(XCoords, xCoords);

  xCoordsGrid->vars[directions::X1] = xCoords[directions::X1];
  xCoordsGrid->vars[directions::X2] = xCoords[directions::X2];
  xCoordsGrid->vars[directions::X3] = xCoords[directions::X3];
}
