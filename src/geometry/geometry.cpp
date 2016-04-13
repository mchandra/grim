#include "geometry.hpp"

geometry::geometry(const int metric,
                   const double blackHoleSpin,
                   const double hSlope,
                   const coordinatesGrid &XCoordsGrid
                  )
{
  location  = XCoordsGrid.location;
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
  for (int mu=0; mu<NDIM; mu++)
  {
    for (int nu=0; nu<NDIM; nu++)
    {
      gCov[mu][nu].eval();
      gCon[mu][nu].eval();
    }
  }

  af::sync();
}

void geometry::computeConnectionCoeffs()
{
  array gammaDownDownDown[NDIM][NDIM][NDIM];

  for (int mu=0; mu<NDIM; mu++)
  {
    for (int nu=0; nu<NDIM; nu++)
    {
      for (int lamda = 0; lamda<NDIM; lamda++)
      {
      	gammaDownDownDown[mu][nu][lamda] = zero;
        computeGammaDownDownDown(mu,nu,lamda,
                                 gammaDownDownDown[mu][nu][lamda]
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

      /* theta = pi*X2 + 0.5*(1 - H_SLOPE)*sin(2*pi*X2) 
         => dtheta/dX2 = pi + pi*(1 - H_SLOPE)*cos(2*pi*X2) */
      array dtheta_dX2 =  M_PI 
                        + M_PI*(1 - hSlope)
                              *af::cos(2*M_PI*XCoords[directions::X2]);

      array sigma =  r*r + af::pow(blackHoleSpin * af::cos(theta), 2.);

      /* -(1 - 2*r/sigma) dt^2 */
      gCov[0][0] = -(1. - 2.*r/sigma);     

      /* (4*r/sigma * dr/dX1) dt dX1 */
      gCov[0][1] = (2.*r/sigma) * dr_dX1; 

      /* (0) dt dX2 */
      gCov[0][2] = 0.; 

      /* -(4*a*r*sin(theta)^2/sigma) dt dphi */
      gCov[0][3] = -(2.*blackHoleSpin*r*af::pow(af::sin(theta), 2.)/sigma);

      /* (4*r/sigma * dr/dX1) dX1 dt */
      gCov[1][0] = gCov[0][1];

      /* ( (1 + 2*r/sigma)*dr/dX1*dr/dX1) dX1 dX1 */
      gCov[1][1] = (1. + 2*r/sigma) * dr_dX1 * dr_dX1;
  
      /* (0) dX1 dX2 */
      gCov[1][2] = 0.;

      /* -(2*a*(1 + 2.*r/sigma)*sin(theta)^2*dr/dX1) dX1 dphi */
      gCov[1][3] =
        -blackHoleSpin*(1. + 2.*r/sigma)*af::pow(af::sin(theta), 2.)*dr_dX1;

      /* (0) dX2 dt */
      gCov[2][0] = gCov[0][2];

      /* (0) dX2 dX1 */
      gCov[2][1] = gCov[1][2];

      /* (sigma*dtheta/dX2*dtheta/dX2) dX2 dX2 */
      gCov[2][2] = sigma*dtheta_dX2*dtheta_dX2;

      /* (0) dX2 dphi */
      gCov[2][3] = 0.;

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
      
      xCoords[directions::X1] = af::exp(XCoords[directions::X1]);
      xCoords[directions::X2] =  M_PI*XCoords[directions::X2] 
                               + 0.5*(1 - hSlope)
                               * af::sin(2.*M_PI*XCoords[directions::X2]);
  
      xCoords[directions::X3] = 2*M_PI*XCoords[directions::X3];

      break;
  }
}

void geometry::computeGammaDownDownDown(const int eta,
                                        const int mu,
                                        const int nu,
                                        array& out
                                       )
{
  const double GAMMA_EPS=1.e-5;
  array XCoords_4D[NDIM];
  array XEpsilon[NDIM];
  array XEpsilonSpatial[3];
  array gCovEpsilon[NDIM][NDIM];

  XCoords_4D[0] = zero;
  XCoords_4D[1] = this->XCoords[directions::X1];
  XCoords_4D[2] = this->XCoords[directions::X2];
  XCoords_4D[3] = this->XCoords[directions::X3];

  for (int alpha=0; alpha<NDIM; alpha++)
  {
    XEpsilon[alpha]=zero;
      
    for (int beta=0; beta<NDIM; beta++)
    {
      gCovEpsilon[alpha][beta]=zero;
    }
  }
  out = zero;

  /* First, take care of d(g_eta_mu)/dX^nu. To compute this numerically we
   * first do XEpsilon^alpha[nu] = X^alpha[nu] + EPS and then compute the
   * metric corresponding to XEpsilon^alpha coordinates and add this to *ans
   * Now we do XEpsilon^alpha[nu] = X^alpha[nu] - EPS and then compute the
   * metric corresponding to XEpsilon^alpha coordinates and then subtract it
   * from *ans. And this entire thing is divided by 2.*EPS, corresponding to a
   * centered difference around X^alpha. By doing the computations seperately
   * for +EPS and -EPS, we cut storage requirements by half.*/

  /* Handle +EPS first */
  for (int alpha=0; alpha<NDIM; alpha++)
  {
    XEpsilon[alpha]=XCoords_4D[alpha];
  }
  XEpsilon[nu] += GAMMA_EPS;
  /* setgCovInXCoords takes in spatial XCoords. Hence the following piece of code */
  for (int dir=directions::X1; dir<=directions::X3; dir++)
  {
    XEpsilonSpatial[dir] = XEpsilon[dir+1];
  }
  setgCovInXCoords(XEpsilonSpatial,gCovEpsilon);
  out += 0.5*gCovEpsilon[eta][mu]/(2.*GAMMA_EPS); 

  /* Now do -EPS */
  for (int alpha=0; alpha<NDIM; alpha++)
  {
    XEpsilon[alpha] = XCoords_4D[alpha];
  }
  XEpsilon[nu] -= GAMMA_EPS;
  /* setgCovInXCoords takes in spatial XCoords. Hence the following piece of code */
  for (int dir=directions::X1; dir<=directions::X3; dir++)
  {
    XEpsilonSpatial[dir] = XEpsilon[dir+1];
  }
  setgCovInXCoords(XEpsilonSpatial,gCovEpsilon);
  out -= 0.5*gCovEpsilon[eta][mu]/(2.*GAMMA_EPS);
  /* End of d(g_eta_mu)/dX^nu */

  /* Now, d(g_eta_nu)/dX^mu */
  /* +EPS first */
  for (int alpha=0; alpha<NDIM; alpha++)
  {
    XEpsilon[alpha] = XCoords_4D[alpha];
  }
  XEpsilon[mu] += GAMMA_EPS;
  /* setgCovInXCoords takes in spatial XCoords. Hence the following piece of code */
  for (int dir=directions::X1; dir<=directions::X3; dir++)
  {
    XEpsilonSpatial[dir] = XEpsilon[dir+1];
  }
  setgCovInXCoords(XEpsilonSpatial,gCovEpsilon);
  out += 0.5*(gCovEpsilon[eta][nu])/(2.*GAMMA_EPS);

  /* Now do -EPS */
  for (int alpha=0; alpha<NDIM; alpha++)
  {
    XEpsilon[alpha] = XCoords_4D[alpha];
  }
  XEpsilon[mu] -= GAMMA_EPS;
  /* setgCovInXCoords takes in spatial XCoords. Hence the following piece of code */
  for (int dir=directions::X1; dir<=directions::X3; dir++)
  {
    XEpsilonSpatial[dir] = XEpsilon[dir+1];
  }
  setgCovInXCoords(XEpsilonSpatial,gCovEpsilon);
  out -= 0.5*gCovEpsilon[eta][nu]/(2.*GAMMA_EPS);
  /* End of d(g_eta_nu)/dX^mu */

  /* Finally, d(g_mu_nu)/dX^eta */
  /* +EPS first */
  for (int alpha=0; alpha<NDIM; alpha++)
  {
    XEpsilon[alpha] = XCoords_4D[alpha];
  }
  XEpsilon[eta] += GAMMA_EPS;
  /* setgCovInXCoords takes in spatial XCoords. Hence the following piece of code */
  for (int dir=directions::X1; dir<=directions::X3; dir++)
  {
    XEpsilonSpatial[dir] = XEpsilon[dir+1];
  }
  setgCovInXCoords(XEpsilonSpatial,gCovEpsilon);
  out -= 0.5*(gCovEpsilon[mu][nu])/(2.*GAMMA_EPS);

  /* Now do -EPS */
  for (int alpha=0; alpha<NDIM; alpha++)
  {
    XEpsilon[alpha] = XCoords_4D[alpha];
  }
  XEpsilon[eta] -= GAMMA_EPS;
  /* setgCovInXCoords takes in spatial XCoords. Hence the following piece of code */
  for (int dir=directions::X1; dir<=directions::X3; dir++)
  {
    XEpsilonSpatial[dir] = XEpsilon[dir+1];
  }
  setgCovInXCoords(XEpsilonSpatial,gCovEpsilon);
  out += 0.5*gCovEpsilon[mu][nu]/(2.*GAMMA_EPS);
  /* End of d(g_mu_nu)/dX^eta */

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
