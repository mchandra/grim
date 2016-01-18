#include "geometry.hpp"

//void geometry::getHostPtrTo(std::string str)
//{
//  if (str=="xCoords")
//  {
//    xCoordsHostPtr = 
//
//
//  }
//
//}

geometry::geometry(grid &XCoords)
{
  this->XCoords = &XCoords;

  /* Allocate space */
  g          = zero;
  array gDet = zero;
  for (int mu=0; mu<NDIM; mu++)
  {
    for (int nu=0; nu<NDIM; nu++)
    {
      gCov[mu][nu] = zero;
      gCon[mu][nu] = zero;

      for (int lamda = 0; lamda<NDIM; lamda++)
      {
      	gammaUpDownDown[lamda][mu][nu] = zero;
      }
    }
  }
  
  setgCovInXCoords(XCoords, gCov);

  setgDetAndgConFromgCov(gCov, gDet, gCon);
  g = af::sqrt(-gDet);

  alpha = 1./af::sqrt(-gCon[0][0]);
}

void geometry::setXCoords(grid &indices, int location, grid &XCoords)
{
  array indicesX1 = indices.vars[directions::X1];
  array indicesX2 = indices.vars[directions::X2];
  array indicesX3 = indices.vars[directions::X3];

  double dX1 = indices.dX1;
  double dX2 = indices.dX2;
  double dX3 = indices.dX3;

  switch (location)
  {
    case CENTER:
      XCoords.vars[directions::X1] = params::X1Start + (indicesX1 + 0.5)*dX1;
      XCoords.vars[directions::X2] = params::X2Start + (indicesX2 + 0.5)*dX2;
      XCoords.vars[directions::X3] = params::X3Start + (indicesX3 + 0.5)*dX3;

      break;

    case LEFT:
      XCoords.vars[directions::X1] = params::X1Start + (indicesX1      )*dX1;
      XCoords.vars[directions::X2] = params::X2Start + (indicesX2 + 0.5)*dX2;
      XCoords.vars[directions::X3] = params::X3Start + (indicesX3 + 0.5)*dX3;

      break;

    case RIGHT:
      XCoords.vars[directions::X1] = params::X1Start + (indicesX1 +   1)*dX1;
      XCoords.vars[directions::X2] = params::X2Start + (indicesX2 + 0.5)*dX2;
      XCoords.vars[directions::X3] = params::X3Start + (indicesX3 + 0.5)*dX3;

      break;

    case BOTTOM:
      XCoords.vars[directions::X1] = params::X1Start + (indicesX1 + 0.5)*dX1;
      XCoords.vars[directions::X2] = params::X2Start + (indicesX2      )*dX2;
      XCoords.vars[directions::X3] = params::X3Start + (indicesX3 + 0.5)*dX3;

      break;

    case TOP:
      XCoords.vars[directions::X1] = params::X1Start + (indicesX1 + 0.5)*dX1;
      XCoords.vars[directions::X2] = params::X2Start + (indicesX2 +   1)*dX2;
      XCoords.vars[directions::X3] = params::X3Start + (indicesX3 + 0.5)*dX3;

      break;

    case FRONT:
      XCoords.vars[directions::X1] = params::X1Start + (indicesX1 + 0.5)*dX1;
      XCoords.vars[directions::X2] = params::X2Start + (indicesX2 + 0.5)*dX2;
      XCoords.vars[directions::X3] = params::X3Start + (indicesX3 +   1)*dX3;

      break;

    case BACK:
      XCoords.vars[directions::X1] = params::X1Start + (indicesX1 + 0.5)*dX1;
      XCoords.vars[directions::X2] = params::X2Start + (indicesX2 + 0.5)*dX2;
      XCoords.vars[directions::X3] = params::X3Start + (indicesX3      )*dX3;

      break;
  }

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

void geometry::setgCovInXCoords(array gCov[NDIM][NDIM])
{
  switch (params::metric)
  {
    case metrics::MINKOWSKI:
      
      gCov[0][0] = -1.;
      gCov[1][1] = 1.;
      gCov[2][2] = 1.;
      gCov[3][3] = 1.;

      break;

    case metrics::MODIFIED_KERR_SCHILD:
      /* x^mu = {t, r, theta, phi}, X^mu = {t, X1, X2, phi} */

      array xCoords[NDIM];

      XCoordsToxCoords(XCoords, xCoords);

      /* Easier to read with (r, theta) than (x[1], x[2]) */
      array r     = xCoords[directions::X1];
      array theta = xCoords[directions::X2];

      /* r = exp(X1) => dr/dX = exp(X1) = r */
      array dr_dX1 = r;

      /* theta = pi*X2 + 0.5*(1 - H_SLOPE)*sin(2*pi*X2) 
         => dtheta/dX2 = pi + pi*(1 - H_SLOPE)*cos(2*pi*X2) */
      array dtheta_dX2 = M_PI + M_PI*(1 - params::hSlope)*af::cos(2*M_PI*XCoords[directions::X2]);

      array sigma =  r*r + af::pow2(params::blackHoleSpin * af::cos(theta) );

      /* -(1 - 2*r/sigma) dt^2 */
      gCov[0][0] = -(1. - 2.*r/sigma);     

      /* (4*r/sigma * dr/dX1) dt dX1 */
      gCov[0][1] = (2.*r/sigma) * dr_dX1; 

      /* (0) dt dX2 */
      gCov[0][2] = 0.; 

      /* -(4*a*r*sin(theta)^2/sigma) dt dphi */
      gCov[0][3] = -(2.*params::blackHoleSpin*r*af::pow2(af::sin(theta))/sigma);

      /* (4*r/sigma * dr/dX1) dX1 dt */
      gCov[1][0] = gCov[0][1];

      /* ( (1 + 2*r/sigma)*dr/dX1*dr/dX1) dX1 dX1 */
      gCov[1][1] = (1. + 2*r/sigma) * dr_dX1 * dr_dX1;
  
      /* (0) dX1 dX2 */
      gCov[1][2] = 0.;

      /* -(2*a*(1 + 2.*r/sigma)*sin(theta)^2*dr/dX1) dX1 dphi */
      gCov[1][3] =
        -params::blackHoleSpin*(1. + 2.*r/sigma)*af::pow2(af::sin(theta))*dr_dX1;

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
      gCov[3][3] = (  af::pow2(af::sin(theta))
                    * (sigma +   af::pow2(params::blackHoleSpin*af::sin(theta))
                               * (1. + 2*r/sigma)
                      )
                   );
      
      break;
  }


}

void geometry::XCoordsToxCoords(const array XCoords[NDIM], 
                                array xCoords[NDIM]
                               )
{
  switch (params::metric)
  {
    case metrics::MINKOWSKI:
      
      for (int mu=0; mu<NDIM; mu++)
      {
        xCoords[mu] = XCoords[mu];
      }
      break;

    case metrics::MODIFIED_KERR_SCHILD:
      
      xCoords[0] = XCoords[0];
      xCoords[1] = af::exp(XCoords[1]);
      xCoords[2] =   M_PI*XCoords[2] 
                   + 0.5*(1 - params::hSlope)
                   * af::sin(2.*M_PI*XCoords[2]);
  
      xCoords[3] = XCoords[3];
      break;
  }
}


geometry::~geometry()
{
}


