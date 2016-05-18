// Original HARM coordinate map R(X1), Theta(X2)
// from Charles Gammie

array GammieRadius(const array X1)
{
  array res = af::exp(X1);
  res.eval();
  return res;
}

array GammieTheta(const array X2)
{
  array res = M_PI*X2+ 0.5*(1 - params::hSlope)*af::sin(2.*M_PI*X2);
  res.eval();
  return res;
}


// Modify Theta without cylindrification, but
// setting a constant spacing on the horizon
array Ftr( const array x );
array ThetaNoCyl(const array X1, const array X2)
{
  if(params::DerefineThetaHorizon)
    {
      // This leads to Theta(X2) as in HARM at large distances
      // and Theta ~ X2 for small radii:
      // X2=cst lines are straight lines in the r-z plane
      // The lines are equally spaced on the excision surface
      // but use increase resolution in the equatorial plane
      // at large radii, as in HARM
      array cTh0 = af::cos(M_PI*X2);
      array cTh1 = af::cos(GammieTheta(X2));
      array cTh = cTh1 + exp(params::X1Start)/af::exp(X1)*(cTh0-cTh1);
      array res = af::acos(cTh);
      res.eval();
      return res;
    }
  else
    {
      // Use standard HARM choice
      array res = GammieTheta(X2);
      res.eval();
      return res;
    }
}


/**************************************
// The following code is largely copied
// from Sasha Tchekhovskoy's version of HARM,
// adapted to arrayFire (and without the generic mapping)
***************************************/

//smooth step function:
// Ftr = 0 if x < 0, Ftr = 1 if x > 1 and smoothly interps. in btw.
array Ftr( const array x )
{
  array condition1 = (x<=0);
  array condition2 = (x>=1);
  array res = condition2;
  res = res + (condition1-1)*(condition2-1)*
    (af::cos(x*5*M_PI) + 70*af::sin(((2*x-1)*M_PI)/2.) + 5*af::sin(((2*x-1)*3*M_PI)/2.) + 64.)/128.;
  res.eval();
  return res;
}

array Ftrgenlin( const array x, const array xa, const array xb, const array ya, const array yb )
{
  array res = (x*ya)/xa + (-((x*ya)/xa) + ((x - xb)*(1 - yb))/(1 - xb) + yb)*Ftr((x - xa)/(xb - xa));
  res.eval();
  return res;
}

//goes from ya to yb as x goes from xa to xb
array Ftrgen(const array x, const array xa, const array xb, const array ya, const array yb )
{
  array res = ya + (yb-ya)*Ftr( (x-xa)/(xb-xa) );
  res.eval();
  return res;
}

array Fangle(const array x )
{
  array condition1 = (x<=-1);
  array condition2 = (x>=1);
  array res = condition2*x + (condition1-1)*(condition2-1)*
    (x + 1. + (-140*af::sin(((x+1)*M_PI)/2.) + (10*af::sin(((x+1)*3*M_PI)/2.))/3. + (2*af::sin(((x+1)*5*M_PI)/2.))/5.)/(64.*M_PI))/2.;
  res.eval();
  return res;
}

array limlin( const array x, const array x0, const array dx, const array y0 )
{
  array res = y0 - dx * Fangle((x-x0)*(-1.)/dx);
  res.eval();
  return res;
}

array minlin(const array x,const array x0,const array dx,const array y0 )
{
  array res = y0 + dx * Fangle((x-x0)/dx);
  res.eval();
  return res;
}

array mins(const array f1,const array f2,const array df )
{
  return limlin(f1, f2, df, f2);
}

array maxs(const array f1,const array f2,const array df )
{
  return (mins(-f1, -f2, df)*(-1.) );
}

//=mins if dir < 0
//=maxs if dir >= 0
array minmaxs(const array f1,const array f2,const array df,const array dir )
{
  array condition = (dir>=0);
  array res = condition*maxs(f1, f2, df)
    +(1-condition)*mins(f1, f2, df);
  res.eval();
  return res;
}


array sinth0( array X0[3], array X[3])
{
  array Xc0[3];
  for(int j=0;j<3;j++)
    Xc0[j]=X[j];
  Xc0[directions::X2] = X0[directions::X2];
  array res = GammieRadius(X0[directions::X1]) * af::sin(ThetaNoCyl(X0[directions::X1],X0[directions::X2])) / GammieRadius(Xc0[directions::X1]);
  res.eval();
  return res;
}

array sinth1in( array X0[3], array X[3])
{
  array X0c[3];
  for(int j=0;j<3;j++)
    X0c[j]=X0[j];
  X0c[directions::X2] = X[directions::X2];
  array res = GammieRadius(X0[directions::X1]) * af::sin(ThetaNoCyl(X0c[directions::X1],X0c[directions::X2])) / GammieRadius(X[directions::X1]); 
  res.eval();
  return res;
}


array th2in( array X0[3], array X[3])
{
  array Xc0[3];
  array Xcmid[3];
  for(int j=0;j<3;j++)
    Xc0[j]=X[j];
  Xc0[directions::X2] = X0[directions::X2];
  for(int j=0;j<3;j++)
    Xcmid[j]=X[j];
  Xcmid[directions::X2] = 0.5;
  
  array th0 = af::asin( sinth0(X0, X) );
  th0.eval();
  array res = (ThetaNoCyl(X[directions::X1],X[directions::X2]) - ThetaNoCyl(Xc0[directions::X1],Xc0[directions::X2]))/
    (ThetaNoCyl(Xcmid[directions::X1],Xcmid[directions::X2]) - ThetaNoCyl(Xc0[directions::X1],Xc0[directions::X2])) *
    (ThetaNoCyl(Xcmid[directions::X1],Xcmid[directions::X2])-th0) + th0;
  res.eval();
  return res;
}

array func1( array X0[3], array X[3])
{
  array res = af::sin(ThetaNoCyl(X[directions::X1],X[directions::X2]));
  res.eval();
  return res;
}

array func2( array X0[3], array X[3])
{
  array Xca[3];
  for(int j=0;j<3;j++)
    Xca[j]=X[j];
  Xca[directions::X2] = 0.;
  
  array sth1in = sinth1in( X0, X);
  array sth2in = af::sin( th2in(X0, X));
  
  array sth1inaxis = sinth1in( X0, Xca);
  array sth2inaxis = af::sin( th2in(X0, Xca) );
  
  array func2 = minmaxs( sth1in, sth2in, af::abs(sth2inaxis-sth1inaxis)+1.e-20, X[directions::X1] - X0[directions::X1] );
  func2.eval();
  return func2;
}
