#include "physics.h"

void fixedQuadIntegration(const struct fluidElement *elem,
                          const struct geometry *geom,
                          REAL scaleFactor,
                          REAL *moments)
{
  /* Using roots of the Legendre polynomaials for integration from [-1, 1].
     We get the quad points and the weights from 
     scipy.special.orthogonal.p_roots() */
  REAL quadPts[NUM_QUAD] = \
      {-9.98909991e-01,  -9.94261260e-01,  -9.85915992e-01,
       -9.73903368e-01,  -9.58267849e-01,  -9.39067544e-01,
       -9.16373862e-01,  -8.90271218e-01,  -8.60856711e-01,
       -8.28239764e-01,  -7.92541712e-01,  -7.53895354e-01,
       -7.12444458e-01,  -6.68343221e-01,  -6.21755705e-01,
       -5.72855216e-01,  -5.21823669e-01,  -4.68850904e-01,
       -4.14133983e-01,  -3.57876457e-01,  -3.00287606e-01,
       -2.41581666e-01,  -1.81977027e-01,  -1.21695421e-01,
       -6.09611002e-02,   3.98821425e-16,   6.09611002e-02,
        1.21695421e-01,   1.81977027e-01,   2.41581666e-01,
        3.00287606e-01,   3.57876457e-01,   4.14133983e-01,
        4.68850904e-01,   5.21823669e-01,   5.72855216e-01,
        6.21755705e-01,   6.68343221e-01,   7.12444458e-01,
        7.53895354e-01,   7.92541712e-01,   8.28239764e-01,
        8.60856711e-01,   8.90271218e-01,   9.16373862e-01,
        9.39067544e-01,   9.58267849e-01,   9.73903368e-01,
        9.85915992e-01,   9.94261260e-01,   9.98909991e-01};

  REAL weights[NUM_QUAD] = \
      {0.00279681,  0.00650034,  0.01018519,  0.01383263,  0.01742871,
       0.02095999,  0.0244133 ,  0.0277758 ,  0.03103497,  0.03417869,
       0.03719527,  0.04007348,  0.04280261,  0.04537251,  0.04777363,
       0.04999702,  0.05203442,  0.05387825,  0.05552165,  0.05695851,
       0.05818347,  0.05919199,  0.05998032,  0.06054551,  0.06088546,
       0.06099892,  0.06088546,  0.06054551,  0.05998032,  0.05919199,
       0.05818347,  0.05695851,  0.05552165,  0.05387825,  0.05203442,
       0.04999702,  0.04777363,  0.04537251,  0.04280261,  0.04007348,
       0.03719527,  0.03417869,  0.03103497,  0.0277758 ,  0.0244133 ,
       0.02095999,  0.01742871,  0.01383263,  0.01018519,  0.00650034,
       0.00279681};
  
  #if (GYROAVERAGING)
    /* Using roots of the shifte Legendre polynomaials for integration 
       from [0, 1] (i.e. [0, inf) in the perpendicular direction).
       We get the quad points and the weights from 
       scipy.special.orthogonal.ps_roots() */

    REAL quadPtsPerp[NUM_QUAD_PERP] = \
      {0.00222152,  0.01166804,  0.02851271,  0.052504  ,  0.08327869,
       0.12037037,  0.16321682,  0.21116853,  0.26349863,  0.31941385,
       0.37806656,  0.43856765,  0.5       ,  0.56143235,  0.62193344,
       0.68058615,  0.73650137,  0.78883147,  0.83678318,  0.87962963,
       0.91672131,  0.947496  ,  0.97148729,  0.98833196,  0.99777848};
  
    REAL weightsPerp[NUM_QUAD_PERP] = \
      {0.0056969 ,  0.01317749,  0.02046958,  0.02745235,  0.03401917,
       0.04007035,  0.04551413,  0.05026797,  0.05425981,  0.05742913,
       0.05972788,  0.06112122,  0.06158803,  0.06112122,  0.05972788,
       0.05742913,  0.05425981,  0.05026797,  0.04551413,  0.04007035,
       0.03401917,  0.02745235,  0.02046958,  0.01317749,  0.0056969};

  #endif

  REAL momentsInOrthTetrad[NUM_ALL_COMPONENTS];

  for (int mu=0; mu<NDIM; mu++)
  {
    momentsInOrthTetrad[N_UP(mu)] = 0.;

    for (int nu=0; nu<NDIM; nu++)
    {
      momentsInOrthTetrad[T_UP_UP(mu, nu)] = 0.;
    
      #if (REAPER_MOMENTS==15)
        for (int lambda=0; lambda<NDIM; lambda++)
        {
          momentsInOrthTetrad[M_UP_UP_UP(mu, nu, lambda)] = 0.;
        }
      #endif
    }
  }
//#pragma omp parallel for 
#if (!GYROAVERAGING)
  for (int iQuad=0; iQuad<NUM_QUAD; iQuad++)
  {
    for (int jQuad=0; jQuad<NUM_QUAD; jQuad++)
    {
      for (int kQuad=0; kQuad<NUM_QUAD; kQuad++)
      {

        REAL t[NDIM], pUpHat[NDIM], pDownHat[NDIM];

        t[1] = quadPts[iQuad];
        t[2] = quadPts[jQuad];
        t[3] = quadPts[kQuad];
      
        for (int i=1; i<NDIM; i++)
        {
          pDownHat[i] = scaleFactor * P_DOWN_FROM_T_FULL(t[i]);
        }

        REAL f;

        computefAndPUpHatUsingOrthTetradPDownHatSpatial
                                  (
                                    pDownHat, geom, elem,
                                    pUpHat, &f
                                  );

        REAL jacobian =   scaleFactor * JACOBIAN_FULL(t[1]) 
                        * scaleFactor * JACOBIAN_FULL(t[2])
                        * scaleFactor * JACOBIAN_FULL(t[3]);
        
        REAL weight =  weights[iQuad]*weights[jQuad]*weights[kQuad]
                     * jacobian/pUpHat[0] * f; 

        for (int mu=0; mu<NDIM; mu++)
        {
          momentsInOrthTetrad[N_UP(mu)] += weight * pUpHat[mu];

          for (int nu=0; nu<NDIM; nu++)
          {
            momentsInOrthTetrad[T_UP_UP(mu, nu)] += 
              weight * pUpHat[mu] * pUpHat[nu];

            #if (REAPER_MOMENTS == 15)
              for (int lambda=0; lambda<NDIM; lambda++)
              {
                momentsInOrthTetrad[M_UP_UP_UP(mu, nu, lambda)] += 
                  weight * pUpHat[mu] * pUpHat[nu] * pUpHat[lambda];
              }
            #endif 

	  }
        }
        
      }
    }
  }
#else
  for (int iQuad=0; iQuad<NUM_QUAD; iQuad++)
  {
    for (int jQuad=0; jQuad<NUM_QUAD_PERP; jQuad++)
    {
      REAL t[NDIM], pUpHat[NDIM], pDownHat[NDIM];

      t[1] = quadPts[iQuad];
      t[2] = quadPtsPerp[jQuad];

      pDownHat[1] = scaleFactor * P_DOWN_FROM_T_FULL(t[1]);
      pDownHat[2] = scaleFactor * P_DOWN_FROM_T_PLUS_HALF(t[2]);
      pDownHat[3] = pDownHat[2];

      REAL f;
      computefAndPUpHatUsingOrthTetradPDownHatSpatial
                                (
                                  pDownHat, geom, elem,
                                  pUpHat, &f
                                );

      REAL jacobian =   scaleFactor * JACOBIAN_FULL(t[1]) 
                      * scaleFactor * JACOBIAN_PLUS_HALF(t[2]);
      
      /* Instead of cartesian coordinates, now integrating in cylindrical coordinates */
      REAL weight =  weights[iQuad]*weightsPerp[jQuad]
                   * jacobian/pUpHat[0] * pUpHat[2] * 2.*M_PI * f; 
      
      for (int mu=0; mu<NDIM; mu++)
      {
        momentsInOrthTetrad[N_UP(mu)] += weight * pUpHat[mu];

        for (int nu=0; nu<NDIM; nu++)
        {
          momentsInOrthTetrad[T_UP_UP(mu, nu)] += 
             weight * pUpHat[mu] * pUpHat[nu];

          #if (REAPER_MOMENTS == 15) 
            for (int lambda=0; lambda<NDIM; lambda++)
            {
              momentsInOrthTetrad[M_UP_UP_UP(mu, nu, lambda)] += 
                weight * pUpHat[mu] * pUpHat[nu] * pUpHat[lambda];
            }
          #endif 

        }
      }
    }
  }
#endif
  for (int mu=0; mu<NDIM; mu++)
  {
    moments[N_UP(mu)] = 0;

    for (int alpha=0; alpha<NDIM; alpha++)
    {
      moments[N_UP(mu)] +=  momentsInOrthTetrad[N_UP(alpha)]
                          * elem->eDownHatUpNoHat[alpha][mu];
    }

    for (int nu=0; nu<NDIM; nu++)
    {
      moments[T_UP_UP(mu, nu)] = 0.;

      for (int alpha=0; alpha<NDIM; alpha++)
      {
        for (int beta=0; beta<NDIM; beta++)
        {
          moments[T_UP_UP(mu, nu)] +=  
                              elem->eDownHatUpNoHat[alpha][mu]
                            * elem->eDownHatUpNoHat[beta][nu]
                            * momentsInOrthTetrad[T_UP_UP(alpha, beta)];
        }
      }

      #if (REAPER_MOMENTS == 15)
      for (int lambda=0; lambda<NDIM; lambda++)
      {
        moments[M_UP_UP_UP(mu, nu, lambda)] = 0.;
        
        for (int alpha=0; alpha<NDIM; alpha++)
        {
          for (int beta=0; beta<NDIM; beta++)
          {
            for (int gamma=0; gamma<NDIM; gamma++)
            {
              moments[M_UP_UP_UP(mu, nu, lambda)] +=
                             elem->eDownHatUpNoHat[alpha][mu]
                           * elem->eDownHatUpNoHat[beta][nu]
                           * elem->eDownHatUpNoHat[gamma][lambda]
                           * momentsInOrthTetrad[M_UP_UP_UP(alpha, beta, gamma)];
            }
          }
        }
      }
      #endif
    }
  }
}

void computefAndPUpHatUsingOrthTetradPDownHatSpatial
(
  REAL pDownHat[NDIM],
  const struct geometry* geom,
  const struct fluidElement* elem,
  REAL pUpHat[NDIM],
  REAL *f
)
{
  /* Calculate pDownHat[0] from normalization of four-momentum */
  pDownHat[0] = -sqrt(1 + pDownHat[1]*pDownHat[1] 
                        + pDownHat[2]*pDownHat[2]
                        + pDownHat[3]*pDownHat[3]);
  
  /* Calculate pUpHat[mu] to evaluate the distribution function */
  pUpHat[0] = -pDownHat[0];
  pUpHat[1] =  pDownHat[1];
  pUpHat[2] =  pDownHat[2];
  pUpHat[3] =  pDownHat[3];
 
  #if (REAPER_MOMENTS==5) 
    /* Maxwell-Juettner in the tetrad frame:
     *
     * f = \alpha exp(a_\hat{0} p^\hat{0}) 
     *
     * where \alpha = \rho/(4 \pi \theta K_2(1/\theta)) 
     *          a_0 = -1 / \theta*/

    *f = elem->primVars[ALPHA] * exp(elem->primVars[A0]*pUpHat[0]);

  #elif (REAPER_MOMENTS==15)
    /* 15 moment Israel-Stewart ansatz in a tetrad frame:
     *
     * f = \alpha exp(a_\hat{\mu} p^\hat{\mu} + b_{\hat{\mu} \hat{\nu}} p^{\hat{\mu} \hat{\nu}} )
     * 
     * We now make a choice that a_{\hat{i}} = 0 in the co-moving tetrad frame. 
     * This corresponds to the Eckart frame. On the other hand if we choose b_{\hat{i}} = 0, 
     * this is the Landau frame. The Eckart frame identifies the 4-velocity with the
     * number flux vector. In this frame, there is a finite heat flux because of non-zero
     * b_{\hat{i}}. On the other hand, the Landau frame identifies the 4-velocity with the heat 
     * flow. 
     *
     * We define bScalar = b_{\hat{\mu} \hat{\nu}} p^{\hat{\mu} \hat{\nu}} */

    REAL bDownmuDownnu[NDIM][NDIM];
    #if (GYROAVERAGING)
      bDownmuDownnu[0][0] = elem->primVars[B00];
      bDownmuDownnu[0][1] = elem->primVars[B01];
      bDownmuDownnu[0][2] = elem->primVars[B02];
      bDownmuDownnu[0][3] = elem->primVars[B02];
      bDownmuDownnu[1][0] = elem->primVars[B01];
      bDownmuDownnu[1][1] = elem->primVars[B11];
      bDownmuDownnu[1][2] = elem->primVars[B12];
      bDownmuDownnu[1][3] = elem->primVars[B12];
      bDownmuDownnu[2][0] = elem->primVars[B02];
      bDownmuDownnu[2][1] = elem->primVars[B12];
      bDownmuDownnu[2][2] = elem->primVars[B22];
      bDownmuDownnu[2][3] = elem->primVars[B22];
      bDownmuDownnu[3][0] = elem->primVars[B02];
      bDownmuDownnu[3][1] = elem->primVars[B12];
      bDownmuDownnu[3][2] = elem->primVars[B22];
      bDownmuDownnu[3][3] = elem->primVars[B22];
    #else 
      bDownmuDownnu[0][0] = elem->primVars[B00];
      bDownmuDownnu[0][1] = elem->primVars[B01];
      bDownmuDownnu[0][2] = elem->primVars[B02];
      bDownmuDownnu[0][3] = elem->primVars[B03];
      bDownmuDownnu[1][0] = elem->primVars[B01];
      bDownmuDownnu[1][1] = elem->primVars[B11];
      bDownmuDownnu[1][2] = elem->primVars[B12];
      bDownmuDownnu[1][3] = elem->primVars[B13];
      bDownmuDownnu[2][0] = elem->primVars[B02];
      bDownmuDownnu[2][1] = elem->primVars[B12];
      bDownmuDownnu[2][2] = elem->primVars[B22];
      bDownmuDownnu[2][3] = elem->primVars[B23];
      bDownmuDownnu[3][0] = elem->primVars[B03];
      bDownmuDownnu[3][1] = elem->primVars[B13];
      bDownmuDownnu[3][2] = elem->primVars[B23];
      bDownmuDownnu[3][3] = elem->primVars[B33];
    #endif

    REAL bScalar = 0.;

    for (int mu=0; mu<NDIM; mu++)
    {
      for (int nu=0; nu<NDIM; nu++)
      {
        bScalar += bDownmuDownnu[mu][nu]*pUpHat[mu]*pUpHat[nu];
      }
    }
    
    /* Diagonal components */
    //bScalar +=    elem->primVars[B00]*pUpHat[0]*pUpHat[0];
    //bScalar +=    elem->primVars[B11]*pUpHat[1]*pUpHat[1];
    //bScalar +=    elem->primVars[B22]*pUpHat[2]*pUpHat[2];
    //bScalar +=    elem->primVars[B33]*pUpHat[3]*pUpHat[3];

    /* Off-diagonal components */
    //bScalar += 2.*elem->primVars[B01]*pUpHat[0]*pUpHat[1];
    //bScalar += 2.*elem->primVars[B02]*pUpHat[0]*pUpHat[2];
    //bScalar += 2.*elem->primVars[B03]*pUpHat[0]*pUpHat[3];
    //bScalar += 2.*elem->primVars[B12]*pUpHat[1]*pUpHat[2];
    //bScalar += 2.*elem->primVars[B13]*pUpHat[1]*pUpHat[3];
    //bScalar += 2.*elem->primVars[B23]*pUpHat[2]*pUpHat[3];

    *f = elem->primVars[ALPHA] * exp(elem->primVars[A0]*pUpHat[0] + bScalar);

    // Now with bscalar as a perturbation
    //*f = elem->primVars[ALPHA] * exp(elem->primVars[A0]*pUpHat[0])*(1. + 0.*bScalar);

  #endif    
  return;
}
