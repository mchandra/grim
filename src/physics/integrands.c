#include "physics.h"

#if (REAPER)
void fixedQuadIntegration(const struct fluidElement *elem,
                          const struct geometry *geom,
                          REAL scaleFactor,
                          REAL *moments,
                          REAL *collisionIntegrals)
{
  /* Using roots of the Legendre polynomaials for integration from [-1, 1].
     We get the quad points and the weights from 
     scipy.special.orthogonal.p_roots() */
  #if (NUM_QUAD==51)
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

  #elif (NUM_QUAD==21)
    REAL quadPts[NUM_QUAD] = \
        {-9.93752171e-01,  -9.67226839e-01,  -9.20099334e-01,
         -8.53363365e-01,  -7.68439963e-01,  -6.67138804e-01,
         -5.51618836e-01,  -4.24342120e-01,  -2.88021317e-01,
         -1.45561854e-01,   1.98918497e-16,   1.45561854e-01,
          2.88021317e-01,   4.24342120e-01,   5.51618836e-01,
          6.67138804e-01,   7.68439963e-01,   8.53363365e-01,
          9.20099334e-01,   9.67226839e-01,   9.93752171e-01};

    REAL weights[NUM_QUAD] = \
        {0.01601723,  0.03695379,  0.05713443,  0.07610011,  0.09344442,
         0.1087973 ,  0.12183142,  0.13226894,  0.13988739,  0.1445244 ,
         0.14608113,  0.1445244 ,  0.13988739,  0.13226894,  0.12183142,
         0.1087973 ,  0.09344442,  0.07610011,  0.05713443,  0.03695379,
         0.01601723};
  #endif

  REAL momentsInOrthTetrad[NUM_ALL_COMPONENTS];

  #if (REAPER_MOMENTS == 15 || REAPER_MOMENTS == 35)
    REAL collisionIntegralsInOrthTetrad[NUM_ALL_COLLISION_INTEGRALS];

    collisionIntegralsInOrthTetrad[0] = 0.;
    collisionIntegralsInOrthTetrad[1] = 0.;
    collisionIntegralsInOrthTetrad[2] = 0.;
    collisionIntegralsInOrthTetrad[3] = 0.;
    collisionIntegralsInOrthTetrad[4] = 0.;
  #endif

  for (int mu=0; mu<NDIM; mu++)
  {
    momentsInOrthTetrad[N_UP(mu)] = 0.;

    for (int nu=0; nu<NDIM; nu++)
    {
      momentsInOrthTetrad[T_UP_UP(mu, nu)] = 0.;
    
      #if (REAPER_MOMENTS==15)
        
        collisionIntegralsInOrthTetrad[COLLISION_INTEGRAL_2(mu, nu)] = 0.;

        for (int lambda=0; lambda<NDIM; lambda++)
        {
          momentsInOrthTetrad[M_UP_UP_UP(mu, nu, lambda)] = 0.;
        }
        
      #elif (REAPER_MOMENTS==35)

        collisionIntegralsInOrthTetrad[COLLISION_INTEGRAL_2(mu, nu)] = 0.;

        for (int lambda=0; lambda<NDIM; lambda++)
        {
          momentsInOrthTetrad[M_UP_UP_UP(mu, nu, lambda)] = 0.;

          collisionIntegralsInOrthTetrad[COLLISION_INTEGRAL_3(mu, nu, lambda)] = 0.;
          
          for (int eta=0; eta<NDIM; eta++)
          {
            momentsInOrthTetrad[R_UP_UP_UP_UP(mu, nu, lambda, eta)] = 0.;
          }
        }
        
      #endif
    }
  }

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

        REAL f, collisionOperator;
        computefAndPUpHat(pDownHat, geom, elem,
                          pUpHat, &f, &collisionOperator
                         );

        REAL jacobian =   scaleFactor * JACOBIAN_FULL(t[1]) 
                        * scaleFactor * JACOBIAN_FULL(t[2])
                        * scaleFactor * JACOBIAN_FULL(t[3]);
        
        REAL weight =  weights[iQuad]*weights[jQuad]*weights[kQuad]
                     * jacobian/pUpHat[0]; 
        
        #if (REAPER_MOMENTS == 15 || REAPER_MOMENTS == 35)
        collisionIntegralsInOrthTetrad[0] += weight * collisionOperator;
        collisionIntegralsInOrthTetrad[1] += weight * collisionOperator * pUpHat[0];
        collisionIntegralsInOrthTetrad[2] += weight * collisionOperator * pUpHat[1];
        collisionIntegralsInOrthTetrad[3] += weight * collisionOperator * pUpHat[2];
        collisionIntegralsInOrthTetrad[4] += weight * collisionOperator * pUpHat[3];  
        #endif

        for (int mu=0; mu<NDIM; mu++)
        {
          momentsInOrthTetrad[N_UP(mu)] += weight * f * pUpHat[mu];

          for (int nu=0; nu<NDIM; nu++)
          {
            momentsInOrthTetrad[T_UP_UP(mu, nu)] += 
              weight * f * pUpHat[mu] * pUpHat[nu];

            #if (REAPER_MOMENTS == 15)

              collisionIntegralsInOrthTetrad[COLLISION_INTEGRAL_2(mu, nu)] +=
                weight * collisionOperator * pUpHat[mu] * pUpHat[nu];

              for (int lambda=0; lambda<NDIM; lambda++)
              {
                momentsInOrthTetrad[M_UP_UP_UP(mu, nu, lambda)] += 
                  weight * f * pUpHat[mu] * pUpHat[nu] * pUpHat[lambda];
              }

            #elif (REAPER_MOMENTS==35)

              collisionIntegralsInOrthTetrad[COLLISION_INTEGRAL_2(mu, nu)] +=
                weight * collisionOperator * pUpHat[mu] * pUpHat[nu];

              for (int lambda=0; lambda<NDIM; lambda++)
              {
                momentsInOrthTetrad[M_UP_UP_UP(mu, nu, lambda)] += 
                  weight * f * pUpHat[mu] * pUpHat[nu] * pUpHat[lambda];

                collisionIntegralsInOrthTetrad[COLLISION_INTEGRAL_3(mu, nu, lambda)] += 
                  weight * collisionOperator * pUpHat[mu] * pUpHat[nu] * pUpHat[lambda];

                for (int eta=0; eta<NDIM; eta++)
                {
                  momentsInOrthTetrad[R_UP_UP_UP_UP(mu, nu, lambda, eta)] += 
                    weight * f * pUpHat[mu] * pUpHat[nu] * pUpHat[lambda] * pUpHat[eta];
                }

              }
            #endif 

	        }
        }
        
      }
    }
  }


  /* Integrals have all been performed in the tetrad frame. Now transform to the
   * coordinate frame */

  collisionIntegrals[0] = collisionIntegralsInOrthTetrad[0];

  for (int mu=0; mu<NDIM; mu++)
  {
    moments[N_UP(mu)] = 0;

    collisionIntegrals[COLLISION_INTEGRAL_1(mu)] = 0;

    for (int alpha=0; alpha<NDIM; alpha++)
    {
      moments[N_UP(mu)] +=  momentsInOrthTetrad[N_UP(alpha)]
                          * elem->eDownHatUpNoHat[alpha][mu];

      collisionIntegrals[COLLISION_INTEGRAL_1(mu)] +=
          collisionIntegralsInOrthTetrad[COLLISION_INTEGRAL_1(alpha)]
        * elem->eDownHatUpNoHat[alpha][mu];
    }

    for (int nu=0; nu<NDIM; nu++)
    {
      moments[T_UP_UP(mu, nu)] = 0.;

      collisionIntegrals[COLLISION_INTEGRAL_2(mu, nu)] = 0;

      for (int alpha=0; alpha<NDIM; alpha++)
      {
        for (int beta=0; beta<NDIM; beta++)
        {
          moments[T_UP_UP(mu, nu)] +=  
                              elem->eDownHatUpNoHat[alpha][mu]
                            * elem->eDownHatUpNoHat[beta][nu]
                            * momentsInOrthTetrad[T_UP_UP(alpha, beta)];

          collisionIntegrals[COLLISION_INTEGRAL_2(mu, nu)] +=
              elem->eDownHatUpNoHat[alpha][mu]
            * elem->eDownHatUpNoHat[beta][nu]
            * collisionIntegralsInOrthTetrad[COLLISION_INTEGRAL_2(alpha, beta)];
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
      #elif (REAPER_MOMENTS == 35)
      for (int lambda=0; lambda<NDIM; lambda++)
      {
        moments[M_UP_UP_UP(mu, nu, lambda)] = 0.;

        collisionIntegrals[COLLISION_INTEGRAL_3(mu, nu, lambda)] = 0;
        
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
              
            collisionIntegrals[COLLISION_INTEGRAL_3(mu, nu, lambda)] +=
                elem->eDownHatUpNoHat[alpha][mu]
              * elem->eDownHatUpNoHat[beta][nu]
              * elem->eDownHatUpNoHat[gamma][lambda]
              * collisionIntegralsInOrthTetrad[COLLISION_INTEGRAL_3(alpha, beta, gamma)];

            }
          }
        }

        for (int eta=0; eta<NDIM; eta++)
        {
          moments[R_UP_UP_UP_UP(mu, nu, lambda, eta)] = 0.;

          for (int alpha=0; alpha<NDIM; alpha++)
          {
            for (int beta=0; beta<NDIM; beta++)
            {
              for (int gamma=0; gamma<NDIM; gamma++)
              {
                for (int delta=0; delta<NDIM; delta++)
                {
                  moments[R_UP_UP_UP_UP(mu, nu, lambda, delta)] += 
                    elem->eDownHatUpNoHat[alpha][mu]
                  * elem->eDownHatUpNoHat[beta][nu]
                  * elem->eDownHatUpNoHat[gamma][lambda]
                  * elem->eDownHatUpNoHat[delta][eta]
                  * momentsInOrthTetrad[R_UP_UP_UP_UP(alpha, beta, gamma, delta)];
                }
              }
            }
          }
        }
      }
      #endif
    }
  }

}

void computefAndPUpHat
(
  REAL pDownHat[NDIM],
  const struct geometry* geom,
  const struct fluidElement* elem,
  REAL pUpHat[NDIM],
  REAL *f,
  REAL *collisionOperator
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

  #elif (REAPER_MOMENTS==15 || REAPER_MOMENTS==35)
    /* Maximum entropy closure in a tetrad frame:
     *
     * f = \alpha exp(aScalar + bScalar + cScalar)
     *
     * where 
     *
     * aScalar = a_\hat{\mu} p^\hat{\mu} 
     * bScalar = b_{\hat{\mu} \hat{\nu}} p^{\hat{\mu}} p^{\hat{\nu}}
     * cScalar = c_{\hat{\mu} \hat{\nu} \hat{\lambda}} p^{\hat{\mu} \hat{\nu} \hat{\lambda}})
     * 
     * We work in the Landau frame i.e. 3 variables (a_1, a_2, a_3) in the above
     * ansatz for f are solved by imposing that f in our tetrad frame has no
     * heat flux i.e.
     *
     * \int d^3p/p^{\hat{0}} p^{\hat{0}} p^{\hat{i} f = 0 */

    REAL bDownDown[NDIM][NDIM];
    /* Note: b_{00} = 0 because it is a redundant variable */
    #if (GYROAVERAGING)
      bDownDown[0][0] = 0.;
      bDownDown[0][1] = elem->primVars[B01];
      bDownDown[0][2] = elem->primVars[B02];
      bDownDown[0][3] = elem->primVars[B02];
      bDownDown[1][0] = elem->primVars[B01];
      bDownDown[1][1] = elem->primVars[B11];
      bDownDown[1][2] = elem->primVars[B12];
      bDownDown[1][3] = elem->primVars[B12];
      bDownDown[2][0] = elem->primVars[B02];
      bDownDown[2][1] = elem->primVars[B12];
      bDownDown[2][2] = elem->primVars[B22];
      bDownDown[2][3] = elem->primVars[B22];
      bDownDown[3][0] = elem->primVars[B02];
      bDownDown[3][1] = elem->primVars[B12];
      bDownDown[3][2] = elem->primVars[B22];
      bDownDown[3][3] = elem->primVars[B22];
    #else 
      bDownDown[0][0] = 0.;
      bDownDown[0][1] = elem->primVars[B01];
      bDownDown[0][2] = elem->primVars[B02];
      bDownDown[0][3] = elem->primVars[B03];
      bDownDown[1][0] = elem->primVars[B01];
      bDownDown[1][1] = elem->primVars[B11];
      bDownDown[1][2] = elem->primVars[B12];
      bDownDown[1][3] = elem->primVars[B13];
      bDownDown[2][0] = elem->primVars[B02];
      bDownDown[2][1] = elem->primVars[B12];
      bDownDown[2][2] = elem->primVars[B22];
      bDownDown[2][3] = elem->primVars[B23];
      bDownDown[3][0] = elem->primVars[B03];
      bDownDown[3][1] = elem->primVars[B13];
      bDownDown[3][2] = elem->primVars[B23];
      bDownDown[3][3] = elem->primVars[B33];
    #endif

    #if (REAPER_MOMENTS==35)
    /* The following 4 variables are redundant. For example, C_001 simply adds to
     * a_1, C_111, C_221 and C_331 because p0^2 ~ 1 + p1^2 + p2^2 + p3^2 */
    
      REAL cDownDownDown[NDIM][NDIM][NDIM];

      cDownDownDown[0][0][0] = 0.;
      cDownDownDown[0][0][1] = 0.;
      cDownDownDown[0][0][2] = 0.;
      cDownDownDown[0][0][3] = 0.;
    
      cDownDownDown[0][1][0] = 0.;
      cDownDownDown[0][1][1] = elem->primVars[C011];
      cDownDownDown[0][1][2] = elem->primVars[C012];
      cDownDownDown[0][1][3] = elem->primVars[C013];

      cDownDownDown[0][2][0] = 0.;
      cDownDownDown[0][2][1] = elem->primVars[C012];
      cDownDownDown[0][2][2] = elem->primVars[C022];
      cDownDownDown[0][2][3] = elem->primVars[C023];

      cDownDownDown[0][3][0] = 0.;
      cDownDownDown[0][3][1] = elem->primVars[C013];
      cDownDownDown[0][3][2] = elem->primVars[C023];
      cDownDownDown[0][3][3] = elem->primVars[C033];

      cDownDownDown[1][0][0] = 0.;
      cDownDownDown[1][0][1] = elem->primVars[C011];
      cDownDownDown[1][0][2] = elem->primVars[C012];
      cDownDownDown[1][0][3] = elem->primVars[C013];

      cDownDownDown[1][1][0] = elem->primVars[C011];
      cDownDownDown[1][1][1] = elem->primVars[C111];
      cDownDownDown[1][1][2] = elem->primVars[C112];
      cDownDownDown[1][1][3] = elem->primVars[C113];

      cDownDownDown[1][2][0] = elem->primVars[C012];
      cDownDownDown[1][2][1] = elem->primVars[C112];
      cDownDownDown[1][2][2] = elem->primVars[C122];
      cDownDownDown[1][2][3] = elem->primVars[C123];

      cDownDownDown[1][3][0] = elem->primVars[C013];
      cDownDownDown[1][3][1] = elem->primVars[C113];
      cDownDownDown[1][3][2] = elem->primVars[C123];
      cDownDownDown[1][3][3] = elem->primVars[C133];

      cDownDownDown[2][0][0] = 0.;
      cDownDownDown[2][0][1] = elem->primVars[C012];
      cDownDownDown[2][0][2] = elem->primVars[C022];
      cDownDownDown[2][0][3] = elem->primVars[C023];

      cDownDownDown[2][1][0] = elem->primVars[C012];
      cDownDownDown[2][1][1] = elem->primVars[C112];
      cDownDownDown[2][1][2] = elem->primVars[C122];
      cDownDownDown[2][1][3] = elem->primVars[C123];

      cDownDownDown[2][2][0] = elem->primVars[C022];
      cDownDownDown[2][2][1] = elem->primVars[C122];
      cDownDownDown[2][2][2] = elem->primVars[C222];
      cDownDownDown[2][2][3] = elem->primVars[C223];

      cDownDownDown[2][3][0] = elem->primVars[C023];
      cDownDownDown[2][3][1] = elem->primVars[C123];
      cDownDownDown[2][3][2] = elem->primVars[C223];
      cDownDownDown[2][3][3] = elem->primVars[C233];

      cDownDownDown[3][0][0] = 0.;
      cDownDownDown[3][0][1] = elem->primVars[C013];
      cDownDownDown[3][0][2] = elem->primVars[C023];
      cDownDownDown[3][0][3] = elem->primVars[C033];

      cDownDownDown[3][1][0] = elem->primVars[C013];
      cDownDownDown[3][1][1] = elem->primVars[C113];
      cDownDownDown[3][1][2] = elem->primVars[C123];
      cDownDownDown[3][1][3] = elem->primVars[C133];

      cDownDownDown[3][2][0] = elem->primVars[C023];
      cDownDownDown[3][2][1] = elem->primVars[C123];
      cDownDownDown[3][2][2] = elem->primVars[C223];
      cDownDownDown[3][2][3] = elem->primVars[C233];

      cDownDownDown[3][3][0] = elem->primVars[C033];
      cDownDownDown[3][3][1] = elem->primVars[C133];
      cDownDownDown[3][3][2] = elem->primVars[C233];
      cDownDownDown[3][3][3] = elem->primVars[C333];
    #endif

    REAL bScalar = 0.;
    REAL cScalar = 0.;

    for (int mu=0; mu<NDIM; mu++)
    {
      for (int nu=0; nu<NDIM; nu++)
      {
        bScalar += bDownDown[mu][nu]*pUpHat[mu]*pUpHat[nu];

        #if (REAPER_MOMENTS==35)
          for (int lambda=0; lambda<NDIM; lambda++)
          {
            cScalar += 
              cDownDownDown[mu][nu][lambda]*pUpHat[mu]*pUpHat[nu]*pUpHat[lambda];
          }
        #endif
      }
    }
    
    *f =   elem->primVars[ALPHA] 
         * exp(  elem->primVars[A0]*pUpHat[0] + bScalar + cScalar
               + elem->primVars[F0_A1]*pUpHat[1] 
               + elem->primVars[F0_A2]*pUpHat[2] 
               + elem->primVars[F0_A3]*pUpHat[3] 
              );
               
    REAL f0 =  elem->primVars[F0_ALPHA]
             * exp(elem->primVars[F0_A0]*pUpHat[0]);

    REAL tau = 0.01;

    /* Anderson-Witting collision operator */
    *collisionOperator = - pUpHat[0] * (*f - f0)/tau;

  #endif    

  return;
}

#endif /* Only compile when REAPER==ON */
