#include "physics.h"

void fixedQuadIntegration5Moments(const struct fluidElement *elem,
                                  const struct geometry *geom,
                                  REAL scaleFactor,
                                  REAL *moments)
{
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

  REAL momentsInOrthTetrad[NUM_ALL_COMPONENTS];

  for (int mu=0; mu<NDIM; mu++)
  {
    momentsInOrthTetrad[N_UP(mu)] = 0.;

    for (int nu=0; nu<NDIM; nu++)
    {
      momentsInOrthTetrad[T_UP_UP(mu, nu)] = 0.;
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
          }
        }
        
      }
    }
  }


  for (int mu=0; mu<NDIM; mu++)
  {
    moments[N_UP(mu)] = 0;

    for (int nu=0; nu<NDIM; nu++)
    {
      moments[N_UP(mu)] +=  momentsInOrthTetrad[N_UP(nu)]
                           *elem->eDownHatUpNoHat[nu][mu];
    }
  }

  for (int mu=0; mu<NDIM; mu++)
  {
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
  
  *f = elem->primVars[ALPHA] * exp(pUpHat[0] * elem->primVars[A0]);
    
  return;
}
