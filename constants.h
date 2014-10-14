#define COMPUTE_DIM 2
#define NDIM 4
#define N1 128
#define N2 128
#define NG 2

#define REAL double
#define TILE_SIZE_X1 8
#define TILE_SIZE_X2 8

#define MKS         (0)
#define MINKOWSKI   (1)
#define GEOMETRY    (MKS)

#if (GEOMETRY==MKS)
  #define R0 0.
  #define R_IN .8*(1. + sqrt(1. - A_SPIN*A_SPIN))
  #define R_OUT 40.
  #define X1_START log(R_IN - R0)
  #define X2_START 1e-3
  #define DX1 (log((R_OUT - R0)/(R_IN - R0))/(REAL)N1)
  #define DX2 ((1.-2.*X2_START)/(REAL)N2)
#elif (GEOMETRY==MINKOWSKI)
  #define X1_START (0.)
  #define X2_START (0.)
  #define X1_END   (1.)
  #define X2_END   (1.)
  #define DX1 ((X1_END - X1_START)/(REAL)N1)
  #define DX2 ((X2_END - X2_START)/(REAL)N2)
#endif

/* Boundary mnemonics */
#define OUTFLOW   (0)
#define MIRROR    (1)
#define DIRICHLET (2)
#define PERIODIC  (3)

#define LEFT_BOUNDARY   (OUTFLOW)
#define RIGHT_BOUNDARY  (OUTFLOW)
#define BOTTOM_BOUNDARY (MIRROR)
#define TOP_BOUNDARY    (MIRROR)
#define INFLOW_CHECK    (1)

#define A_SPIN 0.9375
#define M 1.
#define R_MIN 6.
#define R_MAX 12.
#define H_SLOPE 0.3
#define DT 0.005
#define DT_DUMP 1.
#define KAPPA 1e-3
#define BETA 1e2
#define ADIABATIC_INDEX (4/3.)
#define RHO_MIN (1e-4)
#define U_MIN (1e-5)
#define RHO_MIN_LIMIT (1e-15)
#define U_MIN_LIMIT (1e-15)
#define GAMMA_MAX (5.)
#define TAU_R_SAFETY_FACTOR (1.2)
#define PHI (5.)
#define CONDUCTION (1)
#define RESTART (0)

#define EPS (1e-5)

#define RHO 0
#define UU 1
#define U1 2
#define U2 3
#define U3 4
#define B1 5
#define B2 6
#define B3 7

#if (CONDUCTION)
  #define FF 8
  #define DOF 9
#else
  #define DOF 8
#endif


/* iTile, jTile have ranges [-NG, TILE_SIZE+NG) */
#define INDEX_LOCAL(iTile,jTile,var) (iTile+NG + \
                                      (TILE_SIZE_X1+2*NG)*(jTile+NG + \
                                      (TILE_SIZE_X2+2*NG)*(var)))

/* i, j have ranges [0, N1), [0, N2) */
#define INDEX_GLOBAL(i,j,var) (var + DOF*((i)+(N1)*(j)))

#define INDEX_GLOBAL_WITH_NG(i,j,var) (var + DOF*((i+NG)+(N1+2*NG)*(j+NG)))

#define i_TO_X1_CENTER(i) (X1_START + (i + 0.5)*DX1)
#define j_TO_X2_CENTER(j) (X2_START + (j + 0.5)*DX2)
#define i_TO_X1_FACE(i) (X1_START + (i)*DX1)
#define j_TO_X2_FACE(j) (X2_START + (j)*DX2)

void BLCoords(REAL* r, REAL* theta,
              const REAL X1, const REAL X2)
{
#if (GEOMETRY==MKS)
    *r = exp(X1) + R0;
    *theta = M_PI*(X2) + ((1 - H_SLOPE)/2.)*sin(2.*M_PI*(X2));
#elif (GEOMETRY==MINKOWSKI)
    *r = X1; *theta = X2;
#endif
}


void gCovCalc(REAL gcov[NDIM][NDIM],
              const REAL X1, const REAL X2)
{
#if (GEOMETRY==MKS)
    REAL r, theta;
    BLCoords(&r, &theta, X1, X2);

    REAL rhoSqr = r*r + A_SPIN*A_SPIN*cos(theta)*cos(theta);
    REAL rFactor = r - R0;
    REAL hFactor = M_PI + (1. - H_SLOPE)*M_PI*cos(2.*M_PI*(X2));

    gcov[0][0] = -1. + 2.*r/rhoSqr;
    gcov[0][1] = (2.*r/rhoSqr)*rFactor;
    gcov[0][2] = 0.;
    gcov[0][3] = -2.*A_SPIN*r*sin(theta)*sin(theta)/rhoSqr;

    gcov[1][0] = gcov[0][1];
    gcov[1][1] = (1. + 2.*r/rhoSqr)*rFactor*rFactor;
    gcov[1][2] = 0.;
    gcov[1][3] = -A_SPIN*sin(theta)*sin(theta)*(1.+
                2.*r/rhoSqr)*rFactor;

    gcov[2][0] = gcov[0][2];
    gcov[2][1] = gcov[1][2];
    gcov[2][2] = rhoSqr*hFactor*hFactor;
    gcov[2][3] = 0.;


    gcov[3][0] = gcov[0][3];
    gcov[3][1] = gcov[1][3];
    gcov[3][2] = gcov[2][3];
    gcov[3][3] = sin(theta)*sin(theta)*(rhoSqr +
                            A_SPIN*A_SPIN*sin(theta)*sin(theta)*\
                            (1. + 2.*r/rhoSqr) );

#elif (GEOMETRY==MINKOWSKI)

    gcov[0][0] = -1.;
    gcov[0][1] = 0.;
    gcov[0][2] = 0.;
    gcov[0][3] = 0.;

    gcov[1][0] = 0.;
    gcov[1][1] = 1.;
    gcov[1][2] = 0.;
    gcov[1][3] = 0.;

    gcov[2][0] = 0.;
    gcov[2][1] = 0.;
    gcov[2][2] = 1.;
    gcov[2][3] = 0.;

    gcov[3][0] = 0.;
    gcov[3][1] = 0.;
    gcov[3][2] = 0.;
    gcov[3][3] = 1.;

#endif
}

void gConCalc(REAL gcon[NDIM][NDIM],
              const REAL gcov[NDIM][NDIM],
              const REAL gdet)
{

    gcon[0][0] = 
        (gcov[1][1]*gcov[2][2]*gcov[3][3] + 
         gcov[1][2]*gcov[2][3]*gcov[3][1] + 
         gcov[1][3]*gcov[2][1]*gcov[3][2] - 
         gcov[1][1]*gcov[2][3]*gcov[3][2] - 
         gcov[1][2]*gcov[2][1]*gcov[3][3] - 
         gcov[1][3]*gcov[2][2]*gcov[3][1])/gdet;

    gcon[0][1] = 
        (gcov[0][1]*gcov[2][3]*gcov[3][2] + 
         gcov[0][2]*gcov[2][1]*gcov[3][3] + 
         gcov[0][3]*gcov[2][2]*gcov[3][1] - 
         gcov[0][1]*gcov[2][2]*gcov[3][3] - 
         gcov[0][2]*gcov[2][3]*gcov[3][1] - 
         gcov[0][3]*gcov[2][1]*gcov[3][2])/gdet;

    gcon[0][2] = 
        (gcov[0][1]*gcov[1][2]*gcov[3][3] + 
         gcov[0][2]*gcov[1][3]*gcov[3][1] + 
         gcov[0][3]*gcov[1][1]*gcov[3][2] - 
         gcov[0][1]*gcov[1][3]*gcov[3][2] - 
         gcov[0][2]*gcov[1][1]*gcov[3][3] - 
         gcov[0][3]*gcov[1][2]*gcov[3][1])/gdet;

    gcon[0][3] = 
        (gcov[0][1]*gcov[1][3]*gcov[2][2] + 
         gcov[0][2]*gcov[1][1]*gcov[2][3] + 
         gcov[0][3]*gcov[1][2]*gcov[2][1] - 
         gcov[0][1]*gcov[1][2]*gcov[2][3] - 
         gcov[0][2]*gcov[1][3]*gcov[2][1] - 
         gcov[0][3]*gcov[1][1]*gcov[2][2])/gdet;

    gcon[1][0] = gcon[0][1];
    
    gcon[1][1] = 
        (gcov[0][0]*gcov[2][2]*gcov[3][3] + 
         gcov[0][2]*gcov[2][3]*gcov[3][0] + 
         gcov[0][3]*gcov[2][0]*gcov[3][2] - 
         gcov[0][0]*gcov[2][3]*gcov[3][2] - 
         gcov[0][2]*gcov[2][0]*gcov[3][3] - 
         gcov[0][3]*gcov[2][2]*gcov[3][0])/gdet;

    gcon[1][2] = 
        (gcov[0][0]*gcov[1][3]*gcov[3][2] + 
         gcov[0][2]*gcov[1][0]*gcov[3][3] + 
         gcov[0][3]*gcov[1][2]*gcov[3][0] - 
         gcov[0][0]*gcov[1][2]*gcov[3][3] - 
         gcov[0][2]*gcov[1][3]*gcov[3][0] - 
         gcov[0][3]*gcov[1][0]*gcov[3][2])/gdet;

    gcon[1][3] = 
        (gcov[0][0]*gcov[1][2]*gcov[2][3] + 
         gcov[0][2]*gcov[1][3]*gcov[2][0] + 
         gcov[0][3]*gcov[1][0]*gcov[2][2] - 
         gcov[0][0]*gcov[1][3]*gcov[2][2] - 
         gcov[0][2]*gcov[1][0]*gcov[2][3] - 
         gcov[0][3]*gcov[1][2]*gcov[2][0])/gdet;

    gcon[2][0] = gcon[0][2];
    gcon[2][1] = gcon[1][2];

    gcon[2][2] =
        (gcov[0][0]*gcov[1][1]*gcov[3][3] + 
         gcov[0][1]*gcov[1][3]*gcov[3][0] + 
         gcov[0][3]*gcov[1][0]*gcov[3][1] - 
         gcov[0][0]*gcov[1][3]*gcov[3][1] - 
         gcov[0][1]*gcov[1][0]*gcov[3][3] - 
         gcov[0][3]*gcov[1][1]*gcov[3][0])/gdet;

    gcon[2][3] =
        (gcov[0][0]*gcov[1][3]*gcov[2][1] + 
         gcov[0][1]*gcov[1][0]*gcov[2][3] + 
         gcov[0][3]*gcov[1][1]*gcov[2][0] - 
         gcov[0][0]*gcov[1][1]*gcov[2][3] - 
         gcov[0][1]*gcov[1][3]*gcov[2][0] - 
         gcov[0][3]*gcov[1][0]*gcov[2][1])/gdet;

    gcon[3][0] = gcon[0][3];
    gcon[3][1] = gcon[1][3];
    gcon[3][2] = gcon[2][3];

    gcon[3][3] =
        (gcov[0][0]*gcov[1][1]*gcov[2][2] + 
         gcov[0][1]*gcov[1][2]*gcov[2][0] + 
         gcov[0][2]*gcov[1][0]*gcov[2][1] - 
         gcov[0][0]*gcov[1][2]*gcov[2][1] - 
         gcov[0][1]*gcov[1][0]*gcov[2][2] - 
         gcov[0][2]*gcov[1][1]*gcov[2][0])/gdet;

}

void gDetCalc(REAL* gdet,
              const REAL gcov[NDIM][NDIM])
{
    *gdet = 
        gcov[0][0]*gcov[1][1]*gcov[2][2]*gcov[3][3] + 
        gcov[0][0]*gcov[1][2]*gcov[2][3]*gcov[3][1] + 
        gcov[0][0]*gcov[1][3]*gcov[2][1]*gcov[3][2] + 
        gcov[0][1]*gcov[1][0]*gcov[2][3]*gcov[3][2] + 
        gcov[0][1]*gcov[1][2]*gcov[2][0]*gcov[3][3] + 
        gcov[0][1]*gcov[1][3]*gcov[2][2]*gcov[3][0] + 
        gcov[0][2]*gcov[1][0]*gcov[2][1]*gcov[3][3] + 
        gcov[0][2]*gcov[1][1]*gcov[2][3]*gcov[3][0] + 
        gcov[0][2]*gcov[1][3]*gcov[2][0]*gcov[3][1] + 
        gcov[0][3]*gcov[1][0]*gcov[2][2]*gcov[3][1] + 
        gcov[0][3]*gcov[1][1]*gcov[2][0]*gcov[3][2] + 
        gcov[0][3]*gcov[1][2]*gcov[2][1]*gcov[3][0] - 
        gcov[0][0]*gcov[1][1]*gcov[2][3]*gcov[3][2] - 
        gcov[0][0]*gcov[1][2]*gcov[2][1]*gcov[3][3] - 
        gcov[0][0]*gcov[1][3]*gcov[2][2]*gcov[3][1] - 
        gcov[0][1]*gcov[1][0]*gcov[2][2]*gcov[3][3] - 
        gcov[0][1]*gcov[1][2]*gcov[2][3]*gcov[3][0] - 
        gcov[0][1]*gcov[1][3]*gcov[2][0]*gcov[3][2] - 
        gcov[0][2]*gcov[1][0]*gcov[2][3]*gcov[3][1] - 
        gcov[0][2]*gcov[1][1]*gcov[2][0]*gcov[3][3] - 
        gcov[0][2]*gcov[1][3]*gcov[2][1]*gcov[3][0] - 
        gcov[0][3]*gcov[1][0]*gcov[2][1]*gcov[3][2] - 
        gcov[0][3]*gcov[1][1]*gcov[2][2]*gcov[3][0] - 
        gcov[0][3]*gcov[1][2]*gcov[2][0]*gcov[3][1];
}

void alphaCalc(REAL* alpha,
               const REAL gcon[NDIM][NDIM])
{
    *alpha = 1./sqrt(-gcon[0][0]);
}

void gammaCalc(REAL* gamma,
               const REAL var[DOF],
               const REAL gcov[NDIM][NDIM])
{
    *gamma = 
        sqrt(1 + gcov[1][1]*var[U1]*var[U1] + 
                 gcov[2][2]*var[U2]*var[U2] + 
                 gcov[3][3]*var[U3]*var[U3] + 
              2*(gcov[1][2]*var[U1]*var[U2] + 
                 gcov[1][3]*var[U1]*var[U3] + 
                 gcov[2][3]*var[U2]*var[U3]));
}

void dgammaCalc_dt(REAL* dgamma_dt,
                   const REAL gamma,
                   const REAL var[DOF],
                   const REAL dvar_dt[DOF],
                   const REAL gcov[NDIM][NDIM])
{
    *dgamma_dt = 
        ((gcov[1][1]*var[U1]*dvar_dt[U1] + 
          gcov[2][2]*var[U2]*dvar_dt[U2] + 
          gcov[3][3]*var[U3]*dvar_dt[U3])+ 
         (gcov[1][2]*dvar_dt[U1]*var[U2] + 
          gcov[1][2]*var[U1]*dvar_dt[U2] + 
          gcov[1][3]*dvar_dt[U1]*var[U3] + 
          gcov[1][3]*var[U1]*dvar_dt[U3] + 
          gcov[2][3]*dvar_dt[U2]*var[U3] + 
          gcov[2][3]*var[U2]*dvar_dt[U3]))/gamma;
}

void uconCalc(REAL ucon[NDIM],
              const REAL gamma,
              const REAL alpha,
              const REAL var[DOF],
              const REAL gcon[NDIM][NDIM])
{
    ucon[0] = gamma/alpha;
    ucon[1] = var[U1] - gamma*gcon[0][1]*alpha;
    ucon[2] = var[U2] - gamma*gcon[0][2]*alpha;
    ucon[3] = var[U3] - gamma*gcon[0][3]*alpha;
}

void duconCalc_dt(REAL ducon_dt[NDIM],
                  const REAL dgamma_dt,
                  const REAL alpha,
                  const REAL dvar_dt[DOF],
                  const REAL gcon[NDIM][NDIM])
{
    ducon_dt[0] = dgamma_dt/alpha;
    ducon_dt[1] = dvar_dt[U1] - dgamma_dt*gcon[0][1]*alpha;
    ducon_dt[2] = dvar_dt[U2] - dgamma_dt*gcon[0][2]*alpha;
    ducon_dt[3] = dvar_dt[U3] - dgamma_dt*gcon[0][3]*alpha;
}

void covFromCon(REAL cov[NDIM],
                const REAL con[NDIM],
                const REAL gcov[NDIM][NDIM])
{
    cov[0] = gcov[0][0]*con[0] + gcov[0][1]*con[1] +
             gcov[0][2]*con[2] + gcov[0][3]*con[3];

    cov[1] = gcov[1][0]*con[0] + gcov[1][1]*con[1] +
             gcov[1][2]*con[2] + gcov[1][3]*con[3];

    cov[2] = gcov[2][0]*con[0] + gcov[2][1]*con[1] +
             gcov[2][2]*con[2] + gcov[2][3]*con[3];

    cov[3] = gcov[3][0]*con[0] + gcov[3][1]*con[1] +
             gcov[3][2]*con[2] + gcov[3][3]*con[3];
}

void conFromCov(REAL con[NDIM],
                const REAL cov[NDIM],
                const REAL gcon[NDIM][NDIM])
{
    con[0] = gcon[0][0]*cov[0] + gcon[0][1]*cov[1] +
             gcon[0][2]*cov[2] + gcon[0][3]*cov[3];

    con[1] = gcon[1][0]*cov[0] + gcon[1][1]*cov[1] +
             gcon[1][2]*cov[2] + gcon[1][3]*cov[3];

    con[2] = gcon[2][0]*cov[0] + gcon[2][1]*cov[1] +
             gcon[2][2]*cov[2] + gcon[2][3]*cov[3];

    con[3] = gcon[3][0]*cov[0] + gcon[3][1]*cov[1] +
             gcon[3][2]*cov[2] + gcon[3][3]*cov[3];

}

void conDotCov(REAL* ans,
               const REAL con[NDIM],
               const REAL cov[NDIM])
{
    *ans = con[0]*cov[0] + con[1]*cov[1] +
           con[2]*cov[2] + con[3]*cov[3];
}

void bconCalc(REAL bcon[NDIM],
              const REAL var[DOF],
              const REAL ucon[NDIM],
              const REAL ucov[NDIM])
{
    bcon[0] = var[B1]*ucov[1]+ var[B2]*ucov[2]+ var[B3]*ucov[3];
    
    bcon[1] = (var[B1] + bcon[0]*ucon[1])/ucon[0];
    bcon[2] = (var[B2] + bcon[0]*ucon[2])/ucon[0];
    bcon[3] = (var[B3] + bcon[0]*ucon[3])/ucon[0];
}

void dbconCalc_dt(REAL dbcon_dt[NDIM],
                  const REAL ucon[NDIM],
                  const REAL ducon_dt[NDIM],
                  const REAL ucov[NDIM],
                  const REAL ducov_dt[NDIM],
                  const REAL bcon[NDIM],
                  const REAL var[DOF],
                  const REAL dvar_dt[DOF])
{

    dbcon_dt[0] = 
        (dvar_dt[B1]*ucov[1] + dvar_dt[B2]*ucov[2] + dvar_dt[B3]*ucov[3] +
         var[B1]*ducov_dt[1] + var[B2]*ducov_dt[2] + var[B3]*ducov_dt[3]);

    dbcon_dt[1] =
        (-(var[B1] + bcon[0]*ucon[1])*ducon_dt[0]/(ucon[0]*ucon[0]) +
         (dvar_dt[B1] + bcon[0]*ducon_dt[1] + dbcon_dt[0]*ucon[1])/ucon[0]);

    dbcon_dt[2] =
        (-(var[B2] + bcon[0]*ucon[2])*ducon_dt[0]/(ucon[0]*ucon[0]) +
         (dvar_dt[B2] + bcon[0]*ducon_dt[2] + dbcon_dt[0]*ucon[2])/ucon[0]);

    dbcon_dt[3] =
        (-(var[B3] + bcon[0]*ucon[3])*ducon_dt[0]/(ucon[0]*ucon[0]) +
         (dvar_dt[B3] + bcon[0]*ducon_dt[3] + dbcon_dt[0]*ucon[3])/ucon[0]);
}

void bSqrCalc(REAL* bsqr,
              const REAL bcon[NDIM],
              const REAL bcov[NDIM])
{
    *bsqr = bcon[0]*bcov[0] + bcon[1]*bcov[1] +
            bcon[2]*bcov[2] + bcon[3]*bcov[3];

    if (*bsqr < 1e-20) *bsqr = 1e-20;
}

REAL SlopeLim(REAL y1, REAL y2, REAL y3)
{
  REAL Dqm = 2. * (y2 - y1);
	REAL Dqp = 2. * (y3 - y2);
	REAL Dqc = 0.5 * (y3 - y1);
	REAL s = Dqm * Dqp;
	if (s <= 0.) {
		return 0.;
    }
	else {
		if (fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
			return (Dqm);
		else if (fabs(Dqp) < fabs(Dqc))
			return (Dqp);
		else
			return (Dqc);
	}
}

#ifdef OPENCL
void mhdCalc(REAL mhd[NDIM][NDIM],
             const REAL var[DOF],
             const REAL ucon[NDIM],
             const REAL ucov[NDIM],
             const REAL bcon[NDIM],
             const REAL bcov[NDIM],
             const REAL dvars_dt[DOF],
             REAL const __local *primTile)
{
    REAL P = (ADIABATIC_INDEX - 1.)*var[UU];
    REAL bsqr;
    bSqrCalc(&bsqr, bcon, bcov);

#define DELTA(mu, nu) (mu==nu ? 1 : 0)

    for (int mu=0; mu<NDIM; mu++)
        for (int nu=0; nu<NDIM; nu++) {
#if(CONDUCTION)
            mhd[mu][nu] = (var[RHO] + var[UU] + P + bsqr)*ucon[mu]*ucov[nu] +
                          (P + 0.5*bsqr)*DELTA(mu, nu) - bcon[mu]*bcov[nu] +
                          (var[FF]*bcon[mu]/sqrt(bsqr))*ucov[nu] + 
                           ucon[mu]*(var[FF]*bcov[nu]/sqrt(bsqr));

#else
            mhd[mu][nu] = (var[RHO] + var[UU] + P + bsqr)*ucon[mu]*ucov[nu] +
                          (P + 0.5*bsqr)*DELTA(mu, nu) - bcon[mu]*bcov[nu];
#endif /* Conduction */
        }

#undef DELTA
}
#endif /* OPENCL */

#ifdef OPENCL
void addSources(REAL dU_dt[DOF],
                REAL ucon[NDIM],
                REAL ucov[NDIM],
                REAL bcon[NDIM],
                REAL bcov[NDIM],
                REAL gcon[NDIM][NDIM],
                REAL gcov[NDIM][NDIM],
                REAL mhd[NDIM][NDIM],
                REAL primVars[DOF],
                REAL g, REAL gdet, REAL alpha,
                REAL X1, REAL X2,
                REAL const __local *primTile,
                REAL dvars_dt[DOF],
                REAL ducon_dt[NDIM],
                const int i, const int j,
                const int iTile, const int jTile)

{
    REAL gcovh[NDIM][NDIM], gcovl[NDIM][NDIM];
    REAL conntmp[NDIM][NDIM][NDIM], conn[NDIM][NDIM][NDIM];
    REAL Xl[NDIM], Xh[NDIM];

    for (int k = 0; k < NDIM; k++) {
        Xl[0] = 0.; Xl[1] = X1; Xl[2] = X2; Xl[3] = 0.;
        Xh[0] = 0.; Xh[1] = X1; Xh[2] = X2; Xh[3] = 0.;
        Xl[k] = Xl[k] - EPS;
        Xh[k] = Xh[k] + EPS;
        gCovCalc(gcovh, Xh[1], Xh[2]);
        gCovCalc(gcovl, Xl[1], Xl[2]);
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j < NDIM; j++) {
				conn[i][j][k] = (gcovh[i][j] - gcovl[i][j])/(Xh[k] - Xl[k]);
            }
        }
    }

	/* now rearrange to find \Gamma_{ijk} */
	for (int i = 0; i < NDIM; i++)
		for (int j = 0; j < NDIM; j++)
			for (int k = 0; k < NDIM; k++)
			    conntmp[i][j][k] =
				    0.5 * (conn[j][i][k] + conn[k][i][j] -
					   conn[k][j][i]);

	/* finally, raise index */
	for (int i = 0; i < NDIM; i++)
		for (int j = 0; j < NDIM; j++)
			for (int k = 0; k < NDIM; k++) {
				conn[i][j][k] = 0.;
				for (int l = 0; l < NDIM; l++)
					conn[i][j][k] += gcon[i][l]*conntmp[l][j][k];
			}
    
    for (int j=0; j<NDIM; j++)
        for (int k=0; k<NDIM; k++) {
            dU_dt[UU] = dU_dt[UU] - g*(mhd[j][k]*conn[k][0][j]);
            dU_dt[U1] = dU_dt[U1] - g*(mhd[j][k]*conn[k][1][j]);
            dU_dt[U2] = dU_dt[U2] - g*(mhd[j][k]*conn[k][2][j]);
            dU_dt[U3] = dU_dt[U3] - g*(mhd[j][k]*conn[k][3][j]);

        }

#if(CONDUCTION)
    REAL left, center, right;
    REAL ducon, dT_dt, dT_dX1, dT_dX2, cs, F0, gamma;

    REAL uconCenter[NDIM];
    uconCenter[0] = ucon[0];
    uconCenter[1] = ucon[1];
    uconCenter[2] = ucon[2];
    uconCenter[3] = ucon[3];

    /* Compute -F*(g*u^\mu),\mu - (1)
               = -F*( d(g*u^0)/dt + d(g*u^1)/dX1 + d(g*u^2)/dX2)
     */

    center = g*uconCenter[1];

    X1 = i_TO_X1_CENTER(i+1); X2 = j_TO_X2_CENTER(j);
    for (int var=0; var<DOF; var++) {
        primVars[var] = primTile[INDEX_LOCAL(iTile+1, jTile, var)];
    }
    gCovCalc(gcov, X1, X2);
    gDetCalc(&gdet, gcov);
    gConCalc(gcon, gcov, gdet);
    gammaCalc(&gamma, primVars, gcov);
    alphaCalc(&alpha, gcon);
    uconCalc(ucon, gamma, alpha, primVars, gcon);

    right = sqrt(-gdet)*ucon[1];

    X1 = i_TO_X1_CENTER(i-1); X2 = j_TO_X2_CENTER(j);
    for (int var=0; var<DOF; var++) {
        primVars[var] = primTile[INDEX_LOCAL(iTile-1, jTile, var)];
    }
    gCovCalc(gcov, X1, X2);
    gDetCalc(&gdet, gcov);
    gConCalc(gcon, gcov, gdet);
    gammaCalc(&gamma, primVars, gcov);
    alphaCalc(&alpha, gcon);
    uconCalc(ucon, gamma, alpha, primVars, gcon);

    left = sqrt(-gdet)*ucon[1];

    ducon = SlopeLim(left, center, right);

    dU_dt[FF] += -primTile[INDEX_LOCAL(iTile, jTile, FF)]*ducon/DX1;

    center = g*uconCenter[2];

    X1 = i_TO_X1_CENTER(i); X2 = j_TO_X2_CENTER(j+1);
    for (int var=0; var<DOF; var++) {
        primVars[var] = primTile[INDEX_LOCAL(iTile, jTile+1, var)];
    }
    gCovCalc(gcov, X1, X2);
    gDetCalc(&gdet, gcov);
    gConCalc(gcon, gcov, gdet);
    gammaCalc(&gamma, primVars, gcov);
    alphaCalc(&alpha, gcon);
    uconCalc(ucon, gamma, alpha, primVars, gcon);

    right = sqrt(-gdet)*ucon[2];

    X1 = i_TO_X1_CENTER(i); X2 = j_TO_X2_CENTER(j-1);
    for (int var=0; var<DOF; var++) {
        primVars[var] = primTile[INDEX_LOCAL(iTile, jTile-1, var)];
    }
    gCovCalc(gcov, X1, X2);
    gDetCalc(&gdet, gcov);
    gConCalc(gcon, gcov, gdet);
    gammaCalc(&gamma, primVars, gcov);
    alphaCalc(&alpha, gcon);
    uconCalc(ucon, gamma, alpha, primVars, gcon);

    left = sqrt(-gdet)*ucon[2];

    ducon = SlopeLim(left, center, right);

    dU_dt[FF] += -primTile[INDEX_LOCAL(iTile, jTile, FF)]*ducon/DX2;

    dU_dt[FF] += -primTile[INDEX_LOCAL(iTile, jTile, FF)]*(g*ducon_dt[0]);

    dT_dt = (ADIABATIC_INDEX-1.)*
            ((dvars_dt[UU]/primTile[INDEX_LOCAL(iTile, jTile, RHO)]) -
              (primTile[INDEX_LOCAL(iTile, jTile, UU)]*dvars_dt[RHO]
              /pow(primTile[INDEX_LOCAL(iTile, jTile, RHO)], 2.))
            );

    dT_dX1 = (ADIABATIC_INDEX-1.)*
             SlopeLim(primTile[INDEX_LOCAL(iTile-1, jTile, UU)]/
                      primTile[INDEX_LOCAL(iTile-1, jTile, RHO)],
                      primTile[INDEX_LOCAL(iTile, jTile, UU)]/
                      primTile[INDEX_LOCAL(iTile, jTile, RHO)],
                      primTile[INDEX_LOCAL(iTile+1, jTile, UU)]/
                      primTile[INDEX_LOCAL(iTile+1, jTile, RHO)] )/DX1;

    dT_dX2 = (ADIABATIC_INDEX-1.)*
             SlopeLim(primTile[INDEX_LOCAL(iTile, jTile-1, UU)]/
                      primTile[INDEX_LOCAL(iTile, jTile-1, RHO)],
                      primTile[INDEX_LOCAL(iTile, jTile, UU)]/
                      primTile[INDEX_LOCAL(iTile, jTile, RHO)],
                      primTile[INDEX_LOCAL(iTile, jTile+1, UU)]/
                      primTile[INDEX_LOCAL(iTile, jTile+1, RHO)] )/DX2;

    REAL acon[NDIM], acov[NDIM];

    /* Compute a^\mu = u^\mu\del_\nu u^\mu 
                     = u^0 d(u^\mu)/dt + u^1 d(u^\mu)/dX1 + u^2 d(u^\mu)/dX2
     */

    for (int mu=0; mu<NDIM; mu++)
    {
        acon[mu] = uconCenter[0]*ducon_dt[mu];

        center = uconCenter[mu];

        X1 = i_TO_X1_CENTER(i+1); X2 = j_TO_X2_CENTER(j);
        for (int var=0; var<DOF; var++)
        {
          primVars[var] = primTile[INDEX_LOCAL(iTile+1, jTile, var)];
        }
        gCovCalc(gcov, X1, X2);
        gDetCalc(&gdet, gcov);
        gConCalc(gcon, gcov, gdet);
        gammaCalc(&gamma, primVars, gcov);
        alphaCalc(&alpha, gcon);
        uconCalc(ucon, gamma, alpha, primVars, gcon);

        right = ucon[mu];

        X1 = i_TO_X1_CENTER(i-1); X2 = j_TO_X2_CENTER(j);
        for (int var=0; var<DOF; var++)
        {
          primVars[var] = primTile[INDEX_LOCAL(iTile-1, jTile, var)];
        }
        gCovCalc(gcov, X1, X2);
        gDetCalc(&gdet, gcov);
        gConCalc(gcon, gcov, gdet);
        gammaCalc(&gamma, primVars, gcov);
        alphaCalc(&alpha, gcon);
        uconCalc(ucon, gamma, alpha, primVars, gcon);

        left = ucon[mu];

        ducon = SlopeLim(left, center, right);

        acon[mu] += uconCenter[1]*ducon/DX1;

        X1 = i_TO_X1_CENTER(i); X2 = j_TO_X2_CENTER(j+1);
        for (int var=0; var<DOF; var++)
        {
          primVars[var] = primTile[INDEX_LOCAL(iTile, jTile+1, var)];
        }
        gCovCalc(gcov, X1, X2);
        gDetCalc(&gdet, gcov);
        gConCalc(gcon, gcov, gdet);
        gammaCalc(&gamma, primVars, gcov);
        alphaCalc(&alpha, gcon);
        uconCalc(ucon, gamma, alpha, primVars, gcon);

        right = ucon[mu];

        X1 = i_TO_X1_CENTER(i); X2 = j_TO_X2_CENTER(j-1);
        for (int var=0; var<DOF; var++)
        {
          primVars[var] = primTile[INDEX_LOCAL(iTile, jTile-1, var)];
        }
        gCovCalc(gcov, X1, X2);
        gDetCalc(&gdet, gcov);
        gConCalc(gcon, gcov, gdet);
        gammaCalc(&gamma, primVars, gcov);
        alphaCalc(&alpha, gcon);
        uconCalc(ucon, gamma, alpha, primVars, gcon);

        left = ucon[mu];

        ducon = SlopeLim(left, center, right);

        acon[mu] += uconCenter[2]*ducon/DX2;

        for (int alpha=0; alpha<NDIM; alpha++)
        {
          for (int beta=0; beta<NDIM; beta++)
          {
            acon[mu] += conn[mu][alpha][beta]
                        *uconCenter[alpha]*uconCenter[beta];
          }
        }
    }

    for (int mu=0; mu<NDIM; mu++)
    {
      acov[mu] = 0.;
      for (int nu=0; nu<NDIM; nu++)
      {
        acov[mu] += gcov[mu][nu]*acon[nu];
      }
    }

    REAL qEckartCon[NDIM], kappa;
    REAL dT[NDIM];

    /* Set value for kappa using pitch angle scattering estimates from Sharma,
     * Quataert and Stone, 2008.
     * kappa = .2 (GMr)^0.5
     */
    REAL r, theta;
    X1 = i_TO_X1_CENTER(i); X2 = j_TO_X2_CENTER(j);
    BLCoords(&r, &theta, X1, X2);
    kappa = 1.*sqrt(r)*primTile[INDEX_LOCAL(iTile, jTile, RHO)];

    dT[0] = dT_dt; dT[1] = dT_dX1; dT[2] = dT_dX2; dT[3] = 0.;

    REAL T = (primTile[INDEX_LOCAL(iTile, jTile, UU)]/
              primTile[INDEX_LOCAL(iTile, jTile, RHO)])*
             (ADIABATIC_INDEX-1.);

    for (int mu=0; mu<NDIM; mu++)
    {
      qEckartCon[mu] = 0.;
      for (int nu=0; nu<NDIM; nu++)
      {
        qEckartCon[mu] += -kappa*(uconCenter[mu]*uconCenter[nu] + gcon[mu][nu])
                                *(dT[nu] + T*acov[nu]);
      }                   
    }

    REAL bDotq = bcov[0]*qEckartCon[0] + bcov[1]*qEckartCon[1] +
                 bcov[2]*qEckartCon[2] + bcov[3]*qEckartCon[3];

    REAL bsqr;
    bSqrCalc(&bsqr, bcon, bcov);

    cs = sqrt(ADIABATIC_INDEX*(ADIABATIC_INDEX - 1)*\
              primTile[INDEX_LOCAL(iTile, jTile, UU)]/
              (primTile[INDEX_LOCAL(iTile, jTile, RHO)] +
               primTile[INDEX_LOCAL(iTile, jTile, UU)]) );
    
    REAL qsat = PHI*(primTile[INDEX_LOCAL(iTile, jTile, RHO)] +
                     primTile[INDEX_LOCAL(iTile, jTile, UU)])*pow(cs,
                     3.);
    
//    F0 = qsat/sqrt(bsqr)*tanh(bDotq/sqrt(bsqr)/qsat);
    F0 = bDotq/sqrt(bsqr);

//    REAL tau_r = TAU_R_SAFETY_FACTOR*(kappa*T)/
//                 (primTile[INDEX_LOCAL(iTile, jTile, RHO)] +
//                  ADIABATIC_INDEX*primTile[INDEX_LOCAL(iTile, jTile, UU)]
//                 ) + DT;

    REAL density = primTile[INDEX_LOCAL(iTile, jTile, RHO)];
    REAL pressure = 
    (ADIABATIC_INDEX-1.)*primTile[INDEX_LOCAL(iTile, jTile, UU)];
    REAL temperature = pressure/density;
    REAL thermal_velocity = sqrt(temperature);
    REAL tau_r = TAU_R_SAFETY_FACTOR*kappa/density/pow(thermal_velocity, 2.);

    if (tau_r < DT) tau_r = DT;
    
    dU_dt[FF] += g*(primTile[INDEX_LOCAL(iTile, jTile, FF)] - F0)/tau_r;

#endif /* Conduction */

}
#endif /*OPENCL code*/

void ComputeFluxAndU(REAL flux[DOF],
                     REAL U[DOF],
                     const REAL ucon[NDIM],
                     const REAL ucov[NDIM],
                     const REAL bcon[NDIM],
                     const REAL bcov[NDIM],
                     const REAL gcon[NDIM][NDIM],
                     const REAL gcov[NDIM][NDIM],
                     const REAL mhd[NDIM][NDIM],
                     const REAL var[DOF],
                     const REAL g,
                     const int dir)
{
    flux[RHO] = g*var[RHO]*ucon[dir];

    flux[UU] = g*mhd[dir][0];
    flux[U1] = g*mhd[dir][1];
    flux[U2] = g*mhd[dir][2];
    flux[U3] = g*mhd[dir][3];

    flux[B1] = g*(bcon[1]*ucon[dir] - bcon[dir]*ucon[1]);
    flux[B2] = g*(bcon[2]*ucon[dir] - bcon[dir]*ucon[2]);
    flux[B3] = g*(bcon[3]*ucon[dir] - bcon[dir]*ucon[3]);

    U[RHO] = g*var[RHO]*ucon[0];

    U[UU] = g*mhd[0][0];
    U[U1] = g*mhd[0][1];
    U[U2] = g*mhd[0][2];
    U[U3] = g*mhd[0][3];

    U[B1] = g*(bcon[1]*ucon[0] - bcon[0]*ucon[1]);
    U[B2] = g*(bcon[2]*ucon[0] - bcon[0]*ucon[2]);
    U[B3] = g*(bcon[3]*ucon[0] - bcon[0]*ucon[3]);
    
#if(CONDUCTION)
    flux[FF] = g*(ucon[dir]*var[FF]);
    U[FF] = g*(ucon[0]*var[FF]);
#endif /* Conduction */
}

void ComputedU_dt(REAL dU_dt[DOF],
                  const REAL ucon[NDIM],
                  const REAL ducon_dt[NDIM],
                  const REAL ucov[NDIM],
                  const REAL ducov_dt[NDIM],
                  const REAL bcon[NDIM],
                  const REAL dbcon_dt[NDIM],
                  const REAL bcov[NDIM],
                  const REAL dbcov_dt[NDIM],
                  const REAL gcon[NDIM][NDIM],
                  const REAL gcov[NDIM][NDIM],
                  const REAL var[DOF],
                  const REAL dvar_dt[DOF],
                  const REAL gamma,
                  const REAL dgamma_dt,
                  const REAL alpha,
                  const REAL g)
{
    REAL P, dP_dt, bsqr, dbsqr_dt, tmp1, dtmp1_dt, tmp2, dtmp2_dt;

    bSqrCalc(&bsqr, bcon, bcov);

    dbsqr_dt = bcon[0]*dbcov_dt[0] + dbcon_dt[0]*bcov[0] +
               bcon[1]*dbcov_dt[1] + dbcon_dt[1]*bcov[1] +
               bcon[2]*dbcov_dt[2] + dbcon_dt[2]*bcov[2] +
               bcon[3]*dbcov_dt[3] + dbcon_dt[3]*bcov[3];

    REAL qcon[NDIM], qcov[NDIM];
    REAL dqcon_dt[NDIM], dqcov_dt[NDIM];

    for (int mu=0; mu<NDIM; mu++)
    {
#if (CONDUCTION)
      qcon[mu] = var[FF]*bcon[mu]/sqrt(bsqr);

      dqcon_dt[mu] =    dvar_dt[FF]*bcon[mu]/sqrt(bsqr)
                      + var[FF]/sqrt(bsqr)*dbcon_dt[mu]
                      - 0.5*var[FF]*bcon[mu]/pow(bsqr, 3./2)*dbsqr_dt;
#else
      qcon[mu] = 0.;
      dqcon_dt[mu] = 0.;
#endif
    }

    for (int mu=0; mu<NDIM; mu++)
    {
      qcov[mu] = 0.;
      dqcov_dt[mu] = 0.;
      for (int nu=0; nu<NDIM; nu++)
      {
        qcov[mu] += gcov[mu][nu]*qcon[nu];

        dqcov_dt[mu] += gcov[mu][nu]*dqcon_dt[nu];
      }
    }

    P = (ADIABATIC_INDEX-1.)*var[UU];
    dP_dt = (ADIABATIC_INDEX-1.)*dvar_dt[UU];

    tmp1 = P + var[RHO] + var[UU] + bsqr;
    dtmp1_dt = dP_dt + dvar_dt[RHO] + dvar_dt[UU] + dbsqr_dt;

    tmp2 = P + 0.5*bsqr;
    dtmp2_dt = dP_dt + 0.5*dbsqr_dt;

    dU_dt[RHO] = g*(dvar_dt[RHO]*ucon[0] + var[RHO]*ducon_dt[0]);

    dU_dt[UU] = g*(dtmp1_dt*ucon[0]*ucov[0] +
                   tmp1*(ducon_dt[0]*ucov[0] + ucon[0]*ducov_dt[0]) +
                   dtmp2_dt - dbcon_dt[0]*bcov[0] - bcon[0]*dbcov_dt[0]
                   + dqcon_dt[0]*ucov[0] + qcon[0]*ducov_dt[0]
                   + ducon_dt[0]*qcov[0] + ucon[0]*dqcov_dt[0]);

    dU_dt[U1] = g*(dtmp1_dt*ucon[0]*ucov[1] +
                   tmp1*(ducon_dt[0]*ucov[1] + ucon[0]*ducov_dt[1]) -
                   dbcon_dt[0]*bcov[1] - bcon[0]*dbcov_dt[1]
                   + dqcon_dt[0]*ucov[1] + qcon[0]*ducov_dt[1]
                   + ducon_dt[0]*qcov[1] + ucon[0]*dqcov_dt[1]);

    dU_dt[U2] = g*(dtmp1_dt*ucon[0]*ucov[2] +
                   tmp1*(ducon_dt[0]*ucov[2] + ucon[0]*ducov_dt[2]) -
                   dbcon_dt[0]*bcov[2] - bcon[0]*dbcov_dt[2]
                   + dqcon_dt[0]*ucov[2] + qcon[0]*ducov_dt[2]
                   + ducon_dt[0]*qcov[2] + ucon[0]*dqcov_dt[2]);

    dU_dt[U3] = g*(dtmp1_dt*ucon[0]*ucov[3] +
                   tmp1*(ducon_dt[0]*ucov[3] + ucon[0]*ducov_dt[3]) -
                   dbcon_dt[0]*bcov[3] - bcon[0]*dbcov_dt[3] 
                   + dqcon_dt[0]*ucov[3] + qcon[0]*ducov_dt[3]
                   + ducon_dt[0]*qcov[3] + ucon[0]*dqcov_dt[3]);
                   

    dU_dt[B1] = g*dvar_dt[B1];
    dU_dt[B2] = g*dvar_dt[B2];
    dU_dt[B3] = g*dvar_dt[B3];

#if(CONDUCTION)

    dU_dt[FF] = g*(dvar_dt[FF]*ucon[0] + ducon_dt[0]*var[FF]);

#endif /* Conduction */


}
                  
void VChar(REAL* vmin, REAL* vmax,
           const REAL ucon[DOF], const REAL ucov[DOF],
           const REAL bsqr,
           const REAL gcon[NDIM][NDIM],
           const REAL var[DOF], const int dir)
{
    REAL Acov[NDIM], Acon[NDIM], Bcov[NDIM], Bcon[NDIM];
    REAL Asqr, Bsqr, vasqr, cssqr, cmsqr, Adotu, Bdotu, AdotB;
    REAL A, B, C, discr;
    REAL vp, vm;

    vasqr = bsqr/(bsqr + var[RHO] + ADIABATIC_INDEX*var[UU]);
    cssqr = (ADIABATIC_INDEX)*(ADIABATIC_INDEX-1)*var[UU]/(var[RHO] +
             ADIABATIC_INDEX*var[UU]);

    cmsqr = cssqr + vasqr - cssqr*vasqr;
    
    for (int mu=0; mu<NDIM; mu++)
        Acov[mu] = 0.;
    Acov[dir] = 1.;
    conFromCov(Acon, Acov, gcon);

    for (int mu=0; mu<NDIM; mu++)
        Bcov[mu] = 0.;
    Bcov[0] = 1.;
    conFromCov(Bcon, Bcov, gcon);

    conDotCov(&Asqr, Acon, Acov);
    conDotCov(&Bsqr, Bcon, Bcov);
    conDotCov(&Adotu, ucon, Acov);
    conDotCov(&Bdotu, ucon, Bcov);
    conDotCov(&AdotB, Acon, Bcov);

    A = (Bdotu*Bdotu) - (Bsqr + Bdotu*Bdotu)*cmsqr;
    B = 2.*(Adotu*Bdotu - (AdotB + Adotu*Bdotu)*cmsqr);
    C = Adotu*Adotu - (Asqr + Adotu*Adotu)*cmsqr;

    discr = sqrt(B*B - 4.*A*C);

    vp = -(-B + discr)/(2.*A);
    vm = -(-B - discr)/(2.*A);

    if (vp>vm) {
        *vmax = vp;
        *vmin = vm;
    } else {
        *vmax = vm;
        *vmin = vp;
    }
}

void setFloor(REAL var[DOF], 
              const REAL gamma,
              REAL X1, REAL X2)
{
    REAL r, theta;
    BLCoords(&r, &theta, X1, X2);

    REAL rhoFloor = RHO_MIN_LIMIT;

    REAL uFloor = U_MIN_LIMIT;
    
    if (var[RHO] < rhoFloor)
        var[RHO] = rhoFloor;

    if (var[UU] < uFloor)
        var[UU] = uFloor;

//    if (gamma > GAMMA_MAX) {
//        REAL f = sqrt((GAMMA_MAX*GAMMA_MAX - 1.) /
//				      (gamma * gamma - 1.));
//
//        var[U1] = f*var[U1];
//        var[U2] = f*var[U2];
//        var[U3] = f*var[U3];
//
//    }
}

void setFloorInit(REAL var[DOF], 
              const REAL gamma,
              REAL X1, REAL X2)
{
    REAL r, theta;
    BLCoords(&r, &theta, X1, X2);

    REAL rhoFloor = RHO_MIN*pow(r, -1.5);
    REAL uFloor = U_MIN*pow(r, -2.5);

    if (rhoFloor < RHO_MIN_LIMIT)
        rhoFloor = RHO_MIN_LIMIT;

    if (uFloor < U_MIN_LIMIT)
        uFloor = U_MIN_LIMIT;
    
    if (var[RHO] < rhoFloor)
        var[RHO] = rhoFloor;

    if (var[UU] < uFloor)
        var[UU] = uFloor;

    if (gamma > GAMMA_MAX) {
        REAL f = sqrt((GAMMA_MAX*GAMMA_MAX - 1.) /
				      (gamma * gamma - 1.));

        var[U1] = f*var[U1];
        var[U2] = f*var[U2];
        var[U3] = f*var[U3];

    }
}
