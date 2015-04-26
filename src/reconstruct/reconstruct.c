#include "reconstruct.h"

/* Taken from HARM */
void slopeLim(const REAL left[ARRAY_ARGS NUM_FLUXES],
              const REAL mid[ARRAY_ARGS NUM_FLUXES],
              const REAL right[ARRAY_ARGS NUM_FLUXES],
              REAL ans[ARRAY_ARGS NUM_FLUXES])
{
  for (int var=0; var<NUM_FLUXES; var++)
  {
    /* Monotonized Central slope limiter */
#if (RECONSTRUCTION == MONOTONIZED_CENTRAL)
    REAL Dqm = 2. * (mid[var] - left[var]);
	  REAL Dqp = 2. * (right[var] - mid[var]);
	  REAL Dqc = 0.5 * (right[var] - left[var]);
	  REAL s = Dqm * Dqp;
	  if (s <= 0.) 
    {
		  ans[var] = 0.;
    }
	  else
    {
		  if (fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
      {
			  ans[var] = Dqm;
      }
		  else if (fabs(Dqp) < fabs(Dqc))
      {
			  ans[var] = Dqp;
      }
		  else
      {
			  ans[var] = Dqc;
      }
	  }
#elif (RECONSTRUCTION == MIN_MOD)
		REAL Dqm = (mid[var] - left[var]);
		REAL Dqp = (right[var] - mid[var]);
		REAL s = Dqm*Dqp;
		if(s <= 0.) ans[var] = 0.;
		else if(fabs(Dqm) < fabs(Dqp)) ans[var] = Dqm;
		else ans[var] = Dqp;
#endif
  }

}

REAL slopeLimitedDerivative(const REAL left, const REAL mid, const REAL right)
{
  REAL Dqm = 2. * (mid - left);
  REAL Dqp = 2. * (right - mid);
  REAL Dqc = .5 * (right - left);
  REAL s = Dqm * Dqp;
  if (s <= 0.) 
  {
	  return 0.;
  }
  else
  {
	  if (fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
    {
		  return Dqm;
    }
	  else if (fabs(Dqp) < fabs(Dqc))
    {
		  return Dqp;
    }
	  else
    {
		  return Dqc;
    }
  }
}

REAL MP5_Reconstruct(const REAL Fjm2,
                     const REAL Fjm1,
                     const REAL Fj,
                     const REAL Fjp1,
                     const REAL Fjp2)
{
  double f, d2, d2p, d2m; 
  double dMMm, dMMp;
  double scrh1,scrh2, fmin, fmax; 
  double fAV, fMD, fLC, fUL, fMP;
  static double alpha = 4.0, epsm = 1.e-12;

  f  = 2.0*Fjm2 - 13.0*Fjm1 + 47.0*Fj + 27.0*Fjp1 - 3.0*Fjp2;
  f /= 60.0;   

  fMP = Fj + MINMOD(Fjp1-Fj,alpha*(Fj-Fjm1));

  if ((f - Fj)*(f - fMP) <= epsm) return f;

  d2m = Fjm2 + Fj - 2.0*Fjm1;    /* -- Eq. (2.19) -- */
  d2  = Fjm1 + Fjp1 - 2.0*Fj;
  d2p = Fj + Fjp2 - 2.0*Fjp1;    /* -- Eq. (2.19) -- */

  scrh1 = MINMOD(4.0*d2 - d2p, 4.0*d2p - d2);
  scrh2 = MINMOD(d2, d2p);
  dMMp  = MINMOD(scrh1,scrh2);   /* -- Eq. (2.27) -- */

  scrh1 = MINMOD(4.0*d2m - d2, 4.0*d2 - d2m);
  scrh2 = MINMOD(d2, d2m);
  dMMm  = MINMOD(scrh1,scrh2);   /* -- Eq. (2.27) -- */

  fUL = Fj + alpha*(Fj - Fjm1);   /* -- Eq. (2.8) -- */
  fAV = 0.5*(Fj + Fjp1);        /* -- Eq. (2.16) -- */
  fMD = fAV - 0.5*dMMp; /* -- Eq. (2.28) -- */
  fLC = 0.5*(3.0*Fj - Fjm1) + 4.0/3.0*dMMm;  /* -- Eq. (2.29) -- */

  scrh1 = MY_MIN(Fj, Fjp1); scrh1 = MY_MIN(scrh1, fMD);
  scrh2 = MY_MIN(Fj, fUL);    scrh2 = MY_MIN(scrh2, fLC);
  fmin  = MY_MAX(scrh1, scrh2);  /* -- Eq. (2.24a) -- */

  scrh1 = MY_MAX(Fj, Fjp1); scrh1 = MY_MAX(scrh1, fMD);
  scrh2 = MY_MAX(Fj, fUL);    scrh2 = MY_MAX(scrh2, fLC);
  fmax  = MY_MIN(scrh1, scrh2);  /* -- Eq. (2.24b) -- */

  f = MEDIAN(f, fmin, fmax); /* -- Eq. (2.26) -- */
  return f;
}

void reconstruct(const REAL primTile[ARRAY_ARGS TILE_SIZE],
                 const int dir,
                 const int iTile, const int jTile,
                 const int X1Start, const int X2Start,
                 const int X1Size, const int X2Size,
                 REAL primVarsLeft[ARRAY_ARGS TILE_SIZE],
                 REAL primVarsRight[ARRAY_ARGS TILE_SIZE])
{
  if (dir==X1)
  {
    LOOP_INSIDE_TILE(-2, TILE_SIZE_X1+2, -2, TILE_SIZE_X2+2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  X1Start, X2Start, 
                  X1Size, X2Size, 
                  &zone);

      #if (RECONSTRUCTION == MONOTONIZED_CENTRAL || \
           RECONSTRUCTION == MIN_MOD)
        REAL slope[NUM_FLUXES];
  
        slopeLim(&primTile[INDEX_TILE_MINUS_ONE_X1(&zone, 0)],
                 &primTile[INDEX_TILE(&zone, 0)],
                 &primTile[INDEX_TILE_PLUS_ONE_X1(&zone, 0)],
                 slope);

        for (int var=0; var<NUM_FLUXES; var++)
        {
          /* Left Edge */
          primVarsLeft[INDEX_TILE(&zone, var)] =
            primTile[INDEX_TILE(&zone, var)] - 0.5*slope[var];

          /* Right Edge */
          primVarsRight[INDEX_TILE(&zone, var)] =
            primTile[INDEX_TILE(&zone, var)] + 0.5*slope[var];
        }
      #elif (RECONSTRUCTION == MP5)
        for (int var=0; var<NUM_FLUXES; var++)
        {
          primVarsRight[INDEX_TILE(&zone, var)] = 
            MP5_Reconstruct(
                primTile[INDEX_TILE_MINUS_TWO_X1(&zone, var)],
                primTile[INDEX_TILE_MINUS_ONE_X1(&zone, var)],
                primTile[INDEX_TILE(&zone, var)],
                primTile[INDEX_TILE_PLUS_ONE_X1(&zone, var)],
                primTile[INDEX_TILE_PLUS_TWO_X1(&zone, var)]);

          primVarsLeft[INDEX_TILE(&zone, var)] = 
            MP5_Reconstruct(
                primTile[INDEX_TILE_PLUS_TWO_X1(&zone, var)],
                primTile[INDEX_TILE_PLUS_ONE_X1(&zone, var)],
                primTile[INDEX_TILE(&zone, var)],
                primTile[INDEX_TILE_MINUS_ONE_X1(&zone, var)],
                primTile[INDEX_TILE_MINUS_TWO_X1(&zone, var)]);
        }
      #endif

    }

  }
#if (COMPUTE_DIM==2)
  else if (dir==X2)
  {
    LOOP_INSIDE_TILE(-2, TILE_SIZE_X1+2, -2, TILE_SIZE_X2+2)
    {
      struct gridZone zone;
      setGridZone(iTile, jTile,
                  iInTile, jInTile,
                  X1Start, X2Start, 
                  X1Size, X2Size, 
                  &zone);
    
      #if (RECONSTRUCTION == MONOTONIZED_CENTRAL || \
           RECONSTRUCTION == MIN_MOD)
        REAL slope[NUM_FLUXES];
  
        slopeLim(&primTile[INDEX_TILE_MINUS_ONE_X2(&zone, 0)],
                 &primTile[INDEX_TILE(&zone, 0)],
                 &primTile[INDEX_TILE_PLUS_ONE_X2(&zone, 0)],
                 slope);

        for (int var=0; var<NUM_FLUXES; var++)
        {
          /* Left Edge */
          primVarsLeft[INDEX_TILE(&zone, var)] =
            primTile[INDEX_TILE(&zone, var)] - 0.5*slope[var];
  
          /* Right Edge */
          primVarsRight[INDEX_TILE(&zone, var)] =
            primTile[INDEX_TILE(&zone, var)] + 0.5*slope[var];
        } 

      #elif (RECONSTRUCTION == MP5)
        for (int var=0; var<NUM_FLUXES; var++)
        {
          primVarsRight[INDEX_TILE(&zone, var)] = 
            MP5_Reconstruct(
                primTile[INDEX_TILE_MINUS_TWO_X2(&zone, var)],
                primTile[INDEX_TILE_MINUS_ONE_X2(&zone, var)],
                primTile[INDEX_TILE(&zone, var)],
                primTile[INDEX_TILE_PLUS_ONE_X2(&zone, var)],
                primTile[INDEX_TILE_PLUS_TWO_X2(&zone, var)]);

          primVarsLeft[INDEX_TILE(&zone, var)] = 
            MP5_Reconstruct(
                primTile[INDEX_TILE_PLUS_TWO_X2(&zone, var)],
                primTile[INDEX_TILE_PLUS_ONE_X2(&zone, var)],
                primTile[INDEX_TILE(&zone, var)],
                primTile[INDEX_TILE_MINUS_ONE_X2(&zone, var)],
                primTile[INDEX_TILE_MINUS_TWO_X2(&zone, var)]);

        }
      #endif

    }

  }
  #endif

}
