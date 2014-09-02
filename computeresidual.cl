#include "constants.h"
#include "reconstruct.cl"
#include "phys.cl"
#include "riemannsolver.cl"

__kernel void ComputeResidual(__global const REAL* restrict prim, 
                              __global const REAL* restrict dprim_dt,
                              __global REAL* restrict F,
                              __global REAL* restrict fluxX1,
                              __global REAL* restrict fluxX2,
                              __global const REAL* restrict primBoundaries)
{
    int i = get_global_id(0);
    int j = get_global_id(1);

    int iTile = get_local_id(0);
    int jTile = get_local_id(1);

    // Finite Volume variables
    REAL primEdge[DOF];
    REAL fluxL[DOF], uL[DOF], fluxR[DOF], uR[DOF];
    REAL fluxX1L[DOF], fluxX1R[DOF], fluxX2L[DOF], fluxX2R[DOF];
    REAL cminL, cminR, cmaxL, cmaxR;

    // Geometry variables
    REAL X1, X2;
    REAL gcon[NDIM][NDIM], gcov[NDIM][NDIM], gdet;
    REAL alpha, gamma, g;

    // Physics variables
    REAL dU_dt[DOF], dvars_dt[DOF];
    REAL ucon[NDIM], ucov[NDIM], bcon[NDIM], bcov[NDIM];
    REAL mhd[NDIM][NDIM], bsqr;

    // Time derivative of physics variables
    REAL ducon_dt[NDIM], ducov_dt[NDIM], dbcon_dt[NDIM], dbcov_dt[NDIM];
    REAL dgamma_dt;

    __local REAL primTile[(TILE_SIZE_X1+2*NG)*(TILE_SIZE_X2+2*NG)*DOF];

    if (i>=0 && i<N1 && j>=0 && j<N2) {
        for (int var=0; var<DOF; var++) {
            primTile[INDEX_LOCAL(iTile,jTile,var)] =
            prim[INDEX_GLOBAL(i,j,var)];
        }
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    if (i>=0 && i<N1 && j>=0 && j<N2) {
        for (int var=0; var<DOF; var++) {
            primEdge[var] = primTile[INDEX_LOCAL(iTile,jTile,var)];
            dvars_dt[var] = dprim_dt[INDEX_GLOBAL(i,j,var)];
        }

        X1 = i_TO_X1_CENTER(i); X2 = j_TO_X2_CENTER(j);
        gCovCalc(gcov, X1, X2);
        gDetCalc(&gdet, gcov);
        gConCalc(gcon, gcov, gdet);
        g = sqrt(-gdet);
    
        gammaCalc(&gamma, primEdge, gcov);
        setFloor(primEdge, gamma, X1, X2);
        gammaCalc(&gamma, primEdge, gcov);
        alphaCalc(&alpha, gcon);
        uconCalc(ucon, gamma, alpha, primEdge, gcon);
        covFromCon(ucov, ucon, gcov);
        bconCalc(bcon, primEdge, ucon, ucov);
        covFromCon(bcov, bcon, gcov);
    
        mhdCalc(mhd, primEdge, ucon, ucov, bcon, bcov, dvars_dt, primTile);
    
        dgammaCalc_dt(&dgamma_dt, gamma, primEdge, dvars_dt, gcov);
        duconCalc_dt(ducon_dt, dgamma_dt, alpha, dvars_dt, gcon);
        covFromCon(ducov_dt, ducon_dt, gcov);
    
        dbconCalc_dt(dbcon_dt, ucon, ducon_dt, ucov, ducov_dt, bcon, 
                     primEdge, dvars_dt);
        covFromCon(dbcov_dt, dbcon_dt, gcov);
    
        ComputedU_dt(dU_dt,
                     ucon, ducon_dt, ucov, ducov_dt,
                     bcon, dbcon_dt, bcov, dbcov_dt,
                     gcon, gcov, primEdge, dvars_dt,
                     gamma, dgamma_dt, alpha, g);
    }

    /* Set boundaries */
    for (int var=0; var<DOF; var++) {

        if (iTile==0) {
            if (i>=NG) {
                for (int iNg=-NG; iNg<0; iNg++) {
                    primTile[INDEX_LOCAL(iNg,jTile,var)] =
                        prim[INDEX_GLOBAL(i+iNg,j,var)];
                }
            } else {
                for (int iNg=-NG; iNg<0; iNg++) {
                  #if (LEFT_BOUNDARY==OUTFLOW)
                    primTile[INDEX_LOCAL(iNg,jTile,var)] =
                    primTile[INDEX_LOCAL(0,jTile,var)];
                  #elif (LEFT_BOUNDARY==DIRICHLET)
                    primTile[INDEX_LOCAL(iNg,jTile,var)] = 
                    primBoundaries[INDEX_GLOBAL_WITH_NG(i+iNg,j,var)];
                  #endif
                }
            }
        }
    
        if (iTile==TILE_SIZE_X1-1) {
            if (i<=N1-NG) {
                for (int iNg=0; iNg<NG; iNg++) {
                    primTile[INDEX_LOCAL(TILE_SIZE_X1+iNg,jTile,var)] =
                        prim[INDEX_GLOBAL(i+iNg+1,j,var)];
                }
            } else {
                for (int iNg=0; iNg<NG; iNg++) {
                  #if (RIGHT_BOUNDARY==OUTFLOW)
                    primTile[INDEX_LOCAL(TILE_SIZE_X1+iNg,jTile,var)] =
                    primTile[INDEX_LOCAL(TILE_SIZE_X1-1,jTile,var)];
                  #elif (RIGHT_BOUNDARY==DIRICHLET)
                    primTile[INDEX_LOCAL(TILE_SIZE_X1+iNg,jTile,var)] =
                    primBoundaries[INDEX_GLOBAL_WITH_NG(i+iNg+1,j,var)];
                  #endif
                }
            }
        }
       
        if (jTile==0) {
            if (j>=NG) {
                for (int jNg=-NG; jNg<0; jNg++) {
                    primTile[INDEX_LOCAL(iTile,jNg,var)] =
                        prim[INDEX_GLOBAL(i,j+jNg,var)];
                }
            } else {
                for (int jNg=-NG; jNg<0; jNg++) {
                  #if (BOTTOM_BOUNDARY==MIRROR)
                    primTile[INDEX_LOCAL(iTile,jNg,var)] =
                    primTile[INDEX_LOCAL(iTile,-jNg-1,var)];
                  #elif (BOTTOM_BOUNDARY==DIRICHLET)
                    primTile[INDEX_LOCAL(iTile,jNg,var)] =
                    primBoundaries[INDEX_GLOBAL_WITH_NG(i,j+jNg,var)];
                  #endif
                }
            }
        }
    
        if (jTile==TILE_SIZE_X2-1) {
            if (j<=N2-NG) {
                for (int jNg=0; jNg<NG; jNg++) {
                    primTile[INDEX_LOCAL(iTile,TILE_SIZE_X2+jNg,var)] =
                        prim[INDEX_GLOBAL(i,j+jNg+1,var)];
                }
            } else {
                for (int jNg=0; jNg<NG; jNg++) {
                  #if (TOP_BOUNDARY==MIRROR)
                    primTile[INDEX_LOCAL(iTile,TILE_SIZE_X2+jNg,var)] =
                    primTile[INDEX_LOCAL(iTile,TILE_SIZE_X2-1-jNg,var)];
                  #elif (TOP_BOUNDARY==DIRICHLET)
                    primTile[INDEX_LOCAL(iTile,TILE_SIZE_X2+jNg,var)] =
                    primBoundaries[INDEX_GLOBAL_WITH_NG(i,j+jNg+1,var)];
                  #endif
                }
            }
        }
    
    }

#if (BOTTOM_BOUNDARY==MIRROR)
    if (jTile==0 && j<NG) {
        for (int jNg=-NG; jNg<0; jNg++) {
            primTile[INDEX_LOCAL(iTile,jNg,U2)] =
           -primTile[INDEX_LOCAL(iTile,jNg,U2)];
            primTile[INDEX_LOCAL(iTile,jNg,B2)] =
           -primTile[INDEX_LOCAL(iTile,jNg,B2)];
        }
    }
#endif

#if (TOP_BOUNDARY==MIRROR)
    if (jTile==TILE_SIZE_X2-1 && j>N2-NG) {
        for (int jNg=0; jNg<NG; jNg++) {
            primTile[INDEX_LOCAL(iTile,TILE_SIZE_X2+jNg,U2)] =
           -primTile[INDEX_LOCAL(iTile,TILE_SIZE_X2+jNg,U2)];
            primTile[INDEX_LOCAL(iTile,TILE_SIZE_X2+jNg,B2)] =
           -primTile[INDEX_LOCAL(iTile,TILE_SIZE_X2+jNg,B2)];
        }
    }
#endif

    barrier(CLK_LOCAL_MEM_FENCE);

#if (INFLOW_CHECK)
    if (i>=N1-1 && iTile==TILE_SIZE_X1-1) {
        for (int iNg=0; iNg<NG; iNg++) {
            for (int var=0; var<DOF; var++) {
                primEdge[var] = 
                primTile[INDEX_LOCAL(TILE_SIZE_X1+iNg,jTile,var)];
            }

            X1 = i_TO_X1_CENTER(i+1+iNg); X2 = j_TO_X2_CENTER(j);
            
            gCovCalc(gcov, X1, X2);
            gDetCalc(&gdet, gcov);
            gConCalc(gcon, gcov, gdet);
            g = sqrt(-gdet);

            gammaCalc(&gamma, primEdge, gcov);
            setFloor(primEdge, gamma, X1, X2);
            gammaCalc(&gamma, primEdge, gcov);
            alphaCalc(&alpha, gcon);
            uconCalc(ucon, gamma, alpha, primEdge, gcon);

            if (ucon[1] < 0.) {

                primEdge[U1] /= gamma;
                primEdge[U2] /= gamma;
                primEdge[U3] /= gamma;

                alpha = 1./sqrt(-gcon[0][0]);
                
                primEdge[U1] = gcon[0][1]*alpha;

                double vSqr = gcov[1][1]*primEdge[U1]*primEdge[U1] +
                              gcov[1][2]*primEdge[U1]*primEdge[U2] + 
                              gcov[1][3]*primEdge[U1]*primEdge[U3] +
                              gcov[2][1]*primEdge[U2]*primEdge[U1] +
                              gcov[2][2]*primEdge[U2]*primEdge[U2] +
                              gcov[2][3]*primEdge[U2]*primEdge[U3] +
                              gcov[3][1]*primEdge[U3]*primEdge[U1] +
                              gcov[3][2]*primEdge[U3]*primEdge[U2] +
                              gcov[3][3]*primEdge[U3]*primEdge[U3];
                
                if (fabs(vSqr) < 1e-13) vSqr = 1e-13;
                if (fabs(vSqr) >= 1.) vSqr = 1. - 1./(GAMMA_MAX*GAMMA_MAX);

                gamma = 1./sqrt(1. - vSqr);
                primEdge[U1] *= gamma;
                primEdge[U2] *= gamma;
                primEdge[U3] *= gamma;
            }
            for (int var=0; var<DOF; var++)
                primTile[INDEX_LOCAL(TILE_SIZE_X1+iNg,jTile,var)] = primEdge[var];
        }
    }

    if (i<=0 && iTile==0) {
        for (int iNg=-NG; iNg<0; iNg++) {
            for (int var=0; var<DOF; var++) {
                primEdge[var] = 
                primTile[INDEX_LOCAL(iNg,jTile,var)];
            }

            X1 = i_TO_X1_CENTER(iNg); X2 = j_TO_X2_CENTER(j);
            
            gCovCalc(gcov, X1, X2);
            gDetCalc(&gdet, gcov);
            gConCalc(gcon, gcov, gdet);
            g = sqrt(-gdet);

            gammaCalc(&gamma, primEdge, gcov);
            setFloor(primEdge, gamma, X1, X2);
            gammaCalc(&gamma, primEdge, gcov);
            alphaCalc(&alpha, gcon);
            uconCalc(ucon, gamma, alpha, primEdge, gcon);

            if (ucon[1] > 0.) {

                primEdge[U1] /= gamma;
                primEdge[U2] /= gamma;
                primEdge[U3] /= gamma;

                alpha = 1./sqrt(-gcon[0][0]);
                
                primEdge[U1] = gcon[0][1]*alpha;

                double vSqr = gcov[1][1]*primEdge[U1]*primEdge[U1] +
                              gcov[1][2]*primEdge[U1]*primEdge[U2] + 
                              gcov[1][3]*primEdge[U1]*primEdge[U3] +
                              gcov[2][1]*primEdge[U2]*primEdge[U1] +
                              gcov[2][2]*primEdge[U2]*primEdge[U2] +
                              gcov[2][3]*primEdge[U2]*primEdge[U3] +
                              gcov[3][1]*primEdge[U3]*primEdge[U1] +
                              gcov[3][2]*primEdge[U3]*primEdge[U2] +
                              gcov[3][3]*primEdge[U3]*primEdge[U3];
                
                if (fabs(vSqr) < 1e-13) vSqr = 1e-13;
                if (fabs(vSqr) >= 1.) vSqr = 1. - 1./(GAMMA_MAX*GAMMA_MAX);

                gamma = 1./sqrt(1. - vSqr);
                primEdge[U1] *= gamma;
                primEdge[U2] *= gamma;
                    primEdge[U3] *= gamma;
            }
            for (int var=0; var<DOF; var++) {
                primTile[INDEX_LOCAL(iNg,jTile,var)] = primEdge[var];
            }
        }
    }
#endif
    
/*
     (i,j) 
     _____
    |     |
    |  o  | X1 = i_TO_X1_CENTER(i)
    |_____| X2 = j_TO_X2_CENTER(j)
     _____
    |     |
    |o    | X1 = i_TO_X1_FACE(i)
    |_____| X2 = j_TO_X2_CENTER(j)

     _____
    |     |
    |     | X1 = i_TO_X1_CENTER(i)
    |__o__| X2 = j_TO_X2_FACE(j)

     _____
    |     |
    |     | X1 = i_TO_X1_FACE(i)
    |o____| X2 = j_TO_X2_FACE(j)

*/
    for (int var=0; var<DOF; var++) {
        primEdge[var] = primTile[INDEX_LOCAL(iTile,jTile,var)];
    }

    X1 = i_TO_X1_CENTER(i); X2 = j_TO_X2_CENTER(j);
    gCovCalc(gcov, X1, X2);
    gDetCalc(&gdet, gcov);
    gConCalc(gcon, gcov, gdet);
    g = sqrt(-gdet);
    
    gammaCalc(&gamma, primEdge, gcov);
    setFloor(primEdge, gamma, X1, X2);
    gammaCalc(&gamma, primEdge, gcov);
    alphaCalc(&alpha, gcon);
    uconCalc(ucon, gamma, alpha, primEdge, gcon);
    covFromCon(ucov, ucon, gcov);
    bconCalc(bcon, primEdge, ucon, ucov);
    covFromCon(bcov, bcon, gcov);
    
    mhdCalc(mhd, primEdge, ucon, ucov, bcon, bcov, dvars_dt, primTile);

    addSources(dU_dt, 
               ucon, ucov, bcon, bcov, 
               gcon, gcov, mhd, primEdge, g, gdet, alpha,
               X1, X2, primTile,
               dvars_dt, ducon_dt,
               i, j, iTile, jTile);

/*
      i-2   i-1    i    i+1   i+2
     _____ _____ _____ _____ _____
    |  _  |  _  |  _  |     |     |
    |  o  |  o->|  o  |  o  |  o  |
    |_____|_____|_____|_____|_____|
*/

    ReconstructX1(primTile, iTile-1, jTile, primEdge, RIGHT);

    X1 = i_TO_X1_FACE(i); X2 = j_TO_X2_CENTER(j);
    gCovCalc(gcov, X1, X2);
    gDetCalc(&gdet, gcov);
    gConCalc(gcon, gcov, gdet);
    g = sqrt(-gdet);

    gammaCalc(&gamma, primEdge, gcov);
    setFloor(primEdge, gamma, X1, X2);
    gammaCalc(&gamma, primEdge, gcov);
    alphaCalc(&alpha, gcon);
    uconCalc(ucon, gamma, alpha, primEdge, gcon);
    covFromCon(ucov, ucon, gcov);
    bconCalc(bcon, primEdge, ucon, ucov);
    covFromCon(bcov, bcon, gcov);

    mhdCalc(mhd, primEdge, ucon, ucov, bcon, bcov, dvars_dt, primTile); 

    ComputeFluxAndU(fluxL, uL,
                    ucon, ucov,
                    bcon, bcov,
                    gcon, gcov,
                    mhd, primEdge, g, 1);

    conDotCov(&bsqr, bcon, bcov);
    VChar(&cminL, &cmaxL, ucon, ucov, bsqr, gcon, primEdge, 1); 

/*
      i-2   i-1    i    i+1   i+2
     _____ _____ _____ _____ _____
    |     |  _  |  _  |  _  |     |
    |  o  |  o  |<-o  |  o  |  o  |
    |_____|_____|_____|_____|_____|
*/

    ReconstructX1(primTile, iTile, jTile, primEdge, LEFT);

    X1 = i_TO_X1_FACE(i); X2 = j_TO_X2_CENTER(j);
    gCovCalc(gcov, X1, X2);
    gDetCalc(&gdet, gcov);
    gConCalc(gcon, gcov, gdet);
    g = sqrt(-gdet);

    gammaCalc(&gamma, primEdge, gcov);
    setFloor(primEdge, gamma, X1, X2);
    gammaCalc(&gamma, primEdge, gcov);
    alphaCalc(&alpha, gcon);
    uconCalc(ucon, gamma, alpha, primEdge, gcon);
    covFromCon(ucov, ucon, gcov);
    bconCalc(bcon, primEdge, ucon, ucov);
    covFromCon(bcov, bcon, gcov);

    mhdCalc(mhd, primEdge, ucon, ucov, bcon, bcov, dvars_dt, primTile);

    ComputeFluxAndU(fluxR, uR,
                    ucon, ucov,
                    bcon, bcov,
                    gcon, gcov,
                    mhd, primEdge, g, 1);

    conDotCov(&bsqr, bcon, bcov);
    VChar(&cminR, &cmaxR, ucon, ucov, bsqr, gcon, primEdge, 1); 

    RiemannSolver(fluxX1L, fluxL, fluxR, uL, uR,
                  cminL, cmaxL, cminR, cmaxR);


/*
      i-2   i-1    i    i+1   i+2
     _____ _____ _____ _____ _____
    |     |  _  |  _  |  _  |     |
    |  o  |  o  |  o->|  o  |  o  |
    |_____|_____|_____|_____|_____|
*/

    ReconstructX1(primTile, iTile, jTile, primEdge, RIGHT);

    X1 = i_TO_X1_FACE(i+1); X2 = j_TO_X2_CENTER(j);
    gCovCalc(gcov, X1, X2);
    gDetCalc(&gdet, gcov);
    gConCalc(gcon, gcov, gdet);
    g = sqrt(-gdet);

    gammaCalc(&gamma, primEdge, gcov);
    setFloor(primEdge, gamma, X1, X2);
    gammaCalc(&gamma, primEdge, gcov);
    alphaCalc(&alpha, gcon);
    uconCalc(ucon, gamma, alpha, primEdge, gcon);
    covFromCon(ucov, ucon, gcov);
    bconCalc(bcon, primEdge, ucon, ucov);
    covFromCon(bcov, bcon, gcov);

    mhdCalc(mhd, primEdge, ucon, ucov, bcon, bcov, dvars_dt, primTile);

    ComputeFluxAndU(fluxL, uL,
                    ucon, ucov,
                    bcon, bcov,
                    gcon, gcov,
                    mhd, primEdge, g, 1);

    conDotCov(&bsqr, bcon, bcov);
    VChar(&cminL, &cmaxL, ucon, ucov, bsqr, gcon, primEdge, 1); 

/*
      i-2   i-1    i    i+1   i+2
     _____ _____ _____ _____ _____
    |     |     |  _  |  _  |  _  |
    |  o  |  o  |  o  |<-o  |  o  |
    |_____|_____|_____|_____|_____|
*/

    ReconstructX1(primTile, iTile+1, jTile, primEdge, LEFT);

    X1 = i_TO_X1_FACE(i+1); X2 = j_TO_X2_CENTER(j);
    gCovCalc(gcov, X1, X2);
    gDetCalc(&gdet, gcov);
    gConCalc(gcon, gcov, gdet);
    g = sqrt(-gdet);

    gammaCalc(&gamma, primEdge, gcov);
    setFloor(primEdge, gamma, X1, X2);
    gammaCalc(&gamma, primEdge, gcov);
    alphaCalc(&alpha, gcon);
    uconCalc(ucon, gamma, alpha, primEdge, gcon);
    covFromCon(ucov, ucon, gcov);
    bconCalc(bcon, primEdge, ucon, ucov);
    covFromCon(bcov, bcon, gcov);

    mhdCalc(mhd, primEdge, ucon, ucov, bcon, bcov, dvars_dt, primTile);

    ComputeFluxAndU(fluxR, uR,
                    ucon, ucov,
                    bcon, bcov,
                    gcon, gcov,
                    mhd, primEdge, g, 1);

    conDotCov(&bsqr, bcon, bcov);
    VChar(&cminR, &cmaxR, ucon, ucov, bsqr, gcon, primEdge, 1); 

    RiemannSolver(fluxX1R, fluxL, fluxR, uL, uR,
                  cminL, cmaxL, cminR, cmaxR);

/*
                 _____ 
                |     |
                |  o  |  j+2
                |_____| 
                |     |
                |  o  |  j+1 
                |_____| 
                |     |
                | |o| |   j
                |_____|
                |  ^  |
                | |o| |  j-1
                |_____| 
                |     |
                | |o| |  j-2
                |_____|

*/

    ReconstructX2(primTile, iTile, jTile-1, primEdge, UP);

    X1 = i_TO_X1_CENTER(i); X2 = j_TO_X2_FACE(j);
    gCovCalc(gcov, X1, X2);
    gDetCalc(&gdet, gcov);
    gConCalc(gcon, gcov, gdet);
    g = sqrt(-gdet);

    gammaCalc(&gamma, primEdge, gcov);
    setFloor(primEdge, gamma, X1, X2);
    gammaCalc(&gamma, primEdge, gcov);
    alphaCalc(&alpha, gcon);
    uconCalc(ucon, gamma, alpha, primEdge, gcon);
    covFromCon(ucov, ucon, gcov);
    bconCalc(bcon, primEdge, ucon, ucov);
    covFromCon(bcov, bcon, gcov);

    mhdCalc(mhd, primEdge, ucon, ucov, bcon, bcov, dvars_dt, primTile);

    ComputeFluxAndU(fluxL, uL,
                    ucon, ucov,
                    bcon, bcov,
                    gcon, gcov,
                    mhd, primEdge, g, 2);

    conDotCov(&bsqr, bcon, bcov);
    VChar(&cminL, &cmaxL, ucon, ucov, bsqr, gcon, primEdge, 2); 

/*
                 _____ 
                |     |
                |  o  |  j+2
                |_____| 
                |     |
                | |o| |  j+1 
                |_____| 
                |     |
                | |o| |   j
                |__v__|
                |     |
                | |o| |  j-1
                |_____| 
                |     |
                |  o  |  j-2
                |_____|

*/

    ReconstructX2(primTile, iTile, jTile, primEdge, DOWN);

    X1 = i_TO_X1_CENTER(i); X2 = j_TO_X2_FACE(j);
    gCovCalc(gcov, X1, X2);
    gDetCalc(&gdet, gcov);
    gConCalc(gcon, gcov, gdet);
    g = sqrt(-gdet);

    gammaCalc(&gamma, primEdge, gcov);
    setFloor(primEdge, gamma, X1, X2);
    gammaCalc(&gamma, primEdge, gcov);
    alphaCalc(&alpha, gcon);
    uconCalc(ucon, gamma, alpha, primEdge, gcon);
    covFromCon(ucov, ucon, gcov);
    bconCalc(bcon, primEdge, ucon, ucov);
    covFromCon(bcov, bcon, gcov);

    mhdCalc(mhd, primEdge, ucon, ucov, bcon, bcov, dvars_dt, primTile);

    ComputeFluxAndU(fluxR, uR,
                    ucon, ucov,
                    bcon, bcov,
                    gcon, gcov,
                    mhd, primEdge, g, 2);

    conDotCov(&bsqr, bcon, bcov);
    VChar(&cminR, &cmaxR, ucon, ucov, bsqr, gcon, primEdge, 2); 

    RiemannSolver(fluxX2L, fluxL, fluxR, uL, uR,
                  cminL, cmaxL, cminR, cmaxR);


/*
                 _____ 
                |     |
                |  o  |  j+2
                |_____| 
                |     |
                | |o| |  j+1 
                |_____| 
                |  ^  |
                | |o| |   j
                |_____|
                |     |
                | |o| |  j-1
                |_____| 
                |     |
                |  o  |  j-2
                |_____|

*/

    ReconstructX2(primTile, iTile, jTile, primEdge, UP);

    X1 = i_TO_X1_CENTER(i); X2 = j_TO_X2_FACE(j+1);
    gCovCalc(gcov, X1, X2);
    gDetCalc(&gdet, gcov);
    gConCalc(gcon, gcov, gdet);
    g = sqrt(-gdet);

    gammaCalc(&gamma, primEdge, gcov);
    setFloor(primEdge, gamma, X1, X2);
    gammaCalc(&gamma, primEdge, gcov);
    alphaCalc(&alpha, gcon);
    uconCalc(ucon, gamma, alpha, primEdge, gcon);
    covFromCon(ucov, ucon, gcov);
    bconCalc(bcon, primEdge, ucon, ucov);
    covFromCon(bcov, bcon, gcov);

    mhdCalc(mhd, primEdge, ucon, ucov, bcon, bcov, dvars_dt, primTile);

    ComputeFluxAndU(fluxL, uL,
                    ucon, ucov,
                    bcon, bcov,
                    gcon, gcov,
                    mhd, primEdge, g, 2);

    conDotCov(&bsqr, bcon, bcov);
    VChar(&cminL, &cmaxL, ucon, ucov, bsqr, gcon, primEdge, 2); 

/*
                 _____ 
                |     |
                | |o| |  j+2
                |_____| 
                |     |
                | |o| |  j+1 
                |__V__| 
                |     |
                | |o| |   j
                |_____|
                |     |
                |  o  |  j-1
                |_____| 
                |     |
                |  o  |  j-2
                |_____|

*/

    ReconstructX2(primTile, iTile, jTile+1, primEdge, DOWN);

    X1 = i_TO_X1_CENTER(i); X2 = j_TO_X2_FACE(j+1);
    gCovCalc(gcov, X1, X2);
    gDetCalc(&gdet, gcov);
    gConCalc(gcon, gcov, gdet);
    g = sqrt(-gdet);

    gammaCalc(&gamma, primEdge, gcov);
    setFloor(primEdge, gamma, X1, X2);
    gammaCalc(&gamma, primEdge, gcov);
    alphaCalc(&alpha, gcon);
    uconCalc(ucon, gamma, alpha, primEdge, gcon);
    covFromCon(ucov, ucon, gcov);
    bconCalc(bcon, primEdge, ucon, ucov);
    covFromCon(bcov, bcon, gcov);

    mhdCalc(mhd, primEdge, ucon, ucov, bcon, bcov, dvars_dt, primTile);

    ComputeFluxAndU(fluxR, uR,
                    ucon, ucov,
                    bcon, bcov,
                    gcon, gcov,
                    mhd, primEdge, g, 2);

    conDotCov(&bsqr, bcon, bcov);
    VChar(&cminR, &cmaxR, ucon, ucov, bsqr, gcon, primEdge, 2); 

    RiemannSolver(fluxX2R, fluxL, fluxR, uL, uR,
                  cminL, cmaxL, cminR, cmaxR);


    for (int var=0; var<DOF; var++) {
        
        if (i>=0 && i<N1 && j>=0 && j<N2) {
            F[INDEX_GLOBAL(i,j,var)] = dU_dt[var];
        }

        fluxX1[INDEX_GLOBAL_WITH_NG(i,j,var)] = fluxX1L[var];
        fluxX2[INDEX_GLOBAL_WITH_NG(i,j,var)] = fluxX2L[var];

        if (i==N1-1) {
            fluxX1[INDEX_GLOBAL_WITH_NG(i+1, j, var)] = fluxX1R[var];
        }

        if (j==N2-1) {
            fluxX2[INDEX_GLOBAL_WITH_NG(i, j+1, var)] = fluxX2R[var];
        }
    }

}
