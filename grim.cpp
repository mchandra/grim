#include "grim.h"

struct data
{
  REAL primBoundaries[(N1 + 2*NG)*(N2 + 2*NG)*DOF];
};

PetscErrorCode SetCoordinates(DM dmda)
{
    DM coordDM;
    Vec coordVec;

    int X1Start, X2Start;
    int X1Size, X2Size;

    DMDAGetCorners(dmda, 
                   &X1Start, &X2Start, NULL,
                   &X1Size, &X2Size, NULL);

    DMDASetUniformCoordinates(dmda,
                              0., 1., 0., 1., 0., 1.);
    DMGetCoordinateDM(dmda, &coordDM);
    DMGetCoordinatesLocal(dmda, &coordVec);

    PetscScalar X1, X2, r, theta;
    DMDACoor2d **coord;
    
    DMDAVecGetArray(coordDM, coordVec, &coord);

    for (int j=X2Start; j<X2Start+X2Size; j++) {
        for (int i=X1Start; i<X1Start+X1Size; i++) {
    
            X1 = i_TO_X1_CENTER(i);
            X2 = j_TO_X2_CENTER(j);

            BLCoords(&r, &theta, X1, X2);

            coord[j][i].x = r*PetscSinScalar(theta);
            coord[j][i].y = r*PetscCosScalar(theta);
        }
    }
    DMDAVecRestoreArray(coordDM, coordVec, &coord);

    DMSetCoordinatesLocal(dmda, coordVec);

    return (0);
}

int main(int argc, char **argv)
{
    TS ts;
    SNES snes;
    Vec soln;
    DM dmda;

    int X1Start, X2Start;
    int X1Size, X2Size;

    PetscInitialize(&argc, &argv, PETSC_NULL, help);

    DMDACreate2d(PETSC_COMM_WORLD, 
                 DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                 DMDA_STENCIL_BOX,
                 N1, N2,
                 PETSC_DECIDE, PETSC_DECIDE,
                 DOF, NG, PETSC_NULL, PETSC_NULL, &dmda);

    DMDAGetCorners(dmda, 
                   &X1Start, &X2Start, NULL,
                   &X1Size, &X2Size, NULL);

    SetCoordinates(dmda);
//    DMDASetUniformCoordinates(dmda, 
//                              DX1/2., 1. - DX1/2.,
//                              DX2/2., 1. - DX2/2., 
//                              0., 0.);


    DMDASetFieldName(dmda, RHO, "Density");
    DMDASetFieldName(dmda, UU, "Internal Energy");
    DMDASetFieldName(dmda, U1, "Ur");
    DMDASetFieldName(dmda, U2, "Utheta");
    DMDASetFieldName(dmda, U3, "Uphi");
    DMDASetFieldName(dmda, B1, "Br");
    DMDASetFieldName(dmda, B2, "Btheta");
    DMDASetFieldName(dmda, B3, "Bphi");    

    DMCreateGlobalVector(dmda, &soln);

    TSCreate(PETSC_COMM_WORLD, &ts);
    TSSetDM(ts, dmda);

    clErr = cl::Platform::get(&platforms);
    CheckCLErrors(clErr, "cl::Platform::get");

    clErr = platforms.at(0).getDevices(CL_DEVICE_TYPE_CPU, &devices);
    CheckCLErrors(clErr, "cl::Platform::getDevices");

    context = cl::Context(devices, NULL, NULL, NULL, &clErr);
    CheckCLErrors(clErr, "cl::Context::Context");

    queue = cl::CommandQueue(context, devices.at(0), 0, &clErr);
    CheckCLErrors(clErr, "cl::CommandQueue::CommandQueue");

    std::ifstream sourceFile("computeresidual.cl");
    std::string sourceCode((std::istreambuf_iterator<char>(sourceFile)),
                            std::istreambuf_iterator<char>());
    cl::Program::Sources source(1, std::make_pair(sourceCode.c_str(),
                                sourceCode.length()+1));
    
    program = cl::Program(context, source, &clErr);
    CheckCLErrors(clErr, "cl::Program::Program");

    std::string BuildOptions("\
                              -D X1_SIZE=" +
                             std::to_string(X1Size) +
                             " -D X2_SIZE=" + 
                             std::to_string(X2Size) +
                             " -D TOTAL_X1_SIZE=" + 
                             std::to_string(X1Size+2*NG) + 
                             " -D TOTAL_X2_SIZE=" +
                             std::to_string(X2Size+2*NG) +
                             " -DOPENCL" + 
                             " -DLEFT_BOUNDARY=" + 
                             std::to_string(LEFT_BOUNDARY) +
                             " -DRIGHT_BOUNDARY=" + 
                             std::to_string(RIGHT_BOUNDARY) +
                             " -DBOTTOM_BOUNDARY=" + 
                             std::to_string(BOTTOM_BOUNDARY) +
                             " -DTOP_BOUNDARY=" + 
                             std::to_string(TOP_BOUNDARY) +
                             " -DINFLOW_CHECK=" + 
                             std::to_string(INFLOW_CHECK));

    PetscScalar start = std::clock();
    clErr = program.build(devices, BuildOptions.c_str(), NULL, NULL);
    const char *buildlog = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(
                                                devices.at(0),
                                                &clErr).c_str();
    PetscPrintf(PETSC_COMM_WORLD, "%s\n", buildlog);
    CheckCLErrors(clErr, "cl::Program::build");
    PetscScalar end = std::clock();

    PetscScalar time = (end - start)/(PetscScalar)CLOCKS_PER_SEC;
    PetscPrintf(PETSC_COMM_WORLD, 
                "Time taken for kernel compilation = %f\n", time);

    kernel = cl::Kernel(program, "ComputeResidual", &clErr);
    CheckCLErrors(clErr, "cl::Kernel::Kernel");

    cl_ulong localMemSize = kernel.getWorkGroupInfo<CL_KERNEL_LOCAL_MEM_SIZE>(
                                        devices.at(0), &clErr);
    cl_ulong privateMemSize = kernel.getWorkGroupInfo<CL_KERNEL_PRIVATE_MEM_SIZE>(
                                        devices.at(0), &clErr);
    printf("Local memory used = %llu\n", (unsigned long long)localMemSize);
    printf("Private memory used = %llu\n", (unsigned long long)privateMemSize);

    struct data tsData;
    InitialConditionMTITest(ts, soln, &tsData);

    PetscViewer viewer;
#if(RESTART)
    PetscViewerHDF5Open(PETSC_COMM_WORLD,"restartfile.h5",
                        FILE_MODE_READ, &viewer);
    PetscObjectSetName((PetscObject) soln,"soln");
    VecLoad(soln, viewer);
    PetscViewerDestroy(&viewer);
#endif /* Restart from "init.h5" */

    Monitor(ts, 0., 0., soln, NULL);

    TSSetIFunction(ts, PETSC_NULL, ComputeResidual, &tsData);

    TSSetSolution(ts, soln);
    TSMonitorSet(ts, Monitor, NULL, NULL);
    TSGetSNES(ts, &snes);
    SNESMonitorSet(snes, SNESMonitor, NULL, NULL);
    TSSetType(ts, TSTHETA);
    TSSetFromOptions(ts);

    TSSolve(ts, soln);

    TSDestroy(&ts);
    VecDestroy(&soln);
    DMDestroy(&dmda);

    PetscFinalize();
    return(0);
}

PetscErrorCode ComputeResidual(TS ts,
                               PetscScalar t,
                               Vec Prim, Vec dPrim_dt,
                               Vec F, void *ptr)
{
    struct data *tsData = (struct data*)ptr;
    
    PetscScalar *prim, *dprim_dt, *f;
    VecGetArray(Prim, &prim);
    VecGetArray(dPrim_dt, &dprim_dt);
    VecGetArray(F, &f);

    REAL fluxX1[(N1 + 2*NG)*(N2 + 2*NG)*DOF];
    REAL fluxX2[(N1 + 2*NG)*(N2 + 2*NG)*DOF];

    cl::Buffer primBuffer, dprimBuffer_dt, fbuffer, primBoundariesBuffer;
    PetscInt size = DOF*N1*N2*sizeof(PetscScalar);
    PetscInt sizeWithNG = DOF*(N1 + 2*NG)*(N2 + 2*NG)*sizeof(PetscScalar);


    cl::Buffer fluxX1Buffer, fluxX2Buffer;

    primBuffer = cl::Buffer(context,
                            CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
                            size, &(prim[0]), &clErr);
    primBoundariesBuffer = cl::Buffer(context,
                                      CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
                                      sizeWithNG, &(tsData->primBoundaries[0]),
                                      &clErr);
    dprimBuffer_dt = cl::Buffer(context,
                                CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
                                size, &(dprim_dt[0]), &clErr);
    fbuffer = cl::Buffer(context,
                         CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY,
                         size, &(f[0]), &clErr);

    fluxX1Buffer = cl::Buffer(context,
                              CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY,
                              sizeWithNG, &(fluxX1[0]), &clErr);
    fluxX2Buffer = cl::Buffer(context,
                              CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY,
                              sizeWithNG, &(fluxX2[0]), &clErr);



    clErr = kernel.setArg(0, primBuffer);
    clErr = kernel.setArg(1, dprimBuffer_dt);
    clErr = kernel.setArg(2, fbuffer);
    clErr = kernel.setArg(3, fluxX1Buffer);
    clErr = kernel.setArg(4, fluxX2Buffer);
    clErr = kernel.setArg(5, primBoundariesBuffer);

    cl::NDRange global(N1, N2);
    cl::NDRange local(TILE_SIZE_X1, TILE_SIZE_X2);
    clErr = queue.enqueueNDRangeKernel(kernel,
                                       cl::NullRange,
                                       global, local,
                                       NULL, NULL);

    f = (PetscScalar*)queue.enqueueMapBuffer(fbuffer,
                                             CL_FALSE,
                                             CL_MAP_READ,
                                             0, size,
                                             NULL, NULL, &clErr);

    clErr = queue.finish();

    REAL emf[N2+2*NG][N1+2*NG];

    for (int j=0; j<N2+NG; j++) {
        for (int i=0; i<N1+NG; i++) {
            emf[j][i] = 
                        0.25*(fluxX1[INDEX_GLOBAL_WITH_NG(i,j,B2)] +
                              fluxX1[INDEX_GLOBAL_WITH_NG(i,j-1,B2)] -
                              fluxX2[INDEX_GLOBAL_WITH_NG(i,j,B1)] -
                              fluxX2[INDEX_GLOBAL_WITH_NG(i-1,j,B1)]);
        }
    }

    for (int j=1; j<N2-1; j++) {
        for (int i=1; i<N1-1; i++) {
            
            fluxX1[INDEX_GLOBAL_WITH_NG(i,j,B1)] = 0.;

            fluxX1[INDEX_GLOBAL_WITH_NG(i,j,B2)] = 
                                    0.5*(emf[j][i] + emf[j+1][i]);

            fluxX2[INDEX_GLOBAL_WITH_NG(i,j,B1)] =
                                    -0.5*(emf[j][i] + emf[j][i+1]);

            fluxX2[INDEX_GLOBAL_WITH_NG(i,j,B2)] = 0.;

        }
    }

    for (int j=0; j<N2; j++) {
        for (int i=0; i<N1; i++) {
            for (int var=0; var<DOF; var++) {

                f[INDEX_GLOBAL(i,j,var)] = f[INDEX_GLOBAL(i,j,var)] +
                                (fluxX1[INDEX_GLOBAL_WITH_NG(i+1,j,var)] -
                                 fluxX1[INDEX_GLOBAL_WITH_NG(i,j,var)])/DX1 +
                                (fluxX2[INDEX_GLOBAL_WITH_NG(i,j+1,var)] -
                                 fluxX2[INDEX_GLOBAL_WITH_NG(i,j,var)])/DX2;
            }
        }
    }

    VecRestoreArray(Prim, &prim);
    VecRestoreArray(dPrim_dt, &dprim_dt);
    VecRestoreArray(F, &f);

    return(0.);
}


void Benchmark(TS ts, Vec Prim)
{
    PetscInt NIter = 1;
    std::clock_t start, end;
    PetscScalar time;

    PetscScalar t = 0.;
    Vec dPrim_dt, F;
    VecDuplicate(Prim, &dPrim_dt);
    VecDuplicate(Prim, &F);
    VecSet(dPrim_dt, 0.);
    VecSet(F, 0.);

    start = std::clock();
    for (int n=0; n < NIter; n++) {
        ComputeResidual(ts, t, Prim, dPrim_dt, F, NULL);
    }
    end = std::clock();

    time = (end - start)/(PetscScalar)CLOCKS_PER_SEC;
    PetscPrintf(PETSC_COMM_WORLD, "Time taken for %d iterations = %f\n", NIter, time);

    PetscViewer viewer;
    PetscViewerHDF5Open(PETSC_COMM_WORLD,"InitialResidual.h5",
                        FILE_MODE_WRITE, &viewer);
    PetscObjectSetName((PetscObject) F, "F");
    VecView(F, viewer);
    PetscViewerDestroy(&viewer);

    VecDestroy(&dPrim_dt);
    VecDestroy(&F);
}

void InitialConditionTest(TS ts, Vec X)
{
    DM dmda;
    TSGetDM(ts, &dmda);
    PetscScalar ***x;

    DMDAVecGetArrayDOF(dmda, X, &x);

    for (PetscInt j=0; j<N2; j++) {
        for (PetscInt i=0; i<N1; i++) {
            PetscScalar xcoord = DX1/2. + i*DX1;
            PetscScalar ycoord = DX2/2. + j*DX2;

            PetscScalar xcenter = 0.5;
            PetscScalar ycenter = 0.5;

            PetscScalar r = PetscSqrtScalar((xcoord-xcenter)*(xcoord-xcenter)+
                                            (ycoord-ycenter)*(ycoord-ycenter));
            
//            x[j][i][RHO] = 1. + exp(-r*r/0.01);
//
//            x[j][i][UU] = 1./(ADIABATIC_INDEX-1);
//            x[j][i][U1] = 4.95;
//            x[j][i][U2] = 4.95;
//            x[j][i][U3] = 0.;
//            x[j][i][B1] = 0.;
//            x[j][i][B2] = 0.;
//            x[j][i][B3] = 0.;

            x[j][i][RHO] = 25./(36.*M_PI);
            x[j][i][UU] = 5./(12.*M_PI*(ADIABATIC_INDEX - 1.));
            PetscScalar vx = -0.5*PetscSinScalar(2*M_PI*ycoord);
            PetscScalar vy = 0.5*PetscSinScalar(2*M_PI*xcoord);
            PetscScalar vz = 0.;
            PetscScalar gamma = 1./PetscSqrtScalar(1 - vx*vx - vy*vy - vz*vz);
            x[j][i][U1] = gamma*vx;
            x[j][i][U2] = gamma*vy;
            x[j][i][U3] = gamma*vz;
            x[j][i][B1] =
                -1./PetscSqrtScalar(4*M_PI)*PetscSinScalar(2*M_PI*ycoord);
            x[j][i][B2] =
                1./PetscSqrtScalar(4*M_PI)*PetscSinScalar(4*M_PI*xcoord);
            x[j][i][B3] = 0.0;

//                  Komissarov strong cylindrical explosion
//              PetscScalar R = 1.;
//              PetscScalar rho_outside = 1e-4, rho_inside = 1e-2;
//              PetscScalar u_outside = 3e-5/(5./3 - 1); 
//              PetscScalar u_inside = 1./(5./3 - 1);
//              PetscScalar alpha = 20.;
//              PetscScalar norm_factor = (1. - tanh(-R))/2.;
//
//              x[j][i][RHO] = ((rho_inside - rho_outside)*
//                              (1. - tanh(pow(r, alpha)-R))
//                              /2./norm_factor + rho_outside);
//
//              x[j][i][UU] = ((u_inside - u_outside)*
//                            (1. - tanh(pow(r, alpha)-R))
//                            /2./norm_factor + u_outside);
//
//              PetscScalar vx = 0.;
//              PetscScalar vy = 0.;
//              PetscScalar vz = 0.;
//              PetscScalar gamma = 1.;
//              x[j][i][U1] = gamma*vx;
//              x[j][i][U2] = gamma*vy;
//              x[j][i][U3] = gamma*vz;
//              x[j][i][B1] = .1;
//              x[j][i][B2] = 0.;
//              x[j][i][B3] = 0.0;

        }
    }

    DMDAVecRestoreArrayDOF(dmda, X, &x);
}

struct rootFinderParams
{
  REAL C1, C2, r;
};

REAL inputFunctionForRootFinder(REAL T, void *ptr)
{
  struct rootFinderParams *params = 
    (struct rootFinderParams *)ptr;

  REAL C1 = params->C1;
  REAL C2 = params->C2;
  REAL r  = params->r;

  REAL n = 1./(ADIABATIC_INDEX-1.);

  return pow(1. + (1. + n)*T, 2.)*(1 - 2./r + C1*C1/(pow(r, 4.) * pow(T, 2.*n)) ) - C2;
}

REAL inputDerForRootFinder(REAL T, void *ptr)
{
  struct rootFinderParams *params = 
    (struct rootFinderParams *)ptr;

  REAL C1 = params->C1;
  REAL r  = params->r;

  REAL n = 1./(ADIABATIC_INDEX-1.);

  REAL ans;
  ans  =  2.*(1 + n)*(1. + (1.+n)*T)*(1 - 2./r 
                                     + C1*C1/(pow(r, 4.) * pow(T, 2.*n))
                                    )
         + pow(1 + (1 + n)*T, 2.)*(-2*n*C1*C1/(pow(r, 4.) * pow(T, 2*n+1)) );

  return ans;
}

void inputFunctionAndDerForRootFinder(REAL T, void *ptr, REAL *y, REAL *dy)
{
  struct rootFinderParams *params = 
    (struct rootFinderParams *)ptr;

  REAL C1 = params->C1;
  REAL C2 = params->C2;
  REAL r  = params->r;

  REAL n = 1./(ADIABATIC_INDEX-1.);

  *y = pow(1. + (1. + n)*T, 2.)*(1 - 2./r + C1*C1/(pow(r, 4.) * pow(T, 2.*n)) ) - C2;

  *dy = 2.*(1 + n)*(1. + (1.+n)*T)*(1 - 2./r 
                                    + C1*C1/(pow(r, 4.) * pow(T, 2.*n))
                                   )
        + pow(1 + (1 + n)*T, 2.)*(-2*n*C1*C1/(pow(r, 4.) * pow(T, 2*n+1)) );
}

void InitialConditionMTITest(TS ts, Vec Prim, struct data *tsData)
{
  DM dmda;

  int X1Start, X2Start;
  int X1Size, X2Size;

  TSGetDM(ts, &dmda);

  DMDAGetCorners(dmda, 
                 &X1Start, &X2Start, NULL,
                 &X1Size, &X2Size, NULL);

  Vec localPrim;
  REAL ***prim;
  DMGetLocalVector(dmda, &localPrim);
  DMDAVecGetArrayDOF(dmda, localPrim, &prim);

  REAL r_c = 8.;
  REAL u_c = -sqrt(1./(2.*r_c));
  REAL V_c = sqrt(u_c*u_c/(1. - (3*u_c*u_c)) );

  REAL n = 1./(ADIABATIC_INDEX - 1.);

  REAL T_c = (n*V_c*V_c)/( (1. + n)*(1. - (n*V_c*V_c)) );

  REAL C1 = pow(T_c, n) * u_c * r_c*r_c;

  REAL C2 = pow(1 + (1 + n)*T_c, 2.)*(1. - 2./r_c + u_c*u_c);

  FILE *rCoords;
  rCoords = fopen("rCoords.txt", "w");

  for (int i=X1Start; i<X1Start+X1Size; i++) {
    REAL X1 = i_TO_X1_CENTER(i);
    REAL X2 = j_TO_X2_CENTER(0);

    REAL r, theta;

    BLCoords(&r, &theta, X1, X2);

    fprintf(rCoords, "%f\n", r);
  }
  fclose(rCoords);

  for (int j=X2Start-NG; j<X2Start+X2Size+NG; j++) {
    for (int i=X1Start-NG; i<X1Start+X1Size+NG; i++) {
  
      REAL X1 = i_TO_X1_CENTER(i);
      REAL X2 = j_TO_X2_CENTER(j);

      REAL r, theta;

      BLCoords(&r, &theta, X1, X2);

      gsl_function_fdf fdf;
      struct rootFinderParams params = {C1, C2, r};

      fdf.f = &inputFunctionForRootFinder;
      fdf.df = &inputDerForRootFinder;
      fdf.fdf = &inputFunctionAndDerForRootFinder;
      fdf.params = &params;

      const gsl_root_fdfsolver_type *solverType;
      gsl_root_fdfsolver *solver;

      solverType = gsl_root_fdfsolver_newton;
      solver = gsl_root_fdfsolver_alloc(solverType);

      gsl_root_fdfsolver_set(solver, &fdf, T_c);

      int status = GSL_CONTINUE;
      int maxIter = 10000, iter=0;

      REAL T = T_c, T0;
      do
        {
          printf("iter = %d, T = %f\n", iter, T);
          iter++;
          status = gsl_root_fdfsolver_iterate(solver);
          T0 = T;
          T = gsl_root_fdfsolver_root(solver);
          status = gsl_root_test_delta (T, T0, 0, 1e-10);

          if (status == GSL_SUCCESS)
            printf ("Converged for r = %f\n", r);

        }
      while (status == GSL_CONTINUE && iter < maxIter);

      gsl_root_fdfsolver_free(solver);
      /* End of root finding. We now have T */

      prim[j][i][RHO] = pow(T, n);
      prim[j][i][UU] = T*prim[j][i][RHO]/(ADIABATIC_INDEX-1.);

      REAL uConBL[NDIM], uConMKS[NDIM];
      uConBL[1] = C1/(pow(T, n)*r*r);
      uConBL[2] = 0.;
      uConBL[3] = 0.;

      transformBLtoMKS(uConBL, uConMKS, X1, X2, r, theta);

      REAL gcov[NDIM][NDIM], gcon[NDIM][NDIM];
      REAL gdet, alpha;
      gCovCalc(gcov, X1, X2);
      gDetCalc(&gdet, gcov);
      gConCalc(gcon, gcov, gdet);
      alphaCalc(&alpha, gcon);

      prim[j][i][U1] = uConMKS[1] + 
                       alpha*alpha*uConMKS[0]*gcon[0][1];
      prim[j][i][U2] = uConMKS[2] + 
                       alpha*alpha*uConMKS[0]*gcon[0][2];
      prim[j][i][U3] = uConMKS[3] +
                       alpha*alpha*uConMKS[0]*gcon[0][3];

      prim[j][i][B1] = 0.;
      prim[j][i][B2] = 0.;
      prim[j][i][B3] = 0.;

      for (int var=0; var<DOF; var++) {
        tsData->primBoundaries[INDEX_GLOBAL_WITH_NG(i,j,var)] = prim[j][i][var];
      }
    }
  }

  DMLocalToGlobalBegin(dmda, localPrim, INSERT_VALUES, Prim);
  DMLocalToGlobalEnd(dmda, localPrim, INSERT_VALUES, Prim);

  DMDAVecRestoreArrayDOF(dmda, localPrim, &prim);
  DMRestoreLocalVector(dmda, &localPrim);

}

void transformBLtoMKS(REAL uconBL[NDIM], REAL uconMKS[NDIM], 
                      REAL X1, REAL X2, REAL r, REAL theta)
{
  REAL gcovBL[NDIM][NDIM], gconBL[NDIM][NDIM], transform[NDIM][NDIM];

  for (int ii=0; ii<NDIM; ii++)
      for (int jj=0; jj<NDIM; jj++) {
          gcovBL[ii][jj] = 0.;
          gconBL[ii][jj] = 0.;
          transform[ii][jj] = 0.;
      }

  REAL DD = 1. - 2./r + A_SPIN*A_SPIN/(r*r);
  REAL mu = 1 + A_SPIN*A_SPIN*cos(theta)*cos(theta)/(r*r);

  gcovBL[0][0] = -(1. - 2./(r*mu));
  gcovBL[0][3] = -2.*A_SPIN*sin(theta)*sin(theta)/(r*mu);
  gcovBL[3][0] = gcovBL[0][3];
  gcovBL[1][1] = mu/DD;
  gcovBL[2][2] = r*r*mu;
  gcovBL[3][3] = r*r*sin(theta)*sin(theta)*\
                 (1. + A_SPIN*A_SPIN/(r*r) + 
                  2.*A_SPIN*A_SPIN*sin(theta)*sin(theta)/\
                  (r*r*r*mu));

  gconBL[0][0] = -1. -2.*(1 + A_SPIN*A_SPIN/(r*r))/(r*DD*mu);
  gconBL[0][3] = -2.*A_SPIN/(r*r*r*DD*mu);
  gconBL[3][0] = gconBL[0][3];
  gconBL[1][1] = DD/mu;
  gconBL[2][2] = 1./(r*r*mu);
  gconBL[3][3] = (1. - 2./(r*mu))/(r*r*sin(theta)*sin(theta)*DD);

  transform[0][0] = 1.;
  transform[1][1] = 1.;
  transform[2][2] = 1.;
  transform[3][3] = 1.;
  transform[0][1] = 2.*r/(r*r* - 2.*r + A_SPIN*A_SPIN);
  transform[3][1] = A_SPIN/(r*r - 2.*r + A_SPIN*A_SPIN);

  REAL AA = gcovBL[0][0];
  REAL BB = 2.*(gcovBL[0][1]*uconBL[1] + 
                gcovBL[0][2]*uconBL[2] +
                gcovBL[0][3]*uconBL[3]);
  REAL CC = 1. + gcovBL[1][1]*uconBL[1]*uconBL[1] +
                 gcovBL[2][2]*uconBL[2]*uconBL[2] +
                 gcovBL[3][3]*uconBL[3]*uconBL[3] +
             2.*(gcovBL[1][2]*uconBL[1]*uconBL[2] +
                 gcovBL[1][3]*uconBL[1]*uconBL[3] +
                 gcovBL[2][3]*uconBL[2]*uconBL[3]);

  REAL discriminent = BB*BB - 4.*AA*CC;
  uconBL[0] = -(BB + sqrt(discriminent))/(2.*AA);

  REAL uconKS[NDIM];
  uconKS[0] = transform[0][0]*uconBL[0] + 
              transform[0][1]*uconBL[1] +
              transform[0][2]*uconBL[2] +
              transform[0][3]*uconBL[3];

  uconKS[1] = transform[1][0]*uconBL[0] + 
              transform[1][1]*uconBL[1] +
              transform[1][2]*uconBL[2] +
              transform[1][3]*uconBL[3];

  uconKS[2] = transform[2][0]*uconBL[0] + 
              transform[2][1]*uconBL[1] +
              transform[2][2]*uconBL[2] +
              transform[2][3]*uconBL[3];

  uconKS[3] = transform[3][0]*uconBL[0] + 
              transform[3][1]*uconBL[1] +
              transform[3][2]*uconBL[2] +
              transform[3][3]*uconBL[3];

  PetscScalar rFactor, hFactor;
  rFactor = r - R0;
  hFactor = M_PI + (1. - H_SLOPE)*M_PI*cos(2.*M_PI*X2);
  uconMKS[0] = uconKS[0];
  uconMKS[1] = uconKS[1]*(1./rFactor);
  uconMKS[2] = uconKS[2]*(1./hFactor);
  uconMKS[3] = uconKS[3];

}
void InitialCondition(TS ts, Vec Prim)
{
    DM dmda;
    PetscScalar ***prim, r, theta;
    PetscScalar X1, X2;

    PetscScalar l, delta, sigma, A, lnOfh;
    PetscScalar thetaIn, deltaIn, sigmaIn, AIn;
    PetscScalar uPhiTmp, uconBL[NDIM], uconKS[NDIM], uconKSPrime[NDIM]; 
    PetscScalar gcovBL[NDIM][NDIM], gconBL[NDIM][NDIM], transform[NDIM][NDIM];
    PetscScalar rhoMax = 0., uMax = 0., bsqrMax = 0., rhoAv=0.;
    PetscScalar q, norm, betaActual;
    PetscScalar AA, BB, CC, DD, discriminent, mu;

    int X1Start, X2Start;
    int X1Size, X2Size;

    PetscScalar randNum;
    PetscRandom randNumGen;
    PetscRandomCreate(PETSC_COMM_SELF, &randNumGen);
    PetscRandomSetType(randNumGen, PETSCRAND48);

    TSGetDM(ts, &dmda);

    DMDAGetCorners(dmda, 
                   &X1Start, &X2Start, NULL,
                   &X1Size, &X2Size, NULL);

    PetscScalar AVector[X2Size+2*NG][X1Size+2*NG];

    Vec localPrim;
    DMGetLocalVector(dmda, &localPrim);
    DMDAVecGetArrayDOF(dmda, localPrim, &prim);

#if(CONDUCTION)
    for (int j=X2Start-2; j<X2Start+X2Size+2; j++)
        for (int i=X1Start-2; i<X1Start+X1Size+2; i++) {
            prim[j][i][FF] = 0.;
        }
#endif /* Conduction */

	l =  ( ( (pow(A_SPIN, 2) - 2.*A_SPIN*sqrt(R_MAX) + pow(R_MAX, 2.)) *
		     ( (-2.*A_SPIN*R_MAX*(pow(A_SPIN, 2.) - 2.*A_SPIN*sqrt(R_MAX) +
		        pow(R_MAX, 2.))) / sqrt(2.*A_SPIN*sqrt(R_MAX) + (-3. +
                R_MAX)*R_MAX) +
		       ((A_SPIN + (-2. + R_MAX)*sqrt(R_MAX))*(pow(R_MAX, 3) +
				pow(A_SPIN, 2)*(2. + R_MAX))) / sqrt(1 + (2.*A_SPIN) /
			   pow(R_MAX, 1.5) - 3./R_MAX) ) ) / \
           (pow(R_MAX, 3)*sqrt(2.*A_SPIN*sqrt(R_MAX) +
           (-3.+R_MAX)*R_MAX)*(pow(A_SPIN,2) + (-2.+R_MAX)*R_MAX)) );

    for (int j=X2Start-2; j<X2Start+X2Size+2; j++)
        for (int i=X1Start-2; i<X1Start+X1Size+2; i++) {
            
            X1 = i_TO_X1_CENTER(i);
            X2 = j_TO_X2_CENTER(j);

            BLCoords(&r, &theta, X1, X2);

            delta = r*r - 2*M*r + A_SPIN*A_SPIN;
            sigma = r*r + A_SPIN*A_SPIN*cos(theta)*cos(theta);
            A = (r*r + A_SPIN*A_SPIN)*(r*r + A_SPIN*A_SPIN) - 
                delta*A_SPIN*A_SPIN*sin(theta)*sin(theta);

            thetaIn = M_PI/2.;

            deltaIn = R_MIN*R_MIN - 2*M*R_MIN + A_SPIN*A_SPIN;
            sigmaIn = R_MIN*R_MIN + A_SPIN*A_SPIN*cos(thetaIn)*cos(thetaIn);
            AIn = (R_MIN*R_MIN + A_SPIN*A_SPIN)*(R_MIN*R_MIN + A_SPIN*A_SPIN) - 
                  deltaIn*A_SPIN*A_SPIN*sin(thetaIn)*sin(thetaIn);

            if (r >=R_MIN) {

			    lnOfh = 0.5*log((1. + sqrt(1. + 4.*(l*l*sigma*sigma)*delta/\
                                (A*sin(theta)*A*sin(theta))))/(sigma*delta/A))-
			            0.5*sqrt(1. + 4.*(l*l*sigma*sigma)*delta /
					             (A*A*sin(theta)*sin(theta))) -2.*A_SPIN*r*l/A-
			           (0.5*log((1. + sqrt(1. +
				                           4.*(l*l*sigmaIn*sigmaIn)*deltaIn/ \
				        (AIn*AIn*sin(thetaIn)*sin(thetaIn)))) /
				        (sigmaIn * deltaIn / AIn)) - 0.5 * sqrt(1. +
					    4.*(l*l*sigmaIn*sigmaIn)*deltaIn/ \
                        (AIn*AIn*sin(thetaIn)*sin(thetaIn))) - 
                        2.*A_SPIN*R_MIN*l/AIn);
		    } else {
			    lnOfh = 1.;
            }

            if (lnOfh <0. || r < R_MIN) {

                prim[j][i][RHO] = 1e-7*RHO_MIN;
                prim[j][i][UU] = 1e-7*U_MIN;
                prim[j][i][U1] = 0.;
                prim[j][i][U2] = 0.;
                prim[j][i][U3] = 0.;

            } else {

                prim[j][i][RHO] = pow((exp(lnOfh) - 1.)*
                    (ADIABATIC_INDEX-1.)/(KAPPA*ADIABATIC_INDEX),
                    1./(ADIABATIC_INDEX-1.));

                if (prim[j][i][RHO] > rhoMax)
                    rhoMax = prim[j][i][RHO];

                PetscRandomGetValue(randNumGen, &randNum);
                prim[j][i][UU] = KAPPA*pow(prim[j][i][RHO], ADIABATIC_INDEX)/\
                              (ADIABATIC_INDEX-1.)*(1. + 4e-2*(randNum-0.5));

                if (prim[j][i][UU] > uMax && r > R_MIN)
                    uMax = prim[j][i][UU];

                uPhiTmp = sqrt((-1. + sqrt(1. + 4.*l*l*(sigma*sigma*delta/\
                               (A*A*sin(theta)*sin(theta))))) / 2.);
			
                uconBL[1] = 0.;
                uconBL[2] = 0.;
                uconBL[3] = 2.*A_SPIN*r*sqrt(1. + uPhiTmp*uPhiTmp)/ \
                            sqrt(A*sigma*delta) + 
                            sqrt(sigma/A)*uPhiTmp/sin(theta);
                
                // Transform uconBoyerLinquist to uconKerrSchild and set to prim

                for (int ii=0; ii<NDIM; ii++)
                    for (int jj=0; jj<NDIM; jj++) {
                        gcovBL[ii][jj] = 0.;
                        gconBL[ii][jj] = 0.;
                        transform[ii][jj] = 0.;
                    }

                DD = 1. - 2./r + A_SPIN*A_SPIN/(r*r);
                mu = 1 + A_SPIN*A_SPIN*cos(theta)*cos(theta)/(r*r);

                gcovBL[0][0] = -(1. - 2./(r*mu));
                gcovBL[0][3] = -2.*A_SPIN*sin(theta)*sin(theta)/(r*mu);
                gcovBL[3][0] = gcovBL[0][3];
                gcovBL[1][1] = mu/DD;
                gcovBL[2][2] = r*r*mu;
                gcovBL[3][3] = r*r*sin(theta)*sin(theta)*\
                               (1. + A_SPIN*A_SPIN/(r*r) + 
                                2.*A_SPIN*A_SPIN*sin(theta)*sin(theta)/\
                                (r*r*r*mu));

                gconBL[0][0] = -1. -2.*(1 + A_SPIN*A_SPIN/(r*r))/(r*DD*mu);
                gconBL[0][3] = -2.*A_SPIN/(r*r*r*DD*mu);
                gconBL[3][0] = gconBL[0][3];
                gconBL[1][1] = DD/mu;
                gconBL[2][2] = 1./(r*r*mu);
                gconBL[3][3] = (1. - 2./(r*mu))/(r*r*sin(theta)*sin(theta)*DD);

                transform[0][0] = 1.;
                transform[1][1] = 1.;
                transform[2][2] = 1.;
                transform[3][3] = 1.;
                transform[0][1] = 2.*r/(r*r* - 2.*r + A_SPIN*A_SPIN);
                transform[3][1] = A_SPIN/(r*r - 2.*r + A_SPIN*A_SPIN);

                AA = gcovBL[0][0];
                BB = 2.*(gcovBL[0][1]*uconBL[1] + 
                         gcovBL[0][2]*uconBL[2] +
                         gcovBL[0][3]*uconBL[3]);
                CC = 1. + gcovBL[1][1]*uconBL[1]*uconBL[1] +
                          gcovBL[2][2]*uconBL[2]*uconBL[2] +
                          gcovBL[3][3]*uconBL[3]*uconBL[3] +
                      2.*(gcovBL[1][2]*uconBL[1]*uconBL[2] +
                          gcovBL[1][3]*uconBL[1]*uconBL[3] +
                          gcovBL[2][3]*uconBL[2]*uconBL[3]);

                discriminent = BB*BB - 4.*AA*CC;
                uconBL[0] = -(BB + sqrt(discriminent))/(2.*AA);

                uconKS[0] = transform[0][0]*uconBL[0] + 
                            transform[0][1]*uconBL[1] +
                            transform[0][2]*uconBL[2] +
                            transform[0][3]*uconBL[3];

                uconKS[1] = transform[1][0]*uconBL[0] + 
                            transform[1][1]*uconBL[1] +
                            transform[1][2]*uconBL[2] +
                            transform[1][3]*uconBL[3];

                uconKS[2] = transform[2][0]*uconBL[0] + 
                            transform[2][1]*uconBL[1] +
                            transform[2][2]*uconBL[2] +
                            transform[2][3]*uconBL[3];

                uconKS[3] = transform[3][0]*uconBL[0] + 
                            transform[3][1]*uconBL[1] +
                            transform[3][2]*uconBL[2] +
                            transform[3][3]*uconBL[3];

                PetscScalar rFactor, hFactor;
                rFactor = r - R0;
                hFactor = M_PI + (1. - H_SLOPE)*M_PI*cos(2.*M_PI*X2);
                uconKSPrime[0] = uconKS[0];
                uconKSPrime[1] = uconKS[1]*(1./rFactor);
                uconKSPrime[2] = uconKS[2]*(1./hFactor);
                uconKSPrime[3] = uconKS[3];

                PetscScalar gcov[NDIM][NDIM], gcon[NDIM][NDIM];
                PetscScalar gdet, alpha;
                gCovCalc(gcov, X1, X2);
                gDetCalc(&gdet, gcov);
                gConCalc(gcon, gcov, gdet);
                alphaCalc(&alpha, gcon);

                prim[j][i][U1] = uconKSPrime[1] + 
                                 alpha*alpha*uconKSPrime[0]*gcon[0][1];
                prim[j][i][U2] = uconKSPrime[2] + 
                                 alpha*alpha*uconKSPrime[0]*gcon[0][2];
                prim[j][i][U3] = uconKSPrime[3] +
                                 alpha*alpha*uconKSPrime[0]*gcon[0][3];

            }

            prim[j][i][B1] = 0.;
            prim[j][i][B2] = 0.;
            prim[j][i][B3] = 0.;

        }

    for (int j=X2Start-2; j<X2Start+X2Size+2; j++) {
        for (int i=X1Start-2; i<X1Start+X1Size+2; i++) {        
            prim[j][i][RHO] = prim[j][i][RHO]/rhoMax;
            prim[j][i][UU] = prim[j][i][UU]/rhoMax;

            AVector[j+NG][i+NG] = 0.;
        }
    }

    uMax = uMax/rhoMax;

    for (int j=X2Start-2; j<X2Start+X2Size+2; j++) {
        for (int i=X1Start-2; i<X1Start+X1Size+2; i++) {        

            REAL X1 = i_TO_X1_CENTER(i);
            REAL X2 = j_TO_X2_CENTER(j);
            REAL gcov[NDIM][NDIM];
            gCovCalc(gcov, X1, X2);

            REAL vars[DOF];
            for (int var=0; var<DOF; var++)
                vars[var] = prim[j][i][var];

            REAL gamma;
            gammaCalc(&gamma, vars, gcov);

            setFloorInit(vars, gamma, X1, X2);

            for (int var=0; var<DOF; var++)
                prim[j][i][var] = vars[var];
        }
    }

    for (int j=X2Start-1; j<X2Start+X2Size+1; j++) {
        for (int i=X1Start-1; i<X1Start+X1Size+1; i++) {        
            rhoAv = 0.25*(prim[j][i][RHO] + prim[j][i-1][RHO] +
                          prim[j-1][i][RHO] + prim[j-1][i-1][RHO]);
            
            q = rhoAv - 0.2;

            if (q > 0.)
                AVector[j+NG][i+NG] = q;

        }
    }

    for (int j=X2Start; j<X2Start+X2Size; j++) {
        for (int i=X1Start; i<X1Start+X1Size; i++) {
            
            X1 = i_TO_X1_CENTER(i);
            X2 = j_TO_X2_CENTER(j);

            PetscScalar gcov[NDIM][NDIM], gcon[NDIM][NDIM];
            PetscScalar g, gdet, alpha;
            gCovCalc(gcov, X1, X2);
            gDetCalc(&gdet, gcov);
            gConCalc(gcon, gcov, gdet);
            alphaCalc(&alpha, gcon);
            g = sqrt(-gdet);

            prim[j][i][B1] = -(AVector[j+NG][i+NG] - AVector[j+NG+1][i+NG] +
                               AVector[j+NG][i+1+NG] - AVector[j+1+NG][i+1+NG])/\
                              (2.*DX2*g);

            prim[j][i][B2] = (AVector[j+NG][i+NG] + AVector[j+1+NG][i+NG] -
                            AVector[j+NG][i+1+NG] - AVector[j+1+NG][i+1+NG])/\
                            (2.*DX1*g);

            prim[j][i][B3] = 0.;

            PetscScalar gamma, ucon[NDIM], ucov[NDIM], bcon[NDIM], bcov[NDIM];
            PetscScalar bsqr, var[DOF];

            for (int n=0; n<DOF; n++)
                var[n] = prim[j][i][n];

            gammaCalc(&gamma, var, gcov);
            uconCalc(ucon, gamma, alpha, var, gcon);
            covFromCon(ucov, ucon, gcov);
            bconCalc(bcon, var, ucon, ucov);
            covFromCon(bcov, bcon, gcov);
            bSqrCalc(&bsqr, bcon, bcov);

            if (bsqr > bsqrMax)
                bsqrMax = bsqr;

        }
    }
    betaActual = (ADIABATIC_INDEX-1.)*uMax/(0.5*bsqrMax);
    norm = sqrt(betaActual/BETA);

    for (int j=X2Start; j<X2Start+X2Size; j++) {
        for (int i=X1Start; i<X1Start+X1Size; i++) {
            prim[j][i][B1] = prim[j][i][B1]*norm;
            prim[j][i][B2] = prim[j][i][B2]*norm;
        }
    }

    for (int j=X2Start-2; j<X2Start+X2Size+2; j++) {
        for (int i=X1Start-2; i<X1Start+X1Size+2; i++) {        

            REAL X1 = i_TO_X1_CENTER(i);
            REAL X2 = j_TO_X2_CENTER(j);
            REAL gcov[NDIM][NDIM];
            gCovCalc(gcov, X1, X2);

            REAL vars[DOF];
            for (int var=0; var<DOF; var++)
                vars[var] = prim[j][i][var];

            REAL gamma;
            gammaCalc(&gamma, vars, gcov);

            setFloorInit(vars, gamma, X1, X2);

            prim[j][i][RHO] = (vars[RHO]);
            prim[j][i][UU] = (vars[UU]);

            for (int var=2; var<DOF; var++)
                prim[j][i][var] = vars[var];

        }
    }

    DMLocalToGlobalBegin(dmda, localPrim, INSERT_VALUES, Prim);
    DMLocalToGlobalEnd(dmda, localPrim, INSERT_VALUES, Prim);

    DMDAVecRestoreArrayDOF(dmda, localPrim, &prim);
    DMRestoreLocalVector(dmda, &localPrim);

    PetscRandomDestroy(&randNumGen);
}

PetscErrorCode Monitor(TS ts, 
                       PetscInt step,
                       PetscReal t,
                       Vec Prim,
                       void *ptr)
{
    DM dmda;
    TSGetDM(ts, &dmda);
    PetscScalar ***prim;
    PetscScalar gcov[NDIM][NDIM];

    DMDAVecGetArrayDOF(dmda, Prim, &prim);

    for (int j=0; j<N2; j++) {
        for (int i=0; i<N1; i++) {        

            REAL X1 = i_TO_X1_CENTER(i);
            REAL X2 = j_TO_X2_CENTER(j);
            gCovCalc(gcov, X1, X2);

            REAL vars[DOF];
            for (int var=0; var<DOF; var++)
                vars[var] = prim[j][i][var];

            REAL gamma;
            gammaCalc(&gamma, vars, gcov);

            setFloorInit(vars, gamma, X1, X2);

            for (int var=0; var<DOF; var++)
                prim[j][i][var] = vars[var];
        }
    }

    DMDAVecRestoreArrayDOF(dmda, Prim, &prim);

    REAL dt, dtDump;
    TSGetTimeStep(ts, &dt);
    
    static PetscInt counter = 0;
    static PetscScalar tDump = 0.;
    dtDump = 1.;

    SNES snes;
    TSGetSNES(ts, &snes);
    Vec Res;
    SNESGetFunction(snes, &Res, NULL, NULL);

    if (t >= tDump) {
        printf("Dumping data..\n");
        char filename[50];
        char errname[50];
        sprintf(filename, "plot%04d.h5", counter);
        sprintf(errname, "residual%04d.h5", counter);

        PetscViewer viewer;
        PetscViewerHDF5Open(PETSC_COMM_WORLD, filename,
                            FILE_MODE_WRITE, &viewer);
        PetscObjectSetName((PetscObject) Prim, "soln");
        VecView(Prim, viewer);
        PetscViewerDestroy(&viewer);

        PetscViewerHDF5Open(PETSC_COMM_WORLD, errname,
                            FILE_MODE_WRITE, &viewer);
        PetscObjectSetName((PetscObject) Res, "Err");
        VecView(Res, viewer);
        PetscViewerDestroy(&viewer);

        tDump = tDump + dtDump;
        counter++;
    }

    REAL time;
    TSGetTime(ts, &time);
//    if (time < 1.)
//        TSSetTimeStep(ts, 0.001);
//    else
        TSSetTimeStep(ts, DT);

    return(0.);
}

PetscErrorCode SNESMonitor(SNES snes, PetscInt its, PetscReal norm, void *ptr)
{
    if (its==0) {
        Vec Res;
        SNESGetFunction(snes, &Res, NULL, NULL);

        char errname[50];
        sprintf(errname, "SNESresidual.h5");

        PetscViewer viewer;
        PetscViewerHDF5Open(PETSC_COMM_WORLD, errname,
                            FILE_MODE_WRITE, &viewer);
        PetscObjectSetName((PetscObject) Res, "Err");
        VecView(Res, viewer);
        PetscViewerDestroy(&viewer);
    }

    return(0);
}
