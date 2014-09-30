#include "grim.h"

struct data
{
  REAL primBoundaries[(N1 + 2*NG)*(N2 + 2*NG)*DOF];
};

static struct data tsData;
static REAL fluxX1[(N1 + 2*NG)*(N2 + 2*NG)*DOF];
static REAL fluxX2[(N1 + 2*NG)*(N2 + 2*NG)*DOF];
static REAL emf[N2+2*NG][N1+2*NG];
static REAL AVector[N2+2*NG][N1+2*NG];

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

    char buildOptions[2000];

    sprintf(buildOptions,"\
                          -D X1_SIZE=%d \
                          -D X2_SIZE=%d \
                          -D TOTAL_X1_SIZE=%d \
                          -D TOTAL_X2_SIZE=%d \
                          -DOPENCL \
                          -DLEFT_BOUNDARY=%d \
                          -DRIGHT_BOUNDARY=%d \
                          -DBOTTOM_BOUNDARY=%d \
                          -DTOP_BOUNDARY=%d \
                          -DINFLOW_CHECK=%d",
                          X1Size, X2Size, X1Size+2*NG, X2Size+2*NG, 
                          LEFT_BOUNDARY, RIGHT_BOUNDARY,
                          BOTTOM_BOUNDARY, TOP_BOUNDARY, INFLOW_CHECK);

    std::string BuildOptions(buildOptions);

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

    InitialConditionMTITest(ts, soln, &tsData);
//    InitialConditionLinearModes(ts, soln);

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

//    REAL fluxX1[(N1 + 2*NG)*(N2 + 2*NG)*DOF];
//    REAL fluxX2[(N1 + 2*NG)*(N2 + 2*NG)*DOF];

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

//    REAL emf[N2+2*NG][N1+2*NG];

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

void InitialConditionLinearModes(TS ts, Vec Prim)
{
    DM dmda;
    TSGetDM(ts, &dmda);
    PetscScalar ***prim;

    DMDAVecGetArrayDOF(dmda, Prim, &prim);

    FILE *X1Coords;
    X1Coords = fopen("X1Coords.txt", "w");
      
    for (PetscInt i=0; i<N1; i++) 
    {
      fprintf(X1Coords, "%d   %.18f\n", i, i_TO_X1_CENTER(i) );
    } 
    fclose(X1Coords);

    for (PetscInt j=0; j<N2; j++) 
    {
      for (PetscInt i=0; i<N1; i++) 
      {
        REAL X1 = i_TO_X1_CENTER(i);
        REAL X2 = j_TO_X2_CENTER(j);


        /* delta ~ A exp(i k x) and A = a + i b. 
           Therefore delta ~ a cos(k*X1) - b sin(k*X1) */

        /* 1D tests */
//        REAL rho0 = 1.;
//        REAL u0 = 100.;
//        REAL u10 = 0.;
//        REAL u20 = 0.;
//        REAL u30 = 0.;
//        REAL B10 = 1e-5;
//        REAL B20 = 0.;
//        REAL B30 = 0.;
//        REAL phi0 = 0.;
//
//        REAL k = 2*M_PI;
//        REAL amplitude = 1e-3;
        
        /* Eigenvalue = 19.9123 */
//        REAL delta_rho = -0.00234936399163*sin(k*X1);
//        REAL delta_u =  0.002273900078*sin(k*X1);
//        REAL delta_u1 = -0.00744570296682*cos(k*X1);
//        REAL delta_u2 = 0.;
//        REAL delta_u3 = 0.;
//        REAL delta_B1 = 0.;
//        REAL delta_B2 = 0.;
//        REAL delta_B3 = 0.;
//        REAL delta_phi = 0.999966935141*cos(k*X1);

        /* Eigenvalue = -4.99 - 13.12i */
//        REAL delta_rho = -0.00311030349*cos(k*X1) + 0.00118312960098*sin(k*X1);
//        REAL delta_u =  0.00320815217935*cos(k*X1) - 0.0013948219362*sin(k*X1);
//        REAL delta_u1 = -0.00744109969005*cos(k*X1) - 2.71345707727e-06*sin(k*X1);
//        REAL delta_u2 = 0.;
//        REAL delta_u3 = 0.;
//        REAL delta_B1 = 0.;
//        REAL delta_B2 = 0.;
//        REAL delta_B3 = 0.;
//        REAL delta_phi = 0.999960658464*cos(k*X1);
//
//        prim[j][i][RHO] = rho0 + amplitude*delta_rho;
//        prim[j][i][UU] = u0 + amplitude*delta_u;
//        prim[j][i][U1] = u10 + amplitude*delta_u1;
//        prim[j][i][U2] = u20 + amplitude*delta_u2;
//        prim[j][i][U3] = u30 + amplitude*delta_u3;
//        prim[j][i][B1] = B10 + amplitude*delta_B1;
//        prim[j][i][B2] = B20 + amplitude*delta_B2;
//        prim[j][i][B3] = B30 + amplitude*delta_B3;
//        prim[j][i][FF] = phi0 + amplitude*delta_phi;


        /* 2D tests */
        REAL rho0 = 1.;
        REAL u0 = 100.;
        REAL u10 = 0.;
        REAL u20 = 0.;
        REAL u30 = 0.;
        REAL B10 = 0.01;
        REAL B20 = 0.02;
        REAL B30 = 0.;
        REAL phi0 = 0.;

        REAL k1 = 2*M_PI;
        REAL k2 = 2*M_PI;
        REAL amplitude = 1e-3;

        REAL complex delta_rho = 0.00320842119877 - 0.000881561535595*I;
        REAL complex delta_u = -0.00338288038189 + 0.00107837165346*I;
        REAL complex delta_u1 = -0.00332654066543 - 1.60219221159e-06*I;
        REAL complex delta_u2 = -0.00665554256072 - 1.60219221159e-06*I;
        REAL complex delta_u3 = 0.;
        REAL complex delta_B1 = 6.49351798214e-09 - 7.32543863021e-09*I;
        REAL complex delta_B2 = -6.49351798214e-09 + 7.32543863021e-09*I;
        REAL complex delta_B3 = 0.;
        REAL complex delta_phi = 0.999960479215;

        REAL complex mode = cexp(I*(k1*X1 + k2*X2) );

        prim[j][i][RHO] = rho0 + amplitude*creal(delta_rho*mode);
        prim[j][i][UU] = u0 + amplitude*creal(delta_u*mode);
        prim[j][i][U1] = u10 + amplitude*creal(delta_u1*mode);
        prim[j][i][U2] = u20 + amplitude*creal(delta_u2*mode);
        prim[j][i][U3] = u30 + amplitude*creal(delta_u3*mode);
        prim[j][i][B1] = B10 + amplitude*creal(delta_B1*mode);
        prim[j][i][B2] = B20 + amplitude*creal(delta_B2*mode);
        prim[j][i][B3] = B30 + amplitude*creal(delta_B3*mode);
        prim[j][i][FF] = phi0 + amplitude*creal(delta_phi*mode);

      }
    }


    DMDAVecRestoreArrayDOF(dmda, Prim, &prim);
}

void InitialConditionMTITest(TS ts, Vec Prim, struct data *tsData)
{
  DM dmda;

  int X1Start, X2Start;
  int X1Size, X2Size;

  TSGetDM(ts, &dmda);

  gsl_rng *gslRand;
  gslRand = gsl_rng_alloc(gsl_rng_mt19937);     /* use Mersenne twister */
  gsl_rng_set(gslRand, 1.);

  DMDAGetCorners(dmda, 
                 &X1Start, &X2Start, NULL,
                 &X1Size, &X2Size, NULL);

  Vec localPrim;
  REAL ***prim;
  DMGetLocalVector(dmda, &localPrim);
  DMDAVecGetArrayDOF(dmda, localPrim, &prim);

  FILE *bondiSolnRHO, *bondiSolnUU, *bondiSolnU1, *bondiSolnRCoords;
  bondiSolnRHO = fopen("bondi_soln_rho.txt", "r");
  bondiSolnUU = fopen("bondi_soln_u.txt", "r");
  bondiSolnU1 = fopen("bondi_soln_ur.txt", "r");
  bondiSolnRCoords = fopen("bondi_soln_rCoords.txt", "r");

  char *rhoLine = NULL, *uLine = NULL, *rLine = NULL, *urLine = NULL;
  size_t rhoLen=0; ssize_t rhoRead;
  size_t uLen=0; ssize_t uRead;
  size_t urLen=0; ssize_t urRead;
  size_t rLen=0; ssize_t rRead;

  REAL rho[N1+2*NG], uu[N1+2*NG], ur[N1+2*NG], rCoords[N1+2*NG];

  for (int i=X1Start-NG; i<X1Start+X1Size+NG; i++) {
    rhoRead = getline(&rhoLine, &rhoLen, bondiSolnRHO);
    uRead = getline(&uLine, &uLen, bondiSolnUU);
    urRead = getline(&urLine, &urLen, bondiSolnU1);
    rRead = getline(&rLine, &rLen, bondiSolnRCoords);

    rho[i+NG] = atof(rhoLine);
    uu[i+NG] = atof(uLine);
    ur[i+NG] = atof(urLine);
    rCoords[i+NG] = atof(rLine);
  }

  free(rhoLine); free(uLine); free(urLine);
  fclose(bondiSolnRHO);
  fclose(bondiSolnUU);
  fclose(bondiSolnU1);
  fclose(bondiSolnRCoords);

  for (int j=X2Start-NG; j<X2Start+X2Size+NG; j++) {
    for (int i=X1Start-NG; i<X1Start+X1Size+NG; i++) {
  
      REAL X1 = i_TO_X1_CENTER(i);
      REAL X2 = j_TO_X2_CENTER(j);

      REAL r, theta;

      BLCoords(&r, &theta, X1, X2);

      if (abs(r - rCoords[i+NG])>1e-15)
      {
        PetscPrintf(PETSC_COMM_SELF, "r = %f, rCoords = %f, DX1 = %f\n", r,
        rCoords[i+NG], DX1);
        PetscPrintf(PETSC_COMM_SELF, "Mismatch in rCoords! Check r coords in python script\n");
        exit(1);
      }
  
      prim[j][i][RHO] = rho[i+NG];
      prim[j][i][UU] = uu[i+NG]*(1. + 4e-2 * (gsl_rng_uniform(gslRand)- 0.5));


      REAL uConBL[NDIM];
      uConBL[1] = ur[i+NG];
      uConBL[2] = 0.;
      uConBL[3] = 0.;


      /* Initial conditions given by Ben */
      REAL gcov[NDIM][NDIM], gcon[NDIM][NDIM];
      REAL gdet, alpha;
      gCovCalc(gcov, X1, X2);
      gDetCalc(&gdet, gcov);
      gConCalc(gcon, gcov, gdet);
      alphaCalc(&alpha, gcon);

	    REAL a = gcov[1][1];
	    REAL b = gcon[0][1];
	    REAL c = gcon[0][0];
	    REAL v1 = (c*uConBL[1]/r - sqrt(-a*b*b*b*b -
                 a*b*b*c*uConBL[1]*uConBL[1]/(r*r) - b*b*c))/(a*b*b + c);

      prim[j][i][U1] = v1;
      prim[j][i][U2] = 0.;
      prim[j][i][U3] = 0.;

      /* Monopolar magnetic field */
//      REAL qB = 0.001;
//      prim[j][i][B1] = qB/(r*r*r);
//      prim[j][i][B2] = 0.;
//      prim[j][i][B3] = 0.;

      /* Vertical magnetic field */
      AVector[j+NG][i+NG] = 0.000000000001*0.5*r*sin(theta);

#if (CONDUCTION)
      prim[j][i][FF] = 0.;
#endif

    }
  }

  for (int j=X2Start; j<X2Start+X2Size; j++) {
    for (int i=X1Start; i<X1Start+X1Size; i++) {

      REAL X1 = i_TO_X1_CENTER(i);
      REAL X2 = j_TO_X2_CENTER(j);

      REAL gcov[NDIM][NDIM], gcon[NDIM][NDIM];
      REAL gdet, alpha;
      gCovCalc(gcov, X1, X2);
      gDetCalc(&gdet, gcov);
      gConCalc(gcon, gcov, gdet);
      alphaCalc(&alpha, gcon);

      prim[j][i][B1] = -(AVector[j+NG][i+NG] - AVector[j+NG+1][i+NG] +
                         AVector[j+NG][i+1+NG] - AVector[j+1+NG][i+1+NG])/\
                        (2.*DX2*sqrt(-gdet));

      prim[j][i][B2] = (AVector[j+NG][i+NG] + AVector[j+1+NG][i+NG] -
                            AVector[j+NG][i+1+NG] - AVector[j+1+NG][i+1+NG])/\
                            (2.*DX1*sqrt(-gdet));

      prim[j][i][B3] = 0.;
    }
  }

  for (int j=X2Start-NG; j<X2Start+X2Size+NG; j++) {
    for (int i=X1Start-NG; i<X1Start+X1Size+NG; i++) {

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

void InitialConditionAtmosphereTest(TS ts, Vec Prim, struct data *tsData)
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

  FILE *atmosphereSolnRHO, *atmosphereSolnUU, *atmosphereSolnRCoords;
  atmosphereSolnRHO = fopen("atmosphere_soln_rho.txt", "r");
  atmosphereSolnUU = fopen("atmosphere_soln_u.txt", "r");
  atmosphereSolnRCoords = fopen("atmosphere_soln_rCoords.txt", "r");

  char *rhoLine = NULL, *uLine = NULL, *rLine = NULL;
  size_t rhoLen=0; ssize_t rhoRead;
  size_t uLen=0; ssize_t uRead;
  size_t rLen=0; ssize_t rRead;

  REAL rho[N1+2*NG], uu[N1+2*NG], rCoords[N1+2*NG];

  for (int i=X1Start-NG; i<X1Start+X1Size+NG; i++) {
    rhoRead = getline(&rhoLine, &rhoLen, atmosphereSolnRHO);
    uRead = getline(&uLine, &uLen, atmosphereSolnUU);
    rRead = getline(&rLine, &rLen, atmosphereSolnRCoords);

    rho[i+NG] = atof(rhoLine);
    uu[i+NG] = atof(uLine);
    rCoords[i+NG] = atof(rLine);
  }

  free(rhoLine); free(uLine);
  fclose(atmosphereSolnRHO);
  fclose(atmosphereSolnUU);
  fclose(atmosphereSolnRCoords);

  for (int j=X2Start-NG; j<X2Start+X2Size+NG; j++) {
    for (int i=X1Start-NG; i<X1Start+X1Size+NG; i++) {

      REAL X1 = i_TO_X1_CENTER(i);
      REAL X2 = j_TO_X2_CENTER(j);
      REAL r, theta;
      BLCoords(&r, &theta, X1, X2);

      if (abs(r - rCoords[i+NG])>1e-15)
      {
        PetscPrintf(PETSC_COMM_SELF, "r = %f, rCoords = %f, DX1 = %f\n", r,
        rCoords[i+NG], DX1);
        PetscPrintf(PETSC_COMM_SELF, "Mismatch in rCoords! Check r coords in python script\n");
        exit(1);
      }
  
      REAL gcov[NDIM][NDIM], gcon[NDIM][NDIM], alpha, gdet;
      gCovCalc(gcov, X1, X2);
      gDetCalc(&gdet, gcov);
      gConCalc(gcon, gcov, gdet);
      alphaCalc(&alpha, gcon);
    
      prim[j][i][RHO] = rho[i+NG];
      prim[j][i][UU] = uu[i+NG];

      REAL uConBL[NDIM];
      uConBL[0] = 1./sqrt(-gcov[0][0]); uConBL[1] = 0.;
      uConBL[2] = 0.; uConBL[3] = 0.;

	    REAL a = gcov[1][1];
	    REAL b = gcon[0][1];
	    REAL c = gcon[0][0];
	    REAL v1 = (c*uConBL[1]/r - sqrt(-a*b*b*b*b -
                 a*b*b*c*uConBL[1]*uConBL[1]/(r*r) - b*b*c))/(a*b*b + c);

      prim[j][i][U1] = v1;
      prim[j][i][U2] = 0.;
      prim[j][i][U3] = 0.;

      /* Monopolar magnetic field */
      REAL qB = 0.0016;
      prim[j][i][B1] = qB/(r*r*r);
      prim[j][i][B2] = 0.;
      prim[j][i][B3] = 0.;

      //prim[j][i][FF] = 0.;

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


#if (GEOMETRY==MKS)
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
#endif

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
    dtDump = DT_DUMP;

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
