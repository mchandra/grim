#ifndef _GRIM_H
#define _GRIM_H

//#define VIENNACL_BUILD_INFO
//#define VIENNACL_WITH_OPENCL

#include <petsc.h>
#include <petscdmda.h>
//#include <petscviennacl.h>
#include <CL/cl.hpp>
#include <fstream>
#include <ctime>
#include <petscviewerhdf5.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include "constants.h"

static const char help[] = 
    "GRIM -- General Relativistic Implicit Magnetohydrodynamics";

cl_int clErr;
std::vector<cl::Platform> platforms;
std::vector<cl::Device> devices;
cl::Context context;
cl::CommandQueue queue;
cl::Program program;
cl::Kernel kernel;

//viennacl::ocl::program program;

PetscErrorCode CheckCLErrors(cl_int clErr, std::string errMsg)
{
    if (clErr!=CL_SUCCESS) {
        SETERRQ(PETSC_COMM_SELF, 1, errMsg.c_str());
    }
    return(0.);
}
extern void InitialConditionTest(TS ts, Vec prim);
extern void InitialCondition(TS ts, Vec prim);
extern void Benchmark(TS ts, Vec prim);
extern PetscErrorCode ComputeResidual(TS ts,
                               PetscScalar t,
                               Vec Prim, Vec dPrim_dt,
                               Vec F, void *ptr);
extern PetscErrorCode Monitor(TS ts, 
                              PetscInt step,
                              PetscReal time,
                              Vec Prim,
                              void *ptr);

extern PetscErrorCode SNESMonitor(SNES snes,
                                  PetscInt its,
                                  PetscReal norm,
                                  void *ptr);

void transformBLtoMKS(REAL uConBL[NDIM], REAL uConMKS[NDIM], 
                      REAL X1, REAL X2, REAL r, REAL theta);

void InitialConditionMTITest(TS ts, Vec Prim, struct data *tsData);
#endif
