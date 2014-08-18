ALL: grim

PETSC_DIR=/home/mc/Downloads/petsc_optimized
#PETSC_DIR=/home/mc/Downloads/petsc_debug

CFLAGS = -std=c++0x -lOpenCL -O3 -lgsl -lgslcblas

FFLAGS =

CPPFLAGS = -std=c++0x -lOpenCL -O3 -lgsl -lgslcblas

FPPFLAGS =

include ${PETSC_DIR}/conf/variables

include ${PETSC_DIR}/conf/rules

grim: grim.o
	-${CLINKER} -o grim grim.o ${PETSC_LIB}
	${RM} grim.o
