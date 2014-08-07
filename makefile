ALL: grim

PETSC_DIR=/home/manic/petsc_optimized
#PETSC_DIR=/home/mc/Downloads/petsc_debug

CFLAGS = -std=c++0x -lOpenCL -O3

FFLAGS =

CPPFLAGS = -std=c++0x -lOpenCL -O3

FPPFLAGS =

include ${PETSC_DIR}/conf/variables

include ${PETSC_DIR}/conf/rules

grim: grim.o
	-${CLINKER} -o grim grim.o ${PETSC_LIB}
	${RM} grim.o
