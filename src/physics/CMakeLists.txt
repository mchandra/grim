add_library(physics physics.c physics.h conduction.c integrands.c tetrad.c)
target_link_libraries(physics geometry problem)
