#ifndef GRIM_GRIDZONE_H_
#define GRIM_GRIDZONE_H_

#include <math.h>
#include "../inputs.h"
#include "../geometry/geometry.h"

#define CENTER           (0)
#define FACE_X1          (1)
#define FACE_X2          (2)
#define CORNER           (3)
#define FACE_X1_PLUS_ONE (4)
#define FACE_X2_PLUS_ONE (5)

/* Entire domain is divided into three parts:
 * 1) Global -- The full domain without ghost zones. Reading the global domain
 *              will give the part of the domain that the current processor can
 *              read. Using this means that the associated array does not need
 *              any communication. Used for arrays associated with global Petsc
 *              Vecs.
 *
 * 2) Local -- The part of the full domain that the current processor can read
 *             with ghost zones. Arrays marked local need communication. Used
 *             for arrays associated with local Petsc Vecs.
 *
 * 3) Tile -- Small chunk of the local/global array that can fit inside the
 *            cache. Used for arrays created on our own.
 */
#if (COMPUTE_DIM==1)

/* When GYROAVERAGING==ON, NUM_FLUXES > DOF, so we just use NUM_FLUXES to allocate
 * static storage */
#define TILE_SIZE ((TILE_SIZE_X1+2*NG)*NUM_FLUXES)
#define INDEX_TILE(zone,var) (var + NUM_FLUXES*((zone)->iInTile + NG) )

#define INDEX_TILE_PLUS_ONE_X1(zone,var) (var + NUM_FLUXES*((zone)->iInTile + 1 + NG) )

#define INDEX_TILE_PLUS_TWO_X1(zone,var) (var + NUM_FLUXES*((zone)->iInTile + 2 + NG ) )

#define INDEX_TILE_MINUS_ONE_X1(zone,var) (var + NUM_FLUXES*((zone)->iInTile - 1 + NG) )

#define INDEX_TILE_MINUS_TWO_X1(zone,var) (var + NUM_FLUXES*((zone)->iInTile - 2 + NG) )

#define INDEX_TILE_MANUAL(iInTile,jInTile,var) (var + NUM_FLUXES*(iInTile + NG) )

#define INDEX_PETSC(ptr,zone,var) (ptr[(zone)->i][var])
#define INDEX_PETSC_MANUAL(ptr,i,zone,var) (ptr[i][var])

#elif (COMPUTE_DIM==2)

#define TILE_SIZE ((TILE_SIZE_X1+2*NG)*(TILE_SIZE_X2+2*NG)*NUM_FLUXES)
#define INDEX_TILE(zone,var) (var + NUM_FLUXES*(\
                                          (zone)->iInTile + NG \
                                         +(TILE_SIZE_X1+2*NG)\
                                         *((zone)->jInTile + NG)\
                                        )\
                             )

#define INDEX_TILE_PLUS_ONE_X1(zone,var) (var + NUM_FLUXES*(\
                                                      (zone)->iInTile + 1 + NG \
                                                     +(TILE_SIZE_X1+2*NG)\
                                                     *((zone)->jInTile + NG)\
                                                    )\
                                         )

#define INDEX_TILE_PLUS_TWO_X1(zone,var) (var + NUM_FLUXES*(\
                                                      (zone)->iInTile + 2 + NG \
                                                     +(TILE_SIZE_X1+2*NG)\
                                                     *((zone)->jInTile + NG)\
                                                    )\
                                         )

#define INDEX_TILE_MINUS_ONE_X1(zone,var) (var + NUM_FLUXES*(\
                                                      (zone)->iInTile - 1 + NG \
                                                      +(TILE_SIZE_X1+2*NG)\
                                                      *((zone)->jInTile + NG)\
                                                     )\
                                          )

#define INDEX_TILE_MINUS_TWO_X1(zone,var) (var + NUM_FLUXES*(\
                                                      (zone)->iInTile - 2 + NG \
                                                      +(TILE_SIZE_X1 + 2*NG)\
                                                      *((zone)->jInTile + NG)\
                                                     )\
                                          )

#define INDEX_TILE_PLUS_ONE_X2(zone,var) (var + NUM_FLUXES*(\
                                                      (zone)->iInTile + NG \
                                                     +(TILE_SIZE_X1 + 2*NG)\
                                                     *((zone)->jInTile + 1 + NG)\
                                                    )\
                                         )

#define INDEX_TILE_PLUS_TWO_X2(zone,var) (var + NUM_FLUXES*(\
                                                      (zone)->iInTile + NG \
                                                     +(TILE_SIZE_X1 + 2*NG)\
                                                     *((zone)->jInTile + 2 + NG)\
                                                    )\
                                         )

#define INDEX_TILE_MINUS_ONE_X2(zone,var) (var + NUM_FLUXES*(\
                                                      (zone)->iInTile + NG \
                                                      +(TILE_SIZE_X1 + 2*NG)\
                                                      *((zone)->jInTile - 1 + NG)\
                                                     )\
                                          )

#define INDEX_TILE_MINUS_TWO_X2(zone,var) (var + NUM_FLUXES*(\
                                                      (zone)->iInTile + NG \
                                                      +(TILE_SIZE_X1 + 2*NG)\
                                                      *((zone)->jInTile - 2 + NG)\
                                                     )\
                                          )

#define INDEX_TILE_MANUAL(iInTile,jInTile,var) \
                                          (var + NUM_FLUXES*(\
                                                      iInTile + NG \
                                                      +(TILE_SIZE_X1 + 2*NG)\
                                                      *(jInTile + NG)\
                                                     )\
                                          )

#define INDEX_PETSC(ptr,zone,var) (ptr[(zone)->j][(zone)->i][var])
#define INDEX_PETSC_MANUAL(ptr,i,j,var) (ptr[j][i][var])

#endif

#if (COMPUTE_DIM==1)

  #define LOOP_OVER_TILES(X1Size, X2Size) \
    for (int iTile=0, jTile=0; iTile<(X1Size)/TILE_SIZE_X1; iTile++)

  #define LOOP_INSIDE_TILE(iStart, iEnd, jStart, jEnd) \
    for (int iInTile=(iStart), jInTile=0; iInTile<(iEnd); iInTile++)

  #define ARRAY(ptr) REAL **ptr

#elif (COMPUTE_DIM==2)

  #define LOOP_OVER_TILES(X1Size, X2Size) \
    for (int jTile=0; jTile<(X2Size)/TILE_SIZE_X2; jTile++) \
      for (int iTile=0; iTile<(X1Size)/TILE_SIZE_X1; iTile++)

  #define LOOP_INSIDE_TILE(iStart, iEnd, jStart, jEnd) \
    for (int jInTile=(jStart); jInTile<(jEnd); jInTile++) \
      for (int iInTile=(iStart); iInTile<(iEnd); iInTile++)

  #define ARRAY(ptr) REAL ***ptr

#endif /* LOOP_INSIDE_TILE for different COMPUTE_DIM */

struct gridZone
{
  int i, j;
  int iTile, jTile; 
  int iInTile, jInTile;
  int X1Size, X2Size;

  REAL dX1, dX2;
};

void getXCoords(const struct gridZone zone[ARRAY_ARGS 1],
                const int location,
                REAL X[ARRAY_ARGS NDIM]);

void setGridZone(const int iTile,
                 const int jTile,
                 const int iInTile,
                 const int jInTile,
                 const int X1Start,
                 const int X2Start,
                 const int X1Size,
                 const int X2Size,
                 struct gridZone zone[ARRAY_ARGS 1]);

#endif /* GRIM_GRIDZONE_H_ */
