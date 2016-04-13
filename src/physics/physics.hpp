#ifndef GRIM_PHYSICS_H_
#define GRIM_PHYSICS_H_

#include "../params.hpp"
#include "../grid/grid.hpp"
#include "../geometry/geometry.hpp"
#include "../reconstruction/reconstruction.hpp"

inline int DELTA(int const &mu, int const &nu)
{
  return (mu==nu ? 1 : 0);
}

class fluidElement
{
  public:
    array one;

    /* fluidElement parameters */
    array tau, chi_emhd, nu_emhd;
    
    array rho, u, u1, u2, u3, B1, B2, B3;
    array pressure, temperature;
    array qTilde, deltaPTilde;
    array q, deltaP;
  
    array gammaLorentzFactor, uCon[NDIM], uCov[NDIM];
    array bSqr, bCon[NDIM], bCov[NDIM];
    array soundSpeed;
    
    array NUp[NDIM];
    array TUpDown[NDIM][NDIM];

    array gradT[NDIM];
    array dtuCov[NDIM];
    array graduCov[NDIM][NDIM];
    array divuCov;
    array deltaP0;
    array q0;
    array bNorm;

    fluidElement(const grid &prim,
                 const grid &magneticFields,
                 const geometry &geom,
                 int &numReads,
                 int &numWrites
                );
    void set(const grid &prim,
             const grid &magneticFields,
             const geometry &geom,
             int &numReads,
             int &numWrites
            );
    void setFluidElementParameters(const geometry &geom);
    void computeFluidFluxes(const geometry &geom, 
                            const int direction,
                            grid &flux,
                            int &numReads,
                            int &numWrites
                           );                                

    void computeMagneticFluxes(const geometry &geom, 
                               const int dir,
                               grid &flux,
                               int &numReads,
                               int &numWrites
                              );

    void computeMinMaxCharSpeeds(const geometry &geom,
			                           const int dir,
                    			       array &MinSpeed,
                        				 array &MaxSpeed,
                                 int &numReads,
                                 int &numWrites
                                );

    void computeTimeDerivSources(const geometry &geom,
		  		                       const fluidElement &elemOld,
                            		 const fluidElement &elemNew,
                             		 const double dt,
                                 grid &sources,
                                 int &numReads,
                                 int &numWrites
				                        );

    void computeImplicitSources(const geometry &geom,
	  			                      grid &sources,
                                array &tauDamp,
                                int &numReads,
                                int &numWrites
				                       );

    void computeExplicitSources(const geometry &geom,
                      		      grid &sources,
                                int &numReads,
                                int &numWrites
                               );
    void computeEMHDGradients(const geometry &geom,
                              const double dX[3],
                              int &numReads,
                              int &numWrites
                             );                   
};

class riemannSolver
{
  public:
    fluidElement *elemFace;

    grid *fluidFluxLeft, *fluidFluxRight;
    grid *fluidConsLeft, *fluidConsRight;
    grid *magneticFluxLeft, *magneticFluxRight;
    grid *magneticConsLeft, *magneticConsRight;

    array minSpeedLeft,  maxSpeedLeft;
    array minSpeedRight, maxSpeedRight;

    riemannSolver(const grid &prim, 
                  const grid &magneticFields,
                  const geometry &geom
                 );
    ~riemannSolver();

    void HLLFluxFormula(const array &minSpeed,
                        const array &maxSpeed,
                        const grid &fluxLeft,
                        const grid &fluxRight,
                        const grid &consLeft,
                        const grid &consRight,
                        const int dir,
                        grid &flux
                       );

    void LLFFluxFormula(const array &minSpeed,
                        const array &maxSpeed,
                        const grid &fluxLeft,
                        const grid &fluxRight,
                        const grid &consLeft,
                        const grid &consRight,
                        const int dir,
                        grid &flux
                       );

    void solve(const grid &primLeft,
               const grid &primRight,
               const grid &magneticFieldsLeft,
               const grid &magneticFieldsRight,
               const geometry &geomLeft,
               const geometry &geomRight,
               const int dir,
               grid &fluidFluxes,
               grid &magneticFluxes,
               int &numReads,
               int &numWrites
              );
};

#endif /* GRIM_PHYSICS_H_ */
