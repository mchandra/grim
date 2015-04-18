#include "physics.h"

void setTetrad(const struct geometry *geom,
               struct fluidElement *elem)
{
    /* Make a tetrad with elem.uCon as $\overheadarrow{e}_\hat{0}$ */
    makeTetrad(elem->uCon, elem->bCon, geom, elem->eDownHatUpNoHat, elem->eDownNoHatUpHat);
}

void makeTetrad(const REAL eDown0Hat[NDIM],
                const REAL eDown1Hat[NDIM],
                const struct geometry *geom,
                REAL eDownMuHatUpNu[NDIM][NDIM],
                REAL eDownMuUpNuHat[NDIM][NDIM])
{
    /* eDown0Hat -- $\overheadarrow{e}_\hat{0}$, the zeroth tetrad vector.
     * This is given as a contravarient input in the coordinate basis i.e. 
     * if the input eDown0Hat = {a^0, a^1, a^2, a^3} then,
     *
     * $\overheadarroa{e}_\hat{0} = a^mu $\overheadarrow{e}_mu$ 
     *
     * where $\overheadarrow{e}_{\mu}$ is the coordinate basis.
     *
     * eDownMuHatUpNu -- $e^\nu_\hat{\mu}$, the four tetrad
     * vectors in the coordinate basis. This is the nu'th coordinate component
     * of the mu'th tetrad.
     *
     * $\overheadarrow{e}_\hat{\mu}$ =   eDownMuHatUpNu[mu][nu]
     *                                 * $\overheadarrow{e}_{\nu}$
     *
     * eDownMuUpNuHat -- $e^\hat{\nu}_\mu$, the four coordinate tetrads in the
     * orthonormal basis. This is the nu'th tetrad component of the mu'th
     * coordinate basis.
     *
     * $\overheadarrow{e}_{\mu} =   eDownMuUpNuHat[mu][nu]
     *                            * $\overheadarrow{e}_\hat{\nu}$
     */


    /* First set the given input as the zeroth tetrad in the coordinate basis.
     * This means \hat{\mu} = 0 (zeroth tetrad vector) and we loop over the
     * components of this vector in the coordinate basis. */
#pragma ivdep
    for (int nu=0; nu<NDIM; nu++)
    {
        eDownMuHatUpNu[0][nu] = eDown0Hat[nu];
    }
    normalize(geom, eDownMuHatUpNu[0]);
    /*Done with the zeroth vector (\overhead{e}_\hat{0}) of the tetrad*/

    /*Now onto the first vector of the tetrad (\overhead{e}_\hat{1})*/

    /* Pick some trial vector. Here we let the vector be (0, 1, 0, 0) in the
     * coordinate basis */
    for (int nu=0; nu<NDIM; nu++)
    {
        #if (!GYROAVERAGING)
        eDownMuHatUpNu[1][nu] = DELTA(nu, 1);
        #else
        eDownMuHatUpNu[1][nu] = eDown1Hat[nu];
        #endif
    }
    normalize(geom, eDownMuHatUpNu[1]);
    /* Now make this trial vector orthogonal to (\overhead{e}_\hat{0}) by
     * subtracting out the part of the trial vector parallel to
     * (\overhead{e}_\hat{0}) */
    orthogonalize(eDownMuHatUpNu[0], geom, eDownMuHatUpNu[1]);
    
    /* Now normalize it */
    normalize(geom, eDownMuHatUpNu[1]);
    /* Done with the first vector (\overhead{e}_\hat{1}) of the tetrad*/

    /*Now onto the second vector of the tetrad (\overhead{e}_\hat{2})*/

    /* Pick some trial vector. Here we let the vector be (0, 0, 1, 0) in the
     * coordinate basis */
    for (int nu=0; nu<NDIM; nu++)
    {
        eDownMuHatUpNu[2][nu] = DELTA(nu, 2);
    }
    /* Now make this trial vector orthogonal to (\overhead{e}_\hat{0}) and
     * (\overhead{e}_\hat{1}) by subtracting out the part of the trial vector
     * parallel to (\overhead{e}_\hat{0}) and (\overhead{e}_\hat{1})  */
    orthogonalize(eDownMuHatUpNu[0], geom, eDownMuHatUpNu[2]);
    orthogonalize(eDownMuHatUpNu[1], geom, eDownMuHatUpNu[2]);
    
    /* Now normalize it */
    normalize(geom, eDownMuHatUpNu[2]);
    /* Done with the second vector (\overhead{e}_\hat{2}) of the tetrad*/

    /*Now onto the third vector of the tetrad (\overhead{e}_\hat{3})*/

    /* Pick some trial vector. Here we let the vector be (0, 0, 0, 1) in the
     * coordinate basis */
    for (int nu=0; nu<NDIM; nu++)
    {
        eDownMuHatUpNu[3][nu] = DELTA(nu, 3);
    }
    /* Now make this trial vector orthogonal to (\overhead{e}_\hat{0}),
     * (\overhead{e}_\hat{1}) and (\overhead{e}_\hat{2}) by subtracting out the
     * part of the trial vector parallel to (\overhead{e}_\hat{0}),
     * (\overhead{e}_\hat{1}) and (\overhead{e}_\hat{2}) */
    orthogonalize(eDownMuHatUpNu[0], geom, eDownMuHatUpNu[3]);
    orthogonalize(eDownMuHatUpNu[1], geom, eDownMuHatUpNu[3]);
    orthogonalize(eDownMuHatUpNu[2], geom, eDownMuHatUpNu[3]);
    
    /* Now normalize it */
    normalize(geom, eDownMuHatUpNu[3]);
    /* Done with the third vector (\overhead{e}_\hat{3}) of the tetrad*/

    /* Now fill in the matrix to transform a vector written down in coordinate
     * tetrad to an orthonormal tetrad */

    /* First we lower the coordinate index of eDownMuHatUpNu. Use the metric
     * tensor in the coordinate tetrad */

    /* Loop over tetrad index */
    for (int nu=0; nu<NDIM; nu++)
    {
        /* Loop over coordinate index */
        for (int mu=0; mu<NDIM; mu++)
        {
            eDownMuUpNuHat[mu][nu] = 0.;
            
            for (int alpha=0; alpha<NDIM; alpha++)
            {
                /* Lower the coordinate index in eDownMuHatUpNu */
                eDownMuUpNuHat[mu][nu] +=   eDownMuHatUpNu[nu][alpha]
                                          * geom->gCov[alpha][mu];

            }
        }
    }
    /* Now we have both orthonormal tetrad index and coordinate index down. Need
     * to raise the tetrad index. Use the metric tensor in the orthonormal
     * tetrad. Same as Minkowski! */

    /* Loop over coordinate index */
    for (int mu=0; mu<NDIM; mu++)
    {
        /* Only need to change the 0 component since the metric tensor in
         * orthonormal tetrad is Minkowski */
        eDownMuUpNuHat[mu][0] = -eDownMuUpNuHat[mu][0];
    }

    /* All done! */
}

void normalize(const struct geometry *geom, 
               REAL vecCon[NDIM])
{
    REAL vecCov[NDIM];
    conToCov(vecCon, geom, vecCov);
    
    REAL norm = covDotCon(vecCov, vecCon);

    /* Take the sqrt of the abs value cause it could be either a timelike or a
     * spacelike vector */
    norm = sqrt(fabs(norm));

    for (int mu=0; mu<NDIM; mu++)
    {
        vecCon[mu] /= norm;
    }
}

void orthogonalize(const REAL inputVecCon[NDIM],
                   const struct geometry *geom,
                   REAL trialVecCon[NDIM])
{
    /* Gram-Schmidt orthogonalization. 
     *
     * trialVec =   trialVec 
     *            - <inputVec, trialVec>/<inputVec, inputVec> * inputVec
     */
    REAL inputVecCov[NDIM];
    conToCov(inputVecCon, geom, inputVecCov);

    REAL inputVecDotInputVec = covDotCon(inputVecCon, inputVecCov);
    REAL trialVecDotInputVec = covDotCon(trialVecCon, inputVecCov);

#pragma ivdep
    for (int mu=0; mu<NDIM; mu++)
    {
        trialVecCon[mu] -=   trialVecDotInputVec * inputVecCon[mu] \
                           / inputVecDotInputVec;
    }
}
