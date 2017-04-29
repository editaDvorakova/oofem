/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef degbeam_h
#define degbeam_h

#include "../sm/Elements/nlstructuralelement.h"
#include "spatiallocalizer.h"
#include "load.h"


#define _IFT_DegeneratedBeam_Name "degbeam"

namespace oofem {


/**
 * This class implements an quad element based on Mixed Interpolation of Tensorial Components (MITC).
 * This element is a shell element suitable for both thin and thick shells.
 * The element has 24 DOFs (u,v,w-displacements and three rotations) in each node
 *
 * Tasks:
 * - calculating its B,D matrices and dV.
 */

class DegeneratedBeam : public NLStructuralElement
{
protected:
    /// Element geometry approximation
    static FEILin interp_lin;
    int nPointsX, nPointsY, nPointsZ;


public:

    DegeneratedBeam(int n, Domain *d);
    virtual ~MITC4Shell() { }

    virtual FEInterpolation *giveInterpolation() const;
    virtual FEInterpolation *giveInterpolation(DofIDItem id) const;

    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);


protected:

    virtual void computeGaussPoints();
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer);
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);
    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);


private:
/*    void giveNodeCoordinates(double &x1, double &x2, double &x3, double &x4,
                             double &y1, double &y2, double &y3, double &y4,
                             double &z1, double &z2, double &z3, double &z4);
    void giveDirectorVectors(FloatArray &V1, FloatArray &V2, FloatArray &V3, FloatArray &V4);
    void giveLocalDirectorVectors(FloatArray &V1, FloatArray &V2, FloatArray &V3, FloatArray &V4);
    void giveThickness(double &a1, double &a2, double &a3, double &a4);
    void giveJacobian(FloatArray lcoords, FloatMatrix &jacobianMatrix);
    void giveLocalCoordinates(FloatArray &answer, FloatArray &global);
    const FloatMatrix *computeGtoLRotationMatrix();
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    void giveCharacteristicTensor(FloatMatrix &answer, CharTensor type, GaussPoint *gp, TimeStep *tStep);
    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
    virtual bool computeLocalCoordinates(FloatArray &answer, const FloatArray &coords);
    virtual double computeVolumeAround(GaussPoint *gp);
    void computeLocalBaseVectors(FloatArray &e1, FloatArray &e2, FloatArray &e3);
    void givedNdx(FloatArray &hkx, FloatArray &hky, FloatArray coords);

    void giveMidplaneIPValue(FloatArray &answer, int gpXY, InternalStateType type, TimeStep *tStep);

    // definition & identification
    virtual const char *giveClassName() const { return "MITC4Shell"; }
    virtual const char *giveInputRecordName() const { return _IFT_MITC4Shell_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int computeNumberOfDofs() { return 24; }
    virtual int computeNumberOfGlobalDofs() { return 24; }
    virtual integrationDomain giveIntegrationDomain() const { return _3dDegShell; }
    virtual MaterialMode giveMaterialMode() { return _3dDegeneratedShell; }

*/
  
};
} // end namespace oofem
#endif // degbeam_h
