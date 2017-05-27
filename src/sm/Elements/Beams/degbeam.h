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



#define _IFT_DegeneratedBeam3d_Name "degeneratedbeam3d"
// use the definition at beam3d/mitc4 ???
#define _IFT_DegeneratedBeam3d_zaxis "zaxis" 
#define _IFT_DegeneratedBeam3d_yaxis "yaxis" //optional
/*#define _IFT_DegeneratedBeam3d_nipX "nipx"
#define _IFT_DegeneratedBeam3d_nipY "nipy"
#define _IFT_DegeneratedBeam3d_nipZ "nipz"
*/
#define _IFT_DegeneratedBeam3d_directorType "directortype"

namespace oofem {
    class FEI1dLin;

#ifndef __CHARTENSOR
 #define __CHARTENSOR
enum CharTensor {
    LocalStrainTensor,
    GlobalStrainTensor,
    LocalForceTensor,
    GlobalForceTensor,
};
#endif


/// Comment or uncomment the following line to force full or reduced integration
// #define DegeneratedBeam3d_reducedShearIntegration

/**
 * This class implements an quad element based on Mixed Interpolation of Tensorial Components (MITC).
 * This element is a shell element suitable for both thin and thick shells.
 * The element has 24 DOFs (u,v,w-displacements and three rotations) in each node
 *
 * Tasks:
 * - calculating its B,D matrices and dV.
 */

class DegeneratedBeam3d : public NLStructuralElement
{
protected:
    /// Element geometry approximation
    static FEI1dLin interpolation;
    int nPointsX, nPointsY, nPointsZ;
    int directorType, nDofMans, nGaussPoints;
    FloatArray yaxis, zaxis;
    FloatMatrix GtoLRotationMatrix;


public:

    DegeneratedBeam3d(int n, Domain *d);
    virtual ~DegeneratedBeam3d() { }

    virtual FEInterpolation *giveInterpolation() const;
    virtual FEInterpolation *giveInterpolation(DofIDItem id) const;

    //    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);


protected:

    virtual void computeGaussPoints();
    //    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer);
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);
    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);


private:

    virtual integrationDomain giveIntegrationDomain() const { return _Cube; }
    virtual MaterialMode giveMaterialMode() { return _3dMat; }


    void giveLocalDirectorVectors(FloatMatrix &Vs, FloatMatrix &Vt);
    virtual void giveDirectorVectors(FloatMatrix &Vs, FloatMatrix &Vt);
    void giveJacobian(FloatArray lcoords, FloatMatrix &jacobianMatrix);
    void giveLocalCoordinates(FloatArray &answer, FloatArray &global);
    virtual double computeVolumeAround(GaussPoint *gp);
    void giveWidth(FloatArray &a);
    void giveThickness(FloatArray &b);
    void computeLocalBaseVectors(FloatArray &e1, FloatArray &e2, FloatArray &e3);
    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    void giveDofManDofIDMask(int inode, IntArray &answer) const { answer = { D_u, D_v, D_w, R_u, R_v, R_w }; }

    // definition & identification
    virtual const char *giveClassName() const { return "DegeneratedBeam3d"; }
    virtual const char *giveInputRecordName() const { return _IFT_DegeneratedBeam3d_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int computeNumberOfDofs() { return (numberOfDofMans*6); }
    virtual int computeNumberOfGlobalDofs() { return (numberOfDofMans*6); }


    const FloatMatrix * computeGtoLRotationMatrix();
    void giveCharacteristicTensor(FloatMatrix &answer, CharTensor type, GaussPoint *gp, TimeStep *tStep);


    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);

/*    void giveNodeCoordinates(double &x1, double &x2, double &x3, double &x4,
                             double &y1, double &y2, double &y3, double &y4,
                             double &z1, double &z2, double &z3, double &z4);
    void giveLocalDirectorVectors(FloatArray &V1, FloatArray &V2, FloatArray &V3, FloatArray &V4);
    void giveLocalCoordinates(FloatArray &answer, FloatArray &global);
    const FloatMatrix *computeGtoLRotationMatrix();
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    void giveCharacteristicTensor(FloatMatrix &answer, CharTensor type, GaussPoint *gp, TimeStep *tStep);
    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
    virtual bool computeLocalCoordinates(FloatArray &answer, const FloatArray &coords);
    void givedNdx(FloatArray &hkx, FloatArray &hky, FloatArray coords);

    void giveMidplaneIPValue(FloatArray &answer, int gpXY, InternalStateType type, TimeStep *tStep);

*/
  
};
} // end namespace oofem
#endif // degbeam_h
