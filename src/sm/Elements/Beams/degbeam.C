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
 *               Copyright (C) 1993 - 2017   Borek Patzak
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

#include "../sm/Elements/Beams/degbeam.h"
#include "../sm/Materials/structuralms.h"
#include "../sm/Materials/structuralmaterial.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "fei1dlin.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "../sm/CrossSections/variablecrosssection.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "load.h"
#include "boundaryload.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(DegeneratedBeam3d);

FEI1dLin DegeneratedBeam3d :: interp_lin(1);

DegeneratedBeam3d :: DegeneratedBeam3d(int n, Domain *aDomain) :
    NLStructuralElement(n, aDomain)
{
    nDofMans = 2;
    nPointsX = 2;
    nPointsY = 1;
    nPointsZ = 1;

    nGaussPoints = nPointsX * nPointsY * nPointsZ;
}

FEInterpolation *
DegeneratedBeam3d :: giveInterpolation() const { return & interp_lin; }


FEInterpolation *
DegeneratedBeam3d :: giveInterpolation(DofIDItem id) const
{
    return & interp_lin;
}

void
DegeneratedBeam3d :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_v, D_w, R_u, R_v, R_w
    };
}

IRResultType
DegeneratedBeam3d :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    IR_GIVE_OPTIONAL_FIELD(ir, nPointsX, _IFT_DegeneratedBeam3d_nipX); // move to element?
    IR_GIVE_OPTIONAL_FIELD(ir, nPointsT, _IFT_DegeneratedBeam3d_nipY);
    IR_GIVE_OPTIONAL_FIELD(ir, nPointsZ, _IFT_DegeneratedBeam3d_nipZ);
    IR_GIVE_OPTIONAL_FIELD(ir, nlGeometry, _IFT_NLStructuralElement_nlgeoflag);

    directorType = 0; // default
    IR_GIVE_OPTIONAL_FIELD(ir, directorType, _IFT_DegeneratedBeam3d_directorType); // move to element?

    return this->NLStructuralElement :: initializeFrom(ir);
}

void
DegeneratedBeam3d :: computeGaussPoints()
// Sets up the array containing the eight Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(1);
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 10) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], nPointsX, nPointsY, nPointsZ, this);
    }
}



void
DegeneratedBeam3d ::  giveDirectorVectors(FloatMatrix &Vs, FloatMatrix &Vt)
{
    Vs.resize(3 ,nDofMans);
    Vt.resize(3 ,nDofMans);
    if ( directorType == 0 ) { // normal to the midplane
        FloatArray e1, e2, e3;
        this->computeLocalBaseVectors(e1, e2, e3);
	for (int i =1; i<=nDofMans; i++) {
	    Vt.setColumn(e2,i);
	    Vs.setColumn(e3,i);
	}
    } else if ( directorType == 1 ) {          // nodal average
	OOFEM_ERROR("DirectorType 1 not implemented yet.");
    } else if ( directorType == 2 ) {          // specified at crosssection
	OOFEM_ERROR("DirectorType 1 not implemented yet.");
    } else {
        OOFEM_ERROR("Unsupported directorType");
    }

    return;
}


void
DegeneratedBeam3d :: computeNmatrixAt(const FloatArray &lCoords, FloatMatrix &answer)
// Returns the [6x24] displacement interpolation matrix {N} of the receiver,
// evaluated at gp.
// Zeroes in rows 4, 5, 6.
{
    FloatArray h(nDofMans);
    interp_lin.evalN( h, iLocCoord,  FEIElementGeometryWrapper(this) );

    FloatArray a, b;
    this->giveThickness(a);
    this->giveHeight(b);

    FloatMatrix Vs, Vt;
    this->giveLocalDirectorVectors(Vs, Vt);

    for (int i=0; i<3; i++){ // number of spatial coordinate
	for (int j=0; j<nDofMans; j++){ 
	    answer.at(i+1, j*6+i+1) = h.at(j+1);
	    answer.at(i+1, j*6+4) = 0; // warping
	    answer.at(i+1, j*6+5) = lCoords.at(2)/2 * a.at(j+1) * h.at(j+1) * Vt.at(i+1,j+1); 
	    answer.at(i+1, j*6+6) = lCoords.at(3)/2 * b.at(j+1) * h.at(j+1) * Vs.at(i+1,j+1);
	}
    }
}

void
DegeneratedBeam3d :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    StructuralCrossSection *scs = dynamic_cast< StructuralCrossSection * >( this->giveCrossSection() );
    //SimpleCrossSection *cs = dynamic_cast< SimpleCrossSection * >( this->giveCrossSection() );
    scs->give3dDegeneratedBeam3dStiffMtrx(answer, rMode, gp, tStep);
}



double
DegeneratedBeam3d :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    double detJ, weight;
    FloatMatrix jacobianMatrix(3, 3);

    FloatArray lcoords(3);
    lcoords.at(1) = gp->giveNaturalCoordinate(1);
    lcoords.at(2) = gp->giveNaturalCoordinate(2);
    lcoords.at(3) = gp->giveNaturalCoordinate(3);

    weight = gp->giveWeight();

    this->giveJacobian(lcoords, jacobianMatrix);

    detJ = jacobianMatrix.giveDeterminant();
    return detJ * weight;
}


void
DegeneratedBeam3d :: giveJacobian(FloatArray lcoords, FloatMatrix &jacobianMatrix)
// Returns the jacobianMatrix
{
    FloatArray h(nDofMans);
    FloatMatrix dh(nDofMans,2);
    FloatMatrix dhdx(nDofMans,2);

    jacobianMatrix.resize(3, 3);

    // get local director vector
    FloatMatrix Vs, Vt;
    this->giveLocalDirectorVectors(Vs, Vt);

    // get thickness
    FloatArray a, b;
    this->giveThickness(a);
    this->giveHeight(b);

    interp_lin.evalN( h, lcoords,  FEIElementGeometryWrapper(this) );
    interp_lin.giveDerivatives(dh, lcoords);

    double t = lcoords.at(2);
    double s = lcoords.at(3);
    double dt = 1; // assuming linear interpolation
    double ds = 1;

    // FloatArray x - x-coordinates, z_i = 0, y_i = 0
    for (int j=0; j<nDofMans; j++){
	jacobianMatrix.at(1, 1) += dh.at(j+1,1) * x.at(j+1);
	for (int i=0; i<3; i++){
	    jacobianMatrix.at(1, i+1) += t/2 * a.at(j+1) * dh.at(j+1,1) * Vt.at(i+1,j+1) + s/2 * b.at(j+1) * dh.at(j+1,1) * Vs.at(i+1,j+1); // dx_i/dr
	    jacobianMatrix.at(2, i+1) += dt/2 * a.at(j+1) * h.at(j+1) * Vt.at(i+1,j+1); // dx_i/dt
	    jacobianMatrix.at(3, i+1) += ds/2 * b.at(j+1) * h.at(j+1) * Vs.at(i+1,j+1); // dx_i/ds
	}
    }
    
}



void
DegeneratedBeam3d :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [6x20] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
{
    answer.resize(6, 6*nDofMans);
    FloatArray h(nDofMans);
    FloatMatrix jacobianMatrix(3, 3);
    FloatMatrix dh(nDofMans, 1);
    FloatMatrix dhdx(nDofMans, 1);

    jacobianMatrix.resize(3, 3);

    // get local director vector
    FloatMatrix Vs, Vt;
    this->giveLocalDirectorVectors(Vs, Vt);

    // get thickness
    FloatArray a, b;
    this->giveThickness(a);
    this->giveHeight(b);

    interp_lin.evalN( h, lcoords,  FEIElementGeometryWrapper(this) );
    interp_lin.giveDerivatives(dh, lcoords);

    // get gp coordinates
    double r = gp->giveNaturalCoordinate(1);
    double s = gp->giveNaturalCoordinate(2);
    double t = gp->giveNaturalCoordinate(3);

    double dt = 1; // assuming linear interpolation
    double ds = 1;

    FloatArray lcoords(3);
    lcoords.at(1) = r1;
    lcoords.at(2) = r2;
    lcoords.at(3) = r3;

    // Jacobian Matrix
    this->giveJacobian(lcoords, jacobianMatrix);
    FloatMatrix invJ;
    invJ.beInverseOf(jacobianMatrix);

    dhdx = dh.times(invJ(1,1));


    for (int j=0; j<nDofMans; j++) {
	// du/dx u: V1 = fi2*V3 - fi3*V2 
	// du/dr dr/dx -> inv(1,1)
	answer.at(1, j*6+1) = dh.at(j+1,1) *invJ(1,1);
	answer.at(1, j*6+5) = t/2 * a.at(j+1) * dh.at(1,j+1) *invJ(1,1) * Vt.at(3,j+1) + s/2 * b.at(j+1) * dh.at(1,j+1) *invJ(1,1) * Vs.at(3,j+1);
	answer.at(1, j*6+6) =-t/2 * a.at(j+1) * dh.at(1,j+1) *invJ(1,1) * Vt.at(2,j+1) - s/2 * b.at(j+1) * dh.at(1,j+1) *invJ(1,1) * Vs.at(2,j+1); 
	// du/dt dt/dx -> inv(1,2)
	answer.at(1, j*6+5) += dt/2 * invJ(1,2) *  a.at(j+1) * h.at(j+1) * Vt.at(3,j+1);
	answer.at(1, j*6+6) +=-dt/2 * invJ(1,2) *  a.at(j+1) * h.at(j+1) * Vt.at(2,j+1);
	// du/ds ds/dx -> inv(1,3)
	answer.at(1, j*6+5) += ds/2 * invJ(1,3) *  b.at(j+1) * h.at(j+1) * Vs.at(3,j+1); 
	answer.at(1, j*6+6) +=-ds/2 * invJ(1,3) *  b.at(j+1) * h.at(j+1) * Vs.at(2,j+1); 

	// dv/dy v: V2 = fi3*V1 - fi1*V3
	// dv/dr dr/dy -> inv(2,1) 
	answer.at(2, j*6+2) = dh.at(j+1,1) *invJ(2,1);
	answer.at(2, j*6+6) = t/2 * a.at(j+1) * dh.at(1,j+1) *invJ(2,1) * Vt.at(1,j+1) + s/2 * b.at(j+1) * dh.at(1,j+1) *invJ(2,1) * Vs.at(1,j+1);
	answer.at(2, j*6+4) =-t/2 * a.at(j+1) * dh.at(1,j+1) *invJ(2,1) * Vt.at(3,j+1) - s/2 * b.at(j+1) * dh.at(1,j+1) *invJ(2,1) * Vs.at(3,j+1);
	// dv/dt dt/dy -> inv(2,2)
	answer.at(2, j*6+6) += dt/2 * invJ(2,2) *  a.at(j+1) * h.at(j+1) * Vt.at(1,j+1);
	answer.at(2, j*6+4) +=-dt/2 * invJ(2,2) *  a.at(j+1) * h.at(j+1) * Vt.at(3,j+1);
	// dv/ds ds/dy -> inv(2,3)
	answer.at(2, j*6+6) += ds/2 * invJ(2,3) *  b.at(j+1) * h.at(j+1) * Vs.at(1,j+1); 
	answer.at(2, j*6+4) +=-ds/2 * invJ(2,3) *  b.at(j+1) * h.at(j+1) * Vs.at(3,j+1); 

	// dw/dz w: V3 = fi1*V2 - fi2*V1
	// dw/dr dr/dz -> inv(3,1) 
	answer.at(3, j*6+3) = dh.at(j+1,1) *invJ(3,1);
	answer.at(3, j*6+4) = t/2 * a.at(j+1) * dh.at(1,j+1) *invJ(3,1) * Vt.at(2,j+1) + s/2 * b.at(j+1) * dh.at(1,j+1) *invJ(3,1) * Vs.at(2,j+1);
	answer.at(3, j*6+5) =-t/2 * a.at(j+1) * dh.at(1,j+1) *invJ(3,1) * Vt.at(1,j+1) - s/2 * b.at(j+1) * dh.at(1,j+1) *invJ(3,1) * Vs.at(1,j+1);
	// dw/dt dt/dz -> inv(3,2)
	answer.at(3, j*6+4) += dt/2 * invJ(3,2) *  a.at(j+1) * h.at(j+1) * Vt.at(2,j+1);
	answer.at(3, j*6+5) +=-dt/2 * invJ(3,2) *  a.at(j+1) * h.at(j+1) * Vt.at(1,j+1);
	// dw/ds ds/dz -> inv(3,3)
	answer.at(3, j*6+4) += ds/2 * invJ(3,3) *  b.at(j+1) * h.at(j+1) * Vs.at(2,j+1); 
	answer.at(3, j*6+5) +=-ds/2 * invJ(3,3) *  b.at(j+1) * h.at(j+1) * Vs.at(1,j+1); 


	// dv/dz v: V2 = fi3*V1 - fi1*V3
	// dv/dr dr/dz -> inv(3,1) 
	answer.at(4, j*6+2) = dh.at(j+1,1) *invJ(1,1);
	answer.at(4, j*6+6) = t/2 * a.at(j+1) * dh.at(1,j+1) *invJ(3,1) * Vt.at(1,j+1) + s/2 * b.at(j+1) * dh.at(1,j+1) *invJ(3,1) * Vs.at(1,j+1);
	answer.at(4, j*6+4) =-t/2 * a.at(j+1) * dh.at(1,j+1) *invJ(3,1) * Vt.at(3,j+1) - s/2 * b.at(j+1) * dh.at(1,j+1) *invJ(3,1) * Vs.at(3,j+1);
	// dv/dt dt/dz -> inv(3,2)
	answer.at(4, j*6+6) += dt/2 * invJ(3,2) *  a.at(j+1) * h.at(j+1) * Vt.at(1,j+1);
	answer.at(4, j*6+4) +=-dt/2 * invJ(3,2) *  a.at(j+1) * h.at(j+1) * Vt.at(3,j+1);
	// dv/ds ds/dz -> inv(3,3)
	answer.at(4, j*6+6) += ds/2 * invJ(3,3) *  b.at(j+1) * h.at(j+1) * Vs.at(1,j+1); 
	answer.at(4, j*6+4) +=-ds/2 * invJ(3,3) *  b.at(j+1) * h.at(j+1) * Vs.at(3,j+1); 

	// dw/dy w: V3 = fi1*V2 - fi2*V1
	// dw/dr dr/dy -> inv(2,1) 
	answer.at(4, j*6+3) += dh.at(j+1,1) *invJ(2,1);
	answer.at(4, j*6+4) += t/2 * a.at(j+1) * dh.at(1,j+1) *invJ(2,1) * Vt.at(2,j+1) + s/2 * b.at(j+1) * dh.at(1,j+1) *invJ(2,1) * Vs.at(2,j+1);
	answer.at(4, j*6+5) +=-t/2 * a.at(j+1) * dh.at(1,j+1) *invJ(2,1) * Vt.at(1,j+1) - s/2 * b.at(j+1) * dh.at(1,j+1) *invJ(2,1) * Vs.at(1,j+1);
	// dw/dt dt/dy -> inv(2,2)
	answer.at(4, j*6+4) += dt/2 * invJ(2,2) *  a.at(j+1) * h.at(j+1) * Vt.at(2,j+1);
	answer.at(4, j*6+5) +=-dt/2 * invJ(2,2) *  a.at(j+1) * h.at(j+1) * Vt.at(1,j+1);
	// dw/ds ds/dy -> inv(2,3)
	answer.at(4, j*6+4) += ds/2 * invJ(2,3) *  b.at(j+1) * h.at(j+1) * Vs.at(2,j+1); 
	answer.at(4, j*6+5) +=-ds/2 * invJ(2,3) *  b.at(j+1) * h.at(j+1) * Vs.at(1,j+1); 


	// du/dz u: V1 = fi2*V3 - fi3*V2 
	// du/dr dr/dz -> inv(3,1)
	answer.at(5, j*6+1) = dh.at(j+1,1) *invJ(3,1);
	answer.at(5, j*6+5) = t/2 * a.at(j+1) * dh.at(1,j+1) *invJ(3,1) * Vt.at(3,j+1) + s/2 * b.at(j+1) * dh.at(1,j+1) *invJ(3,1) * Vs.at(3,j+1);
	answer.at(5, j*6+6) =-t/2 * a.at(j+1) * dh.at(1,j+1) *invJ(3,1) * Vt.at(2,j+1) - s/2 * b.at(j+1) * dh.at(1,j+1) *invJ(3,1) * Vs.at(2,j+1); 
	// du/dt dt/dz -> inv(3,2)
	answer.at(5, j*6+5) += dt/2 * invJ(3,2) *  a.at(j+1) * h.at(j+1) * Vt.at(3,j+1);
	answer.at(5, j*6+6) +=-dt/2 * invJ(3,2) *  a.at(j+1) * h.at(j+1) * Vt.at(2,j+1);
	// du/ds ds/dz -> inv(3,3)
	answer.at(5, j*6+5) += ds/2 * invJ(3,3) *  b.at(j+1) * h.at(j+1) * Vs.at(3,j+1); 
	answer.at(5, j*6+6) +=-ds/2 * invJ(3,3) *  b.at(j+1) * h.at(j+1) * Vs.at(2,j+1); 

	// dw/dx w: V3 = fi1*V2 - fi2*V1
	// dw/dr dr/dx -> inv(1,1) 
	answer.at(5, j*6+3) += dh.at(j+1,1) *invJ(1,1);
	answer.at(5, j*6+4) += t/2 * a.at(j+1) * dh.at(1,j+1) *invJ(1,1) * Vt.at(2,j+1) + s/2 * b.at(j+1) * dh.at(1,j+1) *invJ(1,1) * Vs.at(2,j+1);
	answer.at(5, j*6+5) +=-t/2 * a.at(j+1) * dh.at(1,j+1) *invJ(1,1) * Vt.at(1,j+1) - s/2 * b.at(j+1) * dh.at(1,j+1) *invJ(1,1) * Vs.at(1,j+1);
	// dw/dt dt/dx -> inv(1,2)
	answer.at(5, j*6+4) += dt/2 * invJ(1,2) *  a.at(j+1) * h.at(j+1) * Vt.at(2,j+1);
	answer.at(5, j*6+5) +=-dt/2 * invJ(1,2) *  a.at(j+1) * h.at(j+1) * Vt.at(1,j+1);
	// dw/ds ds/dx -> inv(1,3)
	answer.at(5, j*6+4) += ds/2 * invJ(1,3) *  b.at(j+1) * h.at(j+1) * Vs.at(2,j+1); 
	answer.at(5, j*6+5) +=-ds/2 * invJ(1,3) *  b.at(j+1) * h.at(j+1) * Vs.at(1,j+1); 


	// du/dy u: V1 = fi2*V3 - fi3*V2 
	// du/dr dr/dy -> inv(2,1)
	answer.at(6, j*6+1) = dh.at(j+1,1) *invJ(2,1);
	answer.at(6, j*6+5) = t/2 * a.at(j+1) * dh.at(1,j+1) *invJ(2,1) * Vt.at(3,j+1) + s/2 * b.at(j+1) * dh.at(1,j+1) *invJ(2,1) * Vs.at(3,j+1);
	answer.at(6, j*6+6) =-t/2 * a.at(j+1) * dh.at(1,j+1) *invJ(2,1) * Vt.at(2,j+1) - s/2 * b.at(j+1) * dh.at(1,j+1) *invJ(2,1) * Vs.at(2,j+1); 
	// du/dt dt/dy -> inv(2,2)
	answer.at(6, j*6+5) += dt/2 * invJ(2,2) *  a.at(j+1) * h.at(j+1) * Vt.at(3,j+1);
	answer.at(6, j*6+6) +=-dt/2 * invJ(2,2) *  a.at(j+1) * h.at(j+1) * Vt.at(2,j+1);
	// du/ds ds/dy -> inv(2,3)
	answer.at(6, j*6+5) += ds/2 * invJ(2,3) *  b.at(j+1) * h.at(j+1) * Vs.at(3,j+1); 
	answer.at(6, j*6+6) +=-ds/2 * invJ(2,3) *  b.at(j+1) * h.at(j+1) * Vs.at(2,j+1); 

	// dv/dx v: V2 = fi3*V1 - fi1*V3
	// dv/dr dr/dx -> inv(1,1) 
	answer.at(6, j*6+2) += dh.at(j+1,1) *invJ(1,1);
	answer.at(6, j*6+6) += t/2 * a.at(j+1) * dh.at(1,j+1) *invJ(1,1) * Vt.at(1,j+1) + s/2 * b.at(j+1) * dh.at(1,j+1) *invJ(1,1) * Vs.at(1,j+1);
	answer.at(6, j*6+4) +=-t/2 * a.at(j+1) * dh.at(1,j+1) *invJ(1,1) * Vt.at(3,j+1) - s/2 * b.at(j+1) * dh.at(1,j+1) *invJ(1,1) * Vs.at(3,j+1);
	// dv/dt dt/dx -> inv(1,2)
	answer.at(6, j*6+6) += dt/2 * invJ(1,2) *  a.at(j+1) * h.at(j+1) * Vt.at(1,j+1);
	answer.at(6, j*6+4) +=-dt/2 * invJ(1,2) *  a.at(j+1) * h.at(j+1) * Vt.at(3,j+1);
	// dv/ds ds/dx -> inv(1,3)
	answer.at(6, j*6+6) += ds/2 * invJ(1,3) *  b.at(j+1) * h.at(j+1) * Vs.at(1,j+1); 
	answer.at(6, j*6+4) +=-ds/2 * invJ(1,3) *  b.at(j+1) * h.at(j+1) * Vs.at(3,j+1); 



    }



}
