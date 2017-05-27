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

FEI1dLin DegeneratedBeam3d :: interpolation(1);

DegeneratedBeam3d :: DegeneratedBeam3d(int n, Domain *aDomain) :
    NLStructuralElement(n, aDomain)
{
    //  numberOfDofMans = 2;
    //  nGaussPoints = 8;


    //    nPointsX = 2;
    //    nPointsY = 1;
    //    nPointsZ = 1;

    //    nGaussPoints = nPointsX * nPointsY * nPointsZ;
}

FEInterpolation *
DegeneratedBeam3d :: giveInterpolation() const { return & interpolation; }


FEInterpolation *
DegeneratedBeam3d :: giveInterpolation(DofIDItem id) const
{
    return & interpolation;
}


IRResultType
DegeneratedBeam3d :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro
    /*
    IR_GIVE_OPTIONAL_FIELD(ir, nPointsX, _IFT_DegeneratedBeam3d_nipX); // move to element?
    IR_GIVE_OPTIONAL_FIELD(ir, nPointsY, _IFT_DegeneratedBeam3d_nipY);
    IR_GIVE_OPTIONAL_FIELD(ir, nPointsZ, _IFT_DegeneratedBeam3d_nipZ);
    */
    IR_GIVE_OPTIONAL_FIELD(ir, nlGeometry, _IFT_NLStructuralElement_nlgeoflag);
    

    IR_GIVE_OPTIONAL_FIELD(ir, this->yaxis, _IFT_DegeneratedBeam3d_yaxis);
    if ( ir->hasField(_IFT_DegeneratedBeam3d_zaxis) ) {
	IR_GIVE_FIELD(ir, this->zaxis, _IFT_DegeneratedBeam3d_zaxis);
    }
    // to do : refnode, dofstocondense, ... (see beam3d.C:570)

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
	this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], nGaussPoints, this);
    }
}



void
DegeneratedBeam3d ::  giveDirectorVectors(FloatMatrix &Vs, FloatMatrix &Vt)
{
    Vs.resize(3 ,numberOfDofMans);
    Vt.resize(3 ,numberOfDofMans);
    if ( directorType == 0 ) { // normal to the midplane
        FloatArray e1, e2, e3;
        this->computeLocalBaseVectors(e1, e2, e3);
	for (int i =1; i<=numberOfDofMans; i++) {
	    Vt.setColumn(e2,i);
	    Vs.setColumn(e3,i);
	}
    } else if ( directorType == 1 ) {          // nodal average
	OOFEM_ERROR("DirectorType 1 not implemented yet.");
    } else if ( directorType == 2 ) {          // specified at crosssection
	OOFEM_ERROR("DirectorType 2 not implemented yet.");
    } else {
        OOFEM_ERROR("Unsupported directorType");
    }

    return;
}

void
DegeneratedBeam3d ::  giveLocalDirectorVectors(FloatMatrix &Vs, FloatMatrix &Vt)
{
    FloatMatrix Vsg, Vtg;
    this->giveDirectorVectors(Vsg, Vtg);
    this->computeGtoLRotationMatrix();

    Vs.resize(3, numberOfDofMans);
    Vt.resize(3, numberOfDofMans);

    FloatArray global, local;
    for (int i=1; i<=numberOfDofMans; i++){
	global.beColumnOf(Vsg,i);
	local.beProductOf(GtoLRotationMatrix, global);
	Vs.setColumn(local, i);

	global.beColumnOf(Vtg,i);
	local.beProductOf(GtoLRotationMatrix, global);
	Vt.setColumn(local, i);
    }
}

void
DegeneratedBeam3d :: giveWidth(FloatArray &a)
{
    a.resize(numberOfDofMans);
    FloatArray *c1;
	
    for (int i=1; i<=numberOfDofMans; i++) {
	c1 = this->giveNode(i)->giveCoordinates();
	a.at(i) = this->giveCrossSection()->give(CS_Width, * c1, this, false);
    }
}

void
DegeneratedBeam3d :: giveThickness(FloatArray &b)
{
    b.resize(numberOfDofMans);
    FloatArray *c1;
	
    for (int i=1; i<=numberOfDofMans; i++) {
	c1 = this->giveNode(i)->giveCoordinates();
	b.at(i) = this->giveCrossSection()->give(CS_Thickness, * c1, this, false); // todo
    }
}

void
DegeneratedBeam3d :: computeNmatrixAt(const FloatArray &lCoords, FloatMatrix &answer)
// Returns the [6x24] displacement interpolation matrix {N} of the receiver,
// evaluated at gp.
// Zeroes in rows 4, 5, 6.
{
    answer.resize(3, numberOfDofMans*6);
    FloatArray h(numberOfDofMans);
    interpolation.evalN( h, lCoords,  FEIElementGeometryWrapper(this) );

    FloatArray a, b;
    this->giveWidth(a);
    this->giveThickness(b);

    FloatMatrix Vs, Vt;
    this->giveLocalDirectorVectors(Vs, Vt);

    for (int i=0; i<3; i++){ // number of spatial coordinate
	for (int j=0; j<numberOfDofMans; j++){ 
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
    scs->giveStiffnessMatrix_3d(answer, rMode, gp, tStep);
    //    scs->give3dDegeneratedBeam3dStiffMtrx(answer, rMode, gp, tStep);
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
    FloatArray h(numberOfDofMans);
    FloatMatrix dh(numberOfDofMans,2);
    FloatMatrix dhdx(numberOfDofMans,2);

    jacobianMatrix.resize(3, 3);

    // get local director vector
    FloatMatrix Vs, Vt;
    this->giveLocalDirectorVectors(Vs, Vt);

    // get thickness
    FloatArray a, b;
    this->giveWidth(a);
    this->giveThickness(b);

    interpolation.evalN( h, lcoords,  FEIElementGeometryWrapper(this) );
    interpolation.evaldNdx( dh, lcoords,  FEIElementGeometryWrapper(this) );
    // interpolation.giveDerivatives(dh, lcoords);

    double t = lcoords.at(2);
    double s = lcoords.at(3);
    double dt = 1; // assuming linear interpolation
    double ds = 1;

    // FloatArray x - x-coordinates, z_i = 0, y_i = 0
    FloatArray x = {0, 1};
    for (int j=0; j<numberOfDofMans; j++){
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
    answer.resize(6, 6*numberOfDofMans);
    FloatArray h(numberOfDofMans);
    FloatMatrix jacobianMatrix(3, 3);
    FloatMatrix dh(numberOfDofMans, 1);
    FloatMatrix dhdx(numberOfDofMans, 1);

    jacobianMatrix.resize(3, 3);

    // get local director vector
    FloatMatrix Vs, Vt;
    this->giveLocalDirectorVectors(Vs, Vt);

    // get thickness
    FloatArray a, b;
    this->giveWidth(a);
    this->giveThickness(b);

    FloatArray lcoords = gp -> giveNaturalCoordinates();

    interpolation.evalN( h, lcoords,  FEIElementGeometryWrapper(this) );
    interpolation.evaldNdx( dh, lcoords,  FEIElementGeometryWrapper(this) );
    // interpolation.giveDerivatives(dh, lcoords);

    // get gp coordinates
    // double r = lcoords.at(1);
    double s = lcoords.at(2);
    double t = lcoords.at(3);

    double dt = 1; // assuming linear interpolation
    double ds = 1;

    // Jacobian Matrix
    this->giveJacobian(lcoords, jacobianMatrix);
    FloatMatrix invJ;
    invJ.beInverseOf(jacobianMatrix);

    dhdx = dh;
    dhdx.times(invJ(1,1));


    for (int j=0; j<numberOfDofMans; j++) {
	// du/dx u: V1 = fi2*V3 - fi3*V2 
	// du/dr dr/dx -> inv(1,1)
	answer.at(1, j*6+1) = dh.at(j+1,1) *invJ.at(1,1);
	answer.at(1, j*6+5) = t/2 * a.at(j+1) * dh.at(j+1,1) *invJ.at(1,1) * Vt.at(3,j+1) + s/2 * b.at(j+1) * dh.at(j+1,1) *invJ.at(1,1) * Vs.at(3,j+1);
	answer.at(1, j*6+6) =-t/2 * a.at(j+1) * dh.at(j+1,1) *invJ.at(1,1) * Vt.at(2,j+1) - s/2 * b.at(j+1) * dh.at(j+1,1) *invJ.at(1,1) * Vs.at(2,j+1); 
	// du/dt dt/dx -> inv(1,2)
	answer.at(1, j*6+5) += dt/2 * invJ.at(1,2) *  a.at(j+1) * h.at(j+1) * Vt.at(3,j+1);
	answer.at(1, j*6+6) +=-dt/2 * invJ.at(1,2) *  a.at(j+1) * h.at(j+1) * Vt.at(2,j+1);
	// du/ds ds/dx -> inv(1,3)
	answer.at(1, j*6+5) += ds/2 * invJ.at(1,3) *  b.at(j+1) * h.at(j+1) * Vs.at(3,j+1); 
	answer.at(1, j*6+6) +=-ds/2 * invJ.at(1,3) *  b.at(j+1) * h.at(j+1) * Vs.at(2,j+1); 

	// dv/dy v: V2 = fi3*V1 - fi1*V3
	// dv/dr dr/dy -> inv(2,1) 
	answer.at(2, j*6+2) = dh.at(j+1,1) *invJ.at(2,1);
	answer.at(2, j*6+6) = t/2 * a.at(j+1) * dh.at(j+1,1) *invJ.at(2,1) * Vt.at(1,j+1) + s/2 * b.at(j+1) * dh.at(j+1,1) *invJ.at(2,1) * Vs.at(1,j+1);
	answer.at(2, j*6+4) =-t/2 * a.at(j+1) * dh.at(j+1,1) *invJ.at(2,1) * Vt.at(3,j+1) - s/2 * b.at(j+1) * dh.at(j+1,1) *invJ.at(2,1) * Vs.at(3,j+1);
	// dv/dt dt/dy -> inv(2,2)
	answer.at(2, j*6+6) += dt/2 * invJ.at(2,2) *  a.at(j+1) * h.at(j+1) * Vt.at(1,j+1);
	answer.at(2, j*6+4) +=-dt/2 * invJ.at(2,2) *  a.at(j+1) * h.at(j+1) * Vt.at(3,j+1);
	// dv/ds ds/dy -> inv(2,3)
	answer.at(2, j*6+6) += ds/2 * invJ.at(2,3) *  b.at(j+1) * h.at(j+1) * Vs.at(1,j+1); 
	answer.at(2, j*6+4) +=-ds/2 * invJ.at(2,3) *  b.at(j+1) * h.at(j+1) * Vs.at(3,j+1); 

	// dw/dz w: V3 = fi1*V2 - fi2*V1
	// dw/dr dr/dz -> inv(3,1) 
	answer.at(3, j*6+3) = dh.at(j+1,1) *invJ.at(3,1);
	answer.at(3, j*6+4) = t/2 * a.at(j+1) * dh.at(j+1,1) *invJ.at(3,1) * Vt.at(2,j+1) + s/2 * b.at(j+1) * dh.at(j+1,1) *invJ.at(3,1) * Vs.at(2,j+1);
	answer.at(3, j*6+5) =-t/2 * a.at(j+1) * dh.at(j+1,1) *invJ.at(3,1) * Vt.at(1,j+1) - s/2 * b.at(j+1) * dh.at(j+1,1) *invJ.at(3,1) * Vs.at(1,j+1);
	// dw/dt dt/dz -> inv(3,2)
	answer.at(3, j*6+4) += dt/2 * invJ.at(3,2) *  a.at(j+1) * h.at(j+1) * Vt.at(2,j+1);
	answer.at(3, j*6+5) +=-dt/2 * invJ.at(3,2) *  a.at(j+1) * h.at(j+1) * Vt.at(1,j+1);
	// dw/ds ds/dz -> inv(3,3)
	answer.at(3, j*6+4) += ds/2 * invJ.at(3,3) *  b.at(j+1) * h.at(j+1) * Vs.at(2,j+1); 
	answer.at(3, j*6+5) +=-ds/2 * invJ.at(3,3) *  b.at(j+1) * h.at(j+1) * Vs.at(1,j+1); 


	// dv/dz v: V2 = fi3*V1 - fi1*V3
	// dv/dr dr/dz -> inv(3,1) 
	answer.at(4, j*6+2) = dh.at(j+1,1) *invJ.at(1,1);
	answer.at(4, j*6+6) = t/2 * a.at(j+1) * dh.at(j+1,1) *invJ.at(3,1) * Vt.at(1,j+1) + s/2 * b.at(j+1) * dh.at(j+1,1) *invJ.at(3,1) * Vs.at(1,j+1);
	answer.at(4, j*6+4) =-t/2 * a.at(j+1) * dh.at(j+1,1) *invJ.at(3,1) * Vt.at(3,j+1) - s/2 * b.at(j+1) * dh.at(j+1,1) *invJ.at(3,1) * Vs.at(3,j+1);
	// dv/dt dt/dz -> inv(3,2)
	answer.at(4, j*6+6) += dt/2 * invJ.at(3,2) *  a.at(j+1) * h.at(j+1) * Vt.at(1,j+1);
	answer.at(4, j*6+4) +=-dt/2 * invJ.at(3,2) *  a.at(j+1) * h.at(j+1) * Vt.at(3,j+1);
	// dv/ds ds/dz -> inv(3,3)
	answer.at(4, j*6+6) += ds/2 * invJ.at(3,3) *  b.at(j+1) * h.at(j+1) * Vs.at(1,j+1); 
	answer.at(4, j*6+4) +=-ds/2 * invJ.at(3,3) *  b.at(j+1) * h.at(j+1) * Vs.at(3,j+1); 

	// dw/dy w: V3 = fi1*V2 - fi2*V1
	// dw/dr dr/dy -> inv(2,1) 
	answer.at(4, j*6+3) += dh.at(j+1,1) *invJ.at(2,1);
	answer.at(4, j*6+4) += t/2 * a.at(j+1) * dh.at(j+1,1) *invJ.at(2,1) * Vt.at(2,j+1) + s/2 * b.at(j+1) * dh.at(j+1,1) *invJ.at(2,1) * Vs.at(2,j+1);
	answer.at(4, j*6+5) +=-t/2 * a.at(j+1) * dh.at(j+1,1) *invJ.at(2,1) * Vt.at(1,j+1) - s/2 * b.at(j+1) * dh.at(j+1,1) *invJ.at(2,1) * Vs.at(1,j+1);
	// dw/dt dt/dy -> inv(2,2)
	answer.at(4, j*6+4) += dt/2 * invJ.at(2,2) *  a.at(j+1) * h.at(j+1) * Vt.at(2,j+1);
	answer.at(4, j*6+5) +=-dt/2 * invJ.at(2,2) *  a.at(j+1) * h.at(j+1) * Vt.at(1,j+1);
	// dw/ds ds/dy -> inv(2,3)
	answer.at(4, j*6+4) += ds/2 * invJ.at(2,3) *  b.at(j+1) * h.at(j+1) * Vs.at(2,j+1); 
	answer.at(4, j*6+5) +=-ds/2 * invJ.at(2,3) *  b.at(j+1) * h.at(j+1) * Vs.at(1,j+1); 

    }


    ////////////////////////////////////////// reduced integration


#ifdef  DegeneratedBeam3d_reducedShearIntegration

    FloatArray redCoords;
    redCoords = lcoords;
    redCoords.at(1) = 0.0;
    // Jacobian Matrix
    this->giveJacobian(redCoords, jacobianMatrix);
    invJ.beInverseOf(jacobianMatrix);

    interpolation.evalN( h, redCoords,  FEIElementGeometryWrapper(this) );
    interpolation.evaldNdx( dh, redCoords,  FEIElementGeometryWrapper(this) );

    dhdx = dh;
    dhdx.times(invJ(1,1));
#endif

    for (int j=0; j<numberOfDofMans; j++) {

	// du/dz u: V1 = fi2*V3 - fi3*V2 
	// du/dr dr/dz -> inv(3,1)
	answer.at(5, j*6+1) = dh.at(j+1,1) *invJ.at(3,1);
	answer.at(5, j*6+5) = t/2 * a.at(j+1) * dh.at(j+1,1) *invJ.at(3,1) * Vt.at(3,j+1) + s/2 * b.at(j+1) * dh.at(j+1,1) *invJ.at(3,1) * Vs.at(3,j+1);
	answer.at(5, j*6+6) =-t/2 * a.at(j+1) * dh.at(j+1,1) *invJ.at(3,1) * Vt.at(2,j+1) - s/2 * b.at(j+1) * dh.at(j+1,1) *invJ.at(3,1) * Vs.at(2,j+1); 
	// du/dt dt/dz -> inv(3,2)
	answer.at(5, j*6+5) += dt/2 * invJ.at(3,2) *  a.at(j+1) * h.at(j+1) * Vt.at(3,j+1);
	answer.at(5, j*6+6) +=-dt/2 * invJ.at(3,2) *  a.at(j+1) * h.at(j+1) * Vt.at(2,j+1);
	// du/ds ds/dz -> inv(3,3)
	answer.at(5, j*6+5) += ds/2 * invJ.at(3,3) *  b.at(j+1) * h.at(j+1) * Vs.at(3,j+1); 
	answer.at(5, j*6+6) +=-ds/2 * invJ.at(3,3) *  b.at(j+1) * h.at(j+1) * Vs.at(2,j+1); 

	// dw/dx w: V3 = fi1*V2 - fi2*V1
	// dw/dr dr/dx -> inv(1,1) 
	answer.at(5, j*6+3) += dh.at(j+1,1) *invJ.at(1,1);
	answer.at(5, j*6+4) += t/2 * a.at(j+1) * dh.at(j+1,1) *invJ.at(1,1) * Vt.at(2,j+1) + s/2 * b.at(j+1) * dh.at(j+1,1) *invJ.at(1,1) * Vs.at(2,j+1);
	answer.at(5, j*6+5) +=-t/2 * a.at(j+1) * dh.at(j+1,1) *invJ.at(1,1) * Vt.at(1,j+1) - s/2 * b.at(j+1) * dh.at(j+1,1) *invJ.at(1,1) * Vs.at(1,j+1);
	// dw/dt dt/dx -> inv(1,2)
	answer.at(5, j*6+4) += dt/2 * invJ.at(1,2) *  a.at(j+1) * h.at(j+1) * Vt.at(2,j+1);
	answer.at(5, j*6+5) +=-dt/2 * invJ.at(1,2) *  a.at(j+1) * h.at(j+1) * Vt.at(1,j+1);
	// dw/ds ds/dx -> inv(1,3)
	answer.at(5, j*6+4) += ds/2 * invJ.at(1,3) *  b.at(j+1) * h.at(j+1) * Vs.at(2,j+1); 
	answer.at(5, j*6+5) +=-ds/2 * invJ.at(1,3) *  b.at(j+1) * h.at(j+1) * Vs.at(1,j+1); 


	// du/dy u: V1 = fi2*V3 - fi3*V2 
	// du/dr dr/dy -> inv(2,1)
	answer.at(6, j*6+1) = dh.at(j+1,1) *invJ.at(2,1);
	answer.at(6, j*6+5) = t/2 * a.at(j+1) * dh.at(j+1,1) *invJ.at(2,1) * Vt.at(3,j+1) + s/2 * b.at(j+1) * dh.at(j+1,1) *invJ.at(2,1) * Vs.at(3,j+1);
	answer.at(6, j*6+6) =-t/2 * a.at(j+1) * dh.at(j+1,1) *invJ.at(2,1) * Vt.at(2,j+1) - s/2 * b.at(j+1) * dh.at(j+1,1) *invJ.at(2,1) * Vs.at(2,j+1); 
	// du/dt dt/dy -> inv(2,2)
	answer.at(6, j*6+5) += dt/2 * invJ.at(2,2) *  a.at(j+1) * h.at(j+1) * Vt.at(3,j+1);
	answer.at(6, j*6+6) +=-dt/2 * invJ.at(2,2) *  a.at(j+1) * h.at(j+1) * Vt.at(2,j+1);
	// du/ds ds/dy -> inv(2,3)
	answer.at(6, j*6+5) += ds/2 * invJ.at(2,3) *  b.at(j+1) * h.at(j+1) * Vs.at(3,j+1); 
	answer.at(6, j*6+6) +=-ds/2 * invJ.at(2,3) *  b.at(j+1) * h.at(j+1) * Vs.at(2,j+1); 

	// dv/dx v: V2 = fi3*V1 - fi1*V3
	// dv/dr dr/dx -> inv(1,1) 
	answer.at(6, j*6+2) += dh.at(j+1,1) *invJ.at(1,1);
	answer.at(6, j*6+6) += t/2 * a.at(j+1) * dh.at(j+1,1) *invJ.at(1,1) * Vt.at(1,j+1) + s/2 * b.at(j+1) * dh.at(j+1,1) *invJ.at(1,1) * Vs.at(1,j+1);
	answer.at(6, j*6+4) +=-t/2 * a.at(j+1) * dh.at(j+1,1) *invJ.at(1,1) * Vt.at(3,j+1) - s/2 * b.at(j+1) * dh.at(j+1,1) *invJ.at(1,1) * Vs.at(3,j+1);
	// dv/dt dt/dx -> inv(1,2)
	answer.at(6, j*6+6) += dt/2 * invJ.at(1,2) *  a.at(j+1) * h.at(j+1) * Vt.at(1,j+1);
	answer.at(6, j*6+4) +=-dt/2 * invJ.at(1,2) *  a.at(j+1) * h.at(j+1) * Vt.at(3,j+1);
	// dv/ds ds/dx -> inv(1,3)
	answer.at(6, j*6+6) += ds/2 * invJ.at(1,3) *  b.at(j+1) * h.at(j+1) * Vs.at(1,j+1); 
	answer.at(6, j*6+4) +=-ds/2 * invJ.at(1,3) *  b.at(j+1) * h.at(j+1) * Vs.at(3,j+1); 
    }


}


const FloatMatrix *
DegeneratedBeam3d :: computeGtoLRotationMatrix()
// Returns the rotation matrix of the receiver of the size [3,3]
// coords(local) = T * coords(global)
//
// local coordinate (described by vector triplet e1',e2',e3') is defined as follows:
//
// e1'    : [N4-N3]    Ni - means i - th node
// help   : [N2-N3]
// e3'    : e1' x help
// e2'    : e3' x e1'
{
    if ( !GtoLRotationMatrix.isNotEmpty() ) {
        FloatArray e1, e2, e3;


	// e1, e2, e3 - not neccessarily perpendicular - does it matter?
        this->computeLocalBaseVectors(e1, e2, e3); 

        GtoLRotationMatrix.resize(3, 3);

        for ( int i = 1; i <= 3; i++ ) {
            GtoLRotationMatrix.at(1, i) = e1.at(i);
            GtoLRotationMatrix.at(2, i) = e2.at(i);
            GtoLRotationMatrix.at(3, i) = e3.at(i);
        }
    }

    return & GtoLRotationMatrix;
}


void
DegeneratedBeam3d :: computeLocalBaseVectors(FloatArray &e1, FloatArray &e2, FloatArray &e3)
{
    // exact when linear approximation is used - what about higher order approximation?

    FloatArray help;
    FloatArray coordA, coordB;

    // compute A - (node2+node3)/2
    e1.beDifferenceOf(*this->giveNode(numberOfDofMans)->giveCoordinates(), *this->giveNode(1)->giveCoordinates());
    e1.normalize();

    e3 = this->zaxis;
    e3.normalize();
    
    if (this->yaxis.giveSize() == 0)
	this->yaxis.beVectorProductOf(e3,e1);

    e2 = this->yaxis;
    e2.normalize();
}



    /*
bool
DegeneratedBeam3d :: computeGtoLRotationMatrix(FloatMatrix &answer)
// Returns the rotation matrix of the receiver of the size [24,24]
// r(local) = T * r(global)
// for one node (r written transposed): {u,v,w,alpha,beta} = T * {u,v,w,r1,r2,r3}

{
    this->computeGtoLRotationMatrix();

    answer.resize(24, 24);
    answer.zero();

    FloatMatrix LtoDir1, LtoDir2, LtoDir3, LtoDir4;
    this->computeLToDirectorRotationMatrix(LtoDir1, LtoDir2, LtoDir3, LtoDir4);

    for ( int i = 0; i <= 3; i++ ) {
        answer.setSubMatrix(GtoLRotationMatrix, i * 6 + 1, i * 6 + 1);
    }
    ////
    FloatMatrix help;

    help.beProductOf(LtoDir1, GtoLRotationMatrix);
    answer.setSubMatrix(help, 4, 4);

    help.beProductOf(LtoDir2, GtoLRotationMatrix);
    answer.setSubMatrix(help, 10, 10);

    help.beProductOf(LtoDir3, GtoLRotationMatrix);
    answer.setSubMatrix(help, 16, 16);

    help.beProductOf(LtoDir4, GtoLRotationMatrix);
    answer.setSubMatrix(help, 22, 22);

    return 1;
}
    */

void
DegeneratedBeam3d :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveRealStress_3d(answer, gp, strain, tStep);

    //    this->giveStructuralCrossSection()->giveRealStress_3dDegeneratedBeam(answer, gp, strain, tStep);
}

void
DegeneratedBeam3d :: giveCharacteristicTensor(FloatMatrix &answer, CharTensor type, GaussPoint *gp, TimeStep *tStep)
// returns characteristic tensor of the receiver at given gp and tStep
{
    answer.resize(3, 3);
    answer.zero();
    this->computeGtoLRotationMatrix();
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveStructuralCrossSection()->giveMaterial(gp) );

    if ( ( type == GlobalForceTensor ) ) {
        FloatArray stress, localStress, localStrain;
        this->computeStrainVector(localStrain, gp, tStep);
        this->computeStressVector(localStress, localStrain, gp, tStep);
        mat->transformStressVectorTo(stress,  GtoLRotationMatrix, localStress, false);

        answer.at(1, 1) = stress.at(1);
        answer.at(2, 2) = stress.at(2);
        answer.at(3, 3) = stress.at(3);
        answer.at(1, 2) = stress.at(4);
        answer.at(2, 1) = stress.at(4);
        answer.at(2, 3) = stress.at(5);
        answer.at(3, 2) = stress.at(5);
        answer.at(1, 3) = stress.at(6);
        answer.at(3, 1) = stress.at(6);
    } else if ( ( type == GlobalStrainTensor ) ) {
        FloatArray strain, localStrain;
        this->computeStrainVector(localStrain, gp, tStep);
        mat->transformStrainVectorTo(strain,  GtoLRotationMatrix, localStrain, false);

        answer.at(1, 1) = strain.at(1);
        answer.at(2, 2) = strain.at(2);
        answer.at(3, 3) = strain.at(3);
        answer.at(2, 3) = strain.at(4) / 2.;
        answer.at(3, 2) = strain.at(4) / 2.;
        answer.at(1, 3) = strain.at(5) / 2.;
        answer.at(3, 1) = strain.at(5) / 2.;
        answer.at(1, 2) = strain.at(6) / 2.;
        answer.at(2, 1) = strain.at(6) / 2.;
    } else {
        OOFEM_ERROR("unsupported tensor mode");
        exit(1);
    }
}

void
DegeneratedBeam3d :: printOutputAt(FILE *file, TimeStep *tStep)
// Performs end-of-step operations.
{
    FloatArray v;
    // GaussPoint *gp;

    fprintf(file, "element %d (%8d):\n", this->giveLabel(), number);
    /*
    for ( int i = 0; i < nPointsX; i++ ) {
        for ( int j = 0; j < nPointsY; j++ ) {
        fprintf(file, "  GP %d :", i + j + 1);
	    
	    for ( int k = 0; k < nPointsZ; k++ ) {
		int GPnum = nPointsZ * nPointsY * i + nPointsZ * j + k; // ???
		// gp = integrationRulesArray [ 0 ]->getIntegrationPoint( GPnum );
		gp = integrationRulesArray [ 0 ]->getIntegrationPoint( 0 ); // !!!
		fprintf(file, "\n          GP %d.%d :", i + 1, j + 1);

		this->giveIPValue(v, gp, IST_StrainTensor, tStep);
		fprintf(file, "    strains    ");
		for ( auto &val : v ) {
		    fprintf(file, " %.4e", val);
		}
		
		this->giveIPValue(v, gp, IST_StressTensor, tStep);
		fprintf(file, "\n                      stresses   ");
		for ( auto &val : v ) {
		    fprintf(file, " %.4e", val);
		}
	    }
        }
        fprintf(file, "\n");
	}*/
}

int
DegeneratedBeam3d :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    FloatMatrix globTensor;
    CharTensor cht;

    answer.resize(6);

    if (  type == IST_StrainTensor ) {
        cht = GlobalStrainTensor;

        this->giveCharacteristicTensor(globTensor, cht, gp, tStep);

        answer.at(1) = globTensor.at(1, 1); //xx
        answer.at(2) = globTensor.at(2, 2); //yy
        answer.at(3) = globTensor.at(3, 3); //zz
        answer.at(4) = 2 * globTensor.at(2, 3); //yz
        answer.at(5) = 2 * globTensor.at(1, 3); //xz
        answer.at(6) = 2 * globTensor.at(1, 2); //xy

        return 1;
    } else if ( type == IST_StressTensor ) {
        cht = GlobalForceTensor;

        this->giveCharacteristicTensor(globTensor, cht, gp, tStep);

        answer.at(1) = globTensor.at(1, 1); //xx
        answer.at(2) = globTensor.at(2, 2); //yy
        answer.at(3) = globTensor.at(3, 3); //zz
        answer.at(4) = globTensor.at(2, 3); //yz
        answer.at(5) = globTensor.at(1, 3); //xz
        answer.at(6) = globTensor.at(1, 2); //xy

        return 1;
    } else {
        return NLStructuralElement :: giveIPValue(answer, gp, type, tStep);
    }
}


void
DegeneratedBeam3d :: giveLocalCoordinates(FloatArray &answer, FloatArray &global)
// Returns global x coordinates given in global vector
// transformed into local coordinate system of the
// receiver
{
    FloatArray offset;

    this->computeGtoLRotationMatrix();

    offset = global;
    offset.subtract( * this->giveNode(1)->giveCoordinates() );
    answer.beProductOf(GtoLRotationMatrix, offset);
}




    /*
bool
MITC4Shell :: computeLocalCoordinates(FloatArray &answer, const FloatArray &coords)
//converts global coordinates to local planar area coordinates,
//does not return a coordinate in the thickness direction, but
//does check that the point is in the element thickness
{
    // rotate the input point Coordinate System into the element CS
    FloatArray inputCoords_ElCS;
    std :: vector< FloatArray >lc(3);
    FloatArray llc;
    this->giveLocalCoordinates( inputCoords_ElCS, const_cast< FloatArray & >(coords) );
    for ( int _i = 0; _i < 4; _i++ ) {
        this->giveLocalCoordinates( lc [ _i ], * this->giveNode(_i + 1)->giveCoordinates() );
    }
    bool inplane = interpolation.global2local( llc, inputCoords_ElCS, FEIVertexListGeometryWrapper(lc) ) > 0;
    answer.resize(2);
    answer.at(1) = inputCoords_ElCS.at(1);
    answer.at(2) = inputCoords_ElCS.at(2);
    GaussPoint _gp(NULL, 1, answer, 2.0, _2dPlate);
    // now check if the third local coordinate is within the thickness of element
    bool outofplane = ( fabs( inputCoords_ElCS.at(3) ) <= this->giveCrossSection()->give(CS_Thickness, & _gp) / 2. );

    return inplane && outofplane;
}
    */

int
DegeneratedBeam3d :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    answer.resize(3);
    answer.zero();

    FloatMatrix n;
    this->computeNmatrixAt(lcoords, n);

    for ( int i = 1; i <= 3; i++ ) {
	for (int j = 1; j<=numberOfDofMans; j++ )
	    answer.at(i) += this->giveNode(j)->giveCoordinate(i) * n.at(i,i+(j-1)*6);
    }
    return true;
}

    /*
void
DegeneratedBeam3d :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    // works when linear
    if ( iEdge == 1 ) { // edge between nodes 1 2
        answer = {
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12
        };
    } else {
        OOFEM_ERROR("wrong edge number");
    }
}
    */

}
