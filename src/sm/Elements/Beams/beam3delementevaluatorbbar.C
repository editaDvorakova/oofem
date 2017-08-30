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

#include "../sm/Elements/Beams/beam3delementevaluatorbbar.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "domain.h"
#include "node.h"
#include "load.h"
#include "boundaryload.h"
#include "element.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "matresponsemode.h"
#include "crosssection.h"
#include "mathfem.h"
#include "iga/iga.h"
#include "iga/feibspline.h"

namespace oofem {
    /*void Beam3dElementEvaluatorBbar :: computeBMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
    FloatArray ders, patchDers;  // row
    FloatMatrix help;

    Element *elem = this->giveElement();
    // IntegrationRule *iRule = gp->giveIntegrationRule();

    IntArray irlocnum;
    IntArray mask, nodeDofIDMask;
    int dofmandof;

    // get number of dofs in node
    elem->giveDofManDofIDMask(1, nodeDofIDMask);
    dofmandof = nodeDofIDMask.giveSize();

    elem->giveInterpolation()->giveKnotSpanBasisFuncMask(*gp->giveIntegrationRule()->giveKnotSpan(), mask);
    irlocnum.resize( dofmandof*mask.giveSize() );
    for ( int i = 1; i <= mask.giveSize(); i++ ) {
	for ( int j = 1; j <= dofmandof; j++ ) {
	    irlocnum.at((i-1)*dofmandof + j) = dofmandof * ( mask.at(i) - 1 ) + j;
	}
    }
    // this->giveIntegrationElementLocalCodeNumbers(irlocnum, elem, iRule);

    BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (gp->giveElement()->giveInterpolation());

    //  IGAElement *igaelem = dynamic_cast<IGAElement*> (gp->giveElement());
    //  BSplineInterpolation *bBarInterp = dynamic_cast<BSplineInterpolation*> (igaelem->giveBbarInterpolation());
  
    //  const double *const* kv = bBarInterp->giveKnotVector();
    // kv = interp->giveKnotVector();

    FloatMatrix d; // check : derivatives according to ksi?!?!?!
    interp->evaldNdx( d, gp->giveNaturalCoordinates(),
	FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

    // double J;
    // J = interp->giveTransformationJacobian( gp->giveNaturalCoordinates(),
    //		FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

    return;
}*/


const FloatMatrix *
Beam3dElementEvaluatorBbar :: computeBbarMassMatrix()
{
    if ( !invM.isNotEmpty() ) {
	int numberOfIntegrationRules;
	FloatMatrix temp, bj, d, dbj;
	IGAElement *elem = dynamic_cast<IGAElement*>(this->giveElement());
	BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (this->giveElement()->giveInterpolation());
	BSplineInterpolation *bBarInterp = dynamic_cast<BSplineInterpolation*> (elem->giveBbarInterpolation());
	int mSize = interp->giveNumberOfControlPoints(1) - 1; //lower order basis
	IntArray mask;
	FloatMatrix M, locM;
	FloatArray n;
	const IntArray* knotSpan;
	IntArray newKnotSpan[1];

  
	M.resize(mSize, mSize);
	M.zero();

	numberOfIntegrationRules = elem->giveNumberOfIntegrationRules();
	// loop over individual integration rules
	for ( int ir = 0; ir < numberOfIntegrationRules; ir++ ) {
#ifdef __PARALLEL_MODE
	    if ( this->giveElement()->giveKnotSpanParallelMode(ir) == Element_remote ) {
		continue;
	    }
	    //fprintf (stderr, "[%d] Computing element.knotspan %d.%d\n", elem->giveDomain()->giveEngngModel()->giveRank(), elem->giveNumber(), ir);
#endif
	    int nBF = mask.giveSize();
	    locM.clear();
	    locM.resize(nBF,nBF);

	    IntegrationRule *iRule;
	    iRule = elem->giveIntegrationRule(ir);

	    // change knotspan number for lower approximation
	    knotSpan = iRule->giveKnotSpan();
	    newKnotSpan->resize(1);
	    newKnotSpan->at(1) = knotSpan->at(1)-1;

	    bBarInterp->giveKnotSpanBasisFuncMask( *newKnotSpan, mask);
	    // loop over individual integration points
	    for ( GaussPoint *gp: *iRule ) {
		double dV = this->computeVolumeAround(gp);
		bBarInterp->evalN( n, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), newKnotSpan ) );
		temp.beDyadicProductOf(n, n);
		locM.add(dV, temp);
	    }
	    M.assemble(locM, mask, mask);
	}	    
	invM.beInverseOf(M);
	
    }
    return &invM;
}


void Beam3dElementEvaluatorBbar :: computeBendingPartOfBMatrixAt(FloatMatrix &Bb, GaussPoint *gp)
{
    double J, kappa, tau;
    FloatMatrix d;
    FloatArray n;
    
    BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (gp->giveElement()->giveInterpolation());
      
    interp->evaldNdx( d, gp->giveNaturalCoordinates(),
	    FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

    interp->evalN( n, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

    J = interp->giveTransformationJacobian( gp->giveNaturalCoordinates(),
	    FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
      
    kappa = this->giveCurvature( gp->giveNaturalCoordinates(),
			      FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

    tau = this->giveTorsion( gp->giveNaturalCoordinates(),
			      FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

    Bb.resize(6, 6*d.giveNumberOfRows());
    for ( int i = 1; i <= d.giveNumberOfRows(); i++ ) {

	Bb.at(4, i*6 - 2) = d.at(i, 1) / J;
	Bb.at(4, i*6 - 1) = -kappa * n.at(i);
      
	Bb.at(5, i*6 - 2) = kappa * n.at(i);
	Bb.at(5, i*6 - 1) = d.at(i, 1) / J;
	Bb.at(5, i*6) = -tau * n.at(i);

	Bb.at(6, i*6 - 1) = tau * n.at(i);
	Bb.at(6, i*6) = d.at(i, 1) / J;
    }
    return;
}
    

void Beam3dElementEvaluatorBbar :: computeKBbarMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    //KBbar
    FloatMatrix Bm, Bsn, Bsb;
    // compute invM
    this->computeBbarMassMatrix();
    
    
    int numberOfIntegrationRules;
    IGAElement *elem = dynamic_cast<IGAElement*>(this->giveElement());
    BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (this->giveElement()->giveInterpolation());
    BSplineInterpolation *bBarInterp = dynamic_cast<BSplineInterpolation*> (elem->giveBbarInterpolation());

    FloatMatrix mat;
    FloatMatrix Km, Ksn, Ksb;
    double J, iJ, coeff, kappa, tau;
    FloatMatrix Bm1, Bm2, Bsn1, Bsn2, Bsn3, Bsn6, Bsb2, Bsb3, Bsb5;
    FloatMatrix locBm1, locBm2, locBsn1, locBsn2, locBsn3, locBsn6, locBsb2, locBsb3, locBsb5;
    FloatMatrix D, bBarD;
    FloatMatrix nn, nd, temp, help;
    FloatArray N, bBarN, D1;
    IntArray mask, bBarMask;
    int nCP = interp->giveNumberOfControlPoints(1);
    int nCPbBar = nCP - 1; //lower order basis
    
    Bm.resize(nCPbBar, nCP*6);
    Bm1.resize(nCPbBar, nCP);
    Bm2.resize(nCPbBar, nCP);
    
    Bsn.resize(nCPbBar, nCP*6);
    Bsn1.resize(nCPbBar, nCP);
    Bsn2.resize(nCPbBar, nCP);
    Bsn3.resize(nCPbBar, nCP);
    Bsn6.resize(nCPbBar, nCP);

    Bsb.resize(nCPbBar, nCP*6);
    Bsb2.resize(nCPbBar, nCP);
    Bsb3.resize(nCPbBar, nCP);
    Bsb5.resize(nCPbBar, nCP);

    numberOfIntegrationRules = elem->giveNumberOfIntegrationRules();
    // loop over individual integration rules
    for ( int ir = 0; ir < numberOfIntegrationRules; ir++ ) {
#ifdef __PARALLEL_MODE
	if ( this->giveElement()->giveKnotSpanParallelMode(ir) == Element_remote ) {
	    continue;
	}
	//fprintf (stderr, "[%d] Computing element.knotspan %d.%d\n", elem->giveDomain()->giveEngngModel()->giveRank(), elem->giveNumber(), ir);
#endif
	
	IntegrationRule *iRule;
	iRule = elem->giveIntegrationRule(ir);
	
	const IntArray* knotSpan;
	IntArray newKnotSpan[1];
	// change knotspan number for lower approximation
	knotSpan = iRule->giveKnotSpan();
	newKnotSpan->resize(1);
	newKnotSpan->at(1) = knotSpan->at(1)-1;
	
	bBarInterp->giveKnotSpanBasisFuncMask( *newKnotSpan, bBarMask);
	interp->giveKnotSpanBasisFuncMask( *knotSpan, mask);

	int nBF = mask.giveSize();
	int nBFbBar = bBarMask.giveSize();
	locBm1.clear();
	locBm1.resize(nBFbBar,nBF);
	locBm2.clear();
	locBm2.resize(nBFbBar,nBF);
	locBsn1.clear();
	locBsn1.resize(nBFbBar,nBF);
	locBsn2.clear();
	locBsn2.resize(nBFbBar,nBF);
	locBsn3.clear();
	locBsn3.resize(nBFbBar,nBF);
	locBsn6.clear();
	locBsn6.resize(nBFbBar,nBF);
	locBsb2.clear();
	locBsb2.resize(nBFbBar,nBF);
	locBsb3.clear();
	locBsb3.resize(nBFbBar,nBF);
	locBsb5.clear();
	locBsb5.resize(nBFbBar,nBF);
	// loop over individual integration points
	for ( GaussPoint *gp: *iRule ) {
	double dV = this->computeVolumeAround(gp);
	    bBarInterp->evalN( bBarN, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), newKnotSpan ) );
	interp->evalN( N, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), knotSpan ) );
	    
	interp->evaldNdx( D, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), knotSpan ) );
	   
	J = interp->giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
	iJ = 1./J;

	kappa = this->giveCurvature( gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

	tau = this->giveTorsion( gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

	this->computeConstitutiveMatrixAt(mat, rMode, gp, tStep);

	D1.beColumnOf(D, 1);
	D1.times(iJ);
	    
	nn.beDyadicProductOf(bBarN, N);
	nd.beDyadicProductOf(bBarN, D1);

	coeff = dV * sqrt(mat.at(1,1));
	locBm1.add(coeff, nd);
	coeff = -kappa*dV * sqrt(mat.at(1,1));
	locBm2.add(coeff, nn);
	coeff = kappa*dV * sqrt(mat.at(2,2));
	locBsn1.add(coeff, nn);
	coeff = dV * sqrt(mat.at(2,2));
	locBsn2.add(coeff, nd);
	coeff = -tau*dV * sqrt(mat.at(2,2));
	locBsn3.add(coeff, nn);
	coeff = -dV * sqrt(mat.at(2,2));
	locBsn6.add(coeff, nn);
	coeff = tau*dV * sqrt(mat.at(3,3));
	locBsb2.add(coeff, nn);
	coeff = dV * sqrt(mat.at(3,3));
	locBsb3.add(coeff, nd);
	coeff = dV * sqrt(mat.at(3,3));
	locBsb5.add(coeff, nn);

    }
	Bm1.assemble(locBm1, bBarMask, mask);
	Bm2.assemble(locBm2, bBarMask, mask);
	Bsn1.assemble(locBsn1, bBarMask, mask);
	Bsn2.assemble(locBsn2, bBarMask, mask);
        Bsn3.assemble(locBsn3,bBarMask, mask);
	Bsn6.assemble(locBsn6, bBarMask, mask);
	Bsb2.assemble(locBsb2, bBarMask, mask);
	Bsb3.assemble(locBsb3, bBarMask, mask);
	Bsb5.assemble(locBsb5, bBarMask, mask);
    }

IntArray mask1;
    IntArray mask2;
    IntArray mask3;
    IntArray mask4;
    IntArray mask5;
    IntArray mask6;
    mask1.resize(nCP);
    mask2.resize(nCP);
    mask3.resize(nCP);
    mask4.resize(nCP);
    mask5.resize(nCP);
    mask6.resize(nCP);
    for (int i = 0; i < nCP; i++ ) {
      mask1(i) = 6 * i + 1;
      mask2(i) = 6 * i + 2;
      mask3(i) = 6 * i + 3;
      mask4(i) = 6 * i + 4;
      mask5(i) = 6 * i + 5;
      mask6(i) = 6 * i + 6;
    }
    IntArray maskRow;
    maskRow.enumerate(nCP-1);

    Bm.assemble(Bm1, maskRow, mask1);
    Bm.assemble(Bm2, maskRow, mask2);

    Bsn.assemble(Bsn1, maskRow, mask1);
    Bsn.assemble(Bsn2, maskRow, mask2);
    Bsn.assemble(Bsn3, maskRow, mask3);
    Bsn.assemble(Bsn6, maskRow, mask6);

    Bsb.assemble(Bsb2, maskRow, mask2);
    Bsb.assemble(Bsb3, maskRow, mask3);
    Bsb.assemble(Bsb5, maskRow, mask5);

	help.beTProductOf(Bm, invM);
	Km.beProductOf(help, Bm);

	help.beTProductOf(Bsn, invM);
	Ksn.beProductOf(help, Bsn);

	help.beTProductOf(Bsb, invM);
	Ksb.beProductOf(help, Bsb);

	answer = Km;
	answer.add(Ksn);
	answer.add(Ksb);

    }
 


	/*

void Beam3dElementEvaluatorBbar :: computeKMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{

    IGAElement *elem = dynamic_cast<IGAElement*>(this->giveElement());
    int ndofs = elem->computeNumberOfDofs();
    answer.resize(ndofs, ndofs);
    answer.zero();

    FloatMatrix KBbar, Kb;
    // bBar
    this->computeKBbarMatrix(KBbar, rMode, tStep);
    // bending
    this->computeKbMatrix(Kb, rMode, tStep);
    
    // element stiffness matrix
    answer.add(KBbar);
    answer.add(Kb);
}
	*/

void Beam3dElementEvaluatorBbar :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    FloatMatrix KBbar;
    this->computeKBbarMatrix(KBbar, rMode, tStep);

    int numberOfIntegrationRules;
    FloatMatrix temp, bj, d, dbj;
    IGAElement *elem = dynamic_cast<IGAElement*>(this->giveElement());
    StructuralCrossSection *cs = static_cast< StructuralCrossSection * >( elem->giveCrossSection() );
    int ndofs = elem->computeNumberOfDofs();
    bool matStiffSymmFlag = cs->isCharacteristicMtrxSymmetric(rMode);
    IntArray irlocnum;

    answer.resize(ndofs, ndofs);
    answer.zero();

    FloatMatrix m2;
    m2.resize(ndofs, ndofs);
    m2.zero();

    FloatMatrix *m = & answer;
    if ( elem->giveInterpolation()->hasSubPatchFormulation() ) {
        m = & temp;
    }

  

  
    numberOfIntegrationRules = elem->giveNumberOfIntegrationRules();
    // loop over individual integration rules
    for ( int ir = 0; ir < numberOfIntegrationRules; ir++ ) {
#ifdef __PARALLEL_MODE
        if ( this->giveElement()->giveKnotSpanParallelMode(ir) == Element_remote ) {
            continue;
        }
        //fprintf (stderr, "[%d] Computing element.knotspan %d.%d\n", elem->giveDomain()->giveEngngModel()->giveRank(), elem->giveNumber(), ir);
#endif
        m->clear();
        IntegrationRule *iRule;
	iRule = elem->giveIntegrationRule(ir);
	// loop over individual integration points
	for ( GaussPoint *gp: *iRule ) {
	    double dV = this->computeVolumeAround(gp);
	    //this->computeNmMatrixAt(nm, gp);
	    if ( matStiffSymmFlag ) {
		m->plusProductSymmUpper(bj, bj, dV);
	    } else {
		m->plusProductUnsym(bj, bj, dV);
	    }
	    this->computeBendingPartOfBMatrixAt(bj, gp);
	    this->computeConstitutiveMatrixAt(d, rMode, gp, tStep);
	    dbj.beProductOf(d, bj);
	    if ( matStiffSymmFlag ) {
		m->plusProductSymmUpper(bj, dbj, dV);
	    } else {
		m->plusProductUnsym(bj, dbj, dV);
	    }
	}
	
        if ( matStiffSymmFlag ) {
            m->symmetrized();
        }

        // localize irule contribution into element matrix
        if ( this->giveIntegrationElementLocalCodeNumbers(irlocnum, elem, iRule) ) {
            answer.assemble(* m, irlocnum);
        }
    } // end loop over irules
    
    answer.add(KBbar);
}



} // end namespace oofem
