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

#include "../sm/Elements/Beams/beam2delementevaluatordsg.h"
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
void Beam2dElementEvaluatorDsg :: computeNMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
    FloatArray N, help;
    IntArray mask;
    Element *elem = this->giveElement();
    BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (this->giveElement()->giveInterpolation());
    int nCP = interp->giveNumberOfControlPoints(1);
    interp->evalN( N, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

    help.resize(nCP);
    help.zero();
    elem->giveInterpolation()->giveKnotSpanBasisFuncMask(* gp->giveIntegrationRule()->giveKnotSpan(), mask);
		
    help.assemble(N, mask);

    answer.beNMatrixOf(help, 3); 
}

const FloatMatrix *
Beam2dElementEvaluatorDsg :: computeDSGMatrix()
{
    if ( !BDSG.isNotEmpty() ) {
	
	BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (this->giveElement()->giveInterpolation());
	Element *elem = this->giveElement();
	int nCP;
	nCP = interp->giveNumberOfControlPoints(1);

	BDSG.resize(2*nCP, 3*nCP);
	BDSG.zero();

	FloatArray n;
	FloatMatrix d, help;
	double J, dr;
	IntArray mask;
	IntArray knotSpan[1];
	double wgp;

	for ( int icp = 2; icp <= nCP; icp++ ) // loop over collocation points 
	    { 		
		// to optimize: copy all knotSpans(irules) <= ir-2
		giveKnotSpanAt(*knotSpan, collocationPts.at(icp));
		// KNOTSPAN
		for ( int ir = 0; ir <= knotSpan[0].at(1)-interp->giveDegree(); ir++ ) { // whole knotSpans until the actual one
#ifdef __PARALLEL_MODE
		    if ( this->giveElement()->giveKnotSpanParallelMode(ir) == Element_remote ) {
			continue;
		    }
		    //fprintf (stderr, "[%d] Computing element.knotspan %d.%d\n", elem->giveDomain()->giveEngngModel()->giveRank(), elem->giveNumber(), ir);
#endif
		    IntegrationRule *iRule;
		    if ( ir < knotSpan[0].at(1)-interp->giveDegree() )
			iRule = elem->giveIntegrationRule(ir);
		    else {
			int numberOfGaussPoints =  elem->giveIntegrationRule(0)->giveNumberOfIntegrationPoints(); 
			numberOfGaussPoints = 16;

			const FloatArray *knotValues; 
			knotValues = interp -> giveKnotValues(1);
			double du = collocationPts.at(icp) - knotValues->at(knotSpan[0](0)-interp->giveDegree()+1);

			IGAIntegrationElement *ie = new IGAIntegrationElement(1, this->giveElement(), *knotSpan);
			ie->SetUpPointsOnLine(numberOfGaussPoints, _2dBeam);

			// remap local subelement gp coordinates into knot span coordinates and update integration weight
			for ( GaussPoint *gp: *ie ) {
			    const FloatArray &gpcoords = gp->giveNaturalCoordinates();
			    FloatArray newgpcoords;
			    newgpcoords.resize(1);
			    newgpcoords.at(1) = knotValues->at(knotSpan[0](0)-interp->giveDegree()+1) + du * ( gpcoords.at(1) / 2.0 + 0.5 );
			    gp->setNaturalCoordinates(newgpcoords);
			    gp->setWeight(gp->giveWeight() / 2.0 * du );
			}
			iRule = ie;
		    }
			
		       
		    // loop over individual integration points
		    for ( GaussPoint *gp: *iRule ) {
			wgp = gp->giveWeight();

			interp->evaldNdx( d, gp->giveNaturalCoordinates(),
				  FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
			
			interp->evalN( n, gp->giveNaturalCoordinates(),
				       FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
			J = interp->giveTransformationJacobian( gp->giveNaturalCoordinates(),
					FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
			//	iJ = 1./J;
			
			dr = this->giveCurvature(gp->giveNaturalCoordinates(),
			    FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
			
			elem->giveInterpolation()->giveKnotSpanBasisFuncMask(* gp->giveIntegrationRule()->giveKnotSpan(), mask);
			
			for ( int i = 1; i <= d.giveNumberOfRows(); i++ ) {
			    // BmDSG
			    BDSG.at(icp, (mask(0) + i-1)*3 - 2) += d.at(i, 1) * wgp;// * iJ; // mask goes from 0, otherwise i must go from zero
			    BDSG.at(icp, (mask(0) + i-1)*3 - 1) += -n.at(i) * dr * J * wgp;
			    
			    // BsDSG
			    BDSG.at(nCP+icp, (mask(0) + i-1)*3 - 2) += n.at(i) * dr * J * wgp;
			    BDSG.at(nCP+icp, (mask(0) + i-1)*3 - 1) += d.at(i, 1) * wgp;// * iJ;
			    BDSG.at(nCP+icp, (mask(0) + i-1)*3 ) += -n.at(i) * J * wgp;
			}
		    }
		} // end loop over irules
		
	    }
    }
    return & BDSG;
}
  
void Beam2dElementEvaluatorDsg :: giveKnotSpanAt(IntArray &knotSpan, double lcoord)
{
    /*
    BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (this->giveElement()->giveInterpolation());
    const double *const *knotVector; 
    knotVector = interp -> giveKnotVector();
   
    int count = 2;
    for (int i = 2; knotVector[0][i]<lcoord ; i++)
	{
	    if (knotVector[0][i]>lcoord)
		    break;
	    count++;
	}
    knotSpan.resize(1);
    knotSpan(0) = count;
 
    return;
*/
    BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (this->giveElement()->giveInterpolation());
    const FloatArray *knotValues; 
    knotValues = interp -> giveKnotValues(1);
   
    int count = 1;
    for (int i = 2; i <= knotValues[0].giveSize(); i++)
	{
	    if (knotValues[0].at(i)>=lcoord)
		    break;
	    count++;
	}
    knotSpan.resize(1);
    knotSpan(0) = count + interp->giveDegree()-1; // CHANGE THIS ! ! !
 
    return;

}

const FloatMatrix * Beam2dElementEvaluatorDsg :: computeInvAMatrix()
{
    if ( !invA.isNotEmpty() ) {
        FloatMatrix A;;
	BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (this->giveElement()->giveInterpolation());
	Element *elem = this->giveElement();
	this->computeCollocationPoints();     // collocationPts
	
	int nCP;
	nCP = interp->giveNumberOfControlPoints(1);

        A.resize(nCP, nCP);
	A.zero();

	IntArray knotSpan;
    
	FloatArray n;
	FloatArray clpCoords;
	clpCoords.resize(1);
	IntArray mask;
	
        for ( int i = 0; i < nCP; i++ ) {
	    clpCoords(0) = collocationPts(i);
	    giveKnotSpanAt(knotSpan, collocationPts(i));
	    interp->evalN( n, clpCoords, FEIIGAElementGeometryWrapper( elem, &knotSpan ) );
	    elem->giveInterpolation()->giveKnotSpanBasisFuncMask( knotSpan, mask);
	 
	    A.copySubVectorRow(n, i+1, mask.at(1));
	}
	invA.beInverseOf(A);
    }
  
  return & invA;
}
  
const FloatArray * Beam2dElementEvaluatorDsg :: computeCollocationPoints()
{
    if ( !collocationPts.isNotEmpty() ) {
	
	BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (this->giveElement()->giveInterpolation());
 
	int nCP;
	nCP = interp->giveNumberOfControlPoints(1);
	
	int degree;
	degree = interp->giveDegree();

	const double *const *knotVector; 
	knotVector = interp -> giveKnotVector();
        collocationPts.resize(nCP);

	for (int i = 0; i < nCP; i++){
	    double sum = 0;
	    for (int j = 0; j < degree; j++){
		sum = sum + knotVector[0][1+i+j];
	    }
	    collocationPts(i) = sum/degree;
	}
    }
  
  return & collocationPts;
}
  
void Beam2dElementEvaluatorDsg :: computeBendingPartOfBMatrix(FloatArray &Bb, GaussPoint *gp)
{
    double J;
    FloatMatrix d;
    
    BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (gp->giveElement()->giveInterpolation());
      
    interp->evaldNdx( d, gp->giveNaturalCoordinates(),
	    FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

    J = interp->giveTransformationJacobian( gp->giveNaturalCoordinates(),
	    FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
      
    Bb.resize(3*d.giveNumberOfRows());
    for ( int i = 1; i <= d.giveNumberOfRows(); i++ ) {
	Bb.at( i*3 ) = d.at(i, 1) * (1./J);
    }
    return;
}
    
int Beam2dElementEvaluatorDsg :: giveIntegrationElementLocalCodeNumbers(IntArray &answer, Element *elem, IntegrationRule *ie)
{
    IntArray mask, nodeDofIDMask, nodalArray;
    int dofmandof;

    // get number of dofs in node
    elem->giveDofManDofIDMask(1, nodeDofIDMask);
    dofmandof = nodeDofIDMask.giveSize();

    nodalArray.resize( dofmandof );

    answer.clear();
    for (int i=1; i <= elem->giveNumberOfNodes(); i++) {
	for ( int j = 1; j <= dofmandof; j++ ) {
	    // nodalArray.at(j) = dofmandof * (elem->giveNode(i)->giveNumber()-1) + j;
	    nodalArray.at(j) = (i-1)*3+j;
	}
	answer.followedBy(nodalArray);	    
    }
    return 1;
}

void Beam2dElementEvaluatorDsg :: computeBMatrixAt(FloatMatrix &answer, GaussPoint *gp)
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
  
    int nCP;
    nCP = interp->giveNumberOfControlPoints(1);


    this->computeInvAMatrix();     // invA
    this->computeDSGMatrix();    // BDSG = [BmDSG; BsDSG]
    

    FloatMatrix d; // check : derivatives according to ksi?!?!?!
    interp->evaldNdx( d, gp->giveNaturalCoordinates(),
	FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

    double J;
    J = interp->giveTransformationJacobian( gp->giveNaturalCoordinates(),
		FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

    // the first derivatives
    ders.beColumnOf(d, 1);
    ders.times(1./J); // times iJ

    patchDers.resize(nCP);
    patchDers.zero();
    patchDers.assemble(ders, mask);
    

    FloatArray Bs;
    FloatMatrix BsDSG;
    BsDSG.beSubMatrixOf (BDSG, nCP+1, 2*nCP, 1, BDSG.giveNumberOfColumns());
    help.beProductOf(invA, BsDSG);
    Bs.beTProductOf(help, patchDers);

    FloatArray Bm;
    FloatMatrix BmDSG;
    BmDSG.beSubMatrixOf (BDSG, 1, nCP, 1, BDSG.giveNumberOfColumns());
    help.beProductOf(invA, BmDSG);
    Bm.beTProductOf(help, patchDers);

    FloatArray Bb, patchBb;
    this->computeBendingPartOfBMatrix(Bb, gp);
    patchBb.resize(nCP*3);
    patchBb.zero();
    patchBb.assemble(Bb, irlocnum);

    answer.resize(3,3 * nCP);
    answer.copySubVectorRow (Bm, 1, 1);
    answer.copySubVectorRow (patchBb, 2, 1);
    answer.copySubVectorRow (Bs, 3, 1);
  
    return;
}




void
Beam2dElementEvaluatorDsg :: boundaryEdgeGiveNodes(IntArray& bNodes, int boundary)
{
    if ( boundary != 1 ) {
        OOFEM_ERROR("Beam2DElementEvaluator only has 1 edge (the midline) that supports loads. Attempted to apply load to edge %d", boundary);
    }

    int nNodes = this->giveElement()->giveNumberOfNodes();

    bNodes.resize(nNodes);
    for( int i = 1; i<=nNodes; i++){
	//	bNodes.at(i) = this->giveElement()->giveNode(i)->giveNumber();
	bNodes.at(i) = i;
    }
}

} // end namespace oofem
