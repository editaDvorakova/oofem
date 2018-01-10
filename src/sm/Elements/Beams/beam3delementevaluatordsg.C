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

#include "../sm/Elements/Beams/beam3delementevaluatordsg.h"
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
#include "timer.h"

namespace oofem {
void Beam3dElementEvaluatorDsg :: computeNMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
  // change
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

    answer.beNMatrixOf(help, 6); 
}

IntegrationRule * Beam3dElementEvaluatorDsg :: setIR(int i)
{
    int ir, deg, ks;
    const IntArray *mp;
    
    BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (this->giveElement()->giveInterpolation());
    
    mp = interp->giveKnotMultiplicity(1);
    deg = interp->giveDegree(1);

    ks = deg-1;
    ir = 0;
    while (ks<i)
	{
	    ks += mp->at(ir+2);
	    ir ++;
	}

    return this->giveElement()->giveIntegrationRule(ir);

}

int Beam3dElementEvaluatorDsg ::giveNumberOfIR(int knotSpan)
{
    int irn, deg, ks;
    const IntArray *mp;
    
    BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (this->giveElement()->giveInterpolation());
    
    mp = interp->giveKnotMultiplicity(1);
    deg = interp->giveDegree(1);

    ks = deg;
    irn = 0;
    while (ks<knotSpan)
	{
	    ks += mp->at(irn+1);
	    irn ++;
	}

    return irn;
}



const FloatMatrix *
Beam3dElementEvaluatorDsg :: computeDSGMatrix()
{
    if ( !BDSG.isNotEmpty() ) {
	//	Timer timer;
	//	timer.startTimer();

	BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (this->giveElement()->giveInterpolation());
	Element *elem = this->giveElement();
	int nCP;
	nCP = interp->giveNumberOfControlPoints(1);

	BDSG.resize(3*nCP, 6*nCP);
	BDSG.zero();

	FloatArray n;
	FloatMatrix d, help;
	double J, kappa, tau;
	IntArray mask;
	IntArray knotSpan[1];
	double wgp;
	int irn;

	for ( int icp = 2; icp <= nCP; icp++ ) // loop over collocation points 
	    { 		
		// to optimize: copy all knotSpans(irules) <= ir-2
		giveKnotSpanAt(*knotSpan, collocationPts.at(icp));
		irn = this->giveNumberOfIR(knotSpan[0][0]);
		irn = knotSpan[0].at(1)-interp->giveDegree(1); 
		// KNOTSPAN
		for ( int ir = 0; ir <= irn/*knotSpan[0].at(1)-interp->giveDegree()*/; ir++ ) { // whole knotSpans including the actual one
#ifdef __PARALLEL_MODE
		    if ( this->giveElement()->giveKnotSpanParallelMode(ir) == Element_remote ) {
			continue;
		    }
		    //fprintf (stderr, "[%d] Computing element.knotspan %d.%d\n", elem->giveDomain()->giveEngngModel()->giveRank(), elem->giveNumber(), ir);
#endif
		    IntegrationRule *iRule;
		    if ( ir < irn )// knotSpans until the actual one
			iRule = elem->giveIntegrationRule(ir);
		    // iRule = this->setIR(ir);
		    else {
			int numberOfGaussPoints = elem->giveIntegrationRule(0)->giveNumberOfIntegrationPoints(); 
			// numberOfGaussPoints *= 2;
			const FloatArray *knotValues; 
			knotValues = interp -> giveKnotValues(1);
			
			// !!!!! make sure this is right
			double du = collocationPts.at(icp) - knotValues->at(irn+1);
			//double du = collocationPts.at(icp) - knotValues->at(knotSpan[0](0)-interp->giveDegree()+1);

			IGAIntegrationElement *ie = new IGAIntegrationElement(1, this->giveElement(), *knotSpan);
			ie->SetUpPointsOnLine(numberOfGaussPoints, _3dBeam);

			// remap local subelement gp coordinates into knot span coordinates and update integration weight
			for ( GaussPoint *gp: *ie ) {
			    const FloatArray &gpcoords = gp->giveNaturalCoordinates();
			    FloatArray newgpcoords;
			    newgpcoords.resize(1);
			    newgpcoords.at(1) = knotValues->at(irn+1) + du * ( gpcoords.at(1) / 2.0 + 0.5 );
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
			
			kappa = this->giveCurvature( gp->giveNaturalCoordinates(),
			      FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

			tau = this->giveTorsion( gp->giveNaturalCoordinates(),
			      FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

			
			elem->giveInterpolation()->giveKnotSpanBasisFuncMask(* gp->giveIntegrationRule()->giveKnotSpan(), mask);
			
			for ( int i = 1; i <= d.giveNumberOfRows(); i++ ) {
			    // BmDSG
			    BDSG.at(icp, (mask(0) + i-1)*6 - 5) += d.at(i, 1) * wgp; // mask goes from 0, otherwise i must go from zero
			    BDSG.at(icp, (mask(0) + i-1)*6 - 4) += -n.at(i) * kappa * J * wgp; // check this column 
 
			    // BsDSG
			    BDSG.at(nCP+icp, (mask(0) + i-1)*6 - 5) += n.at(i) * kappa * J * wgp; // check this column
			    BDSG.at(nCP+icp, (mask(0) + i-1)*6 - 4) += d.at(i, 1) * wgp;
			    BDSG.at(nCP+icp, (mask(0) + i-1)*6 - 3) += -n.at(i) * tau * J * wgp;
			    BDSG.at(nCP+icp, (mask(0) + i-1)*6 ) += -n.at(i) * J * wgp;
        
			    // BbDSG
			    BDSG.at(2*nCP+icp, (mask(0) + i-1)*6 - 4) += n.at(i) * tau * J * wgp;
			    BDSG.at(2*nCP+icp, (mask(0) + i-1)*6 - 3) += d.at(i, 1) * wgp;
			    BDSG.at(2*nCP+icp, (mask(0) + i-1)*6 - 1) += n.at(i) * J * wgp; // check this column

	
			}
		    }
		} // end loop over irules
		
	    }
	// timer.stopTimer();
	// OOFEM_LOG_INFO( "BDSG: %.2fs\n", timer.getUtime() );
    }
    return & BDSG;
}
  
void Beam3dElementEvaluatorDsg :: giveKnotSpanAt(IntArray &knotSpan, double lcoord)
{/*
    BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (this->giveElement()->giveInterpolation());
    knotSpan.resize(1);

       if (lcoord == 1) lcoord -= 0.00001;

    //    const double *knt = interp->giveKnotVector()[0];
    knotSpan(0) = interp->findSpan(interp->giveNumberOfControlPoints(1), interp->giveDegree(), lcoord, interp->giveKnotVector()[0]);
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
    knotSpan(0) = count + interp->giveDegree(1)-1; // CHANGE THIS ! ! !
 
    return;

}

const FloatMatrix * Beam3dElementEvaluatorDsg :: computeInvAMatrix()
{
    if ( !invA.isNotEmpty() ) {
	//	Timer timer;
	//	timer.startTimer();

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

	//	timer.stopTimer();
	//	OOFEM_LOG_INFO( "invA: %.2fs\n", timer.getUtime() );
    }
  
  return & invA;
}
  
const FloatArray * Beam3dElementEvaluatorDsg :: computeCollocationPoints()
{
    if ( !collocationPts.isNotEmpty() ) {
	
	BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (this->giveElement()->giveInterpolation());
 
	int nCP;
	nCP = interp->giveNumberOfControlPoints(1);
	
	int degree;
	degree = interp->giveDegree(1);

	const FloatArray *knotVector; 
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
  
void Beam3dElementEvaluatorDsg :: computeBendingPartOfBMatrix(FloatMatrix &Bb, GaussPoint *gp)
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

    Bb.resize(3, 6*d.giveNumberOfRows());
    for ( int i = 1; i <= d.giveNumberOfRows(); i++ ) {

	Bb.at(1, i*6 - 2) = d.at(i, 1) / J;
	Bb.at(1, i*6 - 1) = -kappa * n.at(i);
      
	Bb.at(2, i*6 - 2) = kappa * n.at(i);
	Bb.at(2, i*6 - 1) = d.at(i, 1) / J;
	Bb.at(2, i*6) = -tau * n.at(i);

	Bb.at(3, i*6 - 1) = tau * n.at(i);
	Bb.at(3, i*6) = d.at(i, 1) / J;
    }
    return;
}
    
int Beam3dElementEvaluatorDsg :: giveIntegrationElementLocalCodeNumbers(IntArray &answer, Element *elem, IntegrationRule *ie)
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
	    nodalArray.at(j) = (i-1)*6+j;
	}
	answer.followedBy(nodalArray);	    
    }
    return 1;
    }

void Beam3dElementEvaluatorDsg :: computeBMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
    //    Timer timer;
    //    timer.startTimer();

    FloatArray ders, patchDers;  // row

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
    


    FloatArray help;
    FloatArray Bsn;
    FloatMatrix BsnDSG;
    BsnDSG.beSubMatrixOf (BDSG, nCP+1, 2*nCP, 1, BDSG.giveNumberOfColumns());
    // help.beProductOf(invA, BsnDSG);
    // Bsn.beTProductOf(help, patchDers);
    help.beTProductOf(invA, patchDers);
    Bsn.beTProductOf(BsnDSG, help);

    FloatArray Bm;
    FloatMatrix BmDSG;
    BmDSG.beSubMatrixOf (BDSG, 1, nCP, 1, BDSG.giveNumberOfColumns());
    //  help.beProductOf(invA, BmDSG);
    //  Bm.beTProductOf(help, patchDers);
    help.beTProductOf(invA, patchDers);
    Bm.beTProductOf(BmDSG, help);

    FloatArray Bsb;
    FloatMatrix BsbDSG;
    BsbDSG.beSubMatrixOf (BDSG, 2*nCP+1, 3*nCP, 1, BDSG.giveNumberOfColumns());

    // help.beProductOf(invA, BsbDSG);
    // Bsb.beTProductOf(help, patchDers);
    help.beTProductOf(invA, patchDers);
    Bsb.beTProductOf(BsbDSG, help);



    FloatMatrix Bb, patchBb;
    IntArray rloc = {1,2,3};
    this->computeBendingPartOfBMatrix(Bb, gp);
    patchBb.resize(3,nCP*6);
    patchBb.zero();
    patchBb.assemble(Bb, rloc, irlocnum);


    answer.resize(6,6 * nCP);
    answer.copySubVectorRow (Bm, 1, 1);
    answer.copySubVectorRow (Bsn, 2, 1);
    answer.copySubVectorRow (Bsb, 3, 1);
  
    answer.setSubMatrix (patchBb, 4, 1);

    //    timer.stopTimer();
    //    OOFEM_LOG_INFO( "computeBmatrixAt: elem %d gp %d gpcoords %f  %.4fs\n", this->giveElement()->giveNumber(), gp->giveNumber(), gp->giveNaturalCoordinate(1),timer.getUtime() );
    return;
}





    /*
  void Beam3dElementEvaluatorDsg :: computeDofsGtoLMatrix(FloatMatrix &answer, FloatArray coords,  const FEICellGeometry &cellgeo)
{
    FloatMatrix lcs;
    
    answer.resize(6, 6);
    answer.zero();

    this->giveLocalCoordinateSystem(lcs, coords, cellgeo);
    for ( int i = 1; i <= 3; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            answer.at(i, j) = lcs.at(i, j);
            answer.at(i + 3, j + 3) = lcs.at(i, j);
        }
    }

  }
    */
       
    /*
int
Beam3dElementEvaluatorDsg  :: giveLocalCoordinateSystem(FloatMatrix &answer, FloatArray lcoords, const FEICellGeometry &cellgeo)
// move this function to beam3delementevaluator????
//
// returns a unit vectors of local coordinate system at element
// stored rowwise (mainly used by some materials with ortho and anisotrophy)
//
{
    FloatArray t, n, b, help, N;
    double J, norm2;
    FloatMatrix dx, d;

    BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (this->giveElement()->giveInterpolation());
    interp->evalN(N, lcoords, cellgeo);
    J = interp->evaldNdx(d, lcoords, cellgeo);
 
    this->givedxds(dx, lcoords, cellgeo);

    t.resize(3);
    t.at(1) = dx.at(1,1)/J;
    t.at(2) = dx.at(2,1)/J;
    t.at(3) = dx.at(3,1)/J;
    t.normalize();

    FloatArray gcoords;
    interp->local2global(gcoords, lcoords, cellgeo); 

    help.beColumnOf(dx,1);
    norm2 = help.computeSquaredNorm();


    n.resize(3);
    double d1d2;
    d1d2 = dx.at(1,1)*dx.at(1,2) +  dx.at(2,1)*dx.at(2,2) + dx.at(3,1)*dx.at(3,2);

    n.at(1) =  (dx.at(1,2)* norm2  - dx.at(1,1) * d1d2 )/pow(J,4);
    n.at(2) =  (dx.at(2,2)* norm2  - dx.at(2,1) * d1d2 )/pow(J,4);
    n.at(3) =  (dx.at(3,2)* norm2  - dx.at(3,1) * d1d2 )/pow(J,4);
    n.normalize();

    b.beVectorProductOf(t,n);


    answer.resize(3, 3);
    answer.copySubVectorRow (t, 1, 1);
    answer.copySubVectorRow (n, 2, 1);
    answer.copySubVectorRow (b, 3, 1);

    return 1;
}
    */

} // end namespace oofem
