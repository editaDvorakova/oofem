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
		for ( int ir = 0; ir <= knotSpan[0].at(1)-interp->giveDegree(1); ir++ ) { // whole knotSpans until the actual one
#ifdef __PARALLEL_MODE
		    if ( this->giveElement()->giveKnotSpanParallelMode(ir) == Element_remote ) {
			continue;
		    }
		    //fprintf (stderr, "[%d] Computing element.knotspan %d.%d\n", elem->giveDomain()->giveEngngModel()->giveRank(), elem->giveNumber(), ir);
#endif
		    IntegrationRule *iRule;
		    if ( ir < knotSpan[0].at(1)-interp->giveDegree(1) )
			iRule = elem->giveIntegrationRule(ir);
		    else {
			int numberOfGaussPoints = elem->giveIntegrationRule(0)->giveNumberOfIntegrationPoints(); 
			const FloatArray *knotValues; 
			knotValues = interp -> giveKnotValues(1);
			double du = collocationPts.at(icp) - knotValues->at(knotSpan[0](0)-interp->giveDegree(1)+1);

			IGAIntegrationElement *ie = new IGAIntegrationElement(1, this->giveElement(), *knotSpan);
			ie->SetUpPointsOnLine(numberOfGaussPoints, _2dBeam);

			// remap local subelement gp coordinates into knot span coordinates and update integration weight
			for ( GaussPoint *gp: *ie ) {
			    const FloatArray &gpcoords = gp->giveNaturalCoordinates();
			    FloatArray newgpcoords;
			    newgpcoords.resize(1);
			    newgpcoords.at(1) = knotValues->at(knotSpan[0](0)-interp->giveDegree(1)+1) + du * ( gpcoords.at(1) / 2.0 + 0.5 );
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
			
			dr = interp->givedR( help, gp->giveNaturalCoordinates(),
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
    knotSpan(0) = count + interp->giveDegree(1)-1; // CHANGE THIS ! ! !
 
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


double Beam2dElementEvaluatorDsg :: computeVolumeAround(GaussPoint *gp)
{
  // change
    double determinant, weight, area, volume;
    determinant = fabs( this->giveElement()->giveInterpolation()
                       ->giveTransformationJacobian( gp->giveNaturalCoordinates(),
                                                    FEIIGAElementGeometryWrapper( this->giveElement(),
                                                                                 gp->giveIntegrationRule()->giveKnotSpan() ) ) );
    weight      = gp->giveWeight();
    area   = this->giveElement()->giveCrossSection()->give(CS_Area, gp);
    volume   = determinant * weight * area;

    return volume;
}


void Beam2dElementEvaluatorDsg :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
  // check
    static_cast< StructuralCrossSection * >( this->giveElement()->giveCrossSection() )->giveGeneralizedStress_Beam2d(answer, gp, strain, tStep);
}

void Beam2dElementEvaluatorDsg :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
  // ok
  static_cast< StructuralCrossSection * >( this->giveElement()->giveCrossSection() )->give2dBeamStiffMtrx(answer, rMode, gp, tStep);
}

void Beam2dElementEvaluatorDsg :: computeEdgeLoadVectorAt(FloatArray &answer, Load *load,                                              int edge, TimeStep *tStep, ValueModeType mode)
{
  // change
    answer.clear();
    
    if ( edge != 1 ) {
        OOFEM_ERROR("Beam2DElementEvaluator only has 1 edge (the midline) that supports loads. Attempted to apply load to edge %d", edge);
    }
    
    int numberOfIntegrationRules;
    Element *elem = this->giveElement();
    int ndofs = elem->computeNumberOfDofs();
    IntArray irlocnum;
    FloatMatrix n, T;
    FloatArray f, force;
    FloatArray reducedAnswer;
    FloatArray coords;
    double dV;
    answer.resize(ndofs);
    numberOfIntegrationRules = elem->giveNumberOfIntegrationRules();

    // loop over individual integration rules
    for ( int ir = 0; ir < numberOfIntegrationRules; ir++ ) {
#ifdef __PARALLEL_MODE
      if ( this->giveElement()->giveKnotSpanParallelMode(ir) == Element_remote ) {
	continue;
      }
      //fprintf (stderr, "[%d] Computing element.knotspan %d.%d\n", elem->giveDomain()->giveEngngModel()->giveRank(), elem->giveNumber(), ir);
#endif
      reducedAnswer.clear();
      IntegrationRule *iRule = elem->giveIntegrationRule(ir);
      for ( GaussPoint *gp: *iRule ) {
	dV = this->computeVolumeAround(gp);
	this->computeNMatrixAt(n, gp);

	if ( load->giveFormulationType() == Load :: FT_Entity ) {
	  coords = gp->giveNaturalCoordinates();
	  coords.times(2.0);
	  coords.add(-1.0);
	      
	  load->computeValueAt(force, tStep, coords, mode);
	} else {
	// edita 
	  OOFEM_ERROR("Load :: FT_Global unsupported");
	}
	   
	// transform force
	if ( load->giveCoordSystMode() == Load :: CST_Global ) {
	  // transform from global to element local c.s
	  if ( this->computeLoadGToLRotationMtrx(T, gp) ) {
	    force.rotatedWith(T, 'n');
	  }
	}
	f.beTProductOf(n, force);
	reducedAnswer.add(dV,f);
      }

      // localize irule contribution into element matrix
      if ( this->giveIntegrationElementLocalCodeNumbers(irlocnum, elem, iRule) ) {
	answer.assemble(reducedAnswer, irlocnum);
      }
    } // end loop over irules
}

void 
Beam2dElementEvaluatorDsg :: computeNormal (FloatArray &n, FloatArray c, int knotSpan)
{
  // change
    FloatMatrix d;
    BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (this->giveElement()->giveInterpolation());
    interp->givedR( d, c, FEIIGAElementGeometryWrapper( this->giveElement(), this->giveElement() ->giveIntegrationRule(knotSpan)->giveKnotSpan() ) );
    n = {-d.at(2,1),  d.at(1,1), 0};
    n.normalize();

  return;
}
                                
  int Beam2dElementEvaluatorDsg :: computeLoadGToLRotationMtrx(FloatMatrix &answer, GaussPoint *gp)
{
  // change
  this->computeDofsGtoLMatrix(answer,  gp->giveNaturalCoordinates(), gp->giveIntegrationRule()->giveNumber());

  return 1;
}  


  void Beam2dElementEvaluatorDsg :: computeDofsGtoLMatrix(FloatMatrix &answer, FloatArray coords, int knotSpan)
  {
    // change
    // calculates 3x3 transformation matrix at given point and knotspan
    FloatMatrix d;
    double J;
    BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (this->giveElement()->giveInterpolation());

    answer.resize(3,3);
    answer.zero();

    J = interp->giveTransformationJacobian( coords, FEIIGAElementGeometryWrapper( this->giveElement(), this->giveElement() ->giveIntegrationRule(knotSpan)->giveKnotSpan() ) );
    interp->givedR( d, coords, FEIIGAElementGeometryWrapper( this->giveElement(), this->giveElement() ->giveIntegrationRule(knotSpan)->giveKnotSpan() ) );

    answer.at(1,1) = d.at(1,1)/J;
    answer.at(1,2) = d.at(2,1)/J;
    answer.at(2,1) = -d.at(2,1)/J;
    answer.at(2,2) = d.at(1,1)/J;
    answer.at(3,3) = 1.0;

    return;
  }

       
  bool Beam2dElementEvaluatorDsg :: computeGtoLRotationMatrix(FloatMatrix &answer)
  {
    // change
    FloatMatrix GtoL;
    int ndofs = this->giveElement()->computeNumberOfDofs();
    int numberOfIntegrationRules = this->giveElement()->giveNumberOfIntegrationRules();
    int knotSpan;
    FloatArray coords;
    coords.resize(1);

    answer.resize(ndofs, ndofs);
    answer.beUnitMatrix();

    // the first control point
    knotSpan = 0;
    coords.at(1) = 0.0;
    this->computeDofsGtoLMatrix(GtoL, coords, knotSpan);
    answer.setSubMatrix (GtoL, 1, 1);

    // the last control point
    knotSpan = numberOfIntegrationRules-1;
    coords.at(1) = 1.0;
    this->computeDofsGtoLMatrix(GtoL, coords, knotSpan);
    answer.setSubMatrix (GtoL, ndofs-2, ndofs-2);
	  
  return true;
}  


  void Beam2dElementEvaluatorDsg :: computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode)
  {
    // change

    double dens, dV;
    FloatArray force, reducedAnswer, helpForce, ntf;
    FloatMatrix n, T;

    IntArray irlocnum;
    Element *elem = this->giveElement();

    int ndofs = elem->computeNumberOfDofs();
    answer.resize(ndofs);
    answer.zero();

    if ( ( load->giveBCGeoType() != BodyLoadBGT ) || ( load->giveBCValType() != ForceLoadBVT ) ) {
        OOFEM_ERROR("unknown load type");
    }

    // note: force is assumed to be in global coordinate system.
    load->computeComponentArrayAt(force, tStep, mode);

    int numberOfIntegrationRules;
    numberOfIntegrationRules = elem->giveNumberOfIntegrationRules();


    if ( force.giveSize() ) {
      // loop over individual integration rules
      for ( int ir = 0; ir < numberOfIntegrationRules; ir++ ) {
#ifdef __PARALLEL_MODE
	if ( elem->giveKnotSpanParallelMode(ir) == Element_remote ) {
	  continue;
	}
	//fprintf (stderr, "[%d] Computing element.knotspan %d.%d\n", elem->giveDomain()->giveEngngModel()->giveRank(), elem->giveNumber(), ir);
#endif
	reducedAnswer.clear();
	IntegrationRule *iRule = elem->giveIntegrationRule(ir);
	for ( GaussPoint *gp: *iRule ) {
	  this->computeNMatrixAt(n, gp);
	  dV  = this->computeVolumeAround(gp);
	  dens = this->giveElement()->giveCrossSection()->give('d', gp);

	  if ( this->computeLoadGToLRotationMtrx(T, gp) ) {
	    helpForce.beProductOf(T, force);
	  }
	  ntf.beTProductOf(n, helpForce);
	  reducedAnswer.add(dV * dens, ntf);
	}

	// localize irule contribution into element matrix
	if ( this->giveIntegrationElementLocalCodeNumbers(irlocnum, elem, iRule) ) {
	  answer.assemble(reducedAnswer, irlocnum);
	}
      } // end loop over irules


    } else {
      return;
    }


  }

} // end namespace oofem
