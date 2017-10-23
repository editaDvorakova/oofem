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

#include "../sm/Elements/Beams/beam2delementevaluator.h"
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
void Beam2dElementEvaluator :: computeNMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
    FloatArray N;
    FEInterpolation *interp = gp->giveElement()->giveInterpolation();
    interp->evalN( N, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    answer.beNMatrixOf(N, 3); 
}

void Beam2dElementEvaluator :: computeBMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
  FloatMatrix d;
  FloatArray n;
  double J, iJ, dr;
  
  BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (gp->giveElement()->giveInterpolation());
  
  interp->evalN( n, gp->giveNaturalCoordinates(),
		 FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
  
  
  interp->evaldNdx( d, gp->giveNaturalCoordinates(),
		    FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
  
  J = interp->giveTransformationJacobian( gp->giveNaturalCoordinates(),
					  FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
  iJ = 1./J;
  
  dr = this->giveCurvature(gp->giveNaturalCoordinates(),
		       FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
  answer.resize(3, d.giveNumberOfRows() * 3);
  answer.zero();
  
  for ( int i = 1; i <= d.giveNumberOfRows(); i++ ) {
    answer.at(1, i*3 - 2) = d.at(i, 1) * iJ;
    answer.at(1, i*3 - 1) = -n.at(i) * dr;
      
    answer.at(2, i*3 ) = d.at(i, 1) * iJ;

    answer.at(3, i*3 - 2) = n.at(i) * dr;
    answer.at(3, i*3 - 1) = d.at(i, 1) * iJ;
    answer.at(3, i*3 ) = -n.at(i);
  }
}


double Beam2dElementEvaluator :: computeVolumeAround(GaussPoint *gp)
{
    double determinant, weight, area, volume;
    determinant = fabs( this->giveElement()->giveInterpolation()->giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( this->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) ) );
    weight      = gp->giveWeight();
    area   = this->giveElement()->giveCrossSection()->give(CS_Area, gp);
    volume   = determinant * weight * area;

    return volume;
}


void Beam2dElementEvaluator :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    static_cast< StructuralCrossSection * >( this->giveElement()->giveCrossSection() )->giveGeneralizedStress_Beam2d(answer, gp, strain, tStep);
}

void Beam2dElementEvaluator :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
  static_cast< StructuralCrossSection * >( this->giveElement()->giveCrossSection() )->give2dBeamStiffMtrx(answer, rMode, gp, tStep);
}

void Beam2dElementEvaluator :: computeBoundaryEdgeLoadVector(FloatArray &answer, BoundaryLoad *load, int edge, CharType type, ValueModeType mode, TimeStep *tStep, bool global)
{
    answer.clear();
    
    if ( edge != 1 ) {
        OOFEM_ERROR("Beam2DElementEvaluator only has 1 edge (the midline) that supports loads. Attempted to apply load to edge %d", edge);
    }
    
    if ( type != ExternalForcesVector ) {
        return;
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

    if (global) {
      // Loads from sets expects global c.s.
      this->computeGtoLRotationMatrix(T);
      answer.rotatedWith(T, 't');
    }

}


void
Beam2dElementEvaluator :: boundaryEdgeGiveNodes(IntArray& bNodes, int boundary)
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

void 
Beam2dElementEvaluator :: computeNormal (FloatArray &n, FloatArray c, int knotSpan)
{
  // change? - possibly find better way to compute it 
    FloatMatrix d;
    this->givedxds( d, c, FEIIGAElementGeometryWrapper( this->giveElement(), this->giveElement() ->giveIntegrationRule(knotSpan)->giveKnotSpan() ) );
    n = {-d.at(2,1),  d.at(1,1), 0};
    n.normalize();

    return;
}
                                
  int Beam2dElementEvaluator :: computeLoadGToLRotationMtrx(FloatMatrix &answer, GaussPoint *gp)
{
  
  this->computeDofsGtoLMatrix(answer,  gp->giveNaturalCoordinates(), gp->giveIntegrationRule()->giveNumber());
  /*
    // change? - possibly find better way to compute it
    // same result as computeDofsGtoLMatrix
    FloatMatrix d;
    this->givedxds( d, gp->giveNaturalCoordinates(),FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

    answer.resize(3, 3);
    answer.zero();

    FloatArray t;
    t = {d.at(1,1),  d.at(2,1), 0};
    t.normalize();

    answer.at(1,1) =  t.at(1);
    answer.at(1,2) =  t.at(2);
    answer.at(2,1) = -t.at(2);
    answer.at(2,2) =  t.at(1);
    answer.at(3,3) = 1;
*/
  return 1;
}  


  void Beam2dElementEvaluator :: computeDofsGtoLMatrix(FloatMatrix &answer, FloatArray coords, int knotSpan)
  {
    // same as computeLoadGToLRotationMtrx(...)
    // calculates 3x3 transformation matrix at given point and knotspan
    FloatMatrix d;
    double J;
    BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (this->giveElement()->giveInterpolation());

    answer.resize(3,3);
    answer.zero();

    J = interp->giveTransformationJacobian( coords, FEIIGAElementGeometryWrapper( this->giveElement(), this->giveElement() ->giveIntegrationRule(knotSpan)->giveKnotSpan() ) );
    this->givedxds( d, coords, FEIIGAElementGeometryWrapper( this->giveElement(), this->giveElement() ->giveIntegrationRule(knotSpan)->giveKnotSpan() ) );

    answer.at(1,1) = d.at(1,1)/J;
    answer.at(1,2) = d.at(2,1)/J;
    answer.at(2,1) = -d.at(2,1)/J;
    answer.at(2,2) = d.at(1,1)/J;
    answer.at(3,3) = 1.0;

    return;
  }

       
  bool Beam2dElementEvaluator :: computeGtoLRotationMatrix(FloatMatrix &answer)
  {
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


  void Beam2dElementEvaluator :: computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode)
  {

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

    

    
  void Beam2dElementEvaluator :: givedxds(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
  {
      /* calculates coordinates derivatives with respect to curvilinear coordinate s 
	 answer = [dxds, dxdss; dyds, dydss]
      */
    answer.resize(2, 2);
    answer.zero();
    
    FEIIGAElementGeometryWrapper *gw = ( FEIIGAElementGeometryWrapper * ) & cellgeo;
    double x, y;
    FloatMatrix d;
    int i, k, ind, uind;
    BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (this->giveElement()->giveInterpolation());

    interp->evaldNdx(d, lcoords, cellgeo);
    int fsd = interp->giveFsd(); // = 1 
    IntArray span(fsd);
    const FloatArray *vertexCoordsPtr;
    if ( gw->knotSpan ) {
      span = * gw->knotSpan;
    } else {
      for ( i = 0; i < fsd; i++ ) {
	  span(i) = interp->findSpan(interp->giveNumberOfControlPoints(i), interp->giveDegree(), lcoords(i), interp->giveKnotVector()[i]);
      }
    }
    
    int degree = interp->giveDegree();
    uind = span(0) - degree;
    ind = uind + 1;
    
    for ( k = 0; k <= degree; k++ ) {
	vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind + k);
	x = vertexCoordsPtr->at(1);
	y = vertexCoordsPtr->at(2);
	answer.at(1,1) += x*d.at(k+1,1);
	answer.at(2,1) += y*d.at(k+1,1);
	answer.at(1,2) += x*d.at(k+1,2);
	answer.at(2,2) += y*d.at(k+1,2);
	
	//cnt++;
    }
    // dR = (dxdk*dydkk-dydk*dxdkk)/(J*J*J); 
    return;
  }
    
  double Beam2dElementEvaluator :: giveCurvature(const FloatArray &lcoords, const FEICellGeometry &cellgeo)
  {
    double J, dR;
    FloatMatrix d;
    BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (this->giveElement()->giveInterpolation());
    J = interp->evaldNdx(d, lcoords, cellgeo);

    FloatMatrix dx; 
    this->givedxds(dx, lcoords, cellgeo);

    dR = (dx.at(1,1)*dx.at(2,2)-dx.at(2,1)*dx.at(1,2))/(J*J*J); 
    return dR;
  }
    




void
Beam2dElementEvaluator :: computeInternalForces(FloatMatrix &internalForces, int divisions, TimeStep *tStep)
{
    Element *elem = this->giveElement();
    FloatArray strain, stress;

    internalForces.resize(divisions+1, 8);

    int numberOfIntegrationRules = elem->giveNumberOfIntegrationRules();
    FloatArray gpcoords;
    gpcoords.resize(3);
    gpcoords.zero();

    int count = 1;
    FloatArray u;
    elem->computeVectorOf(VM_Total, tStep, u);
    for (int ir = 0; ir < numberOfIntegrationRules; ir++)
	{
	    IntegrationRule *iRule = elem->giveIntegrationRule(ir);
	    GaussPoint *gp = iRule->getIntegrationPoint(0);
	    for (int j = 0; j < divisions+1; j++)
		{
		    gpcoords.at(1) = (double) j*(1./divisions);
		    gp->setNaturalCoordinates(gpcoords);
		    this->computeStrainVector(strain, gp, tStep, u);
		    this->computeStressVector(stress, strain, gp, tStep);
		    internalForces.copySubVectorRow(stress, count, 1);
		    count++;
		}
	}
}



} // end namespace oofem
