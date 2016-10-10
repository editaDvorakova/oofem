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
  FloatMatrix help;
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
  
  dr = interp->givedR( help, gp->giveNaturalCoordinates(),
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
    determinant = fabs( this->giveElement()->giveInterpolation()
                       ->giveTransformationJacobian( gp->giveNaturalCoordinates(),
                                                    FEIIGAElementGeometryWrapper( this->giveElement(),
                                                                                 gp->giveIntegrationRule()->giveKnotSpan() ) ) );
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

void Beam2dElementEvaluator :: computeEdgeLoadVectorAt(FloatArray &answer, Load *load,                                              int edge, TimeStep *tStep, ValueModeType mode)
{
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
Beam2dElementEvaluator :: computeNormal (FloatArray &n, FloatArray c, int knotSpan)
{
    FloatMatrix d;
    BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (this->giveElement()->giveInterpolation());
    interp->givedR( d, c, FEIIGAElementGeometryWrapper( this->giveElement(), this->giveElement() ->giveIntegrationRule(knotSpan)->giveKnotSpan() ) );
    n = {-d.at(2,1),  d.at(1,1), 0};

  return;
}
                                
  int Beam2dElementEvaluator :: computeLoadGToLRotationMtrx(FloatMatrix &answer, GaussPoint *gp)
{
  
  this->computeDofsGtoLMatrix(answer,  gp->giveNaturalCoordinates(), gp->giveIntegrationRule()->giveNumber());

  return 1;
}  


  void Beam2dElementEvaluator :: computeDofsGtoLMatrix(FloatMatrix &answer, FloatArray coords, int knotSpan)
  {
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
      /*
    answer.clear();
 
    int numberOfIntegrationRules;
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


    }*/



} // end namespace oofem
