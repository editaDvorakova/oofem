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

#include "../sm/Elements/Beams/beam3delementevaluator.h"
#include "../sm/Elements/igaelements.h"
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
void Beam3dElementEvaluator :: computeNMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
    FloatArray N;
    FEInterpolation *interp = gp->giveElement()->giveInterpolation();
    interp->evalN( N, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    answer.beNMatrixOf(N, 6); 
}


void Beam3dElementEvaluator :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
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
	if (elem->giveIntegrationType() == 0) {
	    // FULL INTEGRATION
	    iRule = elem->giveIntegrationRule(ir);
	    // loop over individual integration points
	    for ( GaussPoint *gp: *iRule ) {
		double dV = this->computeVolumeAround(gp);
		this->computeBMatrixAt(bj, gp);
		this->computeConstitutiveMatrixAt(d, rMode, gp, tStep);

		dbj.beProductOf(d, bj);
		if ( matStiffSymmFlag ) {
		    m->plusProductSymmUpper(bj, dbj, dV);
		} else {
		    m->plusProductUnsym(bj, dbj, dV);
		}
	    }
	} else if (elem->giveIntegrationType() == 1) {
	    // FULL INTEGRATION
	    iRule = elem->giveIntegrationRule(ir);
	    // loop over individual integration points
	    for ( GaussPoint *gp: *iRule ) {
		double dV = this->computeVolumeAround(gp);
		this->computeB1MatrixAt(bj, gp);
		this->computeConstitutiveMatrixAt(d, rMode, gp, tStep);

		dbj.beProductOf(d, bj);
		if ( matStiffSymmFlag ) {
		    m->plusProductSymmUpper(bj, dbj, dV);
		} else {
		    m->plusProductUnsym(bj, dbj, dV);
		}
	    }

	    // REDUCED INTEGRATION
	    iRule = elem->giveReducedIntegrationRule(ir);
	    // loop over individual integration points
	    for ( GaussPoint *gp: *iRule ) {
		double dV = this->computeVolumeAround(gp);
		this->computeBMatrixAt(bj, gp);
		this->computeConstitutiveMatrixAt(d, rMode, gp, tStep);

		dbj.beProductOf(d, bj);
		if ( matStiffSymmFlag ) {
		    m2.plusProductSymmUpper(bj, dbj, dV);
		} else {
		    m2.plusProductUnsym(bj, dbj, dV);
		}
	    }
	    m->add(m2);

	} else {
	    OOFEM_ERROR("Unsupported integrationType == %d", elem->giveIntegrationType());
	}
	
        if ( matStiffSymmFlag ) {
            m->symmetrized();
        }

        // localize irule contribution into element matrix
        if ( this->giveIntegrationElementLocalCodeNumbers(irlocnum, elem, iRule) ) {
            answer.assemble(* m, irlocnum);
        }
    } // end loop over irules
}



void Beam3dElementEvaluator :: computeBMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
    FloatMatrix B1, B2;
    this->computeB1MatrixAt(B1, gp);
    this->computeB2MatrixAt(B2, gp);
    
    answer = B1;
    answer.add(B2);

/*
    FloatMatrix d;
    FloatMatrix help;
    FloatArray n;
    double J, iJ, kappa, tau;
  
    BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (gp->giveElement()->giveInterpolation());
  
    interp->evalN( n, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    interp->evaldNdx( d, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );  
    J = interp->giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    iJ = 1./J;
  
    kappa = this->giveCurvature( gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

    tau = this->giveTorsion( gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

    answer.resize(6, d.giveNumberOfRows() * 6);
    answer.zero();
  
    for ( int i = 1; i <= d.giveNumberOfRows(); i++ ) {
	answer.at(1, i*6 - 5) = d.at(i, 1) * iJ;
	answer.at(1, i*6 - 4) = -kappa * n.at(i);
      
	answer.at(2, i*6 - 5) = kappa * n.at(i);
	answer.at(2, i*6 - 4) = d.at(i, 1) * iJ;
	answer.at(2, i*6 - 3) = -tau * n.at(i);
	answer.at(2, i*6) = -n.at(i);

	answer.at(3, i*6 - 4) = tau * n.at(i);
	answer.at(3, i*6 - 3) = d.at(i, 1) * iJ;
	answer.at(3, i*6 - 1) = n.at(i);

	answer.at(4, i*6 - 2) = d.at(i, 1) * iJ;
	answer.at(4, i*6 - 1) = -kappa * n.at(i);
      
	answer.at(5, i*6 - 2) = kappa * n.at(i);
	answer.at(5, i*6 - 1) = d.at(i, 1) * iJ;
	answer.at(5, i*6) = -tau * n.at(i);

	answer.at(6, i*6 - 1) = tau * n.at(i);
	answer.at(6, i*6) = d.at(i, 1) * iJ;

	}*/

}

void Beam3dElementEvaluator :: computeB1MatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
    FloatMatrix d;
    FloatMatrix help;
    FloatArray n;
    double J, iJ, kappa, tau;
  
    BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (gp->giveElement()->giveInterpolation());
  
    interp->evalN( n, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    interp->evaldNdx( d, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );  
    J = interp->giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    iJ = 1./J;
  
    kappa = this->giveCurvature( gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

    tau = this->giveTorsion( gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

    answer.resize(6, d.giveNumberOfRows() * 6);
    answer.zero();
  
    for ( int i = 1; i <= d.giveNumberOfRows(); i++ ) {
	answer.at(1, i*6 - 5) = d.at(i, 1) * iJ;
	answer.at(1, i*6 - 4) = -kappa * n.at(i);
      
	answer.at(2, i*6 - 5) = kappa * n.at(i);
	answer.at(2, i*6 - 4) = d.at(i, 1) * iJ;
	answer.at(2, i*6 - 3) = -tau * n.at(i);
	answer.at(2, i*6) = -n.at(i);

	answer.at(3, i*6 - 4) = tau * n.at(i);
	answer.at(3, i*6 - 3) = d.at(i, 1) * iJ;
	answer.at(3, i*6 - 1) = n.at(i);
	/*
	answer.at(4, i*6 - 2) = d.at(i, 1) * iJ;
	answer.at(4, i*6 - 1) = -kappa * n.at(i);
      
	answer.at(5, i*6 - 2) = kappa * n.at(i);
	answer.at(5, i*6 - 1) = d.at(i, 1) * iJ;
	answer.at(5, i*6) = -tau * n.at(i);

	answer.at(6, i*6 - 1) = tau * n.at(i);
	answer.at(6, i*6) = d.at(i, 1) * iJ;
	*/
    }
}

void Beam3dElementEvaluator :: computeB2MatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
    FloatMatrix d;
    FloatMatrix help;
    FloatArray n;
    double J, iJ, kappa, tau;
  
    BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (gp->giveElement()->giveInterpolation());
  
    interp->evalN( n, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    interp->evaldNdx( d, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );  
    J = interp->giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    iJ = 1./J;
  
    kappa = this->giveCurvature( gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

    tau = this->giveTorsion( gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

    answer.resize(6, d.giveNumberOfRows() * 6);
    answer.zero();
  
    for ( int i = 1; i <= d.giveNumberOfRows(); i++ ) {
	/*
	answer.at(1, i*6 - 5) = d.at(i, 1) * iJ;
	answer.at(1, i*6 - 4) = -kappa * n.at(i);
      
	answer.at(2, i*6 - 5) = kappa * n.at(i);
	answer.at(2, i*6 - 4) = d.at(i, 1) * iJ;
	answer.at(2, i*6 - 3) = -tau * n.at(i);
	answer.at(2, i*6) = -n.at(i);

	answer.at(3, i*6 - 4) = tau * n.at(i);
	answer.at(3, i*6 - 3) = d.at(i, 1) * iJ;
	answer.at(3, i*6 - 1) = n.at(i);
	*/
	answer.at(4, i*6 - 2) = d.at(i, 1) * iJ;
	answer.at(4, i*6 - 1) = -kappa * n.at(i);
      
	answer.at(5, i*6 - 2) = kappa * n.at(i);
	answer.at(5, i*6 - 1) = d.at(i, 1) * iJ;
	answer.at(5, i*6) = -tau * n.at(i);

	answer.at(6, i*6 - 1) = tau * n.at(i);
	answer.at(6, i*6) = d.at(i, 1) * iJ;
    }
}


double Beam3dElementEvaluator :: computeVolumeAround(GaussPoint *gp)
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


void Beam3dElementEvaluator :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    static_cast< StructuralCrossSection * >( this->giveElement()->giveCrossSection() )->giveGeneralizedStress_Beam3d(answer, gp, strain, tStep);
}

void Beam3dElementEvaluator :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    static_cast< StructuralCrossSection * >( this->giveElement()->giveCrossSection() )->give3dBeamStiffMtrx(answer, rMode, gp, tStep);
}

void Beam3dElementEvaluator :: computeBoundaryEdgeLoadVector(FloatArray &answer, BoundaryLoad *load, int edge, CharType type, ValueModeType mode, TimeStep *tStep, bool global)
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


// oofeg support
void 
Beam3dElementEvaluator :: computeNormal (FloatArray &n, FloatArray c, int knotSpan)
{
// edita : 3D not supported
FloatMatrix d;
// BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (this->giveElement()->giveInterpolation());
this->givedxds( d, c, FEIIGAElementGeometryWrapper( this->giveElement(), this->giveElement() ->giveIntegrationRule(knotSpan)->giveKnotSpan() ) );
n = {-d.at(2,1),  d.at(1,1), 0};

return;
}

                                
int Beam3dElementEvaluator :: computeLoadGToLRotationMtrx(FloatMatrix &answer, GaussPoint *gp)
{
  
    // this->computeDofsGtoLMatrix(answer,  gp->giveNaturalCoordinates(), gp->giveIntegrationRule()->giveNumber());

    
    FloatMatrix lcs;

    answer.resize(6, 6);
    answer.zero();

    this->giveLocalCoordinateSystem(lcs, gp->giveNaturalCoordinates(),  FEIIGAElementGeometryWrapper( this->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    for ( int i = 1; i <= 3; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            answer.at(i, j) = lcs.at(i, j);
            answer.at(3 + i, 3 + j) = lcs.at(i, j);
        }
    }

    return 1;
}  


void Beam3dElementEvaluator :: computeDofsGtoLMatrix(FloatMatrix &answer, FloatArray coords, int knotSpan)
{
    // edita : change to 3d
    // calculates 3x3 transformation matrix at given point and knotspan
    FloatMatrix d;
    double J;
    BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (this->giveElement()->giveInterpolation());

    answer.resize(6,6);
    answer.zero();

    J = interp->giveTransformationJacobian( coords, FEIIGAElementGeometryWrapper( this->giveElement(), this->giveElement() ->giveIntegrationRule(knotSpan)->giveKnotSpan() ) );
    this->givedxds( d, coords, FEIIGAElementGeometryWrapper( this->giveElement(), this->giveElement() ->giveIntegrationRule(knotSpan)->giveKnotSpan() ) );

    answer.at(1,1) = d.at(1,1)/J;
    answer.at(1,2) = d.at(2,1)/J;
    answer.at(2,1) = -d.at(2,1)/J;
    answer.at(2,2) = d.at(1,1)/J;
    answer.at(3,3) = 1.0;
    answer.at(4,4) = 1.0;
    answer.at(5,5) = 1.0;
    answer.at(6,6) = 1.0;

    /*
    //this is maybe the same
    FloatMatrix lcs;
    
    answer.resize(6, 6);
    answer.zero();

    this->giveLocalCoordinateSystem(lcs, coords, cellgeo);
    for ( int i = 1; i <= 3; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            answer.at(i, j) = lcs.at(i, j);
            answer.at(i + 3, j + 3) = lcs.at(i, j);
        }
    }*/
    return;
}

 

   
  void Beam3dElementEvaluator :: computeDofsGtoLMatrix(FloatMatrix &answer, FloatArray coords,  const FEICellGeometry &cellgeo)
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
     
int
Beam3dElementEvaluator  :: giveLocalCoordinateSystem(FloatMatrix &answer, FloatArray lcoords, const FEICellGeometry &cellgeo)
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
    if (n.computeNorm() < 10e-6){
	IGAElement *elem = dynamic_cast<IGAElement*> (this->giveElement()); 
	elem->giveZaxis();
	FloatArray z = elem->giveZaxis();
	if (z.giveSize()){
	    //FloatArray zaxis = {0,0,1}; // read from input file
	    n.beVectorProductOf(z,t);
	}else{
	    OOFEM_ERROR("zaxis must be specified in case of straight beam");
	}
    }
    n.normalize();

    b.beVectorProductOf(t,n);


    answer.resize(3, 3);
    answer.copySubVectorRow (t, 1, 1);
    answer.copySubVectorRow (n, 2, 1);
    answer.copySubVectorRow (b, 3, 1);

    return 1;
}
    
bool Beam3dElementEvaluator :: computeGtoLRotationMatrix(FloatMatrix &answer)
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
    this->computeDofsGtoLMatrix(GtoL, coords, FEIIGAElementGeometryWrapper( this->giveElement(), this->giveElement() ->giveIntegrationRule(knotSpan)->giveKnotSpan() )); // edita : check multiplicity
    // this->computeDofsGtoLMatrix(GtoL, coords, knotSpan);
    answer.setSubMatrix (GtoL, 1, 1);

    // the last control point
    knotSpan = numberOfIntegrationRules-1;
    coords.at(1) = 1.0;
    this->computeDofsGtoLMatrix(GtoL, coords, FEIIGAElementGeometryWrapper( this->giveElement(), this->giveElement() ->giveIntegrationRule(knotSpan)->giveKnotSpan() )); // edita : check multiplicity
    // this->computeDofsGtoLMatrix(GtoL, coords, knotSpan);
    answer.setSubMatrix (GtoL, ndofs-5, ndofs-5);
	  
    return true;
}  


void Beam3dElementEvaluator :: computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode)
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

void
Beam3dElementEvaluator :: boundaryEdgeGiveNodes(IntArray& bNodes, int boundary)
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
    
void Beam3dElementEvaluator :: givedxds(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    /* calculates coordinates derivatives with respect to curvilinear coordinate s 
       answer = [dxds, dxdss; dyds, dydss]
    */
    answer.resize(3, 3);
    answer.zero();
    
    FEIIGAElementGeometryWrapper *gw = ( FEIIGAElementGeometryWrapper * ) & cellgeo;
    double x, y, z;
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
	z = vertexCoordsPtr->at(3);
	answer.at(1,1) += x*d.at(k+1,1);
	answer.at(2,1) += y*d.at(k+1,1);
	answer.at(3,1) += z*d.at(k+1,1);
	answer.at(1,2) += x*d.at(k+1,2);
	answer.at(2,2) += y*d.at(k+1,2);
	answer.at(3,2) += z*d.at(k+1,2);
	answer.at(1,3) += x*d.at(k+1,3);
	answer.at(2,3) += y*d.at(k+1,3);
	answer.at(3,3) += z*d.at(k+1,3);
	
	//cnt++;
    }
    return;
}
    
void Beam3dElementEvaluator :: giveTangent(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double J;
    FloatMatrix d;
    BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (this->giveElement()->giveInterpolation());
    J = interp->evaldNdx(d, lcoords, cellgeo);

    FloatMatrix dx; 
    this->givedxds(dx, lcoords, cellgeo);

    answer.resize(3);
    answer.at(1) = dx.at(1,1)/J;
    answer.at(2) = dx.at(2,1)/J;
    answer.at(3) = dx.at(3,1)/J;
    answer.normalize();
 
    return;
}
    
void Beam3dElementEvaluator :: giveNormal(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double J;
    FloatMatrix d;
    BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (this->giveElement()->giveInterpolation());
    J = interp->evaldNdx(d, lcoords, cellgeo);

    FloatMatrix dx; 
    this->givedxds(dx, lcoords, cellgeo);

    answer.resize(3);
    answer.at(1) = dx.at(1,2)/J;
    answer.at(2) = dx.at(2,2)/J;
    answer.at(3) = dx.at(3,2)/J;
    answer.normalize();

    return;
}

void Beam3dElementEvaluator :: giveBinormal(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray t, n;
    this->giveTangent(t, lcoords, cellgeo);
    this->giveNormal(n, lcoords, cellgeo);

    answer.beVectorProductOf(t, n);

    return;
}

double Beam3dElementEvaluator :: giveCurvature(const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double kappa, den, num;
    FloatMatrix d;
    FloatArray d1, d2, help;
    BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (this->giveElement()->giveInterpolation());
    interp->evaldNdx(d, lcoords, cellgeo);

    FloatMatrix dx; 
    this->givedxds(dx, lcoords, cellgeo);

    d1.beColumnOf(dx,1);
    d2.beColumnOf(dx,2);
    help.beVectorProductOf(d1,d2);
    num = help.computeNorm();
    den = d1.computeSquaredNorm() *  d1.computeNorm();
    kappa = num/den;

    //  kappa = sqrt ( pow(dx.at(3,2)*dx.at(2,1)-dx.at(2,2)*dx.at(3,1),2) + pow(dx.at(1,2)*dx.at(3,1)-dx.at(3,2)*dx.at(1,1),2) + pow(dx.at(2,2)*dx.at(1,1)-dx.at(1,2)*dx.at(2,1),2) ) / (J*J*J);
    return kappa;
}
    
    
double Beam3dElementEvaluator :: giveTorsion(const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double tau, den, num;
    FloatMatrix dx, dxhelp; 
    FloatArray d1, d2, d3, help, n0;
    BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (this->giveElement()->giveInterpolation());
    interp->evaldNdx(dx, lcoords, cellgeo);
    interp->evalN(n0, lcoords, cellgeo);
    this->givedxds(dx, lcoords, cellgeo);

    num = dx.giveDeterminant(); // mixed product of d1, d2, d3

    d1.beColumnOf(dx,1);
    d2.beColumnOf(dx,2);
    d3.beColumnOf(dx,3);

    help.beVectorProductOf(d1,d2);
    den = help.computeSquaredNorm();
    if (num == 0)// change the calculation of tau
	tau = 0;
    else
	tau = num/den;


    // tau = ( dx.at(1,3)*(dx.at(3,2)*dx.at(2,1)-dx.at(2,2)*dx.at(3,1)) + dx.at(2,3)*(dx.at(1,2)*dx.at(3,1)-dx.at(3,2)*dx.at(1,1)) + dx.at(3,3)*(dx.at(2,2)*dx.at(1,1)-dx.at(1,2)*dx.at(2,1)) ) / ( pow(dx.at(3,2)*dx.at(2,1)-dx.at(2,2)*dx.at(3,1),2) + pow(dx.at(1,2)*dx.at(3,1)-dx.at(3,2)*dx.at(1,1),2) + pow(dx.at(2,2)*dx.at(1,1)-dx.at(1,2)*dx.at(2,1),2) );
    return tau;
}



} // end namespace oofem
