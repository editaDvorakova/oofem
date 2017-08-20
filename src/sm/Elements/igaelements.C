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

#include "../sm/Elements/igaelements.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "domain.h"
#include "node.h"
#include "element.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "matresponsemode.h"
#include "crosssection.h"
#include "mathfem.h"
#include "pointload.h"
#include "load.h"
#include "iga/iga.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "oofegutils.h"
#endif

namespace oofem {
REGISTER_Element(BsplinePlaneStressElement);
REGISTER_Element(NURBSPlaneStressElement);
REGISTER_Element(TSplinePlaneStressElement);
REGISTER_Element(NURBSSpace3dElement);
REGISTER_Element(NURBSBeam2dElement);
REGISTER_Element(NURBSBeam3dElement);
REGISTER_Element(NURBSBeam2dElementDsg);
REGISTER_Element(NURBSBeam3dElementDsg);


  BsplinePlaneStressElement :: BsplinePlaneStressElement(int n, Domain *aDomain) : IGAElement(n, aDomain), PlaneStressStructuralElementEvaluator(), interpolation(2, 2) { }


IRResultType BsplinePlaneStressElement :: initializeFrom(InputRecord *ir)
{
    //PlaneStressStructuralElementEvaluator::initializeFrom(ir);
    return IGAElement :: initializeFrom(ir);
}


int BsplinePlaneStressElement :: checkConsistency()
{
    BSplineInterpolation *interpol = static_cast< BSplineInterpolation * >( this->giveInterpolation() );
    if ( giveNumberOfDofManagers() != interpol->giveNumberOfControlPoints(1) * interpol->giveNumberOfControlPoints(2) ) {
        OOFEM_WARNING("number of control points mismatch");
        return 0;
    }
    return 1;
}



  NURBSPlaneStressElement :: NURBSPlaneStressElement(int n, Domain *aDomain) : IGAElement(n, aDomain), PlaneStressStructuralElementEvaluator(), interpolation(2, 2) { }


IRResultType NURBSPlaneStressElement :: initializeFrom(InputRecord *ir)
{
    //PlaneStressStructuralElementEvaluator::initializeFrom(ir);
    return IGAElement :: initializeFrom(ir);
}

int NURBSPlaneStressElement :: checkConsistency()
{
    NURBSInterpolation *interpol = static_cast< NURBSInterpolation * >( this->giveInterpolation() );
    if ( giveNumberOfDofManagers() != interpol->giveNumberOfControlPoints(1) * interpol->giveNumberOfControlPoints(2) ) {
        OOFEM_WARNING("number of control points mismatch");
        return 0;
    }
    return 1;
}


  TSplinePlaneStressElement :: TSplinePlaneStressElement(int n, Domain *aDomain) : IGATSplineElement(n, aDomain), PlaneStressStructuralElementEvaluator(), interpolation(2, 2) { }




  NURBSSpace3dElement :: NURBSSpace3dElement(int n, Domain *aDomain) : IGAElement(n, aDomain), Space3dStructuralElementEvaluator(), interpolation(3, 3) { }


IRResultType NURBSSpace3dElement :: initializeFrom(InputRecord *ir)
{
    //PlaneStressStructuralElementEvaluator::initializeFrom(ir);
    return IGAElement :: initializeFrom(ir);
}


int NURBSSpace3dElement :: checkConsistency()
{
    NURBSInterpolation *interpol = static_cast< NURBSInterpolation * >( this->giveInterpolation() );
    if ( giveNumberOfDofManagers() != interpol->giveNumberOfControlPoints(1) * interpol->giveNumberOfControlPoints(2) * interpol->giveNumberOfControlPoints(3) ) {
        OOFEM_WARNING("number of control points mismatch");
        return 0;
    }
    return 1;
}


  NURBSBeam2dElementDsg :: NURBSBeam2dElementDsg(int n, Domain *aDomain) : IGAElement(n, aDomain), Beam2dElementEvaluatorDsg(), interpolation(3, 1) { }
  // interpolation(nsd, fsd)

IRResultType NURBSBeam2dElementDsg :: initializeFrom(InputRecord *ir)
{
    return IGAElement :: initializeFrom(ir);
}


int NURBSBeam2dElementDsg :: checkConsistency()
{
    NURBSInterpolationLine2d *interpol = static_cast< NURBSInterpolationLine2d * >( this->giveInterpolation() );
    if ( giveNumberOfDofManagers() != interpol->giveNumberOfControlPoints(1) ) {
        OOFEM_WARNING("number of control points mismatch");
        return 0;
    }
    return 1;
}

void NURBSBeam2dElementDsg :: computeLoadVector(FloatArray &answer, BodyLoad *load, CharType type, ValueModeType mode, TimeStep *tStep)
{
    if ( type != ExternalForcesVector ) {
        answer.clear();
        return;
    }
    // Just a wrapper for the deadweight body load computations:
    PointLoad *p = dynamic_cast< PointLoad * >(load);
    if ( p ) {
        FloatArray lcoords;
        if ( this->computeLocalCoordinates(lcoords, p->giveCoordinates()) ) {
            this->computePointLoadVectorAt(answer, load, tStep, mode);
        }
    } else {
        ///@todo This assumption of dead-weight loads needs to be lifted. We can have other body loads, such as
        this->computeBodyLoadVectorAt(answer, load, tStep, mode);
    }
}

NURBSBeam3dElementDsg :: NURBSBeam3dElementDsg(int n, Domain *aDomain) : IGAElement(n, aDomain), Beam3dElementEvaluatorDsg(), interpolation(3, 1) { }
  // interpolation(nsd, fsd)

IRResultType NURBSBeam3dElementDsg :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                    // Required by IR_GIVE_FIELD macro

    this->zaxis.clear();
    if ( ir->hasField(_IFT_NURBSBeam3dElementDsg_zaxis) ) {
        IR_GIVE_FIELD(ir, this->zaxis, _IFT_NURBSBeam3dElementDsg_zaxis);
    } /*else if ( ir->hasField(_IFT_Beam3d_refnode) ) {
        IR_GIVE_FIELD(ir, referenceNode, _IFT_Beam3d_refnode);
        if ( referenceNode == 0 ) {
            OOFEM_WARNING("wrong reference node specified. Using default orientation.");
        }
    } else if ( ir->hasField(_IFT_Beam3d_refangle) ) {
        IR_GIVE_FIELD(ir, referenceAngle, _IFT_Beam3d_refangle);
    } else {
        OOFEM_WARNING("y-axis, reference node or angle not set");
        return IRRT_NOTFOUND;
	}*/
    return IGAElement :: initializeFrom(ir);
}


int NURBSBeam3dElementDsg :: checkConsistency()
{
    NURBSInterpolationLine2d *interpol = static_cast< NURBSInterpolationLine2d * >( this->giveInterpolation() );
    if ( giveNumberOfDofManagers() != interpol->giveNumberOfControlPoints(1) ) {
        OOFEM_WARNING("number of control points mismatch");
        return 0;
    }
    return 1;
}

void NURBSBeam3dElementDsg :: computeLoadVector(FloatArray &answer, BodyLoad *load, CharType type, ValueModeType mode, TimeStep *tStep)
{
    if ( type != ExternalForcesVector ) {
        answer.clear();
        return;
    }
    // Just a wrapper for the deadweight body load computations:
    PointLoad *p = dynamic_cast< PointLoad * >(load);
    if ( p ) {
        FloatArray lcoords;
        if ( this->computeLocalCoordinates(lcoords, p->giveCoordinates()) ) {
            this->computePointLoadVectorAt(answer, load, tStep, mode);
        }
    } else {
        ///@todo This assumption of dead-weight loads needs to be lifted. We can have other body loads, such as
        this->computeBodyLoadVectorAt(answer, load, tStep, mode);
    }
}


  NURBSBeam2dElement :: NURBSBeam2dElement(int n, Domain *aDomain) : IGAElement(n, aDomain), Beam2dElementEvaluator(), interpolation(3, 1) { }
  // interpolation(nsd, fsd)

IRResultType NURBSBeam2dElement :: initializeFrom(InputRecord *ir)
{
    return IGAElement :: initializeFrom(ir);
}


int NURBSBeam2dElement :: checkConsistency()
{
    NURBSInterpolationLine2d *interpol = static_cast< NURBSInterpolationLine2d * >( this->giveInterpolation() );
    if ( giveNumberOfDofManagers() != interpol->giveNumberOfControlPoints(1) ) {
        OOFEM_WARNING("number of control points mismatch");
        return 0;
    }
    return 1;
}

void NURBSBeam2dElement :: computeLoadVector(FloatArray &answer, BodyLoad *load, CharType type, ValueModeType mode, TimeStep *tStep)
{
    if ( type != ExternalForcesVector ) {
        answer.clear();
        return;
    }
    // Just a wrapper for the deadweight body load computations:
    PointLoad *p = dynamic_cast< PointLoad * >(load);
    if ( p ) {
        FloatArray lcoords;
        if ( this->computeLocalCoordinates(lcoords, p->giveCoordinates()) ) {
            this->computePointLoadVectorAt(answer, load, tStep, mode);
        }
    } else {
        ///@todo This assumption of dead-weight loads needs to be lifted. We can have other body loads, such as
        this->computeBodyLoadVectorAt(answer, load, tStep, mode);
    }
}


  NURBSBeam3dElement :: NURBSBeam3dElement(int n, Domain *aDomain) : IGAElement(n, aDomain), Beam3dElementEvaluator(), interpolation(3, 1) { }
  // interpolation(nsd, fsd)

IRResultType NURBSBeam3dElement :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                    // Required by IR_GIVE_FIELD macro

    this->integrationType = 0; 
    if ( ir->hasField(_IFT_NURBSBeam3dElement_integrationType) ) {
        IR_GIVE_FIELD(ir, this->integrationType, _IFT_NURBSBeam3dElement_integrationType);
    }

    this->zaxis.clear();
    if ( ir->hasField(_IFT_NURBSBeam3dElement_zaxis) ) {
        IR_GIVE_FIELD(ir, this->zaxis, _IFT_NURBSBeam3dElement_zaxis);
    } /*else if ( ir->hasField(_IFT_Beam3d_refnode) ) {
        IR_GIVE_FIELD(ir, referenceNode, _IFT_Beam3d_refnode);
        if ( referenceNode == 0 ) {
            OOFEM_WARNING("wrong reference node specified. Using default orientation.");
        }
    } else if ( ir->hasField(_IFT_Beam3d_refangle) ) {
        IR_GIVE_FIELD(ir, referenceAngle, _IFT_Beam3d_refangle);
    } else {
        OOFEM_WARNING("y-axis, reference node or angle not set");
        return IRRT_NOTFOUND;
	}*/
    result = IGAElement :: initializeFrom(ir);

    int indx = 0;
    //numberOfGaussPoints = 1;
    double du;
    FloatArray newgpcoords;
    IntArray knotSpan;

    if (integrationType == 1) {
        int numberOfKnotSpansU = this->giveInterpolation()->giveNumberOfKnotSpans(1);
        const IntArray *knotMultiplicityU = this->giveInterpolation()->giveKnotMultiplicity(1);
        const FloatArray *knotValuesU = this->giveInterpolation()->giveKnotValues(1);

        newgpcoords.resize(1);
        knotSpan.resize(1);

	IntArray reducedNumberOfGaussPoints;
	this->giveRINumberOfGaussPoints(reducedNumberOfGaussPoints, numberOfKnotSpansU);
	//	int reducedNumberOfGaussPoints = numberOfGaussPoints-2;
        int numberOfIntegrationRules = numberOfKnotSpansU;
        reducedIntegrationRulesArray.resize( numberOfIntegrationRules );
	
	knotSpan.at(1) = -1;
	for ( int ui = 1; ui <= numberOfKnotSpansU; ui++ ) {
	  du = knotValuesU->at(ui + 1) - knotValuesU->at(ui);
	  knotSpan.at(1) += knotMultiplicityU->at(ui);
	  
	  reducedIntegrationRulesArray [ indx ].reset( new IGAIntegrationElement(indx, this, knotSpan) );
	  reducedIntegrationRulesArray [ indx ]->SetUpPointsOnLine(reducedNumberOfGaussPoints(indx), _2dBeam); // HUHU 

	  // remap local subelement gp coordinates into knot span coordinates and update integration weight
	  for ( GaussPoint *gp: *reducedIntegrationRulesArray [ indx ] ) {
	    const FloatArray &gpcoords = gp->giveNaturalCoordinates();
	    
	    newgpcoords.at(1) = knotValuesU->at(ui) + du * ( gpcoords.at(1) / 2.0 + 0.5 );
	    gp->setNaturalCoordinates(newgpcoords);
	    gp->setWeight(gp->giveWeight() / 2.0 * du );
	  }
	  indx++;
        }
    }

    return result;
}


int NURBSBeam3dElement :: checkConsistency()
{
    NURBSInterpolationLine3d *interpol = static_cast< NURBSInterpolationLine3d * >( this->giveInterpolation() );
    if ( giveNumberOfDofManagers() != interpol->giveNumberOfControlPoints(1) ) {
        OOFEM_WARNING("number of control points mismatch");
        return 0;
    }
    return 1;
}

void NURBSBeam3dElement :: giveRINumberOfGaussPoints(IntArray &answer, int nkns) {
    answer.resize(nkns);
    answer.zero();

    NURBSInterpolationLine3d *interpol = static_cast< NURBSInterpolationLine3d * >( this->giveInterpolation() );
    int degree = interpol->giveDegree();
    if (degree == 80){
	
    } else {
	for (int i = 1; i<=nkns; i++){
	    answer.at(i) = numberOfGaussPoints-2;
	}
    }
}


void NURBSBeam3dElement :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if (integrationType == 0) {
	if ( integrationRulesArray.size() == 0 ) {
	    // the gauss point is used only when methods from crosssection and/or material
	    // classes are requested
	    integrationRulesArray.resize(1);
	    integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 3) );
	    this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], this->numberOfGaussPoints, this);
	}
    } else if (integrationType == 1) {
	//	OOFEM_ERROR("Reduced integration not supported yet");
	if ( integrationRulesArray.size() == 0 ) {
	    // the gauss point is used only when methods from crosssection and/or material
	    // classes are requested
	    integrationRulesArray.resize(1);
	    integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 3) );
	    this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], this->numberOfGaussPoints-7, this);
	}
    } else {
	OOFEM_ERROR("IntegrationType %d not supported", integrationType);
    }
}

void NURBSBeam3dElement :: computeLoadVector(FloatArray &answer, BodyLoad *load, CharType type, ValueModeType mode, TimeStep *tStep)
{
    if ( type != ExternalForcesVector ) {
        answer.clear();
        return;
    }
    // Just a wrapper for the deadweight body load computations:
    PointLoad *p = dynamic_cast< PointLoad * >(load);
    if ( p ) {
        FloatArray lcoords;
        if ( this->computeLocalCoordinates(lcoords, p->giveCoordinates()) ) {
            this->computePointLoadVectorAt(answer, load, tStep, mode);
        }
    } else {
        ///@todo This assumption of dead-weight loads needs to be lifted. We can have other body loads, such as
        this->computeBodyLoadVectorAt(answer, load, tStep, mode);
    }
}

// HUHU should be implemented by IGA element (it is the same for Bspline NURBS and TSpline)
// however in such a case it should be generalized in terms of appropriately multiplying
// nseq for those integration elements which span more than just a single knot span
// the reason is to ensure compatible division to quads over which scalar quantity is interpolated
// bilinearly !!!

#ifdef __OOFEG

//#define COMPUTE_STRESS
 #define COMPUTE_STRAIN

void BsplinePlaneStressElement :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
    int indx;
    WCRec p [ 4 ];
    GraphicObj *go;
    double s [ 4 ];
    FEInterpolation *interp = this->giveInterpolation();
    FloatArray val;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    EASValsSetFillStyle(FILL_SOLID);
    const double *const *knotVector = interp->giveKnotVector();
    const IntArray *span;
    int j, nsd = this->giveNsd();
    FloatArray c [ 4 ], cg [ 4 ], u;
    IntArray sign [ 4 ];

    if ( nsd == 2 ) {
        for ( j = 0; j < 4; j++ ) {
            c [ j ].resize(2);
            cg [ j ].resize(2);
            sign [ j ].resize(2);
        }
        sign [ 0 ].at(1) = 1;
        sign [ 0 ].at(2) = 1;
        sign [ 1 ].at(1) = -1;
        sign [ 1 ].at(2) = 1;
        sign [ 2 ].at(1) = -1;
        sign [ 2 ].at(2) = -1;
        sign [ 3 ].at(1) = 1;
        sign [ 3 ].at(2) = -1;
    } else {
        OOFEM_ERROR("not implemented for nsd = %d", nsd);
    }

    indx = gc.giveIntVarIndx();

    StructuralElementEvaluator :: computeVectorOf(VM_Total, tStep, u);

    // loop over individual integration rules (i.e., knot spans)
    for ( auto &iRule: integrationRulesArray ) {
        span = iRule->giveKnotSpan();
        if ( nsd == 2 ) {
            // divide span locally to get finer geometry rep.
            int i, j, k, nseg = 4;
            double du = ( knotVector [ 0 ] [ span->at(1) + 1 ] - knotVector [ 0 ] [ span->at(1) ] ) / nseg;
            double dv = ( knotVector [ 1 ] [ span->at(2) + 1 ] - knotVector [ 1 ] [ span->at(2) ] ) / nseg;
            for ( i = 1; i <= nseg; i++ ) {
                for ( j = 1; j <= nseg; j++ ) {
                    c [ 0 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                    c [ 0 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                    c [ 1 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                    c [ 1 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                    c [ 2 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                    c [ 2 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;
                    c [ 3 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                    c [ 3 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;

                    for ( k = 0; k < 4; k++ ) {
                        interp->local2global( cg [ k ], c [ k ], FEIIGAElementGeometryWrapper( this, iRule->giveKnotSpan() ) );
                        p [ k ].x = ( FPNum ) cg [ k ].at(1);
                        p [ k ].y = ( FPNum ) cg [ k ].at(2);
                        p [ k ].z = 0.;

 #ifdef QUARTER_PLATE_WITH_HOLE_SINGLE_PATCH
                        //move sampling gp out of boundary to overcome degeneracy on quarter plate with hole modelled by single patch
                        if ( c [ k ].at(1) > 0.99999 && c [ k ].at(2) > 0.495 && c [ k ].at(2) < 0.505 ) {
                            c [ k ].at(1) += sign [ k ].at(1) * du / 10.0;
                            c [ k ].at(2) += sign [ k ].at(2) * dv / 10.0;
                            c [ k ].at(3) += sign [ k ].at(3) * dw / 10.0;
                        }
 #endif

                        // create a dummy ip's
                        GaussPoint gp(iRule.get(), 999, c [ k ], 1.0, _PlaneStress);
 #ifdef COMPUTE_STRAIN
                        this->computeStrainVector(val, & gp, tStep, u);
 #endif
 #ifdef COMPUTE_STRESS
                        FloatArray strain;
                        this->computeStrainVector(strain, & gp, tStep, u);
                        this->computeStressVector(val, strain, & gp, tStep);
 #endif
                        s [ k ] = val.at(indx);
                    }

                    if ( ( std::isnan(s [ 0 ]) ) || ( std::isnan(s [ 1 ]) ) || ( std::isnan(s [ 2 ]) ) || ( std::isnan(s [ 3 ]) ) ) {
                        continue;
                    }

                    if ( ( fabs(s [ 0 ]) > 1.e5 ) || ( fabs(s [ 1 ]) > 1.e5 ) || ( fabs(s [ 2 ]) > 1.e5 ) || ( fabs(s [ 3 ]) > 1.e5 ) ) {
                        continue;
                    }

                    //printf ("QWD: %e %e %e %e\n", s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);

                    go =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
                    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
                    EGAttachObject(go, ( EObjectP ) this);
                    EMAddGraphicsToModel(ESIModel(), go);
                }
            }
        }
    } // end loop over knot spans (irules)
}


void NURBSPlaneStressElement :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
    int indx;
    WCRec p [ 4 ];
    GraphicObj *go;
    double s [ 4 ];
    FEInterpolation *interp = this->giveInterpolation();
    FloatArray val, u;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    EASValsSetFillStyle(FILL_SOLID);
    EASValsSetEdgeFlag(true);
    const double *const *knotVector = interp->giveKnotVector();
    const IntArray *span;
    int j, nsd = this->giveNsd();
    FloatArray c [ 4 ], cg [ 4 ];
    IntArray sign [ 4 ];

    if ( nsd == 2 ) {
        for ( j = 0; j < 4; j++ ) {
            c [ j ].resize(2);
            cg [ j ].resize(2);
            sign [ j ].resize(2);
        }
        sign [ 0 ].at(1) = 1;
        sign [ 0 ].at(2) = 1;
        sign [ 1 ].at(1) = -1;
        sign [ 1 ].at(2) = 1;
        sign [ 2 ].at(1) = -1;
        sign [ 2 ].at(2) = -1;
        sign [ 3 ].at(1) = 1;
        sign [ 3 ].at(2) = -1;
    } else {
        OOFEM_ERROR("not implemented for nsd = %d", nsd);
    }

    indx = gc.giveIntVarIndx();

    StructuralElementEvaluator :: computeVectorOf(VM_Total, tStep, u);

    //double maxs=-1.0e10, mins=1.0e10;

    // loop over individual integration rules (i.e., knot spans)
    for ( auto &iRule: integrationRulesArray ) {
        span = iRule->giveKnotSpan();
        if ( nsd == 2 ) {
            // divide span locally to get finer geometry rep.
            int i, j, k, nseg = 8;
            double du = ( knotVector [ 0 ] [ span->at(1) + 1 ] - knotVector [ 0 ] [ span->at(1) ] ) / nseg;
            double dv = ( knotVector [ 1 ] [ span->at(2) + 1 ] - knotVector [ 1 ] [ span->at(2) ] ) / nseg;
            for ( i = 1; i <= nseg; i++ ) {
                for ( j = 1; j <= nseg; j++ ) {
                    c [ 0 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                    c [ 0 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                    c [ 1 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                    c [ 1 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                    c [ 2 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                    c [ 2 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;
                    c [ 3 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                    c [ 3 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;

                    for ( k = 0; k < 4; k++ ) {
                        interp->local2global( cg [ k ], c [ k ], FEIIGAElementGeometryWrapper( this, iRule->giveKnotSpan() ) );
                        p [ k ].x = ( FPNum ) cg [ k ].at(1);
                        p [ k ].y = ( FPNum ) cg [ k ].at(2);
                        p [ k ].z = 0.;

                        //move sampling gp out of boundary to overcome degeneracy on quarter plate with hole modelled by single patch
                        if ( c [ k ].at(1) > 0.99999 && c [ k ].at(2) > 0.495 && c [ k ].at(2) < 0.505 ) {
                            c [ k ].at(1) += sign [ k ].at(1) * du / 10.0;
                            c [ k ].at(2) += sign [ k ].at(2) * dv / 10.0;
                        }

                        // create a dummy ip's
                        GaussPoint gp(iRule.get(), 999, c [ k ], 1.0, _PlaneStress);
 #ifdef COMPUTE_STRAIN
                        this->computeStrainVector(val, & gp, tStep, u);
 #endif
 #ifdef COMPUTE_STRESS
                        FloatArray strain;
                        this->computeStrainVector(strain, & gp, tStep, u);
                        this->computeStressVector(val, strain, & gp, tStep);
 #endif
                        s [ k ] = val.at(indx);

 #ifdef QUARTER_PLATE_WITH_HOLE
                        if ( indx == 2 ) {
                            if ( cg [ k ].at(1) <= 1.0 + 1.0e-10 && cg [ k ].at(2) <= 1.0e-10 ) {
                                fprintf(stderr, "A: syy = %e\n", s [ k ]);
                            }
                        }
                        if ( indx == 1 ) {
                            if ( cg [ k ].at(1) <= 1.0e-10 && cg [ k ].at(2) <= 1.0 + 1.0e-10 ) {
                                fprintf(stderr, "B: sxx = %e\n", s [ k ]);
                            }
                        }

                        double x, y, r, phi, rate, E, G, kap, ny;
                        x = cg [ k ].at(1);
                        y = cg [ k ].at(2);
                        if ( x < 1.0e-10 ) {
                            phi = M_PI / 2.0;
                            r = y;
                        } else {
                            phi = atan(y / x);
                            r = x / cos(phi);
                        }

  #ifdef STRESS
                        // exact stresses quarter plate with hole s0=1 a=1
                        rate = 1.0 / r;
                        rate *= rate;
                        if ( indx == 1 ) {
                            s [ k ] = 0.5 * ( 2.0 + 3.0 * rate * rate * cos(4.0 * phi) - rate * ( 3 * cos(2.0 * phi) + 2.0 * cos(4.0 * phi) ) );
                        }
                        if ( indx == 2 ) {
                            s [ k ] = 0.5 * ( -3.0 * rate * rate * cos(4.0 * phi) - rate * ( cos(2.0 * phi) - 2.0 * cos(4.0 * phi) ) );
                        }
                        if ( indx == 3 ) {
                            s [ k ] = 0.5 * ( 3.0 * rate * rate * sin(4.0 * phi) - rate * ( sin(2.0 * phi) + 2.0 * sin(4.0 * phi) ) );
                        }

                        if ( indx == 2 ) {
                            if ( cg [ k ].at(1) <= 1.0 + 1.0e-10 && cg [ k ].at(2) <= 1.0e-10 ) {
                                fprintf(stderr, "A: syy = %e\n", s [ k ]);
                            }
                        }
                        if ( indx == 1 ) {
                            if ( cg [ k ].at(1) <= 1.0e-10 && cg [ k ].at(2) <= 1.0 + 1.0e-10 ) {
                                fprintf(stderr, "B: sxx = %e\n", s [ k ]);
                            }
                        }
  #endif
  #ifdef DISPL
                        // exact displ quarter plate with hole s0=1 a=1
                        E = 15000.0;
                        ny = 0.25;
                        G = E / ( 2.0 * ( 1.0 + ny ) );
                        kap = ( 3.0 - ny ) / ( 1.0 + ny );
                        rate = 1.0 / r;
                        if ( indx == 1 ) {
                            s [ k ] = ( r * ( kap + 1.0 ) * cos(phi) + 2.0 * rate * ( ( 1.0 + kap ) * cos(phi) + cos(3.0 * phi) ) - 2.0 * rate * rate * rate * cos(3.0 * phi) ) / ( 8.0 * G );
                        }
                        if ( indx == 2 ) {
                            s [ k ] = ( r * ( kap - 3.0 ) * sin(phi) + 2.0 * rate * ( ( 1.0 - kap ) * sin(phi) + sin(3.0 * phi) ) - 2.0 * rate * rate * rate * sin(3.0 * phi) ) / ( 8.0 * G );
                        }

                        if ( indx == 1 ) {
                            if ( cg [ k ].at(1) <= 1.0 + 1.0e-10 && cg [ k ].at(2) <= 1.0e-10 ) {
                                fprintf(stderr, "A: ux = %e\n", s [ k ]);
                            }
                        }
                        if ( indx == 2 ) {
                            if ( cg [ k ].at(1) <= 1.0e-10 && cg [ k ].at(2) <= 1.0 + 1.0e-10 ) {
                                fprintf(stderr, "B: uy = %e\n", s [ k ]);
                            }
                        }
  #endif

                        if ( s [ k ] < mins ) {
                            mins = s [ k ];
                        }
                        if ( s [ k ] > maxs ) {
                            maxs = s [ k ];
                        }
 #endif
                    }

                    if ( ( std::isnan(s [ 0 ]) ) || ( std::isnan(s [ 1 ]) ) || ( std::isnan(s [ 2 ]) ) || ( std::isnan(s [ 3 ]) ) ) {
                        continue;
                    }

                    if ( ( fabs(s [ 0 ]) > 1.e5 ) || ( fabs(s [ 1 ]) > 1.e5 ) || ( fabs(s [ 2 ]) > 1.e5 ) || ( fabs(s [ 3 ]) > 1.e5 ) ) {
                        continue;
                    }

                    //printf ("QWD: %e %e %e %e\n", s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);

                    go =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
                    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
                    EGAttachObject(go, ( EObjectP ) this);
                    EMAddGraphicsToModel(ESIModel(), go);
                }
            }
        }
    } // end loop over knot spans (irules)

 #ifdef QUARTER_PLATE_WITH_HOLE
    fprintf(stderr, "%d: %e %e %e %e\n", indx, mins, maxs, ( 10.0 * mins + maxs ) / 11.0, ( 10.0 * maxs + mins ) / 11.0);
 #endif
}


// refinement of integration elements should be generalized in terms of appropriately multiplying
// nseq for those integration elements which span more tham just a single knot span
// the reason is to ensure compatible division to quads over which scalar quantity is interpolated
// bilinearly !!!

void TSplinePlaneStressElement :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
    int indx;
    WCRec p [ 4 ];
    GraphicObj *go;
    double s [ 4 ];
    FEInterpolation *interp = this->giveInterpolation();
    FloatArray val;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    EASValsSetFillStyle(FILL_SOLID);
    EASValsSetEdgeFlag(true);
    const double *const *knotVector = interp->giveKnotVector();
    const IntArray *span;
    int j, nsd = this->giveNsd();
    FloatArray c [ 4 ], cg [ 4 ], u;
    IntArray sign [ 4 ];

    if ( nsd == 2 ) {
        for ( j = 0; j < 4; j++ ) {
            c [ j ].resize(2);
            cg [ j ].resize(2);
            sign [ j ].resize(2);
        }
        sign [ 0 ].at(1) = 1;
        sign [ 0 ].at(2) = 1;
        sign [ 1 ].at(1) = -1;
        sign [ 1 ].at(2) = 1;
        sign [ 2 ].at(1) = -1;
        sign [ 2 ].at(2) = -1;
        sign [ 3 ].at(1) = 1;
        sign [ 3 ].at(2) = -1;
    } else {
        OOFEM_ERROR("not implemented for nsd = %d", nsd);
    }

    indx = gc.giveIntVarIndx();

    StructuralElementEvaluator :: computeVectorOf(VM_Total, tStep, u);

    // loop over individual integration rules (i.e., knot spans)
    for ( auto &iRule: integrationRulesArray ) {
        span = iRule->giveKnotSpan();
        if ( nsd == 2 ) {
            // divide span locally to get finer geometry rep.
            int i, j, k, nseg = 4;
            double du = ( knotVector [ 0 ] [ span->at(1) + 1 ] - knotVector [ 0 ] [ span->at(1) ] ) / nseg;
            double dv = ( knotVector [ 1 ] [ span->at(2) + 1 ] - knotVector [ 1 ] [ span->at(2) ] ) / nseg;
            for ( i = 1; i <= nseg; i++ ) {
                for ( j = 1; j <= nseg; j++ ) {
                    c [ 0 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                    c [ 0 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                    c [ 1 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                    c [ 1 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                    c [ 2 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                    c [ 2 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;
                    c [ 3 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                    c [ 3 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;

                    for ( k = 0; k < 4; k++ ) {
                        interp->local2global( cg [ k ], c [ k ], FEIIGAElementGeometryWrapper( this, iRule->giveKnotSpan() ) );
                        p [ k ].x = ( FPNum ) cg [ k ].at(1);
                        p [ k ].y = ( FPNum ) cg [ k ].at(2);
                        p [ k ].z = 0.;

 #ifdef QUARTER_PLATE_WITH_HOLE_SINGLE_PATCH
                        //move sampling gp out of boundary to overcome degeneracy on quarter plate with hole modelled by single patch
                        if ( c [ k ].at(1) > 0.99999 && c [ k ].at(2) > 0.495 && c [ k ].at(2) < 0.505 ) {
                            c [ k ].at(1) += sign [ k ].at(1) * du / 10.0;
                            c [ k ].at(2) += sign [ k ].at(2) * dv / 10.0;
                        }
 #endif

                        // create a dummy ip's
                        GaussPoint gp(iRule.get(), 999, c [ k ], 1.0, _PlaneStress);
 #ifdef COMPUTE_STRAIN
                        this->computeStrainVector(val, & gp, tStep, u);
 #endif
 #ifdef COMPUTE_STRESS
                        FloatArray strain;
                        this->computeStrainVector(strain, & gp, tStep, u);
                        this->computeStressVector(val, strain, & gp, tStep);
 #endif
                        s [ k ] = val.at(indx);
                    }

                    if ( ( std::isnan(s [ 0 ]) ) || ( std::isnan(s [ 1 ]) ) || ( std::isnan(s [ 2 ]) ) || ( std::isnan(s [ 3 ]) ) ) {
                        continue;
                    }

                    if ( ( fabs(s [ 0 ]) > 1.e5 ) || ( fabs(s [ 1 ]) > 1.e5 ) || ( fabs(s [ 2 ]) > 1.e5 ) || ( fabs(s [ 3 ]) > 1.e5 ) ) {
                        continue;
                    }

                    //printf ("QWD: %e %e %e %e\n", s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);

                    go =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
                    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
                    EGAttachObject(go, ( EObjectP ) this);
                    EMAddGraphicsToModel(ESIModel(), go);
                }
            }
        }
    } // end loop over knot spans (irules)
}



void NURBSSpace3dElement :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
    int indx;
    WCRec p [ 8 ];
    GraphicObj *go;
    double s [ 8 ];
    FEInterpolation *interp = this->giveInterpolation();
    FloatArray val, u;
    //int huhu = 0;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    EASValsSetFillStyle(FILL_SOLID);
    EASValsSetEdgeFlag(true);
    const double *const *knotVector = interp->giveKnotVector();
    const IntArray *span;
    int j, nsd = this->giveNsd();
    FloatArray c [ 8 ], cg [ 8 ];
    IntArray sign [ 8 ];

    if ( nsd == 3 ) {
        for ( j = 0; j < 8; j++ ) {
            c [ j ].resize(3);
            cg [ j ].resize(3);
            sign [ j ].resize(3);
        }

        sign [ 0 ].at(1) = 1;
        sign [ 0 ].at(2) = 1;
        sign [ 0 ].at(3) = 1;
        sign [ 1 ].at(1) = -1;
        sign [ 1 ].at(2) = 1;
        sign [ 1 ].at(3) = 1;
        sign [ 2 ].at(1) = -1;
        sign [ 2 ].at(2) = -1;
        sign [ 2 ].at(3) = 1;
        sign [ 3 ].at(1) = 1;
        sign [ 3 ].at(2) = -1;
        sign [ 3 ].at(3) = 1;
        sign [ 4 ].at(1) = 1;
        sign [ 4 ].at(2) = 1;
        sign [ 4 ].at(3) = -1;
        sign [ 5 ].at(1) = -1;
        sign [ 5 ].at(2) = 1;
        sign [ 5 ].at(3) = -1;
        sign [ 6 ].at(1) = -1;
        sign [ 6 ].at(2) = -1;
        sign [ 6 ].at(3) = -1;
        sign [ 7 ].at(1) = 1;
        sign [ 7 ].at(2) = -1;
        sign [ 7 ].at(3) = -1;
    } else {
        OOFEM_ERROR("not implemented for nsd = %d", nsd);
    }

    indx = gc.giveIntVarIndx();

 #ifdef SPHERICAL_CS
    if ( gc.giveIntVarType() == IST_StrainTensor ) {
        huhu = 1;
    }
 #endif
 #ifdef MISSES_STRESS
    if ( gc.giveIntVarType() == IST_ErrorIndicatorLevel ) {
        huhu = 2;
    }
 #endif

    StructuralElementEvaluator :: computeVectorOf(VM_Total, tStep, u);

    //double maxs=-1.0e10, mins=1.0e10;

    // loop over individual integration rules (i.e., knot spans)
    for ( auto &iRule: integrationRulesArray ) {
        span = iRule->giveKnotSpan();
        if ( nsd == 3 ) {
            // divide span locally to get finer geometry rep.
            int i, j, k, m, nseg = 8;
            double du = ( knotVector [ 0 ] [ span->at(1) + 1 ] - knotVector [ 0 ] [ span->at(1) ] ) / nseg;
            double dv = ( knotVector [ 1 ] [ span->at(2) + 1 ] - knotVector [ 1 ] [ span->at(2) ] ) / nseg;
            double dw = ( knotVector [ 2 ] [ span->at(3) + 1 ] - knotVector [ 2 ] [ span->at(3) ] ) / nseg;

 #ifdef DRAW_VISIBLE_CONTOUR
            WCRec pp [ 4 ];
            double ss [ 4 ];
            int kk, ii, nd [ 6 ] [ 4 ] = { { 0, 1, 2, 3 }, { 4, 5, 6, 7 }, { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 3, 0, 4, 7 } };
 #endif

            for ( i = 1; i <= nseg; i++ ) {
                for ( j = 1; j <= nseg; j++ ) {
                    for ( m = 1; m <= nseg; m++ ) {
                        c [ 0 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                        c [ 0 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                        c [ 0 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dw * ( m - 1 );
                        c [ 1 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                        c [ 1 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                        c [ 1 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dw * ( m - 1 );
                        c [ 2 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                        c [ 2 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;
                        c [ 2 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dw * ( m - 1 );
                        c [ 3 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                        c [ 3 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;
                        c [ 3 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dw * ( m - 1 );
                        c [ 4 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                        c [ 4 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                        c [ 4 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dw * m;
                        c [ 5 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                        c [ 5 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                        c [ 5 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dw * m;
                        c [ 6 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                        c [ 6 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;
                        c [ 6 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dw * m;
                        c [ 7 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                        c [ 7 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;
                        c [ 7 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dw * m;

                        for ( k = 0; k < 8; k++ ) {
                            interp->local2global( cg [ k ], c [ k ], FEIIGAElementGeometryWrapper( this, iRule->giveKnotSpan() ) );
                            p [ k ].x = ( FPNum ) cg [ k ].at(1);
                            p [ k ].y = ( FPNum ) cg [ k ].at(2);
                            p [ k ].z = ( FPNum ) cg [ k ].at(3);
                        }

 #ifdef DRAW_VISIBLE_CONTOUR
  #ifdef SPHERE_WITH_HOLE

                        // check whether drawing side visible in particular view !!!
                        int haha = 0;
                        for ( kk = 0; kk < 6; kk++ ) {
                            double xx = 0.0, yy = 0.0, zz = 0.0, rr, r;
                            for ( ii = 0; ii < 4; ii++ ) {
                                xx += ( pp [ ii ].x = p [ nd [ kk ] [ ii ] ].x );
                                yy += ( pp [ ii ].y = p [ nd [ kk ] [ ii ] ].y );
                                zz += ( pp [ ii ].z = p [ nd [ kk ] [ ii ] ].z );
                                ss [ ii ] = s [ nd [ kk ] [ ii ] ];
                            }
                            xx /= 4.0;
                            yy /= 4.0;
                            zz /= 4.0;
                            rr = xx * xx + yy * yy;
                            r = rr + zz * zz;

                            if ( zz < 2.0001 /* || xx < 0.0001 */ || yy < 0.0001 || rr < 1.0 || r < 25.0 || r > 5.99 * 5.99 ) {
                                haha = 1;
                            }
                        }
                        if ( haha == 0 ) {
                            continue;
                        }
  #endif
 #endif

                        for ( k = 0; k < 8; k++ ) {
                            // create a dummy ip's
                            GaussPoint gp(iRule.get(), 999, c [ k ], 1.0, _3dMat);
 #ifdef COMPUTE_STRAIN
                            this->computeStrainVector(val, & gp, tStep, u);
 #endif
 #ifdef COMPUTE_STRESS
                            FloatArray strain;
                            this->computeStrainVector(strain, & gp, tStep, u);
                            this->computeStressVector(val, strain, & gp, tStep);
 #endif
                            s [ k ] = val.at(indx);

 #ifdef MISSES_STRESS
                            if ( huhu == 2 ) {
                                double vonMisses;
                                vonMisses = sqrt( ( ( val.at(1) - val.at(2) ) * ( val.at(1) - val.at(2) ) + ( val.at(2) - val.at(3) ) * ( val.at(2) - val.at(3) ) + ( val.at(1) - val.at(3) ) * ( val.at(1) - val.at(3) ) + 6.0 * ( val.at(4) * val.at(4) + val.at(5) * val.at(5) + val.at(6) * val.at(6) ) ) / 2.0 );
                                s [ k ] = vonMisses;
                            }
 #endif
 #ifdef SPHERICAL_CS
                            if ( huhu == 1 ) {
                                double x, y, z, r, r2, rr;
                                FloatMatrix sigma(3, 3), T(3, 3), product(3, 3), result(3, 3);

                                sigma.at(1, 1) = val.at(1);
                                sigma.at(1, 2) = val.at(6);
                                sigma.at(1, 3) = val.at(5);
                                sigma.at(2, 1) = val.at(6);
                                sigma.at(2, 2) = val.at(2);
                                sigma.at(2, 3) = val.at(4);
                                sigma.at(3, 1) = val.at(5);
                                sigma.at(3, 2) = val.at(4);
                                sigma.at(3, 3) = val.at(3);

                                x = cg [ k ].at(1);
                                y = cg [ k ].at(2);
                                z = cg [ k ].at(3);
                                r2 = x * x + y * y + z * z;
                                r = sqrt(r2);
                                rr = sqrt(x * x + y * y);

                                T.at(1, 1) = -z * z * y / rr / r2 - y * rr / r2;
                                T.at(1, 2) = z * z * x / rr / r2 + x * rr / r2;
                                T.at(1, 3) = 0.0;
                                T.at(2, 1) = -z * x / rr / r;
                                T.at(2, 2) = -z * y / rr / r;
                                T.at(2, 3) = rr / r;
                                T.at(3, 1) = x / r;
                                T.at(3, 2) = y / r;
                                T.at(3, 3) = z / r;

                                product.beProductOf(T, sigma);
                                result.beProductTOf(product, T);

                                val.at(1) = result.at(1, 1);
                                val.at(6) = result.at(1, 2);
                                val.at(5) = result.at(1, 3);
                                val.at(6) = result.at(2, 1);
                                val.at(2) = result.at(2, 2);
                                val.at(4) = result.at(2, 3);
                                val.at(5) = result.at(3, 1);
                                val.at(4) = result.at(3, 2);
                                val.at(3) = result.at(3, 3);

                                s [ k ] = val.at(indx);
                            }
 #endif
 #ifdef SPHERE_WITH_HOLE
                            if ( s [ k ] < mins ) {
                                mins = s [ k ];
                            }
                            if ( s [ k ] > maxs ) {
                                maxs = s [ k ];
                            }
 #endif
                        }

                        if ( ( std::isnan(s [ 0 ]) ) || ( std::isnan(s [ 1 ]) ) || ( std::isnan(s [ 2 ]) ) || ( std::isnan(s [ 3 ]) ) ) {
                            continue;
                        }

                        if ( ( std::isnan(s [ 4 ]) ) || ( std::isnan(s [ 5 ]) ) || ( std::isnan(s [ 6 ]) ) || ( std::isnan(s [ 7 ]) ) ) {
                            continue;
                        }

                        if ( ( fabs(s [ 0 ]) > 1.e5 ) || ( fabs(s [ 1 ]) > 1.e5 ) || ( fabs(s [ 2 ]) > 1.e5 ) || ( fabs(s [ 3 ]) > 1.e5 ) ) {
                            continue;
                        }

                        if ( ( fabs(s [ 4 ]) > 1.e5 ) || ( fabs(s [ 5 ]) > 1.e5 ) || ( fabs(s [ 6 ]) > 1.e5 ) || ( fabs(s [ 7 ]) > 1.e5 ) ) {
                            continue;
                        }

                        //printf ("HWD: %e %e %e %e %e %e %e %e\n", s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ], s [ 4 ], s [ 5 ], s [ 6 ], s [ 7 ]);

 #ifdef DRAW_VISIBLE_CONTOUR
  #ifdef SPHERE_WITH_HOLE
                        for ( kk = 0; kk < 6; kk++ ) {
                            double xx = 0.0, yy = 0.0, zz = 0.0, rr, r;
                            for ( ii = 0; ii < 4; ii++ ) {
                                xx += ( pp [ ii ].x = p [ nd [ kk ] [ ii ] ].x );
                                yy += ( pp [ ii ].y = p [ nd [ kk ] [ ii ] ].y );
                                zz += ( pp [ ii ].z = p [ nd [ kk ] [ ii ] ].z );
                                ss [ ii ] = s [ nd [ kk ] [ ii ] ];
                            }
                            xx /= 4.0;
                            yy /= 4.0;
                            zz /= 4.0;
                            rr = xx * xx + yy * yy;
                            r = rr + zz * zz;

                            if ( zz < 2.0001 /* || xx < 0.0001 */ || yy < 0.0001 || rr < 1.0 || r < 25.0 || r > 5.99 * 5.99 ) {
                                go = CreateQuadWD3D(pp, ss [ 0 ], ss [ 1 ], ss [ 2 ], ss [ 3 ]);
                                EGWithMaskChangeAttributes(LAYER_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK, go);
                                EMAddGraphicsToModel(ESIModel(), go);
                            }
                        }
  #endif
 #else
                        go = CreateHexahedronWD(p, s);
                        EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
                        EGAttachObject(go, ( EObjectP ) this);
                        EMAddGraphicsToModel(ESIModel(), go);
 #endif
                    }
                }
            }
        }
    } // end loop over knot spans (irules)

 #ifdef SPHERE_WITH_HOLE
    fprintf(stderr, "%d %e %e %e %e\n", indx, mins, maxs, ( 10.0 * mins + maxs ) / 11.0, ( 10.0 * maxs + mins ) / 11.0);
 #endif
}


void
NURBSBeam2dElementDsg::drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    WCRec p [ 2 ];
    GraphicObj *go;
    FEInterpolation *interp = this->giveInterpolation();

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );

    const double *const *knotVector = interp->giveKnotVector();
    const IntArray *span;
    FloatArray c [ 2 ], cg [ 2 ], u;

    c [ 0 ].resize(2);
    c [ 1 ].resize(2);
    cg[ 0 ].resize(2);
    cg[ 1 ].resize(2);

    // loop over individual integration rules (i.e., knot spans)
    for ( auto &iRule: integrationRulesArray ) {
      int i, k, nseg = 20;
        span = iRule->giveKnotSpan();
	double du = ( knotVector [ 0 ] [ span->at(1) + 1 ] - knotVector [ 0 ] [ span->at(1) ] ) / nseg;
	for ( i = 1; i <= nseg; i++ ) {
	  c [ 0 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
	  c [ 1 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
	  for ( k = 0; k < 2; k++ ) {
	    interp->local2global( cg [ k ], c [ k ], FEIIGAElementGeometryWrapper( this, iRule->giveKnotSpan() ) );
	    p [ k ].x = ( FPNum ) cg [ k ].at(1);
	    p [ k ].y = ( FPNum ) cg [ k ].at(2);
	    p [ k ].z = 0.;
	  }
	  go = CreateLine3D(p);
	  EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
	  EGAttachObject(go, ( EObjectP ) this);
	  EMAddGraphicsToModel(ESIModel(), go);
	}
    }
}


void 
NURBSBeam2dElementDsg::drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType ut) 
{

    WCRec p [ 2 ];
    GraphicObj *go;
    int i, j, k, n, nseg;
    FloatArray u;
    FloatMatrix N;
    FloatArray ur;
    IntArray lc;
    FloatArray normal;
    FEInterpolation *interp = this->giveInterpolation();
    double defScale = gc.getDefScale();

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    EASValsSetLineStyle(SOLID_STYLE);
    EASValsSetFillStyle(FILL_SOLID);
    EASValsSetLineWidth(0);

    nseg = 40;

    const double *const *knotVector = interp->giveKnotVector();
    const IntArray *span;

    this->StructuralElementEvaluator::computeVectorOf(VM_Total, tStep, u);

    FloatArray c [ 2 ], cg [ 2 ], d, dg;
    dg.resize(2);
    double du;

    for ( j = 0; j < 2; j++ ) {
      c [ j ].resize(2);
      cg [ j ].resize(2);
    }

    // loop over individual integration rules (i.e., knot spans)
    for ( int ir = 0; ir < this->giveNumberOfIntegrationRules(); ir++ ) {
      IntegrationRule *iRule = this->giveIntegrationRule(ir);
      span = iRule->giveKnotSpan();
      // divide span locally to get finer geometry rep.
      du = ( knotVector [ 0 ] [ span->at(1) + 1 ] - knotVector [ 0 ] [ span->at(1) ] ) / nseg;
      for ( i = 1; i <= nseg; i++ ) {
	c [ 0 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
	c [ 1 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;

	for ( k = 0; k < 2; k++ ) {
	  // create a dummy ip's
	  GaussPoint gp(iRule, 999, c [ k ], 1.0, _2dBeam);

	  // compute displacements at gp
	  this->computeNMatrixAt(N, & gp);

	  // get local code numbers corresponding to ir
	  this->giveIntegrationElementLocalCodeNumbers(lc, this, gp.giveIntegrationRule());
	  ur.resize( N.giveNumberOfColumns() );
	  for ( n = 1; n <= lc.giveSize(); n++ ) {
	    ur.at(n) = u.at( lc.at(n) );
	  }

	  // interpolate displacements
	  d.beProductOf(N, ur);

	  interp->local2global( cg [ k ], c [ k ], FEIIGAElementGeometryWrapper( this, iRule->giveKnotSpan() ) );

	  this->computeNormal (normal, c[k], ir);

	  dg.at(1) = d.at(1) * normal.at(2) + d.at(2) * normal.at(1);
	  dg.at(2) = -d.at(1) * normal.at(1) + d.at(2) * normal.at(2);
	  // interp->local2global( dg , d , FEIIGAElementGeometryWrapper( this, iRule->giveKnotSpan() ) );
	  // dg.beVectorProductOf(normal, d);

	  p [ k ].x = ( FPNum ) ( cg [ k ].at(1) + dg.at(1) * defScale );
	  p [ k ].y = ( FPNum ) ( cg [ k ].at(2) + dg.at(2) * defScale );
	  p [ k ].z = 0.;
	}

	go =  CreateLine3D(p);
	EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
	EGAttachObject(go, ( EObjectP ) this);
	EMAddGraphicsToModel(ESIModel(), go);
      }
    }                 // end loop over knot spans (irules)
}
  

void 
NURBSBeam2dElementDsg::drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
  WCRec p [ 2 ];
  WCRec p2 [ 2 ];
  int indx;
    GraphicObj *go;
    int i, j, k, nseg;
    FloatArray u;
    FloatMatrix N;
    FloatArray ur, d, n, val;
    IntArray lc;
    FEInterpolation *interp = this->giveInterpolation();
    // double defScale = gc.getDefScale();

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_VARPLOT_PATTERN_LAYER);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    EASValsSetLineStyle(SOLID_STYLE);
    EASValsSetFillStyle(FILL_SOLID);
    EASValsSetLineWidth(0);

    nseg = 10;

    const double *const *knotVector = interp->giveKnotVector();
    const IntArray *span;
    double scale = gc.getDefScale();


    this->StructuralElementEvaluator::computeVectorOf(VM_Total, tStep, u);

    indx = gc.giveIntVarIndx();

    FloatArray c [ 2 ], cg [ 2 ];
    double du;
    double s [ 2 ];

    for ( j = 0; j < 2; j++ ) {
      c [ j ].resize(2);
      cg [ j ].resize(2);
    }

    // loop over individual integration rules (i.e., knot spans)
    for ( int ir = 0; ir < this->giveNumberOfIntegrationRules(); ir++ ) {
      IntegrationRule *iRule = this->giveIntegrationRule(ir);
      span = iRule->giveKnotSpan();
      // divide span locally to get finer geometry rep.
      du = ( knotVector [ 0 ] [ span->at(1) + 1 ] - knotVector [ 0 ] [ span->at(1) ] ) / nseg;
      for ( i = 1; i <= nseg; i++ ) {
	c [ 0 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
	c [ 1 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;

	for ( k = 0; k < 2; k++ ) {
	  // create a dummy ip's
	  GaussPoint gp(iRule, 999, c [ k ], 1.0, _2dBeam);
	  FloatArray strain;
	  if ((gc.giveIntVarType() == IST_StrainTensor)){
	    // copy to array of size 1x6
	    this->computeStrainVector(val, & gp, tStep, u);
	  }
	    else if (gc.giveIntVarType() == IST_StressTensor) {
		// copy to array of size 1x6
		this->computeStrainVector(strain, & gp, tStep, u);
		this->computeStressVector(val, strain, & gp, tStep);
	    } else{
		OOFEM_ERROR ("Internal Variable Type not recognized");
	    }
	  s [ k ] = val.at(indx);
	  fprintf (stderr, " %e",  s [ k ]);
	  
	  interp->local2global( cg [ k ], c [ k ], FEIIGAElementGeometryWrapper( this, iRule->giveKnotSpan() ) );
	  this->computeNormal (n, c[k], ir);
	  //n = {0, 1, 0};
	  p [ k ].x = ( FPNum ) ( cg [ k ].at(1) + n.at(1) * s[k] * scale );
	  p [ k ].y = ( FPNum ) ( cg [ k ].at(2) + n.at(2) * s[k] * scale );
	  p [ k ].z = 0.;

	  p2 [ k ].x = ( FPNum ) ( cg [ 0 ].at(1) + (1-k) * n.at(1) * s[k] * scale );
	  p2 [ k ].y = ( FPNum ) ( cg [ 0 ].at(2) + (1-k) * n.at(2) * s[k] * scale );
	  p2 [ k ].z = 0.;

	}

	go =  CreateLine3D(p);
	EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
	EGAttachObject(go, ( EObjectP ) this);
	EMAddGraphicsToModel(ESIModel(), go);

	go =  CreateLine3D(p2);
	EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
	EGAttachObject(go, ( EObjectP ) this);
	EMAddGraphicsToModel(ESIModel(), go);
      }
    }                // end loop over knot spans (irules)

	    this->computeNormal (n, c[1], this->giveNumberOfIntegrationRules()-1);
	  p2 [ 0 ].x = ( FPNum ) ( cg [ 1 ].at(1) + n.at(1) * s[1] * scale );
	  p2 [ 0 ].y = ( FPNum ) ( cg [ 1 ].at(2) + n.at(2) * s[1] * scale );

	  p2 [ 1 ].x = ( FPNum ) ( cg [ 1 ].at(1) );
	  p2 [ 1 ].y = ( FPNum ) ( cg [ 1 ].at(2) );

	go =  CreateLine3D(p2);
	EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
	EGAttachObject(go, ( EObjectP ) this);
	EMAddGraphicsToModel(ESIModel(), go);
 
}


void
NURBSBeam3dElementDsg::drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    WCRec p [ 2 ];
    GraphicObj *go;
    FEInterpolation *interp = this->giveInterpolation();

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );

    const double *const *knotVector = interp->giveKnotVector();
    const IntArray *span;
    FloatArray c [ 2 ], cg [ 2 ], u;

    c [ 0 ].resize(3);
    c [ 1 ].resize(3);
    cg[ 0 ].resize(3);
    cg[ 1 ].resize(3);

    // loop over individual integration rules (i.e., knot spans)
    for ( auto &iRule: integrationRulesArray ) {
      int i, k, nseg = 20;
        span = iRule->giveKnotSpan();
	double du = ( knotVector [ 0 ] [ span->at(1) + 1 ] - knotVector [ 0 ] [ span->at(1) ] ) / nseg;
	for ( i = 1; i <= nseg; i++ ) {
	  c [ 0 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
	  c [ 1 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
	  for ( k = 0; k < 2; k++ ) {
	    interp->local2global( cg [ k ], c [ k ], FEIIGAElementGeometryWrapper( this, iRule->giveKnotSpan() ) );
	    p [ k ].x = ( FPNum ) cg [ k ].at(1);
	    p [ k ].y = ( FPNum ) cg [ k ].at(2);
	    p [ k ].z = ( FPNum ) cg [ k ].at(3);
	  }
	  go = CreateLine3D(p);
	  EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
	  EGAttachObject(go, ( EObjectP ) this);
	  EMAddGraphicsToModel(ESIModel(), go);
	}
    }
}


void 
NURBSBeam3dElementDsg::drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType ut) 
{

    WCRec p [ 2 ];
    GraphicObj *go;
    int i, j, k, n, nseg;
    FloatArray u;
    FloatMatrix N;
    FloatArray ur;
    IntArray lc;
    FloatArray normal;
    FEInterpolation *interp = this->giveInterpolation();
    double defScale = gc.getDefScale();

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    EASValsSetLineStyle(SOLID_STYLE);
    EASValsSetFillStyle(FILL_SOLID);
    EASValsSetLineWidth(0);

    nseg = 80;

    const double *const *knotVector = interp->giveKnotVector();
    const IntArray *span;

    this->StructuralElementEvaluator::computeVectorOf(VM_Total, tStep, u);

    FloatArray c [ 2 ], cg [ 2 ], d, dg;
    dg.resize(3);
    double du;

    for ( j = 0; j < 2; j++ ) {
      c [ j ].resize(3);
      cg [ j ].resize(3);
    }

    // loop over individual integration rules (i.e., knot spans)
    for ( int ir = 0; ir < this->giveNumberOfIntegrationRules(); ir++ ) {
      IntegrationRule *iRule = this->giveIntegrationRule(ir);
      span = iRule->giveKnotSpan();
      // divide span locally to get finer geometry rep.
      du = ( knotVector [ 0 ] [ span->at(1) + 1 ] - knotVector [ 0 ] [ span->at(1) ] ) / nseg;
      for ( i = 1; i <= nseg; i++ ) {
	c [ 0 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
	c [ 1 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;

	for ( k = 0; k < 2; k++ ) {
	  // create a dummy ip's
	  GaussPoint gp(iRule, 999, c [ k ], 1.0, _3dBeam);

	  // compute displacements at gp
	  this->computeNMatrixAt(N, & gp);

	  // get local code numbers corresponding to ir
	  this->giveIntegrationElementLocalCodeNumbers(lc, this, gp.giveIntegrationRule());
	  ur.resize( N.giveNumberOfColumns() );
	  for ( n = 1; n <= lc.giveSize(); n++ ) {
	    ur.at(n) = u.at( lc.at(n) );
	  }

	  // interpolate displacements
	  d.beProductOf(N, ur);

	  interp->local2global( cg [ k ], c [ k ], FEIIGAElementGeometryWrapper( this, iRule->giveKnotSpan() ) );

	  this->computeNormal (normal, c[k], iRule->giveKnotSpan()->at(1)  );

	  // edita : todo - transformation of dipls to global  
	  dg = d;
	  // dg.at(1) = d.at(1) * normal.at(2) + d.at(2) * normal.at(1);
	  // dg.at(2) = d.at(1) * normal.at(1) + d.at(2) * normal.at(2);
	  // interp->local2global( dg , d , FEIIGAElementGeometryWrapper( this, iRule->giveKnotSpan() ) );
	  // dg.beVectorProductOf(normal, d);

	  p [ k ].x = ( FPNum ) ( cg [ k ].at(1) + dg.at(1) * defScale );
	  p [ k ].y = ( FPNum ) ( cg [ k ].at(2) + dg.at(2) * defScale );
	  p [ k ].z = ( FPNum ) ( cg [ k ].at(3) + dg.at(3) * defScale );
	}

	go =  CreateLine3D(p);
	EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
	EGAttachObject(go, ( EObjectP ) this);
	EMAddGraphicsToModel(ESIModel(), go);
      }
    }                 // end loop over knot spans (irules)
}
  

void 
NURBSBeam3dElementDsg::drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
  WCRec p [ 2 ];
  WCRec p2 [ 2 ];
  int indx;
    GraphicObj *go;
    int i, j, k, nseg;
    FloatArray u;
    FloatMatrix N;
    FloatArray ur, d, n, val;
    IntArray lc;
    FEInterpolation *interp = this->giveInterpolation();
    // double defScale = gc.getDefScale();

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_VARPLOT_PATTERN_LAYER);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    EASValsSetLineStyle(SOLID_STYLE);
    EASValsSetFillStyle(FILL_SOLID);
    EASValsSetLineWidth(0);

    nseg = 10;

    const double *const *knotVector = interp->giveKnotVector();
    const IntArray *span;
    double scale = gc.getDefScale();


    this->StructuralElementEvaluator::computeVectorOf(VM_Total, tStep, u);

    indx = gc.giveIntVarIndx();

    FloatArray c [ 2 ], cg [ 2 ];
    double du;
    double s [ 2 ];

    for ( j = 0; j < 2; j++ ) {
      c [ j ].resize(3);
      cg [ j ].resize(3);
    }

    // loop over individual integration rules (i.e., knot spans)
    for ( int ir = 0; ir < this->giveNumberOfIntegrationRules(); ir++ ) {
      IntegrationRule *iRule = this->giveIntegrationRule(ir);
      span = iRule->giveKnotSpan();
      // divide span locally to get finer geometry rep.
      du = ( knotVector [ 0 ] [ span->at(1) + 1 ] - knotVector [ 0 ] [ span->at(1) ] ) / nseg;
      for ( i = 1; i <= nseg; i++ ) {
	c [ 0 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
	c [ 1 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;

	for ( k = 0; k < 2; k++ ) {
	  // create a dummy ip's
	  GaussPoint gp(iRule, 999, c [ k ], 1.0, _3dBeam);
	  FloatArray strain;
	  if ((gc.giveIntVarType() == IST_StrainTensor)){
	    // copy to array of size 1x6
	    this->computeStrainVector(val, & gp, tStep, u);
	  }
	    else if (gc.giveIntVarType() == IST_StressTensor) {
		// copy to array of size 1x6
		this->computeStrainVector(strain, & gp, tStep, u);
		this->computeStressVector(val, strain, & gp, tStep);
	    } else{
		OOFEM_ERROR ("Internal Variable Type not recognized");
	    }
	  s [ k ] = val.at(indx);
	  fprintf (stderr, " %e",  s [ k ]);
	  
	  interp->local2global( cg [ k ], c [ k ], FEIIGAElementGeometryWrapper( this, iRule->giveKnotSpan() ) );
	  this->computeNormal (n, c[k], iRule->giveKnotSpan()->at(1) );
	  //n = {0, 1, 0};
	  p [ k ].x = ( FPNum ) ( cg [ k ].at(1) + n.at(1) * s[k] * scale );
	  p [ k ].y = ( FPNum ) ( cg [ k ].at(2) + n.at(2) * s[k] * scale );
	  p [ k ].z = ( FPNum ) ( cg [ k ].at(3) + n.at(3) * s[k] * scale );

	  p2 [ k ].x = ( FPNum ) ( cg [ 0 ].at(1) + (1-k) * n.at(1) * s[k] * scale );
	  p2 [ k ].y = ( FPNum ) ( cg [ 0 ].at(2) + (1-k) * n.at(2) * s[k] * scale );
	  p2 [ k ].z = ( FPNum ) ( cg [ 0 ].at(3) + (1-k) * n.at(3) * s[k] * scale );

	}

	go =  CreateLine3D(p);
	EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
	EGAttachObject(go, ( EObjectP ) this);
	EMAddGraphicsToModel(ESIModel(), go);

	go =  CreateLine3D(p2);
	EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
	EGAttachObject(go, ( EObjectP ) this);
	EMAddGraphicsToModel(ESIModel(), go);
      }
    }                // end loop over knot spans (irules)

	    this->computeNormal (n, c[1], this->giveNumberOfIntegrationRules()-1);
	  p2 [ 0 ].x = ( FPNum ) ( cg [ 1 ].at(1) + n.at(1) * s[1] * scale );
	  p2 [ 0 ].y = ( FPNum ) ( cg [ 1 ].at(2) + n.at(2) * s[1] * scale );
	  p2 [ 0 ].z = ( FPNum ) ( cg [ 1 ].at(3) + n.at(3) * s[1] * scale );

	  p2 [ 1 ].x = ( FPNum ) ( cg [ 1 ].at(1) );
	  p2 [ 1 ].y = ( FPNum ) ( cg [ 1 ].at(2) );
	  p2 [ 1 ].z = ( FPNum ) ( cg [ 1 ].at(3) );

	go =  CreateLine3D(p2);
	EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
	EGAttachObject(go, ( EObjectP ) this);
	EMAddGraphicsToModel(ESIModel(), go);
 
}

#endif
} // end namespace oofem
