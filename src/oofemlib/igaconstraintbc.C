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

#include "igaconstraintbc.h"
#include "classfactory.h"
#include "masterdof.h"
#include "floatmatrix.h"
#include "sparsemtrx.h"
#include "unknownnumberingscheme.h"
#include "function.h"
#include "timestep.h"
#include "datastream.h"
#include "contextioerr.h"
#include "node.h"
#include "domain.h"
#include "element.h"
#include "iga/iga.h"
#include "iga/feibspline.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "../sm/Elements/igaelements.h"


namespace oofem {
REGISTER_BoundaryCondition(IGAConstraintBC);

IGAConstraintBC :: IGAConstraintBC(int n, Domain *d) : ActiveBoundaryCondition(n, d),
    md( new Node(0, domain) )
{
    // this is internal lagrange multiplier used to enforce the receiver constrain
    // this allocates a new equation related to this constraint
    this->md->appendDof( new MasterDof( this->md.get(), ( DofIDItem ) ( d->giveNextFreeDofID() ) ) );
    this->lhsType.clear();
    this->rhsType.clear();
}


IGAConstraintBC :: ~IGAConstraintBC()
{
}


IRResultType IGAConstraintBC :: initializeFrom(InputRecord *ir)
{
    IRResultType result;
    rhsTf = 0;

    IR_GIVE_FIELD(ir, weights, _IFT_IGAConstraintBC_weights);
    IR_GIVE_FIELD(ir, rhs, _IFT_IGAConstraintBC_rhs);
    IR_GIVE_FIELD(ir, elements, _IFT_IGAConstraintBC_elements);
    if ( weights.giveSize() != elements.giveSize() ) {
        OOFEM_WARNING("Size mismatch, weights %d and elements %d", weights.giveSize(), elements.giveSize());
        return IRRT_BAD_FORMAT;
    }
    IR_GIVE_FIELD(ir, locArray, _IFT_IGAConstraintBC_locarray);
    if ( weights.giveSize() != elements.giveSize() ) {
        OOFEM_WARNING("Size mismatch, weights %d and locarray %d", weights.giveSize(), locArray.giveSize());
        return IRRT_BAD_FORMAT;
    }
    IR_GIVE_FIELD(ir, dofs, _IFT_IGAConstraintBC_dofs);
    IR_GIVE_OPTIONAL_FIELD(ir, weightsTf, _IFT_IGAConstraintBC_weightsfuncs);
    IR_GIVE_OPTIONAL_FIELD(ir, rhsTf, _IFT_IGAConstraintBC_rhsfuncs);

    coordSys = {1, 0, 0, 0, 1, 0};
    IR_GIVE_OPTIONAL_FIELD(ir, coordSys, _IFT_IGAConstraintBC_cs);

    IR_GIVE_FIELD(ir, lhsType, _IFT_IGAConstraintBC_lhstype);
    IR_GIVE_FIELD(ir, rhsType, _IFT_IGAConstraintBC_rhstype);

    return ActiveBoundaryCondition :: initializeFrom(ir);
}


void IGAConstraintBC :: giveLocArray(const UnknownNumberingScheme &r_s,  IntArray &locr, int &lambda_eq)
{
    int size = this->weights.giveSize();
    IntArray dofLoc;
    dofLoc.resize(1);

    locr.clear();
    // assemble location array
    for ( int i = 1; i <= size; i++ ) { // elements
	IGAElement *elem = dynamic_cast<IGAElement*> (this->domain->giveElement(elements.at(i)));
	int nnodes = elem->giveNumberOfDofManagers();
	

	BSplineInterpolation *interp = dynamic_cast<BSplineInterpolation*> (elem->giveInterpolation());
	IntArray knotSpan(1);
	int nCP = interp->giveNumberOfControlPoints(1);
	const FloatArray * knotVector = interp->giveKnotVector();
	int p = interp->giveDegree(1);
	knotSpan.at(1) = interp->findSpan(nCP, p, locArray.at(i), *knotVector); 
	IntArray mask;
	interp->giveKnotSpanBasisFuncMask(knotSpan, mask);
	int maskSize = mask.giveSize();
	
	for (int nn = mask.at(1); nn <= mask.at(maskSize); nn++) { // nodes

	    for (int j = 1; j <= 6; j++) { // dofs (6 dofs in node)
		Dof *idof = this->domain->giveDofManager( elem->giveNode(nn)->giveNumber())->giveDofWithID( j );
		dofLoc.at(1) = r_s.giveDofEquationNumber(idof);
		locr.followedBy(dofLoc);
	    }
	}
    }
    lambda_eq = r_s.giveDofEquationNumber( *md->begin() ); // E: first lambda
    
}
    
void IGAConstraintBC :: assemble(SparseMtrx &answer, TimeStep *tStep,
                                    CharType type, const UnknownNumberingScheme &r_s,
                                    const UnknownNumberingScheme &c_s)
{
    if ( !this->lhsType.contains( ( int ) type ) ) {
        return;
    }
 
    IntArray lambdaeq(1);
    IntArray locr;
    this->giveLocArray( r_s, locr, lambdaeq.at(1) );

    FloatArray contrib;//, contribt;
    contrib.clear();

    if ( this->isImposed(tStep) ) {
	int size = this->weights.giveSize();
        for ( int i = 1; i <= size; i++ ) { // loop over elements
            double factor = 1.; // timefunction
            if ( weightsTf.giveSize() ) {
                factor = domain->giveFunction( weightsTf.at(i) )->evaluateAtTime( tStep->giveIntrinsicTime() );
            }
	    IGAElement *elem = dynamic_cast<IGAElement*> (this->domain->giveElement(elements.at(i)));
	    FloatArray coords(1);
	    coords.at(1) = this->locArray.at(i);

	    // knotspan ?!?!?!
	    IntArray knotSpan(1);
	    BSplineInterpolation *bspInterp = dynamic_cast<BSplineInterpolation*> ( elem->giveInterpolation());
	    int nCP = bspInterp->giveNumberOfControlPoints(1);
	    const FloatArray* knotVector = bspInterp->giveKnotVector();
	    int p = bspInterp->giveDegree(1);
	    knotSpan.at(1) = bspInterp->findSpan(nCP, p, locArray.at(i), *knotVector); 

	    IntArray mask;
	    bspInterp->giveKnotSpanBasisFuncMask(knotSpan, mask);
	    int ir = mask.at(1) - 1;

 	    IntegrationRule *iRule = elem->giveIntegrationRule(ir);
	    GaussPoint *gp = iRule->getIntegrationPoint(0);
	    gp->setNaturalCoordinates(coords);
	   
	    FloatArray localN;
	    FEInterpolation *interp = elem->giveInterpolation();
	    bspInterp->evalN( localN, coords, FEIIGAElementGeometryWrapper( elem, elem->giveIntegrationRule(ir)->giveKnotSpan() ) );
	    bspInterp->evalN( localN, coords, FEIIGAElementGeometryWrapper( elem, &knotSpan ) );

	    FloatMatrix nMatrix;
	    nMatrix.beNMatrixOf(localN, 6);

	    // FloatMatrix T;
	    // elem->computeDofsGtoLMatrix(T, coords, FEIIGAElementGeometryWrapper( elem, elem->giveIntegrationRule(knotSpan)->giveKnotSpan() ) );
	    FloatMatrix R;
	    elem->computeKnotspanGtoLRotationMatrix(R, knotSpan);
	   
	    FloatMatrix T;
	    elem->computeLtoCSMatrix(T, coordSys, coords, FEIIGAElementGeometryWrapper( elem, &knotSpan));
	  


	    FloatMatrix contribTMatrix;
	    FloatMatrix NR;
	    NR.beProductOf(nMatrix, R);
	    FloatMatrix TNR;
	    TNR.beProductOf(T, NR);
	    contribTMatrix.beTranspositionOf(TNR);

	    FloatArray elemContrib;
	    elemContrib.beColumnOf(contribTMatrix, dofs.at(1));  
	    elemContrib.times(this->weights.at(i) * factor);
	    /*
	    elemContrib.clear();
	    int nSize = localN.giveSize();
	    for (int j = 1; j <= nSize; j++) {
		FloatArray help(6), help2;
		help.zero();
		help.at(dofs.at(1)) = localN.at(j) * this->weights.at(i) * factor;

		help2.beTProductOf(T, help);
		elemContrib.append(help2);	    
		}*/
	    contrib.append(elemContrib);

	}

	FloatMatrix help, helpT;
	help.clear();
	help.addSubVectorCol(contrib, 1, 1);
	helpT.beTranspositionOf(help);

        answer.assemble(lambdaeq, locr, helpT);
        answer.assemble(locr, lambdaeq, help);
    } else {
        // the bc is not imposed at specific time step, however in order to make the equation system regular
        // we initialize the allocated equation to the following form 1*labmda = 0, forcing lagrange multiplier
        // of inactive condition to be zero.
        FloatMatrix help(1, 1);
        help.at(1, 1) = 1.0;
        answer.assemble(lambdaeq, lambdaeq, help);
    }
}

void IGAConstraintBC :: assembleVector(FloatArray &answer, TimeStep *tStep,
                                          CharType type, ValueModeType mode,
                                          const UnknownNumberingScheme &s, FloatArray *eNorms)
{
    IntArray loc, lambdaeq(1);
    FloatArray vec(1);
    double factor = 1.;

    if ( !this->rhsType.contains( ( int ) type ) ) {
        return;
    }

    if ( type == InternalForcesVector ) {
        // compute true residual
        int size = this->weights.giveSize();
        Dof *mdof = *md->begin();
        Dof *idof;

        // assemble location array
        for ( int _i = 1; _i <= size; _i++ ) {
            factor = 1.;
            if ( weightsTf.giveSize() ) {
                factor = domain->giveFunction( weightsTf.at(_i) )->evaluateAtTime( tStep->giveIntrinsicTime() );
            }
            idof = this->domain->giveDofManager( this->elements.at(_i) )->giveDofWithID( this->dofs.at(1) );
            if ( s.giveDofEquationNumber(idof) ) {
                answer.at( s.giveDofEquationNumber(idof) ) += mdof->giveUnknown(mode, tStep) * this->weights.at(_i) * factor;
            }
            if ( s.giveDofEquationNumber( mdof ) ) {
                answer.at( s.giveDofEquationNumber( mdof ) ) += idof->giveUnknown(mode, tStep) * this->weights.at(_i) * factor;
            }
        }
    } else if ( type == ExternalForcesVector ) {
        // use rhs value

        if ( rhsTf ) {
            factor = domain->giveFunction(rhsTf)->evaluateAtTime( tStep->giveIntrinsicTime() );
        }
        this->giveLocArray( s, loc, lambdaeq.at(1) );
        vec.at(1) = rhs * factor;
        answer.assemble(vec, lambdaeq);
    }
}

void IGAConstraintBC :: giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    rows.resize(3);
    cols.resize(3);

    IntArray loc, lambdaeq(1);
    this->giveLocArray( r_s, loc, lambdaeq.at(1) );
    // column block
    rows [ 0 ] = loc;
    cols [ 0 ] = lambdaeq;
    // row block
    cols [ 1 ] = loc;
    rows [ 1 ] = lambdaeq;
    // diagonal enry (some sparse mtrx implementation requaire this)
    rows [ 2 ] = lambdaeq;
    cols [ 2 ] = lambdaeq;
}


contextIOResultType
IGAConstraintBC :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( mode & CM_Definition ) {
        if ( ( iores = weights.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = weightsTf.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = elements.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = dofs.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( !stream.write(rhs) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(rhsTf) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( ( iores = lhsType.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = rhsType.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    if ( ( iores = md->saveContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
IGAConstraintBC :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( mode & CM_Definition ) {
        if ( ( iores = weights.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = weightsTf.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = elements.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = dofs.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( !stream.read(rhs) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(rhsTf) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( ( iores = lhsType.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = rhsType.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    if ( ( iores = md->restoreContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}
} //end of oofem namespace
