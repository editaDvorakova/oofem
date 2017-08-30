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

#ifndef beam3delementevaluatorbbar_h
#define beam3delementevaluatorbbar_h

#include "../sm/Elements/Beams/beam3delementevaluator.h"
#include "../sm/Elements/structuralelementevaluator.h"

namespace oofem {
/**
 * General purpose Beam structural element evaluator.
 */
class Beam3dElementEvaluatorBbar : public Beam3dElementEvaluator
{
public:
    Beam3dElementEvaluatorBbar() : Beam3dElementEvaluator() { }

protected:
    FloatMatrix invM;
   
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep); 
	//{ StructuralElementEvaluator :: computeStiffnessMatrix(answer, rMode, tStep); }
    /**
     * Assembles the strain-displacement matrix of the receiver at given integration point.
     * In case of IGAElements, B is assumed to contain only contribution from nonzero interpolation functions.
     */

    virtual const FloatMatrix* computeBbarMassMatrix();
    // virtual void computeBMatrixAt( FloatMatrix &answer, GaussPoint *gp);


    virtual void computeKBbarMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);

    //  virtual void computeKMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    
    
    virtual void computeBendingPartOfBMatrixAt(FloatMatrix &Bb, GaussPoint *gp);

    //    int giveIntegrationElementLocalCodeNumbers(IntArray &answer, Element *elem, IntegrationRule *ie);

    // void giveKnotSpanAt(IntArray &knotSpan, double lcoord);
    // virtual int giveIntegrationElementLocalCodeNumbers(IntArray &answer, Element *elem, IntegrationRule *ie);


    
    // transformation
    // virtual void giveLocalCoordinates( FloatArray &lcoords, FloatArray &gcoords );

    //int computeIFGToLRotationMtrx(FloatMatrix &answer);

    // void computeDofsGtoLMatrix(FloatMatrix &answer, FloatArray coords,  const FEICellGeometry &cellgeo);
    //void computeLToDirectorRotationMatrix(FloatMatrix &answer1, FloatMatrix &answer2, FloatMatrix &answer3, FloatMatrix &answer4);
    //int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp);
    // IntegrationRule *setIR(int i) ;
    // int giveNumberOfIR(int knotSpan);

    // int giveLocalCoordinateSystem(FloatMatrix &answer, FloatArray lcoords, const FEICellGeometry &cellgeo);
   
    // draw    

    // void computeNormal (FloatArray &n, FloatArray c, int knotSpan);
}; // end of Beam3dElementEvaluator definition
} // end namespace oofem
#endif //beam3delementevaluator_h
