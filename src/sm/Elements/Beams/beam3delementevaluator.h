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

#ifndef beam3delementevaluator_h
#define beam3delementevaluator_h

#include "../sm/Elements/structuralelementevaluator.h"

namespace oofem {
/**
 * General purpose Beam structural element evaluator.
 */
class Beam3dElementEvaluator : public StructuralElementEvaluator
{
public:
    Beam3dElementEvaluator() : StructuralElementEvaluator() { }

protected:
    /**
     * Assemble interpolation matrix at given IP.
     * In case of IGAElements, N is assumed to contain only nonzero interpolation functions.
     */
    virtual void computeNMatrixAt(FloatMatrix &answer, GaussPoint *gp);
    /**
     * Assembles the strain-displacement matrix of the receiver at given integration point.
     * In case of IGAElements, B is assumed to contain only contribution from nonzero interpolation functions.
     */
    virtual void computeBMatrixAt( FloatMatrix &answer, GaussPoint *gp);
    virtual void computeB1MatrixAt( FloatMatrix &answer, GaussPoint *gp);
    virtual void computeB2MatrixAt( FloatMatrix &answer, GaussPoint *gp);
    //virtual void computeBMatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx, int upperIndx);
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual double computeVolumeAround(GaussPoint *gp);
    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);
    void giveDofManDofIDMask(int inode, IntArray &answer) const {
	answer = {D_u, D_v, D_w, R_u, R_v, R_w};
    }

  //TO DO: edita
    
    // transformation
    // virtual void giveLocalCoordinates( FloatArray &lcoords, FloatArray &gcoords );
    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);
    //int computeIFGToLRotationMtrx(FloatMatrix &answer);

    int giveLocalCoordinateSystem(FloatMatrix &answer, FloatArray lcoords, const FEICellGeometry &cellgeo);
    void computeDofsGtoLMatrix(FloatMatrix &answer, FloatArray coords, int knotSpan);
    void computeDofsGtoLMatrix(FloatMatrix &answer, FloatArray coords,  const FEICellGeometry &cellgeo);
    int computeLoadGToLRotationMtrx(FloatMatrix &answer, GaussPoint *gp);
    //void computeLToDirectorRotationMatrix(FloatMatrix &answer1, FloatMatrix &answer2, FloatMatrix &answer3, FloatMatrix &answer4);
    //int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp);
    
     
    void boundaryEdgeGiveNodes(IntArray& bNodes, int boundary);

    void givedxds(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    double giveCurvature(const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    double giveTorsion(const FloatArray &lcoords, const FEICellGeometry &cellgeo);

    //loading

    void  computeBoundaryEdgeLoadVector(FloatArray &answer, BoundaryLoad *load, int edge, CharType type, ValueModeType mode, TimeStep *tStep, bool global);
    void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode);
    void giveTangent(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    void giveNormal(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    void giveBinormal(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);


    // draw
    void computeNormal (FloatArray &n, FloatArray c, int knotSpan);
}; // end of Beam3dElementEvaluator definition
} // end namespace oofem
#endif //beam2delementevaluator_h
