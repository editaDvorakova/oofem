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

#ifndef igaelements_h
#define igaelements_h

#include "iga/iga.h"
#include "iga/feibspline.h"
#include "iga/feinurbs.h"
#include "iga/feinurbsline2d.h"
#include "iga/feitspline.h"
#include "../sm/Elements/PlaneStress/planestresselementevaluator.h"
#include "../sm/Elements/Beams/beam2delementevaluator.h"
#include "../sm/Elements/Beams/beam2delementevaluatordsg.h"
#include "../sm/Elements/3D/space3delementevaluator.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "matresponsemode.h"
#include "mathfem.h"

#define _IFT_BsplinePlaneStressElement_Name "bsplineplanestresselement"
#define _IFT_NURBSPlaneStressElement_Name "nurbsplanestresselement"
#define _IFT_TSplinePlaneStressElement_Name "tsplineplanestresselement"
#define _IFT_NURBSSpace3dElement_Name "nurbs3delement"
#define _IFT_NURBSBeam2dElement_Name "nurbsbeam2delement"
#define _IFT_NURBSBeam2dElementDsg_Name "nurbsbeam2delementdsg"

namespace oofem {
class BsplinePlaneStressElement : public IGAElement, public PlaneStressStructuralElementEvaluator
{
protected:
    BSplineInterpolation interpolation;

public:
    BsplinePlaneStressElement(int n, Domain * aDomain);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int checkConsistency();

    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep) {
        PlaneStressStructuralElementEvaluator :: giveCharacteristicMatrix(answer, mtrx, tStep);
    }
    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *t) {
        PlaneStressStructuralElementEvaluator :: giveCharacteristicVector(answer, type, mode, t);
    }

    virtual FEInterpolation *giveInterpolation() const { return const_cast< BSplineInterpolation * >(& this->interpolation); }
    virtual Element *giveElement() { return this; }
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const {
        PlaneStressStructuralElementEvaluator :: giveDofManDofIDMask(inode, answer);
    }
    virtual int computeNumberOfDofs() { return numberOfDofMans * 2; }
    virtual void updateInternalState(TimeStep *tStep) { PlaneStressStructuralElementEvaluator :: updateInternalState(tStep); }
    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_BsplinePlaneStressElement_Name; }
    virtual const char *giveClassName() const { return "BsplinePlaneStressElement"; }

#ifdef __OOFEG
    // Graphics output
    virtual void drawScalar(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType ut) {
        drawIGAPatchDeformedGeometry(this, this, gc, tStep, ut);
    }
#endif

protected:
    virtual int giveNsd() { return 2; }
    virtual int giveFsd() { return 2; }
};


class NURBSPlaneStressElement : public IGAElement, public PlaneStressStructuralElementEvaluator
{
protected:
    NURBSInterpolation interpolation;

public:
    NURBSPlaneStressElement(int n, Domain * aDomain);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int checkConsistency();

    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep) {
        PlaneStressStructuralElementEvaluator :: giveCharacteristicMatrix(answer, mtrx, tStep);
    }
    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *t) {
        PlaneStressStructuralElementEvaluator :: giveCharacteristicVector(answer, type, mode, t);
    }

    virtual FEInterpolation *giveInterpolation() const { return const_cast< NURBSInterpolation * >(& this->interpolation); }
    virtual Element *giveElement() { return this; }
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const {
        PlaneStressStructuralElementEvaluator :: giveDofManDofIDMask(inode, answer);
    }
    virtual int computeNumberOfDofs() { return numberOfDofMans * 2; }
    virtual void updateInternalState(TimeStep *tStep) { PlaneStressStructuralElementEvaluator :: updateInternalState(tStep); }
    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_NURBSPlaneStressElement_Name; }
    virtual const char *giveClassName() const { return "NURBSPlaneStressElement"; }
#ifdef __OOFEG
    //
    // Graphics output
    //
    virtual void  drawScalar(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType ut) {
        drawIGAPatchDeformedGeometry(this, this, gc, tStep, ut);
    }

#endif

protected:
    virtual int giveNsd() { return 2; }
    virtual int giveFsd() { return 2; }
};



class TSplinePlaneStressElement : public IGATSplineElement, public PlaneStressStructuralElementEvaluator
{
protected:
    TSplineInterpolation interpolation;

public:
    TSplinePlaneStressElement(int n, Domain * aDomain);

    virtual IRResultType initializeFrom(InputRecord *ir) {
        IGATSplineElement :: initializeFrom(ir);
        //PlaneStressStructuralElementEvaluator::initializeFrom(ir);
        return IRRT_OK;
    }

    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep) {
        PlaneStressStructuralElementEvaluator :: giveCharacteristicMatrix(answer, mtrx, tStep);
    }
    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *t) {
        PlaneStressStructuralElementEvaluator :: giveCharacteristicVector(answer, type, mode, t);
    }

    virtual FEInterpolation *giveInterpolation() const { return const_cast< TSplineInterpolation * >(& this->interpolation); }
    virtual Element *giveElement() { return this; }
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const {
        PlaneStressStructuralElementEvaluator :: giveDofManDofIDMask(inode, answer);
    }
    virtual int computeNumberOfDofs() { return numberOfDofMans * 2; }
    virtual void updateInternalState(TimeStep *tStep) { PlaneStressStructuralElementEvaluator :: updateInternalState(tStep); }
    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_TSplinePlaneStressElement_Name; }
    virtual const char *giveClassName() const { return "TSplinePlaneStressElement"; }
#ifdef __OOFEG
    // Graphics output
    virtual void  drawScalar(oofegGraphicContext &gc, TimeStep *tStep);
#endif

protected:
    virtual int giveNsd() { return 2; }
    virtual int giveFsd() { return 2; }

};


class NURBSSpace3dElement : public IGAElement, public Space3dStructuralElementEvaluator
{
protected:
    NURBSInterpolation interpolation;

public:
    NURBSSpace3dElement(int n, Domain * aDomain);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int checkConsistency();

    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep) {
        Space3dStructuralElementEvaluator :: giveCharacteristicMatrix(answer, mtrx, tStep);
    }
    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *t) {
        Space3dStructuralElementEvaluator :: giveCharacteristicVector(answer, type, mode, t);
    }

    virtual FEInterpolation *giveInterpolation() const { return const_cast< NURBSInterpolation * >(& this->interpolation); }
    virtual Element *giveElement() { return this; }

    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const {
        Space3dStructuralElementEvaluator :: giveDofManDofIDMask(inode, answer);
    }
    virtual int computeNumberOfDofs() { return numberOfDofMans * 3; }
    virtual void updateInternalState(TimeStep *tStep) { Space3dStructuralElementEvaluator :: updateInternalState(tStep); }
    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_NURBSSpace3dElement_Name; }
    virtual const char *giveClassName() const { return "NURBSSpace3dElement"; }
#ifdef __OOFEG
    // Graphics output
    virtual void drawScalar(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType ut) {
        drawIGAPatchDeformedGeometry(this, this, gc, tStep, ut);
    }
#endif

protected:
    virtual int giveNsd() { return 3; }
    virtual int giveFsd() { return 3; }
};

class NURBSBeam2dElement : public IGAElement, public Beam2dElementEvaluator
{
protected:
    NURBSInterpolationLine2d interpolation;

public:
    NURBSBeam2dElement(int n, Domain * aDomain);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int checkConsistency();

    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep) {
        Beam2dElementEvaluator :: giveCharacteristicMatrix(answer, mtrx, tStep);
    }
    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *t) {
        Beam2dElementEvaluator :: giveCharacteristicVector(answer, type, mode, t);
    }

    virtual FEInterpolation *giveInterpolation() const { return const_cast< NURBSInterpolationLine2d * >(& this->interpolation); }
    virtual Element *giveElement() { return this; }
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const {
        Beam2dElementEvaluator :: giveDofManDofIDMask(inode, answer);
    }
    virtual int computeNumberOfDofs() { return numberOfDofMans * 3; } 
    virtual void updateInternalState(TimeStep *tStep) { Beam2dElementEvaluator :: updateInternalState(tStep); }
    // transformation
     virtual bool computeGtoLRotationMatrix(FloatMatrix &answer) {
      return (Beam2dElementEvaluator :: computeGtoLRotationMatrix(answer));
      }
    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_NURBSBeam2dElement_Name; }
    virtual const char *giveClassName() const { return "NURBSBeam2dElement"; }
    

    // edita - OOFEG ?
    void computeNormal (FloatArray &n, FloatArray c, int knotSpan) { Beam2dElementEvaluator :: computeNormal ( n, c, knotSpan); }
#ifdef __OOFEG
    //
    // Graphics output
    //
    //virtual void  drawScalar(oofegGraphicContext &gc, TimeStep *tStep);
    /*
    virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType ut) ;
    virtual void drawScalar(oofegGraphicContext &gc, TimeStep *tStep);
    */
#endif

protected:
    virtual int giveNsd() { return 3; }
    virtual int giveFsd() { return 1; }
};


class NURBSBeam2dElementDsg : public IGAElement, public Beam2dElementEvaluatorDsg
{
protected:
    NURBSInterpolationLine2d interpolation;

public:
    NURBSBeam2dElementDsg(int n, Domain * aDomain);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int checkConsistency();

    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep) {
        Beam2dElementEvaluatorDsg :: giveCharacteristicMatrix(answer, mtrx, tStep);
    }
    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *t) {
        Beam2dElementEvaluatorDsg :: giveCharacteristicVector(answer, type, mode, t);
    }

    virtual FEInterpolation *giveInterpolation() const { return const_cast< NURBSInterpolationLine2d * >(& this->interpolation); }
    virtual Element *giveElement() { return this; }
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const {
        Beam2dElementEvaluatorDsg :: giveDofManDofIDMask(inode, answer);
    }
    virtual int computeNumberOfDofs() { return numberOfDofMans * 3; } 
    virtual void updateInternalState(TimeStep *tStep) { Beam2dElementEvaluatorDsg :: updateInternalState(tStep); }
    // transformation
     virtual bool computeGtoLRotationMatrix(FloatMatrix &answer) {
      return (Beam2dElementEvaluatorDsg :: computeGtoLRotationMatrix(answer));
      }
    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_NURBSBeam2dElementDsg_Name; }
    virtual const char *giveClassName() const { return "NURBSBeam2dElementDsg"; }
    

    // edita - OOFEG ?
    void computeNormal (FloatArray &n, FloatArray c, int knotSpan) { Beam2dElementEvaluatorDsg :: computeNormal ( n, c, knotSpan); }
#ifdef __OOFEG
    //
    // Graphics output
    //
    //virtual void  drawScalar(oofegGraphicContext &gc, TimeStep *tStep);

    virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType ut) ;
    virtual void drawScalar(oofegGraphicContext &gc, TimeStep *tStep);

#endif

protected:
    virtual int giveNsd() { return 3; }
    virtual int giveFsd() { return 1; }
};
} // end namespace oofem
#endif //igaelements_h
