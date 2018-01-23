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
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef feinurbsline2d_h
#define feinurbsline2d_h

#include "feibspline.h"
#include "floatarray.h"
#include "floatmatrix.h"
//#include "mathfem.h"

namespace oofem {
/**
 * Interpolation class for NURBS.
 */
class NURBSInterpolationLine2d : public BSplineInterpolation
{
public:
    NURBSInterpolationLine2d(int nsd, int fsd) : BSplineInterpolation(nsd, fsd) { }
    virtual ~NURBSInterpolationLine2d();

    virtual void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual int  global2local(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
    {
        OOFEM_ERROR("NURBSInterpolation :: global2local - Not yet implemented.");
        return 0;
    }
    virtual void giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    // double  givedY(const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual const char *giveClassName() const { return "NURBSInterpolationLine2d"; }

    //virtual double givedR(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
};
} // end namespace oofem
#endif // feinurbs_h
