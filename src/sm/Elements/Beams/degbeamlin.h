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

#ifndef degbeamlin_h
#define degbeamlin_h

#include "../sm/Elements/nlstructuralelement.h"
#include "spatiallocalizer.h"
#include "load.h"
#include "degbeam.h"



#define _IFT_DegeneratedBeamLin3d_Name "degeneratedbeamlin3d"

namespace oofem {
    class FEI1dLin;

/// Comment or uncomment the following line to force full or reduced integration
// #define DegeneratedBeam3d_reducedShearIntegration

/**
 * This class implements an quad element based on Mixed Interpolation of Tensorial Components (MITC).
 * This element is a shell element suitable for both thin and thick shells.
 * The element has 24 DOFs (u,v,w-displacements and three rotations) in each node
 *
 * Tasks:
 * - calculating its B,D matrices and dV.
 */

class DegeneratedBeamLin3d : public DegeneratedBeam3d
{
protected:
    /// Element geometry approximation
    static FEI1dLin interpolation;


public:

    DegeneratedBeamLin3d(int n, Domain *d);
    virtual ~DegeneratedBeamLin3d() { }


private:

    // definition & identification
    virtual const char *giveClassName() const { return "DegeneratedBeamLin3d"; }
    virtual const char *giveInputRecordName() const { return _IFT_DegeneratedBeamLin3d_Name; }

};
} // end namespace oofem
#endif // degbeam_h
