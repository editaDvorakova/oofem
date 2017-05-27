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

#include "../sm/Elements/Beams/degbeam.h"
#include "../sm/Elements/Beams/degbeamquad.h"
#include "../sm/Materials/structuralms.h"
#include "../sm/Materials/structuralmaterial.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "fei1dquad.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "../sm/CrossSections/variablecrosssection.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "load.h"
#include "boundaryload.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(DegeneratedBeamQuad3d);

FEI1dQuad DegeneratedBeamQuad3d :: interpolation(1);

DegeneratedBeamQuad3d :: DegeneratedBeamQuad3d(int n, Domain *aDomain) :
    DegeneratedBeam3d(n, aDomain)
{
    numberOfDofMans = 3;
    nGaussPoints = 8;
}

FEInterpolation *
DegeneratedBeamQuad3d :: giveInterpolation() const { return & interpolation; }


FEInterpolation *
DegeneratedBeamQuad3d :: giveInterpolation(DofIDItem id) const
{
    return & interpolation;
}

}
