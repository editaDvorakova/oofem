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

#ifndef igaexportmodule_h
#define igaexportmodule_h

#include "exportmodule.h"
#include "intarray.h"
#include "nodalrecoverymodel.h"
#include "internalstatevaluetype.h"
#include "unknowntype.h"
#include "valuemodetype.h"

#define _IFT_IGAExportModule_Name "iga"
#define _IFT_IGAExportModule_divisions "divisions"

namespace oofem {
class DofManager;

/**
 * Represents VTK (Visualization Toolkit) export module. It uses VTK file format, Unstructured grid dataset.
 * There is built in support for Region By Region output, taking care about possible nonsmooth character of
 * some internal variables at region boundaries. This, however, is rather complication and since application
 * of VTK is naturally in 3D, the corresponding sections are commented out.
 */
class OOFEM_EXPORT IGAExportModule : public ExportModule
{   
protected:
  // number of segments
  int divisions = {10};

public:
  
    /// Constructor. Creates empty Output Manager with number n.
    IGAExportModule(int n, EngngModel * e);
    /// Destructor
    virtual ~IGAExportModule();
    /// Initializes receiver according to object description stored in input record.
    virtual IRResultType initializeFrom(InputRecord *ir);
    /**
     * Writes the output. Abstract service.
     * @param tStep Time step.
     * @param forcedOutput If true, no testTimeStepOutput should be done.
     */
    virtual void doOutput(TimeStep *tStep, bool forcedOutput = false);
  
    virtual void terminate() { }
    /// Returns class name of the receiver.
    virtual const char *giveClassName() const { return "IGAExportModule"; }

    FILE * giveOutputStream(const char *extension);
protected:
    /**
     * Tests if given time step output is required.
     * @param tStep Time step to check.
     * @return True if output required.
     */
    bool testTimeStepOutput(TimeStep *tStep);
    /**
     * Test if domain output is required.
     * @return True if required.
     */
    bool testDomainOutput(int n);

};
} // end namespace oofem
#endif // igaexportmodule_h

