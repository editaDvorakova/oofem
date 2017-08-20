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

#include "ghexportmodule.h"
#include "timestep.h"
#include "engngm.h"
#include "domain.h"

#include "iga/iga.h"
#include "../sm/Elements/igaelements.h"
#include "gausspoint.h"
#include "node.h"
#include "materialinterface.h"
#include "mathfem.h"
#include "cltypes.h"
#include "element.h"
#include "material.h"
#include "classfactory.h"
#include "crosssection.h"
#include "dof.h"

#include <cstdarg>



namespace oofem {
REGISTER_ExportModule(GHExportModule)


  GHExportModule :: GHExportModule(int n, EngngModel *e) : ExportModule(n,e)
{ }


GHExportModule :: ~GHExportModule()
{ }

IRResultType
GHExportModule :: initializeFrom(InputRecord *ir)
{
    IRResultType result;  // Required by IR_GIVE_FIELD macro
    IR_GIVE_OPTIONAL_FIELD(ir, divisions, _IFT_GHExportModule_divisions);

    return ExportModule :: initializeFrom(ir);
}

bool
  GHExportModule :: testTimeStepOutput(TimeStep *tStep)
{
  return true;
}

bool
GHExportModule :: testDomainOutput(int n)
{
  return true;
}


FILE *
GHExportModule :: giveOutputStream(const char *extension)
{
    std :: string fileName;
    FILE *answer;
    fileName = this->emodel->giveOutputBaseFileName()+ "." + extension + ".gh.out";
    if ( ( answer = fopen(fileName.c_str(), "w") ) == NULL ) {
        OOFEM_ERROR("failed to open file %s", fileName.c_str());
    }
    return answer;
}


void
GHExportModule :: doOutput(TimeStep *tStep, bool forcedOutput)
{
    if ( !( testTimeStepOutput(tStep) || forcedOutput ) ) {
        return;
    }


    // DISPLACEMENTS
    FILE *stream = this->giveOutputStream("d");
    Domain *d  = emodel->giveDomain(1);

    int nnodes = d->giveNumberOfDofManagers();
    FloatArray answer;
    for ( int inode = 1; inode <= nnodes; inode++ ) { 
      d->giveNode( inode )->giveUnknownVectorOfType(answer, DisplacementVector, VM_Total, tStep);
        for ( int i = 1; i <= answer.giveSize(); i++ ) {
            fprintf( stream, "%e ", answer.at(i) );
        }

        fprintf(stream, "\n");
    }
    fclose(stream);

    // FORCES
    stream = this->giveOutputStream("f");
    d  = emodel->giveDomain(1);
    
    FloatMatrix internalForces;
    int nelem = d->giveNumberOfElements();
    for ( int ielem = 1; ielem <= nelem; ielem++ ) {
	Element *elem = d->giveElement(ielem);

	NURBSBeam2dElementDsg *beamelem = dynamic_cast<NURBSBeam2dElementDsg*>( d->giveElement(ielem) );
	beamelem->computeInternalForces(internalForces, divisions, tStep);

      for ( int i = 1; i <= internalForces.giveNumberOfRows(); i++ ) {
	for (int j = 1; j <= internalForces.giveNumberOfColumns(); j++) {  
	  fprintf( stream, "%e ", internalForces.at(i, j) );
	}
        fprintf(stream, "\n");
      }
    }
    fclose(stream);

}

}
