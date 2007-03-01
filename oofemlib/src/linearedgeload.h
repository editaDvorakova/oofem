/* $Header: /home/cvs/bp/oofem/oofemlib/src/linearedgeload.h,v 1.1 2003/04/06 14:08:24 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2000   Borek Patzak                                       



         Czech Technical University, Faculty of Civil Engineering,
     Department of Structural Mechanics, 166 29 Prague, Czech Republic
                                                                               
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                                                                              
*/

//   ******************************
//   *** CLASS LINEAR EDGE LOAD ***
//   ******************************

#ifndef linearedgeload_h 

#include "boundaryload.h"
#include "gausspnt.h"
#include "flotarry.h"

class TimeStep ;


class LinearEdgeLoad : public BoundaryLoad
{
/*
   This class implements a linear boundary load (force, moment,...) that acts
   on straight segment. 
  A boundary load is usually attribute of one or more elements.

 DESCRIPTION
   The boundary load describes its geometry and values (it is assumed, that user will specify 
 all necessary dofs) on element boundary using isoparametric approximation. 
 Elements can request the order of approximation (for setting up the appropriate 
 integration rule order) and the array of values (for each dof) at specific integration point
 on the boundary. 

 Elements must take care, on which boundary the load acts on (side number, ...).

 REMARK
   This class is not restricted to structural problems. For example, in ther-
   mal analysis, a boundary load load would be a  heat source.
*/
protected:
  
 /// Coordinates of start and end point
 FloatArray startCoords, endCoords;
 BL_FormulationType formulation;


public:
  LinearEdgeLoad  (int i,Domain* d) : BoundaryLoad(i,d) {}           // constructor
  
  int          giveApproxOrder () {return 1;}
  IRResultType initializeFrom (InputRecord* ir);
 /** Setups the input record string of receiver
  @param str string to be filled by input record
  @param keyword print record keyword (default true)
  */
 virtual int giveInputRecordString(std::string &str, bool keyword = true);
 bcGeomType giveBCGeoType () const {return EdgeLoadBGT;}
 BL_FormulationType   giveFormulationType () {return formulation;}

  classType    giveClassID() const {return LinearEdgeLoadClass;}
  const char*        giveClassName () const {return "LinearEdgeLoad";}

protected:
  void computeNArray (FloatArray& answer, FloatArray& coords) const;

} ;


#define linearedgeload_h
#endif






