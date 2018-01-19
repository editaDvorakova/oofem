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

#include "feinurbsline3d.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "iga.h"

namespace oofem {
// optimized version of A4.4 for d=1
#define OPTIMIZED_VERSION_A4dot4


NURBSInterpolationLine3d :: ~NURBSInterpolationLine3d() { }


void NURBSInterpolationLine3d :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FEIIGAElementGeometryWrapper *gw = ( FEIIGAElementGeometryWrapper * ) & cellgeo;
    IntArray span(fsd);
    double sum = 0.0, val;
    int count, c = 1, i, k, ind, uind;
#ifdef HAVE_VARIABLE_ARRAY_SIZE
    FloatArray N [ fsd ];
#else
    FloatArray *N = new FloatArray [ fsd ];
#endif

    if ( gw->knotSpan ) {
        span = * gw->knotSpan;
    } else {
        for ( i = 0; i < fsd; i++ ) {
            span(i) = this->findSpan(numberOfControlPoints [ i ], degree [ i ], lcoords(i), knotVector [ i ]);
        }
    }

    for ( i = 0; i < fsd; i++ ) {
        this->basisFuns(N [ i ], span(i), lcoords(i), degree [ i ], knotVector [ i ]);
    }

    count = giveNumberOfKnotSpanBasisFunctions(span);
    answer.resize(count);

    if ( fsd == 1 ) {
        uind = span(0) - degree [ 0 ];
        ind = uind + 1;
        for ( k = 0; k <= degree [ 0 ]; k++ ) {
            answer.at(c++) = val = N [ 0 ](k) * cellgeo.giveVertexCoordinates(ind + k)->at(nsd+1);       // Nu*w
            sum += val;
        }
    } else {
        OOFEM_ERROR("evalN not implemented for fsd = %d", fsd);
    }

    while ( count ) {
        answer.at(count--) /= sum;
    }

#ifndef HAVE_VARIABLE_ARRAY_SIZE
    delete [] N;
#endif
}


     
 void NURBSInterpolationLine3d :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    //
    // Based on Algorithm A4.4 (p. 137) for d=1
    //
    FEIIGAElementGeometryWrapper *gw = ( FEIIGAElementGeometryWrapper * ) & cellgeo;
    const FloatArray *vertexCoordsPtr;
    IntArray span(fsd);
    double w, weight;
    int i, k, ind, uind;
    
    jacobianMatrix.resize(3, 3);
    
#ifdef HAVE_VARIABLE_ARRAY_SIZE
    FloatMatrix ders [ fsd ];
#else
    FloatMatrix *ders = new FloatMatrix [ fsd ];
#endif

    if ( gw->knotSpan ) {
        span = * gw->knotSpan;
    } else {
        for ( i = 0; i < fsd; i++ ) {
            span(i) = this->findSpan(numberOfControlPoints [ i ], degree [ i ], lcoords(i), knotVector [ i ]);
        }
    }

    for ( i = 0; i < fsd; i++ ) {
        this->dersBasisFuns(1, lcoords(i), span(i), degree [ i ], knotVector [ i ], ders [ i ]);
    }



    FloatArray *Aders = new FloatArray [ 3 ];
    FloatArray wders;          // 0th and 1st derivatives in w direction on BSpline

    for ( i = 0; i < 3; i++ ) {
        Aders [ i ].resize(fsd + 1);
        Aders [ i ].zero();
    }

    wders.resize(fsd + 1);
    wders.zero();

    if ( fsd == 1 ) {
        // calculate values and derivatives of nonrational Bspline curve with weights at first (Aders, wders)
        uind = span(0) - degree [ 0 ];
        ind = uind + 1;
        for ( k = 0; k <= degree [ 0 ]; k++ ) {
            vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind + k);
            w = vertexCoordsPtr->at(nsd+1);
            Aders [ 0 ](0) += ders [ 0 ](0, k) * vertexCoordsPtr->at(1) * w;   // xw=sum(Nu*x*w)
	    Aders [ 1 ](0) += ders [ 0 ](0, k) * vertexCoordsPtr->at(2) * w;   // yw=sum(Nu*y*w)
	    Aders [ 2 ](0) += ders [ 0 ](0, k) * vertexCoordsPtr->at(3) * w;   // yw=sum(Nu*z*w)
            wders(0)    += ders [ 0 ](0, k) * w;                               // w=sum(Nu*w)

            Aders [ 0 ](1) += ders [ 0 ](1, k) * vertexCoordsPtr->at(1) * w;   // dxw/du=sum(dNu/du*x*w)
            Aders [ 1 ](1) += ders [ 0 ](1, k) * vertexCoordsPtr->at(2) * w;   // dyw/du=sum(dNu/du*y*w)
            Aders [ 2 ](1) += ders [ 0 ](1, k) * vertexCoordsPtr->at(3) * w;   // dyw/du=sum(dNu/du*z*w)
            wders(1)    += ders [ 0 ](1, k) * w;                               // dw/du=sum(dNu/du*w)
        }

        weight = wders(0);

        // calculation of jacobian matrix according to Eq 4.7
        jacobianMatrix(0, 0) = ( Aders [ 0 ](1) - wders(1) * Aders [ 0 ](0) / weight ) / weight; // dx/du
        jacobianMatrix(1, 1) = ( Aders [ 1 ](1) - wders(1) * Aders [ 1 ](0) / weight ) / weight; // dy/du
        jacobianMatrix(2, 2) = ( Aders [ 2 ](1) - wders(1) * Aders [ 2 ](0) / weight ) / weight; // dz/du
    }  else {
      OOFEM_ERROR("giveTransformationJacobianMatrix not implemented for fsd = %d", fsd);
    }

    delete [] ders;
    delete [] Aders;


}


double
NURBSInterpolationLine3d :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo){

  FloatMatrix jacobian;
  this->giveJacobianMatrixAt(jacobian, lcoords, cellgeo);
    
  double Jacob = sqrt(jacobian(0,0)*jacobian(0,0) + jacobian(1,1)*jacobian(1,1)  + jacobian(2,2)*jacobian(2,2));
  
  if ( fabs(Jacob) < 1.0e-10 ) {
    OOFEM_ERROR("giveTransformationJacobianMatrix - zero Jacobian");
  }
  
  
  return Jacob;
}
  

void NURBSInterpolationLine3d :: lowerDegree()
{
    // change degree
    this->degree[0] = this->degree[0]-1;

    int kvSize = 0;
    int n = this->knotMultiplicity[0].giveSize();
    for (int i = 0; i < n; i++) {
	kvSize += this->knotMultiplicity[0](i);
    }

    // change knotVector
    for (int i = 0; i < kvSize-2; i++) {
	knotVector[0][i] = knotVector[0][i+1]; 
    }

    // change multiplicity
    knotMultiplicity[0](0) = knotMultiplicity[0](0) -1;
    knotMultiplicity[0](n-1) = knotMultiplicity[0](n-1) -1;

    
}

void NURBSInterpolationLine3d :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    /* Based on SurfacePoint A4.3 implementation*/
    FEIIGAElementGeometryWrapper *gw = ( FEIIGAElementGeometryWrapper * ) & cellgeo;
    const FloatArray *vertexCoordsPtr;
    IntArray span(fsd);
    double w, weight = 0.0;
    int i, k, ind, uind;
#ifdef HAVE_VARIABLE_ARRAY_SIZE
    FloatArray N [ fsd ];
#else
    FloatArray *N = new FloatArray [ fsd ];
#endif

    if ( gw->knotSpan ) {
        span = * gw->knotSpan;
    } else {
        for ( i = 0; i < fsd; i++ ) {
            span(i) = this->findSpan(numberOfControlPoints [ i ], degree [ i ], lcoords(i), knotVector [ i ]);
        }
    }

    for ( i = 0; i < fsd; i++ ) {
        this->basisFuns(N [ i ], span(i), lcoords(i), degree [ i ], knotVector [ i ]);
    }
    //change for line in 2d, instead of nsd -> 2
    answer.resize(3);
    answer.zero();

    if ( fsd == 1 ) {
        uind = span(0) - degree [ 0 ];
        ind = uind + 1;
        for ( k = 0; k <= degree [ 0 ]; k++ ) {
            vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind + k);
            w = vertexCoordsPtr->at(nsd+1);
            answer(0) += N [ 0 ](k) * vertexCoordsPtr->at(1) * w;       // xw=sum(Nu*x*w)
            answer(1) += N [ 0 ](k) * vertexCoordsPtr->at(2) * w;       // yw=sum(Nu*y*w)
            answer(2) += N [ 0 ](k) * vertexCoordsPtr->at(3) * w;       // zw=sum(Nu*z*w)
            weight    += N [ 0 ](k) * w;                                // w=sum(Nu*w)
        }
    } else {
        OOFEM_ERROR("local2global not implemented for fsd = %d", fsd);
    }

    answer.times(1.0 / weight);

#ifndef HAVE_VARIABLE_ARRAY_SIZE
    delete [] N;
#endif
}



double NURBSInterpolationLine3d :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FEIIGAElementGeometryWrapper *gw = ( FEIIGAElementGeometryWrapper * ) & cellgeo;
    const FloatArray *vertexCoordsPtr;
    FloatMatrix jacobian(3, 3);
    IntArray span(fsd);
    double Jacob = 0., product, w, weight;
    int count, cnt, i, k, ind, uind;
#ifdef HAVE_VARIABLE_ARRAY_SIZE
    FloatMatrix ders [ fsd ];
#else
    FloatMatrix *ders = new FloatMatrix [ fsd ];
#endif

    if ( gw->knotSpan ) {
        span = * gw->knotSpan;
    } else {
        for ( i = 0; i < fsd; i++ ) {
            span(i) = this->findSpan(numberOfControlPoints [ i ], degree [ i ], lcoords(i), knotVector [ i ]);
        }
    }

    for ( i = 0; i < fsd; i++ ) {
        this->dersBasisFuns(3, lcoords(i), span(i), degree [ i ], knotVector [ i ], ders [ i ]);
    }

    count = giveNumberOfKnotSpanBasisFunctions(span);
    answer.resize(count, 3);


    FloatArray *Aders = new FloatArray [ 3 ];
    FloatArray wders;          // 0th and 1st and 2nd derivatives in w direction on BSpline
    for ( i = 0; i < 3; i++ ) {
        Aders [ i ].resize(fsd + 3);
        Aders [ i ].zero();
    }

    wders.resize(fsd + 3);
    wders.zero();

    if ( fsd == 1 ) {
        // calculate values and derivatives of nonrational Bspline curve with weights at first (Aders, wders)
        uind = span(0) - degree [ 0 ];
        ind = uind + 1;
        for ( k = 0; k <= degree [ 0 ]; k++ ) {
            vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind + k);
            w = vertexCoordsPtr->at(nsd+1);

            Aders [ 0 ](0) += ders [ 0 ](0, k) * vertexCoordsPtr->at(1) * w;   // xw=sum(Nu*x*w)
            Aders [ 1 ](0) += ders [ 0 ](0, k) * vertexCoordsPtr->at(2) * w;   // yw=sum(Nu*y*w)
            Aders [ 2 ](0) += ders [ 0 ](0, k) * vertexCoordsPtr->at(3) * w;   // yw=sum(Nu*z*w)
            wders(0)    += ders [ 0 ](0, k) * w;                               // w=sum(Nu*w)

            Aders [ 0 ](1) += ders [ 0 ](1, k) * vertexCoordsPtr->at(1) * w;   // dxw/du=sum(dNu/du*x*w)
            Aders [ 1 ](1) += ders [ 0 ](1, k) * vertexCoordsPtr->at(2) * w;   // dyw/du=sum(dNu/du*y*w)
            Aders [ 2 ](1) += ders [ 0 ](1, k) * vertexCoordsPtr->at(3) * w;   // dzw/du=sum(dNu/du*z*w)
            wders(1)    += ders [ 0 ](1, k) * w;                               // dw/du=sum(dNu/du*w)

            Aders [ 0 ](2) += ders [ 0 ](2, k) * vertexCoordsPtr->at(1) * w;   // ddxw/du=sum(ddNu/du*x*w)
            Aders [ 1 ](2) += ders [ 0 ](2, k) * vertexCoordsPtr->at(2) * w;   // ddyw/du=sum(ddNu/du*y*w)
            Aders [ 2 ](2) += ders [ 0 ](2, k) * vertexCoordsPtr->at(3) * w;   // ddzw/du=sum(ddNu/du*z*w)
            wders(2)    += ders [ 0 ](2, k) * w;                               // ddw/du=sum(ddNu/du*w)

            Aders [ 0 ](3) += ders [ 0 ](3, k) * vertexCoordsPtr->at(1) * w;   // ddxw/du=sum(ddNu/du*x*w)
            Aders [ 1 ](3) += ders [ 0 ](3, k) * vertexCoordsPtr->at(2) * w;   // ddyw/du=sum(ddNu/du*y*w)
            Aders [ 2 ](3) += ders [ 0 ](3, k) * vertexCoordsPtr->at(3) * w;   // ddzw/du=sum(ddNu/du*z*w)
            wders(3)    += ders [ 0 ](3, k) * w;                               // ddw/du=sum(ddNu/du*w)
        }

        weight = wders(0);

        // calculation of jacobian matrix according to Eq 4.7
        jacobian(0, 0) = ( Aders [ 0 ](1) - wders(1) * Aders [ 0 ](0) / weight ) / weight; // dx/du
	jacobian(1,1) =  ( Aders [ 1 ](1) - wders(1) * Aders [ 1 ](0) / weight ) / weight; // dy/du
	jacobian(2,2) =  ( Aders [ 2 ](1) - wders(1) * Aders [ 2 ](0) / weight ) / weight; // dz/du

	Jacob = sqrt(jacobian(0,0)*jacobian(0,0) + jacobian(1,1)*jacobian(1,1) + jacobian(2,2)*jacobian(2,2));
        //Jacob = jacobian.giveDeterminant();

        //calculation of derivatives of NURBS basis functions with respect to local parameters is not covered by NURBS book
        product = weight * weight;
        cnt = 0;
        ind = uind + 1;
        for ( k = 0; k <= degree [ 0 ]; k++ ) {
            w = cellgeo.giveVertexCoordinates(ind + k)->at(nsd+1);
            // [dNu/du*w*sum(Nu*w) - Nu*w*sum(dNu/du*w)] / [J*sum(Nu*w)^2]
            answer(cnt, 0) = (ders [ 0 ](1, k) * w * weight - ders [ 0 ](0, k) * w * wders(1)) / product; 
            answer(cnt, 1) = (ders [ 0 ](2, k) * w * weight - ( 2 * answer(cnt, 0) * weight * wders(1) + ders [ 0 ](0, k) * w * wders(2) )) / product;
            answer(cnt, 2) = (ders [ 0 ](3, k) * w * weight - ( 3 * answer(cnt, 1) * weight * wders(1) + 3 * answer(cnt, 0) * weight * wders(2) + ders [ 0 ](0, k) * w * wders(3) )) / product;

            cnt++;
        }
    } else {
        OOFEM_ERROR("evaldNdx not implemented for fsd = %d", fsd);
    }

 #ifndef HAVE_VARIABLE_ARRAY_SIZE
    delete [] Aders;
 #endif

#ifndef HAVE_VARIABLE_ARRAY_SIZE
    delete [] ders;
#endif
    return Jacob;
}


} // end namespace oofem
