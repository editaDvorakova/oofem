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
 *               Copyright (C) 1993 - 2020   Borek Patzak
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

#include "hyperelasticmaterial1d.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "mathfem.h"
#include "contextioerr.h"
#include "datastream.h"



namespace oofem {
REGISTER_Material(HyperelasticMaterial1d);

HyperelasticMaterial1d::HyperelasticMaterial1d(int n, Domain *d) : StructuralMaterial(n, d)
{ }

FloatArrayF< 1 >
HyperelasticMaterial1d::giveFirstPKStressVector_1d(const FloatArrayF< 1 > &vF, GaussPoint *gp, TimeStep *tStep) const
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    FloatArrayF<1> vP;
    if(this->hyperelasticMaterialType == HEMT_Biot) {
      vP.at(1) = this->E * (vF.at(1) - 1.);
    } else if(this->hyperelasticMaterialType == HEMT_StVenantKirchhoff) {
      vP.at(1) = vF.at(1) * this->E * 0.5 * (vF.at(1) * vF.at(1) - 1.);
    } else if(this->hyperelasticMaterialType == HEMT_IncompressibleNeoHookean) {
      const double lambda = vF.at(1);
      vP.at(1) = (this->E / 3.0) * (lambda - 1.0 / (lambda * lambda));
    } else if(this->hyperelasticMaterialType == HEMT_Logarithmic) {
      const double lambda = vF.at(1);
      vP.at(1) = this->E * log(lambda) / lambda;
    } else {
      OOFEM_ERROR("Unknow material type");
    }
    // update gp
    status->letTempFVectorBe(vF);
    status->letTempPVectorBe(vP);
    //
    return vP;
}



FloatMatrixF< 1, 1 >
HyperelasticMaterial1d::give1dStressStiffnessMatrix_dPdF(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
// returns the 9x9 tangent stiffness matrix - dP/dF
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    double e = this->E;
    if(this->hyperelasticMaterialType == HEMT_Biot) {
      return {e};
    } else if(this->hyperelasticMaterialType == HEMT_StVenantKirchhoff) {
      FloatArray vF(status->giveTempFVector() );
      FloatMatrix d(1,1);
      d.at(1,1)= e * 0.5 * (vF.at(1) * vF.at(1) - 1.) +  e * vF.at(1)*vF.at(1);
      return d;
    } else if(this->hyperelasticMaterialType == HEMT_IncompressibleNeoHookean) {
      const double lambda = status->giveTempFVector().at(1);
      return { (e / 3.0) * (1.0 + 2.0 / (lambda * lambda * lambda)) };
    } else if(this->hyperelasticMaterialType == HEMT_Logarithmic) {
      const double lambda = status->giveTempFVector().at(1);
      return { e * (1.0 - log(lambda)) / (lambda * lambda) };
    } else {
      OOFEM_ERROR("Unknow material type");
    }


}


FloatArrayF< 1 >
HyperelasticMaterial1d::giveFirstPKStress_fromStretch(double lambda) const
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet
{
    FloatArrayF<1> vP;
    if(this->hyperelasticMaterialType == HEMT_Biot) {
      vP.at(1) = this->E * (lambda - 1.);
    } else if(this->hyperelasticMaterialType == HEMT_StVenantKirchhoff) {
      vP.at(1) = lambda * this->E * 0.5 * (lambda * lambda - 1.);
    } else if(this->hyperelasticMaterialType == HEMT_IncompressibleNeoHookean) {
      vP.at(1) = (this->E / 3.0) * (lambda - 1.0 / (lambda * lambda));
    } else if(this->hyperelasticMaterialType == HEMT_Logarithmic) {
      vP.at(1) = this->E * log(lambda) / lambda;
    } else {
      OOFEM_ERROR("Unknow material type");
    }
    return vP;
}



FloatMatrixF< 1, 1 >
HyperelasticMaterial1d::givedPdF_fromStretch(double lambda) const
// returns the 9x9 tangent stiffness matrix - dP/dF
{

    double e = this->E;
    if(this->hyperelasticMaterialType == HEMT_Biot) {
      return {e};
    } else if(this->hyperelasticMaterialType == HEMT_StVenantKirchhoff) {
      FloatMatrix d(1,1);
      d.at(1,1)= e * 0.5 * (lambda * lambda - 1.) +  e * lambda * lambda;
      return d;
    } else if(this->hyperelasticMaterialType == HEMT_IncompressibleNeoHookean) {
      return { (e / 3.0) * (1.0 + 2.0 / (lambda * lambda * lambda)) };
    } else if(this->hyperelasticMaterialType == HEMT_Logarithmic) {
      return { e * (1.0 - log(lambda)) / (lambda * lambda) };
    } else {
      OOFEM_ERROR("Unknow material type");
    }


}




double
HyperelasticMaterial1d::giveLateralStretch(double lambda) const
{
    if(this->hyperelasticMaterialType == HEMT_Biot) {
      return 1.0 - nu * (lambda - 1.0);
    } else if(this->hyperelasticMaterialType == HEMT_StVenantKirchhoff) {
        const double arg = 1.0 - nu * (lambda * lambda - 1.0);
        if (arg < 0.0) {
            OOFEM_WARNING("HyperelasticMaterial1d::giveLateralStretch: negative argument in sqrt, arg=%e, lambda=%e, nu=%e", arg, lambda, nu);
            return 0.0;
        }
	return sqrt(arg);
    } else if(this->hyperelasticMaterialType == HEMT_IncompressibleNeoHookean) {
      return 1.0 / sqrt(lambda);
    } else if(this->hyperelasticMaterialType == HEMT_Logarithmic) {
      return pow(lambda, -nu);
    } else {
      OOFEM_ERROR("Unknow material type");
    }
}

double
HyperelasticMaterial1d::giveLateralStretchDerivative(double lambda) const
{

    if(this->hyperelasticMaterialType == HEMT_Biot) {
      return -nu;
    } else if(this->hyperelasticMaterialType == HEMT_StVenantKirchhoff) {
        const double ls = this->giveLateralStretch(lambda);
        if (fabs(ls) < 1.e-16) {
            OOFEM_WARNING("HyperelasticMaterial1d::giveLateralStretchDerivative: near-zero lateral stretch, ls=%e", ls);
            return 0.0;
        }
	return -nu * lambda / ls;
    } else if(this->hyperelasticMaterialType == HEMT_IncompressibleNeoHookean) {
      return -0.5 / (lambda * sqrt(lambda));
    } else if(this->hyperelasticMaterialType == HEMT_Logarithmic) {
      return -nu * pow(lambda, -nu - 1.0);
    } else {
      OOFEM_ERROR("Unknow material type");
    }

}




std::unique_ptr<MaterialStatus> 
HyperelasticMaterial1d::CreateStatus(GaussPoint *gp) const
{
  return std::make_unique<StructuralMaterialStatus>(gp);

}


void
HyperelasticMaterial1d::initializeFrom(InputRecord &ir)
{
    StructuralMaterial :: initializeFrom(ir);

    int hyperTypeRecord = 0; // default
    IR_GIVE_OPTIONAL_FIELD(ir, hyperTypeRecord, _IFT_HyperelasticMaterial1d_type);
    // specify the type of hyperelasticmaterial
    if ( hyperTypeRecord == 0 ) {
        this->hyperelasticMaterialType = HEMT_Biot;
    } else if ( hyperTypeRecord == 1 ) {
        this->hyperelasticMaterialType = HEMT_StVenantKirchhoff;
    } else if ( hyperTypeRecord == 2 ) {
        this->hyperelasticMaterialType = HEMT_IncompressibleNeoHookean;
    } else if ( hyperTypeRecord == 3 ) {
        this->hyperelasticMaterialType = HEMT_Logarithmic;
    } else {
      throw ValueInputException(ir, _IFT_HyperelasticMaterial1d_type, "Unknown hyperelastic material type");
    }

    IR_GIVE_FIELD(ir, this->E, _IFT_HyperelasticMaterial1d_E);
    IR_GIVE_OPTIONAL_FIELD(ir, this->nu, _IFT_HyperelasticMaterial1d_nu);

}


void
HyperelasticMaterial1d::saveContext(DataStream &stream, ContextMode mode)
{
    StructuralMaterial::saveContext(stream, mode);

    if (mode & CM_Definition) {
        const int materialType = static_cast<int>(hyperelasticMaterialType);
        if (!stream.write(E) || !stream.write(nu) || !stream.write(materialType)) {
            THROW_CIOERR(CIO_IOERR);
        }
    }
}


void
HyperelasticMaterial1d::restoreContext(DataStream &stream, ContextMode mode)
{
    StructuralMaterial::restoreContext(stream, mode);

    if (mode & CM_Definition) {
        int materialType;
        if (!stream.read(E) || !stream.read(nu) || !stream.read(materialType)) {
            THROW_CIOERR(CIO_IOERR);
        }
        if (E <= 0. || materialType < HEMT_Biot || materialType > HEMT_Logarithmic) {
            THROW_CIOERRM(CIO_BADVERSION, "invalid one-dimensional hyperelastic material definition");
        }
        hyperelasticMaterialType = static_cast<HyperElasticMaterialType>(materialType);
    }
}
} // end namespace oofem
