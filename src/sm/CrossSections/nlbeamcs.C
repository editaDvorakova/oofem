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

#include "sm/CrossSections/nlbeamcs.h"
#include "sm/Elements/structuralelement.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "gausspoint.h"
#include "material.h"
#include "floatarray.h"
#include "verbose.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"
#include "mathfem.h"

namespace oofem {
REGISTER_CrossSection(NlBeamSimpleElasticCrossSection);
REGISTER_CrossSection(NlBeamHyperElasticCrossSection);
REGISTER_CrossSection(NlBeamSimpleElastoPlasticCrossSection);

void
NlBeamSimpleElasticCrossSection::initializeFrom(InputRecord &ir)
{
    CrossSection::initializeFrom(ir);
    IR_GIVE_FIELD(ir, Iy, "iy");
    IR_GIVE_FIELD(ir, A, "area");
    IR_GIVE_OPTIONAL_FIELD(ir, gamma, "gamma");
}


void
NlBeamSimpleElasticCrossSection::postInitialize()
{ }


void
NlBeamSimpleElasticCrossSection::saveContext(DataStream &stream, ContextMode mode)
{
    NlBeamCrossSection::saveContext(stream, mode);

    if (mode & CM_Definition) {
        if (!stream.write(Iy) || !stream.write(A) || !stream.write(gamma)) {
            THROW_CIOERR(CIO_IOERR);
        }
    }
}


void
NlBeamSimpleElasticCrossSection::restoreContext(DataStream &stream, ContextMode mode)
{
    NlBeamCrossSection::restoreContext(stream, mode);

    if (mode & CM_Definition) {
        if (!stream.read(Iy) || !stream.read(A) || !stream.read(gamma)) {
            THROW_CIOERR(CIO_IOERR);
        }
        if (Iy <= 0. || A <= 0. || gamma <= 0.) {
            THROW_CIOERRM(CIO_BADVERSION, "invalid nonlinear-beam cross-section definition");
        }
    }
}


double
NlBeamSimpleElasticCrossSection :: computeMomentFromCurvature(double kappa)
{
  return E*Iy*kappa;

}


double
NlBeamSimpleElasticCrossSection :: computeDerivativeOfMomentWrtCurvature(double kappa)
{
  return E*Iy;
}



double
NlBeamSimpleElasticCrossSection :: computeCurvatureFromMoment(double M)
{
  return M/(E*Iy);
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////



void
NlBeamHyperElasticCrossSection :: initializeFrom(InputRecord &ir)
{
    CrossSection::initializeFrom(ir);

    IR_GIVE_FIELD(ir, b, "b");
    IR_GIVE_FIELD(ir, h, "h");
    IR_GIVE_FIELD(ir, materialNumber, _IFT_NlBeamHyperElasticCrossSection_MaterialNumber);

    if (b <= 0.0 || h <= 0.0) {
        OOFEM_ERROR("NlBeamHyperElasticCrossSection: b and h must be positive");
    }
    //
    IR_GIVE_OPTIONAL_FIELD(ir, nHeightIP, "nheightip");

    if (nHeightIP < 1) {
        OOFEM_ERROR("NlBeamHyperElasticCrossSection: nheightip must be positive");
    }

    if (nHeightIP % 2 == 0) {
        nHeightIP += 1;
    }
    //
    IR_GIVE_OPTIONAL_FIELD(ir, lateralStrain, "lateralstrain");

}


HyperelasticMaterial1d *
NlBeamHyperElasticCrossSection :: giveHyperElasticMaterial() const
{
    auto *material = dynamic_cast< HyperelasticMaterial1d * >(this->giveDomain()->giveMaterial(this->materialNumber));
    if ( !material ) {
        OOFEM_ERROR("NlBeamHyperElasticCrossSection: material %d is not HyperelasticMaterial1d", this->materialNumber);
    }
    return material;
}


Material *
NlBeamHyperElasticCrossSection :: giveMaterial(IntegrationPoint *ip) const
{
    return this->giveHyperElasticMaterial();
}


std::tuple<double, double>
NlBeamHyperElasticCrossSection :: compute_N_and_M(double lambda, double kappa)
{
    auto *material = this->giveHyperElasticMaterial();

    // For a fixed-height Biot section the midpoint integration has a closed
    // form.  Besides being cheaper, evaluating it directly avoids tiny
    // spurious N-kappa and M-lambda couplings from summing many cancelling
    // floating-point terms.  This is the exact law represented by the same
    // nHeightIP midpoint fibers, not a continuum-limit approximation.
    if ( !this->lateralStrain && material->isBiot() ) {
        const double area = b * h;
        const double n = this->nHeightIP;
        const double midpointInertia = b * h * h * h
                                      * (1.0 - 1.0 / (n * n)) / 12.0;
        return std::make_tuple(material->giveE() * area * (lambda - 1.0),
                               material->giveE() * midpointInertia * kappa);
    }

    const double weight = b * h / this->nHeightIP;

    std::vector<double> zStar;
    this->computeSectionKinematics(lambda, kappa, zStar);

    double N = 0.0;
    double M = 0.0;
    for ( int k = 0; k < this->nHeightIP; ++k ) {
        const double lambdaK = lambda + kappa * zStar[k];
        const double P1 = material->giveFirstPKStress_fromStretch(lambdaK).at(1);
        N += P1;
        M += zStar[k] * P1;
    }
    return std::make_tuple(weight * N, weight * M);
}



FloatMatrixF<2,2>
NlBeamHyperElasticCrossSection :: compute_dN_and_dM(double lambda, double kappa)
{
    auto *material = this->giveHyperElasticMaterial();

    if ( !this->lateralStrain && material->isBiot() ) {
        const double area = b * h;
        const double n = this->nHeightIP;
        const double midpointInertia = b * h * h * h
                                      * (1.0 - 1.0 / (n * n)) / 12.0;
        return {material->giveE() * area, 0.0,
                0.0, material->giveE() * midpointInertia};
    }

    const double weight = b * h / this->nHeightIP;

    std::vector<double> zStar, dzStar_dLambda, dzStar_dKappa;
    this->computeSectionKinematics(lambda, kappa, zStar, &dzStar_dLambda, &dzStar_dKappa);

    FloatMatrixF<2,2> J;
    J.at(1,1) = 0.0;
    J.at(1,2) = 0.0;
    J.at(2,1) = 0.0;
    J.at(2,2) = 0.0;

    for ( int k = 0; k < this->nHeightIP; ++k ) {
        const double zs = zStar[k];
        const double dzdl = dzStar_dLambda[k];
        const double dzdk = dzStar_dKappa[k];
        const double lambdaK = lambda + kappa * zs;
	//
	auto P  = material->giveFirstPKStress_fromStretch(lambdaK).at(1);
        auto A  = material->givedPdF_fromStretch(lambdaK).at(1,1);
	//
        const double dlambdaK_dLambda = 1.0 + kappa * dzdl;
        const double dlambdaK_dKappa  = zs + kappa * dzdk;

        J.at(1,1) += A * dlambdaK_dLambda;
        J.at(1,2) += A * dlambdaK_dKappa;
        J.at(2,1) += zs * A * dlambdaK_dLambda + P * dzdl;
        J.at(2,2) += zs * A * dlambdaK_dKappa  + P * dzdk;
    }

    J.at(1,1) *= weight;
    J.at(1,2) *= weight;
    J.at(2,1) *= weight;
    J.at(2,2) *= weight;

    return J;
}

  void
NlBeamHyperElasticCrossSection::computeSectionKinematics(double lambda, double kappa, std::vector<double> &zStar,
                                                         std::vector<double> *dzStar_dLambda,
                                                         std::vector<double> *dzStar_dKappa) const
{
    auto *material = this->giveHyperElasticMaterial();
  /*
    double b, h;
    this->giveEquivalentRectangle(b, h);
  */
    const int n = this->nHeightIP;
    const int m = (n - 1) / 2;
    const double dz = h / n;

    zStar.assign(n, 0.0);
    if ( dzStar_dLambda ) {
        dzStar_dLambda->assign(n, 0.0);
    }
    if ( dzStar_dKappa ) {
        dzStar_dKappa->assign(n, 0.0);
    }

    if ( !this->lateralStrain ) {
        for ( int k = 0; k < n; ++k ) {
            zStar[k] = ( (2.0 * (k + 1.0) - 1.0 - n) / (2.0 * n) ) * h;
        }
        return;
    }

    // Section 6.2, equations (273)--(276), with zero-based indexing.
    zStar[m] = 0.0;

    for ( int k = m + 1; k < n; ++k ) {
        const double lamPrev = lambda + kappa * zStar[k - 1];
        zStar[k] = zStar[k - 1] + dz * material->giveLateralStretch(lamPrev);
    }

    for ( int k = m - 1; k >= 0; --k ) {
        const double lamNext = lambda + kappa * zStar[k + 1];
        zStar[k] = zStar[k + 1] - dz * material->giveLateralStretch(lamNext);
    }

    if ( dzStar_dLambda == nullptr && dzStar_dKappa == nullptr ) {
        return;
    }

    auto &zl = *dzStar_dLambda;
    auto &zk = *dzStar_dKappa;
    zl[m] = 0.0;
    zk[m] = 0.0;

    // Section 6.3, differentiated positive-side recursion.
    for ( int k = m + 1; k < n; ++k ) {
        const int p = k - 1;
        const double lamPrev = lambda + kappa * zStar[p];
        const double Fp = material->giveLateralStretchDerivative(lamPrev);
        zl[k] = zl[p] + dz * Fp * (1.0 + kappa * zl[p]);
        zk[k] = zk[p] + dz * Fp * (zStar[p] + kappa * zk[p]);
    }

    // Section 6.3, differentiated negative-side recursion.
    for ( int k = m - 1; k >= 0; --k ) {
        const int p = k + 1;
        const double lamNext = lambda + kappa * zStar[p];
        const double Fp = material->giveLateralStretchDerivative(lamNext);
        zl[k] = zl[p] - dz * Fp * (1.0 + kappa * zl[p]);
        zk[k] = zk[p] - dz * Fp * (zStar[p] + kappa * zk[p]);
    }
}
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////



void
NlBeamSimpleElastoPlasticCrossSection::initializeFrom(InputRecord &ir)
{
    NlBeamSimpleElasticCrossSection::initializeFrom(ir);
    IR_GIVE_FIELD(ir, M_elastic, "mel");
    IR_GIVE_FIELD(ir, M_plastic, "mpl");
    IR_GIVE_OPTIONAL_FIELD(ir, delta, "delta");

    if (M_elastic <= 0. || M_plastic <= M_elastic) {
        throw ValueInputException(ir, "mpl", "requires 0 < mel < mpl");
    }
    if (delta <= 0. || delta > 1.) {
        throw ValueInputException(ir, "delta", "must lie in (0, 1]");
    }
}


void
NlBeamSimpleElastoPlasticCrossSection::saveContext(DataStream &stream, ContextMode mode)
{
    NlBeamSimpleElasticCrossSection::saveContext(stream, mode);

    if (mode & CM_Definition) {
        if (!stream.write(M_elastic) || !stream.write(M_plastic) || !stream.write(delta)) {
            THROW_CIOERR(CIO_IOERR);
        }
    }
}


void
NlBeamSimpleElastoPlasticCrossSection::restoreContext(DataStream &stream, ContextMode mode)
{
    NlBeamSimpleElasticCrossSection::restoreContext(stream, mode);

    if (mode & CM_Definition) {
        if (!stream.read(M_elastic) || !stream.read(M_plastic) || !stream.read(delta)) {
            THROW_CIOERR(CIO_IOERR);
        }
        if (M_elastic <= 0. || M_plastic <= M_elastic || delta <= 0. || delta > 1.) {
            THROW_CIOERRM(CIO_BADVERSION, "invalid nonlinear-beam elastoplastic definition");
        }
    }
}


void
NlBeamSimpleElastoPlasticCrossSection::updateMomentCurvatureParameters()
{
    if (E <= 0. || Iy <= 0.) {
        OOFEM_ERROR("NlBeamSimpleElastoPlasticCrossSection: E and Iy must be positive");
    }
    kappa_elastic = M_elastic/E/Iy;
    delta_kappa = 2. * (M_plastic-M_elastic)/E/Iy/M_PI;
    kappa_nl = kappa_elastic + delta_kappa * sqrt(1./delta-1.);
    M_nl = M_elastic + E *Iy * delta_kappa * atan(sqrt(1./delta-1.));
}

void
NlBeamSimpleElastoPlasticCrossSection::postInitialize()
{
    this->updateMomentCurvatureParameters();
}




double
NlBeamSimpleElastoPlasticCrossSection :: computeMomentFromCurvature(double kappa)
{
  this->updateMomentCurvatureParameters();
  double M;
  if(fabs(kappa) <= kappa_elastic) {
    M = E * Iy * kappa;
  } else if(fabs(kappa) <= kappa_nl) {
    M = signum(kappa) * (M_elastic + E * Iy * delta_kappa * atan((fabs(kappa) - kappa_elastic)/delta_kappa));
  } else {
    M = signum(kappa) * (M_nl + delta * E * Iy * (fabs(kappa) - kappa_nl));
  }
  return M;
}


double
NlBeamSimpleElastoPlasticCrossSection :: computeDerivativeOfMomentWrtCurvature(double kappa)
{
  this->updateMomentCurvatureParameters();
    double dMdk;

  if(fabs(kappa) <= kappa_elastic) {
    dMdk = E * Iy;
  } else if(fabs(kappa) <= kappa_nl) {
    dMdk = E * Iy / ( 1 + (fabs(kappa)-kappa_elastic) / delta_kappa * (fabs(kappa)-kappa_elastic)/delta_kappa );
  } else {
    dMdk = delta * E * Iy;
  }
  return dMdk;
}



double
NlBeamSimpleElastoPlasticCrossSection :: computeCurvatureFromMoment(double M)
{
  this->updateMomentCurvatureParameters();
  double kappa;
  if(fabs(M) <= M_elastic) {
    kappa = M / (E * Iy);
  } else if(fabs(M) <= M_nl) {
    kappa = signum(M) * (kappa_elastic + delta_kappa * tan((fabs(M) - M_elastic) / (E*Iy*delta_kappa)));
  } else {
    kappa = signum(M) * (kappa_nl + (fabs(M) - M_nl) / (E*Iy*delta));
  }
  return kappa;

}








} // end namespace oofem
