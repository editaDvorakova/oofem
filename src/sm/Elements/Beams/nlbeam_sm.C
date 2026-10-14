#include "../sm/Elements/Beams/nlbeam_sm.h"
#include "material.h"
#include "../sm/CrossSections/nlbeamcs.h"
#include "node.h"

#include "floatarray.h"
#include "floatmatrix.h"
#include "floatmatrixf.h"
#include "intarray.h"

#include "engngm.h"
#include "mathfem.h"
#include "classfactory.h"
#include "function.h"
#include "domain.h"

#include <limits>


#include "../sm/Materials/HyperelasticMaterials/hyperelasticmaterial1d.h"


namespace oofem {
REGISTER_Element(NlBeam_Reissner);
REGISTER_Element(NlBeam_Ziegler);
  REGISTER_Element(NlBeam_Hyperelastic);


NlBeam_SM :: NlBeam_SM(int n, Domain *aDomain) : NLStructuralElement(n, aDomain)
{
    numberOfDofMans = 2;
    this->internalForces.resize(3);
    this->fb_init.resize(3);

}

void
NlBeam_SM :: updateYourself(TimeStep *tStep)
// Updates the receiver at end of step.
{
  StructuralElement :: updateYourself(tStep);
  this->fb_init = this->internalForces;
  this->committedKinematics = this->shootingKinematics;
}


void
NlBeam_SM :: initForNewStep()
// Updates the receiver at end of step.
{
  Element :: initForNewStep();
  this->internalForces =   this->fb_init;
  this->shootingKinematics = this->committedKinematics;
}


void
NlBeam_SM :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_w, R_v
    };
}




void
NlBeam_SM :: initializeFrom(InputRecord &ir, int priority)
{
    // first call parent
    NLStructuralElement :: initializeFrom(ir, priority);
    this->computeLength();
    // Numerical parameters
    // 1. number of segments for numerical integration along the beam, default value 100

    IR_GIVE_OPTIONAL_FIELD(ir, NIP, _IFT_NlBeam_SM_NIP);
    /*
    IR_GIVE_FIELD(ir, EA, "ea");
    IR_GIVE_FIELD(ir, EI, "ei");
    */


    //    IR_GIVE_FIELD(ir, GAs, "gas");
    //
    IR_GIVE_OPTIONAL_FIELD(ir, this->num, _IFT_NlBeam_SM_num);
    this->numericalShooting = this->num;
    this->numericalStiffness = this->num;
    IR_GIVE_OPTIONAL_FIELD(ir, this->numericalShooting, _IFT_NlBeam_SM_numShooting);
    IR_GIVE_OPTIONAL_FIELD(ir, this->numericalStiffness, _IFT_NlBeam_SM_numStiffness);
    IR_GIVE_OPTIONAL_FIELD(ir, this->fdCheck, _IFT_NlBeam_SM_fdCheck);
    IR_GIVE_OPTIONAL_FIELD(ir, this->fdStep, _IFT_NlBeam_SM_fdStep);
    IR_GIVE_OPTIONAL_FIELD(ir, this->localLineSearch, _IFT_NlBeam_SM_localLineSearch);
    IR_GIVE_OPTIONAL_FIELD(ir, this->sectionLineSearch, _IFT_NlBeam_SM_sectionLineSearch);
    IR_GIVE_OPTIONAL_FIELD(ir, this->sym, _IFT_NlBeam_SM_sym);
    // relative tolerance for iterations at the beam level
    IR_GIVE_OPTIONAL_FIELD(ir, beam_tol, _IFT_NlBeam_SM_Beam_Tolerance);
    // maximum number of iterations at the beam level
    IR_GIVE_OPTIONAL_FIELD(ir, beam_maxit, _IFT_NlBeam_SM_Beam_MaxIteration);
    IR_GIVE_OPTIONAL_FIELD(ir, section_tol, _IFT_NlBeam_SM_Section_Tolerance);
    IR_GIVE_OPTIONAL_FIELD(ir, section_maxit, _IFT_NlBeam_SM_Section_MaxIteration);
    IR_GIVE_OPTIONAL_FIELD(ir, nsubsteps_init, _IFT_NlBeam_SMBeam_NumberInitialSubsteps);
    IR_GIVE_OPTIONAL_FIELD(ir, nsubsteps_max, _IFT_NlBeam_SMBeam_NumberMaxSubsteps);
    if ( nsubsteps_init < 1 ) {
        OOFEM_ERROR("NlBeam_SM: nsubsteps must be positive");
    }
    if ( nsubsteps_max < nsubsteps_init ) {
        OOFEM_ERROR("NlBeam_SM: maxsubsteps must not be smaller than nsubsteps");
    }
    if ( beam_tol <= 0.0 || beam_maxit < 1 || section_tol <= 0.0 || section_maxit < 1 || fdStep <= 0.0 ) {
        OOFEM_ERROR("NlBeam_SM: tolerances, FD step, and iteration limits must be positive");
    }
    // initialize initial position vector
    FloatArray a = {this->giveNode(1)->giveCoordinate(1), this->giveNode(1)->giveCoordinate(3), 0, this->giveNode(2)->giveCoordinate(1), this->giveNode(2)->giveCoordinate(3), 0};
    X = FloatArrayF<6>(a);
    this->shootingKinematics = X;
    this->committedKinematics = X;
    // set size of the fields
    this->x.resize(NIP+1);
    this->z.resize(NIP+1);
    this->phi.resize(NIP+1);
    //compute alpha
    this->alpha = this->computeAlpha();
    // Applied magnetic field
    double BaA = 0, mu0 = 1, theta = 0,rho = 0, Br = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, BaA, "baa");
    IR_GIVE_OPTIONAL_FIELD(ir, mu0, "mu0");
    IR_GIVE_OPTIONAL_FIELD(ir, theta, "theta");
    Btildax = (BaA/mu0)*cos(theta);
    Btildaz = (BaA/mu0)*sin(theta);
    // Residual magnetic field
    IR_GIVE_OPTIONAL_FIELD(ir, Br, "br");
    IR_GIVE_OPTIONAL_FIELD(ir, rho, "rho");
    auto Brx = Br*cos(rho);
    auto Brz = Br*sin(rho);
     IR_GIVE_OPTIONAL_FIELD(ir, Ba_ltf, "ba_ltf");
    //
    // calculate the loading -- change the name
    double fx = 0, fz = 0;
    ////
    IR_GIVE_OPTIONAL_FIELD(ir, fx, _IFT_NlBeam_SM_fx);
    IR_GIVE_OPTIONAL_FIELD(ir, fx_ltf, _IFT_NlBeam_SM_fxLTF);
    ////
    IR_GIVE_OPTIONAL_FIELD(ir, fz, _IFT_NlBeam_SM_fz);
    IR_GIVE_OPTIONAL_FIELD(ir, fz_ltf, _IFT_NlBeam_SM_fzLTF);
    ////
    M = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, M, _IFT_NlBeam_SM_m);
    IR_GIVE_OPTIONAL_FIELD(ir, m_ltf, _IFT_NlBeam_SM_mLTF);
    //
    Fx.resize(NIP+1);
    Fz.resize(NIP+1);
    FBrx.resize(NIP+1);
    FBrz.resize(NIP+1);
    this->computeLoadingVariables(fx, fz, Brx, Brz);
    ///
    this->vN.resize(NIP+1);
    this->vQ.resize(NIP+1);
    this->vM.resize(NIP+1);
    ///
    this->s.resize(NIP+1);
    this->x.resize(NIP+1);
    this->z.resize(NIP+1);
    //
    this->phi.resize(NIP+1);
    this->vKappa.resize(NIP+1);

}

  void
  NlBeam_SM :: postInitialize()
  {
    this->crossSection = dynamic_cast< NlBeamSimpleElasticCrossSection* >( this->giveCrossSection() );
    if(this->crossSection == nullptr) {
      OOFEM_ERROR("NlBeam_SM::Incompatible cross-section type");
    }
    auto m  = dynamic_cast< HyperelasticMaterial1d* >( this->giveMaterial() );
    if(m == nullptr) {
      OOFEM_ERROR("NlBeam_SM::Incompatible material type");
    }
    this->crossSection->setE(m->giveE());

    EA  = this->crossSection->giveEA();
    GAs = this->crossSection->giveGAs();
  }

double
NlBeam_SM :: computeLength()
// Returns the length of the receiver.
{
    double dx, dy;
    Node *nodeA, *nodeB;

    if ( beamLength == 0. ) {
        nodeA   = this->giveNode(1);
        nodeB   = this->giveNode(2);
        dx      = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
        dy      = nodeB->giveCoordinate(3) - nodeA->giveCoordinate(3);
        beamLength  = sqrt(dx * dx + dy * dy);
    }

    return beamLength;
}



double
NlBeam_SM :: computeAlpha()
{
    double xA, xB, yA, yB;
    Node *nodeA, *nodeB;
    nodeA  = this->giveNode(1);
    nodeB  = this->giveNode(2);
    xA     = nodeA->giveCoordinate(1);
    xB     = nodeB->giveCoordinate(1);
    yA     = nodeA->giveCoordinate(3);
    yB     = nodeB->giveCoordinate(3);
    auto    alpha  = atan2(yB - yA, xB - xA);
    return alpha;
}





void
NlBeam_SM :: computeLoadingVariables(double fx, double fz, double Brx, double Brz)
{
  // assuming only load with constant intesity
  double dxi = beamLength/NIP;
  Fx.at(1) = 0.5 * fx * dxi;
  Fz.at(1) = 0.5 * fz * dxi;

  FBrx.at(1) = 0.5 * Brx * dxi;
  FBrz.at(1) = 0.5 * Brz * dxi;

  for (int i=2; i <= NIP; i++) {
    Fx.at(i) = Fx.at(i-1) + fx * dxi;
    Fz.at(i) = Fz.at(i-1) + fz * dxi;

    FBrx.at(i) = Brx;//FBrx.at(i-1) + Brx * dxi;
    FBrz.at(i) = Brz;//FBrz.at(i-1) + Brz * dxi;
  }

  Fx.at(NIP+1) = Fx.at(NIP) + 0.5 * fx * dxi;
  Fz.at(NIP+1) = Fz.at(NIP) + 0.5 * fz * dxi;

  FBrx.at(NIP+1) = FBrx.at(NIP);// + 0.5 * Brx * dxi;
  FBrz.at(NIP+1) = FBrz.at(NIP);// + 0.5 * Brz * dxi;
}


double
NlBeam_SM :: computeLoadAt(FloatArray &Px, FloatArray &Pz, double &tBax, double &tBaz, TimeStep *tStep)
{
  double m;
  if(fx_ltf != 0) {
    Px = Fx * this->giveDomain()->giveFunction(fx_ltf)->evaluate(tStep, VM_Total);
  } else {
    Px = Fx;
  }

  if(fz_ltf != 0) {
    Pz = Fz * this->giveDomain()->giveFunction(fz_ltf)->evaluate(tStep, VM_Total);
  } else {
    Pz = Fz;
  }

  if(m_ltf != 0) {
    m = M * this->giveDomain()->giveFunction(m_ltf)->evaluate(tStep, VM_Total);
  } else {
    m = M;
  }

  if(Ba_ltf != 0) {
    tBax = Btildax * this->giveDomain()->giveFunction(Ba_ltf)->evaluate(tStep, VM_Total);
    tBaz = Btildaz * this->giveDomain()->giveFunction(Ba_ltf)->evaluate(tStep, VM_Total);
  } else {
    tBax = Btildax;
    tBaz = Btildaz;
  }


  return m;

}


double
NlBeam_SM :: computeShootingError(const FloatArray &residual) const
{
    const double inverseLength = 1.0 / max(this->beamLength, 1.e-12);
    return sqrt(pow(residual.at(1) * inverseLength, 2)
                + pow(residual.at(2) * inverseLength, 2)
                + residual.at(3) * residual.at(3));
}


double
NlBeam_SM ::performLineSearch(const FloatArray &fa, const FloatArray &dforces, const FloatArray &xa, const  FloatArray &xb_target, TimeStep *tStep)
{
    auto [xbBase, GBase, GrBase, MpBase, dMpBase] = this->integrateAlongBeam(fa, xa, tStep);
    const double baseError = this->computeShootingError(xb_target - xbBase);
    double alpha = 1.0;
    for ( int iteration = 0; iteration < 16; ++iteration ) {
        const FloatArray trialForces = fa + alpha * dforces;
        auto [xbTrial, GTrial, GrTrial, MpTrial, dMpTrial] = this->integrateAlongBeam(trialForces, xa, tStep);
        const double trialError = this->computeShootingError(xb_target - xbTrial);
        if ( std::isfinite(trialError) && trialError < baseError ) {
            return alpha;
        }
        alpha *= 0.5;
    }
    return 0.0;
}




 /*
 Find forces and moment at the left end that lead to given displacements and rotations
 at both ends.
  The input displacements and output forces are taken with respect to the global coordinate system.
 Note that the transformation matrix T is affected by angle alpha that specifies the initial beam geometry.
 */
std::tuple<bool,FloatArrayF<3>, double>
 NlBeam_SM :: findLeftEndForces(const FloatArrayF<6> &x, TimeStep *tStep)
 {
   FloatArray dforces;
   int iter = 0;
   const double tolerance = beam_tol;
   FloatArray fa;
   fa = this->internalForces;
   FloatArray xa = x[{0,1,2}];
   FloatArray xb_target = x[{3,4,5}];

   // Use the end kinematics to select the appropriate branch when the
   // hyperelastic section law is non-convex (notably SVK in compression).
   // For a constant-stretch circular arc these estimates are exact.
   const double deltaPhi = xb_target.at(3) - xa.at(3);
   const double halfDeltaPhi = 0.5 * deltaPhi;
   const double sincHalf = fabs(halfDeltaPhi) < 1.e-8
                         ? 1.0 - halfDeltaPhi * halfDeltaPhi / 6.0
                         : sin(halfDeltaPhi) / halfDeltaPhi;
   const double chord = hypot(xb_target.at(1) - xa.at(1), xb_target.at(2) - xa.at(2));
   const double lambdaGuess = chord / (beamLength * fabs(sincHalf));
   if (lambdaGuess > 0.0 && lambdaGuess == lambdaGuess) {
     this->sectionLambdaGuess = lambdaGuess;
   }
   this->sectionKappaGuess = deltaPhi / beamLength;

   auto [xb_shoot, G, Gr, Mp, dMp] = this->integrateAlongBeam(fa, xa, tStep);
   double Mpr = Mp;

   auto computeNumericalJacobian = [this, &xa, tStep](
       const FloatArray &forces, const FloatArrayF<3> &
   ) {
       FloatMatrix numerical(3, 3);
       for ( int i = 1; i <= 3; ++i ) {
           FloatArray plusForces = forces;
           FloatArray minusForces = forces;
           // The section solve has a finite residual floor.  A central
           // difference with a resolvable force increment is substantially
           // less noisy than the former 1.e-7 forward difference, especially
           // in the initially uncoupled transverse directions.
           const double perturbation = this->fdStep * max(1.0, fabs(forces.at(i)));
           plusForces.at(i) += perturbation;
           minusForces.at(i) -= perturbation;
           auto [plusEnd, plusG, plusGr, plusMp, plusDMp]
               = this->integrateAlongBeam(plusForces, xa, tStep);
           auto [minusEnd, minusG, minusGr, minusMp, minusDMp]
               = this->integrateAlongBeam(minusForces, xa, tStep);
           for ( int j = 1; j <= 3; ++j ) {
               numerical.at(j, i) = (plusEnd.at(j) - minusEnd.at(j)) / (2.0 * perturbation);
           }
       }
       this->integrateAlongBeam(forces, xa, tStep);
       return numerical;
   };

   auto checkOrReplaceJacobian = [this, tStep, &fa, &G, &computeNumericalJacobian, &xb_shoot](int localIteration) {
       if ( !this->fdCheck && !this->numericalShooting ) {
           return;
       }
       const FloatMatrix numericalG = computeNumericalJacobian(fa, xb_shoot);
       if ( this->fdCheck ) {
           FloatMatrix difference(G);
           difference.subtract(numericalG);
           const double fdNorm = numericalG.computeFrobeniusNorm();
           double maxAbsError = 0.0;
           for ( int row = 1; row <= 3; ++row ) {
               for ( int col = 1; col <= 3; ++col ) {
                   maxAbsError = max(maxAbsError, fabs(difference.at(row, col)));
               }
           }
           fprintf(stderr,
                   "NLBEAM_FD_G element=%d time=%.16e local_iteration=%d relative_error=%.16e max_abs_error=%.16e\n",
                   this->giveNumber(), tStep->giveTargetTime(), localIteration,
                   difference.computeFrobeniusNorm() / max(fdNorm, 1.e-30), maxAbsError);
           for ( int row = 1; row <= 3; ++row ) {
               fprintf(stderr,
                       "NLBEAM_FD_G_ROW row=%d analytical=(%.16e,%.16e,%.16e) fd=(%.16e,%.16e,%.16e)\n",
                       row, G.at(row, 1), G.at(row, 2), G.at(row, 3),
                       numericalG.at(row, 1), numericalG.at(row, 2), numericalG.at(row, 3));
           }
       }
       if ( this->numericalShooting ) {
           G = numericalG;
       }
   };

   checkOrReplaceJacobian(0);
   bool finiteJacobian = true;
   for ( int row = 1; row <= 3; ++row ) {
       for ( int column = 1; column <= 3; ++column ) {
           finiteJacobian = finiteJacobian && std::isfinite(G.at(row, column));
       }
   }
   if ( !finiteJacobian ) {
       return {false, FloatArrayF<3>(this->internalForces), 0.0};
   }

   ////
   FloatMatrix jacobi(G);
   auto res = xb_target-xb_shoot;
   double error = this->computeShootingError(res);
   //
   while (error > tolerance && iter < beam_maxit && std::isfinite(error) ) {
     iter++;
     jacobi.solveForRhs(res, dforces);
     const double alpha = this->localLineSearch
                        ? performLineSearch(fa, dforces, xa, xb_target, tStep)
                        : 1.0;
     if ( alpha == 0.0 ) {
       return {false, FloatArrayF<3>(this->internalForces), 0.0};
     }
     fa += alpha * dforces;
     std::tie(xb_shoot, G, Gr, Mp, dMp) = this->integrateAlongBeam(fa, xa, tStep);
     Mpr = Mp;
     checkOrReplaceJacobian(iter);
     finiteJacobian = true;
     for ( int row = 1; row <= 3; ++row ) {
       for ( int column = 1; column <= 3; ++column ) {
         finiteJacobian = finiteJacobian && std::isfinite(G.at(row, column));
       }
     }
     if ( !finiteJacobian ) {
       return {false, FloatArrayF<3>(this->internalForces), 0.0};
     }
     jacobi = G;
     res = xb_target - xb_shoot;
     error = this->computeShootingError(res);
   }

   if (error > tolerance || !std::isfinite(error)) {
     // what to do ???
     return {false, FloatArrayF<3>( this->internalForces), 0};
   }
   return {true, FloatArrayF<3>(fa), Mpr};

 }

 /*
 Find end forces and moments that lead to given displacements and rotations
 at both ends.
 The input displacements and output forces are taken with respect to the global coordinate system.
 */
void
NlBeam_SM :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
   // solution vector
   answer.resize(6);
   FloatArray u;
   this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, u);
   auto y = X + u;
   // only the first three entries of f are computed
   auto [converged, fa, Mp] = this->findLeftEndForces(y, tStep);

   if (!converged) {
     // A global Newton/line-search trial is allowed to be rejected, but calls
     // made while evaluating that trial have already changed the element's
     // working shooting state.  The only guaranteed accepted continuation
     // point is the state committed at the end of the previous load step.
     // Roll back before kinematic substepping so a rejected global trial does
     // not become the local continuation origin.
     this->internalForces = this->fb_init;
     this->shootingKinematics = this->committedKinematics;
     const FloatArray startKinematics(this->committedKinematics);
     const FloatArray kinematicIncrement = y - startKinematics;
     constexpr int refinementFactor = 4;
     double progress = 0.0;
     double increment = 1.0 / this->nsubsteps_init;
     const double minimumIncrement = 1.0 / this->nsubsteps_max;
     int acceptedSubsteps = 0;
     converged = false;

     while ( progress < 1.0 - 10.0 * std::numeric_limits<double>::epsilon() ) {
       const double trialProgress = min(1.0, progress + increment);
       const FloatArray yTrial = startKinematics + trialProgress * kinematicIncrement;
       auto [substepConverged, substepForces, substepMp]
           = this->findLeftEndForces(yTrial, tStep);
       if ( substepConverged ) {
         progress = trialProgress;
         acceptedSubsteps++;
         this->internalForces = substepForces;
         this->shootingKinematics = FloatArrayF<6>(yTrial);
         fa = substepForces;
         Mp = substepMp;
         converged = true;
       } else {
         const double refinedIncrement = increment / refinementFactor;
         if ( refinedIncrement < minimumIncrement ) {
           converged = false;
           break;
         }
         increment = refinedIncrement;
         converged = false;
       }
     }

     if ( !converged ) {
       this->internalForces = fb_init;
       this->shootingKinematics = this->committedKinematics;
       OOFEM_ERROR("NlBeam_SM element %d: local shooting failed after %d accepted substeps (maximum equivalent subdivision %d)",
                   this->giveNumber(), acceptedSubsteps, this->nsubsteps_max);
     }
   } else {
     this->internalForces = fa;
     this->shootingKinematics = FloatArrayF<6>(y);
   }


   // internal forces
   answer.at(1) = this->internalForces.at(1);
   answer.at(2) = this->internalForces.at(2);
   answer.at(3) = this->internalForces.at(3);
   //
   double c1 =  y.at(5) - y.at(2);
   double c2 =  y.at(4) - y.at(1);
   FloatArray Px, Pz;
   double tBax, tBaz;
   this->computeLoadAt(Px, Pz, tBax, tBaz,  tStep);
     //
   answer.at(4) = -answer.at(1) - Px.at(NIP + 1);
   answer.at(5) = -answer.at(2) - Pz.at(NIP + 1);
   answer.at(6) = c1 * answer.at(1) - c2 * answer.at(2) - answer.at(3) + Mp;
}


 /*
 Evaluate the tangent stiffness matrix based on given displacements at both ends
 and given end forces at the left end (must be provided, which means that typically findLeftEndForces is run first).
 */
 void
 NlBeam_SM :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
 {
   answer.resize(6,6);
   // Stiffness assembly is not guaranteed to be preceded by residual
   // assembly.  In particular, at the first Newton iteration of a new load
   // step internalForces can still represent the committed displacement from
   // the previous step.  Synchronize the shooting solution before evaluating
   // G and the geometric terms of the tangent.
   FloatArray currentInternalForces;
   this->giveInternalForcesVector(currentInternalForces, tStep, 0);
   // solution vector
   FloatArray u;
   this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, u);
   auto y = X + u;
   FloatArray xa= FloatArray({y.at(1), y.at(2), y.at(3)});
   auto [xb_shoot, G, vGr, Mp, dMp] = this->integrateAlongBeam(this->internalForces, xa, tStep);
   FloatMatrix Ginv;
   Ginv.beInverseOf(G);
   //
   FloatMatrix Gr = {{1,0,0}, {0,1,0}, {vGr.at(1), vGr.at(2), vGr.at(3)}};
   auto K11 = -1. * ( Ginv * Gr );
   //ading sub matrix starting at 1 1
   answer.setSubMatrix(K11, 1, 1);
   //ading sub matrix starting at 1 4
   answer.setSubMatrix(Ginv, 1, 4);
   //ading sub vector starting at row 4, col 1
   FloatArray K1 = {answer.at(1,1),answer.at(1,2),answer.at(1,3),answer.at(1,4),answer.at(1,5),answer.at(1,6)};
   answer.addSubVectorRow(-1. * K1, 4, 1);
   //ading sub vector starting at row 5, col 1
   FloatArray K2 = {answer.at(2,1),answer.at(2,2),answer.at(2,3),answer.at(2,4),answer.at(2,5),answer.at(2,6)};
   answer.addSubVectorRow(-1. * K2, 5, 1);
   // construct the sixth row
   FloatArray K3 =  {answer.at(3,1),answer.at(3,2),answer.at(3,3),answer.at(3,4),answer.at(3,5),answer.at(3,6)};
   FloatArray K6 =   (y.at(5) - y.at(2) + dMp.at(1)) * K1 + (-y.at(4) + y.at(1) + dMp.at(2)) * K2 + (dMp.at(3)-1.) * K3;
   K6 += FloatArray({this->internalForces.at(2), -this->internalForces.at(1), 0, -this->internalForces.at(2), this->internalForces.at(1), 0});
   answer.addSubVectorRow(K6, 6, 1);
   //
   answer.at(6,3) += dMp.at(4);
   /*
   answer.at(6,1) =  answer.at(1,6);
   answer.at(6,2) =  answer.at(2,6);
   answer.at(6,3) =  answer.at(3,6);
   answer.at(6,4) =  answer.at(4,6);
   answer.at(6,5) =  answer.at(5,6);
   */

   if ( this->fdCheck || this->numericalStiffness ) {
     const FloatMatrix analyticalStiffness(answer);
     FloatMatrix numericalTangent;
     this->computeStiffnessMatrix_num(numericalTangent, rMode, tStep);
     if ( this->fdCheck ) {
       FloatMatrix difference(analyticalStiffness);
       difference.subtract(numericalTangent);
       const double fdNorm = numericalTangent.computeFrobeniusNorm();
       double maxAbsError = 0.0;
       int maxRow = 0, maxColumn = 0;
       for ( int row = 1; row <= 6; ++row ) {
         for ( int column = 1; column <= 6; ++column ) {
           const double componentError = fabs(difference.at(row, column));
           if ( componentError > maxAbsError ) {
             maxAbsError = componentError;
             maxRow = row;
             maxColumn = column;
           }
         }
       }
       fprintf(stderr,
               "NLBEAM_FD_K element=%d time=%.16e relative_error=%.16e max_abs_error=%.16e max_row=%d max_column=%d\n",
               this->giveNumber(), tStep->giveTargetTime(),
               difference.computeFrobeniusNorm() / max(fdNorm, 1.e-30),
               maxAbsError, maxRow, maxColumn);
       for ( int row = 1; row <= 6; ++row ) {
         fprintf(stderr,
                 "NLBEAM_FD_K_ROW row=%d analytical=(%.16e,%.16e,%.16e,%.16e,%.16e,%.16e) fd=(%.16e,%.16e,%.16e,%.16e,%.16e,%.16e)\n",
                 row,
                 analyticalStiffness.at(row, 1), analyticalStiffness.at(row, 2), analyticalStiffness.at(row, 3),
                 analyticalStiffness.at(row, 4), analyticalStiffness.at(row, 5), analyticalStiffness.at(row, 6),
                 numericalTangent.at(row, 1), numericalTangent.at(row, 2), numericalTangent.at(row, 3),
                 numericalTangent.at(row, 4), numericalTangent.at(row, 5), numericalTangent.at(row, 6));
       }
     }
     if ( this->numericalStiffness ) {
       answer = numericalTangent;
     }
   }

   if(sym) {
     FloatMatrix at;
     at.beTranspositionOf(answer);
     answer += at;
     answer *= 0.5;
   }
}



 void
NlBeam_SM :: computeStiffnessMatrix_num(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
 {
   (void)rMode;
   FloatArray u, iF, iFPlus, iFMinus;
   this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, u);

   const FloatArray savedInternalForces = this->internalForces;
   const FloatArrayF<6> savedShootingKinematics = this->shootingKinematics;
   const double savedLambdaGuess = this->sectionLambdaGuess;
   const double savedKappaGuess = this->sectionKappaGuess;
   auto restoreState = [this, &savedInternalForces, &savedShootingKinematics,
                        savedLambdaGuess, savedKappaGuess]() {
       this->internalForces = savedInternalForces;
       this->shootingKinematics = savedShootingKinematics;
       this->sectionLambdaGuess = savedLambdaGuess;
       this->sectionKappaGuess = savedKappaGuess;
   };

   this->giveInternalForcesVector_from_u(iF, tStep, u);
   answer.resize(6,6);

   for ( int i = 1; i <= 6; ++i ) {
     const double perturbation = this->fdStep * max(1.0, fabs(u.at(i)));
     FloatArray uPlus = u;
     FloatArray uMinus = u;
     uPlus.at(i) += perturbation;
     uMinus.at(i) -= perturbation;

     restoreState();
     this->giveInternalForcesVector_from_u(iFPlus, tStep, uPlus);
     restoreState();
     this->giveInternalForcesVector_from_u(iFMinus, tStep, uMinus);
     for ( int row = 1; row <= 6; ++row ) {
       answer.at(row, i) = (iFPlus.at(row) - iFMinus.at(row)) / (2.0 * perturbation);
     }
   }
   restoreState();
 }


void
NlBeam_SM :: giveInternalForcesVector_from_u(FloatArray &answer, TimeStep *tStep, const FloatArray &u)
{
   // solution vector
   answer.resize(6);
   auto y = X + u;
   // only the first three entries of f are computed
   auto [converged, fa, Mp] = this->findLeftEndForces(y, tStep);
   if ( !converged ) {
     OOFEM_ERROR("NlBeam_SM element %d: perturbed local shooting failed", this->giveNumber());
   }
   // internal forces
   answer.at(1) = fa.at(1);
   answer.at(2) = fa.at(2);
   answer.at(3) = fa.at(3);
   //
   double c1 =  y.at(5) - y.at(2);
   double c2 =  y.at(4) - y.at(1);

   FloatArray Px, Pz;
   double tBax, tBaz;
   this->computeLoadAt(Px, Pz, tBax, tBaz,  tStep);
     //
   answer.at(4) = -answer.at(1) - Px.at(NIP + 1);
   answer.at(5) = -answer.at(2) - Pz.at(NIP + 1);
   answer.at(6) = c1 * answer.at(1) - c2 * answer.at(2) - answer.at(3)+ Mp;
}





double
NlBeam_SM :: computeMomentFromCurvature(double kappa)
{

  return this->crossSection->computeMomentFromCurvature(kappa);
    //return EI*kappa;
}

double
NlBeam_SM :: computeDerivativeOfMomentWrtCurvature(double kappa)
{
  return this->crossSection->computeDerivativeOfMomentWrtCurvature(kappa);
  //return EI;
}

double
NlBeam_SM :: computeCurvatureFromMoment(double M)
{
  return this->crossSection->computeCurvatureFromMoment(M);
  //return M/EI;
}




Interface *NlBeam_SM :: giveInterface(InterfaceType it)
{
    switch ( it ) {
    case VTKXMLExportModuleElementInterfaceType:
      return static_cast< VTKXMLExportModuleElementInterface * >( this );
    default:
      return StructuralElement :: giveInterface(it);
    }
}

void
NlBeam_SM :: giveCompositeExportData(std :: vector< ExportRegion > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep )
{

    vtkPieces.resize(1);

    int numCells = this->NIP;
    const int numCellNodes  = 2; // linear line
    int nNodes = numCells * numCellNodes;
    //
    vtkPieces.at(0).setNumberOfCells(numCells);
    vtkPieces.at(0).setNumberOfNodes(nNodes);
    //
    int val    = 1;
    int offset = 0;
    IntArray nodes(numCellNodes);
    Node *nodeA = this->giveNode(1);
    int nodeNum = 1;
    double dxi = beamLength/NIP;

    FloatArray nodeCoords(3);
    IntArray connectivity(2);
    for ( int iElement = 1; iElement <= numCells; iElement++ ) {
      for (int iNode = 1; iNode <= numCellNodes; iNode++) {
	double L = (iElement-1) * dxi + (iNode -1) * dxi;
	nodeCoords.at(1) = nodeA->giveCoordinate(1) + L * cos(alpha);
	nodeCoords.at(3) = nodeA->giveCoordinate(3) + L * sin(alpha);
	nodeCoords.at(2) = 0;
	vtkPieces.at(0).setNodeCoords(nodeNum, nodeCoords);
	nodeNum++;
	connectivity.at(iNode) = val++;
      }
      vtkPieces.at(0).setConnectivity(iElement, connectivity);
      offset += 2;
      vtkPieces.at(0).setOffset(iElement, offset);
      vtkPieces.at(0).setCellType(iElement, 3);
    }


    int n = primaryVarsToExport.giveSize();
    vtkPieces [ 0 ].setNumberOfPrimaryVarsToExport(primaryVarsToExport, nNodes);
    for ( int i = 1; i <= n; i++ ) {
        UnknownType utype = ( UnknownType ) primaryVarsToExport.at(i);
        if ( utype == DisplacementVector ) {
	  for ( int nN = 1; nN <= nNodes; nN++ ) {
	    int lN = nN % 2;
	    int iNode;
	    if(lN == 0) {
	      iNode = nN/2 + 1;
	    } else {
	      iNode = (nN + 1) / 2;
	    }
	    double L = (iNode -1) * dxi;
	    FloatArray X(3);
	    X.at(1) = nodeA->giveCoordinate(1) + L * cos(alpha);
	    X.at(3) = nodeA->giveCoordinate(3) + L * sin(alpha);
	    //
	    FloatArray u;
	    u = {this->x.at(iNode) - X.at(1),0,this->z.at(iNode) - X.at(3)};
	    vtkPieces.at(0).setPrimaryVarInNode(utype, nN, u);
	  }
        }
    }

    /////////////////
    int nC = cellVarsToExport.giveSize();
    vtkPieces [ 0 ].setNumberOfCellVarsToExport(cellVarsToExport, numCells);
    for ( int i = 1; i <= nC; i++ ) {
      InternalStateType istype = ( InternalStateType ) cellVarsToExport.at(i);
      if ( istype == IST_BeamForceMomentTensor ) {
	for ( int iC = 1; iC <= numCells; iC++ ) {
	  FloatArray bfmt(6);
	  bfmt.at(1) = vN.at(iC);
	  bfmt.at(2) = vQ.at(iC);
	  bfmt.at(3) = vM.at(iC);
	  vtkPieces.at(0).setCellVar(istype, iC, bfmt);
	}
      } else if ( istype == IST_MaterialNumber ) {
	for (int iC = 1; iC <= numCells; iC++) {
	  vtkPieces.at(0).setCellVar(istype, iC, {1});
	}
      }
    }


}


void
NlBeam_SM :: doMatlabOutput(std::string fileName)
{

  FILE *FID;
    fileName += "_element" + std::to_string(this->giveNumber());
  auto temp = fileName;
  std :: string functionname;
  size_t backslash = temp.rfind("/");
  if (backslash != std :: string :: npos ) {
    functionname = temp.substr(backslash+1, std :: string :: npos);
  } else {
    functionname = temp;
  }

  //
  fileName += ".m";
  if ( ( FID = fopen(fileName.c_str(), "w") ) == NULL ) {
    OOFEM_ERROR("failed to open file %s", fileName.c_str() );
  }
  //
  fprintf( FID, " function [x z phi kappa N Q M]= %s \n", functionname.c_str() );
  //
  fprintf(FID, "x=[");
  for ( double val: this->x ) {
    fprintf( FID, "%5.8f,", val);
  }
  fprintf(FID, "];\n");
  //
  fprintf(FID, "z=[");
  for ( double val: this->z ) {
    fprintf( FID, "%5.8f,", val);
  }
  fprintf(FID, "];\n");
  //
  fprintf(FID, "phi=[");
  for ( double val: this->phi ) {
    fprintf( FID, "%5.8f,", val );
  }
  fprintf(FID, "];\n");
    //
  fprintf(FID, "kappa=[");
  for ( double val: this->vKappa ) {
    fprintf( FID, "%5.8f,", val );
  }
  fprintf(FID, "];\n");
  //

  fprintf(FID, "N=[");
  for ( double val: this->vN ) {
    fprintf( FID, "%5.8f,", val );
  }
  fprintf(FID, "];\n");
  //
  fprintf(FID, "Q=[");
  for ( double val: this->vQ ) {
    fprintf( FID, "%5.8f,", val );
  }
  fprintf(FID, "];\n");
  //
  fprintf(FID, "M=[");
  for ( double val: this->vM ) {
    fprintf( FID, "%5.8f,", val );
  }
  fprintf(FID, "];\n");
  //
  fprintf(FID, "end\n");
  //
  fclose(FID);
}






  /////////////////////////Reissner//////////////////////////////////////////

  std::tuple<FloatArrayF<3>,FloatMatrixF<3,3>, FloatArrayF<3>, double, FloatArrayF<4>>
NlBeam_Reissner :: integrateAlongBeam(const FloatArray &fa, const FloatArray &xa, TimeStep *tStep)
 {
   FloatArray Px, Pz;
   double tBax, tBaz;
   double m = this->computeLoadAt(Px, Pz, tBax, tBaz,  tStep);
   // initial configuration
   this->x.at(1) = xa.at(1);
   this->z.at(1) = xa.at(2);
   this->phi.at(1) = xa.at(3) - alpha;
   //
   double Xab = fa.at(1), Zab = fa.at(2), Mab = fa.at(3);
   // store variables for postprocessing
   this->s.at(1) = 0.;
   this->vM.at(1) = -Mab;
   //
   FloatArray dx(4), dz(4), dphi(4);
   dphi  = {0,0,0,1};
   //
   FloatArray dN(4), dQ(4), dM(4);
   dM = {0 , 0, -1, 0};
   //
   double dxi = beamLength/NIP;
   //
   double M = -Mab;
   double Mp = 0;
   FloatArray dMp(4);
   dMp.zero();
   double Mm = 0;
   FloatArray dMm(4);
   dMm.zero();
   for (int i=2; i <= NIP+1; i++) {
     this->s.at(i) = this->s.at(i-1) + dxi;
     //
     double kappa = this->computeCurvatureFromMoment(M);
     double dMdKappa = this->computeDerivativeOfMomentWrtCurvature(kappa);
     auto dkappa = dM / dMdKappa;
     //
     double phi_mid = this->phi.at(i-1) + 0.5 * kappa * dxi;
     auto dphi_mid = dphi + 0.5 * dkappa * dxi;
     //
     double c = cos(phi_mid);
     double s = sin(phi_mid);
     double N_mid = - c * (Xab + Px.at(i-1)) + s * (Zab + Pz.at(i-1)) + FBrx.at(i) * ( c * tBax - s * tBaz);
     double Q_mid = - s * (Xab + Px.at(i-1)) - c * (Zab + Pz.at(i-1)) + FBrx.at(i) * ( c * tBaz + s * tBax);
     // store variables for postprocessing
     this->vN.at(i-1) = N_mid;
     this->vQ.at(i-1) = Q_mid;
     //
     auto dN_mid = - Q_mid * dphi_mid;
     dN_mid.at(1) += - c;
     dN_mid.at(2) +=   s;
     //
     auto dQ_mid = N_mid * dphi_mid;
     dQ_mid.at(1) += - s;
     dQ_mid.at(2) += - c;
     //
     double lambda_mid = 1 + N_mid / EA;
     auto dLambda_mid = dN_mid / EA;
     //
     double gamma_mid = Q_mid / GAs;
     auto dGamma_mid = dQ_mid / GAs;
     //
     double delta_x = ( c * lambda_mid + s * gamma_mid ) * dxi;
     FloatArray ddx;
     ddx = ( c * dLambda_mid + s * dGamma_mid - ( s * lambda_mid - c * gamma_mid) *dphi_mid ) * dxi;
     //
     double delta_z = ( c * gamma_mid - s * lambda_mid ) * dxi;
     FloatArray ddz;
     ddz =( c * dGamma_mid - s * dLambda_mid - ( s * gamma_mid + c * lambda_mid) * dphi_mid ) * dxi;
     //
     this->x.at(i) = x.at(i-1) + delta_x;
     dx += ddx;
     //
     this->z.at(i) = z.at(i-1) + delta_z;
     dz += ddz;
     //
     Mp += - m * dxi + Px.at(i-1) * delta_z - Pz.at(i-1) * delta_x;
     dMp += Px.at(i-1) * ddz - Pz.at(i-1) * ddx;
     //
     Mm +=  FBrx.at(i) * ( tBaz * delta_x - tBax*delta_z ) + FBrz.at(i) * ( s * tBaz - c * tBax ) * dxi;
     dMm += FBrx.at(i) * ( tBaz * ddx - tBax * ddz ) + FBrz.at(i) * ( c * tBaz + s * tBax ) * dxi * dphi_mid;
     //
     M = -Mab + Xab *  (z.at(i) - z.at(1)) - Zab * (x.at(i) - x.at(1)) + Mp + Mm;
     //
     dM = dMp + dMm + FloatArray({z.at(i)-z.at(1) + Xab * dz.at(1) - Zab * dx.at(1), -x.at(i) + x.at(1) + Xab * dz.at(2) - Zab * dx.at(2), -1. + Xab * dz.at(3) - Zab * dx.at(3), Xab * dz.at(4) - Zab * dx.at(4)});
     //
     //
     kappa = this->computeCurvatureFromMoment(M);
     dMdKappa = this->computeDerivativeOfMomentWrtCurvature(kappa);
     dkappa = dM / dMdKappa;
     this->phi.at(i) = phi_mid + 0.5 * kappa * dxi;
     dphi = dphi_mid + 0.5 * dxi * dkappa;
     // store variables for postprocessing
     this->vKappa.at(i) = kappa;
     this->vM.at(i) = M;
     //

   }
   // right end displacements
   FloatArrayF<3> xb;
   xb.at(1) = x.at(NIP+1); xb.at(2) = z.at(NIP+1); xb.at(3) = phi.at(NIP+1) + alpha;
   // jacobi matrix G
   FloatMatrix G(3,3);
   G.copySubVectorRow(FloatArray({dx.at(1), dx.at(2), dx.at(3)}), 1,1);
   G.copySubVectorRow(FloatArray({dz.at(1), dz.at(2), dz.at(3)}), 2,1);
   G.copySubVectorRow(FloatArray({dphi.at(1), dphi.at(2), dphi.at(3)}), 3,1);
   // 3rd column of Gr
   FloatArrayF<3> Gr;
   Gr.at(1) = dx.at(4);
   Gr.at(2) = dz.at(4);
   Gr.at(3) = dphi.at(4);
   return {xb, FloatMatrixF<3,3> (G), Gr, Mp+Mm, FloatArrayF<4> (dMp+dMm)};

 }





  /////////////////////////Ziegler//////////////////////////////////////////
  std::tuple<FloatArrayF<3>,FloatMatrixF<3,3>, FloatArrayF<3>, double, FloatArrayF<4>>
NlBeam_Ziegler :: integrateAlongBeam(const FloatArray &fa, const FloatArray &xa, TimeStep *tStep)
 {
   FloatArray Px, Pz;
   double tBax, tBaz;
   double m = this->computeLoadAt(Px, Pz, tBax, tBaz,  tStep);
   // initial configuration
   this->x.at(1) = xa.at(1);
   this->z.at(1) = xa.at(2);
   this->phi.at(1) = xa.at(3) - alpha;
   //
   double Xab = fa.at(1), Zab = fa.at(2), Mab = fa.at(3);
   // store variables for postprocessing
   this->s.at(1) = 0.;
   this->vM.at(1) = -Mab;
   //
   FloatArray dx(4), dz(4), dphi(4);
   dphi  = {0,0,0,1};
   //
   FloatArray dN(4), dQ(4), dM(4);
   dM = {0 , 0, -1, 0};
   //
   double dxi = beamLength/NIP;
   //
   double M = -Mab;
   double Mp = 0;
   FloatArray dMp(4);
   dMp.zero();
   double Mm = 0;
   FloatArray dMm(4);
   dMm.zero();
   //
   auto chi_mid = 0.;
   for (int i=2; i <= NIP+1; i++) {
     this->s.at(i) = this->s.at(i-1) + dxi;
     //
     double kappa = this->computeCurvatureFromMoment(M);
     double dMdKappa = this->computeDerivativeOfMomentWrtCurvature(kappa);
     auto dkappa = dM / dMdKappa;
     //
     double phi_mid = phi.at(i-1) + 0.5 * kappa * dxi;
     auto dphi_mid = dphi + 0.5 * dkappa * dxi;
     //
     auto iter = 0;
     double c, s, N_mid, Q_mid, lambda_mid;
     while (iter <= chi_iter) {
       iter++;
       c = cos(phi_mid - chi_mid);
       s = sin(phi_mid - chi_mid);
       N_mid = - c * (Xab + Px.at(i-1)) + s * (Zab + Pz.at(i-1));
       Q_mid = - s * (Xab + Px.at(i-1)) - c * (Zab + Pz.at(i-1));
       //
       lambda_mid = 1 + N_mid / EA;
       //
       auto F = GAs * chi_mid - lambda_mid * Q_mid;
       if(fabs(F) < chi_tol) {
	 break;
       }
       auto dF = GAs + lambda_mid * N_mid - Q_mid * Q_mid / EA;
       chi_mid -= F/dF;
       if(iter == chi_iter){
	 //@todo: cut the step...
	 OOFEM_ERROR("Maximum number of iteraton for calculation of chi reached");
       }
     }
     // store variables for postprocessing
     this->vN.at(i-1) = N_mid;
     this->vQ.at(i-1) = Q_mid;
     // linearization
     auto denom  = (GAs + lambda_mid * N_mid - Q_mid * Q_mid / EA);
     auto dchi_mid = (lambda_mid * N_mid - Q_mid * Q_mid / EA) * dphi_mid;
     dchi_mid.at(1) += (-c * Q_mid / EA - s * lambda_mid);
     dchi_mid.at(2) += (s * Q_mid / EA - c * lambda_mid);
     dchi_mid /= denom;
     //
     auto dN_mid = - Q_mid * (dphi_mid - dchi_mid);
     dN_mid.at(1) += - c;
     dN_mid.at(2) +=   s;
     // not needed ?
     /*     auto dQ_mid = N_mid * (dphi_mid - dchi_mid);
     dQ_mid.at(1) += - s;
     dQ_mid.at(2) += - c;
     */
     //
     auto dLambda_mid = dN_mid / EA;
     //
     double delta_x = c * lambda_mid  * dxi;
     FloatArray ddx;
     ddx = ( c * dLambda_mid - s * lambda_mid * ( dphi_mid - dchi_mid) ) * dxi;
     //
     double delta_z = - s * lambda_mid * dxi;
     FloatArray ddz;
     ddz = (- s * dLambda_mid - c * lambda_mid * ( dphi_mid - dchi_mid ) ) * dxi;
     //
     this->x.at(i) = this->x.at(i-1) + delta_x;
     dx += ddx;
     //
     this->z.at(i) = this->z.at(i-1) + delta_z;
     dz += ddz;
     //
     Mp += - m * dxi + Px.at(i-1) * delta_z - Pz.at(i-1) * delta_x;
     dMp += Px.at(i-1) * ddz - Pz.at(i-1) * ddx;
     //
     M = -Mab + Xab *  (z.at(i) - z.at(1)) - Zab * (x.at(i) - x.at(1)) + Mp;
     //
     dM = dMp + FloatArray({z.at(i)-z.at(1) + Xab * dz.at(1) - Zab * dx.at(1), -x.at(i) + x.at(1) + Xab * dz.at(2) - Zab * dx.at(2), -1. + Xab * dz.at(3) - Zab * dx.at(3), Xab * dz.at(4) - Zab * dx.at(4)});
     //
     //
     kappa = this->computeCurvatureFromMoment(M);
     dMdKappa = this->computeDerivativeOfMomentWrtCurvature(kappa);
     dkappa = dM / dMdKappa;
     this->phi.at(i) = phi_mid + 0.5 * kappa * dxi;
     dphi = dphi_mid + 0.5 * dxi * dkappa;
     // store variables for postprocessing
     vM.at(i) = M;
     vKappa.at(i) = kappa;
   }
   // right end displacements
   FloatArrayF<3> xb;
   xb.at(1) = x.at(NIP+1); xb.at(2) = z.at(NIP+1); xb.at(3) = phi.at(NIP+1) + alpha;
   // jacobi matrix G
   FloatMatrix G(3,3);
   G.copySubVectorRow(FloatArray({dx.at(1), dx.at(2), dx.at(3)}), 1,1);
   G.copySubVectorRow(FloatArray({dz.at(1), dz.at(2), dz.at(3)}), 2,1);
   G.copySubVectorRow(FloatArray({dphi.at(1), dphi.at(2), dphi.at(3)}), 3,1);
   // 3rd column of Gr
   FloatArrayF<3> Gr;
   Gr.at(1) = dx.at(4);
   Gr.at(2) = dz.at(4);
   Gr.at(3) = dphi.at(4);
   return {xb, FloatMatrixF<3,3> (G), Gr, Mp+Mm, FloatArrayF<4> (dMp+dMm)};
 }

 /////////////////////////Kirchhoff Hyperelastic model//////////////////////////////////////////
/* std::tuple<FloatArrayF<3>,FloatMatrixF<3,3>, FloatArrayF<3>, double, FloatArrayF<4>>
NlBeam_Hyperelastic :: integrateAlongBeam(const FloatArray &fa, const FloatArray &xa, TimeStep *tStep)
{
    // initial configuration
    this->x.at(1) = xa.at(1);
    this->z.at(1) = xa.at(2);
    this->phi.at(1) = xa.at(3) - alpha;
    //
    double Xab = fa.at(1), Zab = fa.at(2), Mab = fa.at(3);
    //
    FloatArray kappa(NIP);
    FloatArray M(NIP);
    FloatArray lambda(NIP);
    // @todo: update
    FloatArray dkappa(NIP);
    FloatArray dlambda(NIP);
    FloatArray dM_mid(NIP);
    //
    double dxi = beamLength/NIP;
    // Moment at the left end
    M.at(1) = -Mab;
    // compute lambda_0 and kappa_0, initial guess lambda_0 = 1 and \kappa_0 = 0
    auto [l, k, J] = this->compute_stretch_and_curvature(1,  0,  0, Xab, Zab, M.at(1), dxi);
    // ---- sensitivities (dx, dz, dphi) wrt 3 unknowns at the left end ----
    FloatArray dphi(3), dx(3), dz(3), dM(3);
    // Choose appropriate initial sensitivities
    dphi = {0,0,1};
    //
    dx.zero();
    dz.zero();
    // Initial dM
    dM.at(1) =  this->z.at(1) + Xab * dz.at(1) - Zab * dx.at(1);
    dM.at(2) = -this->x.at(1) + Xab * dz.at(2) - Zab * dx.at(2);
    dM.at(3) = -1.0           + Xab * dz.at(3) - Zab * dx.at(3);
    //
    FloatArray dlambda_i, dlambda_i_1;
    FloatArray dphi_i, dphi_i_1;
    // ---- main integration loop: i = 2...NIP
    for (int i = 2; i <= NIP; i++) {
	double phi_mid = phi.at(i-1) + kappa.at(i-1) * 0.5 * dxi;
	auto dphi_mid = dphi_i_1 + 0.5 * dkappa  * 0.5 * dxi;
	//
	auto M_mid = M.at(i-1) - (Xab*sin(phi(i-1)) + Zab*cos(phi(i-1))) * lambda.at(i-1) * 0.5 * dxi;
	//
        dM_mid= dM - ( Xab * cos(phi.at(i-1)) * dphi * lambda.at(i-1) - Zab * sin(phi.at(i-1)) * dphi * lambda.at(i-1) + Xab * sin(phi.at(i-1)) * dlambda + Zab * cos(phi.at(i-1)) * dlambda ) * 0.5 * dxi;
	dM_mid(1) += - sin(phi.at(i)) * lambda.at(i) * 0.5 * dxi;
	dM_mid(2) += - cos(phi.at(i)) * lambda.at(i) * 0.5 * dxi;
	//
	//	auto [l, k, J] = compute_stretch_and_curvature(lambda.at(i-1), kappa.at(i-1), phi_mid, M_mid, Xab, Zab);
        //
	phi.at(i) = phi_mid + kappa.at(i) * 0.5 * dxi;
	dphi = dphi_mid + dkappa  * 0.5 * dxi;
	//
        double delta_x = ( lambda.at(i-1) * cos(phi.at(i-1)) + lambda.at(i) * cos(phi.at(i)) ) * 0.5 * dxi;
        this->x.at(i+1) = this->x.at(i) + delta_x;
        //
        FloatArray ddx;
	//        ddx = (    dlambda_i_1  * cos(phi.at(i-1))
	//	+  dlambda_i_1  * cos(phi.at(i))
	//	-   lambda.at(i-1)     * sin(phi.at(i-1))   * dphi.at(i-2)
        //        -   lambda.at(i)       * sin(phi.at(i)) * dphi.at(i-2) ) * 0.5 * dxi;

        dx += ddx;
        //
        double delta_z = - ( lambda.at(i-1) * sin(phi.at(i-1)) + lambda.at(i) * sin(phi.at(i)) ) * 0.5 * dxi;
        this->z.at(i+1) = this->z.at(i) + delta_z;
        //
        FloatArray ddz;

        //ddz = - ( dlambda_i    * sin(phi_i)
        //        + dlambda_ip1  * sin(phi_ip1)
        //        + lambda_i     * cos(phi_i)   * dphi_i
        //        + lambda_ip1   * cos(phi_ip1) * dphi ) * 0.5 * dxi;

        dz += ddz;
        //
        M = -Mab + Xab * (this->z.at(i)-this->z.at(1)) - Zab * (this->x.at(i)-this->x.at(1));
        //
        dM.at(1) =  this->z.at(i) + Xab * dz.at(1) - Zab * dx.at(1);
        dM.at(2) = -this->x.at(i) + Xab * dz.at(2) - Zab * dx.at(2);
        dM.at(3) = -1.0              + Xab * dz.at(3) - Zab * dx.at(3);
    }

    // ---- right end displacements xb ----
    FloatArrayF<3> xb;
    xb.at(1) = this->x.at(NIP+1);
    xb.at(2) = this->z.at(NIP+1);
    xb.at(3) = this->phi.at(NIP+1);

    // ---- Jacobian G = [dx(:,end), dz(:,end), dphi(:,end)]' ----
    FloatMatrix G(3,3);
    G.copySubVectorRow(FloatArray({dx.at(1), dx.at(2), dx.at(3)}),   1, 1);
    G.copySubVectorRow(FloatArray({dz.at(1), dz.at(2), dz.at(3)}),   2, 1);
    G.copySubVectorRow(FloatArray({dphi.at(1), dphi.at(2), dphi.at(3)}), 3, 1);

    return { xb, FloatMatrixF<3,3>(G), FloatArrayF<3>(), 0., FloatArrayF<4> () };
}
*/

/*
std::tuple<double, double, FloatMatrixF<2,2>>
NlBeam_Hyperelastic :: compute_stretch_and_curvature(double l, double k, double Xab, double Zab, double M_mid, double phi_mid, double dxi)
{
  auto lambda = l;
  auto kappa = k;
  int iter = 0;
  int maxiter = 10;
  //
  double FN, FM;
  //  auto [FN,FM] = this->cs->compute_N_M(lambda,kappa);
  auto res = {FN, FM-M_mid};
  double error = norm(res);
  double tolerance = 1.e-6;
  FloatMatrix J;

  while ((error>tolerance)&&(iter < maxiter)) {
    iter = iter + 1;
    //auto [FN,FM] = this->cs->compute_N_and_M(lambda,kappa);
    auto phi = phi_mid + kappa * 0.5 * dxi;
    //
    FN = FN + (Xab * cos(phi) - Zab * sin(phi));
    FM = FM + lambda * 0.5 * dxi * (Xab * sin(phi) + Zab * cos(phi) );
    //
    res = {FN, FM - M_mid};
    error = norm(res);
    //
    //J = this->crossSection->compute_dN_dM(lambda, kappa);
    J.at(1,2) +=  - (Xab * sin(phi) + Zab * cos(phi) ) * 0.5 * dxi;
    J.at(2,1) +=  + (Xab * sin(phi) + Zab * cos(phi) ) * 0.5 * dxi;
    J.at(2,2) +=  + (Xab * cos(phi) - Zab * sin(phi) ) * lambda * dxi * dxi * 0.25;

    auto dlambda = - ( FN * J.at(2,2) - ( FM - M_mid ) * J.at(1,2) ) / J.giveDeterminant();
    auto dkappa  = - (-FN * J.at(2,1) + ( FM - M_mid ) * J.at(1,1) ) / J.giveDeterminant();

    lambda += dlambda;
    kappa  += dkappa;
  }

  return std::make_tuple( lambda, kappa, J );
}
*/





 /////////////////////////Kirchhoff Hyperelastic model//////////////////////////////////////////
 std::tuple<FloatArrayF<3>,FloatMatrixF<3,3>, FloatArrayF<3>, double, FloatArrayF<4>>
NlBeam_Hyperelastic :: integrateAlongBeam(const FloatArray &fa, const FloatArray &xa, TimeStep *tStep)
{
    (void)tStep;
    // initial configuration (local beam coordinates)
    this->x.at(1) = xa.at(1);
    this->z.at(1) = xa.at(2);
    this->phi.at(1) = xa.at(3) - alpha;

    const double Xab = fa.at(1);
    const double Zab = fa.at(2);
    const double Mab = fa.at(3);
    const double dxi = beamLength / NIP;

    // Section variables at all trapezoidal points (left + NIP updates)
    FloatArray lambda(NIP + 1), kappa(NIP + 1), Msec(NIP + 1);

    // 1) Initial values: phi0 = phia, x0 = xa, z0 = za, M0 = -Mab
    Msec.at(1) = -Mab;

    // 2) Solve for lambda0 and kappa0 using the modified left-end equations.
    auto [lambda0, kappa0, J0] = this->compute_stretch_and_curvature(this->sectionLambdaGuess, this->sectionKappaGuess,
                                                                     Xab, Zab, Msec.at(1), this->phi.at(1), 0.0);
    lambda.at(1) = lambda0;
    kappa.at(1) = kappa0;

    // Sensitivities with respect to [Xab, Zab, Mab]
    FloatArray dx(4), dz(4);
    dx.zero();
    dz.zero();
    FloatArray dphi = {0.0, 0.0, 0.0, 1.0};
    FloatArray dM = {0.0, 0.0, -1.0, 0.0};

    FloatArray dlambda_prev(4), dkappa_prev(4);
    {
        const double D = J0.at(1,1) * J0.at(2,2) - J0.at(1,2) * J0.at(2,1);
        const double phi_eval = this->phi.at(1); // because dxi = 0 in the initial solve
        const double c = cos(phi_eval);
        const double s = sin(phi_eval);

        FloatArray dRNs(4), dRMs(4);
        dRNs.zero();
        dRMs.zero();
        dRNs.at(1) = c;
        dRNs.at(2) = -s;
        dRNs.at(3) = 0.0;
        // The initial section equilibrium is also a function of the
        // prescribed left-end rotation.  This fourth sensitivity seeds Gr,
        // i.e. the third column of d(x_b)/d(x_a).  Omitting it makes the
        // analytical stiffness use an incomplete base-rotation derivative
        // whenever the resolved transverse end force is nonzero.
        dRNs.at(4) = -Xab * s - Zab * c;

        dRMs.at(3) = 1.0; // -dM0 with dM0 = -dMab

	///////////////
	FloatArray dlambda_temp = (dRMs * J0.at(1,2) - dRNs * J0.at(2,2)) / D;
	//
	FloatArray dkappa_temp = (dRNs * J0.at(2,1) - dRMs * J0.at(1,1)) / D;
	//
	/*
	FloatArray tmp1, tmp2, tmp3, dlambda_temp;
	tmp1 = dRMs * J0.at(1,2);
	tmp2 = dRNs * J0.at(2,2);
	tmp3 = tmp1 - tmp2;
	dlambda_temp = tmp3 / D;


	//        FloatArray dlambda_temp = (dRMs * J0.at(1,2) - dRNs * J0.at(2,2)) / D;
        FloatArray dkappa_temp;
	dkappa_temp  = (dRNs * J0.at(2,1) - dRMs * J0.at(1,1)) / D;
	*/
        dlambda_prev = dlambda_temp;
        dkappa_prev = dkappa_temp;
    }

    // 3) Trapezoidal marching along the beam
    for (int i = 2; i <= NIP + 1; ++i) {
        const int ip = i - 1;

        const double phi_mid = this->phi.at(ip) + 0.5 * dxi * kappa.at(ip);
        const FloatArray dphi_mid = dphi + 0.5 * dxi * dkappa_prev;

        const double force_term_prev = Xab * sin(this->phi.at(ip)) + Zab * cos(this->phi.at(ip));
        const double M_mid = Msec.at(ip) - force_term_prev * 0.5 * dxi * lambda.at(ip);

        FloatArray dM_mid = dM;
        dM_mid -= ((Xab * cos(this->phi.at(ip)) - Zab * sin(this->phi.at(ip))) * lambda.at(ip) * 0.5 * dxi) * dphi;
        dM_mid -= (force_term_prev * 0.5 * dxi) * dlambda_prev;
        dM_mid.at(1) += -sin(this->phi.at(ip)) * lambda.at(ip) * 0.5 * dxi;
        dM_mid.at(2) += -cos(this->phi.at(ip)) * lambda.at(ip) * 0.5 * dxi;

        auto [lambda_i, kappa_i, Ji] = this->compute_stretch_and_curvature(lambda.at(ip), kappa.at(ip), Xab, Zab, M_mid, phi_mid, dxi);
        lambda.at(i) = lambda_i;
        kappa.at(i) = kappa_i;

        const double theta = phi_mid + 0.5 * dxi * kappa_i;
        const double c = cos(theta);
        const double s = sin(theta);

        const double FNphi = -Xab * sin(theta) - Zab * cos(theta);
        const double FNX = c;
        const double FNZ = -s;
        const double FMphi = lambda_i * 0.5 * dxi * (Xab * c - Zab * s);
        const double FMX = lambda_i * 0.5 * dxi * s;
        const double FMZ = lambda_i * 0.5 * dxi * c;
        const double D = Ji.at(1,1) * Ji.at(2,2) - Ji.at(1,2) * Ji.at(2,1);

        FloatArray dRNs(3), dRMs(3);
        dRNs = FNphi * dphi_mid;
        dRNs.at(1) += FNX;
        dRNs.at(2) += FNZ;

        dRMs = FMphi * dphi_mid - dM_mid;
        dRMs.at(1) += FMX;
        dRMs.at(2) += FMZ;

        FloatArray dlambda_i = (dRMs * Ji.at(1,2) - dRNs * Ji.at(2,2)) / D;
        FloatArray dkappa_i  = (dRNs * Ji.at(2,1) - dRMs * Ji.at(1,1)) / D;

        this->phi.at(i) = phi_mid + 0.5 * dxi * kappa_i;
        dphi = dphi_mid + 0.5 * dxi * dkappa_i;

        const double delta_x = (lambda.at(ip) * cos(this->phi.at(ip)) + lambda.at(i) * cos(this->phi.at(i))) * 0.5 * dxi;
        this->x.at(i) = this->x.at(ip) + delta_x;
        FloatArray ddx = (cos(this->phi.at(ip)) * dlambda_prev
                        + cos(this->phi.at(i)) * dlambda_i
                        - lambda.at(ip) * sin(this->phi.at(ip)) * (dphi_mid - 0.5 * dxi * dkappa_prev)
                        - lambda.at(i) * sin(this->phi.at(i)) * dphi) * 0.5 * dxi;
        dx += ddx;

        const double delta_z = -(lambda.at(ip) * sin(this->phi.at(ip)) + lambda.at(i) * sin(this->phi.at(i))) * 0.5 * dxi;
        this->z.at(i) = this->z.at(ip) + delta_z;
        FloatArray ddz = -1.0 * (sin(this->phi.at(ip)) * dlambda_prev
                         + sin(this->phi.at(i)) * dlambda_i
                         + lambda.at(ip) * cos(this->phi.at(ip)) * (dphi_mid - 0.5 * dxi * dkappa_prev)
                         + lambda.at(i) * cos(this->phi.at(i)) * dphi) * 0.5 * dxi;
        dz += ddz;

        Msec.at(i) = -Mab + Xab * (this->z.at(i) - this->z.at(1)) - Zab * (this->x.at(i) - this->x.at(1));
        dM.at(1) = this->z.at(i) - this->z.at(1) + Xab * dz.at(1) - Zab * dx.at(1);
        dM.at(2) = -this->x.at(i) + this->x.at(1) + Xab * dz.at(2) - Zab * dx.at(2);
        dM.at(3) = -1.0 + Xab * dz.at(3) - Zab * dx.at(3);
	dM.at(4) =	  Xab * dz.at(4) - Zab * dx.at(4);

        dlambda_prev = dlambda_i;
        dkappa_prev = dkappa_i;
    }


    FloatArrayF<3> xb;
    xb.at(1) = this->x.at(NIP + 1);
    xb.at(2) = this->z.at(NIP + 1);
    xb.at(3) = this->phi.at(NIP + 1) + alpha;

    FloatMatrix G(3,3);
    G.copySubVectorRow(FloatArray({dx.at(1), dx.at(2), dx.at(3)}), 1, 1);
    G.copySubVectorRow(FloatArray({dz.at(1), dz.at(2), dz.at(3)}), 2, 1);
    G.copySubVectorRow(FloatArray({dphi.at(1), dphi.at(2), dphi.at(3)}), 3, 1);
   // 3rd column of Gr
   FloatArrayF<3> Gr;
   Gr.at(1) = dx.at(4);
   Gr.at(2) = dz.at(4);
   Gr.at(3) = dphi.at(4);
   //return {xb, FloatMatrixF<3,3> (G), Gr, Mp+Mm, FloatArrayF<4> (dMp+dMm)};
   return {xb, FloatMatrixF<3,3>(G), Gr, 0.0, FloatArrayF<4>()};

}

std::tuple<double, double, FloatMatrixF<2,2>>
NlBeam_Hyperelastic :: compute_stretch_and_curvature(double l, double k, double Xab, double Zab, double M_mid, double phi_mid, double dxi)
{
    double lambda = l;
    double kappa = k;
    const int maxiter = this->section_maxit;
    const double tolerance = this->section_tol;
    FloatMatrixF<2,2> J;

    for (int iter = 0; iter < maxiter; ++iter) {
        const double theta = phi_mid + 0.5 * dxi * kappa;
        const double c = cos(theta);
        const double s = sin(theta);

        auto nm = this->crossSection->compute_N_and_M(lambda, kappa);
        const double fN = std::get<0>(nm);
        const double fM = std::get<1>(nm);

        J = this->crossSection->compute_dN_and_dM(lambda, kappa);
	//
        const double RN = fN + Xab * c - Zab * s;
        const double RM = fM + lambda * 0.5 * dxi * (Xab * s + Zab * c) - M_mid;
	//@todo: get h from cross-section
        const double h = this->crossSection->giveDepth();
        const double residual = sqrt( h * h * RN * RN + RM * RM);
        if (residual < tolerance) {
            J.at(1,2) += -(Xab * s + Zab * c) * 0.5 * dxi;
            J.at(2,1) += +(Xab * s + Zab * c) * 0.5 * dxi;
            J.at(2,2) += +(Xab * c - Zab * s) * lambda * 0.25 * dxi * dxi;
            return std::make_tuple(lambda, kappa, J);
        }

        J.at(1,2) += -(Xab * s + Zab * c) * 0.5 * dxi;
        J.at(2,1) += +(Xab * s + Zab * c) * 0.5 * dxi;
        J.at(2,2) += +(Xab * c - Zab * s) * lambda * 0.25 * dxi * dxi;

        const double D = J.at(1,1) * J.at(2,2) - J.at(1,2) * J.at(2,1);
        if (fabs(D) < 1.e-16) {
            const double nan = std::numeric_limits<double>::quiet_NaN();
            return std::make_tuple(nan, nan, J);
        }

        const double dlambda = (J.at(1,2) * RM - J.at(2,2) * RN) / D;
        const double dkappa  = (J.at(2,1) * RN - J.at(1,1) * RM) / D;

        if ( !this->sectionLineSearch ) {
            lambda += dlambda;
            kappa += dkappa;
            if ( !std::isfinite(lambda) || !std::isfinite(kappa) ) {
                const double nan = std::numeric_limits<double>::quiet_NaN();
                return std::make_tuple(nan, nan, J);
            }
            continue;
        }

        double stepLength = 1.0;
        bool accepted = false;
        for ( int lineSearchIteration = 0; lineSearchIteration < 16; ++lineSearchIteration ) {
            const double trialLambda = lambda + stepLength * dlambda;
            const double trialKappa = kappa + stepLength * dkappa;
            const double trialTheta = phi_mid + 0.5 * dxi * trialKappa;
            const double trialC = cos(trialTheta);
            const double trialS = sin(trialTheta);
            const auto trialNm = this->crossSection->compute_N_and_M(trialLambda, trialKappa);
            const double trialRN = std::get<0>(trialNm) + Xab * trialC - Zab * trialS;
            const double trialRM = std::get<1>(trialNm)
                                 + trialLambda * 0.5 * dxi * (Xab * trialS + Zab * trialC)
                                 - M_mid;
            const double trialResidual = sqrt(h * h * trialRN * trialRN + trialRM * trialRM);
            if ( std::isfinite(trialResidual) && trialResidual < residual ) {
                lambda = trialLambda;
                kappa = trialKappa;
                accepted = true;
                break;
            }
            stepLength *= 0.5;
        }
        if ( !accepted ) {
            const double nan = std::numeric_limits<double>::quiet_NaN();
            return std::make_tuple(nan, nan, J);
        }
    }

    const double nan = std::numeric_limits<double>::quiet_NaN();
    return std::make_tuple(nan, nan, J);
}

void
NlBeam_Hyperelastic :: initializeFrom(InputRecord &ir, int priority)
{
    NlBeam_SM :: initializeFrom(ir, priority);
    IR_GIVE_OPTIONAL_FIELD(ir, b, "b");
    IR_GIVE_OPTIONAL_FIELD(ir, h, "h");
    IR_GIVE_OPTIONAL_FIELD(ir, nHeightIP, "nheightip");
    IR_GIVE_OPTIONAL_FIELD(ir, lateralStrain, "lateralstrain");
}

void
NlBeam_Hyperelastic :: postInitialize()
{
    auto cs = dynamic_cast< NlBeamHyperElasticCrossSection* >( this->giveCrossSection() );
    if(cs == nullptr) {
        OOFEM_ERROR("NlBeam_Hyperelastic: Incompatible cross-section type");
    }
    this->crossSection = cs;

    // The hyperelastic material is a cross-section property. Accept a
    // redundant, matching element material for compatibility, but never use
    // it to mutate a shared cross-section.
    if ( this->giveMaterialNumber() != 0 && this->giveMaterialNumber() != cs->giveMaterialNumber() ) {
        OOFEM_ERROR("NlBeam_Hyperelastic: element material %d differs from cross-section material %d",
                    this->giveMaterialNumber(), cs->giveMaterialNumber());
    }
    cs->giveHyperElasticMaterial();
}






} // end namespace oofem
