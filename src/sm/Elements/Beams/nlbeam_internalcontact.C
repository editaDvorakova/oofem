#include "../sm/Elements/Beams/nlbeam_internalcontact.h"
#include "../sm/Materials/HyperelasticMaterials/hyperelasticmaterial1d.h"
#include "material.h"
#include "crosssection.h"
#include "node.h"
#include "../sm/CrossSections/nlbeamcs.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "floatarray.h"
#include "engngm.h"
#include "mathfem.h"
#include "classfactory.h"
#include "contextioerr.h"
#include "datastream.h"

namespace oofem {
REGISTER_Element(NlBeamInternalContact);


  NlBeamInternalContact :: NlBeamInternalContact (int n, Domain *aDomain) : NlBeam_SM(n, aDomain), Jacobi(3,3), Jacobi44(4,4), Kblock(3,3)
{
  this->internalForces.resize(3);
  this->tempInternalForces.resize(3);
}


void
NlBeamInternalContact :: initializeFrom(InputRecord &ir, int priority)
{
    StructuralElement :: initializeFrom(ir, priority);

    IR_GIVE_OPTIONAL_FIELD(ir, TOL_BEAM, _IFT_NlBeam_InternalContact_Beam_Tolerance);
    IR_GIVE_OPTIONAL_FIELD(ir, MAXIT_BEAM, _IFT_NlBeam_InternalContact_Beam_MaxIteration);
    IR_GIVE_OPTIONAL_FIELD(ir, processLockSteps, _IFT_NlBeamInternalContact_ProcessLockSteps);

    IR_GIVE_FIELD(ir, leftSegmentLength, _IFT_NlBeamInternalContact_LeftSegmentLength);
    IR_GIVE_FIELD(ir, rightSegmentLength, _IFT_NlBeamInternalContact_RightSegmentLength);
    IR_GIVE_FIELD(ir, leftActiveSegmentLength, _IFT_NlBeamInternalContact_LeftActiveSegmentLength);
    IR_GIVE_FIELD(ir, rightActiveSegmentLength, _IFT_NlBeamInternalContact_RightActiveSegmentLength);

    int cmode;
    IR_GIVE_FIELD(ir, cmode, _IFT_NlBeamInternalContact_ContactMode);
    IR_GIVE_FIELD(ir, DX, _IFT_NlBeamInternalContact_dx);
    IR_GIVE_FIELD(ir, friction, _IFT_NlBeamInternalContact_Friction);

    if (TOL_BEAM <= 0.) {
        throw ValueInputException(ir, _IFT_NlBeam_InternalContact_Beam_Tolerance, "must be positive");
    }
    if (MAXIT_BEAM < 1) {
        throw ValueInputException(ir, _IFT_NlBeam_InternalContact_Beam_MaxIteration, "must be at least 1");
    }
    if (processLockSteps < 0) {
        throw ValueInputException(ir, _IFT_NlBeamInternalContact_ProcessLockSteps,
                                  "must not be negative");
    }
    if (DX <= 0.) {
        throw ValueInputException(ir, _IFT_NlBeamInternalContact_dx, "must be positive");
    }
    if (friction < 0.) {
        throw ValueInputException(ir, _IFT_NlBeamInternalContact_Friction, "must not be negative");
    }
    if (leftSegmentLength <= 0. || rightSegmentLength <= 0.) {
        throw ValueInputException(ir, _IFT_NlBeamInternalContact_LeftSegmentLength,
                                  "both segment lengths must be positive");
    }
    if (leftActiveSegmentLength < 0. || leftActiveSegmentLength > leftSegmentLength) {
        throw ValueInputException(ir, _IFT_NlBeamInternalContact_LeftActiveSegmentLength,
                                  "must lie between 0 and lsl");
    }
    if (rightActiveSegmentLength < 0. || rightActiveSegmentLength > rightSegmentLength) {
        throw ValueInputException(ir, _IFT_NlBeamInternalContact_RightActiveSegmentLength,
                                  "must lie between 0 and rsl");
    }
    if (cmode < CMT_N || cmode > CMT_CB) {
        throw ValueInputException(ir, _IFT_NlBeamInternalContact_ContactMode,
                                  "must be an integer from 0 to 8");
    }
    contactMode = static_cast<ContactModeType>(cmode);

    matNum = {1, 2};
    IR_GIVE_OPTIONAL_FIELD(ir, matNum, _IFT_NlBeamInternalContact_matNum);
    if (matNum.giveSize() != 2) {
        throw ValueInputException(ir, _IFT_NlBeamInternalContact_matNum,
                                  "must contain exactly two material numbers");
    }

    this->computeLength();
    this->givePitch();
    if (beamLength <= 0.) {
        throw ValueInputException(ir, "nodes", "must define a beam with positive length");
    }

    auto nodeA = this->giveNode(1);
    auto nodeB = this->giveNode(2);
    cosAlpha = (nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1)) / beamLength;
    sinAlpha = (nodeB->giveCoordinate(3) - nodeA->giveCoordinate(3)) / beamLength;

    trialContactMode = contactMode;
    trialLeftActiveSegmentLength = leftActiveSegmentLength;
    trialRightActiveSegmentLength = rightActiveSegmentLength;
    tempInternalForces = internalForces;
    previousStepDisplacements.clear();
    contactDataInitialized = false;
    this->resetIterationBranchState();
}


void
NlBeamInternalContact :: postInitialize()
{
  // Geometry is not part of the element context; rebuild it both after a
  // normal input initialization and after reconstruction from a restart.
  beamLength = 0.;
  pitch = 10.;
  this->computeLength();
  this->givePitch();
  if (beamLength <= 0.) {
    OOFEM_ERROR("NlBeamInternalContact: nodes must define a beam with positive length");
  }
  auto nodeA = this->giveNode(1);
  auto nodeB = this->giveNode(2);
  cosAlpha = (nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1)) / beamLength;
  sinAlpha = (nodeB->giveCoordinate(3) - nodeA->giveCoordinate(3)) / beamLength;

  this->crossSection = dynamic_cast< NlBeamCrossSection* >( this->giveCrossSection() );
  if(this->crossSection == nullptr) {
    OOFEM_ERROR("NlBeamInternalContact: incompatible cross-section type");
  }
  auto m  = dynamic_cast< HyperelasticMaterial1d* >( this->giveMaterial() );
  if(m == nullptr) {
      OOFEM_ERROR("NlBeamInternalContact: incompatible material type");
  }
  this->crossSection->setE(m->giveE());
  EA  = this->crossSection->giveEA();
  EI  = this->crossSection->giveEI();
  if (EA <= 0. || EI <= 0.) {
    OOFEM_ERROR("NlBeamInternalContact: cross-section must provide positive EA and EI");
  }
  contactDataInitialized = true;
}


void
NlBeamInternalContact :: ensureContactDataInitialized()
{
  if (!contactDataInitialized) {
    this->postInitialize();
  }
}


  

double
NlBeamInternalContact :: L2norm(double x, double y, double z)
{
  return sqrt(x*x+y*y+z*z);
}


// ============================================================================================
// FUNCTIONS RELATED TO THE CROSS-SECTION MOMENT-CURVATURE RELATION
// ============================================================================================
double
NlBeamInternalContact :: computeDerMomentFromCurvature(double kappa)
{
  return this->crossSection->computeDerivativeOfMomentWrtCurvature(kappa);
}

double
NlBeamInternalContact :: computeCurvatureFromMoment(double M)
{
  return this->crossSection->computeCurvatureFromMoment(M);
}
// ============================================================================================
// FUNCTIONS RELATED TO THE CONTACT MODEL WITH DRY COULOMB FRICTION
// ============================================================================================

double
NlBeamInternalContact :: evalContactLoadingFunction(double Nc, double Qc, ContactModeType cmode)
{
  double fc;  
  if (cmode == CMT_AA || cmode == CMT_AB || cmode == CMT_AC || cmode==CMT_CA) {
    // right segment above the left one
      fc = fabs(Nc) - friction*Qc; 
  } else {
    // right segment below the left one
      fc = fabs(Nc) + friction*Qc;
  }
  return fc;
}

double
NlBeamInternalContact :: evalDerContactLoadingFunction(double Nc, double dN, double dQ, ContactModeType cmode)
{
  // contribution of normal force increment
  double dfc = signum(Nc)*dN;
  // contribution of shear force increment
  if (cmode==CMT_AA || cmode==CMT_AB || cmode==CMT_AC || cmode==CMT_CA) {
    // right segment above the left one
    dfc -= friction*dQ; 
  } else { // right segment below the left one
      dfc += friction*dQ;
  }
  return dfc;
}


bool
NlBeamInternalContact :: integrateAlongSegment(FloatArray &fab, double deltaPhi, double Lb, FloatArray &ub, FloatMatrix &jac_b, double &Lc, FloatArray &uc, FloatMatrix &jac_c, bool inflection_only)
{
  
  bool inflection_detected = false;
  FloatMatrix jacobi(3,3), jacobi_prev(3,3), jac_s(3,3);
  // initialization at the left end
  double Xab = fab.at(1);
  double Zab = fab.at(2);
  double Mab = fab.at(3);
  double M = -Mab;
  FloatArray dM(3);
  dM.zero();
  dM.at(3) = -1;
  double kappa = computeCurvatureFromMoment(M);
  double dMdkappa = computeDerMomentFromCurvature(kappa);
  FloatArray dkappa;
  dkappa = dM;
  dkappa.times(1./dMdkappa);
  FloatArray u(3), u_prev;
  u_prev = ub;
 
  // basic loop over spatial steps  
  int nstep = ceil(Lb/DX); // the spatial step size is fixed and the segment length does not need to be its integer multiple
  double x = 0;
  for (int istep = 1; istep <= nstep; istep++){
    x += DX;
    //
    u_prev = u;
    //
    jacobi_prev = jacobi;
    // rotation at midstep and its derivatives with respect to the left-end forces
    double phi_mid = u_prev.at(3) + kappa * DX / 2.;
    FloatArray dphi_mid(3);
    for (int j = 1; j <= 3; j++)
      dphi_mid.at(j) = jacobi.at(3,j)+dkappa.at(j)*DX/2.;
    // normal force at midstep and its derivatives with respect to the left-end forces
    double N_mid = -Xab*cos(phi_mid)+Zab*sin(phi_mid);
    FloatArray dN_mid(3);
    dN_mid = dphi_mid;
    dN_mid.times(Xab*sin(phi_mid)+Zab*cos(phi_mid));
    dN_mid.at(1) -= cos(phi_mid);
    dN_mid.at(2) += sin(phi_mid);
    // horizontal displacement at the end of the step
    u.at(1) = u_prev.at(1) + DX * ((1.+N_mid/EA)*cos(phi_mid)-1.);
   // vertical displacement at the end of the step
    u.at(2) = u_prev.at(2) - DX * ( 1. + N_mid / EA ) * sin(phi_mid);
    //update first two rows of jacobi marix
    for (int j=1; j<=3; j++){
      jacobi.at(1,j) += DX*(dN_mid.at(j)/EA)*cos(phi_mid) - DX*(1.+N_mid/EA)*sin(phi_mid)*dphi_mid.at(j);
      jacobi.at(2,j) -= DX*(dN_mid.at(j)/EA)*sin(phi_mid) + DX*(1.+N_mid/EA)*cos(phi_mid)*dphi_mid.at(j);
    }

    
    // bending moment and curvature at the end of the step and their derivatives with respect to the left-end forces
    double M_prev = M; // (store the moment at the beginning of the step, needed for the inflection check)
    M = -Mab+Xab*u.at(2)-Zab*(x+u.at(1));
    for (int j = 1; j <= 3; j++) {
      dM.at(j) = Xab*jacobi.at(2,j)-Zab*jacobi.at(1,j);
    }
    dM.at(1) += u.at(2);
    dM.at(2) += -(x+u.at(1));
    dM.at(3) += -1.;
    kappa = computeCurvatureFromMoment(M);
    dMdkappa = computeDerMomentFromCurvature(kappa);
    dkappa = dM;
    dkappa.times(1./dMdkappa);
    // rotation at the end of the step
    u.at(3) = phi_mid+kappa*DX/2.;  
    // update last row of Jacobi matrix
    for (int j = 1; j <= 3; j++){
      jacobi.at(3,j) = dphi_mid.at(j)+dkappa.at(j)*DX/2.;
    }
     
    // test whether inflection occurs (the first occurence is considered)
    if (!inflection_detected && M*M_prev<=0. && M_prev!=M){
      // inflection point is detected - output parameters Lc, uc and jac_c will be computed
      inflection_detected = true;
      double aux = M_prev / (M_prev-M);
      //@todo:check this
      Lc = x - DX + aux * DX;
      // displacements, rotation and their derivatives at Lc by linear interpolation within the current step
      for (int i = 1; i <= 3; i++){
	uc.at(i) = u_prev.at(i) + aux*(u.at(i)-u_prev.at(i));
	for (int j = 1; j <= 3; j++) {
	  jac_c.at(i,j) = aux*jacobi.at(i,j) + (1.-aux)*jacobi_prev.at(i,j);
	}
      }
      // derivatives with respect to the spatial coordinate
      for (int i = 1; i <= 3; i++) {
	jac_c.at(i,4) = (u.at(i)-u_prev.at(i))/DX;
      }
      // return if the user cares only about the inflection point and not about the segment end
      if (inflection_only)
	return true;
      // change the rotation by deltaPhi (useful for analysis of inverted contact modes)
      u.at(3) += deltaPhi;
    }

    // test whether the end of the segment has been reached
    if (istep==nstep){
      // displacements, rotation and their derivatives at Lb by linear interpolation within the last step
      double aux = Lb/DX - (nstep-1);
      for (int i = 1; i <= 3; i++){
	ub.at(i) = u_prev.at(i) + aux*(u.at(i)-u_prev.at(i));
	for (int j = 1; j <= 3; j++) {
	  jac_b.at(i,j) = aux*jacobi.at(i,j) + (1.-aux)*jacobi_prev.at(i,j);
	}
      }
      // derivatives with respect to the spatial coordinate
      for (int i = 1; i <= 3; i++) {
	jac_b.at(i,4) = (u.at(i)-u_prev.at(i))/DX;
      }
    } 
  } // end of loop over spatial steps
  
  return inflection_detected;
}

void
NlBeamInternalContact :: transform_ub2ua(const FloatArray &ub, FloatArray &answer)
{
  answer.resize(3);
  answer.at(1) = (beamLength+ub.at(1))*cos(ub.at(3)) - ub.at(2)*sin(ub.at(3)) - beamLength;
  answer.at(2) = (beamLength+ub.at(1))*sin(ub.at(3)) + ub.at(2)*cos(ub.at(3));
  answer.at(3) = -ub.at(3);
}

void
NlBeamInternalContact :: transform_fab2fba(const FloatArray &ub, const FloatArray &fab, FloatArray &answer)
{
  answer.resize(3);
  answer.at(1) = fab.at(1)*cos(ub.at(3)) - fab.at(2)*sin(ub.at(3));
  answer.at(2) = fab.at(1)*sin(ub.at(3)) + fab.at(2)*cos(ub.at(3));
  answer.at(3) = fab.at(1)*ub.at(2) - fab.at(2)*(beamLength+ub.at(1)) - fab.at(3);
}


void
NlBeamInternalContact :: transform_fba2fab(const FloatArray &ub, const FloatArray &fba, FloatArray &answer)
{
  FloatArray ua;
  transform_ub2ua(ub, ua);
  transform_fab2fba(ua, fba, answer);
}

bool
NlBeamInternalContact :: checkLeftSegmentLengthAdmissibility(double Lc)
{
  return (Lc>=0 && Lc<=leftSegmentLength);
}

bool
NlBeamInternalContact :: checkRightSegmentLengthAdmissibility(double Lc)
{
  return (Lc>=0 && Lc<=rightSegmentLength);
}

/*
  For the given contact mode (which must be one of the tip modes),
  this method checks whether the given rotation jump at the contact point would be admissible.

  Input:
  deltaPhi ... rotation jump
  cmode ... contact mode (one of the tip modes)
  
  Return value:
  'true' if the rotation jump is admissible
*/
bool
NlBeamInternalContact :: checkRotationAdmissibility(double deltaPhi, ContactModeType cmode)
{
  if (cmode == CMT_AC || cmode == CMT_CB) // right segment must rotate counterclockwise wrt to the left one
    return (sin(deltaPhi)>=0.); 
  else // right segment must rotate clockwise
    return (sin(deltaPhi)<=0.); 
}

/////////////////////////////////////////////////////////////////////
// TIP CONTACT
/////////////////////////////////////////////////////////////////////


NlBeamInternalContact::LocalSolveStatus
NlBeamInternalContact :: findLeftEndForcesLocal_Tip_SoS(const FloatArray &ub_target, FloatArray &fab, double &Lac, double Lb, ContactModeType cmode, int process, double &Nca, double &Qca, double &deltaPhi)
{
  bool converged = false;
  FloatArray res(3), fba(3), uac(3),  dforces(3), ubc(3);
  FloatMatrix jacobi_ac(3,4), jacobi_bc(3,4), jacobi(3,3);
  double Mca, error;
  double phib = ub_target.at(3);
  double c = cos(phib);
  double s = sin(phib);
  // transformation matrix 
  FloatMatrix F(3,3);
  F = {{c, s, ub_target.at(2)},{-s,c, -beamLength-ub_target.at(1) },{0,0,-1}};
  // weights for conversion to dimensionless quantities
  double weight_disp = 1./beamLength;
  double weight_mom = beamLength/EI; 
  double weight_force = beamLength*weight_mom; 

  int iter = 0;
  while (iter++ < MAXIT_BEAM){// main iterative loop of the Newton-Raphson method
    
    // integrate along the left segment
    double L_dummy;
    FloatArray u_dummy(3);
    FloatMatrix jac_dummy(3,4);
    if (process==0){ // sticking - we do not care about the inflection point here (dummy variables used)
      trialProcess = PT_Stick;
      this->integrateAlongSegment(fab, 0., Lac, uac, jacobi_ac, L_dummy, u_dummy, jac_dummy, false);
    } else { // sliding
      trialProcess = PT_Slide;
      double La = 10.*leftSegmentLength; // a very safe estimate - we only want to find the inflection point
      bool inflex = integrateAlongSegment(fab, 0., La, u_dummy, jac_dummy, Lac, uac, jacobi_ac, true);
      if (!inflex) {
	// This is only the current Newton iterate; absence of an inflection
	// point here does not prove physical separation.
	return LocalSolveStatus::NonConverged;
      }
    }
    // internal forces at the contact point
    Nca = -fab.at(1) * cos(uac.at(3))+fab.at(2)*sin(uac.at(3));
    Qca = -fab.at(1) * sin(uac.at(3))-fab.at(2)*cos(uac.at(3));
    Mca = -fab.at(3) + fab.at(1) * uac.at(2) - fab.at(2) * (Lac + uac.at(1));
    // from equilibrium, we calculate the corresponding initial guess of the right-end forces
    transform_fab2fba(ub_target, fab, fba);
    // integrate along the right segment - we do not care about the inflection point here
    integrateAlongSegment(fba, 0., Lb, ubc, jacobi_bc, L_dummy, u_dummy, jac_dummy, false);
    // evaluate the residual
    // the first two components are always displacement differences
    res.at(1) = Lac + uac.at(1) - beamLength - ub_target.at(1) + (Lb+ubc.at(1))*c + ubc.at(2)*s;
    res.at(2) = uac.at(2) - ub_target.at(2) - (Lb+ubc.at(1))*s + ubc.at(2)*c;
    // the third component depends on the type of process
    if (process==0){ // sticking
      res.at(3) = Mca;
      //@todo
      FloatArray re(res);
      re.at(1) *= weight_disp;
      re.at(2) *= weight_disp;
      re.at(3) *= weight_mom;
      error = re.computeNorm();
    } else { // sliding
      res.at(3) = evalContactLoadingFunction(Nca, Qca, cmode);
      FloatArray re(res);
      re.at(1) *= weight_disp;
      re.at(2) *= weight_disp;
      re.at(3) *= weight_force;
      error = re.computeNorm();
    }
    // check the convergence criterion
    if (error<TOL_BEAM){// converged solution found
      // evaluate the rotation jump at the contact point
      deltaPhi = ub_target.at(3) + ubc.at(3) - uac.at(3);
      if (!stiffEvalMode){
	return LocalSolveStatus::Converged;
      }
      converged = true; // return postponed, since Jacobi should be updated
    }
    
    // set up the Jacobi matrix (derivatives of the residual wrt to the left-end forces)

    // prepare auxiliary variables
    FloatArray du(3), dw(3), dN(3), dQ(3);
    for (int j = 1; j <= 3; j++){
      for (int i = 1; i <= 3; i++){
	du.at(j) += jacobi_bc.at(1,i)*F.at(i,j);
	dw.at(j) += jacobi_bc.at(2,i)*F.at(i,j);
      }
    }
    if (process==1){ // sliding  
      // prepare derivatives of internal forces at the contact point
      for ( int j = 1; j <= 3; j++){
	//@todo:check
	dN.at(j) = (fab.at(1)*sin(uac.at(3))+fab.at(2)*cos(uac.at(3))) * jacobi_ac.at(3,j);
 	dQ.at(j) = (-fab.at(1)*cos(uac.at(3))+fab.at(2)*sin(uac.at(3))) * jacobi_ac.at(3,j);
      }
      dN.at(1) += -cos(uac.at(3));
      dN.at(2) +=  sin(uac.at(3));
      dQ.at(1) += -sin(uac.at(3));
      dQ.at(2) += -cos(uac.at(3));
    }
    // evaluate the part of Jacobi matrix that has a regular structure (additional terms come later)
    for (int j = 1; j <= 3; j++){
      jacobi.at(1,j) = jacobi_ac.at(1,j) + du.at(j)*c + dw.at(j)*s;
      jacobi.at(2,j) = jacobi_ac.at(2,j) - du.at(j)*s + dw.at(j)*c;
      if (process==0) // sticking
	jacobi.at(3,j) = fab.at(1)*jacobi_ac.at(2,j)-fab.at(2)*jacobi_ac.at(1,j); 
      else // sliding
	jacobi.at(3,j) = evalDerContactLoadingFunction(Nca, dN.at(j), dQ.at(j), cmode);
    }
    if (process==0){ // additional terms for sticking (correction of the third row)
      jacobi.at(3,1) += uac.at(2);
      jacobi.at(3,2) += -(Lac+uac.at(1));
      jacobi.at(3,3) += -1.;
    } else { // additional terms for sliding (correction of the first two rows)
      // this piece of code takes into account the fact that changes of end forces lead to changes in the position of inflection point
      double dMdx = fab.at(1)*jacobi_ac.at(2,4)-fab.at(2)*(1.+jacobi_ac.at(1,4)); 
      FloatArray dMc(3), dLac(3);
      dMc.at(1) = uac.at(2);
      dMc.at(2) = -(Lac+uac.at(1));
      dMc.at(3) = -1.;
      for (int j = 1; j <= 3; j++){
	dMc.at(j) += fab.at(1)*jacobi_ac.at(2,j)-fab.at(2)*jacobi_ac.at(1,j);
	dLac.at(j) = -dMc.at(j)/dMdx; // derivatives of the position of inflection point wrt the left-end forces
      }
      // The Coulomb residual is evaluated at the moving inflection point.
      // Its total derivative therefore contains the spatial variation of the
      // contact forces in addition to their fixed-position derivatives above.
      const double dNdx = -Qca * jacobi_ac.at(3,4);
      const double dQdx =  Nca * jacobi_ac.at(3,4);
      const double dfcdx = evalDerContactLoadingFunction(Nca, dNdx, dQdx, cmode);
      for (int j = 1; j <= 3; j++){
	jacobi.at(1,j) += (1.+jacobi_ac.at(1,4))*dLac.at(j);
	jacobi.at(2,j) += jacobi_ac.at(2,4)*dLac.at(j);
	jacobi.at(3,j) += dfcdx*dLac.at(j);
      }  
    }

    if (converged && stiffEvalMode){ // prepare a block of the stiffness matrix
      FloatMatrix iJ;
      iJ.beInverseOf(jacobi);
      auto Xphi = -s*fab.at(1)-c*fab.at(2);
      auto Zphi =  c*fab.at(1)-s*fab.at(2);
      FloatMatrix R(2,3);
      R.at(1,1) = 1. + (c*jacobi_bc.at(1,3)+s*jacobi_bc.at(2,3))*fab.at(2);
      R.at(1,2) =     -(c*jacobi_bc.at(1,3)+s*jacobi_bc.at(2,3))*fab.at(1);
      R.at(1,3) = (Lb+ubc.at(1))*s - ubc.at(2)*c - (c*jacobi_bc.at(1,1)+s*jacobi_bc.at(2,1))*Xphi - (c*jacobi_bc.at(1,2)+s*jacobi_bc.at(2,2))*Zphi; 
      R.at(2,1) =      (c*jacobi_bc.at(2,3)-s*jacobi_bc.at(1,3))*fab.at(2);
      R.at(2,2) = 1. - (c*jacobi_bc.at(2,3)-s*jacobi_bc.at(1,3))*fab.at(1);
      R.at(2,3) = (Lb+ubc.at(1))*c + ubc.at(2)*s - (c*jacobi_bc.at(2,1)-s*jacobi_bc.at(1,1))*Xphi - (c*jacobi_bc.at(2,2)-s*jacobi_bc.at(1,2))*Zphi;
      for (int j=1; j<=3; j++)
	for (int i=1; i<=3; i++)
	  Kblock.at(i,j) = iJ.at(i,1)*R.at(1,j) + iJ.at(i,2)*R.at(2,j);
      return LocalSolveStatus::Converged;
      
    }
    
    // solve the linearized problem
    jacobi.solveForRhs(res, dforces);
    // correct the left-end forces
    fab.subtract(dforces); 
  } // end of while loop
  return LocalSolveStatus::NonConverged;
}


double
NlBeamInternalContact :: shiftToIntervalFromMinusPiToPi(double x)
{
  double answer = x - 2.*M_PI*floor((x+M_PI)/(2.* M_PI));
  return answer;
}


ContactModeType
NlBeamInternalContact :: suggestedMode(double deltaPhi, bool seglength_ok, ContactModeType current_cmode)
{
  // if the tip in contact has slipped beyond the other tip, it is possible that the other tip gets activated
  if (!seglength_ok){
    switch (current_cmode){
    case CMT_AC:
      return CMT_CB;
    case CMT_BC:
      return CMT_CA;
    case CMT_CA:
      return CMT_BC;
    case CMT_CB:
      return CMT_AC;
    default:
      break;
    }
  }
  // now we know that the tip in contact has NOT slipped beyond the other tip
  // the rotation jump is first shifted into the interval from -PI to PI 
  deltaPhi = shiftToIntervalFromMinusPiToPi(deltaPhi);
  // if the rotation jump is admissible, the current_cmode should be the correct one,
  // otherwise we can expect a change into a smooth mode, depending on which boundary
  // of the admissible interval has been crossed
  switch (current_cmode){
  case CMT_AC:
    if (deltaPhi>=0.)
      return CMT_AC;
    else if (deltaPhi>-M_PI/2.)
      return CMT_AA;
    else
      return CMT_AB;
  case CMT_BC:
    if (deltaPhi<=0.)
      return CMT_BC;
    else if (deltaPhi<M_PI/2.)
      return CMT_BB;
    else
      return CMT_BA;
  case CMT_CA: 
    if (deltaPhi<=0.)
      return CMT_CA;
    else if (deltaPhi<M_PI/2.)
      return CMT_AA;
    else
      return CMT_BA;
  case CMT_CB: 
    if (deltaPhi>=0.)
      return CMT_CB;
    if (deltaPhi>-M_PI/2.)
      return CMT_BB;
    else
      return CMT_AB;
  default:
    break;
  }
  return current_cmode;
}


NlBeamInternalContact::LocalSolveStatus
NlBeamInternalContact :: findLeftEndForcesLocal_Tip(FloatArray &ub_target, FloatArray &fab,
                                                     ContactModeType cmode,
                                                     ContactModeType &resultMode)
{
  double Nca, Qca, deltaPhi, fc, Lac;
  bool rotation_ok;
  LocalSolveStatus status;
  resultMode = cmode;
  // store the initial guess for possible restart with sliding assumption
  FloatArray f_init(3), fba(3), ua_target(3);
  
  bool tipIsOnRightSegment;
  switch(cmode){
  case CMT_AC: case CMT_BC:
    tipIsOnRightSegment = true;
    break;
  case CMT_CA: case CMT_CB:
    tipIsOnRightSegment = false;
    transform_ub2ua(ub_target, ua_target);
    transform_fab2fba(ub_target, fab, fba);
    break;
  default:
    OOFEM_ERROR("NlBeamInternalContact: invalid tip-contact mode %d", cmode);
  }

  if (Process!=PT_Slide){// if it is known that the process is sliding, proceed directly to that case
  // sticking is assumed first
    
    trialProcess = PT_Stick;
    if (tipIsOnRightSegment){
      Lac = trialLeftActiveSegmentLength;
      f_init = fab; // fab stored for possible restart with sliding
      status = findLeftEndForcesLocal_Tip_SoS(ub_target, fab, Lac, rightSegmentLength, cmode, 0, Nca, Qca, deltaPhi);
    } else {
      Lac = trialRightActiveSegmentLength;
      // from equilibrium, we calculate the corresponding initial guess of the right-end forces
      f_init = fba; // fba stored for possible restart with sliding
      status = findLeftEndForcesLocal_Tip_SoS(ua_target, fba, Lac, leftSegmentLength, cmode, 0, Nca, Qca, deltaPhi);
      deltaPhi = -deltaPhi; // compensating for swapped segments
      transform_fba2fab(ub_target, fba, fab);
    }
    if (status != LocalSolveStatus::Converged) {
      // give sliding a chance - make sure that restart takes place with the original guess of fab or fba
      if (tipIsOnRightSegment) {
	fab = f_init;
      } else {
	fba = f_init;
      }
    } else { // a solution for sticking has been found but its admissibility still needs to be checked
      // A unilateral tip contact cannot transmit tension. Test this before
      // the Coulomb function; with the wrong normal-force sign its
      // absolute-value equation has no admissible sliding root.
      const bool compressiveNormal =
        (cmode == CMT_AC || cmode == CMT_CA) ? Qca >= 0. : Qca <= 0.;
      if (!compressiveNormal) {
        trialLeftActiveSegmentLength = leftSegmentLength;
        trialRightActiveSegmentLength = rightSegmentLength;
        resultMode = CMT_N;
        return LocalSolveStatus::Inadmissible;
      }
      fc = evalContactLoadingFunction(Nca, Qca, cmode);
      rotation_ok = checkRotationAdmissibility(deltaPhi, cmode);
      if (fc<=0. && rotation_ok){ // the solution with sticking assumed is admissible
	if (tipIsOnRightSegment){
	  trialRightActiveSegmentLength = rightSegmentLength;
	} else {
	  trialLeftActiveSegmentLength = leftSegmentLength;
	}
	resultMode = cmode;
	return LocalSolveStatus::Converged;
      }
    }
  }

  // sliding is assumed
  trialProcess = PT_Slide;
  if (tipIsOnRightSegment) {
    Lac = trialLeftActiveSegmentLength;
    status = findLeftEndForcesLocal_Tip_SoS(ub_target, fab, Lac, rightSegmentLength, cmode, 1, Nca, Qca, deltaPhi);
  }  else {
    Lac = trialRightActiveSegmentLength;
    status = findLeftEndForcesLocal_Tip_SoS(ua_target, fba, Lac, leftSegmentLength, cmode, 1, Nca, Qca, deltaPhi);
    deltaPhi = -deltaPhi; // compensating for swapped segments
  }
  if (status == LocalSolveStatus::NonConverged) {
    return status;
  }
  if (status == LocalSolveStatus::Inadmissible) {
    trialLeftActiveSegmentLength = leftSegmentLength;
    trialRightActiveSegmentLength = rightSegmentLength;
    resultMode = CMT_N;
    return status;
  }
  // a solution for sliding has been found but its admissibility still needs to be checked
  const bool compressiveNormal =
    (cmode == CMT_AC || cmode == CMT_CA) ? Qca >= 0. : Qca <= 0.;
  if (!compressiveNormal) {
    trialLeftActiveSegmentLength = leftSegmentLength;
    trialRightActiveSegmentLength = rightSegmentLength;
    resultMode = CMT_N;
    return LocalSolveStatus::Inadmissible;
  }
  rotation_ok = checkRotationAdmissibility(deltaPhi, cmode);
  bool seglength_ok;
  if (tipIsOnRightSegment) {
    seglength_ok = checkLeftSegmentLengthAdmissibility(Lac);
  } else {
    seglength_ok = checkRightSegmentLengthAdmissibility(Lac);
  }
  
  if (seglength_ok){ 
    if (tipIsOnRightSegment){
      trialLeftActiveSegmentLength = Lac;
      trialRightActiveSegmentLength = rightSegmentLength;
    } else {
      transform_fba2fab(ub_target, fba, fab); 
      trialLeftActiveSegmentLength = leftSegmentLength;
      trialRightActiveSegmentLength = Lac;
    }
    if (rotation_ok) {
      resultMode = cmode;
      return LocalSolveStatus::Converged; // the sliding solution is admissible
    }
  } else {
    trialLeftActiveSegmentLength = leftSegmentLength;
    trialRightActiveSegmentLength = rightSegmentLength;
  }
  resultMode = suggestedMode(deltaPhi, seglength_ok, cmode);
  return LocalSolveStatus::Inadmissible;
}

/////////////////////////////////////////////////////////////////////
// SMOOTH CONTACT
/////////////////////////////////////////////////////////////////////


NlBeamInternalContact::LocalSolveStatus
NlBeamInternalContact :: findLeftEndForcesLocal_Smooth_Rolling(FloatArray &ub_target, FloatArray &fab, double deltaPhi, double &Lac, double &Lbc, double &Nc, double &Qc)
{

  FloatArray res(3), dforces(3), ub_loc(3), uc(3);
  FloatMatrix jac_c_dummy(3,4), jac(3,3), jac_b(3,4);
 
  double L = Lac + Lbc;
  // weight for conversion to dimensionless quantities
  double weight_disp = 1./beamLength;

  int iter = 0;
  while (iter++<MAXIT_BEAM){
    // integrate along beam of length L
     bool inflection_found = integrateAlongSegment(fab, deltaPhi, L, ub_loc, jac_b, Lac, uc, jac_c_dummy, false);
     Jacobi.beSubMatrixOf(jac_b,1,3,1,3);
     
     Nc = -fab.at(1)*cos(uc.at(3))+fab.at(2)*sin(uc.at(3));
     Qc = -fab.at(1)*sin(uc.at(3))-fab.at(2)*cos(uc.at(3));

     // evaluate the residual and check the convergence criterion
    for (int i = 1; i <= 3; i++)
      res.at(i) = ub_target.at(i)-ub_loc.at(i); 
    // take into account the fact that, due to previous sliding,
    // the current beam length may differ from the initial distance between joints
    // (displacement ub_target[0] is taken with respect to the initial length of the beam)
    res.at(1) += beamLength - L;
    double error = L2norm(weight_disp*res.at(1), weight_disp*res.at(2), res.at(3));
    if (error<TOL_BEAM){// converged
      if (!inflection_found)
	return LocalSolveStatus::Inadmissible;
      else{
	Lbc = L - Lac;
	return LocalSolveStatus::Converged;
      }
    }
    // compute the iterative correction of left-end forces
    jac = Jacobi;
    jac.solveForRhs(res,dforces);
    fab.add(dforces);
  }
  // no convergence in MAXIT_BEAM iterations
  return LocalSolveStatus::NonConverged;
}



NlBeamInternalContact::LocalSolveStatus
NlBeamInternalContact :: findLeftEndForcesLocal_Smooth_Sliding(FloatArray &ub_target, FloatArray &fab, double deltaPhi, double &Lac, double &Lbc, ContactModeType cmode)
{
  bool converged = false;
  FloatArray res(4), dforces(4), ub_loc(3), uc(3);
  FloatMatrix jac_b(3,4), jac_c(3,4), jac(4,4);
  double Nc, Qc;
  FloatArray dN(3), dQ(3);
  double L = Lac + Lbc;
  // weights for conversion to dimensionless quantities
  double weight_disp = 1./beamLength;
  double weight_force = beamLength*beamLength/EI; 

  int iter = 0;
  while (iter++<MAXIT_BEAM){
    bool inflection_found = integrateAlongSegment(fab, deltaPhi, L, ub_loc, jac_b, Lac, uc, jac_c,false); 
    if (!inflection_found) {
      // The iterate is not yet a solution, so this is a numerical failure,
      // not evidence that the converged configuration has no contact.
      return LocalSolveStatus::NonConverged;
    }
    // internal forces at the contact point
    Nc = -fab.at(1)*cos(uc.at(3))+fab.at(2)*sin(uc.at(3));
    Qc = -fab.at(1)*sin(uc.at(3))-fab.at(2)*cos(uc.at(3));
    // evaluate the residual and check the convergence criterion
    res = ub_loc;
    res.subtract(ub_target);
    res.resizeWithValues(4);
    res.at(4) = evalContactLoadingFunction(Nc, Qc, cmode);
    // take into account the fact that, due to sliding,
    // the current beam length may differ from the initial distance between joints
    // (displacement ub_target[0] is taken with respect to the initial length of the beam)
    res.at(1) += L - beamLength;
    double error = L2norm(weight_disp*sqrt(res.at(1)*res.at(1)+res.at(2)*res.at(2)), res.at(3), weight_force*res.at(4));
   
    if (error < TOL_BEAM) {
      // the actual return is a bit postponed, to make update of Jacobi possible if needed
      converged = true; 
    }

    // prepare the Jacobi matrix for the set of 4 equations
    if (!converged || stiffEvalMode){
      // evaluate the first three rows
      for(int i = 1; i <= 3; i++) {
	for (int j = 1; j <= 4; j++) {
	  Jacobi44.at(i,j) = jac_b.at(i,j);
	}
      }
      // correct one entry by a term that reflects the line res[0] += L - beamLength;
      Jacobi44.at(1,4) += 1.;
      // derivatives of internal forces at the contact point
      for (int j = 1; j <= 3; j++){
	dN.at(j) = (fab.at(1)*sin(uc.at(3))+fab.at(2)*cos(uc.at(3))) * jac_c.at(3,j);
	dQ.at(j) = (-fab.at(1)*cos(uc.at(3))+fab.at(2)*sin(uc.at(3))) * jac_c.at(3,j);
      }
      dN.at(1) += -cos(uc.at(3));
      dN.at(2) += sin(uc.at(3));
      dQ.at(1) += -sin(uc.at(3));
      dQ.at(2) += -cos(uc.at(3));
      // evaluate the fourth row of Jacobi matrix 
      for (int j = 1; j <= 3; j++) {
	Jacobi44.at(4,j) = evalDerContactLoadingFunction(Nc, dN.at(j), dQ.at(j), cmode);
      }
      //Jacobi44.at(4,4) = 0.; // the total length does not affect the inflection point
    }
    // this is the postponed return after convergence
    if (converged){
      Lbc = L - Lac;
      return LocalSolveStatus::Converged;
    }
   // compute the iterative correction of left-end forces and of active beam length
    jac = Jacobi44;
    jac.solveForRhs(res,dforces);
    for(int i = 1; i <= 3; i++) {
      fab.at(i) -= dforces.at(i);
    }
    
    L -= dforces.at(4);
  }
  // no convergence in MAXIT_BEAM iterations
  return LocalSolveStatus::NonConverged;
}


NlBeamInternalContact::LocalSolveStatus
NlBeamInternalContact :: findLeftEndForcesLocal_Smooth(FloatArray &ub_target, FloatArray &fab,
                                                        ContactModeType cmode,
                                                        ContactModeType &resultMode)
{
  double Nc, Qc;
  // store the initial guess for possible restart with sliding assumption
  FloatArray fab_init(3);
  double deltaPhi, Lac, Lbc;
  LocalSolveStatus status;
  LocalSolveStatus rollingStatus = LocalSolveStatus::NonConverged;
  fab_init = fab;
  resultMode = cmode;

   deltaPhi = 0.;
   if (cmode==CMT_AB) {
     deltaPhi = M_PI;
   } else if (cmode==CMT_BA) {
     deltaPhi = -M_PI;
   }
   
  if (Process!=PT_Slide){// if it is known that the process is sliding, proceed directly to that case
    // first assume that no additional sliding occurs during this step (i.e., fix Lac+Lbc and iterate on fab)
    trialProcess = PT_Roll;  
    Lac = trialLeftActiveSegmentLength;
    Lbc = trialRightActiveSegmentLength;
    rollingStatus = findLeftEndForcesLocal_Smooth_Rolling(ub_target, fab, deltaPhi, Lac, Lbc, Nc, Qc);
    if (rollingStatus != LocalSolveStatus::Converged) {
      // give sliding a chance - make sure that the solution that assumes sliding will later start from the original guess of left-end forces 
      fab = fab_init;
    } else { // rolling solution has been found, now let us see whether it is admissible
      if (Lac > leftSegmentLength) {
        trialLeftActiveSegmentLength = leftSegmentLength;
        trialRightActiveSegmentLength = Lbc;
        resultMode = cmode == CMT_AA ? CMT_CA : CMT_CB;
        return LocalSolveStatus::Inadmissible;
      }
      if (Lbc > rightSegmentLength) {
        trialLeftActiveSegmentLength = Lac;
        trialRightActiveSegmentLength = rightSegmentLength;
        resultMode = cmode == CMT_AA ? CMT_AC : CMT_BC;
        return LocalSolveStatus::Inadmissible;
      }
      // A unilateral contact cannot transmit tension. Check this before the
      // Coulomb function: for a tensile normal force, its absolute-value
      // equation has no physically admissible sliding root and Newton can
      // oscillate around zero tangential force.
      const bool compressiveNormal =
        (cmode == CMT_AA || cmode == CMT_AB) ? Qc >= 0. : Qc <= 0.;
      if (!compressiveNormal) {
        resultMode = CMT_N;
        return LocalSolveStatus::Inadmissible;
      }
      double fc = evalContactLoadingFunction(Nc, Qc, cmode);
      if (fc<=0.){ // admissible solution
	trialLeftActiveSegmentLength = Lac;
	trialRightActiveSegmentLength = Lbc;
	resultMode = cmode;
	return LocalSolveStatus::Converged;
      }
    }
  }
  // now assume that sliding occurs during this step (i.e., iterate on fab as well as on Lac+Lbc)
  trialProcess = PT_Slide;  
  Lac = trialLeftActiveSegmentLength;
  Lbc = trialRightActiveSegmentLength;
  status = findLeftEndForcesLocal_Smooth_Sliding(ub_target, fab, deltaPhi, Lac, Lbc, cmode);
  if (status == LocalSolveStatus::NonConverged) {
    if (rollingStatus == LocalSolveStatus::Inadmissible) {
      resultMode = CMT_N;
      return LocalSolveStatus::Inadmissible;
    }
    return status;
  }
  if (status == LocalSolveStatus::Converged) { // the active segment lengths change due to sliding
    if ((Lac+Lbc)>leftSegmentLength+rightSegmentLength) {
      resultMode = CMT_N; // most likely jumps to no contact
      return LocalSolveStatus::Inadmissible;
    }
    else if (Lac > leftSegmentLength){
      trialLeftActiveSegmentLength = leftSegmentLength;
      trialRightActiveSegmentLength = Lbc;
      if (cmode == CMT_AA)
	resultMode = CMT_CA; // most likely changes into tip contact
      else
	resultMode = CMT_CB;
      return LocalSolveStatus::Inadmissible;
    }
    else if (Lbc > rightSegmentLength) {
      trialLeftActiveSegmentLength = Lac;
      trialRightActiveSegmentLength = rightSegmentLength;
      if (cmode == CMT_AA)
	resultMode = CMT_AC; // most likely changes into tip contact
      else
	resultMode = CMT_BC;
      return LocalSolveStatus::Inadmissible;
    }
    // admissible solution
    trialLeftActiveSegmentLength = Lac;
    trialRightActiveSegmentLength = Lbc;
    resultMode = cmode;
    return LocalSolveStatus::Converged;
  }
  resultMode = CMT_N; // no inflection point: most likely smoothly loses contact
  return LocalSolveStatus::Inadmissible;
}


ContactModeType
NlBeamInternalContact :: predictContactMode(FloatArray &ub)
{
  // assuming that both segments are straight, find their intersection
  
  double phib = ub.at(3);
  const double sinPhib = sin(phib);
  // Near the reference configuration, roundoff can leave an infinitesimal
  // relative rotation. Treat the segments as parallel before evaluating the
  // otherwise ill-conditioned line-intersection formula below.
  if (fabs(sinPhib) <= TOL_BEAM)
    return CMT_N;
  double Lbc = -ub.at(2)/sinPhib;
  if (Lbc<0. || Lbc>rightSegmentLength) {
    return CMT_N;
  }
  double Lac = beamLength+ub.at(1)-Lbc*cos(phib);
  if (Lac < 0. || Lac > leftSegmentLength) {
    return CMT_N;
  }
  if (rightSegmentLength - Lbc < leftSegmentLength - Lac) { // right tip probably touches the left segment first 
    if (ub.at(2)>0.) {
      return CMT_BC;
    } else {
      return CMT_AC;
    }
  }
  // left tip probably touches the right segment first
  if (ub.at(2) > 0.) {
    return CMT_CA;
  }
  return CMT_CB;
}


double
NlBeamInternalContact :: findTipContactTime(FloatArray &u_prev, FloatArray &du, double segLength)
{
  double t;
  // simple linearized estimate, used as the initial guess for the iterative scheme
  double db = du.at(2) + segLength*(sin(u_prev.at(3)+du.at(3))-sin(u_prev.at(3))); 
  if (db!=0.){
    t = (-u_prev.at(2)-segLength*sin(u_prev.at(3)))/db;
  } else {
    t = 0.5;
  }
  // iterative solution of nonlinear equation
  for (int iter = 1; iter <= MAXIT_CONTACT_TIME; iter++){ 
    double f = u_prev.at(2) + t*du.at(2) + segLength*sin(u_prev.at(3)+t*du.at(3));
    if (fabs(f) < TOL_CONTACT_TIME * segLength) {
	return t;// converged solution
    }
    double fp = du.at(2) + segLength*cos(u_prev.at(3)+t*du.at(3))*du.at(3);
    if (fp!=0.) {
	t -= f/fp;
    }
  }
  // converged solution not found
  return -1;
}


ContactModeType
NlBeamInternalContact :: predictContactMode(FloatArray ub, FloatArray ub_prev)
{
  // This routine detects a new contact event along the increment from the
  // committed configuration. With no relative increment there is no event.
  // In particular, do not let findTipContactTime() turn a stationary,
  // near-zero gap into an artificial impact at its fallback time t = 0.5.
  if (ub.distance(ub_prev) == 0.) {
    return CMT_N;
  }
  // assuming that both segments are straight,
  //        check at which time 'tb' and location 'xb' the right tip hits the left segment
  // (the 'time' is a dimensionless parameter running from 0 at the beginning to 1 at the end of the step)
  double xa, xb = 0;
  FloatArray du(3);
  double tol = 0;
  du = ub;
  du.subtract(ub_prev);
  double tb = findTipContactTime(ub_prev, du, rightSegmentLength);
  if (tb>=0.-tol && tb<=1.+tol) {
    // find location of the hit event
    xb = beamLength + ub_prev.at(1) + tb*du.at(1) - rightSegmentLength*cos(ub_prev.at(3)+tb*du.at(3));
  }
  // now check at which time 'ta' and location 'xa' the left tip hits the right segment
  FloatArray ua(3), ua_prev(3);
  transform_ub2ua(ub, ua);
  transform_ub2ua(ub_prev, ua_prev);
  du = ua;
  du.subtract(ua_prev);
  double ta = findTipContactTime(ua_prev, du, leftSegmentLength);
  if (ta >= 0. && ta <= 1.) {
    // find location of the hit event
    xa = beamLength + ua_prev.at(1) + ta*du.at(1) - leftSegmentLength*cos(ua_prev.at(3)+ta*du.at(3));
  }
  // determine which hit event occurs at an admissible time and location
  bool left_hits = ( ta >= 0. && ta <= 1. && xa >= 0. && xa <= rightSegmentLength);
  bool right_hits = (tb >= 0. && tb <= 1. && xb >= 0. && xb <= leftSegmentLength);

  if (!left_hits && !right_hits) {
    return CMT_N;
  }
  if (left_hits && right_hits) {
    left_hits = ta < tb;
  }
  if (left_hits) {
    trialLeftActiveSegmentLength = leftSegmentLength;
    trialRightActiveSegmentLength = xa;
    if (ua.at(2)-ua_prev.at(2)+leftSegmentLength*(sin(ua.at(3))-sin(ua_prev.at(3))) > 0.) {
      return CMT_CA;
    } else {
      return CMT_CB;
    }
  }
  trialLeftActiveSegmentLength = xb;
  trialRightActiveSegmentLength = rightSegmentLength;
  if (ub.at(2)-ub_prev.at(2)+rightSegmentLength*(sin(ub.at(3))-sin(ub_prev.at(3))) > 0.) {
    return CMT_AC;
  } else {
    return CMT_BC;
  }
}


bool
NlBeamInternalContact :: findLeftEndForcesLocal(FloatArray &ub, FloatArray &ub_prev,
                                                 FloatArray &fab,
                                                 ContactModeType initialMode)
{
  ContactModeType assumed_cmode = initialMode;
  trialLeftActiveSegmentLength = leftActiveSegmentLength;
  trialRightActiveSegmentLength = rightActiveSegmentLength;
  for (int attempt = 1; attempt <= MAXIT_CONTACT_MODE; attempt++) {
    LocalSolveStatus status = LocalSolveStatus::Inadmissible;
    switch (assumed_cmode){
      // smooth modes
    case CMT_AA:
    case CMT_BB:
    case CMT_AB:
    case CMT_BA:
      status = findLeftEndForcesLocal_Smooth(ub, fab, assumed_cmode, trialContactMode);
      break;
      // tip modes
    case CMT_AC:
    case CMT_BC:
    case CMT_CA:
    case CMT_CB:
      status = findLeftEndForcesLocal_Tip(ub, fab, assumed_cmode, trialContactMode);
      break;
      // no-contact mode
    case CMT_N:
      fab.zero();
      if (contactMode != CMT_N) {
	trialContactMode = predictContactMode(ub);
      } else {
	trialContactMode = predictContactMode(ub, ub_prev);
      }
      status = trialContactMode == CMT_N ? LocalSolveStatus::Converged
                                        : LocalSolveStatus::Inadmissible;
      break;
      // this should never happen
    default:
      OOFEM_ERROR("Mode %d has not been implemented yet in findLeftEndForcesLocal\n",assumed_cmode);
    }
    if (status == LocalSolveStatus::NonConverged) {
      return false;
    }
    if (status == LocalSolveStatus::Converged && trialContactMode == assumed_cmode) {
      Process = trialProcess;
      return true;
    }
    if (trialContactMode == assumed_cmode) {
      return false;
    }
    // solution not accepted, let us try another one with modified assumption which mode occurs
    assumed_cmode = trialContactMode;
  }
  // At first activation from a committed open state, several admissible-set
  // candidates can meet at the nonsmooth zero-force point. If a tangent is
  // requested before the residual has selected a branch, retain the committed
  // open-set tangent instead of rejecting the whole global trial state. The
  // subsequent residual evaluation still performs the full branch search.
  if (stiffEvalMode && initialMode == CMT_N) {
    fab.zero();
    trialContactMode = CMT_N;
    trialProcess = PT_Unknown;
    trialLeftActiveSegmentLength = leftActiveSegmentLength;
    trialRightActiveSegmentLength = rightActiveSegmentLength;
    return true;
  }
  // no admissible solution found after MAXIT_CONTACT_MODE attempts
  return false;
}


void
NlBeamInternalContact :: construct_T(FloatMatrix &T, double phia)
{
  T.resize(3,3);
  T.at(1,1) = T.at(2,2) = cosAlpha*cos(phia) + sinAlpha*sin(phia);
  T.at(1,2) =  sinAlpha *cos(phia) - cosAlpha * sin(phia);
  T.at(2,1) = - T.at(1,2);
  T.at(3,3) = 1.;
}


void
NlBeamInternalContact :: construct_T(FloatMatrix &T, double phia, double alpha)
{
  T.resize(3,3);
  T.at(1,1) = T.at(2,2) = cos(alpha-phia);
  T.at(1,2) = sin(alpha-phia);
  T.at(2,1) = -T.at(1,2);
  T.at(3,3) = 1.;
}

/*
void
NlBeamInternalContact :: construct_T(FloatMatrix &T, double phia, double rot)
{
  T.resize(3,3);
  T.at(1,1) = T.at(2,2) = (cosAlpha*cos(rot) + sinAlpha*sin(rot))*cos(phia) + (sinAlpha*cos(rot)-cosAlpha*sin(rot))*sin(phia);
  T.at(1,2) =  (sinAlpha*cos(rot)-cosAlpha*sin(rot)) *cos(phia) - (cosAlpha*cos(rot)+sinAlpha*sin(rot)) * sin(phia);
  T.at(2,1) = - T.at(1,2);
  T.at(3,3) = 1.;
}
*/

void
NlBeamInternalContact :: construct_Tprime(FloatMatrix &T, double phia)
{
    T.resize(3,3);
  T.at(1,1) = T.at(2,2) = sinAlpha *cos(phia) - cosAlpha*sin(phia);
  T.at(1,2) = -cosAlpha*cos(phia) - sinAlpha*sin(phia);
  T.at(2,1) = -T.at(1,2);
  /*
  T.resize(3,3);
  T.at(1,1) = T.at(2,2) = sin(alpha-phia);
  T.at(1,2) = -cos(alpha-phia);
  T.at(2,1) = -T.at(1,2);
  */

}

void
NlBeamInternalContact :: construct_l(FloatArray &l, double phia)
{
  l.resize(3);
  /*l.at(1) = beamLength * (cos(phia)*cosAlpha - sin(phia) * sinAlpha - cosAlpha);
  l.at(2) = beamLength * (sinAlpha * cos(phia) + cosAlpha * sin(phia) - sinAlpha);
  */
  l.at(1) = beamLength * (cos(phia)-1.);
  l.at(2) = beamLength * sin(phia);
  
}

void
NlBeamInternalContact :: construct_l_loc(FloatArray &l, double phia)
{
  l.resize(3);
  l.at(1) = beamLength * (cos(phia)-1.);
  l.at(2) = beamLength * sin(phia);
}


void
NlBeamInternalContact :: construct_l(FloatArray &l, double phia, double L)
{
  l.resize(3);
  /*  l.at(1) = L * (cos(phia)*cosAlpha - sin(phia) * sinAlpha - cosAlpha);
  l.at(2) = L * (sinAlpha * cos(phia) + cosAlpha * sin(phia) - sinAlpha);
  */
  l.at(1) = L * (cos(phia) - 1.);
  l.at(2) = L *  sin(phia);
  
}



void
NlBeamInternalContact :: construct_l_IC(FloatArray &l, double phia, double L)
{
  l.at(1) = beamLength * cos(phia) - L;
  l.at(2) = beamLength * sin(phia);
}

void
NlBeamInternalContact :: construct_lprime(FloatArray &l, double phia)
{
  l.resize(3);
  /*  l.at(1) = -beamLength * (sin(phia)*cosAlpha + cos(phia) * sinAlpha);
  l.at(2) =  beamLength * (cosAlpha * cos(phia) -  sinAlpha * sin(phia));
  */
  
  l.at(1) = -beamLength * sin(phia);
  l.at(2) = beamLength * cos(phia);
  
}


bool
NlBeamInternalContact :: findLeftEndForces(const FloatArray &u, const FloatArray &u_prev,
                                            FloatArray &fab, ContactModeType initialMode,
                                            ProcessType initialProcess)
{
  FloatArray ub_loc(3), ub_prev_loc(3), fab_loc(3);
  FloatMatrix T(3,3), T_prev(3,3);
  double phia = u.at(3);
  double phia_prev = u_prev.at(3);
  Process = initialProcess;
  trialProcess = PT_Unknown;
  // compute displacements of the right end with respect to the auxiliary coordinate system
  FloatArray l_ub_loc, l_ub_prev_loc;
  construct_l(l_ub_loc, phia);
  construct_l(l_ub_prev_loc, phia_prev);
  construct_T(T, phia); 
  construct_T(T_prev, phia_prev);


  FloatArray u_a,u_b, u_ba, uprev_ba, uprev_a, uprev_b;
  u_a.beSubArrayOf(u, {1,2,3});
  u_b.beSubArrayOf(u, {4,5,6});
  u_ba.beDifferenceOf(u_b, u_a);
  
  uprev_a.beSubArrayOf(u_prev, {1,2,3});
  uprev_b.beSubArrayOf(u_prev, {4,5,6});
  uprev_ba.beDifferenceOf(uprev_b, uprev_a);

  ub_loc.beProductOf(T, u_ba);
  ub_loc.add(l_ub_loc);
  ub_prev_loc.beProductOf(T_prev,uprev_ba);
  ub_prev_loc.add(l_ub_prev_loc);
  fab_loc.beProductOf(T,fab);
  // compute the corresponding left-end forces and moment in the auxiliary coordinate system
  if(!findLeftEndForcesLocal(ub_loc, ub_prev_loc, fab_loc, initialMode)) {
    return false;
  }
  //@todo: comment the following 3 lines???
  /* contactMode = trialContactMode;
  leftActiveSegmentLength = trialLeftActiveSegmentLength;
  rightActiveSegmentLength = trialRightActiveSegmentLength;
  */
  fab.beTProductOf(T,fab_loc);
  //  printf("%g %g %g\n",fab.at(1),fab.at(2),fab.at(3));
  return true;
}


void
NlBeamInternalContact :: givePreviousStepDisplacements(FloatArray &answer,
                                                        const FloatArray &current,
                                                        TimeStep *tStep)
{
  if (previousStepDisplacements.giveSize() == current.giveSize()) {
    answer = previousStepDisplacements;
    return;
  }

  // The first step has no element-owned history yet. Preserve the established
  // OOFEM convention, including nonzero data prescribed at the initial time.
  FloatArray increment;
  this->computeVectorOf({D_u, D_w, R_v}, VM_Incremental, tStep, increment);
  answer = current - increment;
  previousStepDisplacements = answer;
}

bool
NlBeamInternalContact :: hasIterationBranchCache(const FloatArray &u,
                                                    const FloatArray &uPrev,
                                                    TimeStep *tStep) const
{
  return tStep && branchCacheValid &&
         branchCacheStep == tStep->giveNumber() &&
         branchCacheState == tStep->giveSolutionStateCounter() &&
         u.giveSize() == branchCacheDisplacements.giveSize() &&
         uPrev.giveSize() == branchCachePreviousDisplacements.giveSize() &&
         u.distance(branchCacheDisplacements) == 0. &&
         uPrev.distance(branchCachePreviousDisplacements) == 0.;
}


bool
NlBeamInternalContact :: hasExactTangentCache(const FloatArray &u,
                                                const FloatArray &uPrev,
                                                TimeStep *tStep) const
{
  return tStep && tangentCacheValid &&
         tangentCacheStep == tStep->giveNumber() &&
         tangentCacheState == tStep->giveSolutionStateCounter() &&
         u.giveSize() == tangentCacheDisplacements.giveSize() &&
         uPrev.giveSize() == tangentCachePreviousDisplacements.giveSize() &&
         u.distance(tangentCacheDisplacements) == 0. &&
         uPrev.distance(tangentCachePreviousDisplacements) == 0.;
}


void
NlBeamInternalContact :: publishTangentCache(const FloatArray &u,
                                                const FloatArray &uPrev,
                                                const FloatArray &fab,
                                                TimeStep *tStep)
{
  tangentCacheValid = true;
  tangentCacheStep = tStep ? tStep->giveNumber() : -1;
  tangentCacheState = tStep ? tStep->giveSolutionStateCounter() : -1;
  tangentCacheMode = trialContactMode;
  tangentCacheProcess = trialProcess;
  tangentCacheLeftActiveSegmentLength = trialLeftActiveSegmentLength;
  tangentCacheRightActiveSegmentLength = trialRightActiveSegmentLength;
  tangentCacheForces = fab;
  tangentCacheDisplacements = u;
  tangentCachePreviousDisplacements = uPrev;
  tangentCacheJacobi = Jacobi;
  tangentCacheJacobi44 = Jacobi44;
  tangentCacheKblock = Kblock;
}


void
NlBeamInternalContact :: restoreTangentCache(FloatArray &fab)
{
  trialContactMode = tangentCacheMode;
  Process = trialProcess = tangentCacheProcess;
  trialLeftActiveSegmentLength = tangentCacheLeftActiveSegmentLength;
  trialRightActiveSegmentLength = tangentCacheRightActiveSegmentLength;
  fab = tangentCacheForces;
  Jacobi = tangentCacheJacobi;
  Jacobi44 = tangentCacheJacobi44;
  Kblock = tangentCacheKblock;
}


bool
NlBeamInternalContact :: isProcessLockActive(TimeStep *tStep) const
{
  if (lockedProcess == PT_Unknown || processHistoryStep != tStep->giveNumber()) {
    return false;
  }
  const int substep = tStep->giveSubStepNumber();
  return substep >= processLockFirstSubstep && substep <= processLockLastSubstep;
}


void
NlBeamInternalContact :: resetIterationBranchState()
{
  branchCacheStep = -1;
  branchCacheState = -1;
  branchCacheValid = false;
  branchCacheMode = CMT_N;
  branchCacheProcess = PT_Unknown;
  branchCacheDisplacements.clear();
  branchCachePreviousDisplacements.clear();
  processHistoryStep = -1;
  processHistorySubstep = -1;
  processHistoryCount = 0;
  for (int i = 0; i < 4; ++i) {
    processModeHistory[i] = CMT_N;
    processHistory[i] = PT_Unknown;
  }
  lockedContactMode = CMT_N;
  lockedProcess = PT_Unknown;
  processLockFirstSubstep = -1;
  processLockLastSubstep = -1;
}


void
NlBeamInternalContact :: recordIterationBranch(TimeStep *tStep)
{
  if (processLockSteps == 0) {
    return;
  }

  const int step = tStep->giveNumber();
  const int substep = tStep->giveSubStepNumber();
  if (processHistoryStep != step) {
    processHistoryStep = step;
    processHistorySubstep = -1;
    processHistoryCount = 0;
    lockedContactMode = CMT_N;
    lockedProcess = PT_Unknown;
    processLockFirstSubstep = -1;
    processLockLastSubstep = -1;
  }
  if (processHistorySubstep == substep) {
    return;
  }
  processHistorySubstep = substep;

  if (lockedProcess != PT_Unknown && substep > processLockLastSubstep) {
    lockedContactMode = CMT_N;
    lockedProcess = PT_Unknown;
    processLockFirstSubstep = -1;
    processLockLastSubstep = -1;
  }

  if (processHistoryCount < 4) {
    processModeHistory[processHistoryCount] = trialContactMode;
    processHistory[processHistoryCount] = trialProcess;
    ++processHistoryCount;
  } else {
    for (int i = 0; i < 3; ++i) {
      processModeHistory[i] = processModeHistory[i + 1];
      processHistory[i] = processHistory[i + 1];
    }
    processModeHistory[3] = trialContactMode;
    processHistory[3] = trialProcess;
  }

  if (processHistoryCount < 4 || lockedProcess != PT_Unknown) {
    return;
  }

  const bool sameMode = processModeHistory[0] == processModeHistory[1] &&
                        processModeHistory[0] == processModeHistory[2] &&
                        processModeHistory[0] == processModeHistory[3];
  bool hasSlide = false;
  bool hasNonSlidingContact = false;
  int processChanges = 0;
  for (int i = 0; i < 4; ++i) {
    hasSlide = hasSlide || processHistory[i] == PT_Slide;
    hasNonSlidingContact = hasNonSlidingContact ||
                           processHistory[i] == PT_Stick ||
                           processHistory[i] == PT_Roll;
    if (i > 0 && processHistory[i] != processHistory[i - 1]) {
      ++processChanges;
    }
  }
  // Chatter is not necessarily the exact A-B-A-B sequence. For example,
  // slide-slide-stick-slide already contains two branch switches and can
  // repeat indefinitely. Prefer sliding for a short window; the locked
  // branch is still released immediately if it becomes inadmissible.
  if (sameMode && hasSlide && hasNonSlidingContact && processChanges >= 2 &&
      processHistory[3] == PT_Slide) {
    lockedContactMode = processModeHistory[3];
    lockedProcess = PT_Slide;
    processLockFirstSubstep = substep + 1;
    processLockLastSubstep = substep + processLockSteps;
    OOFEM_LOG_INFO("NlBeamInternalContact element %d: process chatter in step %d; "
                   "retaining sliding for Newton substeps %d-%d\n",
                   this->giveNumber(), step, processLockFirstSubstep,
                   processLockLastSubstep);
  }
}


bool
NlBeamInternalContact :: solveForKinematics(const FloatArray &u, const FloatArray &uPrev,
                                             FloatArray &fab, TimeStep *tStep,
                                             bool evaluateTangent, bool useIterationState)
{
  this->ensureContactDataInitialized();

  const ContactModeType savedTrialContactMode = trialContactMode;
  const ProcessType savedProcess = Process;
  const ProcessType savedTrialProcess = trialProcess;
  const double savedTrialLeftLength = trialLeftActiveSegmentLength;
  const double savedTrialRightLength = trialRightActiveSegmentLength;
  const bool savedStiffEvalMode = stiffEvalMode;
  const bool savedTangentCacheValid = tangentCacheValid;
  const FloatArray savedFab(fab);
  const FloatMatrix savedJacobi(Jacobi);
  const FloatMatrix savedJacobi44(Jacobi44);
  const FloatMatrix savedKblock(Kblock);

  ContactModeType initialMode = contactMode;
  ProcessType initialProcess = PT_Unknown;
  bool usesHint = false;
  const bool usesLock = useIterationState && tStep && this->isProcessLockActive(tStep);
  if (usesLock) {
    initialMode = lockedContactMode;
    initialProcess = lockedProcess;
    usesHint = true;
  } else if (useIterationState && this->hasIterationBranchCache(u, uPrev, tStep)) {
    // Residual and tangent evaluations at the same trial state must use the
    // same active-set hint. This reproduces the ordering assumed by the
    // original standalone implementation and removes call-order dependence.
    initialMode = branchCacheMode;
    initialProcess = branchCacheProcess;
    usesHint = true;
  }

  auto restoreAttemptState = [&]() {
    trialContactMode = savedTrialContactMode;
    Process = savedProcess;
    trialProcess = savedTrialProcess;
    trialLeftActiveSegmentLength = savedTrialLeftLength;
    trialRightActiveSegmentLength = savedTrialRightLength;
    fab = savedFab;
    Jacobi = savedJacobi;
    Jacobi44 = savedJacobi44;
    Kblock = savedKblock;
    tangentCacheValid = savedTangentCacheValid;
  };

  auto attempt = [&](ContactModeType mode, ProcessType process) {
    restoreAttemptState();
    Process = process;
    trialProcess = PT_Unknown;
    stiffEvalMode = evaluateTangent;
    return findLeftEndForces(u, uPrev, fab, mode, process);
  };

  bool success = attempt(initialMode, initialProcess);
  bool triedNoContact = initialMode == CMT_N;
  if (usesLock && success &&
      (trialContactMode != lockedContactMode || trialProcess != lockedProcess)) {
    success = false;
  }
  if (!success && usesHint) {
    if (usesLock) {
      OOFEM_LOG_INFO("NlBeamInternalContact element %d: releasing process lock "
                     "after the locked branch became inadmissible\n",
                     this->giveNumber());
      // Discard the cycle which created this lock.  Otherwise the fallback
      // result can complete the opposite A-B-A-B phase and immediately
      // recreate the same inadmissible lock.
      this->resetIterationBranchState();
    }
    success = attempt(contactMode, PT_Unknown);
    triedNoContact = triedNoContact || contactMode == CMT_N;
  }
  if (!success && !triedNoContact) {
    // At the parallel reference configuration the zero-force contact state is
    // non-unique. A tiny roundoff-level rotation can make the contact branch
    // singular, so accept separation only inside a narrow dimensionless
    // neighbourhood of that configuration.
    FloatArray uA, uB, uBA, ubLocal, localOffset;
    FloatMatrix localTransform(3, 3);
    uA.beSubArrayOf(u, {1, 2, 3});
    uB.beSubArrayOf(u, {4, 5, 6});
    uBA.beDifferenceOf(uB, uA);
    construct_T(localTransform, u.at(3));
    construct_l(localOffset, u.at(3));
    ubLocal.beProductOf(localTransform, uBA);
    ubLocal.add(localOffset);

    double lengthScale = beamLength;
    if (leftSegmentLength > lengthScale) {
      lengthScale = leftSegmentLength;
    }
    if (rightSegmentLength > lengthScale) {
      lengthScale = rightSegmentLength;
    }
    const double referenceDistance =
      L2norm(ubLocal.at(1) / lengthScale, ubLocal.at(2) / lengthScale, ubLocal.at(3));
    const double referenceTolerance = TOL_BEAM > 1.e-10 ? TOL_BEAM : 1.e-10;
    if (referenceDistance <= referenceTolerance) {
      success = attempt(CMT_N, PT_Unknown);
      success = success && trialContactMode == CMT_N;
    }
  }
  stiffEvalMode = savedStiffEvalMode;
  if (success) {
    if (!useIterationState) {
      // Diagnostic/finite-difference probes need the computed forces but must
      // not alter the active set, trial lengths, or cached Jacobians used by
      // the subsequent physical residual and tangent evaluations.
      FloatArray diagnosticForces(fab);
      restoreAttemptState();
      fab = diagnosticForces;
      return true;
    }
    if (evaluateTangent) {
      this->publishTangentCache(u, uPrev, fab, tStep);
    }
    if (trialContactMode == CMT_N) {
      Process = trialProcess = PT_Unknown;
    }
    if (useIterationState && tStep) {
      branchCacheStep = tStep->giveNumber();
      branchCacheState = tStep->giveSolutionStateCounter();
      branchCacheValid = true;
      branchCacheMode = trialContactMode;
      branchCacheProcess = trialProcess;
      branchCacheDisplacements = u;
      branchCachePreviousDisplacements = uPrev;
    }
    return true;
  }

  restoreAttemptState();
  stiffEvalMode = savedStiffEvalMode;
  return false;
}


void
NlBeamInternalContact :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
  FloatArray u;
  this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, u);

  if (!useUpdatedGpRecord) {
    FloatArray uPrev;
    this->givePreviousStepDisplacements(uPrev, u, tStep);

    FloatArray candidateForces = internalForces;
    // Residual and line-search probes may visit inadmissible intermediate
    // states. Select the branch here without publishing a tangent; the matrix
    // assembly builds and caches the tangent transactionally afterwards.
    if (!solveForKinematics(u, uPrev, candidateForces, tStep, false)) {
      OOFEM_ERROR("NlBeamInternalContact element %d: local contact solve failed in step %d",
                  this->giveNumber(), tStep->giveNumber());
    }
    tempInternalForces = candidateForces;
    this->recordIterationBranch(tStep);
  }

  answer.resize(6);
  answer.at(1) = tempInternalForces.at(1);
  answer.at(2) = tempInternalForces.at(2);
  answer.at(3) = tempInternalForces.at(3);
  answer.at(4) = -answer.at(1);
  answer.at(5) = -answer.at(2);
  const double c1 = beamLength * sinAlpha + u.at(5) - u.at(2);
  const double c2 = -beamLength * cosAlpha - u.at(4) + u.at(1);
  answer.at(6) = c1 * answer.at(1) + c2 * answer.at(2) - answer.at(3);
}


void
NlBeamInternalContact :: giveInternalForcesVector_fromU(FloatArray &answer, TimeStep *tStep,
                                                        FloatArray &u, FloatArray &u_prev)
{
  FloatArray candidateForces = internalForces;
  if (!solveForKinematics(u, u_prev, candidateForces, tStep, false, false)) {
    OOFEM_ERROR("NlBeamInternalContact element %d: local contact solve failed in step %d",
                this->giveNumber(), tStep->giveNumber());
  }

  answer.resize(6);
  answer.at(1) = candidateForces.at(1);
  answer.at(2) = candidateForces.at(2);
  answer.at(3) = candidateForces.at(3);
  answer.at(4) = -answer.at(1);
  answer.at(5) = -answer.at(2);
  const double c1 = beamLength * sinAlpha + u.at(5) - u.at(2);
  const double c2 = -beamLength * cosAlpha - u.at(4) + u.at(1);
  answer.at(6) = c1 * answer.at(1) + c2 * answer.at(2) - answer.at(3);
}


void
NlBeamInternalContact :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
  answer.resize(6,6);
  answer.zero();
  FloatArray u(6), u_prev(6);
  this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, u);
  this->givePreviousStepDisplacements(u_prev, u, tStep);
  FloatArray lprime(3), fab_loc(3), ub_loc(3), Tu(3), fab(3);
  FloatMatrix T(3,3), Tprime(3,3), Ginv(3,3), Ginv44(4,4), TtGinv(3,3);
  const bool currentBranchCached = this->hasIterationBranchCache(u, u_prev, tStep);
  if (this->hasExactTangentCache(u, u_prev, tStep)) {
    this->restoreTangentCache(fab);
    predictorTangentAvailable = false;
  } else if (!currentBranchCached && predictorTangentAvailable &&
             tangentCacheValid && tStep->giveSubStepNumber() == 0) {
    // NRSolver intentionally requests one old tangent before the first
    // residual of a new step. Keep the complete cached branch together:
    // forces, active lengths, process, and local Jacobians.
    this->restoreTangentCache(fab);
    predictorTangentAvailable = false;
  } else {
    fab = currentBranchCached ? tempInternalForces : internalForces;
    if (!solveForKinematics(u, u_prev, fab, tStep, true)) {
      OOFEM_ERROR("NlBeamInternalContact element %d: local tangent solve failed in step %d",
                  this->giveNumber(), tStep->giveNumber());
    }
    tempInternalForces = fab;
    predictorTangentAvailable = false;
  }
  tempInternalForces = fab;
  // The step predictor intentionally reuses the previous local contact tangent
  // and maps it through the current trial geometry. This is a quasi-tangent,
  // not the fully assembled tangent from the preceding configuration.
  double phia = u.at(3);
  construct_T(T, phia); 
  construct_Tprime(Tprime, phia); 
  construct_lprime(lprime, phia);
  //  
  fab_loc.beProductOf(T,fab);

  // Relative right-end kinematics in the local frame attached to the left end.
  // This is needed below to map the tangent back after swapping the segments
  // in the reversed tip-contact modes CA and CB.
  FloatArray u_a, u_b, u_ba, l_loc(3);
  u_a.beSubArrayOf(u, {1,2,3});
  u_b.beSubArrayOf(u, {4,5,6});
  u_ba.beDifferenceOf(u_b, u_a);
  construct_l(l_loc, phia);
  ub_loc.beProductOf(T, u_ba);
  ub_loc.add(l_loc);

  // get Jacobi matrix in local coordinates and invert it
  switch (trialContactMode) {
    // smooth modes
  case CMT_AA:
  case CMT_BB:
  case CMT_AB:
  case CMT_BA:
    if (trialProcess == PT_Slide) {
      Ginv44.beInverseOf(Jacobi44);
      Ginv.beSubMatrixOf(Ginv44,1,3,1,3);
    } else if (trialProcess == PT_Roll) {
      Ginv.beInverseOf(Jacobi);
    } else {
      OOFEM_ERROR("NlBeamInternalContact: missing process state for smooth-contact tangent");
    }
    break;
    // tip modes
  case CMT_AC:
  case CMT_BC:
    Ginv.beSubMatrixOf(Kblock,1,3,1,3);
    break;
  case CMT_CA:
  case CMT_CB: {
    // The tip solve swaps the two segments. Kblock is therefore
    // d(fba)/d(ua), where ua is the relative displacement after the swap.
    // Transform it back to d(fab)/d(ub) by the complete chain rule.
    FloatArray ua_loc(3), fba_loc(3);
    transform_ub2ua(ub_loc, ua_loc);
    transform_fab2fba(ub_loc, fab_loc, fba_loc);

    const double c = cos(ub_loc.at(3));
    const double s = sin(ub_loc.at(3));
    FloatMatrix swappedDisplacementJacobian(3,3);
    FloatMatrix inverseForceTransform(3,3);
    FloatMatrix geometricForceDerivative(3,3);
    FloatMatrix transformedKblock(3,3);

    // d ua / d ub
    swappedDisplacementJacobian.at(1,1) = c;
    swappedDisplacementJacobian.at(1,2) = -s;
    swappedDisplacementJacobian.at(1,3) = -ua_loc.at(2);
    swappedDisplacementJacobian.at(2,1) = s;
    swappedDisplacementJacobian.at(2,2) = c;
    swappedDisplacementJacobian.at(2,3) = beamLength + ua_loc.at(1);
    swappedDisplacementJacobian.at(3,3) = -1.;

    // d fab / d fba at fixed ua
    inverseForceTransform.at(1,1) = c;
    inverseForceTransform.at(1,2) = s;
    inverseForceTransform.at(2,1) = -s;
    inverseForceTransform.at(2,2) = c;
    inverseForceTransform.at(3,1) = ua_loc.at(2);
    inverseForceTransform.at(3,2) = -(beamLength + ua_loc.at(1));
    inverseForceTransform.at(3,3) = -1.;

    // Explicit derivative of the force transformation with respect to ua,
    // with fba held fixed.
    geometricForceDerivative.at(1,3) = -fab_loc.at(2);
    geometricForceDerivative.at(2,3) = fab_loc.at(1);
    geometricForceDerivative.at(3,1) = -fba_loc.at(2);
    geometricForceDerivative.at(3,2) = fba_loc.at(1);

    transformedKblock.beProductOf(inverseForceTransform, Kblock);
    transformedKblock.add(geometricForceDerivative);
    Ginv.beProductOf(transformedKblock, swappedDisplacementJacobian);
    break;
  }
    // no-contact mode
  case CMT_N:
    answer.zero();
    return;
    // this should never happen
  default:
    OOFEM_ERROR("NlBeamInternalContact: invalid contact mode %d in tangent", trialContactMode);
  }

  TtGinv.beTProductOf(T, Ginv);

  
  // compute product Ttransposed*Ginverse*T and store it in the upper stiffness block
  FloatMatrix aux;
  aux.beProductOf(TtGinv, T);
  answer.setSubMatrix(aux, 1,1);
  answer.times(-1);
  answer.setSubMatrix(aux, 1,4);

  // compute additional stiffness terms in the third column (associated with phia)
  Tu.beProductOf(Tprime, u_ba);
  FloatArray col3_1,col3_2;
  Tu.add(lprime);
  col3_1.beProductOf(TtGinv,Tu);
  col3_2.beTProductOf(Tprime, fab_loc);
  answer.addSubVectorCol(col3_1, 1,3);
  answer.addSubVectorCol(col3_2, 1,3);
			      
  // construct the sixth row (complete formula)
  double c1 = beamLength*sinAlpha + u.at(5) - u.at(2);
  double c2 = -beamLength*cosAlpha - u.at(4) + u.at(1);
  answer.at(6,1) =  fab.at(2);
  answer.at(6,2) = -fab.at(1);
  answer.at(6,4) = -fab.at(2);
  answer.at(6,5) =  fab.at(1);
  for(int j = 1; j <= 6; j++) {
    answer.at(4,j) = - answer.at(1,j);
    answer.at(5,j) = - answer.at(2,j);
    answer.at(6,j) += c1 * answer.at(1,j) + c2 * answer.at(2,j) - answer.at(3,j);
  }

  //compute the stiffness matrix numerically
  /*  if(tStep->giveNumber() > 132){
    double du = 1.e-6; // small perturbation for numerical evaluation of stiffness
    FloatArray f;
    f = this->internalForces;
    FloatMatrix f_pert(6,6), Knum(6,6);
    
    for (int i=1; i<=6; i++){
      u.at(i) += du;
      this->giveInternalForcesVector_fromU(f, tStep, u, u_prev);
      u.at(i) -= du;
      for (int j=1; j<=6; j++){
	f_pert.at(j,i) = f.at(j);
      }
    }

    f = this->internalForces;
    this->giveInternalForcesVector_fromU(f, tStep, u, u_prev);
    
    
    double stifferr = 0;
    double stiffdiagsum = 0;
    for (int i=1; i<=6; i++){
      for (int j=1; j<=6; j++){
	Knum.at(i,j) = (f_pert.at(i,j)-f.at(i))/du;
	stifferr += fabs(Knum.at(i,j)-answer.at(i,j));
	if (i==j) {
	  stiffdiagsum += answer.at(i,j);
	}
      }
    }
    double ei = stifferr/stiffdiagsum;
    if(tStep->giveNumber() > 230){
      answer = Knum;
    }
    }*/

}

void
NlBeamInternalContact :: initForNewStep()
{
    NlBeam_SM :: initForNewStep();
    trialContactMode = contactMode;
    trialLeftActiveSegmentLength = leftActiveSegmentLength;
    trialRightActiveSegmentLength = rightActiveSegmentLength;
    tempInternalForces = internalForces;
    // Retain the last accepted process together with its Jacobian. It is the
    // predictor tangent if the next step asks for a matrix before a residual.
    trialProcess = Process;
    predictorTangentAvailable = tangentCacheValid;
    stiffEvalMode = false;
    this->resetIterationBranchState();
}


void
NlBeamInternalContact :: saveContext(DataStream &stream, ContextMode mode)
{
    NlBeam_SM :: saveContext(stream, mode);
    if (mode & CM_Definition) {
        if (!stream.write(leftSegmentLength) ||
            !stream.write(rightSegmentLength) ||
            !stream.write(DX) ||
            !stream.write(friction) ||
            !stream.write(TOL_BEAM) ||
            !stream.write(MAXIT_BEAM) ||
            !stream.write(processLockSteps) ||
            !stream.write(static_cast<int>(contactMode)) ||
            !stream.write(leftActiveSegmentLength) ||
            !stream.write(rightActiveSegmentLength)) {
            THROW_CIOERR(CIO_IOERR);
        }
        contextIOResultType result = matNum.storeYourself(stream);
        if (result != CIO_OK) {
            THROW_CIOERR(result);
        }
    }
    if (mode & CM_State) {
        if (!stream.write(static_cast<int>(contactMode)) ||
            !stream.write(leftActiveSegmentLength) ||
            !stream.write(rightActiveSegmentLength)) {
            THROW_CIOERR(CIO_IOERR);
        }
        contextIOResultType result = internalForces.storeYourself(stream);
        if (result != CIO_OK) {
            THROW_CIOERR(result);
        }
        result = previousStepDisplacements.storeYourself(stream);
        if (result != CIO_OK) {
            THROW_CIOERR(result);
        }
    }
}


void
NlBeamInternalContact :: restoreContext(DataStream &stream, ContextMode mode)
{
    NlBeam_SM :: restoreContext(stream, mode);
    if (mode & CM_Definition) {
        int restoredDefinitionMode;
        if (!stream.read(leftSegmentLength) ||
            !stream.read(rightSegmentLength) ||
            !stream.read(DX) ||
            !stream.read(friction) ||
            !stream.read(TOL_BEAM) ||
            !stream.read(MAXIT_BEAM) ||
            !stream.read(processLockSteps) ||
            !stream.read(restoredDefinitionMode) ||
            !stream.read(leftActiveSegmentLength) ||
            !stream.read(rightActiveSegmentLength)) {
            THROW_CIOERR(CIO_IOERR);
        }
        contextIOResultType result = matNum.restoreYourself(stream);
        if (result != CIO_OK) {
            THROW_CIOERR(result);
        }
        if (leftSegmentLength <= 0. || rightSegmentLength <= 0. ||
            DX <= 0. || friction < 0. || TOL_BEAM <= 0. || MAXIT_BEAM < 1 ||
            processLockSteps < 0 ||
            matNum.giveSize() != 2) {
            THROW_CIOERRM(CIO_BADVERSION, "invalid internal-contact definition in restart");
        }
        if (restoredDefinitionMode < CMT_N || restoredDefinitionMode > CMT_CB ||
            leftActiveSegmentLength < 0. || leftActiveSegmentLength > leftSegmentLength ||
            rightActiveSegmentLength < 0. || rightActiveSegmentLength > rightSegmentLength) {
            THROW_CIOERRM(CIO_BADVERSION, "invalid internal-contact initial state in restart");
        }
        contactMode = static_cast<ContactModeType>(restoredDefinitionMode);
        trialContactMode = contactMode;
        trialLeftActiveSegmentLength = leftActiveSegmentLength;
        trialRightActiveSegmentLength = rightActiveSegmentLength;
        contactDataInitialized = false;
        this->resetIterationBranchState();
    }
    if (mode & CM_State) {
        int restoredMode;
        if (!stream.read(restoredMode) ||
            !stream.read(leftActiveSegmentLength) ||
            !stream.read(rightActiveSegmentLength)) {
            THROW_CIOERR(CIO_IOERR);
        }
        contextIOResultType result = internalForces.restoreYourself(stream);
        if (result != CIO_OK) {
            THROW_CIOERR(result);
        }
        result = previousStepDisplacements.restoreYourself(stream);
        if (result != CIO_OK) {
            THROW_CIOERR(result);
        }
        if (internalForces.giveSize() != 3 || previousStepDisplacements.giveSize() != 6) {
            THROW_CIOERRM(CIO_BADVERSION, "invalid internal-contact history in restart");
        }
        if (leftActiveSegmentLength < 0. || leftActiveSegmentLength > leftSegmentLength ||
            rightActiveSegmentLength < 0. || rightActiveSegmentLength > rightSegmentLength) {
            THROW_CIOERRM(CIO_BADVERSION, "invalid active segment lengths in restart");
        }
        if (restoredMode < CMT_N || restoredMode > CMT_CB) {
            THROW_CIOERRM(CIO_BADVERSION, "invalid internal-contact mode in restart");
        }
        contactMode = static_cast<ContactModeType>(restoredMode);
        trialContactMode = contactMode;
        trialLeftActiveSegmentLength = leftActiveSegmentLength;
        trialRightActiveSegmentLength = rightActiveSegmentLength;
        tempInternalForces = internalForces;
        Process = trialProcess = PT_Unknown;
        stiffEvalMode = false;
        tangentCacheValid = false;
        predictorTangentAvailable = false;
        this->resetIterationBranchState();
    }
}


void
NlBeamInternalContact :: updateYourself(TimeStep *tStep)
{
    FloatArray currentDisplacements;
    this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, currentDisplacements);

    NlBeam_SM :: updateYourself(tStep);
    contactMode = trialContactMode;
    leftActiveSegmentLength = trialLeftActiveSegmentLength;
    rightActiveSegmentLength = trialRightActiveSegmentLength;
    this->internalForces = this->tempInternalForces;
    previousStepDisplacements = currentDisplacements;
    Process = trialProcess;
    predictorTangentAvailable = tangentCacheValid;
    stiffEvalMode = false;
}


Interface *NlBeamInternalContact :: giveInterface(InterfaceType it)
{
    switch ( it ) {
    case VTKXMLExportModuleElementInterfaceType:
      return static_cast< VTKXMLExportModuleElementInterface * >( this );
    default:
      return StructuralElement :: giveInterface(it);
    }
}

void
NlBeamInternalContact :: printOutputAt(FILE *file, TimeStep *tStep)
{
  static const char *modeNames[] = {"N", "AA", "AB", "AC", "BA", "BB", "BC", "CA", "CB"};
  static const char *processNames[] = {"unknown", "roll", "stick", "slide"};
  NlBeam_SM :: printOutputAt(file, tStep);
  fprintf(file,
          " Contact Mode %d (%s) Process %s LeftActive %.15e RightActive %.15e\n",
          static_cast<int>(contactMode), modeNames[contactMode], processNames[Process],
          leftActiveSegmentLength, rightActiveSegmentLength);
}


void
NlBeamInternalContact :: giveCompositeExportData(std :: vector< ExportRegion > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep )
{
  this->ensureContactDataInitialized();

  //two segments
  vtkPieces.resize(2);
  
  FloatArray ul, xLeft, xRight, Lleft, Lright;
  FloatMatrix uMatrixLeft, uMatrixRight, T;
  
  FloatArray uab, ua(3),ub(3);
  this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, uab);
  ua.beSubArrayOf(uab,{1,2,3});
  ub.beSubArrayOf(uab,{4,5,6});
  //
  /*To be removed*/
  double phia = uab.at(3);
  FloatArray ub_loc, ub_ua(ub), loc, fab_loc, u0(3);
  ub_ua.subtract(ua);  
  // compute displacements of the right end with respect to the auxiliary coordinate system
  construct_l_loc(loc, phia);
  //construct_T(T, phia,0.0); 
  construct_T(T, phia); 
  ub_loc.beProductOf(T, ub_ua);
  ub_loc.add(loc);
  fab_loc.beProductOf(T, this->internalForces);
  FloatMatrix Tn = {{1,0},{0,1}};  
  this->computeSegmentDisplacements(uMatrixLeft,Lleft, fab_loc, leftActiveSegmentLength, leftSegmentLength, u0, Tn, uab );
  // 
  ////
  FloatMatrix Tb,Ttb;
  FloatArray ur(3);
  ur.at(1) = ub_loc.at(1) + beamLength;
  ur.at(2) = ub_loc.at(2);
  ur.at(3) = ub_loc.at(3);
  double alphab = M_PI + ub_loc.at(3);
  //construct_T(Tb, alphab, 0.0);
  construct_T(Tb, alphab,0.0);
  Ttb.beTranspositionOf(Tb);
  // transformation of forces and moment
  FloatArray f;
  transform_fab2fba(ub_loc, fab_loc, f);
  this->computeSegmentDisplacements(uMatrixRight,Lright, f, rightActiveSegmentLength, rightSegmentLength, ur, Ttb, uab, true );

 
  Node *nodeA = this->giveNode(1);
  //
  const int numCellNodes  = 2; // linear line
  //
  int numCellsLeft = Lleft.giveSize()-1;
  int nNodesLeft = numCellsLeft * numCellNodes;
  vtkPieces.at(0).setNumberOfCells(numCellsLeft);
  vtkPieces.at(0).setNumberOfNodes(nNodesLeft);
  //
  int numCellsRight = Lright.giveSize()-1;
  int nNodesRight = numCellsRight * numCellNodes;
  vtkPieces.at(1).setNumberOfCells(numCellsRight);
  vtkPieces.at(1).setNumberOfNodes(nNodesRight);
  //
  /// Piece 0
  int val    = 1;
  int offset = 0;
  IntArray nodes(numCellNodes);
  FloatArray ug(3);
  int nodeNum = 1;
  FloatArray nodeCoords(3);
  IntArray connectivity(2);
  for ( int iElement = 1; iElement <= numCellsLeft; iElement++ ) {
    for (int nN = 1; nN <= numCellNodes; nN++) {
      int lN = nN % 2;
      int iNode;
      if(lN == 0) {
	iNode = nN/2 + 1;
      } else {
	iNode = (nN + 1) / 2;
      }

      nodeCoords.at(1) = nodeA->giveCoordinate(1);// + Lleft.at(iElement+iNode-1) * cosAlpha;
      nodeCoords.at(2) = 0;
      nodeCoords.at(3) = nodeA->giveCoordinate(3);// + Lleft.at(iElement+iNode-1) * sinAlpha;
      //
      vtkPieces.at(0).setNodeCoords(nodeNum, nodeCoords);
      nodeNum++;
      connectivity.at(iNode) = val++;
    }
    vtkPieces.at(0).setConnectivity(iElement, connectivity);
    offset += 2;
    vtkPieces.at(0).setOffset(iElement, offset);
    vtkPieces.at(0).setCellType(iElement, 3);
  }
  /// Piece 1
  nodeNum = 1;
  val    = 1;
  offset = 0;
  for ( int iElement = 1; iElement <= numCellsRight; iElement++ ) {
    for (int nN = 1; nN <= numCellNodes; nN++) {
      int lN = nN % 2;
      int iNode;
      if(lN == 0) {
	iNode = nN/2 + 1;
      } else {
	iNode = (nN + 1) / 2;
      }

      nodeCoords.at(1) = nodeA->giveCoordinate(1);// - Lright.at(iElement+iNode-1) * cosAlpha;
      nodeCoords.at(2) = 0;
      nodeCoords.at(3) = nodeA->giveCoordinate(3);// - Lright.at(iElement+iNode-1) * sinAlpha;
      //
      vtkPieces.at(1).setNodeCoords(nodeNum, nodeCoords);
      nodeNum++;
      connectivity.at(iNode) = val++;
    }
    vtkPieces.at(1).setConnectivity(iElement, connectivity);
    offset += 2;
    vtkPieces.at(1).setOffset(iElement, offset);
    vtkPieces.at(1).setCellType(iElement, 3);
  }


  

    int n = primaryVarsToExport.giveSize();
    vtkPieces [ 0 ].setNumberOfPrimaryVarsToExport(primaryVarsToExport, nNodesLeft);
    vtkPieces [ 1 ].setNumberOfPrimaryVarsToExport(primaryVarsToExport, nNodesRight);
    for ( int i = 1; i <= n; i++ ) {
        UnknownType utype = ( UnknownType ) primaryVarsToExport.at(i);
        if ( utype == DisplacementVector ) {
	  for ( int nN = 1; nN <= nNodesLeft; nN++ ) {
	    int lN = nN % 2;
	    int iNode;
	    if(lN == 0) {
	      iNode = nN/2 + 1;
	    } else {
	      iNode = (nN + 1) / 2;
	    }
	    /*
	    double L = Lleft.at(iNode);
	    FloatArray l;
	    this->construct_l(l, uab.at(3), L);
	    this->construct_T(T, uab.at(3));
	    FloatArray u_l, u_g;
	    u_l = {uMatrixLeft.at(iNode, 1),uMatrixLeft.at(iNode, 2), 0};
	    u_l.subtract(l);
	    u_g.beTProductOf(T, u_l);
	    ug.at(1) = u_g.at(1) + uab.at(1);
	    ug.at(3) = u_g.at(2) + uab.at(2);
	    ug.at(2) = 0;
	    vtkPieces.at(0).setPrimaryVarInNode(i, nN, ug);
	    */
	    vtkPieces.at(0).setPrimaryVarInNode(utype, nN, {uMatrixLeft.at(iNode, 1),0,uMatrixLeft.at(iNode, 2)});
	  }
	  for ( int nN = 1; nN <= nNodesRight; nN++ ) {
	    int lN = nN % 2;
	    int iNode;
	    if(lN == 0) {
	      iNode = nN/2 + 1;
	    } else {
	      iNode = (nN + 1) / 2;
	    }

	    /*	    double L = beamLength - Lright.at(iNode);
	    FloatArray l;
	    this->construct_l(l, uab.at(3), L);
	    this->construct_T(T, uab.at(3));
	    FloatArray u_l, u_g;
	    u_l = {uMatrixRight.at(iNode, 1),uMatrixRight.at(iNode, 2), 0};
	    u_l.subtract(l);
	    u_g.beTProductOf(T, u_l);
	    ug.at(1) = u_g.at(1) + uab.at(1);
	    ug.at(3) = u_g.at(2) + uab.at(2);
	    ug.at(2) = 0;
	    vtkPieces.at(1).setPrimaryVarInNode(i, nN, ug);
	    */
	    vtkPieces.at(1).setPrimaryVarInNode(utype, nN, {uMatrixRight.at(iNode, 1),0,uMatrixRight.at(iNode, 2)});

	  }
        }
    }


    int nC = cellVarsToExport.giveSize();
    vtkPieces [ 0 ].setNumberOfCellVarsToExport(cellVarsToExport, numCellsLeft);
    vtkPieces [ 1 ].setNumberOfCellVarsToExport(cellVarsToExport, numCellsRight);
    for ( int i = 1; i <= nC; i++ ) {
      InternalStateType type = ( InternalStateType ) cellVarsToExport.at(i);
        if ( type == IST_MaterialNumber ) {
	  for (int iC = 1; iC <= numCellsLeft; iC++) {
	    vtkPieces [ 0 ].setCellVar(type, iC, {static_cast<double>(matNum.at(1))});
	  }
	  for (int iC = 1; iC <= numCellsRight; iC++) {
	    vtkPieces [ 1 ].setCellVar(type, iC, {static_cast<double>(matNum.at(2))});
	  }
	}	
    }
    
}

void
NlBeamInternalContact ::computeSegmentDisplacements(FloatMatrix &uMatrix, FloatArray &Larray, const FloatArray &fab, double Lb, double segmentLength, const FloatArray &u0, const FloatMatrix &T, const FloatArray &uab, bool rightSegment)
{
  FloatArray ub;
  bool inflection_detected = false;
  double aux;
  
  // initialization at the left end
  double Xab = fab.at(1);
  double Zab = fab.at(2);
  double Mab = fab.at(3);
  double M = -Mab;
  double kappa = computeCurvatureFromMoment(M);
  double x = 0.;
  FloatArray u(3), u_prev(3);
  //
  FloatArray l, ul(3);
  FloatMatrix Tg;
  //
  // basic loop over spatial steps
  int nstep = ceil(Lb/DX); // the spatial step size is fixed and the segment length does not need to be its integer multiple
  Larray.resize(nstep+1);
  uMatrix.resize(nstep+1, 2);

  if(rightSegment) {
    this->construct_l(l, uab.at(3), beamLength);
    l.times(0);
    this->construct_T(Tg, uab.at(3));
  } else {
    this->construct_l(l, uab.at(3), 0);
    this->construct_T(Tg, uab.at(3));
  }
  // global displacement for the first point
  FloatArray ug;
  ul.at(1) = u0.at(1);
  ul.at(2) = u0.at(2);
  //
  ul.subtract(l);
  ug.beTProductOf(Tg, ul);
  ug.at(1) += uab.at(1);
  ug.at(2) += uab.at(2);
  uMatrix.at(1,1) = ug.at(1);
  uMatrix.at(1,2) = ug.at(2);
  //
  Larray.at(1) = x;
  //  
  for (int istep=1; istep<=nstep; istep++){
    x += DX;
    u_prev = u;
    // rotation at midstep and its derivatives with respect to the left-end forces
    double phi_mid = u_prev.at(3) + kappa * DX / 2.;
    // normal force at midstep 
    double N_mid = -Xab*cos(phi_mid)+Zab*sin(phi_mid);
    // horizontal displacement at the end of the step
    u.at(1) = u_prev.at(1) + DX * ((1.+N_mid/EA)*cos(phi_mid)-1.);
   // vertical displacement at the end of the step
    u.at(2) = u_prev.at(2) - DX * ( 1. + N_mid / EA ) * sin(phi_mid);
    // bending moment and curvature at the end of the step 
    double M_prev = M; // (store the moment at the beginning of the step, needed for the inflection check)
    M = -Mab+Xab*u.at(2)-Zab*(x+u.at(1));
    kappa = computeCurvatureFromMoment(M);
    // rotation at the end of the step
    u.at(3) = phi_mid+kappa*DX/2.;
    // test whether inflection occurs (the first occurence is considered)
    if (!inflection_detected && M*M_prev<=0. && M_prev!=M){
      // inflection point is detected - output parameter Lc, uc will be computed      inflection_detected = true;
    }

    // test whether the end of the segment has been reached
    if (istep==nstep){
      // displacements and rotation at Lb by linear interpolation within the last step
      aux = Lb/DX - (nstep-1);
      ub.resize(3);
      for (int i=1; i<=3; i++) {
	ub.at(i) = u_prev.at(i) + aux*(u.at(i)-u_prev.at(i));
      }
      /*
      uMatrix.at(istep+1, 1) = u0.at(1) + T.at(1,1)*(Lb+ub.at(1)) + T.at(1,2)*ub.at(2);
      uMatrix.at(istep+1, 2) = u0.at(2) + T.at(2,1)*(Lb+ub.at(1)) + T.at(2,2)*ub.at(2);
      */
      ul.at(1) = u0.at(1) + T.at(1,1)*(Lb+ub.at(1)) + T.at(1,2)*ub.at(2);
      ul.at(2) = u0.at(2) + T.at(2,1)*(Lb+ub.at(1)) + T.at(2,2)*ub.at(2);

      Larray.at(istep+1) = Lb;
      // plot also the straight segment behind the contact point, if it exists
      double Lstraight = segmentLength - Lb;

      if (Lstraight>0.){
	/*	
	uMatrix.at(istep+1, 1) = u0.at(1) + T.at(1,1)*(Lb+ub.at(1)+Lstraight*cos(ub.at(3))) + T.at(1,2)*(ub.at(2)-Lstraight*sin(ub.at(3))) ;
	uMatrix.at(istep+1, 2) = u0.at(2) + T.at(2,1)*(Lb+ub.at(1)+Lstraight*cos(ub.at(3))) + T.at(2,2)*(ub.at(2)-Lstraight*sin(ub.at(3)));
	*/
	ul.at(1) = u0.at(1) + T.at(1,1)*(Lb+ub.at(1)+Lstraight*cos(ub.at(3))) + T.at(1,2)*(ub.at(2)-Lstraight*sin(ub.at(3))) ;
	ul.at(2) = u0.at(2) + T.at(2,1)*(Lb+ub.at(1)+Lstraight*cos(ub.at(3))) + T.at(2,2)*(ub.at(2)-Lstraight*sin(ub.at(3)));
	
	//
	Larray.at(istep-1) = Lb;
      }
    } else {
      /*
      uMatrix.at(istep+1,1) = u0.at(1) + T.at(1,1)*(x+u.at(1)) + T.at(1,2)*u.at(2);
      uMatrix.at(istep+1,2) = u0.at(2) + T.at(2,1)*(x+u.at(1)) + T.at(2,2)*u.at(2);
      */
      FloatArray ll;
      if(rightSegment)
	this->construct_l(ll, u0.at(3), x);
      else {
	this->construct_l(ll, u0.at(3), x);
      }
      ll.times(0.);
      ul.at(1) = u0.at(1) + T.at(1,1)*(x+u.at(1) - ll.at(1)) + T.at(1,2)*(u.at(2) - ll.at(2));
      ul.at(2) = u0.at(2) + T.at(2,1)*(x+u.at(1) - ll.at(1)) + T.at(2,2)*(u.at(2) - ll.at(2));

      
      Larray.at(istep+1) = x;
    }
    if(rightSegment) {
      this->construct_l(l, uab.at(3), beamLength-x);
    } else {
      this->construct_l(l, uab.at(3), x);
    }
    FloatArray ug;
    //ul.subtract(l);
    ug.beTProductOf(Tg, ul);
    ug.at(1) += uab.at(1);
    ug.at(2) += uab.at(2);
    uMatrix.at(istep+1,1) = ug.at(1);
    uMatrix.at(istep+1,2) = ug.at(2);
    /*    if(rightSegment) {
      uMatrix.at(istep+1,1) += beamLength * Tg(1,1);
      uMatrix.at(istep+1,2) += beamLength * Tg(1,2);
      }*/
    //
    
    
    /* ul.subtract(l);
    ug.beTProductOf(Tg, ul);
    uMatrix.at(istep+1,1) = ug.at(1) + uab.at(1);
    uMatrix.at(istep+1,2) = ug.at(2) + uab.at(2);
    */

  } // end of loop over spatial steps
  
}


double
NlBeamInternalContact :: computeLength()
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
NlBeamInternalContact :: givePitch()
// Returns the pitch of the receiver.
{
    double xA, xB, yA, yB;
    Node *nodeA, *nodeB;

    if ( pitch == 10. ) {             // 10. : dummy initialization value
        nodeA  = this->giveNode(1);
        nodeB  = this->giveNode(2);
        xA     = nodeA->giveCoordinate(1);
        xB     = nodeB->giveCoordinate(1);
        yA     = nodeA->giveCoordinate(3);
        yB     = nodeB->giveCoordinate(3);
        pitch  = atan2(yB - yA, xB - xA);
    }

    return pitch;
}


  
}
