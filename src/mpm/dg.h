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

#ifndef dg_h
#define dg_h

#include "engngm.h"
#include "sparselinsystemnm.h"
#include "sparsenonlinsystemnm.h"
#include "sparsemtrx.h"
#include "primaryfield.h"
#include "function.h"
#include "dofdistributedprimaryfield.h"
#include "mpm.h"

///@name Input fields for DGproblem
//@{
#define _IFT_DGProblem_Name "dgproblem"
#define _IFT_DGProblem_initt "initt"
#define _IFT_DGProblem_deltat "deltat"
#define _IFT_DGProblem_deltatfunction "deltatfunction"
#define _IFT_DGProblem_prescribedtimes "prescribedtimes"
#define _IFT_DGProblem_alpha "alpha"
#define _IFT_DGProblem_keepTangent "keeptangent" ///< Fixes the tangent to be reused on each step.
#define _IFT_DGProblem_exportFields "exportfields" ///< Fields to export for staggered problems.
#define _IFT_DGProblem_problemType "ptype" 
//@}

namespace oofem {

/**
 * Callback class for assembling mid point effective tangents. 
 * @todo Need to parametrize individual contributing terms, ther locations and multilication factors.
 */
class ScalarAdvectionLhsAssembler : public MatrixAssembler
{
protected:
    double alpha;
    double deltaT;
    Variable::VariableQuantity q;

public:
    ScalarAdvectionLhsAssembler(double alpha, double deltaT, Variable::VariableQuantity q) ;
    void matrixFromElement(FloatMatrix &mat, Element &element, TimeStep *tStep) const override;
};

/**
 * Callback class for assembling mid point effective tangents. 
 * @todo Need to parametrize individual contributing terms, ther locations and multilication factors.
 */
class ScalarAdvectionRhsAssembler : public MatrixAssembler
{
protected:
    double alpha;
    double deltaT;
    Variable::VariableQuantity q;

public:
    ScalarAdvectionRhsAssembler(double alpha, double deltaT, Variable::VariableQuantity q) ;
    void matrixFromElement(FloatMatrix &mat, Element &element, TimeStep *tStep) const override;
};

/**
 * Callback class for assembling residuals
 */
class ScalarAdvectionResidualAssembler : public MatrixAssembler
{
    protected:
    double alpha;
    double deltaT;
    Variable::VariableQuantity q;

public:
    ScalarAdvectionResidualAssembler(double alpha, double deltaT, Variable::VariableQuantity q) : MatrixAssembler(), alpha(alpha), deltaT(deltaT), q(q) {}
    void matrixFromElement(FloatMatrix &mat, Element &element, TimeStep *tStep) const override;
};


/**
 * This class represents generic Discontinuous Galerkin solver. The problem can be growing/decreasing, signalized by flag "changingproblemsize"
 * in the problem description. The solution is stored in UnknownsField, which can obtain/ project solution from/to DOFs (nodes). If the problem
 * keeps the same equation numbers, solution is taken from UnknownsField without any projection, which is more efficient. See the matlibmanual
 * for solution strategy of balance equations and the solution algorithm.
 *
 * @todo Documentation errors (there is no "UnknownsField" used here).
 */
class DGProblem : public EngngModel
{
protected:

    Variable::VariableQuantity unknownQuantity = Variable::VariableQuantity::VolumeFraction;
    
    LinSystSolverType solverType = ST_Direct;
    SparseMtrxType sparseMtrxType = SMT_Skyline;
    std :: unique_ptr< DofDistributedPrimaryField > field;

    std :: unique_ptr< SparseMtrx > lhsMatrix;
    std :: unique_ptr< SparseMtrx > rhsMatrix;

    FloatArray solution;
    FloatArray internalForces;
    FloatArray eNorm;

    /// Numerical method used to solve the problem
    std :: unique_ptr< SparseLinearSystemNM > nMethod;

    /// Initial time from which the computation runs. Default is zero.
    double initT = 0.;
    /// Length of time step.
    double deltaT = 0.;
    double alpha = 0.;
    /// Associated time function for time step increment.
    int dtFunction = 0;
    /// Specified times where the problem is solved
    FloatArray prescribedTimes;
    bool keepTangent = false, hasTangent = false;
    IntArray exportFields;
    /// identifies what problem to solve (UP, UPV, etc) 
    std::string problemType; 


    // class representing boundary entity (edge, surface, etc.)
    class DGBoundaryEntity
    {
    public:
        IntArray elements; // list of elements sharing entity (size should be 1 or 2)
        IntArray elementBoundaryIDs; // list of boundary IDs of elements sharing entity
        DGBoundaryEntity() : elements(), elementBoundaryIDs() { }
        virtual ~DGBoundaryEntity() { }
        void addElement(int elem, int boundaryID) { elements.followedBy(elem); elementBoundaryIDs.followedBy(boundaryID); }
    };

    // list of boundary entities
    std :: vector< std::unique_ptr<DGBoundaryEntity>> boundaryEntities;


public:
    DGProblem(int i, EngngModel * _master);

    void constructBoundaryEntities ();

    void solveYourselfAt(TimeStep *tStep) override;
    void updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d) override;
    bool newDofHandling() override { return true; }
    void updateSolution(FloatArray &solutionVector, TimeStep *tStep, Domain *d) override;
    void updateInternalRHS(FloatArray &answer, TimeStep *tStep, Domain *d, FloatArray *eNorm) override;
    void updateMatrix(SparseMtrx &mat, TimeStep *tStep, Domain *d) override;
    double giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof) override;
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    virtual void applyIC();

    int requiresUnknownsDictionaryUpdate() override;
    int giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep) override;
    void updateDomainLinks() override;

    Function *giveDtFunction();
    double giveDeltaT(int n);
    double giveDiscreteTime(int iStep);

    TimeStep *giveNextStep() override;
    TimeStep *giveSolutionStepWhenIcApply(bool force = false) override;
    NumericalMethod *giveNumericalMethod(MetaStep *mStep) override;

    void initializeFrom(InputRecord &ir) override;

    bool requiresEquationRenumbering(TimeStep *tStep) override;
    int forceEquationNumbering() override;

    void printOutputAt(FILE *file, TimeStep *tStep) override;
    
    void updateYourself(TimeStep *tStep) override;
    
    int checkConsistency() override;
    FieldPtr giveField (FieldType key, TimeStep *tStep) override;
    // identification
    const char *giveInputRecordName() const { return _IFT_DGProblem_Name; }
    const char *giveClassName() const override { return "DGProblem"; }
    fMode giveFormulation() override { return TL; }

  
  double giveFinalTime() //override
  {
    if(prescribedTimes.giveSize()) {
      return prescribedTimes.at(prescribedTimes.giveSize());
    } else {
      return deltaT * numberOfSteps;
    }
  }

};

} // end namespace oofem
#endif // dg_h
