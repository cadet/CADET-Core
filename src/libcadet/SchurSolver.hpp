// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2014: Eric von Lieres¹, Joel Andersson,
//                         Andreas Puettmann¹, Sebastian Schnittert¹,
//                         Samuel Leweke¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#ifndef SCHURSOLVER_HPP_
#define SCHURSOLVER_HPP_

#include <sundials/sundials_spgmr.h>
#include <idas/idas_impl.h>

#include "SimulatorPImpl.hpp"
#include "JacobianData.hpp"
#include "Timer.hpp"

namespace cadet {

class SchurSolver
{
public:

    // Constructor
    SchurSolver(const SimulatorPImpl & sim);
    // Destructor
    ~SchurSolver();

    void assembleJacCDisc(double alpha);
    void assembleJacPDisc(int pblk, double alpha);

    void factorizeJacCDisc();
    void factorizeJacPDisc(int pblk);

    //  solve_with_jacC()
    //  solves the linear system
    //    J_column * x = rhs
    //  with the given right hand side vector rhs.
    //
    //  The solution x will be returned in rhs.
    //
    //  NOTE: This function uses the discretized system
    //        Jacobian, which is assumed to have been
    //        assembled and factorized before!
    void solveWithJacCDisc(double* rhs);

    //  solve_with_jacP()
    //  solves the linear system
    //    J_particle * x = rhs
    //  with the given right hand side vector rhs
    //  for particle p.
    //
    //  The solution x will be returned in rhs.
    //
    //  NOTE: This function uses the discretized system
    //        Jacobian, which is moreover assumed to
    //        have been factorized before!
    void solveWithJacPDisc(int pblk, double* rhs);

    //  schurComplementTimesVector()
    //  multiplies the Schur complement S by the vector v
    //  and stores the result in z:
    //     z = S * v
    int schurComplementTimesVector(void* userData, N_Vector NV_v, N_Vector NV_z);
    int schurSolve(IDAMem IDA_mem, N_Vector NV_rhs, N_Vector weight,
            N_Vector NV_yCur, N_Vector NV_yDotCur, N_Vector NV_resCur);


    // Setter for parameters set by user
    inline void configure(double schurSafety = 1e-8, int maxRestarts = 0, int maxKrylov = 0, int gramSchmidtType = 1)
    {
        setSchurSafety(schurSafety);
        setMaxRestarts(maxRestarts);
        setMaxKrylov(maxKrylov);
        setGramSchmidtType(gramSchmidtType);
    }

    inline void setMaxKrylov(int maxKrylov)             { _maxKrylov = (maxKrylov > 0) ? maxKrylov : _cc.neq_bnd(); }
    inline int  getMaxKrylov()                    const { return _maxKrylov; }

    inline void setMaxRestarts(int maxRestarts)         {_maxRestarts = maxRestarts; }
    inline int  getMaxRestarts()                  const { return _maxRestarts; }

    inline void setGramSchmidtType(int gramSchmidtType) {_gramSchmidtType = gramSchmidtType; }
    inline int  getGramSchmidtType()              const { return _gramSchmidtType; }

    inline void setSchurSafety(double schurSafety)      { _schurSafety = schurSafety; }
    inline double getSchurSafety()                const { return _schurSafety; }

    inline long int getSpgmrIters()               const { return _spgmrIters; }
    inline void resetSpgmrIters()                       { _spgmrIters = 0; }

    inline const JacobianData& getJacobianData()  const { return _jac; }
    inline       JacobianData& getJacobianData()        { return _jac; }


    // Timer read functions
    inline double timerFact()           const { return _timerFact.getTime(); }
    inline double timerSol()            const { return _timerSol.getTime(); }
    inline double timerSolSpgmr()       const { return _timerSolSpgmr.getTime(); }
    inline double timerSolSpgmrAxb()    const { return _timerSolSpgmrAxb.getTime(); }

    inline double timerFactPar()        const { return _timerFactPar.getTime(); }
    inline double timerSolPar()         const { return _timerSolPar.getTime(); }
    inline double timerSolSpgmrAxbPar() const { return _timerSolSpgmrAxbPar.getTime(); }

private:

    const SimulatorPImpl&       _sim;
    const CadetConstants&       _cc;

    // Parameters set by user
    int _maxKrylov;
    int _maxRestarts;
    int _gramSchmidtType;
    double _schurSafety;        //!< Schur safety factor

    int _wmo;                   //!< Maximum order of weno scheme
    long int _spgmrIters;       //!< counter for SPGMR iterations

    SpgmrMem _spgmrMemBlock;    //!< Pointer to the memory block of the SPGRM linear solver

    JacobianData _jac;

    OmpTimer                _timerFact;             //!< OpenMP timer for jacobian factorization
    OmpTimer                _timerSol;              //!< OpenMP timer for solution of linear equation system
    OmpTimer                _timerSolSpgmr;         //!< OpenMP timer for SPGMR routine
    OmpTimer                _timerSolSpgmrAxb;      //!< OpenMP timer for Matrix times vector routine
    // Timer for parallel regions
    OmpTimer                _timerFactPar;          //!< OpenMP timer for jacobian factorization, parallel part
    OmpTimer                _timerSolPar;           //!< OpenMP timer for solution of linear equation system, parallel part
    OmpTimer                _timerSolSpgmrAxbPar;   //!< OpenMP timer for Matrix times vector routine, parallel part

    const char* SpgmrGetReturnFlagName(int flag) const;

};


} // namespace cadet


#endif /* SCHURSOLVER_HPP_ */
