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

#include <cstdlib>

#include "SchurSolver.hpp"

#include "CadetLogger.hpp"
#include "CadetConvenience.hpp"
#include "AdsorptionModel.hpp"
#include "ChromatographyModel.hpp"
#include "TimeIntegrator.hpp"

namespace cadet {

SchurSolver::SchurSolver(const SimulatorPImpl & sim) :
    _sim(sim),
    _cc(_sim.getCadetConstants()),
    _jac(JacobianData(_cc))
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    // Initialize with default values
    this->configure();
    log::emit<Debug1>() << CURRENT_FUNCTION << ": Configured" << log::endl;

    // ========================================================================================
    //    Setup for the SPGMR iterative linear solver for the schur complement
    // ========================================================================================

    // Create a template vector for the malloc routine of SPGMR
    N_Vector NV_spgmr_tmpl = N_VNew_Serial(_cc.neq_bnd());
    N_VConst(0.0, NV_spgmr_tmpl);

    // Size of allocated memory is either _maxKrylov or _cc.neq_bnd()
    _spgmrMemBlock = SpgmrMalloc(_maxKrylov, NV_spgmr_tmpl);

    N_VDestroy_Serial(NV_spgmr_tmpl);

    log::emit<Debug1>() << CURRENT_FUNCTION << ": SPGMR memory allocated" << log::endl;
    // ========================================================================================

    // Initialize jacobian data object
    _jac.initAllZero();

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


SchurSolver::~SchurSolver()
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    SpgmrFree(_spgmrMemBlock);

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


void SchurSolver::assembleJacCDisc(double alpha)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    double *jacC      = _jac.getJacC()     + _cc.max_wk();
    double *jacC_disc = _jac.getJacCDisc() + _cc.max_wk() + _jac.kl_col();

    for (int i = 0; i < _cc.neq_col(); ++i)
    {
        // Copy the content of jacC to jacC_disc
        for (int j = -_cc.max_wk(); j < _cc.max_wk(); ++j)
            jacC_disc[j] = jacC[j];

        // Add time derivative
        *jacC_disc += alpha;

        // Go to next column of the Jacobian
        jacC      += _jac.ld_jc();
        jacC_disc += _jac.ld_jc_disc();
    }

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}

void SchurSolver::assembleJacPDisc(int pblk, double alpha)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    double *jacP      = _jac.getJacP(pblk);
    double *jacP_disc = _jac.getJacPDisc(pblk) + _jac.kl_par();

    double invBetaP = 1.0 / _sim.getChromatographyModel().getValue<double>(PAR_POROSITY) - 1.0;

    for (int j = 0; j < _cc.npar(); ++j)
    {
        // Mobile phase
        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            for (int r = 0; r < _jac.ld_jp(); ++r)
                jacP_disc[r] = jacP[r];

            jacP_disc[_jac.ku_par()]               += alpha;
            jacP_disc[_jac.ku_par() + _cc.ncomp()] += alpha * invBetaP;

            jacP      += _jac.ld_jp();
            jacP_disc += _jac.ld_jp_disc();
        }

        // Stationary phase
        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            for (int r = 0; r < _jac.ld_jp(); ++r)
                jacP_disc[r] = jacP[r];

            if (_sim.getAdsorptionModel().isDifferential(comp))
                jacP_disc[_jac.ku_par()] += alpha;

            jacP      += _jac.ld_jp();
            jacP_disc += _jac.ld_jp_disc();
        }
    }

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}

void SchurSolver::factorizeJacCDisc()
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    double      *matrix;
    lapackInt_t *piv;

    lapackInt_t n, kl, ku, ldab, flag;

    for (int comp = 0; comp < _cc.ncomp(); ++comp)
    {
        matrix = _jac.getJacCDisc(comp);
        piv    = _jac.getJacCPiv(comp);

        n    = _cc.ncol();
        kl   = _jac.kl_col();
        ku   = _jac.ku_col();
        ldab = _jac.ld_jc_disc();
        flag = 0;

        DGBTRF(&n, &n, &kl, &ku, matrix, &ldab, piv, &flag);

        if (flag != 0)
            log::emit<Error>() << "in " << CURRENT_FUNCTION << "/DGBTRF! flag = " << flag << log::endl;
    }

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}

void SchurSolver::factorizeJacPDisc(int pblk)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    double      *matrix = _jac.getJacPDisc(pblk);
    lapackInt_t *piv    = _jac.getJacPPiv(pblk);

    lapackInt_t n    = _cc.neq_par();
    lapackInt_t kl   = _jac.kl_par();
    lapackInt_t ku   = _jac.ku_par();
    lapackInt_t ldab = _jac.ld_jp_disc();
    lapackInt_t flag = 0;

    DGBTRF(&n, &n, &kl, &ku, matrix, &ldab, piv, &flag);

    if (flag != 0)
        log::emit<Error>() << "in " << CURRENT_FUNCTION << "/DGBTRF! flag = " << flag << log::endl;

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}

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
//
void SchurSolver::solveWithJacCDisc(double *rhs)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    double      *matrix;
    lapackInt_t *piv;

    lapackInt_t n, kl, ku, nrhs, ldab, flag;

    char trans[] = "Transpose";

    for (int comp = 0; comp < _cc.ncomp(); ++comp)
    {
        matrix = _jac.getJacCDisc(comp);
        piv    = _jac.getJacCPiv(comp);

        n    = _cc.ncol();
        kl   = _jac.kl_col();
        ku   = _jac.ku_col();
        nrhs = 1;
        ldab = _jac.ld_jc_disc();
        flag = 0;

        DGBTRS(trans, &n, &kl, &ku, &nrhs, matrix, &ldab, piv, rhs, &n, &flag);

        rhs += _cc.ncol();
        if (flag != 0)
            log::emit<Error>() << "in " << CURRENT_FUNCTION << "/DGBTRS! flag = " << flag << log::endl;
    }

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}

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
//
void SchurSolver::solveWithJacPDisc(int pblk, double *rhs)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    double      *matrix = _jac.getJacPDisc(pblk);
    lapackInt_t *piv    = _jac.getJacPPiv(pblk);

    lapackInt_t n    = _cc.neq_par();
    lapackInt_t kl   = _jac.kl_par();
    lapackInt_t ku   = _jac.ku_par();
    lapackInt_t nrhs = 1;
    lapackInt_t ldab = _jac.ld_jp_disc();
    lapackInt_t flag = 0;

    char trans[] = "Transpose";

    DGBTRS(trans, &n, &kl, &ku, &nrhs, matrix, &ldab, piv, rhs, &n, &flag);

    if (flag != 0)
        log::emit<Error>() << "in " << CURRENT_FUNCTION << "/DGBTRS! flag = " << flag << log::endl;

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}

//  schurComplementTimesVector()
//  multiplies the Schur complement S by the vector v
//  and stores the result in z:
//     z = S * v
//
int SchurSolver::schurComplementTimesVector(void *userData, N_Vector NV_v, N_Vector NV_z)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;
    _timerSolSpgmrAxb.start();

    double *v = NV_DATA_S(NV_v);
    double *z = NV_DATA_S(NV_z);

    N_VScale(1.0, NV_v, NV_z);

    _timerSolSpgmrAxbPar.start();

    #pragma omp parallel for schedule(dynamic)
    for (int pblk = -1; pblk < _cc.npblk(); ++pblk)
    {
        int numel;
        SparseMatrixElement* jdata;

        double* tmp = new double[max(_cc.neq_col(),_cc.neq_par())];

        if (pblk == -1)
        {
            // Column...
            for (int n = 0; n < _cc.neq_col(); ++n)
                tmp[n] = 0.0;

            _jac.sparseMV(_jac.getJacCB(), _jac.numel_jcb(), 1.0, v, tmp);

            solveWithJacCDisc(tmp);

            numel = _jac.numel_jbc();
            jdata = _jac.getJacBC();

        }

        else
        {
            // Particle...
            for (int n = 0; n < _cc.neq_par(); ++n)
                tmp[n] = 0.0;

            _jac.sparseMV(_jac.getJacPB(pblk), _jac.numel_jpb(), 1.0, v, tmp);

            solveWithJacPDisc(pblk, tmp);

            numel = _jac.numel_jbp();
            jdata = _jac.getJacBP(pblk);
        }

        #pragma omp critical
        {
            _timerSolSpgmrAxbPar.stop();
            _jac.sparseMV(jdata, numel, -1, tmp, z);
            _timerSolSpgmrAxbPar.start();
        }

        delete [] tmp;
    }
    _timerSolSpgmrAxbPar.stop();

    _timerSolSpgmrAxb.stop();
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;

    return 0;
}


int SchurSolver::schurSolve(IDAMem IDA_mem, N_Vector NV_rhs, N_Vector weight,
        N_Vector NV_yCur, N_Vector NV_yDotCur, N_Vector NV_resCur)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    SimulatorPImpl* sim = static_cast<SimulatorPImpl*>(IDA_mem->ida_lmem);
    double t     = IDA_mem->ida_tn;
    double alpha = IDA_mem->ida_cj;

    N_Vector NV_tmp = sim->getTimeIntegrator().getNvTemp1();
    double *tmp = NV_DATA_S(NV_tmp);
    double *rhs = NV_DATA_S(NV_rhs);

    // Factorize partial Jacobians only if required!
    if (sim->getTimeIntegrator().factorizeJac() == true)
    {
        _timerFact.start();
        _timerFactPar.start();
        //===============================================================
        // Assemble and factorize discretized system Jacobians
        //===============================================================

        #pragma omp parallel for
        for (int pblk = -1; pblk < _cc.npblk(); ++pblk)
        {
            if (pblk == -1) // column
            {
                assembleJacCDisc(alpha);
            }

            else // particle
            {
                assembleJacPDisc(pblk, alpha);
            }
        }

        factorizeJacCDisc();

        #pragma omp parallel for
        for (int pblk = 0; pblk < _cc.npblk(); ++pblk)
        {
            if (pblk == -1) // column
            {
//                factorizeJacCDisc();
            }

            else // particle
            {
                factorizeJacPDisc(pblk);
            }
        }

/*
        #pragma omp parallel for
        for (int pblk = -1; pblk < _cc.npblk(); ++pblk)
        {
            if (pblk == -1) // column
            {
                assembleJacCDisc(alpha);
                factorizeJacCDisc();
            }

            else // particle
            {
                assembleJacPDisc(pblk, alpha);
                factorizeJacPDisc(pblk);
            }
        }

*/


        //===============================================================
        _timerFactPar.stop();
        _timerFact.stop();
    }


    //===============================================================
    // Part 1: Solve   L * z_hat = rhs
    //===============================================================
    _timerSol.start();

    _timerSolPar.start();

    #pragma omp parallel for
    for (int pblk = -1; pblk < _cc.npblk(); ++pblk)
    {
        if (pblk == -1) // column
            solveWithJacCDisc(_cc.offsetCol(NV_rhs));
        else // particle
            solveWithJacPDisc(pblk, _cc.offsetPar(NV_rhs, pblk));
    }
    _timerSolPar.stop();


    for (int pblk = -1; pblk < _cc.npblk(); ++pblk)
    {
        if (pblk == -1)
            _jac.sparseMV(_jac.getJacBC(), _jac.numel_jbc(), -1.0, _cc.offsetCol(NV_rhs), _cc.offsetBnd(NV_rhs));
        else
            _jac.sparseMV(_jac.getJacBP(pblk), _jac.numel_jbp(), -1.0, _cc.offsetPar(NV_rhs,pblk), _cc.offsetBnd(NV_rhs));
    }
    //===============================================================


    // NOTE: rhs now contains the intermediate solution z_hat !!!


    //===============================================================
    // Initialize temporary storage
    // Column and particle parts <- 0
    for (int i = 0; i < _cc.neq() - _cc.neq_bnd(); ++i)
        tmp[i] = 0.0;

    // Boundary part <- z_b_hat
    for (int i = _cc.neq() - _cc.neq_bnd(); i < _cc.neq(); ++i)
        tmp[i] = rhs[i];
    //===============================================================


    //===============================================================
    // Part 2: Solve   D * z_tilde = z_hat
    //===============================================================

    // column and particle parts remain unchanged!

    // The only thing to be done is the iterative (and approximate)
    // solution of the Schur complement system:
    //     S * z_b_tilde = z_b_hat

    //==== Create init-guess/solution vector by bending pointer =====
    N_Vector NV_x_temp = N_VNewEmpty_Serial(_cc.neq_bnd());
    NV_DATA_S(NV_x_temp) = _cc.offsetBnd(NV_rhs);
    //===============================================================

    //==== Create weight vector by bending pointer ==================
    N_Vector NV_bnd_weight = N_VNewEmpty_Serial(_cc.neq_bnd());
    NV_DATA_S(NV_bnd_weight) = _cc.offsetBnd(weight);
    //===============================================================

    //==== Create right hand side vector by pointer bending =========
    N_Vector NV_b_temp = N_VNewEmpty_Serial(_cc.neq_bnd());
    NV_DATA_S(NV_b_temp) = _cc.offsetBnd(NV_tmp);
    //===============================================================


    double tolerance = _cc.sqrt_neq() * IDA_mem->ida_epsNewt * _schurSafety;

    double res_norm;
    int nli, nps;

    int flag = 0;

    _timerSolSpgmr.start();

    // Here IDA_mem->ida_lmem is a pointer to the current instance of schurSolver.
    // schurComplementTimesVectorWrapper will redirect to it's implementation of
    // schurComplementTimesVector.
    flag = SpgmrSolve(_spgmrMemBlock, IDA_mem->ida_lmem, NV_x_temp, NV_b_temp,
            PREC_NONE, _gramSchmidtType, tolerance, _maxRestarts, NULL,
            NV_bnd_weight, NV_bnd_weight,
            &cadet::TimeIntegrator::schurComplementTimesVectorWrapper,
            NULL, &res_norm, &nli, &nps);

    _timerSolSpgmr.stop();

    if (flag != 0)
        log::emit<Warning>() << "Bad convergence in SpgmrSolve. flag: " << flag << " - "
        << SpgmrGetReturnFlagName(flag) << ". Time = " << t << log::endl;

    _spgmrIters += nli;

    // Free NVector memory space
    N_VDestroy_Serial(NV_b_temp);
    N_VDestroy_Serial(NV_bnd_weight);
    N_VDestroy_Serial(NV_x_temp);

    //===============================================================


    // NOTE: rhs now contains the intermediate solution z_tilde !!!


    //===============================================================
    // Part 3: Solve   U * z = z_tilde
    //===============================================================

    // boundary part remains unchanged!

    _timerSolPar.start();

    #pragma omp parallel for
    for (int pblk = -1; pblk < _cc.npblk(); ++pblk)
    {
        if (pblk == -1) // column
        {
            double *tmp_c = _cc.offsetCol(NV_tmp);
            double *rhs_c = _cc.offsetCol(NV_rhs);

            _jac.sparseMV(_jac.getJacCB(), _jac.numel_jcb(), 1.0, _cc.offsetBnd(NV_rhs), tmp_c);

            solveWithJacCDisc(tmp_c);

            for (int i = 0; i < _cc.neq_col(); ++i)
                rhs_c[i] -= tmp_c[i];
        }

        else // particle
        {
            double *tmp_p = _cc.offsetPar(NV_tmp, pblk);
            double *rhs_p = _cc.offsetPar(NV_rhs, pblk);

            _jac.sparseMV(_jac.getJacPB(pblk), _jac.numel_jpb(), 1.0, _cc.offsetBnd(NV_rhs), tmp_p);

            solveWithJacPDisc(pblk, tmp_p);

            for (int i = 0; i < _cc.neq_par(); ++i)
                rhs_p[i] -= tmp_p[i];
        }
    }
    _timerSolPar.stop();
    //===============================================================

    _timerSol.stop();
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;

    return 0;
}


const char* SchurSolver::SpgmrGetReturnFlagName(int flag) const
{
    switch (flag)
    {
    case  0: return "SPGMR_SUCCESS";            // Converged
    case  1: return "SPGMR_RES_REDUCED";        // Did not converge, but reduced
                                                // norm of residual

    case  2: return "SPGMR_CONV_FAIL";          // Failed to converge
    case  3: return "SPGMR_QRFACT_FAIL";        // QRfact found singular matrix
    case  4: return "SPGMR_PSOLVE_FAIL_REC";    // psolve failed recoverably
    case  5: return "SPGMR_ATIMES_FAIL_REC";    // atimes failed recoverably
    case  6: return "SPGMR_PSET_FAIL_REC";      // pset faild recoverably

    case -1: return "SPGMR_MEM_NULL";           // mem argument is NULL
    case -2: return "SPGMR_ATIMES_FAIL_UNREC";  // atimes returned failure flag
    case -3: return "SPGMR_PSOLVE_FAIL_UNREC";  // psolve failed unrecoverably
    case -4: return "SPGMR_GS_FAIL";            // Gram-Schmidt routine faiuled
    case -5: return "SPGMR_QRSOL_FAIL";         // QRsol found singular R
    case -6: return "SPGMR_PSET_FAIL_UNREC";    // pset failed unrecoverably
    default: return "NO_VALID_FLAG";
    }
}



} // namespace cadet
