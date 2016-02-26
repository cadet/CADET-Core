// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2016: Eric von Lieres¹, Joel Andersson¹,
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

#include <algorithm>

#include "JacobianData.hpp"
#include "CadetLogger.hpp"

namespace cadet {

JacobianData::JacobianData(const CadetConstants& cc) :
    _cc(cc)
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    // Set column jacobian related constants
    _kl_col = _cc.max_wk() - 1;
    _ku_col = _cc.max_wk();
    _ld_jc = _kl_col + 1 + _ku_col;
    _ld_jc_disc = 2 * _kl_col + _ku_col + 1;

    // Set particle jacobian related constants
    _npblk = _cc.npblk();
    _kl_par = 3 * _cc.ncomp();
    _ku_par = 2 * _cc.ncomp();
    _ld_jp = _kl_par + 1 + _ku_par;
    _ld_jp_disc = 2 * _kl_par + _ku_par + 1;

    // Set boundary jacobian blocks related constants
    _numel_jcb = _cc.neq_col();
    _numel_jbc = _cc.neq_col();
    _numel_jpb = _cc.ncomp();
    _numel_jbp = _cc.ncomp();

    // Set AD assembled jacobian constants
    _jacAdDirs = std::max(_ku_col, _ku_par) + std::max(_kl_col, _kl_par) + 1;
    _diagDir   = std::max(_ku_col, _ku_par);

    log::emit<Debug1>() << CURRENT_FUNCTION << ": AD directions needed for Jacobian: " << _jacAdDirs << log::endl;
    log::emit<Debug1>() << CURRENT_FUNCTION << ": Direction corresponding to main diagonal: " << _diagDir << log::endl;

    // Allocate column jacobian block
    // no. diagonals * no. equations in column
    _jacC       = new double[_ld_jc      * _cc.neq_col()];
    _jacC_disc  = new double[_ld_jc_disc * _cc.neq_col()];
    _jacC_piv   = new lapackInt_t[_cc.neq_col()];

    // Allocate particle jacobian blocks
    // no. diagonals * no. equations in particle * no. particle blocks
    _jacP       = new double[_ld_jp      * _cc.neq_par() * _npblk];
    _jacP_disc  = new double[_ld_jp_disc * _cc.neq_par() * _npblk];
    _jacP_piv   = new lapackInt_t[_cc.neq_par() * _npblk];


    // Allocate boundary jacobian blocks
    _jacCB = new SparseMatrixElement[_numel_jcb];
    _jacBC = new SparseMatrixElement[_numel_jbc];
    _jacPB = new SparseMatrixElement[_npblk * _numel_jpb];
    _jacBP = new SparseMatrixElement[_npblk * _numel_jbp];

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}

// destructor. frees all internal memory
JacobianData::~JacobianData()
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    delete [] _jacC;
    delete [] _jacC_disc;
    delete [] _jacC_piv;

    delete [] _jacP;
    delete [] _jacP_disc;
    delete [] _jacP_piv;

    delete [] _jacCB;
    delete [] _jacBC;
    delete [] _jacPB;
    delete [] _jacBP;

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


void JacobianData::initAllZero()
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    double initValueDouble = 0.0;
    int initValueInt = 0;

    // Initialize all column related arrays
    for (int i = 0; i < (_ld_jc      * _cc.neq_col()); ++i)
        _jacC[i] = initValueDouble;
    for (int i = 0; i < (_ld_jc_disc * _cc.neq_col()); ++i)
        _jacC_disc[i] = initValueDouble;
    for (int i = 0; i < (_cc.neq_col()); ++i)
        _jacC_piv[i] = initValueInt;

    // Initialize all particle related arrays
    for (int i = 0; i < (_ld_jp      * _cc.neq_par() * _npblk); ++i)
        _jacP[i] = initValueDouble;
    for (int i = 0; i < (_ld_jp_disc * _cc.neq_par() * _npblk); ++i)
        _jacP_disc[i] = initValueDouble;
    for (int i = 0; i < (_cc.neq_par() * _npblk); ++i)
        _jacP_piv[i] = initValueInt;

    // Initialize sparse matrices
    for (int i = 0; i < (_numel_jcb); ++i)
        _jacCB[i].setElement(initValueInt, initValueInt, initValueDouble);

    for (int i = 0; i < (_numel_jbc); ++i)
        _jacBC[i].setElement(initValueInt, initValueInt, initValueDouble);

    for (int i = 0; i < (_numel_jpb); ++i)
        _jacPB[i].setElement(initValueInt, initValueInt, initValueDouble);

    for (int i = 0; i < (_numel_jbp); ++i)
        _jacBP[i].setElement(initValueInt, initValueInt, initValueDouble);


    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


void JacobianData::initAdForJac(active* yAd, int diagDir) const
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    // Start with diagonal Jacobian element
    int dir = diagDir;

    // Column
    for (int eqc = 0; eqc < _cc.neq_col(); ++eqc)  // Loop over column equations
    {
        yAd[eqc].setADValue(dir, 1.0);

        if (dir == diagDir + _kl_col)
            dir = diagDir - _ku_col;
        else
            dir++;
    }

    // Particles
    for (int pblk = 0; pblk < _cc.npblk(); ++pblk) // Loop over particle blocks
    {
        dir = diagDir;

        for (int eqp = 0; eqp < _cc.neq_par(); ++eqp) // Loop over particle equations
        {
            yAd[_cc.neq_col() + pblk * _cc.neq_par() + eqp].setADValue(dir, 1.0);

            if (dir == diagDir + _kl_par)
                dir = diagDir - _ku_par;
            else
                dir++;
        }
    }

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


void JacobianData::setFromAd(const active* resAd, int diagDir)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    //=============================================================================
    // Copy Jacobian entries from resAd to band Jacobian data structures...
    //=============================================================================
    int dir;

    // Column part
    for (int eqc = 0; eqc < _cc.neq_col(); ++eqc)  // Loop over column equations
    {
        dir = diagDir - _ku_col + eqc % _ld_jc;

        for (int diag = 0; diag < _ld_jc; ++diag)  // Loop over diagonals
        {
            _jacC[eqc * _ld_jc + diag] = resAd[eqc].getADValue(dir);

            if (dir == diagDir + _kl_col)
                dir = diagDir - _ku_col;
            else
                dir++;
        }
    }

    // Particle part
    for (int pblk = 0; pblk < _cc.npblk(); ++pblk)  // Loop over particle blocks
        for (int eqp = 0; eqp < _cc.neq_par(); ++eqp)  // Loop over particle equations
        {
            dir = diagDir - _ku_par + eqp % _ld_jp;

            for (int diag = 0; diag < _ld_jp; ++diag)  // Loop over diagonals
            {
                *(getJacP(pblk) + eqp * _ld_jp + diag) = resAd[_cc.neq_col() + pblk * _cc.neq_par() + eqp].getADValue(dir);

                if (dir == diagDir + _kl_par)
                    dir = diagDir - _ku_par;
                else
                    dir++;
            }
        }
    //=============================================================================

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}

#ifdef VERIFY_ANALYTICAL_JAC
void JacobianData::compareWithAd(const active* resAd, int diagDir)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    //=============================================================================
    // Compare Jacobian entries from resAd with analytically computed Jacobian
    // from band jacobian data structures...
    //=============================================================================
    int dir;
    double diff = 0.0;
    double maxC = 0.0;
    double maxP = 0.0;

    // Column part
    for (int eqc = 0; eqc < _cc.neq_col(); ++eqc)  // Loop over column equations
    {
        dir = diagDir - _ku_col + eqc % _ld_jc;

        for (int diag = 0; diag < _ld_jc; ++diag)  // Loop over diagonals
        {
            diff = abs(_jacC[eqc * _ld_jc + diag] - resAd[eqc].getADValue(dir));
            maxC = (maxC < diff) ? diff : maxC;

            if (dir == diagDir + _kl_col)
                dir = diagDir - _ku_col;
            else
                dir++;
        }
    }
    _maxDiffC = (_maxDiffC < maxC) ? maxC : _maxDiffC;

    // Particle part
    for (int pblk = 0; pblk < _cc.npblk(); ++pblk)  // Loop over particle blocks
        for (int eqp = 0; eqp < _cc.neq_par(); ++eqp)  // Loop over particle equations
        {
            dir = diagDir - _ku_par + eqp % _ld_jp;

            for (int diag = 0; diag < _ld_jp; ++diag)  // Loop over diagonals
            {
                diff  = abs(*(getJacP(pblk) + eqp * _ld_jp + diag) - resAd[_cc.neq_col() + pblk * _cc.neq_par() + eqp].getADValue(dir));
                maxP = (maxP < diff) ? diff : maxP;

                if (dir == diagDir + _kl_par)
                    dir = diagDir - _ku_par;
                else
                    dir++;
            }
        }
    _maxDiffP = (_maxDiffP < maxP) ? maxP : _maxDiffP;
    //=============================================================================

    log::emit<Debug1>() << "JacDiff: Col = " << maxC << " Part = " << maxP << log::endl;
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}
#endif

} // namespace cadet
