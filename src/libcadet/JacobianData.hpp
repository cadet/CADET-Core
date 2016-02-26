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

#ifndef JACOBIANDATA_HPP_
#define JACOBIANDATA_HPP_

#include "active.hpp"

#include "SimulatorPImpl.hpp"
#include "CadetConvenience.hpp"
#include "SparseMatrixElement.hpp"


namespace cadet
{

class JacobianData
{
public:

    // Constructor
    JacobianData(const CadetConstants& cc);
    // Destructor
    ~JacobianData();

    void initAllZero();
    void initAdForJac(active* yAd, int diagDir) const;
    void setFromAd(const active* resAd, int diagDir);

#ifdef VERIFY_ANALYTICAL_JAC
    void compareWithAd(const active* resAd, int diagDir);
#endif

    inline double* getJacC()                        const { return _jacC; }
    inline double* getJacC(int comp)                const { return _jacC + comp * _ld_jc * _cc.ncol(); }
    inline double* getJacC(int comp, int col)       const { return _jacC + comp * _ld_jc * _cc.ncol() + col * _ld_jc; }
    inline double* getJacP(int pblk)                const { return _jacP + pblk * _ld_jp * _cc.neq_par(); }

    inline double* getJacCDisc()                    const { return _jacC_disc; }
    inline double* getJacCDisc(int comp)            const { return _jacC_disc + comp * _ld_jc_disc * _cc.ncol(); }
    inline double* getJacPDisc(int pblk)            const { return _jacP_disc + pblk * _ld_jp_disc * _cc.neq_par(); }

    inline lapackInt_t* getJacCPiv()                const { return _jacC_piv; }
    inline lapackInt_t* getJacCPiv(int comp)        const { return _jacC_piv + comp * _cc.ncol(); }
    inline lapackInt_t* getJacPPiv(int pblk)        const { return _jacP_piv + pblk * _cc.neq_par(); }

    inline SparseMatrixElement* getJacCB()          const { return _jacCB; }
    inline SparseMatrixElement* getJacBC()          const { return _jacBC; }

    inline SparseMatrixElement& getJacCB(int eq)    const { return _jacCB[eq]; }
    inline SparseMatrixElement& getJacBC(int eq)    const { return _jacBC[eq]; }

    inline SparseMatrixElement* getJacPB(int pblk)  const { return _jacPB + pblk * _numel_jpb; }
    inline SparseMatrixElement* getJacBP(int pblk)  const { return _jacBP + pblk * _numel_jbp; }

    inline int kl_col()       const { return _kl_col; }
    inline int ku_col()       const { return _ku_col; }
    inline int ld_jc()        const { return _ld_jc; }
    inline int ld_jc_disc()   const { return _ld_jc_disc; }

    inline int kl_par()       const { return _kl_par; }
    inline int ku_par()       const { return _ku_par; }
    inline int ld_jp()        const { return _ld_jp; }
    inline int ld_jp_disc()   const { return _ld_jp_disc; }

    inline int numel_jcb()    const { return _numel_jcb; }
    inline int numel_jbc()    const { return _numel_jbc; }
    inline int numel_jpb()    const { return _numel_jpb; }
    inline int numel_jbp()    const { return _numel_jbp; }

    inline int jacAdDirs()    const { return _jacAdDirs; }
    inline int diagDir()      const { return _diagDir; }

    inline double maxDiffC()  const { return _maxDiffC; }
    inline double maxDiffP()  const { return _maxDiffP; }

    // Simple multiplication with a sparse matrix b += alpha * A * x
    inline void sparseMV(SparseMatrixElement* A, int numel, double alpha, double* x, double* b) const {
        for (int i = 0; i < numel; ++i)
            b[A[i].row()] += alpha * A[i].value() * x[A[i].col()];
    }

private:

    double* _jacC; // Column level jacobian
    double* _jacP; // Particle level jacobian
    SparseMatrixElement* _jacCB;
    SparseMatrixElement* _jacBC;
    SparseMatrixElement* _jacPB;
    SparseMatrixElement* _jacBP;

    // BDF discretized system Jacobians
    double* _jacC_disc;
    double* _jacP_disc;

    // Pivot element matrices for factorization
    // of discretized system Jacobians
    lapackInt_t* _jacC_piv;
    lapackInt_t* _jacP_piv;


    int _kl_col;       //!< Number of subdiagonals in column block
    int _ku_col;       //!< Number of superdiagonals in column block
    int _ld_jc;        //!< Leading dimension for column level jacobian
    int _ld_jc_disc;   //!< Leading dimension for column level discretized jacobian

    int _npblk;        //!< Number of paticle blocks ( = ncol)
    int _kl_par;       //!< Number of subdiagonals in particle blocks
    int _ku_par;       //!< Number of superdiagonals in particle blocks
    int _ld_jp;        //!< Leading dimension for particle level jacobian
    int _ld_jp_disc;   //!< Leading dimension for particle level discretized jacobian

    int _numel_jcb;    //!< Number of elements in the column-boundary block
    int _numel_jbc;    //!< Number of elements in the boundary-column block
    int _numel_jpb;    //!< Number of elements in the particle-boundary block
    int _numel_jbp;    //!< Number of elements in the boundary-particle block

    int _jacAdDirs;    //!< Number of directions needed for Jacobian assembly
    int _diagDir;      //!< Direction of actives, which corresponds to Jacobian diagonal entries

    double _maxDiffC;   //!< Stores the maximum difference in column jacobian when comparing analytic and AD implementation
    double _maxDiffP;   //!< Stores the maximum difference in particle jacobian when comparing analytic and AD implementation

    const CadetConstants& _cc;

};

} // namespace cadet

#endif  //JACOBIANDATA_HPP_
