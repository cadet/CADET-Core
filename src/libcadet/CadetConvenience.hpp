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

#ifndef CADETCONVENIENCE_HPP_
#define CADETCONVENIENCE_HPP_

#include <string>

#include <nvector/nvector_serial.h>

#include "SparseMatrixElement.hpp"

// nullptr macro should be removed when compiled according to C++11 standard
//#define nullptr NULL


namespace cadet
{

class CadetConstants
{
public:
    /// \brief Set all default constants that depend on ncomp, ncol, npar
    CadetConstants(int ncomp, int ncol, int npar);

    // Inline getter member functions ...
    // ... for basic discretization parameters
    inline int    ncomp()     const { return _ncomp; }
    inline int    ncol()      const { return _ncol; }
    inline int    npar()      const { return _npar; }
    inline int    nbnd()      const { return _nbnd; }

    inline int    nstatec()   const { return _nstatec; }
    inline int    nstateq()   const { return _nstateq; }
    inline int    nstate()    const { return _nstate; }
    inline int    npblk()     const { return _ncol; }

    inline int    neq_col()   const { return _neq_col; }
    inline int    neq_bnd()   const { return _neq_bnd; }
    inline int    neq_par()   const { return _neq_par; }
    inline int    neq()       const { return _neq; }

    inline double sqrt_neq()  const { return _sqrt_neq; }
    inline int    max_wk()    const { return _max_wk; }


    // ========================================================================================
    //    Helper functions to access N_Vector parts (e.g. state vector) easily
    // ========================================================================================
    inline double* offsetCol(N_Vector vec) const {
        return NV_DATA_S(vec);
    }
    inline double* offsetPar(N_Vector vec, int pblk) const {
        return (NV_DATA_S(vec) + _neq_col + (pblk * _neq_par));
    }
    inline double* offsetBnd(N_Vector vec) const {
        return (NV_DATA_S(vec) + _neq_col + (_npblk * _neq_par));
    }
    // ========================================================================================


    // ========================================================================================
    //    Helper functions to access certain fields in a part of a vector
    //       (e.g. the state vector) easily
    // ========================================================================================
    // Mobile phase c in the column
    template <typename T>
    inline T& colC(T* partBegin, const int colCell, const int comp) const {
        return *(partBegin + colCell + (comp * _ncol));
    }
    // const version
    template <typename T>
    inline const T& colC(const T* partBegin, const int colCell, const int comp) const {
        return *(partBegin + colCell + (comp * _ncol));
    }

    // Mobile phase c in the particle
    template <typename T>
    inline T& parC(T* partBegin, const int parCell, const int comp) const {
        return *(partBegin + comp + (parCell * _nstate * _ncomp));
    }
    // const version
    template <typename T>
    inline const T& parC(const T* partBegin, const int parCell, const int comp) const {
        return *(partBegin + comp + (parCell * _nstate * _ncomp));
    }

    // Stationary phase q in the particle
    template <typename T>
    inline T& parQ(T* partBegin, const int parCell, const int comp) const {
        return *(partBegin + _ncomp + comp + (parCell * _nstate * _ncomp));
    }
    // const version
    template <typename T>
    inline const T& parQ(const T* partBegin, const int parCell, const int comp) const {
        return *(partBegin + _ncomp + comp + (parCell * _nstate * _ncomp));
    }

    // Boundary fluxes j on the boundary
    template <typename T>
    inline T& bndJ(T* partBegin, const int bndCell, const int comp) const {
        return *(partBegin + bndCell + (comp * _nbnd));
    }
    // const version
    template <typename T>
    inline const T& bndJ(const T* partBegin, const int bndCell, const int comp) const {
        return *(partBegin + bndCell + (comp * _nbnd));
    }
    // ========================================================================================

private:
    // ---------------------------------------------------------------------------
    //   Basic discretization parameters
    // ---------------------------------------------------------------------------
    int _ncomp;     //!< Number of components
    int _ncol;      //!< Number of spatial discretization on column level
    int _npar;      //!< Number of spatial discretization on particle level
    int _nbnd;      //!< Number of spatial discretization on boundary level

    int _nstatec;   //!< Number of mobile states a species in the bead can be in
    int _nstateq;   //!< Number of bound states a species in the bead can be in (e.g. q_1, q_2 ...)
    int _nstate;    //!< Number of overall states a species in the bead can be in
    int _npblk;     //!< Number of particle blocks

    int _neq_col;   //!< Number of equations, column level
    int _neq_par;   //!< Number of equations, particle interior (for one particle)
    int _neq_bnd;   //!< Number of equations, particle surfaces
    int _neq;       //!< Total number of equations (length of state vectors)

    double _sqrt_neq;   //!< Sqare root of the total number of equations
    const int _max_wk;  //!< Maximum possible order of Weno Scheme - hardcoded to be 3

};

///todo Functions not assigned to a class yet:


} // namespace cadet

#endif /* CADETCONVENIENCE_HPP_ */
