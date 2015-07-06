// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2015: Eric von Lieres¹, Joel Andersson¹,
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

#include <string>
#include <cmath>

#include "CadetConvenience.hpp"

namespace cadet
{

CadetConstants::CadetConstants(int ncomp, int ncol, int npar, int nsec) :
    _max_wk(3)
{
    // Set basic constants
    _ncomp = ncomp;
    _ncol = ncol;
    _npar = npar;
    _nbnd = ncol;
    _nsec = nsec;

    _nstatec = 1;
    _nstateq = 1;
    _nstate = _nstatec + _nstateq;
    _npblk   = ncol;

    // Set number of equations
    _neq_col = _ncol * _ncomp;
    _neq_bnd = _ncol * _ncomp;
    _neq_par = (_nstatec + _nstateq) * _npar * _ncomp;
    _neq = _neq_col + _ncol * _neq_par + _neq_bnd;

    // Misc
    _sqrt_neq = sqrt((double) _neq);

}

}  // namespace cadet
