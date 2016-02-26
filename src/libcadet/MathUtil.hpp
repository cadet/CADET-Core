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

#ifndef MATHUTIL_HPP_
#define MATHUTIL_HPP_

namespace cadet 
{

#if defined(ACTIVE_ADOLC)
    template <typename real_t> real_t sqr(const real_t& x) { return x * x; }
#else
	inline double sqr(const double x) { return x * x; }
#endif

}

#endif /* MATHUTIL_HPP_ */
