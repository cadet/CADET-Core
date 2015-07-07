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

#include <cstddef>
#include <cmath>
#include <algorithm>
#include "active.hpp"


#if defined(ACTIVE_ADOLC) || defined(ACTIVE_SFAD) || defined(ACTIVE_SETFAD)
    ACTIVE_INIT;
#endif

namespace cadet
{
    namespace AD
    {

#if defined(ACTIVE_ADOLC)

        size_t getMaxDirections() { return adtl::ADOLC_numDir; }
        void setMaxDirections(size_t numDir) { adtl::ADOLC_numDir = numDir; }

#elif defined(ACTIVE_SFAD) || defined(ACTIVE_SETFAD)

        size_t getMaxDirections() { return sfad::getGradientSize(); }
        void setMaxDirections(size_t n) { sfad::setGradientSize(n); }

#endif

    } // namespace AD
}  // namespace cadet
