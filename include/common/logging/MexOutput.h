// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2015: Eric von Lieres¹, Joel Andersson,
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

#if !defined(__MEXOutput_h__) && defined(MATLAB_MEX_FILE)
#define __MEXOutput_h__

#include <mex.h>

namespace logging {

    /*! \brief Provides an interface to the MEX output
     */
    class MexOutput {
        public:
            /*! \brief operator that can output a simple character.
             *
             * \param c the character that needs to be outputed
             * \return a reference to itself allowing chaining of
             *         opertor<< calls.
             */
            MexOutput & operator<<(const char c) {
                mexPrintf("%c", c);
                return *this;
            }
    };

} /* logging */

#endif
