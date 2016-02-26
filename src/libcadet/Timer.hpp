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

#ifndef TIMER_HPP_
#define TIMER_HPP_

#ifdef _OPENMP
#include <omp.h>
#endif

namespace cadet {

class OmpTimer
{
public:
    // Constructor
    OmpTimer() : _startTime(0.0), _elapsedTime(0.0) {}

    inline void start() {
#ifdef _OPENMP
        _startTime    = omp_get_wtime();
#endif
        }
    inline void stop()  {
#ifdef _OPENMP
        _elapsedTime += omp_get_wtime() - _startTime;
#endif
        }
    inline double getTime() const {
#ifdef _OPENMP
        return _elapsedTime;
#else
        return -1.0;
#endif
    }

private:

    double _startTime;
    double _elapsedTime;
};

}  // namespace cadet



#endif /* TIMER_HPP_ */
