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

#ifndef SIMPLEPROGRESS_HPP_
#define SIMPLEPROGRESS_HPP_

#include <sstream>
#include <iomanip>

#include "CadetLogger.hpp"

namespace cadet {

class SimpleProgress
{
public:
    SimpleProgress(double min, double max, int nsteps)
    {
        _min = min;
        _max = max;

        _steps = new double[nsteps];

        for (int i = 0; i < nsteps; ++i)
            _steps[i] = min + (i + 1) * (max - min) / nsteps;

        _nextStep = 0;

        std::cout << std::flush;
    }

    ~SimpleProgress() { delete [] _steps; }

    inline void printIfHitNext(double current)
    {
        if (current >= _steps[_nextStep])
        {
            _oss.str("");
            _oss << std::fixed << std::setprecision(0) << (current - _min) / (_max - _min) * 100.0 << "%   ";
            log::emit() << _oss.str();
            std::cout << std::flush;
            _nextStep++;
        }
    }

private:

    std::ostringstream _oss;
    int     _nextStep;
    double  _min;
    double  _max;
    double* _steps;

};

}  // namespace cadet

#endif  // SIMPLEPROGRESS_HPP_
