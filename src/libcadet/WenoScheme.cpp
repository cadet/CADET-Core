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

#include <vector>
#include <algorithm>

#include "WenoScheme.hpp"
#include "TimeIntegrator.hpp"

namespace cadet {


// Initialization of static Weno coefficients
const double WenoScheme::_wenoD2[2]   = { 2.0/3.0, 1.0/3.0 };

const double WenoScheme::_wenoC2[2*2] = { 0.5, -0.5,
                                          0.5,  1.5 };
const double WenoScheme::_wenoJbvv2[2*3*3] =
{0, 2,   0,-2,   0, 0,
 0,-2,   2, 2,  -2, 0,
 0, 0,  -2, 0,   2, 0};

const double WenoScheme::_wenoD3[3]   = { 0.3, 0.6, 0.1 };

const double WenoScheme::_wenoC3[3*3] = { 1.0/3.0, -1.0/6.0,  1.0/3.0,
                                          5.0/6.0,  5.0/6.0, -7.0/6.0,
                                         -1.0/6.0,  1.0/3.0, 11.0/6.0  };
const double WenoScheme::_wenoJbvv3[3*5*5] = // Used to generate Jbv: vec(Jbv) = A*v
{0,0,8.0/3,   0,0,-19.0/3,        0,0,11.0/3,            0,0,0,             0,0,0,
 0,0,-19.0/3, 0,8.0/3,50.0/3,     0,-13.0/3,-31.0/3,     0,5.0/3,0,         0,0,0,
 0,0,11.0/3,  0,-13.0/3,-31.0/3,  20.0/3,26.0/3,20.0/3, -31.0/3,-13.0/3,0,  11.0/3,0,0,
 0,0,0,       0,5.0/3,0,         -31.0/3,-13.0/3,0,      50.0/3,8.0/3,0,   -19.0/3,0,0,
 0,0,0,       0,0,0,              11.0/3,0,0,           -19.0/3,0,0,         8.0/3,0,0};


WenoScheme::WenoScheme(const SimulatorPImpl& sim):
    _sim(sim),
    _cc(sim.getCadetConstants()),
    _weps(std::vector<double>(_cc.ncomp(), 0.0)),
    _samplePointsInletConc(101)
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    // Initialize with default values
    this->configure();
    log::emit<Debug1>() << CURRENT_FUNCTION << ": Configured" << log::endl;

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}

void WenoScheme::setWenoEpsilons(N_Vector NV_y, const std::vector<double>& sectionTimes, int sec)
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    // ==========================================================================
    //    Set WENO epsilon
    // ==========================================================================

    ///todo describe the task that this routine is doing!

    // Get access to the column part of the stae vector
    double* cdata = _cc.offsetCol(NV_y);

    // Sample the inlet concentrations of the current section for the largest value
    double stepSize = (sectionTimes.at(sec + 1) - sectionTimes.at(sec)) / double(_samplePointsInletConc);
    double t;

    std::vector<double> concInlet(_cc.ncomp(), 0.0);

    // Create dummy vectors for the inletConcentration function - we dont need derivatives here!
    int msip = _sim.getTimeIntegrator().getMaxSensInletParams();
    std::vector<bool> dummyInletParamIsSensitive = std::vector<bool>(msip, false);
    std::vector<std::vector<double> > dummyDeriv = std::vector<std::vector<double> >(_cc.ncomp(), std::vector<double>(msip, 0.0));

    for (int comp = 0; comp < _cc.ncomp(); ++comp)  // Iterate over components
    {
        // Find largest concentration in the current state vector
        double maxC = _cc.colC<double>(cdata, 0, comp);

        for (int col = 1; col < _cc.ncol(); ++col)  // Iterate over column cells
            if (_cc.colC<double>(cdata, col, comp) > maxC)
                maxC = _cc.colC<double>(cdata, col, comp);

        for (int step = 0; step < _samplePointsInletConc; ++step)
        {
            t = sectionTimes.at(sec) + stepSize * double(step);

            _sim.getTimeIntegrator().inletConcentration(t, sec, concInlet, dummyInletParamIsSensitive, dummyDeriv);

            if (concInlet.at(comp) > maxC)
                maxC = concInlet.at(comp);
        }

        ///todo rethink this!
        // Set references (i.e. typical concentrations) for the weno epsilons
        _weps.at(comp) = _wenoEps * std::max(pow(maxC,2), 1.0e-6);
    }
    // ==========================================================================

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}

} // namespace cadet
