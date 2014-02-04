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

#include "Simulator.hpp"
#include "Gnuplot.hpp"

class inletParams : public cadet::InletBase
{
public:
    void inletConcentration(double t, int sec,
            std::vector<double>& inletConc,
            const std::vector<bool>& inletParamIsSensitive,
            std::vector<std::vector<double> >& dInletConc_dp)
    {
        inletConc.at(0) = 7.14e-3;  // Always return a constant inlet concentration for comp 0
    }
};

int main(int argc, char** argv)
{
    using namespace cadet;
    int ncomp = 1;
    std::cout << "Creating simulator... ";
    Simulator sim(ncomp, 16, 4, 1, MULTI_COMPONENT_LANGMUIR, GENERAL_RATE_MODEL);
    std::cout << "done." << std::endl;

    sim.setParameterSectionDependent(VELOCITY,          false);
    sim.setParameterSectionDependent(COL_DISPERSION,    false);
    sim.setParameterSectionDependent(PAR_DIFFUSION,     false);
    sim.setParameterSectionDependent(FILM_DIFFUSION,    false);
    sim.setParameterSectionDependent(PAR_SURFDIFFUSION, false);

    sim.setParameterValue(5.75e-4,  VELOCITY);                // Chromatography model parameters
    sim.setParameterValue(0.014,    COL_LENGTH);
    sim.setParameterValue(0.37,     COL_POROSITY);
    sim.setParameterValue(5.75e-8,  COL_DISPERSION);
    sim.setParameterValue(4.5e-5,   PAR_RADIUS);
    sim.setParameterValue(0.75,     PAR_POROSITY);
    sim.setParameterValue(6.07e-11, PAR_DIFFUSION,  0);
    sim.setParameterValue(6.9e-6,   FILM_DIFFUSION, 0);

    sim.setKineticAdsorptionModel(true);                      // Adsorption model parameters
    sim.setParameterValue(1.14,     MCL_KA,   0);
    sim.setParameterValue(0.002,    MCL_KD,   0);
    sim.setParameterValue(4.88,     MCL_QMAX, 0);

    std::vector<double> init(ncomp, 0.0);                     // Create and set initial conditions vector
    std::cout << "Init simulator... ";
    sim.initialize(init, init);
    std::cout << "done." << std::endl;


    inletParams inPar;
    sim.setInletProfile(&inPar);                              // Register class that implements the inlet concentration function

    init.clear();                                             // Make init the section times vector
    init.push_back(0);
    init.push_back(12000);
    sim.setSectionTimes(init);                                // Set section times vector

    std::cout << "Start simulation... ";
    sim.integrate();                                          // Solve the system
    std::cout << "done." << std::endl;

    std::vector<double> solutionTimes;
    sim.getSolutionTimes(solutionTimes);                      // Get times of written solution
    std::vector<double> chromatogram;
    sim.getSolutionColumnOutlet(chromatogram, 0);                     // Get chromatogram for comp 0

#if !(defined(_WIN32) || defined(__WIN32__) || defined(WIN32))
    std::cout << "Plotting solution (type 1 to quit) ... ";
    Gnuplot chromatogramPlot("points");
    chromatogramPlot.plot_xy(solutionTimes, chromatogram, "Component 1");
#endif

    std::cin >> ncomp;
    return 0;
}
