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

#include <iostream>
#include <sstream>
#include <limits>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <algorithm>
#include <functional>

#include "driver/CadetCSDriver.hpp"
#include "hdf5/HDF5Reader.hpp"
#include "hdf5/HDF5Writer.hpp"
#include "xml/XMLReader.hpp"
#include "xml/XMLWriter.hpp"

#ifndef _WIN32
#include "Gnuplot.hpp"                    // Gnuplot interface for C++; handles POSIX-Pipe-communication with Gnuplot
#endif

namespace cadet {


void printHelp()
{
    using namespace logging;
    log::emit() << "Usage: cadet-cs FILE" << log::endl;
    log::emit() << "Simulate the case defined in FILE" << log::endl;
    log::emit() << "FILE can be either an XML file or an HDF5 file in CADET format" << log::endl;
    log::emit() << "Example: cadet-cs input{.h5|.xml}" << log::endl;
    log::emit() << log::endl;
#ifdef CADET_PUBLISH_BRANCH_INFO
    log::emit() << "This is CADET version " << cadet::getLibraryVersion() << " built from commit " << cadet::getLibraryCommitHash() << " on branch " << cadet::getLibraryBranchRefspec() << log::endl;
#else
    log::emit() << "This is CADET version " << cadet::getLibraryVersion() << " built from commit " << cadet::getLibraryCommitHash() << log::endl;
#endif
    log::emit() << "Report bugs to: cadet@fz-juelich.de" << log::endl;
    log::emit() << "CADET homepage: <http://www.cadet-web.de>" << log::endl;
    log::emit() << "Fork CADET on GitHub: <https://github.com/modsim/CADET>" << log::endl;
}



void pressEnterToContinue()
{
    log::emit<Info>() << "Press ENTER to continue... ";
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}



#ifndef _WIN32
template <typename reader_t, typename writer_t>
void plotWithGnuplot(const CadetCS<reader_t, writer_t>& ccs)
{
    // ============================================================================================================
    //    Plot sensitivitites at column outlet for comp 0 w.r.t. parameter X
    // ============================================================================================================
    std::vector<Gnuplot*> plots;
    for (std::size_t sens = 0; sens < ccs.getNSensParams(); ++sens)
    {
        // Gnuplot sneak preview
        Gnuplot* newPlot = new Gnuplot("points");
        newPlot->cmd("set size 1,1\nset origin 0,0");
        newPlot->set_title("Sensitivities w.r.t. " + ccs.getSensModelParamNames().at(sens));
        newPlot->plot_xy(ccs.getSolutionTimes(), ccs.getSensitivityColumnOutlet().at(sens).at(0), "sens comp 0 " + ccs.getSensModelParamNames().at(sens));
        plots.push_back(newPlot);
    }
    // ============================================================================================================

    // Plot chromatogram for comp 0
    Gnuplot chromatogramPlot("points");
    chromatogramPlot.set_style("lp linestyle 2");
    chromatogramPlot.cmd("set size 1,1\nset origin 0,0");
    chromatogramPlot.set_title("Chromatogram");
    chromatogramPlot.plot_xy(ccs.getSolutionTimes(), ccs.getSolutionColumnOutlet().at(0), "comp 0");

    // Pause
    log::emit<Info>() << "Plots created!" << log::endl;
    pressEnterToContinue();

    // Memory cleanup
    for (std::vector<Gnuplot*>::const_iterator it = plots.begin(); it < plots.end(); ++it)
        delete *it;
    // ============================================================================================================
}
#endif



template <typename reader_t, typename writer_t>
void calculateSensStuff(CadetCS<reader_t, writer_t>& cs)
{
    // Calculate G = \int_0^T t c_N dt using numeric integration
    double G_mean;
    double delta_t;
    double G = 0.0;

    for (std::size_t i = 0; i < cs.getSolutionTimes().size() - 1; ++i)
    {
        G_mean  = (cs.getSolutionColumnOutlet().at(0).at(i) + cs.getSolutionColumnOutlet().at(0).at(i + 1)) * 0.5;
        delta_t = cs.getSolutionTimes().at(i + 1) - cs.getSolutionTimes().at(i);
        G      += cs.getSolutionTimes().at(i) * G_mean * delta_t;
    }
    log::emit<Info>() << "G = " << G << log::endl;


    for (std::size_t sens = 0; sens < cs.getNSensParams(); ++sens)
    {
        // Calculate dG / dp using numeric integration
        double sens_mean;
        double delta_t;
        double dGi_dt = 0.0;
        for (std::size_t i = 0; i < cs.getSolutionTimes().size() - 1; ++i)
        {
            sens_mean = (cs.getSensitivityColumnOutlet().at(0).at(0).at(i) + cs.getSensitivityColumnOutlet().at(0).at(0).at(i + 1)) * 0.5;
            delta_t   = cs.getSolutionTimes().at(i + 1) - cs.getSolutionTimes().at(i);
            dGi_dt   += cs.getSolutionTimes().at(i) * sens_mean * delta_t;
        }
        log::emit<Info>() << "dG/d'" << cs.getSensModelParamNames().at(sens) << "' = " << dGi_dt << log::endl;
    }
}


} // namespace cadet



template <typename reader_t, typename writer_t>
void runAndPlot(cadet::CadetCS<reader_t, writer_t>& ccs)
{
    ccs.run();

#ifndef _WIN32
//    cadet::plotWithGnuplot(ccs);
#endif    
}



int main(int argc, char** argv)
{
    using namespace cadet;

    if (argc != 2)
    {
        printHelp();
        return 1;
    }

    // Depending on the file type use different reader in the driver
    char* fileEnding = strrchr(argv[1], '.');
    fileEnding++; // Move to first character after the "."

    try
    {
        if (strcmp(fileEnding, "h5") == 0)
        {
            CadetCS<HDF5Reader, HDF5Writer> cs(argv[1]);
            runAndPlot(cs);
        }
        else if (strcmp(fileEnding, "xml") == 0)
        {
            CadetCS<XMLReader, XMLWriter> cs(argv[1]);
            runAndPlot(cs);
        }
        else
            log::emit<Error>() << "Wrong input file format! Aborting..." << log::endl;
    }
    catch (const cadet::CadetException& e)
    {
        log::emit<Except>() << e.msg() << log::endl;
    }
    catch (const H5::Exception& e)
    {
        log::emit<Except>() << e.getFuncName() << ": " << e.getDetailMsg() << log::endl;
    }

    return 0;
}

