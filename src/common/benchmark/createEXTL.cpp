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

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdlib.h>

#include "CadetEnumeration.hpp"
#include "hdf5/HDF5Writer.hpp"
#include "CadetLogger.hpp"

using namespace cadet;


void printHelp()
{
    log::emit() << "Usage: createEXTL <FILE>" << log::endl;
    log::emit() << "Create an HDF5 input file for externally dependent single-component-langmuir benchmark case." << log::endl;
    log::emit() << "To be used with cadet-cs from the Cromatography Analysis and DEsign Toolkit (CADET)." << log::endl;
    log::emit() << "When called with no arguments, default file name is EXTL.h5." << log::endl;
    log::emit() << "Example: createEXTL newfile.h5" << log::endl;
    log::emit() << log::endl;
    log::emit() << "Report bugs to: cadet@fz-juelich.de" << log::endl;
    log::emit() << "CADET homepage: <http://www.cadet-web.de>" << log::endl;
    log::emit() << "Fork CADET on GitHub: <https://github.com/modsim/CADET>" << log::endl;
}


int main(int argc, char** argv)
{
    std::string fileName("EXTL.h5");  // default output filename
    if (argc == 2) fileName = argv[1];
    if (argc > 2)
    {
        printHelp();
        return 1;
    }

    ChromatographyType  chromatography_type     = GENERAL_RATE_MODEL;

    // Model parameters
    AdsorptionType      adsorption_type         = EXTERNAL_LANGMUIR;
    int                 ncomp                   = 2;

    std::vector<double> init_c(ncomp, 0.0);     // Initialize vector
    std::vector<double> init_q(ncomp, 0.0);     // Initialize vector

                        init_c.at(0)            = 0.0;
                        init_q.at(0)            = 0.0;

                        init_c.at(1)            = 0.0;
                        init_q.at(1)            = 0.0;

    double              col_dispersion          = 6e-7;
    double              col_length              = 0.1;
    double              col_porosity            = 0.4;
    std::vector<double> film_diffusion(ncomp, 0.0);         // Initialize vector
                        film_diffusion.at(0)    = 1.0;
                        film_diffusion.at(1)    = 1.0;
    std::vector<double> par_diffusion(ncomp, 0.0);          // Initialize vector
                        par_diffusion.at(0)     = 1e-6;
                        par_diffusion.at(1)     = 1e-6;
    double              par_porosity            = 0.5;
    double              par_radius              = 4.5e-5;
    std::vector<double> par_surfdiffusion(ncomp, 0.0);      // Initialize vector
                        par_surfdiffusion.at(0) = 0.0;
                        par_surfdiffusion.at(1) = 0.0;
    double              velocity                = 7.3611e-4;

    // Adsorption parameters

    std::vector<double> extl_ka(ncomp, 0.0);       // Initialize vector
    std::vector<double> extl_ka_T(ncomp, 0.0);     // Initialize vector
    std::vector<double> extl_ka_TT(ncomp, 0.0);    // Initialize vector
    std::vector<double> extl_ka_TTT(ncomp, 0.0);   // Initialize vector
    std::vector<double> extl_kd(ncomp, 0.0);       // Initialize vector
    std::vector<double> extl_kd_T(ncomp, 0.0);     // Initialize vector
    std::vector<double> extl_kd_TT(ncomp, 0.0);    // Initialize vector
    std::vector<double> extl_kd_TTT(ncomp, 0.0);   // Initialize vector
    std::vector<double> extl_qmax(ncomp, 0.0);     // Initialize vector
    std::vector<double> extl_qmax_T(ncomp, 0.0);   // Initialize vector
    std::vector<double> extl_qmax_TT(ncomp, 0.0);  // Initialize vector
    std::vector<double> extl_qmax_TTT(ncomp, 0.0); // Initialize vector
    int                 is_kinetic              = 0;

                        extl_ka.at(0)            = 1.0;
                        extl_ka_T.at(0)          = 0.0;
                        extl_ka_TT.at(0)         = 0.0;
                        extl_ka_TTT.at(0)        = 0.0;
                        extl_kd.at(0)            = 4.0821760312; // scaled(308.15K): 0.1665652443; // non-scaled: 4.0821760312;
                        extl_kd_T.at(0)          = -0.0127068336;
                        extl_kd_TT.at(0)         = 0.0;
                        extl_kd_TTT.at(0)        = 0.0;
                        extl_qmax.at(0)          = -10.8607615182;// scaled(308.15K): 2.4776360544; // non-scaled: -10.8607615182;
                        extl_qmax_T.at(0)        = 0.0432854051;
                        extl_qmax_TT.at(0)       = 0.0;
                        extl_qmax_TTT.at(0)      = 0.0;

                        extl_ka.at(1)            = 1.0;
                        extl_ka_T.at(1)          = 0.0;
                        extl_ka_TT.at(1)         = 0.0;
                        extl_ka_TTT.at(1)        = 0.0;
                        extl_kd.at(1)            = 2 *  4.0821760312;
                        extl_kd_T.at(1)          = 2 * -0.0127068336;
                        extl_kd_TT.at(1)         = 0.0;
                        extl_kd_TTT.at(1)        = 0.0;
                        extl_qmax.at(1)          = -10.8607615182;
                        extl_qmax_T.at(1)        = 0.0432854051;
                        extl_qmax_TT.at(1)       = 0.0;
                        extl_qmax_TTT.at(1)      = 0.0;

    // Inlet

    int                 nsec                    = 1;
    std::vector<double> section_times(nsec+1, 0.0);   // Initialize vector
    std::vector<int> section_continuity;
    
                        section_times.at(0)     = 0.0;
                        section_times.at(1)     = 3500.0; // [s]
    std::vector<std::vector<double> > const_coeff(nsec, std::vector<double>(ncomp, 0.0));      // Initialize vector
    std::vector<std::vector<double> > lin_coeff  (nsec, std::vector<double>(ncomp, 0.0));      // Initialize vector
    std::vector<std::vector<double> > quad_coeff (nsec, std::vector<double>(ncomp, 0.0));      // Initialize vector
    std::vector<std::vector<double> > cube_coeff (nsec, std::vector<double>(ncomp, 0.0));      // Initialize vector

    // Sec 1
                        const_coeff.at(0).at(0) = 0.05;                  // component 1
                        lin_coeff  .at(0).at(0) = 0.0;
                        quad_coeff .at(0).at(0) = 0.0;
                        cube_coeff .at(0).at(0) = 0.0;

                        const_coeff.at(0).at(1) = 0.05;                  // component 2

    // Sec 2 ...

    // External
    double              ext_velocity            = 1e-4;    // [m/s]

    // Read external profile from data file
    std::vector<double> ext_profile;
    std::vector<double> ext_prof_delta;
    std::string filename("temp_profile.csv");
    std::string line;
    std::ifstream file(filename.c_str());
    if (file.is_open())
    {
        while (file.good())
        {
            std::getline(file, line);
            if (line.size() > 1) // discard empty lines
            {
                std::replace(line.begin(), line.end(), ',', '.');   // interpret comma or point as decimal point.
                std::size_t pos = line.find(';');                   // find semicolon separator

                ext_profile.push_back(atof(line.substr(pos+1, line.npos).c_str())); // store value
                ext_prof_delta.push_back(atof(line.substr(0, pos).c_str()));        // store delta
            }
        }
        file.close();
    }
    else
    {
        log::emit<Error>() << "Unable to open file: " << filename << log::endl;
        return EXIT_FAILURE;
    }

    // Discretization

    int                 ncol                    = 64;
    int                 npar                    = 4;
    ReconstructionType  reconstruction          = WENO_REC;

    std::vector<double> par_disc_vector(npar + 1, 0.0); // Initialize vector
    ParticleDiscType    par_disc_type           = EQUIDISTANT_PAR;
//                        par_disc_vector.at(0)   = 0.0;
//                        par_disc_vector.at(1)   = 0.5;
//                        par_disc_vector.at(2)   = 1.0;

    // WENO

    int                 boundary_model          = 0;
    double              weno_eps                = 1e-6;
    int                 weno_order              = 3;

    // Schur solver

    int                 gs_type                 = 1;
    int                 max_krylov              = 0;
    int                 max_restarts            = 0;
    double              schur_safety            = 1e-8;

    // Time integrator settings

    double              abstol                  = 1e-8;
    double              reltol                  = 0.0;
    double              abstol_sens_ref         = 1e-5;
    double              init_step_size          = 1e-6;
    int                 max_steps               = 100000;

    // Solver settings

    bool                print_progress          = true;
    bool                print_statistics        = true;
    bool                print_timing            = true;
    bool                print_paramlist         = false;
    bool                print_config            = false;
    bool                use_analytic_jacobian   = true;
    bool                write_at_user_times     = true;
    std::string         log_level               = "INFO";
    std::vector<double> user_solution_times;
    for (int t = 0; t <= 3500; t += 5)
        user_solution_times.push_back(t);
    bool                write_solution_times            = true;
    bool                write_solution_column_outlet    = true;
    bool                write_solution_column_inlet     = true;
    bool                write_solution_all              = true;
    bool                write_sens_column_outlet        = true;
    bool                write_sens_all                  = false;

    // Sensitivity

    int                      nsens                   = 0;
    std::string              sens_method             = "ad1";
    std::vector<std::string> sens_name(nsens, "");    // Initialize vector
    std::vector<int>         sens_comp(nsens, 0);     // Initialize vector
    std::vector<int>         sens_section(nsens, 0);  // Initialize vector
    std::vector<double>      sens_abstol(nsens, 0.0); // Initialize vector
    std::vector<double>      sens_fd_delta(nsens, 0.0); // Initialize vector

    // Sensitivity 1
//                             sens_name    .at(0)     = e2s(EXTL_KD);
//                             sens_comp    .at(0)     = 0;
//                             sens_section .at(0)     = -1;
//                             sens_abstol  .at(0)     = (extl_kd.at(0) > 1e-13) ? abstol / extl_kd.at(0) : abstol_sens_ref;
//                             sens_fd_delta.at(0)     = 1e-3;

    // Sensitivity 2
//                             sens_name    .at(1)     = e2s(EXTL_KD_T);
//                             sens_comp    .at(1)     = 0;
//                             sens_section .at(1)     = -1;
//                             sens_abstol  .at(1)     = (extl_kd_T.at(0) > 1e-13) ? abstol / extl_kd_T.at(0) : abstol_sens_ref;
//                             sens_fd_delta.at(1)     = 1e-3;

    // Sensitivity 3
//                             sens_name    .at(2)     = e2s(COL_DISPERSION);
//                             sens_comp    .at(2)     = -1;                   // does not apply
//                             sens_section .at(2)     = -1;                   // all sections
//                             sens_abstol  .at(2)     = (col_dispersion > 1e-13) ? abstol / col_dispersion : abstol_sens_ref;
//                             sens_fd_delta.at(2)     = 1e-3;

    // Sensitivity 4
//                             sens_name    .at(3)     = "LIN_COEFF";
//                             sens_comp    .at(3)     = 0;
//                             sens_section .at(3)     = 0;
//                             sens_abstol  .at(3)     = (lin_coeff.at(0).at(0) > 1e-13) ? abstol / lin_coeff.at(0).at(0) : abstol_sens_ref;
//                             sens_fd_delta.at(3)     = 1e-3;



    // Writing the actual HDF5 file
    // -----------------------------------------------------------------------------------------------------------------------

    std::ostringstream oss;

    try {
        HDF5Writer h5w;
        h5w.openFile(fileName, "co");
        h5w.extendibleFields(true);  // Make the created vectors extendible, i.e. unlimited maximum size

        h5w.setGroup(e2s(GRP_IN));

            h5w.scalar<std::string> (e2s(CHROMATOGRAPHY_TYPE), e2s(chromatography_type));

        h5w.setGroup(e2s(GRP_IN_MODEL));

            h5w.scalar<std::string> (e2s(ADSORPTION_TYPE),     e2s(adsorption_type));
            h5w.scalar<int>         (e2s(NCOMP),               ncomp);
            h5w.scalar<double>      (e2s(COL_DISPERSION),      col_dispersion);
            h5w.scalar<double>      (e2s(COL_LENGTH),          col_length);
            h5w.scalar<double>      (e2s(COL_POROSITY),        col_porosity);
            h5w.vector<double>      (e2s(FILM_DIFFUSION),      film_diffusion);
            h5w.vector<double>      (e2s(PAR_DIFFUSION),       par_diffusion);
            h5w.scalar<double>      (e2s(PAR_POROSITY),        par_porosity);
            h5w.scalar<double>      (e2s(PAR_RADIUS),          par_radius);
            h5w.vector<double>      (e2s(PAR_SURFDIFFUSION),   par_surfdiffusion);
            h5w.scalar<double>      (e2s(VELOCITY),            velocity);

        h5w.setGroup(e2s(GRP_IN_ADSORPTION));

            h5w.scalar<int>         (e2s(IS_KINETIC),          is_kinetic);
            h5w.vector<double>      (e2s(EXTL_KA),             extl_ka);
            h5w.vector<double>      (e2s(EXTL_KA_T),           extl_ka_T);
            h5w.vector<double>      (e2s(EXTL_KA_TT),          extl_ka_TT);
            h5w.vector<double>      (e2s(EXTL_KA_TTT),         extl_ka_TTT);
            h5w.vector<double>      (e2s(EXTL_KD),             extl_kd);
            h5w.vector<double>      (e2s(EXTL_KD_T),           extl_kd_T);
            h5w.vector<double>      (e2s(EXTL_KD_TT),          extl_kd_TT);
            h5w.vector<double>      (e2s(EXTL_KD_TTT),         extl_kd_TTT);
            h5w.vector<double>      (e2s(EXTL_QMAX),           extl_qmax);
            h5w.vector<double>      (e2s(EXTL_QMAX_T),         extl_qmax_T);
            h5w.vector<double>      (e2s(EXTL_QMAX_TT),        extl_qmax_TT);
            h5w.vector<double>      (e2s(EXTL_QMAX_TTT),       extl_qmax_TTT);

        h5w.setGroup(e2s(GRP_IN_INLET));

            h5w.scalar<int>         (e2s(NSEC),                nsec);
            h5w.vector<double>      (e2s(SECTION_TIMES),       section_times);
            h5w.vector<int>         (e2s(SECTION_CONTINUITY),  section_continuity);

            for (int i = 0; i < nsec; ++i)
            {
                oss.str("");
                oss << e2s(GRP_IN_INLET) << "/sec_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << i;
                h5w.setGroup(oss.str());

                h5w.vector<double>  ("CONST_COEFF",             const_coeff.at(i));
                h5w.vector<double>  ("LIN_COEFF",               lin_coeff.at(i));
                h5w.vector<double>  ("QUAD_COEFF",              quad_coeff.at(i));
                h5w.vector<double>  ("CUBE_COEFF",              cube_coeff.at(i));
            }

        h5w.setGroup(e2s(GRP_IN_EXTERNAL));

            h5w.scalar<double>      (e2s(EXT_VELOCITY),        ext_velocity);
            h5w.vector<double>      (e2s(EXT_PROFILE),         ext_profile);
            h5w.vector<double>      (e2s(EXT_PROF_DELTA),      ext_prof_delta);

        h5w.setGroup(e2s(GRP_IN_DISCRETIZATION));

            h5w.scalar<int>         (e2s(NCOL),                ncol);
            h5w.scalar<int>         (e2s(NPAR),                npar);
            h5w.scalar<std::string> (e2s(PAR_DISC_TYPE),       e2s(par_disc_type));
            h5w.scalar<std::string> (e2s(RECONSTRUCTION),      e2s(reconstruction));

        h5w.setGroup(e2s(GRP_IN_WENO));

            h5w.scalar<int>         (e2s(BOUNDARY_MODEL),      boundary_model);
            h5w.scalar<double>      (e2s(WENO_EPS),            weno_eps);
            h5w.scalar<int>         (e2s(WENO_ORDER),          weno_order);
            h5w.scalar<int>         (e2s(WENO_ORDER),          weno_order);

        h5w.setGroup(e2s(GRP_IN_SOLVER));

            h5w.scalar<int>         (e2s(PRINT_PROGRESS),      print_progress);
            h5w.scalar<int>         (e2s(PRINT_STATISTICS),    print_statistics);
            h5w.scalar<int>         (e2s(PRINT_TIMING),        print_timing);
            h5w.scalar<int>         (e2s(PRINT_PARAMLIST),     print_paramlist);
            h5w.scalar<int>         (e2s(PRINT_CONFIG),        print_config);
            h5w.scalar<int>         (e2s(USE_ANALYTIC_JACOBIAN), use_analytic_jacobian);
            h5w.scalar<int>         (e2s(WRITE_AT_USER_TIMES), write_at_user_times);
            h5w.scalar<std::string> (e2s(LOG_LEVEL),           log_level);
            h5w.vector<double>      (e2s(USER_SOLUTION_TIMES), user_solution_times);
            h5w.scalar<int>         (e2s(WRITE_SOLUTION_TIMES),             write_solution_times);
            h5w.scalar<int>         (e2s(WRITE_SOLUTION_COLUMN_OUTLET),     write_solution_column_outlet);
            h5w.scalar<int>         (e2s(WRITE_SOLUTION_COLUMN_INLET),      write_solution_column_inlet);
            h5w.scalar<int>         (e2s(WRITE_SOLUTION_ALL),               write_solution_all);
            h5w.scalar<int>         (e2s(WRITE_SENS_COLUMN_OUTLET),         write_sens_column_outlet);
            h5w.scalar<int>         (e2s(WRITE_SENS_ALL),                   write_sens_all);

        h5w.setGroup(e2s(GRP_IN_SCHUR));

            h5w.scalar<int>         (e2s(GS_TYPE),             gs_type);
            h5w.scalar<int>         (e2s(MAX_KRYLOV),          max_krylov);
            h5w.scalar<int>         (e2s(MAX_RESTARTS),        max_restarts);
            h5w.scalar<double>      (e2s(SCHUR_SAFETY),        schur_safety);

        h5w.setGroup(e2s(GRP_IN_TIME));

            h5w.scalar<double>      (e2s(ABSTOL),              abstol);
            h5w.scalar<double>      (e2s(RELTOL),              reltol);
            h5w.vector<double>      (e2s(INIT_C),              init_c);
            h5w.vector<double>      (e2s(INIT_Q),              init_q);
            h5w.scalar<double>      (e2s(INIT_STEP_SIZE),      init_step_size);
            h5w.scalar<int>         (e2s(MAX_STEPS),           max_steps);

        h5w.setGroup(e2s(GRP_IN_SENSITIVITY));

            h5w.scalar<int>         (e2s(NSENS),                nsens);
            h5w.scalar<std::string> (e2s(SENS_METHOD),          sens_method);
            for (int i = 0; i < nsens; ++i)
            {
                oss.str("");
                oss << e2s(GRP_IN_SENSITIVITY) << "/param_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << i;
                h5w.setGroup(oss.str());

                h5w.scalar<std::string> (e2s(SENS_NAME),        sens_name.at(i));
                h5w.scalar<int>         (e2s(SENS_COMP),        sens_comp.at(i));
                h5w.scalar<int>         (e2s(SENS_SECTION),     sens_section.at(i));
                h5w.scalar<double>      (e2s(SENS_ABSTOL),      sens_abstol.at(i));
                h5w.scalar<double>      (e2s(SENS_FD_DELTA),    sens_fd_delta.at(i));
            }

        h5w.closeFile();
    } 
    catch (H5::Exception& e)
    {
        log::emit<Except>() << e.getDetailMsg() << log::endl;
    }

    return 0;
}
