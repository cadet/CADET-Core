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

#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>

#include "CadetEnumeration.hpp"
#include "hdf5/HDF5Writer.hpp"
#include "CadetLogger.hpp"

using namespace cadet;


void printHelp()
{
    log::emit() << "Usage: createRegTestC14 <FILE>" << log::endl;
    log::emit() << "Create an HDF5 input file for a benchmark case using 14 components." << log::endl;
    log::emit() << "To be used with cadet-cs from the Cromatography Analysis and DEsign Toolkit (CADET)." << log::endl;
    log::emit() << "When called with no arguments, default file name is regC14.h5." << log::endl;
    log::emit() << "Example: createC14 newfile.h5" << log::endl;
    log::emit() << log::endl;
    log::emit() << "Report bugs to: cadet@fz-juelich.de" << log::endl;
    log::emit() << "CADET homepage: <http://www.cadet-web.de>" << log::endl;
    log::emit() << "Fork CADET on GitHub: <https://github.com/modsim/CADET>" << log::endl;
}


int main(int argc, char** argv)
{
    std::string fileName("regC14.h5");  // default output filename
    if (argc == 2) fileName = argv[1];
    if (argc > 2)
    {
        printHelp();
        return 1;
    }

    ChromatographyType  chromatography_type     = GENERAL_RATE_MODEL;

    // Model

    AdsorptionType      adsorption_type         = MOBILE_PHASE_MODULATORS;
    int                 ncomp                   = 14;
    std::vector<double> film_diffusion(ncomp, 1.0);     // Initialize vector
    std::vector<double> par_diffusion(ncomp, 1e-5);     // Initialize vector
    std::vector<double> par_surfdiffusion(ncomp, 0.0);  // Initialize vector

    double              col_dispersion          = 2.0e-9;
    double              col_length              = 0.15;
    double              col_porosity            = 0.4;

    double              par_porosity            = 0.6;
    double              par_radius              = 2.5e-6;

    double              velocity                = 4.94e-4;
    std::vector<double> init_c(ncomp, 0.0);     // Initialize vector
    std::vector<double> init_q(ncomp, 0.0);     // Initialize vector

                        init_c.at(0)            = 1.0e+2;     // Initial concentration in c for comp 1
//                        init_c.at(1)            = 0.0;      // Initial concentration in c for comp 2
//                        init_c.at(2)            = 0.0;      // Initial concentration in c for comp 3
//                        init_c.at(3)            = 0.0;      // Initial concentration in c for comp 4
//
                        init_q.at(0)            = 160;   // Initial concentration in q for comp 1
//                        init_q.at(1)            = 0.0;      // Initial concentration in q for comp 2
//                        init_q.at(2)            = 0.0;      // Initial concentration in q for comp 3
//                        init_q.at(3)            = 0.0;      // Initial concentration in q for comp 4

    // Adsorption

    std::vector<double> mpm_ka(ncomp,    0.0);     // Initialize vector
    std::vector<double> mpm_kd(ncomp,    1.76e-3); // Initialize vector
    std::vector<double> mpm_qmax(ncomp,  160.0);   // Initialize vector
    std::vector<double> mpm_beta(ncomp,  4.1);     // Initialize vector
    std::vector<double> mpm_gamma(ncomp, 0.0);     // Initialize vector

    int                 is_kinetic              = 1;

                        mpm_ka.at(0)            = 0.0;      // comp 1
                        mpm_ka.at(1)            = 2.73e4;   // comp 2
                        mpm_ka.at(2)            = 1.16e3;   // comp 3
                        mpm_ka.at(3)            = 8.28e3;   // comp 4
                        mpm_ka.at(4)            = 1.13e5;   // comp 5
                        mpm_ka.at(5)            = 8.41e2;   // comp 6
                        mpm_ka.at(6)            = 4.35e5;   // comp 7
                        mpm_ka.at(7)            = 2.62e2;   // comp 8
                        mpm_ka.at(8)            = 1.69e2;   // comp 9
                        mpm_ka.at(9)            = 5.35e4;   // comp 10
                        mpm_ka.at(10)           = 5.62e2;   // comp 11
                        mpm_ka.at(11)           = 2.07e2;   // comp 12
                        mpm_ka.at(12)           = 1.72e4;   // comp 13
                        mpm_ka.at(13)           = 2.71e5;   // comp 14

                        mpm_qmax.at(0)          = 0.0;      // comp 1

                        mpm_beta.at(0)          = 0.0;      // comp 1


    // Inlet

    int                 nsec                    = 2;
    std::vector<double> section_times(nsec+1, 0.0);   // Initialize vector

                        section_times.at(0)     = 0.0;
                        section_times.at(1)     = 1.5231e+001;
                        section_times.at(2)     = 9.1236e+003;

    std::vector<int> section_continuity(nsec-1, 0);

    std::vector<std::vector<double> > const_coeff(nsec, std::vector<double>(ncomp, 0.0));      // Initialize vector
    std::vector<std::vector<double> > lin_coeff  (nsec, std::vector<double>(ncomp, 0.0));      // Initialize vector
    std::vector<std::vector<double> > quad_coeff (nsec, std::vector<double>(ncomp, 0.0));      // Initialize vector
    std::vector<std::vector<double> > cube_coeff (nsec, std::vector<double>(ncomp, 0.0));      // Initialize vector

    // Only values different from 0.0 must be touched!
    // Sec 1
                        const_coeff.at(0).at(0)  = 1.0e+2;      // component 1
                        const_coeff.at(0).at(1)  = 1.0e-2;      // component 2
                        const_coeff.at(0).at(2)  = 1.0e-2;      // component 3
                        const_coeff.at(0).at(3)  = 1.0e-2;      // component 4
                        const_coeff.at(0).at(4)  = 1.0e-2;      // component 5
                        const_coeff.at(0).at(5)  = 1.0e-2;      // component 6
                        const_coeff.at(0).at(6)  = 1.0e-2;      // component 7
                        const_coeff.at(0).at(7)  = 1.0e-2;      // component 8
                        const_coeff.at(0).at(8)  = 1.0e-2;      // component 9
                        const_coeff.at(0).at(9)  = 1.0e-2;      // component 10
                        const_coeff.at(0).at(10) = 1.0e-2;      // component 11
                        const_coeff.at(0).at(11) = 1.0e-2;      // component 12
                        const_coeff.at(0).at(12) = 1.0e-2;      // component 13
                        const_coeff.at(0).at(13) = 1.0e-2;      // component 14

    // Sec 2
//                        const_coeff.at(1).at(0) = 1.0e+2;       // component 1
                        const_coeff.at(1).at(0) = 1.0e+2 + 1.3153e-2 * 1.5231e+001;       // component 1

                        lin_coeff.at(1).at(0)   = 1.3153e-2;    // component 1


    // Discretization

    int                 ncol                    = 128;
    int                 npar                    = 3;
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

    // Solver

    bool                print_progress          = true;
    bool                print_statistics        = true;
    bool                print_timing            = true;
    bool                print_paramlist         = false;
    bool                print_config            = true;
    bool                use_analytic_jacobian   = true;
    bool                write_at_user_times     = true;
    std::string         log_level               = "INFO";
    std::vector<double> user_solution_times;

    for (int t = 0; t <= 15; t += 1)
        user_solution_times.push_back(t);

    user_solution_times.push_back(15.2);
    for (int t = 1; t <= 100; t += 1)
        user_solution_times.push_back(15.2 + 0.001 * static_cast<double>(t));

    user_solution_times.push_back(15.4);
    user_solution_times.push_back(15.5);
    user_solution_times.push_back(15.75);

    for (int t = 16; t <= 9123; t += 1)
        user_solution_times.push_back(t);

    bool                write_solution_times            = true;
    bool                write_solution_column_outlet    = true;
    bool                write_solution_column_inlet     = true;
    bool                write_solution_all              = false;
    bool                write_sens_column_outlet        = true;
    bool                write_sens_all                  = false;

    // Schur solver

    int                 gs_type                 = 1;
    int                 max_krylov              = 0;
    int                 max_restarts            = 0;
    double              schur_safety            = 1e-8;

    // Time integrator

    double              abstol                  = 1e-08;
    double              reltol                  = 0.0;
    double              abstol_sens_ref         = 1e-7;

    double              init_step_size          = 1e-6;
    int                 max_steps               = 10000;

    // Sensitivity

    int                      nsens                   = 4;
    std::string              sens_method             = "ad1";
    std::vector<std::string> sens_name(nsens, "");    // Initialize vector
    std::vector<int>         sens_comp(nsens, 0);     // Initialize vector
    std::vector<int>         sens_section(nsens, 0);  // Initialize vector
    std::vector<double>      sens_abstol(nsens, 0.0); // Initialize vector
    std::vector<double>      sens_fd_delta(nsens, 0.0); // Initialize vector

    // Sensitivity 1
                             sens_name    .at(0)     = e2s(COL_DISPERSION);
                             sens_comp    .at(0)     = -1;                   // does not apply
                             sens_section .at(0)     = -1;                   // all sections
                             sens_abstol  .at(0)     = (col_dispersion > 1e-13) ? abstol / col_dispersion : abstol_sens_ref;
                             sens_fd_delta.at(0)     = 1e-3;

    // Sensitivity 2
                             sens_name    .at(1)     = e2s(MPM_KA);
                             sens_comp    .at(1)     = 1;
                             sens_section .at(1)     = -1;
                             sens_abstol  .at(1)     = (mpm_ka.at(1) > 1e-13) ? abstol / mpm_ka.at(1) : abstol_sens_ref;
                             sens_fd_delta.at(1)     = 1e-3;

    // Sensitivity 3
                             sens_name    .at(2)     = "CONST_COEFF";
                             sens_comp    .at(2)     = 0;
                             sens_section .at(2)     = 0;
                             sens_abstol  .at(2)     = (const_coeff.at(0).at(0) > 1e-13) ? abstol / const_coeff.at(0).at(0) : abstol_sens_ref;
                             sens_fd_delta.at(2)     = 1e-3;

    // Sensitivity 4
                             sens_name    .at(3)     = "LIN_COEFF";
                             sens_comp    .at(3)     = 0;
                             sens_section .at(3)     = 0;
                             sens_abstol  .at(3)     = (lin_coeff.at(0).at(0) > 1e-13) ? abstol / lin_coeff.at(0).at(0) : abstol_sens_ref;
                             sens_fd_delta.at(3)     = 1e-3;



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
            h5w.scalar<int>         (e2s(NCOMP),               ncomp);
            h5w.vector<double>      (e2s(INIT_C),              init_c);
            h5w.vector<double>      (e2s(INIT_Q),              init_q);
            h5w.scalar<std::string> (e2s(ADSORPTION_TYPE),     e2s(adsorption_type));
            h5w.scalar<double>      (e2s(COL_DISPERSION),      col_dispersion);
            h5w.scalar<double>      (e2s(COL_LENGTH),          col_length);
            h5w.scalar<double>      (e2s(COL_POROSITY),        col_porosity);
            h5w.scalar<double>      (e2s(PAR_POROSITY),        par_porosity);
            h5w.scalar<double>      (e2s(PAR_RADIUS),          par_radius);
            h5w.scalar<double>      (e2s(VELOCITY),            velocity);
            h5w.vector<double>      (e2s(FILM_DIFFUSION),      film_diffusion);
            h5w.vector<double>      (e2s(PAR_DIFFUSION),       par_diffusion);
            h5w.vector<double>      (e2s(PAR_SURFDIFFUSION),   par_surfdiffusion);
        
        h5w.setGroup(e2s(GRP_IN_ADSORPTION));

            h5w.scalar<int>         (e2s(IS_KINETIC),          is_kinetic);
            h5w.vector<double>      (e2s(MPM_KA),              mpm_ka);
            h5w.vector<double>      (e2s(MPM_KD),              mpm_kd);
            h5w.vector<double>      (e2s(MPM_QMAX),            mpm_qmax);
            h5w.vector<double>      (e2s(MPM_BETA),            mpm_beta);
            h5w.vector<double>      (e2s(MPM_GAMMA),           mpm_gamma);

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

        h5w.setGroup(e2s(GRP_IN_DISCRETIZATION));

            h5w.scalar<int>         (e2s(NCOL),                ncol);
            h5w.scalar<int>         (e2s(NPAR),                npar);
            h5w.scalar<std::string> (e2s(PAR_DISC_TYPE),       e2s(par_disc_type));
            h5w.vector<double>      (e2s(PAR_DISC_VECTOR),     par_disc_vector);
            h5w.scalar<std::string> (e2s(RECONSTRUCTION),      e2s(reconstruction));

        h5w.setGroup(e2s(GRP_IN_WENO));

            h5w.scalar<int>         (e2s(BOUNDARY_MODEL),      boundary_model);
            h5w.scalar<double>      (e2s(WENO_EPS),            weno_eps);
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
            h5w.scalar<int>         (e2s(WRITE_SOLUTION_TIMES),             write_solution_times);
            h5w.scalar<int>         (e2s(WRITE_SOLUTION_COLUMN_OUTLET),     write_solution_column_outlet);
            h5w.scalar<int>         (e2s(WRITE_SOLUTION_COLUMN_INLET),      write_solution_column_inlet);
            h5w.scalar<int>         (e2s(WRITE_SOLUTION_ALL),               write_solution_all);
            h5w.scalar<int>         (e2s(WRITE_SENS_COLUMN_OUTLET),         write_sens_column_outlet);
            h5w.scalar<int>         (e2s(WRITE_SENS_ALL),                   write_sens_all);
            h5w.vector<double>      (e2s(USER_SOLUTION_TIMES), user_solution_times);

        h5w.setGroup(e2s(GRP_IN_SCHUR));

            h5w.scalar<int>         (e2s(GS_TYPE),             gs_type);
            h5w.scalar<int>         (e2s(MAX_KRYLOV),          max_krylov);
            h5w.scalar<int>         (e2s(MAX_RESTARTS),        max_restarts);
            h5w.scalar<double>      (e2s(SCHUR_SAFETY),        schur_safety);

        h5w.setGroup(e2s(GRP_IN_TIME));

            h5w.scalar<double>      (e2s(ABSTOL),              abstol);
            h5w.scalar<double>      (e2s(RELTOL),              reltol);
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
