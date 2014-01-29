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

#include "Simulator.hpp"
#include "hdf5/HDF5Writer.hpp"
#include "CadetLogger.hpp"

using namespace cadet;


void printHelp()
{
    log::emit() << "Usage: createSCL <FILE>" << log::endl;
    log::emit() << "Create an HDF5 input file for single-component-langmuir benchmark case." << log::endl;
    log::emit() << "To be used with cadet-cs from the Cromatography Analysis and DEsign Toolkit (CADET)." << log::endl;
    log::emit() << "When called with no arguments, default file name is SCL.h5." << log::endl;
    log::emit() << "Example: createSCL newfile.h5" << log::endl;
    log::emit() << log::endl;
    log::emit() << "Report bugs to: cadet@fz-juelich.de" << log::endl;
    log::emit() << "CADET homepage: <http://www.cadet-web.de>" << log::endl;
    log::emit() << "Fork CADET on GitHub: <https://github.com/modsim/CADET>" << log::endl;
}


int main(int argc, char** argv)
{
    std::string fileName("SCL.h5");  // default output filename
    if (argc == 2) fileName = argv[1];
    if (argc > 2)
    {
        printHelp();
        return 1;
    }

    ChromatographyType  chromatography_type     = GENERAL_RATE_MODEL; // check

    // Model

    AdsorptionType      adsorption_type         = MULTI_COMPONENT_LANGMUIR; // check
    int                 ncomp                   = 1; // check
    double              col_dispersion          = 5.75e-8; // check
    double              col_length              = 0.014; // check
    double              col_porosity            = 0.37; // check
    std::vector<double> film_diffusion(ncomp, 0.0);   // Initialize vector
                        film_diffusion.at(0)    = 6.9e-6; // check
    std::vector<double> par_diffusion(ncomp, 0.0);   // Initialize vector
                        par_diffusion.at(0)     = 6.07e-11; // check
    double              par_porosity            = 0.75; // check
    double              par_radius              = 4.5e-5; // check
    std::vector<double> par_surfdiffusion(ncomp, 0.0);   // Initialize vector
                        par_surfdiffusion.at(0) = 0.0; // check
    double              velocity                = 5.75e-4; // check
    
    std::vector<double> init_c(ncomp, 0.0);     // Initialize vector
    std::vector<double> init_q(ncomp, 0.0);     // Initialize vector
                        init_c.at(0)            = 0.0; // check
                        init_q.at(0)            = 0.0; // check

    // Adsorption

    std::vector<double> mcl_ka(ncomp, 0.0);     // Initialize vector
    std::vector<double> mcl_kd(ncomp, 0.0);     // Initialize vector
    std::vector<double> mcl_qmax(ncomp, 0.0);   // Initialize vector
    int                 is_kinetic              = 1; // check
                        mcl_ka.at(0)            = 1.14; // check
                        mcl_kd.at(0)            = 0.002; // check
                        mcl_qmax.at(0)          = 4.88; // check

    // Inlet

    int                 nsec                    = 1; // check
    std::vector<double> section_times(nsec+1, 0.0);   // Initialize vector
                        section_times.at(0)     = 0.0; // check
                        section_times.at(1)     = 10000.0; // check

    std::vector<std::vector<double> > const_coeff(nsec, std::vector<double>(ncomp, 0.0));      // Initialize vector
    std::vector<std::vector<double> > lin_coeff  (nsec, std::vector<double>(ncomp, 0.0));      // Initialize vector
    std::vector<std::vector<double> > quad_coeff (nsec, std::vector<double>(ncomp, 0.0));      // Initialize vector
    std::vector<std::vector<double> > cube_coeff (nsec, std::vector<double>(ncomp, 0.0));      // Initialize vector

    // Sec 1
                        const_coeff.at(0).at(0) = 7.14e-3; // check                  // component 1
                        lin_coeff  .at(0).at(0) = 0.0;
                        quad_coeff .at(0).at(0) = 0.0;
                        cube_coeff .at(0).at(0) = 0.0;

    // Sec 2 ...

    // Discretization

    ReconstructionType  reconstruction          = WENO_REC;
    int                 ncol                    = 64; // check
    int                 npar                    = 16; // check
    std::vector<double> par_disc_vector(npar + 1, 0.0); // Initialize vector
    ParticleDiscType    par_disc_type           = EQUIDISTANT_PAR; // check
//                        par_disc_vector.at(0)   = 0.0;
//                        par_disc_vector.at(1)   = 0.5;
//                        par_disc_vector.at(2)   = 1.0;

    // WENO

    int                 boundary_model          = 0; // check
    double              weno_eps                = 1e-6; // check
    int                 weno_order              = 3; // check

    // Solver

    bool                print_progress          = false; // check
    bool                print_statistics        = true; // check
    bool                print_timing            = true; // check
    bool                print_paramlist         = false; // check
    bool                print_config            = true; // check
    bool                use_analytic_jacobian   = true; // check
    bool                write_at_user_times     = true; // check
    std::string         log_level               = "INFO"; // check
    std::vector<double> user_solution_times;
    for (int t = 0; t <= 10000; t += 10)
        user_solution_times.push_back(t);
    bool                write_solution_times            = true; // check
    bool                write_solution_column_outlet    = true; // check
    bool                write_solution_column_inlet     = true; // check
    bool                write_solution_all              = true; // check
    bool                write_sens_column_outlet        = true; // check
    bool                write_sens_all                  = true; // check


    // Schur solver

    int                 gs_type                 = 1; // check
    int                 max_krylov              = 0; // check
    int                 max_restarts            = 0; // check
    double              schur_safety            = 1e-8; // check

    // Time integrator

    double              abstol                  = 1e-8; // check
    double              reltol                  = 0.0; // check
    double              abstol_sens_ref         = 1e-5; // check
    double              init_step_size          = 1e-6; // check
    int                 max_steps               = 10000; // check

    // Sensitivity

    int                      nsens                   = 0;
    std::string              sens_method             = "ad1";
    std::vector<std::string> sens_name(nsens, "");    // Initialize vector
    std::vector<int>         sens_comp(nsens, 0);     // Initialize vector
    std::vector<int>         sens_section(nsens, 0);  // Initialize vector
    std::vector<double>      sens_abstol(nsens, 0.0); // Initialize vector
    std::vector<double>      sens_fd_delta(nsens, 0.0); // Initialize vector

    // Sensitivity 1
//                             sens_name    .at(0)     = e2s(COL_DISPERSION);
//                             sens_comp    .at(0)     = -1;                   // does not apply
//                             sens_section .at(0)     = -1;                   // all sections
//                             sens_abstol  .at(0)     = (col_dispersion > 1e-13) ? abstol / col_dispersion : abstol_sens_ref;
//                             sens_fd_delta.at(0)     = 1e-3;

    // Sensitivity 2
//                             sens_name    .at(1)     = e2s(MCL_KA);
//                             sens_comp    .at(1)     = 0;
//                             sens_section .at(1)     = -1;
//                             sens_abstol  .at(1)     = (mcl_ka.at(0) > 1e-13) ? abstol / mcl_ka.at(0) : abstol_sens_ref;
//                             sens_fd_delta.at(1)     = 1e-3;

    // Sensitivity 3
//                             sens_name    .at(2)     = "CONST_COEFF";
//                             sens_comp    .at(2)     = 0;
//                             sens_section .at(2)     = 0;
//                             sens_abstol  .at(2)     = (const_coeff.at(0).at(0) > 1e-13) ? abstol / const_coeff.at(0).at(0) : abstol_sens_ref;
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
            h5w.vector<double>      (e2s(INIT_C),              init_c);
            h5w.vector<double>      (e2s(INIT_Q),              init_q);

        h5w.setGroup(e2s(GRP_IN_ADSORPTION));

            h5w.scalar<int>         (e2s(IS_KINETIC),          is_kinetic);
            h5w.vector<double>      (e2s(MCL_KA),              mcl_ka);
            h5w.vector<double>      (e2s(MCL_KD),              mcl_kd);
            h5w.vector<double>      (e2s(MCL_QMAX),            mcl_qmax);

        h5w.setGroup(e2s(GRP_IN_INLET));

            h5w.scalar<int>         (e2s(NSEC),                nsec);
            h5w.vector<double>      (e2s(SECTION_TIMES),       section_times);

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
