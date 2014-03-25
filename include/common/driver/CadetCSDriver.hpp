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

#include "CadetLogger.hpp"
#include "Cadet.hpp"

#ifdef _OPENMP
    #include <omp.h>
#endif

namespace cadet {

class inletParams : public InletBase
{
public:
    inletParams(int ncomp, int nsec) :
        const_coeff(nsec, std::vector<double>(ncomp, 0.0)),
        lin_coeff  (nsec, std::vector<double>(ncomp, 0.0)),
        quad_coeff (nsec, std::vector<double>(ncomp, 0.0)),
        cube_coeff (nsec, std::vector<double>(ncomp, 0.0)),
        _ncomp(ncomp),
        _nsec(nsec)
    {}

    void inletConcentration(double t, int sec,
            std::vector<double>& inletConc,
            const std::vector<bool>& inletParamIsSensitive,
            std::vector<std::vector<double> >& dInletConc_dp)
    {
        // t:                       absolute simulation time
        // sec:                     section for which the boundary condition is needed
        // inletConc:               concentration at the boundary at time t of each componenet (size: NCOMP x 1)
        // inletParamIsSensitive:   blala
        // dInletConc_dp:           derivatives of all c_in w.r.t. all parameters of interest (size: NCOMP x NSENSPAR)

        // This function evaluates a cubic spline, which is a piecewise polynomial given on some intervals
        // called sections. On each section a polynomial of degree 3 is evaluated:
        // 
        //   p_i(t) = CUBE_COEFF[i] * (t - t_i)^3 + QUAD_COEFF[i] * (t - t_i)^2 + LIN_COEFF[i] * (t - t_i) + CONST_COEFF[i],
        //   
        // where p_i is the polynomial on section i given by the interval [t_i, t_{i+1}].

        if (sec < 0)
        {
            // Time is out of all given sections, so return 0
            for (std::size_t comp = 0; comp < inletConc.size(); ++comp)
            {
                inletConc.at(comp) = 0.0;

                // Derivatives w.r.t. each parameter
                dInletConc_dp.at(comp).at(getIndex(CONST_COEFF, comp, sec)) = 0.0;
                dInletConc_dp.at(comp).at(getIndex(LIN_COEFF,   comp, sec)) = 0.0;
                dInletConc_dp.at(comp).at(getIndex(QUAD_COEFF,  comp, sec)) = 0.0;
                dInletConc_dp.at(comp).at(getIndex(CUBE_COEFF,  comp, sec)) = 0.0;
            }
            return;
        }

        t -= sectionTimes.at(sec);

        for (std::size_t comp = 0; comp < inletConc.size(); ++comp)
        {
            // Polynomial of 3rd order
            // Evaluate using Horner scheme
            inletConc.at(comp) = const_coeff.at(sec).at(comp) + t * (lin_coeff.at(sec).at(comp) + t * (quad_coeff.at(sec).at(comp) + t * cube_coeff.at(sec).at(comp)));
            inletConc.at(comp) = std::max(0.0, inletConc.at(comp));

            // Derivatives w.r.t. each parameter
            dInletConc_dp.at(comp).at(getIndex(CONST_COEFF, comp, sec)) = 1;
            dInletConc_dp.at(comp).at(getIndex(LIN_COEFF,   comp, sec)) = t;
            dInletConc_dp.at(comp).at(getIndex(QUAD_COEFF,  comp, sec)) = pow(t,2);
            dInletConc_dp.at(comp).at(getIndex(CUBE_COEFF,  comp, sec)) = pow(t,3);
        }
    }

    int getMaxParams() const { return _ncoeffs * _ncomp * _nsec; }

    int getIndex(std::string coeff, int comp, int sec) const throw (cadet::CadetException)
    {
        int param = 0;
        if     (coeff == CONST_COEFF) param = 0;
        else if(coeff ==   LIN_COEFF) param = 1;
        else if(coeff ==  QUAD_COEFF) param = 2;
        else if(coeff ==  CUBE_COEFF) param = 3;
        else
            throw cadet::CadetException("Unknown coefficient");
        return sec * _ncomp * _ncoeffs + param * _ncomp + comp;
    }

    std::vector<std::vector<double> > const_coeff;
    std::vector<std::vector<double> > lin_coeff;
    std::vector<std::vector<double> > quad_coeff;
    std::vector<std::vector<double> > cube_coeff;
    std::vector<double> sectionTimes;

    static const std::string CONST_COEFF;
    static const std::string   LIN_COEFF;
    static const std::string  QUAD_COEFF;
    static const std::string  CUBE_COEFF;

private:
    int _ncomp;
    int _nsec;
    static const int _ncoeffs;
};

// Initialize static const members of class inletParams
const std::string inletParams::CONST_COEFF = "CONST_COEFF";
const std::string inletParams::LIN_COEFF   = "LIN_COEFF";
const std::string inletParams::QUAD_COEFF  = "QUAD_COEFF";
const std::string inletParams::CUBE_COEFF  = "CUBE_COEFF";

const int inletParams::_ncoeffs = 4;


template <typename reader_t>
class ExternalProfile : public ExternalBase
{
public:
    ExternalProfile() {}

    ExternalProfile(const std::string& filename)
    {
        reader_t reader;
        read(reader, filename, true);
    }

    ExternalProfile(const std::string& filename, reader_t& reader)
    {
        read(reader, filename, false);
    }

    ExternalProfile(reader_t& reader)
    {
        read(reader, "", false);
    }

    void externalProfile(double zc0, double t, double* value)
    {
        // zc0 denotes relative position in the column, i.e. zc0 is in [0,1]
        // 
        // The coordinate system of the external profile begins at the column outlet
        // and points backward to the column inlet
        //
        //         _lc                                                        0
        //      <---|---------------------------------------------------------|
        //           _________________________________________________________
        //          |                                                         |
        //   Inlet  |  =>               Column                            =>  | Outlet
        //          |_________________________________________________________|
        //                                    
        //                                                                                                  
        //          Profile moves to the right with velocity _vT
        //          Thus, after some time, the coordinate system has shifted
        //
        //     _lc + _vtT * t                                             _vtT * t          0
        //      <---|---------------------------------------------------------|-------------|
        //           _________________________________________________________
        //          |                                                         |
        //   Inlet  |  =>               Column                            =>  | Outlet
        //          |_________________________________________________________|
        //                                                                   

        double zc = zc0 * _lc;              // compute non-reduced z position in the column
        double zT = (_lc - zc) + _vT * t;   // z coordinate for ext-profile at time t and position zc

        double zT0 = zT / _lT;              // reduced z coordinate for ext-profile

        profile1(zT0, value);
    }

    void profile1(double zT0, double* value)
    {
        // Organization of the external profile
        // 
        // 
        //  0         delta[1]               delta[1] + delta[2]
        //  |-----------|-----------------------------|-------------|-----------------|----------|
        //   \         / \                           /
        //     delta[1]    -------  delta[2] -------
        //  
        // In order for this layout to work, first position has to be 0, i.e. _tempDelta[0] = 0 is required
        // 
        // The index of the cell containing zT0 is computed by looping through the data until
        // 
        //           delta[1] + ... + delta[idx+1]
        //   zT0 >  -------------------------------
        //                   Length Profile
        //
        // The requested position zT0 is located in the interval 
        //      [ delta[1] + ... + delta[idx], delta[1] + ... + delta[idx] + delta[idx+1] ]
        // with width delta[idx+1].
        // 
        // The column moves over the profile with velocity _vtT
        // 
        // Outlet             Inlet
        //  |------ Column ------|  -> _vT
        // 
        //  |-----------|-----------------------------|-------------|-----------------|----------|
        //  0         delta[1]               delta[1] + delta[2]

        // Use constant extrapolation on both sides of the external profile
        if      (zT0 <= 0) *value = _tempProfile.front();
        else if (zT0 >= 1) *value = _tempProfile.back();

        // In between read external profile, use linear interpolation
        else
        {
            // Find the corresponding index for the range of the external vector we are working on.
            std::size_t idx;
            double  sum = 0;
            for (idx = 0; idx < _tempDelta.size() - 1; ++idx)
            {
                sum += _tempDelta.at(idx + 1) / _lT;
                if (sum > zT0) break;
            }

            // Now idx is the index of the left and idx + 1 is the index of the right data point

            double dz  = _tempDelta.at(idx + 1) / _lT;  // interval width

            double Ti1 = _tempProfile.at(idx);      // value at left data point
            double Ti2 = _tempProfile.at(idx + 1);  // value at right data point

            double z1 = sum - dz; // left data point

            // linear interpolation
            *value = (Ti2 - Ti1) / dz * (zT0 - z1) + Ti1;
        }
    }

private:
    double              _vT;           //!< Velocity of the movement of the external profile in [m/s]
    double              _lc;           //!< Length of the column in [m]
    double              _lT;           //!< Length of external profile in [m]
    std::vector<double> _tempProfile;  //!< External profile
    std::vector<double> _tempDelta;    //!< Delta between two T measurements in [m]

    void read(reader_t& reader, const std::string& filename, bool openAndClose)
    {
        if (openAndClose)
            reader.openFile(filename);

        reader.setGroup(e2s(GRP_IN_MODEL));
        // Read column length to rescale z coordinate.
        _lc = reader.template scalar<double>(e2s(COL_LENGTH));

        reader.setGroup(e2s(GRP_IN_EXTERNAL));

        // Read the external profile and corresponding deltas
        _tempProfile = reader.template vector<double>(e2s(EXT_PROFILE));
        _tempDelta = reader.template vector<double>(e2s(EXT_PROF_DELTA));
        // Sum up to get external profile length
        _lT = 0;
        for (std::size_t i = 0; i < _tempDelta.size(); ++i)
            _lT += _tempDelta.at(i);

        // Velocity is applied to the profile in flow direction
        _vT = reader.template scalar<double>(e2s(EXT_VELOCITY));

        if (openAndClose)
            reader.closeFile();
    }
};



enum SensMethod
{
    ALGORITHMIC_DIFFERENTIATION_1,
    FINITE_DIFFERENCES_1,
    FINITE_DIFFERENCES_2,
    FINITE_DIFFERENCES_4
};


template <typename reader_t, typename writer_t>
class CadetCS
{
public:
    CadetCS(const std::string& filename);
    CadetCS(const std::string& filename, reader_t& reader, writer_t& writer);

    ~CadetCS();

    void run();

    std::size_t getNSensParams() const { return _nsens; }

    const std::vector<std::string>& getSensModelParamNames() const { return _sensNames; }
    
    const std::vector<double>& getSolutionTimes() const { return _solutionTimes; }

    const std::vector<std::vector<std::vector<double> > >& getSensitivityColumnOutlet() const { return _sensitivityColumnOutlet; }

    const std::vector<std::vector<double> >& getSolutionColumnOutlet() const { return _solutionColumnOutlet; }

private:

    std::string _filename;
    reader_t    _reader;
    writer_t    _writer;

    ChromatographyType  _chromType;
    AdsorptionType      _adsType;
    SensMethod          _sensMethod;

    std::size_t _ncomp;
    std::size_t _ncol;
    std::size_t _npar;
    std::size_t _nsens;
    std::size_t _nsec;
    std::size_t _nsim;
    std::size_t _nAddFdRuns;

    bool _writeSolutionTimes;
    bool _writeSolutionColumnOutlet;
    bool _writeSolutionColumnInlet;
    bool _writeSolutionAll;
    bool _writeSensColumnOutlet;
    bool _writeSensAll;

    std::vector<ParameterName>  _sensParamIds;
    std::vector<int>            _sensComps;
    std::vector<int>            _sensSecs;
    std::vector<std::string>    _sensNames;
    std::vector<double>         _sensFdDeltaAbs;    //!< A vector holding the absolute values of the disturbances used in finite difference sensitivity computation

    std::ostringstream _oss;

    ExternalProfile<reader_t>   _externalProfile;
    std::vector<inletParams*>   _inletParams;
    std::vector<Simulator*>     _sim;

    std::vector<double>                             _solutionTimes;             //!< A vector holding the final solution times (timepoints at which a solution was computed)
//    std::vector<double>                             _allSolutions;              //!< A vector holding the solution at all solution times, for all column and particle cells for all components
    std::vector<double>                             _colSolutions;              //!< A vector holding the solution at all solution times, for all components, at all column cells
    std::vector<double>                             _parSolutions;              //!< A vector holding the solution at all solution times, for all particles, for all particle cells, for all phases, for all components
    std::vector<double>                             _bndSolutions;              //!< A vector holding the solution at all solution times, for all components, at all boundary cells
//    std::vector<double>                             _allSensitivities;          //!< A vector holding the sensitivities for all sensitive parameters, at all solution times, for all column and particle cells, for all components
    std::vector<double>                             _colSensitivities;          //!< A vector holding the sensitivities for all sensitive parameters, at all solution times, for all components, at all column cells
    std::vector<double>                             _parSensitivities;          //!< A vector holding the sensitivities for all sensitive parameters, at all solution times, for all particles, for all particle cells, for all phases, for all components
    std::vector<double>                             _bndSensitivities;          //!< A vector holding the sensitivities for all sensitive parameters, at all solution times, for all components, at all boundary cells
    std::vector<std::vector<double> >               _solutionColumnOutlet;      //!< A vector holding the solutions at the column outlet (chromatograms) for all components
    std::vector<std::vector<double> >               _solutionColumnInlet;       //!< A vector holding the solutions at the column inlet (boundary conditions) for all components
    std::vector<std::vector<std::vector<double> > > _sensitivityColumnOutlet;   //!< A vector holding the sensitivities at the column outlet (sensograms) for all components

    void setupInstance();

    void initialize();

    void configureInlet();

    void setParameters();

    void setSensitivities();

    void simulate();

    void extractSolutionData();

    void calculateSensStuff();

    void writeOutput();
};



template <typename reader_t, typename writer_t>
CadetCS<reader_t, writer_t>::CadetCS(const std::string& filename) :
    _filename(filename)
{
    setupInstance();
}



template <typename reader_t, typename writer_t>
CadetCS<reader_t, writer_t>::CadetCS(const std::string& filename, reader_t& reader, writer_t& writer) :
    _filename(filename), _reader(reader), _writer(writer)
{
    setupInstance();
}



template <typename reader_t, typename writer_t>
CadetCS<reader_t, writer_t>::~CadetCS()
{
    for (std::size_t sim = 0; sim < _nsim; sim++)
    {
        delete _inletParams.at(sim);
        delete _sim.at(sim);
    }
}



template <typename reader_t, typename writer_t>
void CadetCS<reader_t, writer_t>::setupInstance()
{
    _reader.openFile(_filename, "r");

    // ============================================================================================================
    //    Check the method for sensitivity computation
    // ============================================================================================================
    _reader.setGroup(e2s(GRP_IN_SENSITIVITY));
    // Read no. of sensitivites to be computed
    _nsens = _reader.template scalar<int>(e2s(NSENS));

    // Read sensitivity computation method
    std::string sensMethod = _reader.template scalar<std::string>(e2s(SENS_METHOD));

    // Algorithmic differnetiation used for sensitivity computation
    if      (sensMethod == "ad1")
    {
        _sensMethod = ALGORITHMIC_DIFFERENTIATION_1;
        _nAddFdRuns = 0;
    }
    // Finite difference method used for sensitivity computation
    else if (sensMethod == "fd1")
    {
        _sensMethod = FINITE_DIFFERENCES_1;
        _nAddFdRuns = 1;
    }
    else if (sensMethod == "fd2")
    {
        _sensMethod = FINITE_DIFFERENCES_2;
        _nAddFdRuns = 2;
    }
    else if (sensMethod == "fd4")
    {
        _sensMethod = FINITE_DIFFERENCES_4;
        _nAddFdRuns = 4;
    }
    else
        throw CadetException("Wrong sensitivity computation method specified!");

    _nsim = 1 + _nAddFdRuns * _nsens; // We always need one simulator for the standard case
    // ============================================================================================================

    // ============================================================================================================
    //    Create as many instances of simulator as needed
    // ============================================================================================================
    _reader.setGroup(e2s(GRP_IN));
    s2e(_reader.template scalar<std::string>(e2s(CHROMATOGRAPHY_TYPE)), _chromType);
    _reader.setGroup(e2s(GRP_IN_MODEL));
    s2e(_reader.template scalar<std::string>(e2s(ADSORPTION_TYPE)), _adsType);

    _ncomp = _reader.template scalar<int>(e2s(NCOMP));

    _reader.setGroup(e2s(GRP_IN_DISCRETIZATION));
    _ncol  = _reader.template scalar<int>(e2s(NCOL));
    _npar  = _reader.template scalar<int>(e2s(NPAR));

    _reader.setGroup(e2s(GRP_IN_INLET));
    _nsec  = _reader.template scalar<int>(e2s(NSEC));

    for (std::size_t sim = 0; sim < _nsim; ++sim)
        _sim.push_back(new Simulator(_ncomp, _ncol, _npar, _nsec, _adsType, _chromType));
    // ============================================================================================================


    // ============================================================================================================
    //    Create as many instances of inletParams as needed
    // ============================================================================================================
    for (std::size_t sim = 0; sim < _nsim; ++sim)
        _inletParams.push_back(new inletParams(_ncomp, _nsec));
    // ============================================================================================================

    // ============================================================================================================
    //    Create an instances of ExternalProfile
    // ============================================================================================================
    _reader.setGroup(e2s(GRP_IN_MODEL));
    if (_reader.exists("external"))
        _externalProfile = ExternalProfile<reader_t>(_reader);
    // ============================================================================================================
}



template <typename reader_t, typename writer_t>
void CadetCS<reader_t, writer_t>::run()
{
    initialize();
    configureInlet();
    setParameters();
    setSensitivities();

    simulate();
    _reader.closeFile();
    extractSolutionData();
    writeOutput();
}



template <typename reader_t, typename writer_t>
void CadetCS<reader_t, writer_t>::initialize()
{
    // ============================================================================================================
    //    Initialize logging framework
    // ============================================================================================================
    _reader.setGroup(e2s(GRP_IN_SOLVER));
#ifndef LOGGING_DISABLE
    std::string logLevel = _reader.template scalar<std::string>(e2s(LOG_LEVEL));
    if      (logLevel == "ERROR")    { log::emit() -= Level::normal; log::emit() -= Level::warning; }
    else if (logLevel == "WARNING")  { log::emit() -= Level::normal; }
    else if (logLevel == "INFO")     { log::emit() += Level::info; }
    else if (logLevel == "DEBUG1")   { log::emit() += Level::info; log::emit() += Level::debug1; }
    else if (logLevel == "DEBUG2")   { log::emit() += Level::info; log::emit() += Level::debug1; log::emit() += Level::debug2; }
    else if (logLevel == "TRACE1")   { log::emit() += Level::info; log::emit() += Level::debug1; log::emit() += Level::debug2; log::emit() += Level::trace1; }
    else if (logLevel == "TRACE2")   { log::emit() += Level::info; log::emit() += Level::debug1; log::emit() += Level::debug2; log::emit() += Level::trace1; log::emit() += Level::trace2; }
    else                             { log::emit() += Level::info; } // For wrong log level specification assume info level
#endif
    log::emit<Info>() << "This is " << Color::green << cadet::getLibraryVersion() << Color::reset << "!" << log::endl << log::endl;
    // ============================================================================================================


    for (std::vector<Simulator*>::iterator sim = _sim.begin(); sim < _sim.end(); ++sim)
    {
        // ============================================================================================================
        //    Configuring the simulators
        // ============================================================================================================
        _reader.setGroup(e2s(GRP_IN_SCHUR));
        (*sim)->configureSchurSolver
        (
                _reader.template scalar<double>(e2s(SCHUR_SAFETY)),
                _reader.template scalar<int>   (e2s(MAX_RESTARTS)),
                _reader.template scalar<int>   (e2s(MAX_KRYLOV)),
                _reader.template scalar<int>   (e2s(GS_TYPE))
        );

        _reader.setGroup(e2s(GRP_IN_TIME));
        (*sim)->configureTimeIntegrator
        (
                _reader.template scalar<double>(e2s(RELTOL)),
                _reader.template scalar<double>(e2s(ABSTOL)),
                _reader.template scalar<double>(e2s(INIT_STEP_SIZE)),
                _reader.template scalar<int>   (e2s(MAX_STEPS))
        );

        _reader.setGroup(e2s(GRP_IN_WENO));
        (*sim)->configureWenoScheme
        (
                _reader.template scalar<int>   (e2s(WENO_ORDER)),
                _reader.template scalar<double>(e2s(WENO_EPS)),
                _reader.template scalar<int>   (e2s(BOUNDARY_MODEL))
        );

        _reader.setGroup(e2s(GRP_IN_SOLVER));
        (*sim)->configurePrinting
        (
                _reader.template scalar<int>(e2s(PRINT_PROGRESS)),
                _reader.template scalar<int>(e2s(PRINT_STATISTICS)),
                _reader.template scalar<int>(e2s(PRINT_TIMING)),
                _reader.template scalar<int>(e2s(PRINT_PARAMLIST)),
                _reader.template scalar<int>(e2s(PRINT_CONFIG))
        );
        log::emit<Debug1>() << "Simulator configured!" << log::endl;
        // ============================================================================================================


        (*sim)->useAnalyticJacobian(_reader.template scalar<int>(e2s(USE_ANALYTIC_JACOBIAN)));

        // ============================================================================================================
        //    Set solution times vector
        // ============================================================================================================
        if (_reader.template scalar<int>(e2s(WRITE_AT_USER_TIMES)))
            (*sim)->setSolutionTimes(_reader.template vector<double>(e2s(USER_SOLUTION_TIMES)));
        log::emit<Debug1>() << "Solution times vector set!" << log::endl;
        // ============================================================================================================


        // ============================================================================================================
        //    Set the particle discretization scheme
        // ============================================================================================================
        _reader.setGroup(e2s(GRP_IN_DISCRETIZATION));

        if (_reader.template scalar<std::string>(e2s(PAR_DISC_TYPE)).compare(e2s(EQUIVOLUME_PAR)) == 0)
        {
            (*sim)->setPartDiscSchemeEqVol();
            log::emit<Info>() << "Using equivolume particle discretization scheme." << log::endl;
        }
        else if (_reader.template scalar<std::string>(e2s(PAR_DISC_TYPE)).compare(e2s(USER_DEFINED_PAR)) == 0)
        {
            (*sim)->setPartDiscSchemeUser(_reader.template vector<double>(e2s(PAR_DISC_VECTOR)));
            log::emit<Info>() << "Using a user defined particle discretization scheme." << log::endl;
        }

        log::emit<Debug1>() << "Particle discretization scheme set!" << log::endl;
        // ============================================================================================================
    }

    // ============================================================================================================
    //    Store information on writing output
    // ============================================================================================================
    _reader.setGroup(e2s(GRP_IN_SOLVER));
    _writeSolutionTimes        = _reader.template scalar<int>(e2s(WRITE_SOLUTION_TIMES));
    _writeSolutionColumnOutlet = _reader.template scalar<int>(e2s(WRITE_SOLUTION_COLUMN_OUTLET));
    _writeSolutionColumnInlet  = _reader.template scalar<int>(e2s(WRITE_SOLUTION_COLUMN_INLET));
    _writeSolutionAll          = _reader.template scalar<int>(e2s(WRITE_SOLUTION_ALL));
    _writeSensColumnOutlet     = _reader.template scalar<int>(e2s(WRITE_SENS_COLUMN_OUTLET));
    _writeSensAll              = _reader.template scalar<int>(e2s(WRITE_SENS_ALL));
    // ============================================================================================================
    
    // ============================================================================================================
    //    Set number of OpenMP threads
    // ============================================================================================================
#ifdef _OPENMP    
    if (_reader.exists(e2s(NTHREADS)))
    {
        const int nThreads = _reader.template scalar<int>(e2s(NTHREADS));
        omp_set_num_threads(nThreads);
    }
#endif
    // ============================================================================================================
}



template <typename reader_t, typename writer_t>
void CadetCS<reader_t, writer_t>::configureInlet()
{
    for (std::vector<Simulator*>::iterator sim = _sim.begin(); sim < _sim.end(); ++sim)
    {
        std::size_t sim_idx = sim - _sim.begin();

        // ============================================================================================================
        //    Set section times vector
        // ============================================================================================================
        _reader.setGroup(e2s(GRP_IN_INLET));
        const std::vector<double>& sectionTimes = _reader.template vector<double>(e2s(SECTION_TIMES));

        std::vector<int> sectionContinuityInt;
        if (_reader.exists(e2s(SECTION_CONTINUITY)))
        {
            sectionContinuityInt = _reader.template vector<int>(e2s(SECTION_CONTINUITY));
            log::emit<Debug1>() << "Section continuity vector present" << log::endl;
        }
        else
        {
            sectionContinuityInt = std::vector<int>(sectionTimes.size() - 2, 0);
            log::emit<Debug1>() << "Section continuity vector not present, assuming discontinuous sections" << log::endl;
        }
        std::vector<bool> sectionContinuity;
        sectionContinuity.insert(sectionContinuity.end(), sectionContinuityInt.begin(), sectionContinuityInt.end());

        _inletParams.at(sim_idx)->sectionTimes.assign(sectionTimes.begin(), sectionTimes.end());
        (*sim)->setSectionTimes(sectionTimes, sectionContinuity);

        log::emit<Debug1>() << "Section times vector set!" << log::endl;
        // ============================================================================================================


        // ============================================================================================================
        //    Set the inlet concentration function and related parameters
        // ============================================================================================================
        for (std::size_t sec = 0; sec < _nsec; ++sec)
        {
            _oss.str("");
            _oss << e2s(GRP_IN_INLET) << "/sec_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << sec;
            _reader.setGroup(_oss.str());
            std::vector<double> const_coeff = _reader.template vector<double>(inletParams::CONST_COEFF);
            std::vector<double> lin_coeff   = _reader.template vector<double>(inletParams::LIN_COEFF);
            std::vector<double> quad_coeff  = _reader.template vector<double>(inletParams::QUAD_COEFF);
            std::vector<double> cube_coeff  = _reader.template vector<double>(inletParams::CUBE_COEFF);
            _inletParams.at(sim_idx)->const_coeff.at(sec).assign(const_coeff.begin(), const_coeff.end());
            _inletParams.at(sim_idx)->lin_coeff  .at(sec).assign(lin_coeff.begin(),   lin_coeff.end());
            _inletParams.at(sim_idx)->quad_coeff .at(sec).assign(quad_coeff.begin(),  quad_coeff.end());
            _inletParams.at(sim_idx)->cube_coeff .at(sec).assign(cube_coeff.begin(),  cube_coeff.end());
        }

        (*sim)->setMaxSensInletParams(_inletParams.at(sim_idx)->getMaxParams());
        (*sim)->setInletProfile(_inletParams.at(sim_idx));

        log::emit<Debug1>() << "Inlet concentration function and related parameters set!" << log::endl;
        // ============================================================================================================
    }
}



template <typename reader_t, typename writer_t>
void CadetCS<reader_t, writer_t>::setParameters()
{
    for (std::vector<Simulator*>::iterator sim = _sim.begin(); sim < _sim.end(); ++sim)
    {
        // ============================================================================================================
        //    Setting all parameters
        // ============================================================================================================

        // Chromatography model parameters
        _reader.setGroup(e2s(GRP_IN_MODEL));

        // Check for section dependent column dispersion
        bool secDep = _reader.isVector(e2s(COL_DISPERSION));
        (*sim)->setParameterSectionDependent(COL_DISPERSION, secDep);
        if (secDep)
        {
            for (std::size_t i = 0; i < _nsec; ++i)
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(COL_DISPERSION), i), COL_DISPERSION, -1, i);
        }
        else
        {
            (*sim)->setParameterValue(_reader.template scalar<double>(e2s(COL_DISPERSION)), COL_DISPERSION);
        }

        // Check for section dependent velocity
        secDep = _reader.isVector(e2s(VELOCITY));
        (*sim)->setParameterSectionDependent(VELOCITY, secDep);
        if (secDep)
        {
            for (int i = 0; i < _nsec; ++i)
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(VELOCITY), i), VELOCITY, -1, i);
        }
        else
        {
            (*sim)->setParameterValue(_reader.template scalar<double>(e2s(VELOCITY)), VELOCITY);
        }

        (*sim)->setParameterValue(_reader.template scalar<double>(e2s(COL_LENGTH)  ), COL_LENGTH);
        (*sim)->setParameterValue(_reader.template scalar<double>(e2s(COL_POROSITY)), COL_POROSITY);
        (*sim)->setParameterValue(_reader.template scalar<double>(e2s(PAR_RADIUS)  ), PAR_RADIUS);
        (*sim)->setParameterValue(_reader.template scalar<double>(e2s(PAR_POROSITY)), PAR_POROSITY);

        // Vectorial parameters
        const std::vector<double> filmDiff = _reader.template vector<double>(e2s(FILM_DIFFUSION));
        if ((filmDiff.size() == _nsec * _ncomp) && (_nsec > 1))
        {
            (*sim)->setParameterSectionDependent(FILM_DIFFUSION, true);
            for (std::size_t i = 0; i < _nsec; ++i)
            {
                for (std::size_t comp = 0; comp < _ncomp; ++comp)
                    (*sim)->setParameterValue(filmDiff[comp + i * _ncomp], FILM_DIFFUSION, comp, i);
            }
        }
        else
        {
            (*sim)->setParameterSectionDependent(FILM_DIFFUSION, false);
            for (std::size_t comp = 0; comp < _ncomp; ++comp)
                (*sim)->setParameterValue(filmDiff[comp], FILM_DIFFUSION, comp);
        }

        const std::vector<double> parDiff = _reader.template vector<double>(e2s(PAR_DIFFUSION));
        if ((parDiff.size() == _nsec * _ncomp) && (_nsec > 1))
        {
            (*sim)->setParameterSectionDependent(PAR_DIFFUSION, true);
            for (std::size_t i = 0; i < _nsec; ++i)
            {
                for (std::size_t comp = 0; comp < _ncomp; ++comp)
                    (*sim)->setParameterValue(parDiff[comp + i * _ncomp], PAR_DIFFUSION, comp, i);
            }
        }
        else
        {
            (*sim)->setParameterSectionDependent(PAR_DIFFUSION, false);
            for (std::size_t comp = 0; comp < _ncomp; ++comp)
                (*sim)->setParameterValue(parDiff[comp], PAR_DIFFUSION, comp);
        }

        const std::vector<double> parSurfDiff = _reader.template vector<double>(e2s(PAR_SURFDIFFUSION));
        if ((parSurfDiff.size() == _nsec * _ncomp) && (_nsec > 1))
        {
            (*sim)->setParameterSectionDependent(PAR_SURFDIFFUSION, true);
            for (std::size_t i = 0; i < _nsec; ++i)
            {
                for (std::size_t comp = 0; comp < _ncomp; ++comp)
                    (*sim)->setParameterValue(parSurfDiff[comp + i * _ncomp], PAR_SURFDIFFUSION, comp, i);
            }
        }
        else
        {
            (*sim)->setParameterSectionDependent(PAR_SURFDIFFUSION, false);
            for (std::size_t comp = 0; comp < _ncomp; ++comp)
                (*sim)->setParameterValue(parSurfDiff[comp], PAR_SURFDIFFUSION, comp);
        }

        // Adsorption model parameters
        _reader.setGroup(e2s(GRP_IN_ADSORPTION));
        (*sim)->setKineticAdsorptionModel(_reader.template scalar<int>(e2s(IS_KINETIC)));

        switch (_adsType)
        {
        case LINEAR:
            for (std::size_t comp = 0; comp < _ncomp; ++comp) // vectorial parameters
            {
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(LIN_KA),   comp), LIN_KA, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(LIN_KD),   comp), LIN_KD, comp);
            }
            break;
        case MULTI_COMPONENT_LANGMUIR:
            for (std::size_t comp = 0; comp < _ncomp; ++comp) // vectorial parameters
            {
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(MCL_KA),   comp), MCL_KA, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(MCL_KD),   comp), MCL_KD, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(MCL_QMAX), comp), MCL_QMAX, comp);
            }
            break;
        case MOBILE_PHASE_MODULATORS:
            for (std::size_t comp = 0; comp < _ncomp; ++comp) // vectorial parameters
            {
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(MPM_KA),    comp), MPM_KA, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(MPM_KD),    comp), MPM_KD, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(MPM_QMAX),  comp), MPM_QMAX, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(MPM_BETA),  comp), MPM_BETA, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(MPM_GAMMA), comp), MPM_GAMMA, comp);
            }
            break;
        case STERIC_MASS_ACTION:
            (*sim)->setParameterValue(_reader.template scalar<double>(e2s(SMA_LAMBDA)), SMA_LAMBDA);
            for (std::size_t comp = 0; comp < _ncomp; ++comp) // vectorial parameters
            {
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(SMA_KA),    comp), SMA_KA, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(SMA_KD),    comp), SMA_KD, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(SMA_NU),    comp), SMA_NU, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(SMA_SIGMA), comp), SMA_SIGMA, comp);
            }
            break;
        case SELF_ASSOCIATION:
            (*sim)->setParameterValue(_reader.template scalar<double>(e2s(SAI_LAMBDA)), SAI_LAMBDA);
            for (std::size_t comp = 0; comp < _ncomp; ++comp) // vectorial parameters
            {
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(SAI_KA1),   comp), SAI_KA1, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(SAI_KA2),   comp), SAI_KA2, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(SAI_KD),    comp), SAI_KD, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(SAI_NU),    comp), SAI_NU, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(SAI_SIGMA), comp), SAI_SIGMA, comp);
            }
            break;
        case EXTERNAL_LANGMUIR:
            for (std::size_t comp = 0; comp < _ncomp; ++comp) // vectorial parameters
            {
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTL_KA),       comp), EXTL_KA, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTL_KA_T),     comp), EXTL_KA_T, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTL_KA_TT),    comp), EXTL_KA_TT, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTL_KA_TTT),   comp), EXTL_KA_TTT, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTL_KD),       comp), EXTL_KD, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTL_KD_T),     comp), EXTL_KD_T, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTL_KD_TT),    comp), EXTL_KD_TT, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTL_KD_TTT),   comp), EXTL_KD_TTT, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTL_QMAX),     comp), EXTL_QMAX, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTL_QMAX_T),   comp), EXTL_QMAX_T, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTL_QMAX_TT),  comp), EXTL_QMAX_TT, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTL_QMAX_TTT), comp), EXTL_QMAX_TTT, comp);
            }
            (*sim)->setExternalProfile(&_externalProfile);
            break;
        case EXTERNAL_STERIC_MASS_ACTION:
            (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTSMA_LAMBDA)),     EXTSMA_LAMBDA);
            (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTSMA_LAMBDA_T)),   EXTSMA_LAMBDA_T);
            (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTSMA_LAMBDA_TT)),  EXTSMA_LAMBDA_TT);
            (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTSMA_LAMBDA_TTT)), EXTSMA_LAMBDA_TTT);
            for (std::size_t comp = 0; comp < _ncomp; ++comp) // vectorial parameters
            {
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTSMA_KA),       comp), EXTSMA_KA, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTSMA_KA_T),     comp), EXTSMA_KA_T, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTSMA_KA_TT),    comp), EXTSMA_KA_TT, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTSMA_KA_TTT),   comp), EXTSMA_KA_TTT, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTSMA_KD),       comp), EXTSMA_KD, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTSMA_KD_T),     comp), EXTSMA_KD_T, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTSMA_KD_TT),    comp), EXTSMA_KD_TT, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTSMA_KD_TTT),   comp), EXTSMA_KD_TTT, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTSMA_NU),       comp), EXTSMA_NU, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTSMA_NU_T),     comp), EXTSMA_NU_T, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTSMA_NU_TT),    comp), EXTSMA_NU_TT, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTSMA_NU_TTT),   comp), EXTSMA_NU_TTT, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTSMA_SIGMA),    comp), EXTSMA_SIGMA, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTSMA_SIGMA_T),  comp), EXTSMA_SIGMA_T, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTSMA_SIGMA_TT), comp), EXTSMA_SIGMA_TT, comp);
                (*sim)->setParameterValue(_reader.template scalar<double>(e2s(EXTSMA_SIGMA_TTT),comp), EXTSMA_SIGMA_TTT, comp);
            }
            (*sim)->setExternalProfile(&_externalProfile);
            break;
        default:
            throw CadetException("Wrong adsorption type specified");
        }

        log::emit<Debug1>() << "All parameters set!" << log::endl;
        // ============================================================================================================
    }
}



template <typename reader_t, typename writer_t>
void CadetCS<reader_t, writer_t>::setSensitivities()
{
    // ============================================================================================================
    //    Sensitivities selection
    // ============================================================================================================
    ParameterName paramId;
    std::size_t sim;
    double paramValue;
    double fdDeltaAbs;
    bool isInletParam;

    for (std::size_t sens = 0; sens < _nsens; ++sens)
    {
        sim = 1 + _nAddFdRuns * sens; // Current simulator index

        _oss.str("");
        _oss << e2s(GRP_IN_SENSITIVITY) << "/param_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << sens;
        _reader.setGroup(_oss.str());

        // Read information on sensitivities from input file
        std::string pname    = _reader.template scalar<std::string>(e2s(SENS_NAME));
        int         comp     = _reader.template scalar<int>        (e2s(SENS_COMP));
        int         sec      = _reader.template scalar<int>        (e2s(SENS_SECTION));
        double      absTolS  = _reader.template scalar<double>     (e2s(SENS_ABSTOL));
        double      fdDelta  = _reader.template scalar<double>     (e2s(SENS_FD_DELTA));

        // Store sensitive parameter for data extraction
        _sensNames.push_back(pname);

        isInletParam = false;
        // If we have an inlet parameter we need to set pname to "INLET_PARAMETER"
        if (pname == inletParams::CONST_COEFF     || pname == inletParams::LIN_COEFF   ||
                pname == inletParams::QUAD_COEFF  || pname == inletParams::CUBE_COEFF)
        {
            isInletParam = true;
            pname = "INLET_PARAMETER";
        }
        // Convert parameter name string to enum
        s2e(pname, paramId);

        switch (_sensMethod)
        {
        // ============================================================================================================
        //    AD Sensitivities
        // ============================================================================================================
        case ALGORITHMIC_DIFFERENTIATION_1:
            if (isInletParam)
            {
                // Figure out the index of the inlet parameter
                comp = _inletParams.at(0)->getIndex(_sensNames.back(), comp, sec);
            }

            // Set the parameter sensitive
            _sim.at(0)->setSensitiveParameter(paramId, absTolS, comp, sec);

            break;
        // ============================================================================================================


        // ============================================================================================================
        //    FD first order aproximation sensitivities (forward difference, 1st order)
        // ============================================================================================================
        case FINITE_DIFFERENCES_1:
            if (isInletParam)
            {
                // Care about inlet parameters
                // If we have an inlet parameter we need to change the parameter value in the right inletParams instance
                if      (_sensNames.back() == inletParams::CONST_COEFF)
                {
                    _inletParams.at(sim)->const_coeff.at(sec).at(comp) *= 1 + fdDelta;
                    fdDeltaAbs = fdDelta * _inletParams.at(sim)->const_coeff.at(sec).at(comp);
                }
                else if (_sensNames.back() == inletParams::LIN_COEFF)
                {
                    _inletParams.at(sim)->lin_coeff  .at(sec).at(comp) *= 1 + fdDelta;
                    fdDeltaAbs = fdDelta * _inletParams.at(sim)->lin_coeff.at(sec).at(comp);
                }
                else if (_sensNames.back() == inletParams::QUAD_COEFF)
                {
                    _inletParams.at(sim)->quad_coeff .at(sec).at(comp) *= 1 + fdDelta;
                    fdDeltaAbs = fdDelta * _inletParams.at(sim)->quad_coeff.at(sec).at(comp);
                }
                else if (_sensNames.back() == inletParams::CUBE_COEFF)
                {
                    _inletParams.at(sim)->cube_coeff .at(sec).at(comp) *= 1 + fdDelta;
                    fdDeltaAbs = fdDelta * _inletParams.at(sim)->cube_coeff.at(sec).at(comp);
                }
            }
            else
            {
                // Care about chromatography and adsorption parameters
                paramValue = _sim.at(sim)->getParameterValue(paramId, comp, sec);
                _sim.at(sim)->setParameterValue(paramValue * (1 + fdDelta), paramId, comp, sec);
                fdDeltaAbs = fdDelta * paramValue;
            }

            log::emit<Debug1>() << "Sensitivity #" << sens << " for FD1 set." << log::endl;
            break;
        // ============================================================================================================


        // ============================================================================================================
        //    FD second order approximation sensitivities (central difference, 2nd order)
        // ============================================================================================================
        case FINITE_DIFFERENCES_2:
            if (isInletParam)
            {
                // Care about inlet parameters
                // If we have an inlet parameter we need to change the parameter value in the right inletParams instance
                if      (_sensNames.back() == inletParams::CONST_COEFF)
                {
                    _inletParams.at(sim + 0)->const_coeff.at(sec).at(comp) *= 1 + fdDelta;
                    _inletParams.at(sim + 1)->const_coeff.at(sec).at(comp) *= 1 - fdDelta;
                    fdDeltaAbs = fdDelta * _inletParams.at(sim)->const_coeff.at(sec).at(comp);
                }
                else if (_sensNames.back() == inletParams::LIN_COEFF)
                {
                    _inletParams.at(sim + 0)->lin_coeff  .at(sec).at(comp) *= 1 + fdDelta;
                    _inletParams.at(sim + 1)->lin_coeff  .at(sec).at(comp) *= 1 - fdDelta;
                    fdDeltaAbs = fdDelta * _inletParams.at(sim)->lin_coeff.at(sec).at(comp);
                }
                else if (_sensNames.back() == inletParams::QUAD_COEFF)
                {
                    _inletParams.at(sim + 0)->quad_coeff .at(sec).at(comp) *= 1 + fdDelta;
                    _inletParams.at(sim + 1)->quad_coeff .at(sec).at(comp) *= 1 - fdDelta;
                    fdDeltaAbs = fdDelta * _inletParams.at(sim)->quad_coeff.at(sec).at(comp);
                }
                else if (_sensNames.back() == inletParams::CUBE_COEFF)
                {
                    _inletParams.at(sim + 0)->cube_coeff .at(sec).at(comp) *= 1 + fdDelta;
                    _inletParams.at(sim + 1)->cube_coeff .at(sec).at(comp) *= 1 - fdDelta;
                    fdDeltaAbs = fdDelta * _inletParams.at(sim)->cube_coeff.at(sec).at(comp);
                }
            }
            else
            {
                // Care about chromatography and adsorption parameters
                paramValue = _sim.at(sim)->getParameterValue(paramId, comp, sec);
                _sim.at(sim + 0)->setParameterValue(paramValue * (1 + fdDelta), paramId, comp, sec);
                _sim.at(sim + 1)->setParameterValue(paramValue * (1 - fdDelta), paramId, comp, sec);
                fdDeltaAbs = fdDelta * paramValue;
            }

            log::emit<Debug1>() << "Sensitivity #" << sens << " for FD2 set." << log::endl;
            break;
        // ============================================================================================================


        // ============================================================================================================
        //    FD fourth order approximation sensitivities (central difference, 4th order)
        // ============================================================================================================
        case FINITE_DIFFERENCES_4:
            if (isInletParam)
            {
                // Care about inlet parameters
                // If we have an inlet parameter we need to change the parameter value in the right inletParams instance
                if      (_sensNames.back() == inletParams::CONST_COEFF)
                {
                    _inletParams.at(sim + 0)->const_coeff.at(sec).at(comp) *= 1 + fdDelta * 2;
                    _inletParams.at(sim + 1)->const_coeff.at(sec).at(comp) *= 1 + fdDelta;
                    _inletParams.at(sim + 2)->const_coeff.at(sec).at(comp) *= 1 - fdDelta;
                    _inletParams.at(sim + 3)->const_coeff.at(sec).at(comp) *= 1 - fdDelta * 2;
                    fdDeltaAbs = fdDelta * _inletParams.at(sim)->const_coeff.at(sec).at(comp);
                }
                else if (_sensNames.back() == inletParams::LIN_COEFF)
                {
                    _inletParams.at(sim + 0)->lin_coeff.at(sec).at(comp) *= 1 + fdDelta * 2;
                    _inletParams.at(sim + 1)->lin_coeff.at(sec).at(comp) *= 1 + fdDelta;
                    _inletParams.at(sim + 2)->lin_coeff.at(sec).at(comp) *= 1 - fdDelta;
                    _inletParams.at(sim + 3)->lin_coeff.at(sec).at(comp) *= 1 - fdDelta * 2;
                    fdDeltaAbs = fdDelta * _inletParams.at(sim)->lin_coeff.at(sec).at(comp);
                }
                else if (_sensNames.back() == inletParams::QUAD_COEFF)
                {
                    _inletParams.at(sim + 0)->quad_coeff.at(sec).at(comp) *= 1 + fdDelta * 2;
                    _inletParams.at(sim + 1)->quad_coeff.at(sec).at(comp) *= 1 + fdDelta;
                    _inletParams.at(sim + 2)->quad_coeff.at(sec).at(comp) *= 1 - fdDelta;
                    _inletParams.at(sim + 3)->quad_coeff.at(sec).at(comp) *= 1 - fdDelta * 2;
                    fdDeltaAbs = fdDelta * _inletParams.at(sim)->quad_coeff.at(sec).at(comp);
                }
                else if (_sensNames.back() == inletParams::CUBE_COEFF)
                {
                    _inletParams.at(sim + 0)->cube_coeff.at(sec).at(comp) *= 1 + fdDelta * 2;
                    _inletParams.at(sim + 1)->cube_coeff.at(sec).at(comp) *= 1 + fdDelta;
                    _inletParams.at(sim + 2)->cube_coeff.at(sec).at(comp) *= 1 - fdDelta;
                    _inletParams.at(sim + 3)->cube_coeff.at(sec).at(comp) *= 1 - fdDelta * 2;
                    fdDeltaAbs = fdDelta * _inletParams.at(sim)->cube_coeff.at(sec).at(comp);
                }
            }
            else
            {
                // Care about chromatography and adsorption parameters
                paramValue = _sim.at(sim)->getParameterValue(paramId, comp, sec);
                _sim.at(sim + 0)->setParameterValue(paramValue * (1 + fdDelta * 2), paramId, comp, sec);
                _sim.at(sim + 1)->setParameterValue(paramValue * (1 + fdDelta)    , paramId, comp, sec);
                _sim.at(sim + 2)->setParameterValue(paramValue * (1 - fdDelta)    , paramId, comp, sec);
                _sim.at(sim + 3)->setParameterValue(paramValue * (1 - fdDelta * 2), paramId, comp, sec);
                fdDeltaAbs = fdDelta * paramValue;
            }

            log::emit<Debug1>() << "Sensitivity #" << sens << " for FD4 set." << log::endl;
            break;
        // ============================================================================================================

        default:
            throw CadetException("Wrong sensitivity comptation method specified!");
        }

        // Store sensitive parameter associated data for later extraction
        _sensParamIds.push_back(paramId);
        _sensComps.push_back(comp);
        _sensSecs.push_back(sec);
        _sensFdDeltaAbs.push_back(fdDeltaAbs);
    }

    log::emit<Debug1>() << "All sensitivities set!" << log::endl;
    // ============================================================================================================
}


template <typename reader_t, typename writer_t>
void CadetCS<reader_t, writer_t>::simulate()
{
    // ============================================================================================================
    //    Create and set initial conditions vector
    // ============================================================================================================
    _reader.setGroup(e2s(GRP_IN_MODEL));
    _sim.at(0)->initialize(_reader.template vector<double>(e2s(INIT_C)), _reader.template vector<double>(e2s(INIT_Q)));

    // ============================================================================================================

    // ============================================================================================================
    //    Some console output
    // ============================================================================================================
    log::emit<Info>() << "File '" << _filename << "' was successfully read." << log::endl;
    log::emit<Debug1>() << "Number of model sensitivities: " << _sim.at(0)->getNSensModelParams() << log::endl;
    log::emit<Debug1>() << "Number of inlet sensitivities: " << _sim.at(0)->getNSensInletParams() << log::endl;
    // ============================================================================================================


    // ============================================================================================================
    //    Solve the main system
    // ============================================================================================================
    _sim.at(0)->integrate();
    // ============================================================================================================



    // Get solution times vector
    _sim.at(0)->getSolutionTimes(_solutionTimes);

    // Set the solution times for finite difference sensitivities,
    // to prevent the solution data vectors from having different sizes!
    if (_sensMethod != ALGORITHMIC_DIFFERENTIATION_1)
    {
        std::size_t run;
        std::size_t sens;
        for (std::vector<Simulator*>::iterator sim = _sim.begin() + 1; sim < _sim.end(); ++sim)
        {
            run = sim - _sim.begin();
            sens = (run - 1) / _nAddFdRuns;
            // ============================================================================================================
            //    Create and set initial conditions vector
            // ============================================================================================================
            (*sim)->setSolutionTimes(_solutionTimes);
            (*sim)->initialize(_reader.template vector<double>(e2s(INIT_C)), _reader.template vector<double>(e2s(INIT_Q)));

            log::emit<Debug1>() << "Initial conditions for sens. " << _sensNames.at(sens) << "[comp " << _sensComps.at(sens)
                    << ", sec " << _sensSecs.at(sens) << "]" <<" in FD run " << run << " set!" << log::endl;
            // ============================================================================================================


            // ============================================================================================================
            //    Solve the finite difference systems
            // ============================================================================================================
            log::emit<Info>() << "FD run #" << run << ": " << _sensNames.at(sens) << "[comp " << _sensComps.at(sens) 
                    << ", sec " << _sensSecs.at(sens) << "]" << log::endl;
            (*sim)->integrate();
            // ============================================================================================================

            std::vector<double> tmp;
            (*sim)->getSolutionTimes(tmp);
        }
    }
}



template <typename reader_t, typename writer_t>
void CadetCS<reader_t, writer_t>::extractSolutionData()
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    std::size_t sim;
    std::size_t i;

    // ============================================================================================================
    //    Extract solution
    // ============================================================================================================
    // Get column, particle and boundary solutions
    _sim.at(0)->getColSolutions(_colSolutions);
    _sim.at(0)->getParSolutions(_parSolutions);
    _sim.at(0)->getBndSolutions(_bndSolutions);


    // Get chromatograms and inlet concentrations
    for (std::size_t comp = 0; comp < _ncomp; ++comp)
    {
        std::vector<double> tmp;
        _solutionColumnOutlet.push_back(tmp);
        _sim.at(0)->getSolutionColumnOutlet(_solutionColumnOutlet.at(comp), comp);

        _solutionColumnInlet.push_back(tmp);
        _sim.at(0)->getSolutionColumnInlet(_solutionColumnInlet.at(comp), comp);
    }
    // ============================================================================================================


    // ============================================================================================================
    //    Extract sensitivities (depending on sensitivity computation method)
    // ============================================================================================================

    // Common things to do
    std::vector<double> tmp;
    std::vector<double> empty;

    std::vector< std::vector<double> > colSolutionsSens;
    std::vector< std::vector<double> > parSolutionsSens;
    std::vector< std::vector<double> > bndSolutionsSens;

    std::vector<double> sensFdDeltaAbsCol;
    std::vector<double> sensFdDeltaAbsPar;
    std::vector<double> sensFdDeltaAbsBnd;

    for (std::size_t sens = 0; sens < _nsens; ++sens)
    {
        std::vector<std::vector<double> > vecvec;
        _sensitivityColumnOutlet.push_back(vecvec);
        for (std::size_t comp = 0; comp < _ncomp; ++comp)
        {
            _sensitivityColumnOutlet.at(sens).push_back(empty);
        }
    }

    // Build up vectors to be used in FD computations
    if (_sensMethod != ALGORITHMIC_DIFFERENTIATION_1)
    {
        // Create a vector holding a vector with the solutions of each simulator instance
        // i.e. the vector contains _nAddFdRuns - 1 solution-vectors
        for (std::size_t addRun = 0; addRun < _nAddFdRuns; ++addRun)
        {
            colSolutionsSens.push_back(empty);
            parSolutionsSens.push_back(empty);
            bndSolutionsSens.push_back(empty);
        }

        for (std::size_t sens = 0; sens < _nsens; ++sens)
        {
            // Extract solution f(x) to the final sensitivitiy vectors (needed for FD1)
            _sim.at(0)->getColSolutions(tmp);
            _colSensitivities.insert(_colSensitivities.end(), tmp.begin(), tmp.end());
            _sim.at(0)->getParSolutions(tmp);
            _parSensitivities.insert(_parSensitivities.end(), tmp.begin(), tmp.end());
            _sim.at(0)->getBndSolutions(tmp);
            _bndSensitivities.insert(_bndSensitivities.end(), tmp.begin(), tmp.end());

            // Extract all other solutions [f(x + Delta), etc.] from remaining simulator instances
            for (std::size_t addRun = 0; addRun < _nAddFdRuns; ++addRun)
            {
                sim = 1 + sens * _nAddFdRuns + addRun;

                _sim.at(sim)->getColSolutions(tmp);
                colSolutionsSens.at(addRun).insert(colSolutionsSens.at(addRun).end(), tmp.begin(), tmp.end());
                _sim.at(sim)->getParSolutions(tmp);
                parSolutionsSens.at(addRun).insert(parSolutionsSens.at(addRun).end(), tmp.begin(), tmp.end());
                _sim.at(sim)->getBndSolutions(tmp);
                bndSolutionsSens.at(addRun).insert(bndSolutionsSens.at(addRun).end(), tmp.begin(), tmp.end());
            }

            // Create vectors containing the finite difference deltas
            tmp.assign(_colSolutions.size(), _sensFdDeltaAbs.at(sens));
            sensFdDeltaAbsCol.insert(sensFdDeltaAbsCol.end(), tmp.begin(), tmp.end());
            tmp.assign(_parSolutions.size(), _sensFdDeltaAbs.at(sens));
            sensFdDeltaAbsPar.insert(sensFdDeltaAbsPar.end(), tmp.begin(), tmp.end());
            tmp.assign(_bndSolutions.size(), _sensFdDeltaAbs.at(sens));
            sensFdDeltaAbsBnd.insert(sensFdDeltaAbsBnd.end(), tmp.begin(), tmp.end());
        }
    }


    switch (_sensMethod)
    {
    case ALGORITHMIC_DIFFERENTIATION_1:

        // Get column, particle and boundary sensitivities
        _sim.at(0)->getColSensitivities(_colSensitivities);
        _sim.at(0)->getParSensitivities(_parSensitivities);
        _sim.at(0)->getBndSensitivities(_bndSensitivities);

        // Get sensograms
        for (std::size_t sens = 0; sens < _nsens; ++sens)
            for (std::size_t comp = 0; comp < _ncomp; ++comp)
                _sim.at(0)->getSensitivityColumnOutlet(_sensitivityColumnOutlet.at(sens).at(comp),
                        comp, _sensParamIds.at(sens), _sensComps.at(sens), _sensSecs.at(sens));
        break;

    case FINITE_DIFFERENCES_1: // Forward difference, 1st order

        for (std::size_t sens = 0; sens < _nsens; ++sens)
        {
            sim = 1 + _nAddFdRuns * sens;

            // Get sensograms
            for (std::size_t comp = 0; comp < _ncomp; ++comp)
            {
                // Get 'f(x + Delta)' from first simulator instance
                _sim.at(sim)->getSolutionColumnOutlet(_sensitivityColumnOutlet.at(sens).at(comp), comp);

                // Compute 'f(x + Delta) - f(x) / Delta' using lambda expression (C++11 standard)
                i = 0;
                std::for_each(_sensitivityColumnOutlet.at(sens).at(comp).begin(), _sensitivityColumnOutlet.at(sens).at(comp).end(),
                    // Here the lambda function starts
                    [&] (double& el) { el = (el - _solutionColumnOutlet.at(comp).at(i)) / _sensFdDeltaAbs.at(sens); ++i; });
                    // [&]: We catch out-of-scope variables by reference; The element we operate on is (double& el)
            }
            log::emit<Debug1>() << "Sensitivities for " << Color::red << _sensNames.at(sens) <<
                    Color::reset << " at column outlet extracted!" << log::endl;
        }

        // Compute column, particle and boundary sensitivities
        // 'f(x + Delta) - f(x) / Delta'
        i = 0;
        std::for_each(_colSensitivities.begin(), _colSensitivities.end(),
                [&] (double& el) { el = (colSolutionsSens.at(0).at(i) - el) / sensFdDeltaAbsCol.at(i); ++i; });
        i = 0;
        std::for_each(_parSensitivities.begin(), _parSensitivities.end(),
                [&] (double& el) { el = (parSolutionsSens.at(0).at(i) - el) / sensFdDeltaAbsPar.at(i); ++i; });
        i = 0;
        std::for_each(_bndSensitivities.begin(), _bndSensitivities.end(),
                [&] (double& el) { el = (bndSolutionsSens.at(0).at(i) - el) / sensFdDeltaAbsBnd.at(i); ++i; });

        log::emit<Debug1>() << "All sensitivities for all parameters extracted!" << log::endl;
        break;

    case FINITE_DIFFERENCES_2:  // Central difference, 2nd order

        for (std::size_t sens = 0; sens < _nsens; ++sens)
        {
            sim = 1 + _nAddFdRuns * sens;

            // Get sensograms
            for (std::size_t comp = 0; comp < _ncomp; ++comp)
            {
                // Get 'f(x + Delta)' from first simulator instance
                _sim.at(sim + 0)->getSolutionColumnOutlet(_sensitivityColumnOutlet.at(sens).at(comp), comp);
                // Get 'f(x - Delta)' from second simulator instance
                _sim.at(sim + 1)->getSolutionColumnOutlet(tmp, comp);

                // Compute 'f(x + Delta) - f(x - Delta) * 0.5 / Delta' using lambda expression (C++11 standard)
                i = 0;
                std::for_each(_sensitivityColumnOutlet.at(sens).at(comp).begin(), _sensitivityColumnOutlet.at(sens).at(comp).end(),
                    [&] (double& el) { el = (el - tmp.at(i)) * 0.5 / _sensFdDeltaAbs.at(sens); ++i; });
            }
            log::emit<Debug1>() << "Sensitivities for " << Color::red << _sensNames.at(sens) <<
                    Color::reset << " at column outlet extracted!" << log::endl;
        }

        // Compute column, particle and boundary sensitivities
        // 'f(x + Delta) - f(x - Delta) * 0.5 / Delta'
        i = 0;
        std::for_each(_colSensitivities.begin(), _colSensitivities.end(),
                [&] (double& el) { el = (colSolutionsSens.at(0).at(i) - colSolutionsSens.at(1).at(i)) * 0.5 / sensFdDeltaAbsCol.at(i); ++i; });
        i = 0;
        std::for_each(_parSensitivities.begin(), _parSensitivities.end(),
                [&] (double& el) { el = (parSolutionsSens.at(0).at(i) - parSolutionsSens.at(1).at(i)) * 0.5 / sensFdDeltaAbsPar.at(i); ++i; });
        i = 0;
        std::for_each(_bndSensitivities.begin(), _bndSensitivities.end(),
                [&] (double& el) { el = (bndSolutionsSens.at(0).at(i) - bndSolutionsSens.at(1).at(i)) * 0.5 / sensFdDeltaAbsBnd.at(i); ++i; });

        log::emit<Debug1>() << "All sensitivities for all parameters extracted!" << log::endl;
        break;

    case FINITE_DIFFERENCES_4:  // Central difference, 4th order

        for (std::size_t sens = 0; sens < _nsens; ++sens)
        {
            sim = 1 + _nAddFdRuns * sens;

            // Get sensograms
            for (std::size_t comp = 0; comp < _ncomp; ++comp)
            {
                // Get 'f(x + Delta * 2)' from first simulator instance
                _sim.at(sim + 0)->getSolutionColumnOutlet(_sensitivityColumnOutlet.at(sens).at(comp), comp);
                // Get 'f(x + Delta)' from second simulator instance
                _sim.at(sim + 1)->getSolutionColumnOutlet(tmp, comp);
                // Get 'f(x - Delta)' from third simulator instance
                std::vector<double> tmp2;
                _sim.at(sim + 2)->getSolutionColumnOutlet(tmp2, comp);
                // Get 'f(x - Delta * 2)' from fourth simulator instance
                std::vector<double> tmp3;
                _sim.at(sim + 3)->getSolutionColumnOutlet(tmp3, comp);

                // Compute '-f(x + Delta * 2) + 8 * f(x + Delta) - 8 * f(x - Delta) + f(x - Delta * 2) / (12 * Delta)'
                // using lambda expression (C++11 standard)
                i = 0;
                std::for_each(_sensitivityColumnOutlet.at(sens).at(comp).begin(), _sensitivityColumnOutlet.at(sens).at(comp).end(),
                    [&] (double& e) { e = (-e + 8 * tmp.at(i) - 8 * tmp2.at(i) + tmp3.at(i)) / (12 * _sensFdDeltaAbs.at(sens)); ++i; });
            }
            log::emit<Debug1>() << "Sensitivities for " << Color::red << _sensNames.at(sens) <<
                    Color::reset << " at column outlet extracted!" << log::endl;
        }

        // Compute column, particle and boundary sensitivities
        // '-f(x + Delta * 2) + 8 * f(x + Delta) - 8 * f(x - Delta) + f(x - Delta * 2) / (12 * Delta)'
        i = 0;
        std::for_each(_colSensitivities.begin(), _colSensitivities.end(),
                [&] (double& el) { el = (-colSolutionsSens.at(0).at(i) + 8 * colSolutionsSens.at(1).at(i) - 8 *
                        colSolutionsSens.at(2).at(i) + colSolutionsSens.at(3).at(i)) / (12 * sensFdDeltaAbsCol.at(i)); ++i; });
        i = 0;
        std::for_each(_parSensitivities.begin(), _parSensitivities.end(),
                [&] (double& el) { el = (-parSolutionsSens.at(0).at(i) + 8 * parSolutionsSens.at(1).at(i) - 8 *
                        parSolutionsSens.at(2).at(i) + parSolutionsSens.at(3).at(i)) / (12 * sensFdDeltaAbsPar.at(i)); ++i; });
        i = 0;
        std::for_each(_bndSensitivities.begin(), _bndSensitivities.end(),
                [&] (double& el) { el = (-bndSolutionsSens.at(0).at(i) + 8 * bndSolutionsSens.at(1).at(i) - 8 *
                        bndSolutionsSens.at(2).at(i) + bndSolutionsSens.at(3).at(i)) / (12 * sensFdDeltaAbsBnd.at(i)); ++i; });

        log::emit<Debug1>() << "All sensitivities for all parameters extracted!" << log::endl;
        break;

    default:
        break;
    }
    // ============================================================================================================
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}



template <typename reader_t, typename writer_t>
void CadetCS<reader_t, writer_t>::writeOutput()
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;
    log::emit<Info>() << "Writing output to file '" << _filename << "' ... ";

    // Create a writer for the output file
    _writer.openFile(_filename, "rw");

    // [HDF5 specific]:
    // Since HDF5 datasets cannot be deleted and also can not be updated in a satisfactory way,
    // the output group is first unlinked from the file and is then rewritten with new content.
    _writer.unlinkGroup(e2s(GRP_OUT));

    // [HDF5 specific]: In the output group not later extension is necessary.
    // Furthermore, chunking would slow down I/O processes a lot.
    _writer.extendibleFields(false);

    // [HDF5 specific]: Switch on compression using deflate filter (gzip)
    _writer.compressFields(true);



    // ============================================================================================================
    //    Write statistical data, timings, coordinates, etc.
    // ============================================================================================================

    // ...

    // ============================================================================================================



    // ============================================================================================================
    //    Write output of the simulation
    // ============================================================================================================
    _writer.setGroup(e2s(GRP_OUT_SOLUTION));

    // Write solution times
    if (_writeSolutionTimes)
        _writer.template vector<double>(e2s(SOLUTION_TIMES), _solutionTimes);

    // Write solution at column outlet (chromatograms) for all components
    if (_writeSolutionColumnOutlet)
    {
        for (std::size_t comp = 0; comp < _ncomp; ++comp)
        {
            _oss.str("");
            _oss << e2s(SOLUTION_COLUMN_OUTLET) << "_COMP_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << comp;
            _writer.template vector<double>(_oss.str(), _solutionColumnOutlet.at(comp));
        }
    }

    // Write solution at column inlet (boundary condition) for all components
    if (_writeSolutionColumnInlet)
    {
        for (std::size_t comp = 0; comp < _ncomp; ++comp)
        {
            _oss.str("");
            _oss << e2s(SOLUTION_COLUMN_INLET) << "_COMP_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << comp;
            _writer.template vector<double>(_oss.str(), _solutionColumnInlet.at(comp));
        }
    }

    // Set rank and dimensions for column, particle and boundary solution tensors
    std::size_t colRank = 3;
    std::vector<std::size_t> colDims(colRank, 0);
    colDims.at(0) = _solutionTimes.size();
    colDims.at(1) = _ncomp;
    colDims.at(2) = _ncol;

    std::size_t parRank = 5;
    std::vector<std::size_t> parDims(parRank, 0);
    parDims.at(0) = _solutionTimes.size();
    parDims.at(1) = _ncol;
    parDims.at(2) = _npar;
    parDims.at(3) = 2;
    parDims.at(4) = _ncomp;

    // Write whole solution
    if (_writeSolutionAll)
    {
        _writer.template tensor<double>(e2s(SOLUTION_COLUMN),   colRank, &colDims[0], _colSolutions);
        _writer.template tensor<double>(e2s(SOLUTION_PARTICLE), parRank, &parDims[0], _parSolutions);
        _writer.template tensor<double>(e2s(SOLUTION_BOUNDARY), colRank, &colDims[0], _bndSolutions);
    }

    // Only write sensitivity data when at least one sensitivity was computed.
    if (_nsens > 0)
    {
        // Write sensitivities at column outlet (sensograms) for all sensitivities and components
        if (_writeSensColumnOutlet)
        {
            for (std::size_t sens = 0; sens < _nsens; ++sens)
            {
                for (std::size_t comp = 0; comp < _ncomp; ++comp)
                {
                    // Set group
                    _oss.str("");
                    _oss << e2s(GRP_OUT_SENS) << "/param_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << sens;
                    _writer.setGroup(_oss.str());

                    // Write data
                    _oss.str("");
                    _oss << e2s(SENS_COLUMN_OUTLET) << "_COMP_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << comp;
                    _writer.template vector<double>(_oss.str(), _sensitivityColumnOutlet.at(sens).at(comp));
                }
            }
        }

        // Increase rank and dimensions for sensitivity tensors
        colRank++;
        colDims.insert(colDims.begin(), _nsens);
        parRank++;
        parDims.insert(parDims.begin(), _nsens);

        //  Write whole sensitivities
        _writer.setGroup(e2s(GRP_OUT_SENS));
        if (_writeSensAll)
        {
            _writer.template tensor<double>(e2s(SENS_COLUMN),   colRank, &colDims[0], _colSensitivities);
            _writer.template tensor<double>(e2s(SENS_PARTICLE), parRank, &parDims[0], _parSensitivities);
            _writer.template tensor<double>(e2s(SENS_BOUNDARY), colRank, &colDims[0], _bndSensitivities);
        }
    }

    // Close the file
    _writer.closeFile();

    log::emit<Info>() << "done." << log::endl;
    // ============================================================================================================
}

} // namespace cadet


