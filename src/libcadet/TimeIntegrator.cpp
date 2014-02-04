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
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cstdlib>

#include <idas/idas.h>
#include <idas/idas_impl.h>
#include <nvector/nvector_serial.h>

#include "TimeIntegrator.hpp"

#include "active.hpp"
#include "AdsorptionModel.hpp"
#include "ChromatographyModel.hpp"
#include "SchurSolver.hpp"
#include "JacobianData.hpp"
#include "WenoScheme.hpp"
#include "CadetConvenience.hpp"
#include "CadetLogger.hpp"
#include "SimpleProgress.hpp"



namespace cadet {


// Constructor
TimeIntegrator::TimeIntegrator(SimulatorPImpl& sim) :
    _writeAtUserTimes(false),
    _psim(sim),
    _cc(sim.getCadetConstants())
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    // Initialize with default values
    this->configure();
    this->configurePrinting();
    log::emit<Debug1>() << CURRENT_FUNCTION << ": Configured" << log::endl;

    // ========================================================================================
    //    Allocation and initialization of N_Vectors
    // ========================================================================================

    // Allocate an initialize state vectors
    _NV_y       = N_VNew_Serial(_cc.neq());
    _NV_yDot    = N_VNew_Serial(_cc.neq());

    N_VConst(0.0, _NV_y);      // _NV_y will be set to fulfill initial conditions later
    N_VConst(0.0, _NV_yDot);   // Initialize _NV_yDot to be zero everywhere (equilibrium)

    // Allocate and initialize temp space vectors
    _NV_temp1  = N_VNew_Serial(_cc.neq());
    _NV_temp2  = N_VNew_Serial(_cc.neq());

    N_VConst(0.0, _NV_temp1);
    N_VConst(0.0, _NV_temp2);

    // Allocate and initialize sensitivities vectors
    _NVp_yS     = N_VCloneVectorArray_Serial(NUMBER_DIRECTIONS, _NV_y);
    _NVp_ySDot  = N_VCloneVectorArray_Serial(NUMBER_DIRECTIONS, _NV_yDot);

    // ========================================================================================
    log::emit<Debug1>() << CURRENT_FUNCTION << ": N_Vector memory allocated" << log::endl;

    _idaMemBlock  = nullptr;

    // Allocate memory for AD computations
    _yAd    = new active[_cc.neq()];
    _resAd  = new active[_cc.neq()];
    // ========================================================================================

    // Misc initializations
    _sec = 0;

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


// Destructor
TimeIntegrator::~TimeIntegrator()
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    // ========================================================================================
    //    Destruction of N_Vectors
    // ========================================================================================
    N_VDestroy_Serial(_NV_y);
    N_VDestroy_Serial(_NV_yDot);

    N_VDestroy_Serial(_NV_temp1);
    N_VDestroy_Serial(_NV_temp2);

    N_VDestroyVectorArray_Serial(_NVp_yS, NUMBER_DIRECTIONS);
    N_VDestroyVectorArray_Serial(_NVp_ySDot, NUMBER_DIRECTIONS);
    // ========================================================================================

    delete [] _yAd;
    delete [] _resAd;

    IDAFree(&_idaMemBlock);

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


int TimeIntegrator::getJacAdDirs() const { return _psim.getSchurSolver().getJacobianData().jacAdDirs(); }
int TimeIntegrator::getDiagDir()   const { return _psim.getSchurSolver().getJacobianData().diagDir(); }



void TimeIntegrator::setInitialConditions(const std::vector<double>& initC, const std::vector<double>& initQ)
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    // ========================================================================================
    // Make shure the initial condition for the salt component in bound phase
    // is equivalent to lambda in case of SMA or SAI kinetics.
    // ========================================================================================
    double lambda;
    std::vector<double> initQ_wrt;
    initQ_wrt.assign(initQ.begin(), initQ.end());

    if (_psim.getAdsorptionModel().getAdsorptionType() == STERIC_MASS_ACTION)
        lambda = _psim.getAdsorptionModel().getValue<double>(SMA_LAMBDA);
    else if (_psim.getAdsorptionModel().getAdsorptionType() == SELF_ASSOCIATION)
        lambda = _psim.getAdsorptionModel().getValue<double>(SAI_LAMBDA);

    if (_psim.getAdsorptionModel().getAdsorptionType() == STERIC_MASS_ACTION ||
            _psim.getAdsorptionModel().getAdsorptionType() == SELF_ASSOCIATION)
        if (lambda != initQ.at(0))
        {
            initQ_wrt.at(0) = lambda;
            log::emit<Warning>() << CURRENT_FUNCTION << Color::yellow <<
                    ": Initial condition for salt component was automatically set "
                    "to parameter value LAMBDA! (must be equal)" << Color::reset << log::endl;
        }
    // ========================================================================================


    // ========================================================================================
    //    Fill state vector with initial conditions
    // ========================================================================================

    double* yc = _cc.offsetCol(_NV_y);
    for (int comp = 0; comp < _cc.ncomp(); ++comp)          // Iterate over components
        for (int col = 0; col < _cc.ncol(); ++col)          // Iterate over column cells
            // Fill in the initial conditions for all components
            _cc.colC<double>(yc, col, comp) = initC.at(comp);

    double* yp;
    for (int pblk = 0; pblk < _cc.npblk(); ++pblk)          // Iterate over particle blocks
    {
        // Skipping column part + act. blocks * particle parts
        yp = _cc.offsetPar(_NV_y, pblk);
        for (int par = 0; par < _cc.npar(); ++par)          // Iterate over particle cells
            for (int comp = 0; comp < _cc.ncomp(); ++comp)  // Iterate over components
            {
                _cc.parC<double>(yp, par, comp) = initC.at(comp);
                _cc.parQ<double>(yp, par, comp) = initQ_wrt.at(comp);
            }
    }

    // Skipping column part + no. particle blocks * particle parts
    double* yb = _cc.offsetBnd(_NV_y);
    for (int comp = 0; comp < _cc.ncomp(); ++comp)          // Iterate over components
        for (int bnd = 0; bnd < _cc.nbnd(); ++bnd)          // Iterate over boundary cells (same # as column cells)
            // Since boundary values are fluxes into the beads, these are all zero in the beginning
            _cc.bndJ<double>(yb, bnd, comp) = 0.0;

    // ========================================================================================

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


void TimeIntegrator::setSectionTimes(const std::vector<double>& sectionTimes) throw (CadetException)
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    _sectionTimes.assign(sectionTimes.begin(), sectionTimes.end());
    _sectionContinuity = std::vector<bool>(_sectionTimes.size() - 1, false);

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


void TimeIntegrator::setSectionTimes(const std::vector<double>& sectionTimes, const std::vector<bool>& sectionContinuity) throw (CadetException)
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    _sectionTimes.assign(sectionTimes.begin(), sectionTimes.end());
    _sectionContinuity.assign(sectionContinuity.begin(), sectionContinuity.end());

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


void TimeIntegrator::setSolutionTimes(const std::vector<double>& solutionTimes) throw (CadetException)
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    _writeAtUserTimes = true;
    _solutionTimes.assign(solutionTimes.begin(), solutionTimes.end());

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}




void TimeIntegrator::initializeIntegrator()
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    // ========================================================================================
    //    IDA initialization part 1
    // ========================================================================================

    // Check if IDAS object was already created
    if (_idaMemBlock != nullptr)
        IDAFree(&_idaMemBlock);

    // IDAS Step 4: Create IDAS object
    _idaMemBlock = IDACreate();

    // IDAS Step 4.1: Specify error handler function
    IDASetErrHandlerFn(_idaMemBlock, &cadet::TimeIntegrator::idaErrorHandlerWrapper, this);

    // IDAS Step 5: Initialize the solver
    IDAInit(_idaMemBlock, &cadet::TimeIntegrator::residualDaeWrapper, _t, _NV_y, _NV_yDot);

    // IDAS Step 6: Specify integration tolerances (SS: scalar/scalar)
    IDASStolerances(_idaMemBlock, _relTol, _absTol);


    // IDAS Step 7.1: Set optional inputs

    // Attach user data structure
    IDASetUserData(_idaMemBlock, &_psim);

    // Set maximum number of steps
    IDASetMaxNumSteps(_idaMemBlock, _maxSteps);

//    // Set maximum step size ///todo make that user choosable
//    IDASetMaxStep(_idaMemBlock, 50.0);

    // Specify the linear solver.
    IDAMem IDA_mem = static_cast<IDAMem>(_idaMemBlock);

    IDA_mem->ida_linit          = NULL;
    IDA_mem->ida_lsetup         = NULL;
    IDA_mem->ida_lsolve         = &cadet::TimeIntegrator::schurSolveWrapper;
    IDA_mem->ida_lperf          = NULL;
    IDA_mem->ida_lfree          = NULL;
    IDA_mem->ida_setupNonNull   = FALSE;
    IDA_mem->ida_lmem           = &_psim;

    // ========================================================================================
    log::emit<Debug1>() << CURRENT_FUNCTION << ": IDA init part 1 done" << log::endl;



    // ========================================================================================
    //    Set differential/algebraic variables
    // ========================================================================================

    // Allocate id (is_differential) vector
    N_Vector NV_id = N_VNew_Serial(_cc.neq());

    // Set the id vector to be 0 everywhere
    // 1 = differential variable, 0 = algebraic variable
    N_VConst(0, NV_id);

    // Set the differential variables to be 1

    // The mass transport equations in the column are differential
    double* idc = _cc.offsetCol(NV_id);
    for (int comp = 0; comp < _cc.ncomp(); ++comp)  // Iterate over components
        for (int col = 0; col < _cc.ncol(); ++col)  // Iterate over column cells
            _cc.colC<double>(idc, col, comp) = 1;

    // Particles need special treatment
    double* idp;
    for (int pblk = 0; pblk < _cc.npblk(); ++pblk)  // Iterate over particle blocks
    {
        idp = _cc.offsetPar(NV_id, pblk);
        for (int par = 0; par < _cc.npar(); ++par)  // Iterate over particle cells
            for (int comp = 0; comp < _cc.ncomp(); ++comp)  // Iterate over components
            {
                // The mass transport equations in the particles are differential
                _cc.parC<double>(idp, par, comp) = 1;
                // The reaction equations in the particle can be differential or algebraic
                _cc.parQ<double>(idp, par, comp) = _psim.getAdsorptionModel().isDifferential(comp);
            }
    }

    // Boundary fluxes are already algebraic (0), need no further attention

    // Tell IDAS the differential variables
    IDASetId(_idaMemBlock, NV_id);

    // Destroy the id vector
    N_VDestroy_Serial(NV_id);

    // ========================================================================================
    log::emit<Debug1>() << CURRENT_FUNCTION << ": Differential/algebraic variables set" << log::endl;


    // ========================================================================================
    // Suppress algebraic variables from error testing
    // ========================================================================================
    IDASetSuppressAlg(_idaMemBlock, TRUE);
    // ========================================================================================



    // ========================================================================================
    // Call setup routine for chromatography model
    // ========================================================================================
    _psim.getChromatographyModel().specialSetup();
    // ========================================================================================



    // ========================================================================================
    // Collect all sensitive parameters and toggle sensitivity computation
    // ========================================================================================
    _sensModelParams = _psim.getSensModelParams();
    _sensInletParams = getSensInletParams();
    _wantSensitivities = (getNSensParams() > 0) ? true : false;
    // ========================================================================================



    // ========================================================================================
    // Allocate space for output storage vectors
    // ========================================================================================
    _storeSizeInc = 200; ///todo make that user selectable or defined by heuristics, first guess was 200
    _nVecInc = 1;

    if (_writeAtUserTimes)
    {
        // Allocate all needed space once
        _stateAllTimes.reserve(_solutionTimes.size() * _cc.neq());
        _sensAllTimes.reserve(getNSensParams() * _solutionTimes.size() * _cc.neq());
        log::emit<Debug1>() << CURRENT_FUNCTION << ": Solution output memory allocated" << log::endl;
    }
    else
        // Allocate space for the first _storeSizeInc timesteps
        incStorageVectorSize(true);
    // ========================================================================================

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


void TimeIntegrator::initializeSensitivities() throw (CadetException)
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    checkSufficientDirections();

    //==========================================================================
    // AD Jacobi Setup
    //==========================================================================
    // Initialize all actives and their AD values to 0.0
    for (int eq = 0; eq < _cc.neq(); ++eq)  // Loop over all equations
    {
        getYAd(eq)   = 0.0;
        getResAd(eq) = 0.0;

        for (int dir = 0; dir < active::getMaxDirections(); ++dir)  // Iterate over all directions
        {
            getYAd(eq)  .setADValue(dir, 0.0);
            getResAd(eq).setADValue(dir, 0.0);
        }
    }

    // Set AD seeds for Jacobian assembly
    _psim.getSchurSolver().getJacobianData().initAdForJac(getYAd(), getDiagDir());
    //==========================================================================
    log::emit<Debug1>() << CURRENT_FUNCTION << ": AD for Jac inititialized" << log::endl;


    //==========================================================================
    // Turn off solution of sensitivity systems
    //   (this will be overridden by a call to IDASensInit! see below... )
    //  In fact, this has only an effect, if at first a computation with sensitivities is performed and then
    //  (without clearing and reallocating internal memory by cs_free/cs_malloc)
    //  another computation without sensitivities is started.
    IDASensToggleOff(_idaMemBlock);
    //==========================================================================


    //==========================================================================
    // SENSITIVITY ANALYSIS SETUP
    //==========================================================================
    if (getNSensParams() > 0)
    {

        double* yS;
        double* ySDot;

        //==============================================
        // Initialize sensitivity vectors with 0.0
        //==============================================
        for (int dir = 0; dir < active::getMaxDirections(); ++dir)  // Iterate over all directions
        {
            yS    = NV_DATA_S(_NVp_yS[dir]);
            ySDot = NV_DATA_S(_NVp_ySDot[dir]);

            for (int eq = 0; eq < _cc.neq(); ++eq)
            {
                yS   [eq] = 0.0;
                ySDot[eq] = 0.0;
            }
        }
        //==============================================
        log::emit<Debug1>() << CURRENT_FUNCTION << ": AD for Sens inititialized" << log::endl;


        //==================================================================================================
        // Set initial solution sensitivity w.r.t. SMA/SAI lambda = 1 for salt concentration in bound phase
        //==================================================================================================
        vpc_it it1 = std::find(_sensModelParams.begin(), _sensModelParams.end(), ParamID(SMA_LAMBDA, -1));  // search for a sensitive SMA_LAMBDA
        vpc_it it2 = std::find(_sensModelParams.begin(), _sensModelParams.end(), ParamID(SAI_LAMBDA, -1));  // search for a sensitive SAI_LAMBDA
        vpc_it it = (it1 != _sensModelParams.end()) ? it1 : it2;

        if (it != _sensModelParams.end())                       // One lambda was set sensitive!
        {
            int posLambda = it - _sensModelParams.begin();      // Index of the current lambda in the _sensModelParams vector

            for (int pblk = 0; pblk < _cc.npblk(); ++pblk)      // Iterate over particle blocks
            {
                yS = _cc.offsetPar(_NVp_yS[posLambda], pblk);
                for (int par = 0; par < _cc.npar(); ++par)      // Iterate over particle cells
                    _cc.parQ<double>(yS, par, 0) = 1.0;
            }

            log::emit<Debug1>() << CURRENT_FUNCTION << ": SMA/SAI_LAMBDA correction" << log::endl;
        }
        //==================================================================================================


        //============================================================================
        // Set parameter AD values = 0.0
        //============================================================================
        //
        // all active parameters are already initialized with AD values = 0.0
        //
        //============================================================================


        //============================================================================
        // Activate sensitivity computation for ith parameter in ith direction
        //============================================================================
        for (vpc_it it = _sensModelParams.begin(); it < _sensModelParams.end(); ++it)
        {
            // note the skip of getJacAdDirs()! sensitivities start after Jacobian entries!
            if (_psim.getChromatographyModel().contains(*it))                                           ///todo not very nice - find a good solution!
                _psim.getChromatographyModel().getParam(*it).setADValue(getJacAdDirs() + (it - _sensModelParams.begin()), 1.0);
            if (_psim.getAdsorptionModel().contains(*it))
                _psim.getAdsorptionModel().getParam(*it).setADValue(getJacAdDirs() + (it - _sensModelParams.begin()), 1.0);
        }
        //============================================================================
        log::emit<Debug1>() << CURRENT_FUNCTION << ": AD for sensitive params activated" << log::endl;



        //============================================================================
        // Initialize IDA sensitivity computation
        //============================================================================
        IDASensInit(_idaMemBlock, getNSensParams(), IDA_STAGGERED,
                &cadet::TimeIntegrator::residualSensWrapper, _NVp_yS, _NVp_ySDot);
        //============================================================================


        //============================================================================
        // Set sensitivity integration tolerances
        //============================================================================
        _relTolS = _relTol;

        _absTolS = new double[getNSensParams()];

        // Collect absolute sensitivities tolerances for model parameters
        for (vpc_it it = _sensModelParams.begin(); it < _sensModelParams.end(); ++it)
        {
            if (_psim.getChromatographyModel().contains(*it))
                _absTolS[it - _sensModelParams.begin()] = _psim.getChromatographyModel().getParam(*it).getAbsTolS();
            if (_psim.getAdsorptionModel().contains(*it))
                _absTolS[it - _sensModelParams.begin()] = _psim.getAdsorptionModel().getParam(*it).getAbsTolS();
        }

        // Collect absolute sensitivities tolerances for inlet parameters
        for (vpc_it it = _sensInletParams.begin(); it < _sensInletParams.end(); ++it)
            _absTolS[getNSensModelParams() + (it - _sensInletParams.begin())] = _inletParamAbsTolS.at(std::get<1>(*it));

        IDASensSStolerances(_idaMemBlock, _relTolS, _absTolS);

        delete [] _absTolS;
        //============================================================================


        //============================================================================
        // Activate sensitivity error control
        //============================================================================
        IDASetSensErrCon(_idaMemBlock, TRUE);
        //============================================================================

    }
    //==========================================================================
    // END SENSITIVITY ANALYSIS SETUP
    //==========================================================================

    // Decrease the number of directions in active type to needed ones
    adtl::setNumDir(getJacAdDirs() + getNSensParams());

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}

///todo write function that changes particle discretization type
///todo write function that tests inputConcentration function on heart and kidneys



void TimeIntegrator::integrate() throw (CadetException)
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;
    _timerAll.start();

#ifdef VERIFY_ANALYTIC_JAC
    _useAnalyticJac = true;
#endif

    vdc_it it;

    if (_printConfig) printConfiguration();

    checkParameters();
    checkTimeConsistency();

    // Decide whether to use user specified solution output times (IDA_NORMAL)
    // or internal integrator steps (IDA_ONE_STEP)
    if (_writeAtUserTimes)
    {
        _itask = IDA_NORMAL;
        log::emit<Info>() << "Solution is written at " << _solutionTimes.size() << " user specified time points." << log::endl;
    }
    else
    {
        _itask = IDA_ONE_STEP;
        log::emit<Info>() << "Solution is written at integration time points." << log::endl;
    }

    log::emit<Info>() << "Starting integration ..."  << log::endl;

    _t = _sectionTimes.at(0);
    const double tEnd = _sectionTimes[_sectionTimes.size() - 1];
    while (_t < tEnd)
    {
        log::emit<Debug1>() << CURRENT_FUNCTION << ": Integration for sec " << _sec << " started" << log::endl;

        _sec = getSection(_t);
        _startTime = _sectionTimes.at(_sec);

        // Determine continuous time slice
        int skip = 0;
        for (size_t i = _sec; i < _sectionTimes.size() - 2; ++i)
        {
            if (!_sectionContinuity.at(i))
                break;

            ++skip;
        }

        _endTime = _sectionTimes.at(_sec + skip + 1);
        _t = _startTime;

//        _psim.getWenoScheme().setWenoEpsilons(_NV_y, _sectionTimes, _sec);  ///todo check this!
        _psim.getSchurSolver().resetSpgmrIters();

        // =========================================================================
        //    IDA initialization part 2
        // =========================================================================
        _stepSize = _initStepSize * (_endTime - _startTime);

        // IDAS Step 7.3: Set the initial step size
        IDASetInitStep(_idaMemBlock, _stepSize);

        // IDAS Step 7.4: Set the stop time
        IDASetStopTime(_idaMemBlock, _endTime);

        // IDAS Step 5.2: Re-initialization of the solver ///todo get to know why!
        IDAReInit(_idaMemBlock, _startTime, _NV_y, _NV_yDot);

        // =========================================================================
        log::emit<Debug1>() << CURRENT_FUNCTION << ": IDA init part 2 done" << log::endl;

        // =========================================================================
        //    Compute consistent initial values
        // =========================================================================
        _psim.getChromatographyModel().calcIC(_t);
        _psim.getChromatographyModel().calcICSens(_t);
        // =========================================================================
        log::emit<Debug1>() << CURRENT_FUNCTION << ": Consistent initial conditions computed" << log::endl;

        // Update Jacobian
        _psim.getChromatographyModel().sectionSetup(_sec, _t);

        // Inititalize the IDA solver flag
        _solverFlag = IDA_SUCCESS;

        if (_writeAtUserTimes)
        {
            // Write initial conditions only if desired by user
            if (_sec == 0 && _solutionTimes.front() == 0.0) writeSolution();

            // Initialize iterator and forward it to the first solution time that lies inside the current section
            it = _solutionTimes.begin();
            while (*it <= _startTime) ++it;
            log::emit<Debug1>() << CURRENT_FUNCTION << ": Iterator forwarded to t = " << *it << log::endl;
        }
        else
        {
            // Always write initial conditions if solutions are written at integration times
            if (_sec == 0) writeSolution();

            // Here _tout - only during the first call to IDASolve - specifies the direction
            // and rough scale of the independent variable, see IDAS Guide p.33
            _tout = _endTime;
        }

        //  A small class that prints out progress of the integration
        if (_printProgress) log::emit<Info>() << "Section " << _sec + 1 << ":   ";
        SimpleProgress sp(_startTime, _endTime, 10);

        // Main loop which integrates the system until reaching the end time of the current section
        // or until an error occures
        while ((_solverFlag == IDA_SUCCESS) || (_solverFlag == IDA_ROOT_RETURN))
        {
            // Update _tout if we write solutions at user specified times
            if (_writeAtUserTimes)
                // Check if user specified times are sufficiently long.
                // otherwise integrate till IDA_TSTOP_RETURN
                _tout = (it != _solutionTimes.end()) ? *it : _sectionTimes.back() + 1;
            log::emit<Debug2>() << CURRENT_FUNCTION << ": _tout = " << _tout << log::endl;

            // IDA Step 11: Advance solution in time
            _solverFlag = IDASolve(_idaMemBlock, _tout, &_t, _NV_y, _NV_yDot, _itask);
            log::emit<Debug2>() << CURRENT_FUNCTION << ": _solverFlag = " << _solverFlag << log::endl;

            // Extract sensitivity information from IDA
            if (_wantSensitivities) IDAGetSens(_idaMemBlock, &_t, _NVp_yS);

            switch (_solverFlag)
            {
            case IDA_SUCCESS:           // _tout was reached
                if (_printProgress) sp.printIfHitNext(_t);  // Print progress
                writeSolution();
                ++it;
                log::emit<Debug2>() << CURRENT_FUNCTION << ": IDA_SUCCESS at t = " << _t << " advancing to " << (*it) << log::endl;
                break;
            case IDA_ROOT_RETURN:       // A root was found
                // Eventually call some routine
                log::emit<Debug2>() << CURRENT_FUNCTION << ": IDA_ROOT_RETURN at t = " << _t << log::endl;
                break;
            case IDA_TSTOP_RETURN:      // Section end time was reached
                if ((!_writeAtUserTimes) && (_endTime == _sectionTimes.back()))
                    // Write a solution for the ultimate endTime in the last section,
                    // when we write at integration times.
                    writeSolution();
//                if ((_writeAtUserTimes) && (_endTime == _sectionTimes.back()) && (_endTime != _solutionTimes.back()))
//                    // When writing at user specifyed times AND the end of the simulation was reached
//                    // AND the ultimate endTime was not in the user supplied times
//                    // we still might want to append the last solution. Then uncomment this section!
//                    writeSolution();
                log::emit<Debug2>() << CURRENT_FUNCTION << ": IDA_TSTOP_RETURN at t = " << _t << log::endl;
                break;
                // ... Exit while loop
            default:                    // An error occured
                log::emit<Debug2>() << CURRENT_FUNCTION << ": ERROR at t = " << _t << ": " << IDAGetReturnFlagName(_solverFlag) << log::endl;
                throw CadetException("Error in IDASolve!"); ///todo might not be necessary
                break;
            } // switch

        } // while

        if (_printProgress) log::emit() << log::endl;
        if (_printStatistics) printIntegratorStats();
    } // for (_sec ...)

    _timerAll.stop();
    if (_printTiming) printTimings();

#ifdef VERIFY_ANALYTICAL_JAC
    std::ostringstream oss;
    oss << scientific << setprecision(2) << _psim.getSchurSolver().getJacobianData().maxDiffC();
    log::emit<Info>() << Color::red << "Jacobian column   part max diff: " << oss.str() << Color::reset << log::endl;
    oss.str("");
    oss << scientific << setprecision(2) << _psim.getSchurSolver().getJacobianData().maxDiffP();
    log::emit<Info>() << Color::red << "Jacobian particle part max diff: " << oss.str() << Color::reset << log::endl;
#endif

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}



void TimeIntegrator::getSolutionTimes(std::vector<double>& userVector) const
{
    userVector.assign(_solutionTimes.begin(), _solutionTimes.end());
}


void TimeIntegrator::getAllSolutions(std::vector<double>& userVector) const
{
    userVector.assign(_stateAllTimes.begin(), _stateAllTimes.end());
}

void TimeIntegrator::getColSolutions(std::vector<double>& userVector) const
{
    userVector.clear();
    for (vdc_it it = _stateAllTimes.begin(); it < _stateAllTimes.end(); it += _cc.neq())
        userVector.insert(userVector.end(), it, it + _cc.neq_col());
}

void TimeIntegrator::getParSolutions(std::vector<double>& userVector) const
{
    userVector.clear();
    for (vdc_it it = _stateAllTimes.begin() + _cc.neq_col(); it < _stateAllTimes.end(); it += _cc.neq())
        userVector.insert(userVector.end(), it, it + (_cc.npblk() * _cc.neq_par()));
}

void TimeIntegrator::getBndSolutions(std::vector<double>& userVector) const
{
    userVector.clear();
    for (vdc_it it = _stateAllTimes.begin() + _cc.neq_col() + (_cc.npblk() * _cc.neq_par());
            it < _stateAllTimes.end(); it += _cc.neq())
        userVector.insert(userVector.end(), it, it + _cc.neq_bnd());
}


void TimeIntegrator::getAllSensitivities(std::vector<double>& userVector) const
{
    userVector.assign(_sensAllTimes.begin(), _sensAllTimes.end());
}

void TimeIntegrator::getColSensitivities(std::vector<double>& userVector) const
{
    userVector.clear();
    for (int sens = 0; sens < getNSensParams(); ++sens)
        for (vdc_it it = _sensAllTimes.begin() + sens * _cc.neq(); it < _sensAllTimes.end(); it += _cc.neq() * getNSensParams())
            userVector.insert(userVector.end(), it, it + _cc.neq_col());
}

void TimeIntegrator::getParSensitivities(std::vector<double>& userVector) const
{
    userVector.clear();
    for (int sens = 0; sens < getNSensParams(); ++sens)
        for (vdc_it it = _sensAllTimes.begin() + sens * _cc.neq() + _cc.neq_col();
                it < _sensAllTimes.end(); it += _cc.neq() * getNSensParams())
            userVector.insert(userVector.end(), it, it + (_cc.npblk() * _cc.neq_par()));
}

void TimeIntegrator::getBndSensitivities(std::vector<double>& userVector) const
{
    userVector.clear();
    for (int sens = 0; sens < getNSensParams(); ++sens)
        for (vdc_it it = _sensAllTimes.begin() + sens * _cc.neq() + _cc.neq_col() + (_cc.npblk() * _cc.neq_par());
                it < _sensAllTimes.end(); it += _cc.neq() * getNSensParams())
            userVector.insert(userVector.end(), it, it + _cc.neq_bnd());
}



void TimeIntegrator::getSolutionColumnOutlet(std::vector<double>& userVector, const int comp) const throw (CadetException)
{
    // Ensure 0 <= comp < ncomp
    if ((comp >= _cc.ncomp()) || (comp < 0)) throw CadetException("Component out of range!");

    // Copy column outlet entries
    userVector.clear();
    for (vdc_it it = _stateAllTimes.begin() + (comp + 1) * _cc.ncol() - 1; it < _stateAllTimes.end(); it += _cc.neq())
        userVector.push_back(*it);
}

void TimeIntegrator::getSolutionColumnInlet(std::vector<double>& userVector, const int comp) const throw (CadetException)
{
    // Ensure 0 <= comp < ncomp
    if ((comp >= _cc.ncomp()) || (comp < 0)) throw CadetException("Component out of range!");

    userVector.clear();
    size_t sec;
    std::vector<double> concInlet(_cc.ncomp(), 0.0);

    // Create dummy vectors for the inletConcentration function - we dont need derivatives here!
    int msip = getMaxSensInletParams();
    std::vector<bool> dummyInletParamIsSensitive = std::vector<bool>(msip, false);
    std::vector<std::vector<double> > dummyDeriv = std::vector<std::vector<double> >(_cc.ncomp(), std::vector<double>(msip, 0.0));

    // Call inlet function to get column inlet entries
    for(vdc_it it = _solutionTimes.begin(); it < _solutionTimes.end(); ++it)
    {
        sec = getSection(*it);

        inletConcentration(*it, sec, concInlet, dummyInletParamIsSensitive, dummyDeriv);
        userVector.push_back(concInlet.at(comp));
    }
}


void TimeIntegrator::getSensitivityColumnOutlet(std::vector<double>& userVector, const ParamID param, const int comp, const int sec) const throw (CadetException)
{
    // Ensure 0 <= comp < ncomp
    if ((comp >= _cc.ncomp()) || (comp < 0)) throw CadetException("Component out of range!");
    // Ensure -1 <= sec < nsec
    if ((sec >= _cc.nsec()) || (sec < -1)) throw CadetException("Section out of range!");

    userVector.clear();
    std::vector<ParamID> sensParams(_sensModelParams.begin(), _sensModelParams.end());
    sensParams.insert(sensParams.end(), _sensInletParams.begin(), _sensInletParams.end());

    ParamID pid = param;
    if (std::get<0>(param) == INLET_PARAMETER)
    {
        // Special treatment of inlet parameters since they use comp as global index and do not have a section
        pid = ParamID(std::get<0>(param), std::get<1>(param), -2);
    }

    vpc_it it1 = std::find(sensParams.begin(), sensParams.end(), pid);
    if (it1 == sensParams.end())
    {
        std::ostringstream ss;
        ss << "getSensitivityColumnOutlet(): Parameter " << e2s(std::get<0>(param)) << "[comp " << std::get<1>(param) << ", sec " << std::get<2>(param) << "] is not sensitive"
            << "; called with [comp " << comp << ", sec " << sec << "]";
        throw CadetException(ss.str());
    }

    int parIdx = it1 - sensParams.begin();
    // iterator.begin() + skip sensitivities + skip components
    for (vdc_it it2 = _sensAllTimes.begin() + parIdx * _cc.neq() + (comp + 1) * _cc.ncol() - 1;
            it2 < _sensAllTimes.end(); it2 += getNSensParams() * _cc.neq())
        userVector.push_back(*it2);
}


int TimeIntegrator::getSection(double t) const
{
    for (size_t i = 0; i < _sectionTimes.size() - 1; ++i)
    {
        if ((t >= _sectionTimes[i]) && (t < _sectionTimes[i+1]))
            return i;
    }

    if (t == _sectionTimes[_sectionTimes.size() - 1])
        return _sectionTimes.size() - 2;
    return -1;
}


template <typename T>
void TimeIntegrator::putStat(const std::string& str, const T& val, const std::string& delim) const
{
    std::ostringstream oss;
    oss << delim.at(0) << std::setw(35) << std::left << str << delim.at(1)
            << std::setw(10) << std::right << val << delim.at(2);
#ifndef BENCHMARK_MODE
    log::emit<Info>() << Color::green << oss.str().c_str() << Color::reset << log::endl;
#else
    std::cout << oss.str() << std::endl;
#endif
}


void TimeIntegrator::putTime(const std::string& str, double val, const std::string& delim) const
{
    std::ostringstream oss;
    oss << delim.at(0) << std::setw(35) << std::left << str << delim.at(1)
            << std::setw(10) << std::right << std::fixed << std::setprecision(3) << val << delim.at(2);
#ifndef BENCHMARK_MODE
    log::emit<Info>() << Color::green << oss.str().c_str() << Color::reset << log::endl;
#else
    std::cout << oss.str() << std::endl;
#endif
}


void TimeIntegrator::printTimings() const
{
    std::ostringstream oss;

    log::emit<Info>() << log::endl;
    putStat<std::string>("===================================","==========","|=|");
#ifdef _OPENMP
    putStat<std::string>("SIMULATOR TIMINGS","OpenMP");
#else
    putStat<std::string>("SIMULATOR TIMINGS","none");
#endif
    putStat<std::string>("===================================","==========","|=|");
    putStat<std::string>("> SERIAL CODE PARTS", "sec");
    putStat<std::string>("-----------------------------------","----------","|+|");
    putTime             ("TOTAL", this->timerAll());
    putStat<std::string>("","");
    putTime             ("Residual DAE", _psim.getChromatographyModel().timerResDae());
    if (_wantSensitivities) putTime("Residual Sens", _psim.getChromatographyModel().timerResSens());
    putTime             ("Factorization", _psim.getSchurSolver().timerFact());
    putTime             ("Solution", _psim.getSchurSolver().timerSol());
    putTime             ("Other", this->timerOther());
    putStat<std::string>("","");
    putTime             ("Solution/SPGMR", _psim.getSchurSolver().timerSolSpgmr());
    putTime             ("Solution/SPGMR/Axb", _psim.getSchurSolver().timerSolSpgmrAxb());
    putStat<std::string>("-----------------------------------","----------","|+|");
#ifdef _OPENMP
    oss << omp_get_max_threads() << " thread";
    if (omp_get_max_threads() > 1) oss << "s";
#else
    oss << "1 thread";
#endif
    putStat<std::string>("> PARALLEL CODE PARTS", oss.str());
    putStat<std::string>("-----------------------------------","----------","|+|");
    putTime             ("TOTAL/PAR", this->timerAllPar());
    putStat<std::string>("","");
    putTime             ("Residual DAE/PAR", _psim.getChromatographyModel().timerResDaePar());
    if (_wantSensitivities) putTime("Residual Sens/PAR", _psim.getChromatographyModel().timerResSensPar());
    putTime             ("Factorization/PAR", _psim.getSchurSolver().timerFactPar());
    putTime             ("Solution/PAR", _psim.getSchurSolver().timerSolPar());
    putTime             ("Solution/SPGMR/Axb/PAR", _psim.getSchurSolver().timerSolSpgmrAxbPar());
    putStat<std::string>("===================================","==========","|=|");

    log::emit<Info>() << log::endl;
}


void TimeIntegrator::printIntegratorStats() const
{
    long int nsteps;        // Cumulative number of internal steps
    long int nrevals;       // Cumulative number of calls to user's residual function
    long int nlinsetups;    // Cumulative number of calls to linear solver setup function
    long int netfails;      // Cumulative number of local error test failures that have occured

    int klast;              // Method order used during the last internal step
    int kcur;               // Method order to be used on the next internal step

    double hinused;         // Actual value of initial step size
    double hlast;           // Step size taken on the last internal step
    double hcur;            // Step size to be attempted on the next internal step
    double tcur;            // Current internal time reached by the solver

    IDAGetIntegratorStats(_idaMemBlock, &nsteps, &nrevals, &nlinsetups, &netfails,
            &klast, &kcur, &hinused, &hlast, &hcur, &tcur);

    std::ostringstream oss;
    oss << "SEC " << _sec + 1;

    log::emit<Info>() << log::endl;
    putStat<std::string>("===================================","==========","|=|");
    putStat<std::string>("INTEGRATOR STATISTICS", oss.str());
    putStat<std::string>("===================================","==========","|=|");
    putStat<long int>   ("Number of time steps", nsteps);
    putStat<long int>   ("Residual evaluations", nrevals);
    putStat<long int>   ("SPGMR iterations", _psim.getSchurSolver().getSpgmrIters());
//  putStat<long int>   ("Linear solver setups", nlinsetups);
    putStat<long int>   ("Error test failures", netfails);
    putStat<std::string>("","");
    putStat<int>        ("Last BDF order", klast);
//  putStat<int>        ("Next BDF order", kcur);
    putStat<std::string>("","");
    putStat<double>     ("Initial step size", hinused);
    putStat<double>     ("Last step size", hlast);
//  putStat<double>     ("Next step size", hcur);
    putStat<double>     ("Internally reached time", tcur);
    if (_wantSensitivities)
    {
        long int nresSevals;    // Number of calls to the sensitivity residual function
        long int nresevalsS;    // Number of calls to the user-supplied residual function for sensitivity
        long int nSetfails;     // Number of sensitivity local error test failures
        long int nlinsetupsS;   // Number of calls to the linear solver setup function for sensitivities

        IDAGetSensStats(_idaMemBlock, &nresSevals, &nresevalsS, &nSetfails, &nlinsetupsS);

        putStat<std::string>("","");
        putStat<long int>("Sens. residual evaluations", nresSevals);
        putStat<long int>("Sens. user residual evaluations", nresevalsS);
        putStat<long int>("Sens. error test failures", nSetfails);
//      putStat<long int>("Sens. linear solver setups", nlinsetupsS);
    }
    putStat<std::string>("===================================","==========","|=|");
    log::emit<Info>() << log::endl;

}


void TimeIntegrator::printConfiguration() const
{
    log::emit<Info>() << log::endl;
    putStat<std::string>("===================================","==========","|=|");
    putStat<std::string>("SIMULATOR CONFIGURATION", "");
    putStat<std::string>("===================================","==========","|=|");
    putStat<bool>       ("Print progress", _printProgress);
    putStat<bool>       ("Print statistics", _printStatistics);
    putStat<bool>       ("Print timing", _printTiming);
    putStat<bool>       ("Print parameter list", _printParamList);
    putStat<std::string>("-----------------------------------","----------","|+|");
    putStat<std::string>("> TIME INTEGRATOR", "");
    putStat<std::string>("-----------------------------------","----------","|+|");
    putStat<double>     ("Relative tolerance", _relTol);
    putStat<double>     ("Absolute tolerance", _absTol);
    putStat<std::string>("","");
    putStat<std::string>("Output times specified", (_writeAtUserTimes) ? "user" : "internal");
    putStat<double>     ("Initial ref. step size", _initStepSize);
    putStat<int>        ("Max. number of steps", _maxSteps);
    if (_wantSensitivities) putStat<int>("Number of sensitivities", getNSensParams());
    putStat<std::string>("-----------------------------------","----------","|+|");
    putStat<std::string>("> SCHUR SOLVER", "");
    putStat<std::string>("-----------------------------------","----------","|+|");
#ifdef VERIFY_ANALYTICAL_JAC
    putStat<std::string>("Jacobian type", "**VERIFY**");
#else
    putStat<std::string>("Jacobian type", (_useAnalyticJac) ? "analytic" : "AD");
#endif
    putStat<double>     ("Schur safety factor", _psim.getSchurSolver().getSchurSafety());
    putStat<int>        ("Max. no. of restarts", _psim.getSchurSolver().getMaxRestarts());
    putStat<int>        ("Max. size of Krylov subspace", _psim.getSchurSolver().getMaxKrylov());
    putStat<std::string>("Gram-Schmidt othogonal. type", (_psim.getSchurSolver().getGramSchmidtType()-1) ? "classic" : "modified");
    putStat<std::string>("-----------------------------------","----------","|+|");
    putStat<std::string>("> WENO SCHEME", "");
    putStat<std::string>("-----------------------------------","----------","|+|");
    putStat<int>        ("Weno order", _psim.getWenoScheme().getWenoOrder());
    putStat<double>     ("Weno epsilon", _psim.getWenoScheme().getWenoEps());
    putStat<int>        ("Boundary model", _psim.getWenoScheme().getBoundaryModel());
    putStat<std::string>("===================================","==========","|=|");
    log::emit<Info>() << log::endl;
}


// Error handler for IDA
void TimeIntegrator::idaErrorHandler(int error_code, const char* module,
        const char* function, char* msg, void* eh_data) throw (CadetException)
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    std::ostringstream oss;
    oss << "In function '" << function << "' of module '" << module
            << "', error code '" << IDAGetReturnFlagName(error_code) << "':" << std::endl << msg;

    if (error_code < 0) throw CadetException(oss.str());
    else log::emit<Except>() << oss.str() << log::endl;

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}



// Wrapper functions for IDA Callback routines
void TimeIntegrator::idaErrorHandlerWrapper(int error_code, const char* module,
        const char* function, char* msg, void* eh_data)
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    TimeIntegrator* ti = static_cast<TimeIntegrator*>(eh_data);
    try { ti->idaErrorHandler(error_code, module, function, msg, nullptr); }
    catch (const CadetException& e) { log::emit<Except>() << e.msg() << log::endl; }

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


int TimeIntegrator::residualDaeWrapper(double t, N_Vector NV_y, N_Vector NV_yDot,
        N_Vector NV_res, void* userData)
{
    SimulatorPImpl* sim = static_cast<SimulatorPImpl*>(userData);
    return sim->getChromatographyModel().residualDae(t, NV_y, NV_yDot, NV_res, nullptr);
}

int TimeIntegrator::residualSensWrapper(int ns, double t, N_Vector NV_y, N_Vector NV_yDot, N_Vector NV_res,
        N_Vector* NV_yS, N_Vector* NV_ySDot, N_Vector* NV_resS,
        void *userData, N_Vector NV_tmp1, N_Vector NV_tmp2, N_Vector NV_tmp3)
{
    SimulatorPImpl* sim = static_cast<SimulatorPImpl*>(userData);
    return sim->getChromatographyModel().residualSens(ns, t, NV_y, NV_yDot, NV_res, NV_yS, NV_ySDot, NV_resS,
            nullptr, NV_tmp1, NV_tmp2, NV_tmp3);
}

int TimeIntegrator::schurSolveWrapper(IDAMem IDA_mem, N_Vector NV_rhs, N_Vector weight,
        N_Vector NV_yCur, N_Vector NV_yDotCur, N_Vector NV_resCur)
{
    SimulatorPImpl* sim = static_cast<SimulatorPImpl*>(IDA_mem->ida_lmem);
    return sim->getSchurSolver().schurSolve(IDA_mem, NV_rhs, weight, NV_yCur, NV_yDotCur, NV_resCur);
}

int TimeIntegrator::schurComplementTimesVectorWrapper(void* userData, N_Vector NV_v, N_Vector NV_z)
{
    SimulatorPImpl* sim = static_cast<SimulatorPImpl*>(userData);
    return sim->getSchurSolver().schurComplementTimesVector(nullptr, NV_v, NV_z);
}



// Private member functions

void TimeIntegrator::checkSufficientDirections() const throw (CadetException)
{
    if (getJacAdDirs() + getNSensParams() > active::getMaxDirections())
    {
        std::ostringstream ss;
        ss << "FATAL ERROR: Not enough directions available for active data types!" << std::endl;
        ss << "Maximum active data type directions:       " << active::getMaxDirections() << std::endl;
        ss << "Directions required for Jacobian assembly: " << getJacAdDirs() << std::endl;
        ss << "Directions required for sensitivities:     " << getNSensParams() << std::endl;
        throw CadetException(ss.str());
    }
}


void TimeIntegrator::checkTimeConsistency() throw (CadetException)
{
    // Ensure that at least one section is defined
    if (_sectionTimes.size() < 2)
        throw CadetException("At least one section has to be specified!");

    // Ensure that all section start times are smaller than their end times
    for (vdc_it it = _sectionTimes.begin(); it < (_sectionTimes.end() - 1); ++it)
        if (*it >= *(it + 1))
            throw CadetException("The start time of each section must be smaller than its end time!");

    // Ensure that at least one solution output time is specified
    if (_writeAtUserTimes &&  !(_solutionTimes.size() > 0 ))
        throw CadetException("At least one solution time for writing output must be specified!");
}


void TimeIntegrator::checkParameters() const  throw (CadetException)
{
    if (_printParamList)
        log::emit<Info>() << log::endl << log::endl
        << Color::green << _psim.getChromatographyModel().info() << _psim.getAdsorptionModel().info()
        << Color::reset << log::endl;

    // Check all parameter values to be inside their predefined limits
    if (!_psim.getAdsorptionModel().checkAll() || !_psim.getChromatographyModel().checkAll())
        throw CadetException("At least one parameter is out of bounds");

    log::emit<Debug1>() << CURRENT_FUNCTION << ": All parameter checked to be inside bounds" << log::endl;
}


void TimeIntegrator::incStorageVectorSize(bool force)
{
    if ((_stateAllTimes.capacity() < _stateAllTimes.size() + _cc.neq()) || (force == true))
    {
        _stateAllTimes.reserve(_nVecInc * _storeSizeInc * _cc.neq());
        _sensAllTimes.reserve(_nVecInc * _storeSizeInc * getNSensParams() * _cc.neq());

        log::emit<Debug1>() << CURRENT_FUNCTION << ": Solution output memory incremented: " << _nVecInc << log::endl;
        ++_nVecInc;
    }
}


void TimeIntegrator::writeSolution()
{
    if (!_writeAtUserTimes)
    {
        // If needed, reallocate space for the storage vectors
        incStorageVectorSize();

        // Store the current time
        _solutionTimes.push_back(_t);
    }

    log::emit<Debug2>() << "Writing solution at time: " << _t << log::endl;

    // Append the current solution (current state state vector) to the storage vector
    double* y = NV_DATA_S(_NV_y);
    _stateAllTimes.insert(_stateAllTimes.end(), y, y + _cc.neq());

    // Append all the current sensitivities (current sensitivities vectors) to the storage vector
    for (int sens = 0; sens < getNSensParams(); ++sens)
    {
        double* yyS = NV_DATA_S(_NVp_yS[sens]);
        _sensAllTimes.insert(_sensAllTimes.end(), yyS, yyS + _cc.neq());
    }
    log::emit<Debug2>() << CURRENT_FUNCTION << ": Solution written to memory" << log::endl;
}


std::vector<ParamID> TimeIntegrator::getSensInletParams() const
{
    std::vector<ParamID> sensInletParams;
    for (int i = 0; i < getMaxSensInletParams(); ++i)
        if (_inletParamIsSensitive.at(i))
            sensInletParams.push_back(ParamID(INLET_PARAMETER, i, -2));
    return sensInletParams;
}


} // namespace cadet
