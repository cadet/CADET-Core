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

#ifndef TIMEINTEGRATOR_HPP_
#define TIMEINTEGRATOR_HPP_

#include <vector>

#include <idas/idas_impl.h>

#include "SimulatorPImpl.hpp"
#include "CadetException.hpp"
#include "CadetConvenience.hpp"
#include "ChromatographyModel.hpp"
#include "SchurSolver.hpp"
#include "Timer.hpp"
#include "active.hpp"


namespace cadet {

class TimeIntegrator {
public:

    // Constructor
    TimeIntegrator(SimulatorPImpl& sim);
    // Destructor
    ~TimeIntegrator();


    void setInitialConditions(const std::vector<double>& initC, const std::vector<double>& initQ);
    void setSectionTimes(const std::vector<double>& sectionTimes) throw (CadetException);
    void setSectionTimes(const std::vector<double>& sectionTimes, const std::vector<bool>& sectionContinuity) throw (CadetException);
    void setSolutionTimes(const std::vector<double>& solutionTimes) throw (CadetException);

    void initializeIntegrator();
    void initializeSensitivities() throw (CadetException);
    void integrate() throw (CadetException);

    void getSolutionTimes(std::vector<double>& userVector) const;

    void getAllSolutions(std::vector<double>& userVector) const;
    void getColSolutions(std::vector<double>& userVector) const;
    void getParSolutions(std::vector<double>& userVector) const;
    void getBndSolutions(std::vector<double>& userVector) const;

    void getAllSensitivities(std::vector<double>& userVector) const;
    void getColSensitivities(std::vector<double>& userVector) const;
    void getParSensitivities(std::vector<double>& userVector) const;
    void getBndSensitivities(std::vector<double>& userVector) const;


    /// \brief Returns the column outlet concentrations (chromatogram) for a given component
    void getSolutionColumnOutlet(std::vector<double>& userVector, const int comp) const throw (CadetException);
    /// \brief Returns the column inlet concentrations (boundary condition) for a given component
    void getSolutionColumnInlet(std::vector<double>& userVector, const int comp) const throw (CadetException);
    /// \brief Returns the sensitivity at the column outlet (sensogram) w.r.t. a given
    ///        parameter and component
    void getSensitivityColumnOutlet(std::vector<double>& userVector, const ParamID param, const int comp, const int sec = -1) const throw (CadetException);

    void printIntegratorStats() const;
    void printConfiguration() const;
    void printTimings() const;

    void idaErrorHandler(int error_code, const char* module,
            const char* function, char* msg, void* eh_data) throw (CadetException);

    // Static wrapper functions for IDAS callback functions
    /// \brief Wrapper callback function that redirects calls from IDAS to idaErrorHandler
    static void idaErrorHandlerWrapper(int error_code, const char* module,
            const char* function, char* msg, void* eh_data);
    /// \brief Wrapper callback function that redirects calls from IDAS to residualDae
    static int residualDaeWrapper(double t, N_Vector NV_y, N_Vector NV_yDot,
            N_Vector NV_res, void* userData);
    /// \brief Wrapper callback function that redirects calls from IDAS to residualSens
    static int residualSensWrapper(int ns, double t, N_Vector NV_y, N_Vector NV_yDot, N_Vector NV_res,
            N_Vector* NV_yS, N_Vector* NV_ySDot, N_Vector* NV_resS,
            void *userData, N_Vector NV_tmp1, N_Vector NV_tmp2, N_Vector NV_tmp3);
    /// \brief Wrapper callback function that redirects calls from IDAS to schurSolve
    static int schurSolveWrapper(IDAMem IDA_mem, N_Vector NV_rhs, N_Vector weight,
            N_Vector NV_yCur, N_Vector NV_yDotCur, N_Vector NV_resCur);
    /// \brief Wrapper callback function that redirects calls from IDAS to schurComplementTimesVector
    static int schurComplementTimesVectorWrapper(void* userData,
            N_Vector NV_v, N_Vector NV_z);


    // Inline getter and setter
    inline N_Vector getNvY()                const { return _NV_y; }
    inline N_Vector getNvYDot()             const { return _NV_yDot; }
    inline double*  getY()                  const { return NV_DATA_S(_NV_y); }
    inline double*  getYDot()               const { return NV_DATA_S(_NV_yDot); }

    inline N_Vector getNvYS(int param)      const { return _NVp_yS[param]; }
    inline N_Vector getNvYSDot(int param)   const { return _NVp_ySDot[param]; }
    inline double*  getYS(int param)        const { return NV_DATA_S(_NVp_yS[param]); }
    inline double*  getYSDot(int param)     const { return NV_DATA_S(_NVp_ySDot[param]); }

    inline N_Vector getNvTemp1()            const { return _NV_temp1; }
    inline N_Vector getNvTemp2()            const { return _NV_temp2; }
    inline double*  getTemp1()              const { return NV_DATA_S(_NV_temp1); }
    inline double*  getTemp2()              const { return NV_DATA_S(_NV_temp2); }

    inline active& getYAd(int pos)          const { return _yAd[pos]; }
    inline active& getResAd(int pos)        const { return _resAd[pos]; }

    inline active* getYAd()                 const { return _yAd; }
    inline active* getResAd()               const { return _resAd; }

    inline int getNSensParams()             const { return _sensModelParams.size() + _sensInletParams.size(); }
    inline int getNSensInletParams()        const { return _sensInletParams.size(); }
    inline int getNSensModelParams()        const { return _sensModelParams.size(); }
    inline int getCurrentSection()          const { return _sec; }
    int getSection(double t)                const;

    inline bool useAnalyticJacobian()       const { return _useAnalyticJac; }
    inline void useAnalyticJacobian(const bool useAnalyticJac) { _useAnalyticJac = useAnalyticJac; }

    inline bool factorizeJac()              const { return _factorizeJac; }
    inline void setFactorizeJac(bool factorizeJac = true) { _factorizeJac = factorizeJac; }

    int getJacAdDirs() const;
    int getDiagDir()   const;

    void inletConcentration(double t, int sec, std::vector<double>& inletConc,
            const std::vector<bool>& inletParamIsSensitive,
            std::vector<std::vector<double> >& dInletConc_dp) const
    {
        _inletBase->inletConcentration(t, sec, inletConc, inletParamIsSensitive, dInletConc_dp);
    }

    void inletConcentration(double t, std::vector<double>& inletConc,
            const std::vector<bool>& inletParamIsSensitive,
            std::vector<std::vector<double> >& dInletConc_dp) const
    {
        _inletBase->inletConcentration(t, getSection(t), inletConc, inletParamIsSensitive, dInletConc_dp);
    }

    inline void setInletProfile(InletBase* inletBase) { _inletBase = inletBase; }

    inline void setInletParamIsSensitive(const int index, double absTolS)
    {
        _inletParamIsSensitive.at(index) = true;
        _inletParamAbsTolS.at(index) = absTolS;
    }
    inline bool getInletParamIsSensitive(const int index)      const { return _inletParamIsSensitive.at(index); }
    inline const std::vector<bool>& getInletParamIsSensitive() const { return _inletParamIsSensitive; }

    inline void setMaxSensInletParams(const int nInPar)
    {
        _inletParamIsSensitive  = std::vector<bool>(nInPar, false);
        _inletParamAbsTolS      = std::vector<double>(nInPar, 0.0);
    }
    inline int  getMaxSensInletParams()                        const { return _inletParamIsSensitive.size(); }

    inline void resetSensInletParams()                               { setMaxSensInletParams(getMaxSensInletParams()); }

    // Setter for parameters set by user
    inline void configure(double relTol = 1e-12, double absTol = 1e-9,
            double initStepSize = 1e-6, int maxSteps = 10000)
    {
        setRelTol(relTol);
        setAbsTol(absTol);
        setInitStepSize(initStepSize);
        setMaxSteps(maxSteps);
        setFactorizeJac();
    }

    inline void configurePrinting(bool printProgress = false, bool printStatistics = false,
            bool printTiming = false, bool printParamList = false, bool printConfig = false)
    {
        setPrintProgress(printProgress);
        setPrintStatistics(printStatistics);
        setPrintTiming(printTiming);
        setPrintParamList(printParamList);
        setPrintConfig(printConfig);
    }

    inline void setInitStepSize(double initStepSize)    { _initStepSize = initStepSize; }
    inline void setMaxSteps(int maxSteps)               { _maxSteps = maxSteps; }
    inline void setAbsTol(double absTol)                { _absTol = absTol; }
    inline void setRelTol(double relTol)                { _relTol = relTol; }

    inline void setPrintProgress(const bool decision)   { _printProgress = decision; }
    inline void setPrintStatistics(const bool decision) { _printStatistics = decision; }
    inline void setPrintTiming(const bool decision)     { _printTiming = decision; }
    inline void setPrintParamList(const bool decision)  { _printParamList = decision; }
    inline void setPrintConfig(const bool decision)     { _printConfig = decision; }
    inline void setOpenMPThreads(const int threads)     { _nThreads = threads; }

    inline double timerAll()   const { return _timerAll.getTime(); }
    inline double timerOther() const {
        return _timerAll.getTime()
                - _psim.getChromatographyModel().timerResDae()
                - _psim.getChromatographyModel().timerResSens()
                - _psim.getSchurSolver().timerFact()
                - _psim.getSchurSolver().timerSol(); }
    inline double timerAllPar() const {
        return _psim.getChromatographyModel().timerResDaePar()
                + _psim.getChromatographyModel().timerResSensPar()
                + _psim.getSchurSolver().timerSolPar()
                + _psim.getSchurSolver().timerSolSpgmrAxbPar(); }

private:

    typedef std::vector<ParamID>::const_iterator    vpc_it;  //!< (v)ector (p)aramID (c)onst_(it)erator
    typedef std::vector<double>::const_iterator     vdc_it;  //!< (v)ector (d)ouble  (c)onst_(it)erator

    // Parameters set by user
    int                     _nThreads;
    int                     _maxSteps;
    double                  _initStepSize;
    double                  _absTol;
    double                  _relTol;

    // Booleans configuring the printing behaviour
    bool                    _printProgress;
    bool                    _printStatistics;
    bool                    _printTiming;
    bool                    _printParamList;
    bool                    _printConfig;

    /// \brief  Contains the start time of the first section, the start/end times of all intermediate sections
    ///         and the end time of the last section, (size = nsec + 1)
    std::vector<double>     _sectionTimes;
    /// \brief  Determines whether the transition from section i to section i+1 is continuous. The solver will
    //          be reset only at discontinuous transitions. The i-th element corresponds to the transition from
    //          _sectionTimes[i+1] to _sectionTimes[i+2]. Therefore size = nsec - 1.
    std::vector<bool>       _sectionContinuity;

    bool                    _writeAtUserTimes;      //!< false: use internal interator time steps for writing output, true: use user specified times
    std::vector<double>     _solutionTimes;         //!< Contains the user specified times for writing solutions to the output

    int                     _nVecInc;               //!< The number of vector increments that were already performed
    /// \brief  The increment (no. of time steps), by which the storage vectors sizes should be increased
    ///         when internal integrator time steps are used for output
    int                     _storeSizeInc;
    std::vector<double>     _stateAllTimes;         //!< Stores state vector at all output times
    std::vector<double>     _sensAllTimes;          //!< Stores all sensitivity vectors at all output times

    void*                   _idaMemBlock;           //!< IDAS memory block

    // Integration data
    double                  _startTime;
    double                  _endTime;
    double                  _t;                     //!< current time
    double                  _tout;
    int                     _sec;                   //!< current section

    double                  _stepSize;
    int                     _itask;                 //!< Determines whether to use user specified solution output times (IDA_NORMAL) or internal integrator steps (IDA_ONE_STEP)

    int                     _solverFlag;            //!< Reports the stop reason of the integrator after a call to IDASolve()


    N_Vector                _NV_y;                  //!< state vector
    N_Vector                _NV_yDot;               //!< state vector time derivative

    // Sensitivity data
    N_Vector*               _NVp_yS;                //!< sensitivities vector
    N_Vector*               _NVp_ySDot;             //!< sensitivities vector time derivative

    std::vector<ParamID>    _sensModelParams;
    std::vector<ParamID>    _sensInletParams;       //!< No of sensitive parameters in inletConcentration
    std::vector<bool>       _inletParamIsSensitive; //!< Stores for which of all possible inlet parameter sensitivities should be computed

    bool                    _wantSensitivities;     //!< true: computation of at least one sensitivity is activated, false: otherwise
    double                  _relTolS;
    double*                 _absTolS;
    std::vector<double>     _inletParamAbsTolS;     //!< Stores the absolute tolerances for the sensitivity computation of inlet parameters

    // Temporary storage
    N_Vector                _NV_temp1;
    N_Vector                _NV_temp2;

    // AD Jacobian stuff
    active*                 _yAd;                   //!< state vector  ///todo might be implemented as std::vector?
    active*                 _resAd;                 //!< residual vector

    bool                    _useAnalyticJac;        //!< Compute jacobians analyticaly when set true, by AD otherwise
    bool                    _factorizeJac;          //!< Determines whether a new factorization of the jacobian blocks is needed in the next call of the schurSolve routine or not

    // Function pointer to inlet concentration function
    InletBase*              _inletBase;

    SimulatorPImpl&         _psim;

    const CadetConstants&   _cc;

    OmpTimer                _timerAll;              //!< OpenMP timer for overall timing

    // Private member functions
    void checkSufficientDirections() const throw (CadetException);
    void checkTimeConsistency() throw (CadetException);
    void checkParameters() const throw (CadetException);
    void incStorageVectorSize(bool force = false);
    void writeSolution();

    std::vector<ParamID> getSensInletParams() const;

    template <typename T>
    void putStat(const std::string& str, const T& val,     const std::string& delim = "|||") const;
    void putTime(const std::string& str, double val, const std::string& delim = "|||") const;

};

} // namespace cadet


#endif // TIMEINTEGRATOR_HPP_
