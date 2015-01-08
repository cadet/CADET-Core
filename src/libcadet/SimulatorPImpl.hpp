// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2015: Eric von Lieres¹, Joel Andersson,
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

#ifndef SIMULATORPIMPL_HPP_
#define SIMULATORPIMPL_HPP_

#include <cstdint>

// Define CADET_VERSION
#ifndef CADET_VERSION
  #define CADET_VERSION_STRING "CADET NO-VERSION-DEFINED"
#else
  #define CADET_VERSION_STRING "CADET " CADET_VERSION
#endif

#ifdef MATLABMEX
    // The Intel MKL shipped with Matlab uses 64 bit integers
    typedef int64_t lapackInt_t;
#else
    typedef int lapackInt_t;
#endif

// Choose the right LAPACK library
#ifdef MATLABMEX
    // This library is built for MEX use, so take Intel MKL supplied by Matlab
    #define USESYSTEMLAPACK 1
#else
    #ifndef _WIN32
        // Linux or MacOS: Use system LAPACK
        #define USESYSTEMLAPACK 1
    #endif
#endif

// Use mylapack version on windows - no lumped rate model possible! ALSO SEE definitions of DGBTRF ...
#if !defined(_WIN32) || defined(USESYSTEMLAPACK)
    // lapack version declarations
    extern "C" void dgbtrf_(lapackInt_t *m, lapackInt_t *n, lapackInt_t *kl, lapackInt_t *ku, double *ab,
            lapackInt_t *ldab, lapackInt_t *ipiv, lapackInt_t *info);

    extern "C" void dgbtrs_(char *trans, lapackInt_t *n, lapackInt_t *kl, lapackInt_t * ku, lapackInt_t *nrhs,
            double *ab, lapackInt_t *ldab, lapackInt_t *ipiv, double *b, lapackInt_t *ldb, lapackInt_t *info);

    extern "C" void dgbmv_(char *trans, lapackInt_t *m, lapackInt_t *n, lapackInt_t *kl, lapackInt_t *ku,
            double *alpha, double *a, lapackInt_t *lda, double *x, lapackInt_t *incx, double *beta,
            double *y, lapackInt_t *incy);
#endif

// mylapack version declarations
void mydgbtrf(int *m, int *n, int *kl, int *ku, double *ab,
        int *ldab, int *ipiv, int *info);

void mydgbtrs(char *trans, int *n, int *kl, int * ku, int *nrhs,
        double *ab, int *ldab, int *ipiv, double *b, int *ldb, int *info);

void mydgbmv(char *trans, int *m, int *n, int *kl, int *ku,
        double *alpha, double *a, int *lda, double *x, int *incx, double *beta,
        double *y, int *incy);


#ifdef USESYSTEMLAPACK
    #define DGBTRF dgbtrf_
    #define DGBTRS dgbtrs_
    #define DGBMV  dgbmv_
#else
    // WARNING: This routines cannot handle NPAR = 1 !!!
    #define DGBTRF mydgbtrf
    #define DGBTRS mydgbtrs
    #define DGBMV  mydgbmv
#endif



#include <string>
#include <vector>
#include <tuple>

#include "Simulator.hpp"        // to have enumeration types available
#include "CadetException.hpp"   // to have exception types available


namespace cadet
{

// Forward declarations
class CadetConstants;
class AdsorptionModel;
class ChromatographyModel;
class TimeIntegrator;
class ParticleDiscretization;
class WenoScheme;
class SchurSolver;
class OmpTimer;

// Tuple of ParameterName, component, and section makes up for the identification type of a parameter
typedef std::tuple<ParameterName, int, int> ParamID;


// Private implementation of class Simulator
class SimulatorPImpl
{
public:

    SimulatorPImpl(int ncomp, int ncol, int npar, int nsec, AdsorptionType adsType, ChromatographyType chromType) throw (CadetException);
    ~SimulatorPImpl();

    void setParameterValue(double value, const ParameterName id, int comp = -1, int sec = -1) throw (CadetException);
    double getParameterValue(const ParameterName id, int comp = -1, int sec = -1) const throw (CadetException);

    void setSensitiveParameter(const ParameterName id, double absTolS = 1e-5, int comp = -1, int sec = -1) throw (CadetException);
    void setParameterSectionDependent(const ParameterName id, bool depends);

    void resetSensParams();

    std::vector<ParamID> getSensModelParams() const;
    std::vector<int> getSensModelParamSecs() const;
    std::vector<int> getSensModelParamComps() const;
    std::vector<std::string> getSensModelParamNames() const;

    int getNSensModelParams() const;
    int getNSensInletParams() const;

    void setPartDiscSchemeUser(const std::vector<double> & cellInterfaces);
    void setPartDiscSchemeEqDist();
    void setPartDiscSchemeEqVol();

    std::vector<double> getParCellCoords() const;

    void initialize(const std::vector<double>& initC, const std::vector<double>& initQ);
    void initialize(const std::vector<double>& initState);
    void initializeWithGivenSensitivities(const std::vector<double>& initState, const std::vector<double>& initSens);
    void integrate();

    // Member functions to access private members through a Simulator pointer/reference (e.g. from inside IDA)
    inline const CadetConstants&             getCadetConstants()         const { return *_cc; }

    inline const AdsorptionModel&            getAdsorptionModel()        const { return *_adsModel; }
    inline       AdsorptionModel&            getAdsorptionModel()              { return *_adsModel; }

    inline const ChromatographyModel&        getChromatographyModel()    const { return *_chromModel; }
    inline       ChromatographyModel&        getChromatographyModel()          { return *_chromModel; }

    inline const TimeIntegrator&             getTimeIntegrator()         const { return *_timeIntegrator; }
    inline       TimeIntegrator&             getTimeIntegrator()               { return *_timeIntegrator; }

    inline const ParticleDiscretization&     getParticleDiscretization() const { return *_parDisc; }

    inline const WenoScheme&                 getWenoScheme()             const { return *_wenoScheme; }
    inline       WenoScheme&                 getWenoScheme()                   { return *_wenoScheme; }

    inline const SchurSolver&                getSchurSolver()            const { return *_schurSolver; }
    inline       SchurSolver&                getSchurSolver()                  { return *_schurSolver; }

private:

    // Pointer to private member objects
    CadetConstants*             _cc;
    AdsorptionModel*            _adsModel;
    ChromatographyModel*        _chromModel;
    TimeIntegrator*             _timeIntegrator;

    ParticleDiscretization*     _parDisc;

    WenoScheme*                 _wenoScheme;
    SchurSolver*                _schurSolver;

};


} // namespace cadet

#endif /* SIMULATORPIMPL_HPP_ */
