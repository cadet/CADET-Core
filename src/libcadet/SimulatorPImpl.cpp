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

#include <string>
#include <sstream>
#include <vector>

#include "SimulatorPImpl.hpp"

#include "CadetConvenience.hpp"
#include "ParticleDiscretization.hpp"
#include "WenoScheme.hpp"
#include "SchurSolver.hpp"
#include "JacobianData.hpp"
#include "AdsorptionModel.hpp"
#include "ChromatographyModel.hpp"
#include "TimeIntegrator.hpp"
#include "Timer.hpp"

#include "CadetException.hpp"
#include "AdsorptionModel_MCL.hpp"
#include "AdsorptionModel_MPM.hpp"
#include "AdsorptionModel_SMA.hpp"
#include "AdsorptionModel_SAI.hpp"
#include "AdsorptionModel_EXTL.hpp"
#include "GeneralRateModel.hpp"
#include "active.hpp"

ACTIVE_INIT;

namespace cadet
{

// ================================================================================================
// Public Simulator interface implementation details
// ================================================================================================

Simulator::Simulator(int ncomp, int ncol, int npar, AdsorptionType ads_type, ChromatographyType chromType) :
    _sim(new SimulatorPImpl(ncomp, ncol, npar, ads_type, chromType))
{
}

Simulator::~Simulator()
{
    delete _sim;
}

void Simulator::setParameterValue(double value, const ParameterName id, int comp)
{
    _sim->setParameterValue(value, id, comp);
}

double Simulator::getParameterValue(const ParameterName id, int comp) const
{
    return _sim->getParameterValue(id, comp);
}

void Simulator::setSensitiveParameter(const ParameterName id, double absTolS, int comp)
{
    _sim->setSensitiveParameter(id, absTolS, comp);
}

void Simulator::setMaxSensInletParams(const int nInPar)
{
    _sim->getTimeIntegrator().setMaxSensInletParams(nInPar);
}

void Simulator::resetSensParams()
{
    _sim->resetSensParams();
}

std::vector<std::string> Simulator::getSensModelParamNames() const
{
    return _sim->getSensModelParamNames();
}

std::vector<int> Simulator::getSensModelParamComps() const
{
    return _sim->getSensModelParamComps();
}

int Simulator::getNSensParams() const
{
    return getNSensModelParams() + getNSensInletParams();
}

int Simulator::getNSensModelParams() const
{
    return _sim->getNSensModelParams();
}

int Simulator::getNSensInletParams() const
{
    return _sim->getNSensInletParams();
}

void Simulator::setPartDiscSchemeUser(const std::vector<double>& cellInterfaces)
{
    _sim->setPartDiscSchemeUser(cellInterfaces);
}

void Simulator::setPartDiscSchemeEqDist()
{
    _sim->setPartDiscSchemeEqDist();
}

void Simulator::setPartDiscSchemeEqVol()
{
    _sim->setPartDiscSchemeEqVol();
}

std::vector<double> Simulator::getParCellCoords() const
{
    return _sim->getParCellCoords();
}

void Simulator::setInletProfile(InletBase* inlet)
{
    _sim->getTimeIntegrator().setInletProfile(inlet);
}

void Simulator::setExternalProfile(ExternalBase* externalProfile)
{
    _sim->getAdsorptionModel().setExternalProfile(externalProfile);
}



void Simulator::setKineticAdsorptionModel(bool isKinetic)
{
    _sim->getAdsorptionModel().setIsKinetic(isKinetic);
}


void Simulator::setSolutionTimes(const std::vector<double>& solutionsTimes)
{
    _sim->getTimeIntegrator().setSolutionTimes(solutionsTimes);
}


void Simulator::setSectionTimes(const std::vector<double>& times)
{
    _sim->getTimeIntegrator().setSectionTimes(times);
}

void Simulator::setSectionTimes(const std::vector<double>& times, const std::vector<bool>& continuity)
{
    _sim->getTimeIntegrator().setSectionTimes(times, continuity);
}

void Simulator::initialize(const std::vector<double>& initC, const std::vector<double>& initQ)
{
    _sim->initialize(initC, initQ);
}

void Simulator::integrate()
{
    _sim->integrate();
}

void Simulator::useAnalyticJacobian(const bool analyticJac)
{
    _sim->getTimeIntegrator().useAnalyticJacobian(analyticJac);
}



void Simulator::getSolutionTimes(std::vector<double>& userVector) const
{
    _sim->getTimeIntegrator().getSolutionTimes(userVector);
}


void Simulator::getAllSolutions(std::vector<double>& userVector) const
{
    _sim->getTimeIntegrator().getAllSolutions(userVector);
}

void Simulator::getColSolutions(std::vector<double>& userVector) const
{
    _sim->getTimeIntegrator().getColSolutions(userVector);
}

void Simulator::getParSolutions(std::vector<double>& userVector) const
{
    _sim->getTimeIntegrator().getParSolutions(userVector);
}

void Simulator::getBndSolutions(std::vector<double>& userVector) const
{
    _sim->getTimeIntegrator().getBndSolutions(userVector);
}


void Simulator::getAllSensitivities(std::vector<double>& userVector) const
{
    _sim->getTimeIntegrator().getAllSensitivities(userVector);
}

void Simulator::getColSensitivities(std::vector<double>& userVector) const
{
    _sim->getTimeIntegrator().getColSensitivities(userVector);
}

void Simulator::getParSensitivities(std::vector<double>& userVector) const
{
    _sim->getTimeIntegrator().getParSensitivities(userVector);
}

void Simulator::getBndSensitivities(std::vector<double>& userVector) const
{
    _sim->getTimeIntegrator().getBndSensitivities(userVector);
}


void Simulator::getSolutionColumnOutlet(std::vector<double>& userVector, const int comp) const
{
    _sim->getTimeIntegrator().getSolutionColumnOutlet(userVector, comp);
}

void Simulator::getSolutionColumnInlet(std::vector<double>& userVector, const int comp) const
{
    _sim->getTimeIntegrator().getSolutionColumnInlet(userVector, comp);
}

void Simulator::getSensitivityColumnOutlet(std::vector<double>& userVector, const int comp, const ParamID& param) const
{
    _sim->getTimeIntegrator().getSensitivityColumnOutlet(userVector, param, comp);
}
void Simulator::getSensitivityColumnOutlet(std::vector<double>& userVector, const int comp, const ParameterName id, const int paramComp) const
{
    getSensitivityColumnOutlet(userVector, comp, ParamID(id, paramComp));
}




void Simulator::configureSchurSolver(double schurSafety, int maxRestarts, int maxKrylov, int gramSchmidtType)
{
    _sim->getSchurSolver().configure(schurSafety, maxRestarts, maxKrylov, gramSchmidtType);
}

void Simulator::configureWenoScheme(int wenoOrder, double wenoEps, int boundaryModel)
{
    _sim->getWenoScheme().configure(wenoOrder, wenoEps, boundaryModel);
}

void Simulator::configureTimeIntegrator(double relTol, double absTol, double initStepSize, int maxSteps)
{
    _sim->getTimeIntegrator().configure(relTol, absTol, initStepSize, maxSteps);
}

void Simulator::configurePrinting(bool printProgress, bool printStatistics, bool printTiming, bool printParamList, bool printConfig)
{
    _sim->getTimeIntegrator().configurePrinting(printProgress, printStatistics, printTiming, printParamList, printConfig);
}

// ================================================================================================




// ================================================================================================
// SimulatorPImpl implementation details
// ================================================================================================

SimulatorPImpl::SimulatorPImpl(int ncomp, int ncol, int npar, AdsorptionType adsType, ChromatographyType chromType) throw (CadetException)
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    // Set the default set of constants
    _cc = new CadetConstants(ncomp, ncol, npar);

    // Create a new adsorption model, according to given adsorption type
    switch (adsType)
    {
    case MULTI_COMPONENT_LANGMUIR:
        _adsModel = new AdsorptionModel_MCL(*this);
        break;
    case MOBILE_PHASE_MODULATORS:
        _adsModel = new AdsorptionModel_MPM(*this);
        break;
    case STERIC_MASS_ACTION:
        _adsModel = new AdsorptionModel_SMA(*this);
        break;
    case SELF_ASSOCIATION:
        _adsModel = new AdsorptionModel_SAI(*this);
        break;
    case EXTERNAL_LANGMUIR:
        _adsModel = new AdsorptionModel_EXTL(*this);
        break;
    default:
        std::ostringstream ss;
        ss << "SimulatorPImpl::SimulatorPImpl(): Unknown isotherm model specified: " << e2s(adsType);
        throw CadetException(ss.str());
        break;
    }

    _parDisc        = new ParticleDiscretization(npar);
    _wenoScheme     = new WenoScheme(*this);
    _schurSolver    = new SchurSolver(*this);
    _timeIntegrator = new TimeIntegrator(*this);

    switch (chromType)
    {
    case GENERAL_RATE_MODEL:
        _chromModel = new GeneralRateModel(*this);
        break;
    case LUMPED_RATE_MODEL:
//        _chromModel = new LumpedRateModel(*this);
        break;
    default:
        std::ostringstream ss;
        ss << "SimulatorPImpl::SimulatorPImpl(): Unknown chromatography model specified: " << e2s(chromType);
        throw CadetException(ss.str());
        break;
    }

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}

SimulatorPImpl::~SimulatorPImpl()
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    delete _cc;
    delete _schurSolver;
    delete _wenoScheme;
    delete _timeIntegrator;
    delete _parDisc;
    delete _chromModel;
    delete _adsModel;

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}

void SimulatorPImpl::setParameterValue(double value, const ParameterName id, int comp) throw (CadetException)
{
    if (_chromModel->contains(id, comp))
    {
        _chromModel->setValue(value, id, comp);
    }
    else if (_adsModel->contains(id, comp))
    {
        _adsModel->setValue(value, id, comp);
    }
    else
    {
        std::ostringstream ss;
        ss << "SimulatorPImpl::setParameterValue(): Parameter does not exist: " << e2s(id) << "[" << comp << "]";
        throw CadetException(ss.str());
    }
}

double SimulatorPImpl::getParameterValue(const ParameterName id, int comp) const throw (CadetException)
{
    if (_chromModel->contains(id, comp))
    {
        return _chromModel->getValue<double>(id, comp);
    }
    else if (_adsModel->contains(id, comp))
    {
        return _adsModel->getValue<double>(id, comp);
    }
    else
    {
        std::ostringstream ss;
        ss << "SimulatorPImpl::getParameterValue(): Parameter does not exist: " << e2s(id) << "[" << comp << "]";
        throw CadetException(ss.str());
    }
}



// the copied from SimulatorPImpl.hpp --> needs to be implemented...


// activate sensitivity computation for parameter identified by pname
void SimulatorPImpl::setSensitiveParameter(const ParameterName id, double absTolS, int comp) throw (CadetException)
{
    if (_chromModel->contains(id, comp))
    {
        _chromModel->setSensitive(id, absTolS, comp);
    }
    else if (_adsModel->contains(id, comp))
    {
        _adsModel->setSensitive(id, absTolS, comp);
    }
    else if (id == INLET_PARAMETER)
    {
        _timeIntegrator->setInletParamIsSensitive(comp, absTolS); // comp is here the global index of the sensitive inlet param
    }
    else
    {
        std::ostringstream ss;
        ss << "SimulatorPImpl::setSensitiveParameter(): Parameter does not exist: " << e2s(id) << "[" << comp << "]";
        throw CadetException(ss.str());
    }
}


// reset to no sensitivity computation
void SimulatorPImpl::resetSensParams()
{
    _chromModel->resetSensParams();
    _adsModel->resetSensParams();
    _timeIntegrator->resetSensInletParams();
}


// get the ParamID of all sensitive parameters
std::vector<ParamID> SimulatorPImpl::getSensModelParams() const
{
    std::vector<ParamID> sensModelParams = _chromModel->getSensParams();
    std::vector<ParamID> tmpVec = _adsModel->getSensParams();
    sensModelParams.insert(sensModelParams.end(), tmpVec.begin(), tmpVec.end());
    return sensModelParams;
}


// get the names of all sensitive parameters
std::vector<std::string> SimulatorPImpl::getSensModelParamNames() const
{
    std::vector<std::string> sensParamNames = _chromModel->getSensParamNames();
    std::vector<std::string> tmpVec = _adsModel->getSensParamNames();
    sensParamNames.insert(sensParamNames.end(), tmpVec.begin(), tmpVec.end());
    return sensParamNames;
}

// get the components of all sensitive parameters
std::vector<int> SimulatorPImpl::getSensModelParamComps() const
{
    std::vector<int> sensParamComps = _chromModel->getSensParamComps();
    std::vector<int> tmpVec = _adsModel->getSensParamComps();
    sensParamComps.insert(sensParamComps.end(), tmpVec.begin(), tmpVec.end());
    return sensParamComps;
}

// get number of activated sensitivities
int SimulatorPImpl::getNSensModelParams() const
{
    return _timeIntegrator->getNSensModelParams();
}

int SimulatorPImpl::getNSensInletParams() const
{
    return _timeIntegrator->getNSensInletParams();
}

// set the particle discretization scheme
void SimulatorPImpl::setPartDiscSchemeUser(const std::vector<double>& cellInterfaces)
{
    _parDisc->setUserdefinedRadialDisc(cellInterfaces);
}

void SimulatorPImpl::setPartDiscSchemeEqDist()
{
    _parDisc->setEquidistantRadialDisc();
}

void SimulatorPImpl::setPartDiscSchemeEqVol()
{
    _parDisc->setEquivolumeRadialDisc();
}

//// get the (dimensionless!) cell centers of column/particle discretization cells
//// dimensionless values need to be multiplied by L or R, respectively, to get absolute values
//std::vector<double> SimulatorPImpl::getColCellCoords()
//{
//    return _colDisc->getColCellCoords();
//}

std::vector<double> SimulatorPImpl::getParCellCoords() const
{
    return _parDisc->getParCellCoords();
}

//// set the column inlet profile function (user provided!)
//// the function inletConcentration() takes the time as an input
//// and returns the concentrations of all ncomp chemical species in c[0..ncomp-1].
//// if sensitivities w.r.t. inflow profile parameters are desired, the number of these parameters must be provided!
//// in this case, also the derivates of the different inlet concentrations w.r.t. to these parameters
//// must be provided by inletConcentration() in dcdp[NCOMP][nSensInflowParams]
//void setColumnInletProfile(void inletConcentration(double t, double *c, double **dcdp = NULL), int nSensInflowParams = 0);
//
//// check if all parameters have been assigned physically meaningful values
//// (might be useful before start of integration...)
//bool checkParameters();
//
//// print parameters to stdout
//void printParameters();
//
//// die alten cs_set_output routinen sind sehr unschoen mit den ganzen mehrfach pointern.
//// ich wuerde sie zuvor allerdings erstmal gerne so lassen. noch keine idee, wie man das schoener loesen kann...
//void setOutput( .... );
//
//// activate printing to stdout if b==true
//// deactivate printing to stdout if b==false
//void setPrintToStdout(bool b);
//
//// set time points of (potential) discontinuities
//// where an integrator restart with adjusted initial conditions might be required
//void setSectionTimes(const std::vector<double> const & const times);

// initialize all the IDA stuff, compute AD directions etc.
void SimulatorPImpl::initialize(const std::vector<double>& initC, const std::vector<double>& initQ)
{
    _timeIntegrator->setInitialConditions(initC, initQ);
    _timeIntegrator->initializeIntegrator();
    _timeIntegrator->initializeSensitivities();
}

// run the integration from tstart till tend
// what about the output time points? pass here (as before) or setup earlier?
void SimulatorPImpl::integrate()
{
    _timeIntegrator->integrate();
}

//// print timing statistics to stdout
//void printTimingStats();

} // namespace cadet
