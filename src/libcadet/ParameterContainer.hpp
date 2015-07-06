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

#ifndef PARAMETERCONTAINER_HPP_
#define PARAMETERCONTAINER_HPP_

#include <map>
#include <string>
#include <sstream>

#include "SimulatorPImpl.hpp"
#include "CadetConvenience.hpp"
#include "CadetException.hpp"
#include "CadetLogger.hpp"
#include "Parameter.hpp"

namespace cadet
{

/// \brief An abstract base class holding a hash map of pointers to parameters
template<typename ParamType>
class ParameterContainer
{
public:

    // Constructor
    ParameterContainer(const SimulatorPImpl& sim) : _sim(sim), _cc(sim.getCadetConstants()) { }

    // Destructor
    ~ParameterContainer() { }

    /// \brief Adds a new parameter to the container
    /// The parameter must have been created before
    void addParam(Parameter<ParamType>& param) { addParam(ParamID(param.getId(), param.getComp(), param.getSec()), &param); }

    /// \brief Adds a new parameter to the container
    /// The parameter must have been created before
    void addParam(const ParameterName id, int comp, int section, Parameter<ParamType>* const ptr) { addParam(ParamID(id, comp, section), ptr); }

    /// \brief Adds a new parameter to the container
    /// The parameter must have been created before
    void addParam(const ParamID& pid, Parameter<ParamType>* const ptr) { _params[pid] = ptr; }

    /// \brief Removes a parameter from the container
    void removeParam(const ParameterName id, int comp, int section) { removeParam(ParamID(id, comp, section)); }

    /// \brief Removes a parameter from the container
    void removeParam(const ParamID& pid) { _params.erase(pid); }

    /// \brief Removes a parameter from the container
    void removeParam(const Parameter<ParamType>& param) { _params.erase(ParamID(param.getId(), param.getComp(), param.getSec())); }

    /// \brief Removes all matching parameters from the container
    void removeParamsOfSection(const ParameterName id, int sec);

    /// \brief Returns a refernce to the parameter named 'pname' describing the component 'comp' in section 'sec'
    /// by default, comp and sec are -1 (parameter does not depend on any component and section)
    virtual const Parameter<ParamType>& getParam(const ParameterName id, int comp = -1, int sec = -1) const throw (CadetException);
    Parameter<ParamType>& getParam(const ParamID& param) throw (CadetException);

    /// \brief Set the value of the parameter named 'pname' describing the component 'comp' in section 'sec'
    /// by default, comp and sec are -1 (parameter does not depend on any component and section)
    void setValue(double value, const ParameterName id, int comp = -1, int sec = -1) throw (CadetException);

    /// \brief Get the value of the parameter named 'pname' describing the component 'comp' in section 'sec'
    /// by default, comp and sec are -1 (parameter does not depend on any component and section)
    template<typename ReturnType>
    const ReturnType getValue(const ParameterName id, int comp = -1, int sec = -1) const throw (CadetException);

    /// \brief Mark the parameter named 'pname' describing the component 'comp' as sensitive
    /// by default, comp and sec are -1 (parameter does not depend on any component and section)
    void setSensitive(const ParameterName id, double absTolS = 1e-5, int comp = -1, int sec = -1) throw (CadetException);

    /// \brief Mark all parameters as non sensitive
    void resetSensParams();

    /// \brief Get the ParamID of all parameters marked sensitive in the map
    const std::vector<ParamID> getSensParams() const;

    /// \brief Get the names of all parameters marked sensitive in the map
    const std::vector<std::string> getSensParamNames() const;

    /// \brief Get the id of all parameters marked sensitive in the map
    const std::vector<ParameterName> getSensParamIds() const;

    /// \brief Get the component numbers of all parameters marked sensitive in the map
    const std::vector<int> getSensParamComps() const;

    /// \brief Get the section numbers of all parameters marked sensitive in the map
    const std::vector<int> getSensParamSecs() const;

    /// \brief Get the number of parameters marked sensitive in the map
    int getNumSens() const;

    /// \brief Check that all parameter values are in between their predefined bounds
    // Returns: false if any wrong value
    //          true  if all parameters are OK
    bool checkAll() const;

    // A string containing info on all hold parameters
    const std::string info() const;

    /// \brief Check if the parameter 'pname' with component 'comp' and section 'sec' is contained in our map
    bool contains(const ParameterName id, int comp = -1, int sec = -1);
    bool contains(const ParamID& param);

protected:
    // Const and Non-const versions of map-iterators
    typedef typename std::map<ParamID, Parameter<ParamType>* >::iterator MapIt;
    typedef typename std::map<ParamID, Parameter<ParamType>* >::const_iterator ConstMapIt;

    // STL map to hold the Parameter objects
    std::map<ParamID, Parameter<ParamType>* > _params;

    const SimulatorPImpl& _sim;
    const CadetConstants& _cc;

}; // class ParameterContainer




// Overloaded operator<< for printing parameter info to streams
template <typename ParamType>
std::ostream& operator<<(std::ostream& out, const ParameterContainer<ParamType>& pc)
{
    out << pc.info();
    return out;
}



// Implementation details


template <typename ParamType>
void ParameterContainer<ParamType>::removeParamsOfSection(const ParameterName id, int sec)
{
    for (auto it = _params.cbegin(); it != _params.cend(); /* no increment */)
    {
        if ((std::get<0>(it->first) == id) && (std::get<2>(it->first) == sec))
            _params.erase(it++);
        else
            ++it;
    }    
}


template <typename ParamType>
const Parameter<ParamType>& ParameterContainer<ParamType>::getParam(const ParameterName id, int comp, int sec) const throw (CadetException)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    ConstMapIt it = _params.find(ParamID(id, comp, sec));

    if (it != _params.end())
    {
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
        return *(it->second);
    }
    else
    {
        std::ostringstream ss;
        ss << "getParam(): Parameter not existent in map: " << e2s(id) << "[comp " << comp << ", sec " << sec << "]";
        throw CadetException(ss.str());
    }
}


template <typename ParamType>
Parameter<ParamType>& ParameterContainer<ParamType>::getParam(const ParamID& param) throw (CadetException)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    MapIt it = _params.find(param);

    if (it != _params.end())
    {
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
        return *(it->second);
    }
    else
    {
        std::ostringstream ss;
        ss << "getParam(): Parameter not existent in map: " << e2s(std::get<0>(param)) << "[comp " << std::get<1>(param) << ", sec " << std::get<2>(param) << "]";
        throw CadetException(ss.str());
    }
}


template <typename ParamType>
void ParameterContainer<ParamType>::setValue(double value, const ParameterName id, int comp, int sec) throw (CadetException)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    MapIt it = _params.find(ParamID(id, comp, sec));

    if (it != _params.end())
        it->second->setValue(value);
    else
    {
        std::ostringstream ss;
        ss << "setValue(): Parameter not existent in map: " << e2s(id) << "[comp " << comp << ", sec" << sec << "]";
        throw CadetException(ss.str());
    }

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


template <typename ParamType> template<typename ReturnType>
const ReturnType ParameterContainer<ParamType>::getValue(const ParameterName id, int comp, int sec) const throw (CadetException)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    ConstMapIt it = _params.find(ParamID(id, comp, sec));

    if (it != _params.end())
    {
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
        return it->second->template getValue<ReturnType> ();
    }
    else
    {
        std::ostringstream ss;
        ss << "getValue(): Parameter not existent in map: " << e2s(id) << "[comp " << comp << ", sec " << sec << "]";
        throw CadetException(ss.str());
    }
}

template<typename ParamType>
void ParameterContainer<ParamType>::setSensitive(const ParameterName id, double absTolS, int comp, int sec) throw (CadetException)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    MapIt it = _params.find(ParamID(id, comp, sec));

    if (it != _params.end())
    {
        it->second->setSensitive();
        it->second->setAbsTolS(absTolS);
    }
    else
    {
        std::ostringstream ss;
        ss << "setSensitive(): Parameter not existent in map: " << e2s(id) << "[comp " << comp << ", sec " << sec << "]";
        throw CadetException(ss.str());
    }

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


template <typename ParamType>
void ParameterContainer<ParamType>::resetSensParams()
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    for (MapIt it = _params.begin(); it != _params.end(); ++it)
    {
        it->second->resetSensitive();
        it->second->setAbsTolS(1e-5);
    }

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


template <typename ParamType>
const std::vector<ParamID> ParameterContainer<ParamType>::getSensParams() const
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    std::vector<ParamID> sensParams;
    for (ConstMapIt it = _params.begin(); it != _params.end(); ++it)
    {
        if (it->second->isSensitive())
            sensParams.push_back(it->first);
    }

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    return sensParams;
}


template <typename ParamType>
const std::vector<std::string> ParameterContainer<ParamType>::getSensParamNames() const
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    std::vector<std::string> sensParamNames;
    for (ConstMapIt it = _params.begin(); it != _params.end(); ++it)
    {
        if (it->second->isSensitive())
            sensParamNames.push_back(it->second->getName());
    }

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    return sensParamNames;
}


template <typename ParamType>
const std::vector<ParameterName> ParameterContainer<ParamType>::getSensParamIds() const
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    std::vector<ParameterName> sensParamIds;
    for (ConstMapIt it = _params.begin(); it != _params.end(); ++it)
    {
        if (it->second->isSensitive())
            sensParamIds.push_back(it->second->getId());
    }

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    return sensParamIds;
}


template <typename ParamType>
const std::vector<int> ParameterContainer<ParamType>::getSensParamComps() const
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    std::vector<int> sensParamComps;
    for (ConstMapIt it = _params.begin(); it != _params.end(); ++it)
    {
        if (it->second->isSensitive())
            sensParamComps.push_back(it->second->getComp());
    }

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    return sensParamComps;
}


template <typename ParamType>
const std::vector<int> ParameterContainer<ParamType>::getSensParamSecs() const
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    std::vector<int> sensParamSecs;
    for (ConstMapIt it = _params.begin(); it != _params.end(); ++it)
    {
        if (it->second->isSensitive())
            sensParamSecs.push_back(it->second->getSec());
    }

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    return sensParamSecs;
}


template <typename ParamType>
int ParameterContainer<ParamType>::getNumSens() const
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    int sensCount = 0;
    for (ConstMapIt it = _params.begin(); it != _params.end(); ++it)
    {
        if (it->second->isSensitive())
            sensCount++;
    }

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    return sensCount;
}


template <typename ParamType>
bool ParameterContainer<ParamType>::checkAll() const
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    bool ok = true;

    for (ConstMapIt it = _params.begin(); it != _params.end(); ++it)
        ok = ok && it->second->check();

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    return ok;
}


template <typename ParamType>
const std::string ParameterContainer<ParamType>::info() const
{
    std::ostringstream oss;
    for (ConstMapIt it = _params.begin(); it != _params.end(); ++it)
        oss << it->second->info() << std::endl;
    return oss.str();
}


template <typename ParamType>
bool ParameterContainer<ParamType>::contains(const ParameterName id, int comp, int sec)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    ConstMapIt it = _params.find(ParamID(id, comp, sec));
    if (it != _params.end())
    {
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
        return true;
    }
    else
    {
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
        return false;
    }
}

template <typename ParamType>
bool ParameterContainer<ParamType>::contains(const ParamID& param)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    ConstMapIt it = _params.find(param);
    if (it != _params.end())
    {
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
        return true;
    }
    else
    {
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
        return false;
    }
}


} // namespace cadet


#endif // PARAMETERCONTAINER_HPP_
