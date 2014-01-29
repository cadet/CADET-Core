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

#ifndef PARAMETER_HPP_
#define PARAMETER_HPP_

#include <string>
#include <ostream>
#include <sstream>
#include <iomanip>

#include "CadetLogger.hpp"
#include "active.hpp"

namespace cadet
{

// A Parameter encapsulates a passive or an active data type,
// which can be used for algorithmic differentiation.
// Therefore, this class provides getADValue() and setADValue().
//
// Parameter objects are *NOT* intended to
// be used in calculations. Hence overloaded operators
// are not provided.
template <typename ParamType>
class Parameter
{
public:

    // Constructor
    Parameter(ParameterName id, const std::string& name, const int comp, const double value, const double absTolS,
              const double lowerbound, const bool lowerbndstrict,
              const double upperbound, const bool upperbndstrict);

    // Copy constructor
    Parameter(const Parameter<ParamType>& p);

    // Destructor
    ~Parameter();

    // Copy assignment operator
    const Parameter<ParamType>& operator=(const Parameter<ParamType>& p);


    // Inline member functions
    inline ParameterName getId() const { return _id; }
    inline const std::string& getName() const { return _name; }
    inline int getComp() const { return _comp; }
    inline double getAbsTolS() const { return _absTolS; }
    inline void setAbsTolS(double absTolS) { _absTolS = absTolS; }
    inline void setSensitive() { _issensitive = true; }
    inline void resetSensitive() { _issensitive = false; }
    inline bool isSensitive() const { return _issensitive; }

    // Check if parameter has a meaningful value
    bool check() const;

    // GetValue() and setValue() need to be implemented for any ParamType
    template <typename ReturnType>
    inline const ReturnType& getValue() const;

    inline void setValue(const double& value);

    // The AD functions should only be implemented for active types
    // (they are meaningless for passives!)
    inline double getADValue(int dir) const;
    inline void setADValue(int dir, double adval);

    // A string containing info on this parameter
    const std::string info() const;

private:

    ParameterName   _id;
    ParamType       _value;
    int             _comp;            //!< chemical component index. -1 if global
    double          _absTolS;         //!< absolute tolerance for sensitivity computation
    double          _lowerbound;
    double          _upperbound;
    bool            _lowerbndstrict;  //!< decides whether the lower bound is strict ( < ) or loose ( <= )
    bool            _upperbndstrict;  //!< decides whether the upper bound is strict ( > ) or loose ( >= )
    bool            _issensitive;
    std::string     _name;

    double*         _returnval;         //!< used only for returning references on double values
                                        // really a dirty hack ... should be replaced by some nicer construct

    inline std::string _type() const;

};



// Overloaded operator<< for printing parameter info to streams
template <typename ParamType>
std::ostream& operator<<(std::ostream& out, const Parameter<ParamType>& p)
{
    out << p.info();
    return out;
}



// ====================================================================================================================
//    IMPLEMENTATION PART
// ====================================================================================================================

// Constructor
template <typename ParamType>
Parameter<ParamType>::Parameter(ParameterName id, const std::string& name, const int comp, const double value, const double absTolS,
        const double lowerbound, const bool lowerbndstrict,
        const double upperbound, const bool upperbndstrict) :
    _id             (id),
    _value          (value),
    _comp           (comp),
    _absTolS        (absTolS),
    _lowerbound     (lowerbound),
    _upperbound     (upperbound),
    _lowerbndstrict (lowerbndstrict),
    _upperbndstrict (upperbndstrict),
    _issensitive    (false),
    _name           (name)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    // allocate memory for returning on the heap
    _returnval  = new double;
    *_returnval = 0.0;

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


// copy constructor
template <typename ParamType>
Parameter<ParamType>::Parameter(const Parameter<ParamType>& p) :
    _id             (p._id),
    _value          (p._value),
    _comp           (p._comp),
    _absTolS        (p._absTolS),
    _lowerbound     (p._lowerbound),
    _upperbound     (p._upperbound),
    _lowerbndstrict (p._lowerbndstrict),
    _upperbndstrict (p._upperbndstrict),
    _issensitive    (p._issensitive),
    _name           (p._name)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    _returnval  = new double;
    *_returnval = *(p._returnval);

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


// destructor
template <typename ParamType>
Parameter<ParamType>::~Parameter()
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    delete _returnval;

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


// copy assignment operator
template <typename ParamType>
const Parameter<ParamType>& Parameter<ParamType>::operator=(const Parameter<ParamType>& p)
{
    _id             = p._id;
    _value          = p._value;
    _comp           = p._comp;
    _absTolS        = p._absTolS;
    _lowerbound     = p._lowerbound;
    _upperbound     = p._upperbound;
    _lowerbndstrict = p._lowerbndstrict;
    _upperbndstrict = p._upperbndstrict;
    _issensitive    = p._issensitive;
    _name           = p._name;
    *_returnval     = *(p._returnval);

    return *this;
}


// check if parameter has a meaningful value
template <typename ParamType>
bool Parameter<ParamType>::check() const
{
    bool ok = true;

    ok = ok && (_lowerbndstrict ? _value > _lowerbound : _value >= _lowerbound);
    ok = ok && (_upperbndstrict ? _value < _upperbound : _value <= _upperbound);

    return ok;
}


// a string containing info on this parameter
template <typename ParamType>
const std::string Parameter<ParamType>::info() const
{
    bool ok = check();

    std::ostringstream info;
    info                 << std::setw(17) << std::left  << _name.substr(0,17);
    info << " | Id: "    << std::setw(3)  << std::right << _id;
    info << " | Comp: "  << std::setw(3)  << _comp;
    info << " | Value: " << std::setw(9)  << getValue<double> ();
    info << " | Type: "  << std::setw(8)  << _type();
    info << " | LB: "    << std::setw(3)  << _lowerbound;
    info << " | UB: "    << std::setw(3)  << _upperbound;
    info << " | Sens: "  << std::setw(1)  << _issensitive;
    info << " | ATolS: " << std::setw(9)  << _absTolS;
    if (ok)
        info << " | in bnds";
    else
        info << " | OUT OF BOUNDS!";
    return info.str();
}


//==============================================================
// specializations for doubles
//==============================================================
template <> template <>
inline const double & Parameter<double>::getValue<double>() const
{
    return _value;
}

template <>
inline void Parameter<double>::setValue(const double& value)
{
    _value = value;
}

template <>
inline std::string Parameter<double>::_type() const
{
    return "double";
}
//==============================================================



//==============================================================
// specializations for actives
//==============================================================

// Constructor specialization
template <>
inline Parameter<active>::Parameter(ParameterName id, const std::string& name, const int comp, const double value, const double absTolS,
        const double lowerbound, const bool lowerbndstrict,
        const double upperbound, const bool upperbndstrict) :
    _id             (id),
    _value          (value),
    _comp           (comp),
    _absTolS        (absTolS),
    _lowerbound     (lowerbound),
    _upperbound     (upperbound),
    _lowerbndstrict (lowerbndstrict),
    _upperbndstrict (upperbndstrict),
    _issensitive    (false),
    _name           (name)
{
    // Initialize all directional derivatives with 0.0
    for (int dir = 0; dir < active::getMaxDirections(); ++dir)
        _value.setADValue(dir, 0.0);

    // allocate memory for returning on the heap
    _returnval  = new double;
    *_returnval = 0.0;
}

template <> template <>
inline const double& Parameter<active>::getValue<double>() const
{
    *_returnval = _value.getValue();
    return *_returnval;
}

template <> template <>
inline const active& Parameter<active>::getValue<active>() const
{
    return _value;
}

template <>
inline void Parameter<active>::setValue(const double& value)
{
    _value.setValue(value);
}

template <>
inline double Parameter<active>::getADValue(int dir) const
{
    return _value.getADValue(dir);
}

template <>
inline void Parameter<active>::setADValue(int dir, double adval)
{
    _value.setADValue(dir, adval);
}

template <>
inline std::string Parameter<active>::_type() const
{
    return "active";
}
//==============================================================


} // namespace cadet

#endif // PARAMETER_HPP_
