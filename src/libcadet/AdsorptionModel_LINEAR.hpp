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

#ifndef ADSORPTIONMODEL_LINEAR_HPP_
#define ADSORPTIONMODEL_LINEAR_HPP_

#include "AdsorptionModel.hpp"

namespace cadet
{

/// \brief Implementation of the linear adsorption model
/// All parameter-related functions are inherited from the ParameterContainer class
class AdsorptionModel_LINEAR : public AdsorptionModel
{
public:

    // Constructor
    AdsorptionModel_LINEAR(const SimulatorPImpl& sim) :
        AdsorptionModel(sim, LINEAR)
    {
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        double inf = std::numeric_limits<double>::infinity();

        this->configure();
        log::emit<Debug1>() << CURRENT_FUNCTION << ": Configured" << log::endl;

        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            addParam(Parameter<active> (LIN_KA,   e2s(LIN_KA),   comp, -1, 0.0, 0.0, 0.0, false, inf, true));
            addParam(Parameter<active> (LIN_KD,   e2s(LIN_KD),   comp, -1, 0.0, 0.0, 0.0, false, inf, true));
        }

        log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    // Destructor
    virtual ~AdsorptionModel_LINEAR()
    {
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    // Public members

    virtual void setIsKinetic(bool isKinetic)
    {
        _isKinetic = isKinetic;
        for (int comp = 0; comp < _cc.ncomp(); ++comp)
            _isDifferential.at(comp) = isKinetic;
    }


    virtual void evaluateResidual(const double t, const double z, const int comp, const active * q, active * res, const active * p) const
        { evaluateResidual<active, active, active>(comp, q, res); }
    virtual void evaluateResidual(const double t, const double z, const int comp, const active * q, active * res, const double * p) const
        { evaluateResidual<active, active, double>(comp, q, res); }
    virtual void evaluateResidual(const double t, const double z, const int comp, const double * q, active * res, const active * p) const
        { evaluateResidual<double, active, active>(comp, q, res); }
    virtual void evaluateResidual(const double t, const double z, const int comp, const double * q, double * res, const double * p) const
        { evaluateResidual<double, double, double>(comp, q, res); }


    virtual void setJacobian(const double t, const double z, const int comp, const double* q, double* jac) const throw (CadetException)
    {
        const double ka = getValue<double>(LIN_KA, comp);
        const double kd = getValue<double>(LIN_KD, comp);

        jac[0]            = kd;   // dres/dq_i
        jac[-_cc.ncomp()] = -ka;  // dres/dc_i
    }

private:

    template<typename StateType, typename ResidType, typename ParamType>
    void evaluateResidual(const int comp, const StateType *q, ResidType *res) const
    {
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        ParamType ka = getValue<ParamType>(LIN_KA, comp);
        ParamType kd = getValue<ParamType>(LIN_KD, comp);

        // Liquid phase concentration
        const StateType* c = q -_cc.ncomp();

        // Residual
        *res = - (ka * (*c) - kd * (*q));

        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }
};

} // namespace cadet

#endif // ADSORPTIONMODEL_LINEAR_HPP_
