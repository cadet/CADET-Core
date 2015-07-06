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

#ifndef ADSORPTIONMODEL_MCL_HPP_
#define ADSORPTIONMODEL_MCL_HPP_

#include "AdsorptionModel.hpp"

namespace cadet
{

/// \brief Implementation of the Multi Component Langmuir adsorption model
/// All parameter-related functions are inherited from the ParameterContainer class
class AdsorptionModel_MCL : public AdsorptionModel
{
public:

    // Constructor
    AdsorptionModel_MCL(const SimulatorPImpl& sim) :
        AdsorptionModel(sim, MULTI_COMPONENT_LANGMUIR)
    {
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        const double inf = std::numeric_limits<double>::infinity();

        this->configure();
        log::emit<Debug1>() << CURRENT_FUNCTION << ": Configured" << log::endl;

        _kA.reserve(_cc.ncomp());
        _kD.reserve(_cc.ncomp());
        _qMax.reserve(_cc.ncomp());

        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            _kA.push_back(Parameter<active> (MCL_KA,   e2s(MCL_KA),   comp, -1, 0.0, 0.0, 0.0, false, inf, true));
            addParam(_kA[comp]);
            _kD.push_back(Parameter<active> (MCL_KD,   e2s(MCL_KD),   comp, -1, 0.0, 0.0, 0.0, false, inf, true));
            addParam(_kD[comp]);
            _qMax.push_back(Parameter<active> (MCL_QMAX, e2s(MCL_QMAX), comp, -1, 0.0, 0.0, 0.0, true,  inf, true));
            addParam(_qMax[comp]);
        }

        log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    // Destructor
    virtual ~AdsorptionModel_MCL()
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
        const double              ka   = _kA[comp].getValue<double>();
        const double              kd   = _kD[comp].getValue<double>();

        // Liquid phase concentration
        const double* c = q -_cc.ncomp();

        double qsum = 1.0;
        for (int j = 0; j < _cc.ncomp(); ++j)
            qsum -= q[-comp + j] / _qMax[j].getValue<double>();

        // Jacobian
        for (int j = 0; j < _cc.ncomp(); ++j)
            jac[-comp + j] = ka * (*c) * _qMax[comp].getValue<double>() / _qMax[j].getValue<double>();  // dres/dq_j

        jac[0]            += kd;                                    // dres/dq_i
        jac[-_cc.ncomp()]  = -ka * _qMax[comp].getValue<double>() * qsum;            // dres/dc
    }

private:

    template<typename StateType, typename ResidType, typename ParamType>
    void evaluateResidual(const int comp, const StateType *q, ResidType *res) const
    {
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        const ParamType              ka   = _kA[comp].getValue<ParamType>();
        const ParamType              kd   = _kD[comp].getValue<ParamType>();

        // Liquid phase concentration
        const StateType* c = q -_cc.ncomp();

        ResidType qsum = 1.0;
        for (int j = 0; j < _cc.ncomp(); ++j)
            qsum -= q[-comp + j] / _qMax[j].getValue<ParamType>();

        // Residual
        *res = - (ka * (*c) * _qMax[comp].getValue<ParamType>() * qsum - kd * *q);

        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    std::vector<Parameter<active>>  _kA;
    std::vector<Parameter<active>>  _kD;
    std::vector<Parameter<active>>  _qMax;
};

} // namespace cadet

#endif // ADSORPTIONMODEL_MCL_HPP_
