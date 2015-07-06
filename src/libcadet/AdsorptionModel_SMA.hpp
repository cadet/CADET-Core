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

#ifndef ADSORPTIONMODEL_SMA_HPP_
#define ADSORPTIONMODEL_SMA_HPP_

#include "AdsorptionModel.hpp"

namespace cadet
{

/// \brief Implementation of the Steric Mass Action adsorption model
/// All parameter-related functions are inherited from the ParameterContainer class
class AdsorptionModel_SMA : public AdsorptionModel
{
public:

    // Constructor
    AdsorptionModel_SMA(const SimulatorPImpl & sim) :
        AdsorptionModel(sim, STERIC_MASS_ACTION),
        _lambda(SMA_LAMBDA, e2s(SMA_LAMBDA), -1, -1, 0.0, 0.0, 0.0, false, std::numeric_limits<double>::infinity(), true)
    {
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        const double inf = std::numeric_limits<double>::infinity();

        this->configure();
        log::emit<Debug1>() << CURRENT_FUNCTION << ": Configured" << log::endl;

        addParam(_lambda);

        _kA.reserve(_cc.ncomp());
        _kD.reserve(_cc.ncomp());
        _nu.reserve(_cc.ncomp());
        _sigma.reserve(_cc.ncomp());

        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            _kA.push_back(Parameter<active> (SMA_KA,    e2s(SMA_KA),    comp, -1, 0.0, 0.0, 0.0, false, inf, true));
            addParam(_kA[comp]);
            _kD.push_back(Parameter<active> (SMA_KD,    e2s(SMA_KD),    comp, -1, 0.0, 0.0, 0.0, false, inf, true));
            addParam(_kD[comp]);
            _nu.push_back(Parameter<active> (SMA_NU,    e2s(SMA_NU),    comp, -1, 0.0, 0.0, 0.0, false, inf, true));
            addParam(_nu[comp]);
            _sigma.push_back(Parameter<active> (SMA_SIGMA, e2s(SMA_SIGMA), comp, -1, 0.0, 0.0, 0.0, false, inf, true));
            addParam(_sigma[comp]);
        }

        log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    // Destructor
    virtual ~AdsorptionModel_SMA()
    {
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    // Public members

    virtual void setIsKinetic(bool isKinetic)
    {
        _isKinetic = isKinetic;
        for (int comp = 1; comp < _cc.ncomp(); ++comp)  // start only at comp 1, since salt-eq. is always non-differential
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
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        const double              ka     = _kA[comp].getValue<double>();
        const double              kd     = _kD[comp].getValue<double>();

        // Liquid phase concentration
        const double* c = q -_cc.ncomp();

        if (comp == 0)  // Salt component
        {
            jac[0] = 1.0;
            for (int j = 1; j < _cc.ncomp(); ++j)
                jac[j] = _nu[j].getValue<double>();
        }
        else  // Protein component
        {
            // Salt concentrations in liquid and solid phase
            const double c0 = c[-comp];
            const double q0 = q[-comp];

            double q0_bar = q0;

            for (int j = 1; j < _cc.ncomp(); ++j)
                q0_bar -= _sigma[j].getValue<double>() * q[-comp + j];

            const double c0_pow_nu     = pow(c0, _nu[comp].getValue<double>());
            const double q0_bar_pow_nu = pow(q0_bar, _nu[comp].getValue<double>());

            // Jacobian
            jac[-_cc.ncomp() - comp] = kd * (*q) * _nu[comp].getValue<double>() * c0_pow_nu / c0;                      // dres_i / dc0
            jac[-_cc.ncomp()] = -ka * q0_bar_pow_nu;                                                // dres_i / dci
            jac[-comp] = -ka * (*c) * _nu[comp].getValue<double>() * q0_bar_pow_nu / q0_bar;                              // dres_i / dq0

            for (int j = 1; j < _cc.ncomp(); ++j)
                jac[-comp + j] = -ka * (*c) * _nu[comp].getValue<double>() * q0_bar_pow_nu / q0_bar * (-_sigma[j].getValue<double>());  // dres_i / dqj

            jac[0] += kd * c0_pow_nu;                                                               // dres_i / dqi
        }

        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }


private:
    // new residual formulation in k_a and k_d needs to be checked!!!
    template<typename StateType, typename ResidType, typename ParamType>
    void evaluateResidual(const int comp, const StateType *q, ResidType *res) const
    {
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        const ParamType              ka     = _kA[comp].getValue<ParamType>();
        const ParamType              kd     = _kD[comp].getValue<ParamType>();

        // Liquid phase concentration
        const StateType* c = q -_cc.ncomp();

        if (comp == 0)
        {
            // Salt component
            *res = *q - _lambda.getValue<ParamType>();

            for (int j = 1; j < _cc.ncomp(); ++j)
                *res += _nu[j].getValue<ParamType>() * q[j];
        }
        else
        {
            // Protein 
            // Salt concentrations in liquid and solid phase
            const StateType c0 = c[-comp];
            const StateType q0 = q[-comp];

            ResidType q0_bar = q0;

            for (int j = 1; j < _cc.ncomp(); ++j)
                q0_bar -= _sigma[j].getValue<ParamType>() * q[-comp + j];

            ResidType c0_pow_nu = pow(c0, _nu[comp].getValue<ParamType>());
            ResidType q0_bar_pow_nu = pow(q0_bar, _nu[comp].getValue<ParamType>());

            // Residual
            *res = kd * (*q) * c0_pow_nu - ka * (*c) * q0_bar_pow_nu;
        }

        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    std::vector<Parameter<active>>  _kA;
    std::vector<Parameter<active>>  _kD;
    std::vector<Parameter<active>>  _nu;
    std::vector<Parameter<active>>  _sigma;

    Parameter<active> _lambda;    
};

} // namespace cadet

#endif // ADSORPTIONMODEL_SMA_HPP_
