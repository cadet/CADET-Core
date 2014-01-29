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
        AdsorptionModel(sim, STERIC_MASS_ACTION)
    {
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        double inf = std::numeric_limits<double>::infinity();

        addParam(Parameter<active> (SMA_LAMBDA, e2s(SMA_LAMBDA), -1, 0.0, 0.0, 0.0, false, inf, true));

        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            addParam(Parameter<active> (SMA_KA,    e2s(SMA_KA),    comp, 0.0, 0.0, 0.0, false, inf, true));
            addParam(Parameter<active> (SMA_KD,    e2s(SMA_KD),    comp, 0.0, 0.0, 0.0, false, inf, true));
            addParam(Parameter<active> (SMA_NU,    e2s(SMA_NU),    comp, 0.0, 0.0, 0.0, false, inf, true));
            addParam(Parameter<active> (SMA_SIGMA, e2s(SMA_SIGMA), comp, 0.0, 0.0, 0.0, false, inf, true));
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

        double              ka     = getValue<double>           (SMA_KA, comp);
        double              kd     = getValue<double>           (SMA_KD, comp);
        std::vector<double> nu     = getValueForAllComp<double> (SMA_NU);
        std::vector<double> sigma  = getValueForAllComp<double> (SMA_SIGMA);

        // Liquid phase concentration
        const double* c = q -_cc.ncomp();

        if (comp == 0)  // Salt component
        {
            jac[0] = 1.0;
            for (int j = 1; j < _cc.ncomp(); ++j)
                jac[j] = nu.at(j);
        }
        else  // Protein component
        {
            // Salt concentrations in liquid and solid phase
            double c0 = c[-comp];
            double q0 = q[-comp];

            double q0_bar = q0;

            for (int j = 1; j < _cc.ncomp(); ++j)
                q0_bar -= sigma.at(j) * q[-comp + j];

            double c0_pow_nu     = pow(c0, nu.at(comp));
            double q0_bar_pow_nu = pow(q0_bar, nu.at(comp));

            // Jacobian
            jac[-_cc.ncomp() - comp] = kd * *q * nu.at(comp) * c0_pow_nu / c0;                      // dres_i / dc0
            jac[-_cc.ncomp()] = -ka * q0_bar_pow_nu;                                                // dres_i / dci
            jac[-comp] = -ka * *c * nu.at(comp) * q0_bar_pow_nu / q0_bar;                              // dres_i / dq0

            for (int j = 1; j < _cc.ncomp(); ++j)
                jac[-comp + j] = -ka * *c * nu.at(comp) * q0_bar_pow_nu / q0_bar * (-sigma.at(j));  // dres_i / dqj

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

        ParamType              lambda = getValue<ParamType>           (SMA_LAMBDA);
        ParamType              ka     = getValue<ParamType>           (SMA_KA, comp);
        ParamType              kd     = getValue<ParamType>           (SMA_KD, comp);
        std::vector<ParamType> nu     = getValueForAllComp<ParamType> (SMA_NU);
        std::vector<ParamType> sigma  = getValueForAllComp<ParamType> (SMA_SIGMA);

        // Liquid phase concentration
        const StateType* c = q -_cc.ncomp();

        if (comp == 0)
        {
            // Salt component
            *res = *q - lambda;

            for (int j = 1; j < _cc.ncomp(); ++j)
                *res += nu.at(j) * q[j];
        }
        else
        {
            // Protein 
            // Salt concentrations in liquid and solid phase
            StateType c0 = c[-comp];
            StateType q0 = q[-comp];

            ResidType q0_bar = q0;

            for (int j = 1; j < _cc.ncomp(); ++j)
                q0_bar -= sigma.at(j) * q[-comp + j];

            ResidType c0_pow_nu = pow(c0, nu.at(comp));
            ResidType q0_bar_pow_nu = pow(q0_bar, nu.at(comp));

            // Residual
            *res = kd * *q * c0_pow_nu - ka * *c * q0_bar_pow_nu;
        }

        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }
};

} // namespace cadet

#endif // ADSORPTIONMODEL_SMA_HPP_
