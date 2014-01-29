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

#ifndef ADSORPTIONMODEL_SAI_HPP_
#define ADSORPTIONMODEL_SAI_HPP_

#include "AdsorptionModel.hpp"

namespace cadet
{

/// \brief Implementation of the Self Association adsorption model
/// All parameter-related functions are inherited from the ParameterContainer class
class AdsorptionModel_SAI : public AdsorptionModel
{
public:

    // Constructor
    AdsorptionModel_SAI(const SimulatorPImpl & sim) :
        AdsorptionModel(sim, SELF_ASSOCIATION)
    {
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        double inf = std::numeric_limits<double>::infinity();

        addParam(Parameter<active> (SAI_LAMBDA, e2s(SAI_LAMBDA), -1, 0.0, 0.0, 0.0, false, inf, true));

        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            addParam(Parameter<active> (SAI_KA1,   e2s(SAI_KA1),   comp, 0.0, 0.0, 0.0, false, inf, true));
            addParam(Parameter<active> (SAI_KA2,   e2s(SAI_KA2),   comp, 0.0, 0.0, 0.0, false, inf, true));
            addParam(Parameter<active> (SAI_KD,    e2s(SAI_KD),    comp, 0.0, 0.0, 0.0, false, inf, true));
            addParam(Parameter<active> (SAI_NU,    e2s(SAI_NU),    comp, 0.0, 0.0, 0.0, false, inf, true));
            addParam(Parameter<active> (SAI_SIGMA, e2s(SAI_SIGMA), comp, 0.0, 0.0, 0.0, false, inf, true));
        }

        log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    // Destructor
    virtual ~AdsorptionModel_SAI()
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

        double              ka1    = getValue<double>           (SAI_KA1, comp);
        double              ka2    = getValue<double>           (SAI_KA2, comp);
        double              kd     = getValue<double>           (SAI_KD, comp);
        std::vector<double> nu     = getValueForAllComp<double> (SAI_NU);
        std::vector<double> sigma  = getValueForAllComp<double> (SAI_SIGMA);

        // Liquid phase concentration
        const double* c = q -_cc.ncomp();

        if (comp == 0) // Salt component
        {
            jac[0] = 1.0;
            for (int j = 1; j < _cc.ncomp(); ++j)
                jac[j] = nu.at(j);
        }
        else // Protein component
        {
            // Salt concentrations in liquid and solid phase
            double c0 = c[-comp];
            double q0 = q[-comp];

            double q0_bar = q0;
            
            for (int j = 1; j < _cc.ncomp(); ++j)
                q0_bar -= sigma.at(j) * q[-comp + j];

            double c0_pow_nu     = pow(c0, nu.at(comp));
            double q0_bar_pow_nu = pow(q0_bar, nu.at(comp));

            double c0_pow_nuMINUS1     = pow(c0, nu.at(comp) - 1);
            double q0_bar_pow_nuMINUS1 = pow(q0_bar, nu.at(comp) - 1);

            // Jacobian entrys
            jac[-_cc.ncomp() - comp] = kd * (*q) * (nu.at(comp) * c0_pow_nuMINUS1);         // dres / dc0

            jac[-comp] = -ka1 * (*c) * nu.at(comp) * q0_bar_pow_nuMINUS1
                    - ka2 * (*c * *c) * nu.at(comp) * q0_bar_pow_nuMINUS1;                  // dres / dq0_bar

            jac[-_cc.ncomp()] = -(ka1 * q0_bar_pow_nu) - (2 * ka2 * q0_bar_pow_nu * (*c));  // dres / dci

            for (int j = 1; j < _cc.ncomp(); ++j)
                jac[-comp + j] = -ka1 * (*c) * nu.at(comp)
                            * q0_bar_pow_nuMINUS1 * (-sigma.at(j))
                            - ka2 * (*c * *c) * nu.at(comp)
                            * q0_bar_pow_nuMINUS1 * (-sigma.at(j));                         // dres / dqi

            jac[0] += kd * c0_pow_nu;
        }
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }


private:
    // new residual formulation in k_a and k_d needs to be checked!!!
    template<typename StateType, typename ResidType, typename ParamType>
    void evaluateResidual(const int comp, const StateType *q, ResidType *res) const
    {
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        ParamType              lambda = getValue<ParamType>           (SAI_LAMBDA);
        ParamType              ka1    = getValue<ParamType>           (SAI_KA1, comp);
        ParamType              ka2    = getValue<ParamType>           (SAI_KA2, comp);
        ParamType              kd     = getValue<ParamType>           (SAI_KD, comp);
        std::vector<ParamType> nu     = getValueForAllComp<ParamType> (SAI_NU);
        std::vector<ParamType> sigma  = getValueForAllComp<ParamType> (SAI_SIGMA);

        // Liquid phase concentration
        const StateType* c = q -_cc.ncomp();

        if (comp == 0) // Salt component
        {
            *res = *q - lambda;

            for (int j = 1; j < _cc.ncomp(); ++j)
                *res += nu.at(j) * q[j];
        }
        else // Protein component
        {
            // Salt concentrations in liquid and solid phase
            StateType c0 = c[-comp];
            StateType q0 = q[-comp];

            ResidType q0_bar = q0;

            for (int j = 1; j < _cc.ncomp(); ++j)
                q0_bar -= sigma.at(j) * q[-comp + j];

            ResidType c0_pow_nu     = pow(c0,     nu.at(comp));
            ResidType q0_bar_pow_nu = pow(q0_bar, nu.at(comp));

            // Residual
            *res = - (ka1 * (*c) * q0_bar_pow_nu) - (ka2 * (*c * *c) * q0_bar_pow_nu) + (kd * c0_pow_nu * *q);
        }

        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }
};

} // namespace cadet

#endif // ADSORPTIONMODEL_SAI_HPP_
