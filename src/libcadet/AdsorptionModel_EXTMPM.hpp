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

#ifndef ADSORPTIONMODEL_EXTMPM_HPP_
#define ADSORPTIONMODEL_EXTMPM_HPP_

#include "AdsorptionModel.hpp"

namespace cadet
{

/// \brief Implementation of the Mobile Phyase Modulators adsorption model with externally dependent parameters
/// All parameter-related functions are inherited from the ParameterContainer class
class AdsorptionModel_EXTMPM : public AdsorptionModel
{
public:

    // Constructor
    AdsorptionModel_EXTMPM(const SimulatorPImpl& sim) :
        AdsorptionModel(sim, EXTERNAL_MOBILE_PHASE_MODULATORS)
    {
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        const double inf = std::numeric_limits<double>::infinity();

        this->configure();
        log::emit<Debug1>() << CURRENT_FUNCTION << ": Configured" << log::endl;

        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            addParam(Parameter<active> (EXTMPM_KA,        e2s(EXTMPM_KA),        comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTMPM_KA_T,      e2s(EXTMPM_KA_T),      comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTMPM_KA_TT,     e2s(EXTMPM_KA_TT),     comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTMPM_KA_TTT,    e2s(EXTMPM_KA_TTT),    comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTMPM_KD,        e2s(EXTMPM_KD),        comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTMPM_KD_T,      e2s(EXTMPM_KD_T),      comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTMPM_KD_TT,     e2s(EXTMPM_KD_TT),     comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTMPM_KD_TTT,    e2s(EXTMPM_KD_TTT),    comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTMPM_QMAX,      e2s(EXTMPM_QMAX),      comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTMPM_QMAX_T,    e2s(EXTMPM_QMAX_T),    comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTMPM_QMAX_TT,   e2s(EXTMPM_QMAX_TT),   comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTMPM_QMAX_TTT,  e2s(EXTMPM_QMAX_TTT),  comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTMPM_BETA,      e2s(EXTMPM_BETA),      comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTMPM_BETA_T,    e2s(EXTMPM_BETA_T),    comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTMPM_BETA_TT,   e2s(EXTMPM_BETA_TT),   comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTMPM_BETA_TTT,  e2s(EXTMPM_BETA_TTT),  comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTMPM_GAMMA,     e2s(EXTMPM_GAMMA),     comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTMPM_GAMMA_T,   e2s(EXTMPM_GAMMA_T),   comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTMPM_GAMMA_TT,  e2s(EXTMPM_GAMMA_TT),  comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTMPM_GAMMA_TTT, e2s(EXTMPM_GAMMA_TTT), comp, -1, 0.0, 0.0, -inf, true, inf, true));
        }

        log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    // Destructor
    virtual ~AdsorptionModel_EXTMPM()
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
        { evaluateResidual<active, active, active>(t, z, comp, q, res); }
    virtual void evaluateResidual(const double t, const double z, const int comp, const active * q, active * res, const double * p) const
        { evaluateResidual<active, active, double>(t, z, comp, q, res); }
    virtual void evaluateResidual(const double t, const double z, const int comp, const double * q, active * res, const active * p) const
        { evaluateResidual<double, active, active>(t, z, comp, q, res); }
    virtual void evaluateResidual(const double t, const double z, const int comp, const double * q, double * res, const double * p) const
        { evaluateResidual<double, double, double>(t, z, comp, q, res); }


    virtual void setJacobian(const double t, const double z, const int comp, const double* q, double* jac) const throw (CadetException)
    {
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        const double              ka        = getValue<double>           (EXTMPM_KA,     comp);
        const double              ka_T      = getValue<double>           (EXTMPM_KA_T,   comp);
        const double              ka_TT     = getValue<double>           (EXTMPM_KA_TT,  comp);
        const double              ka_TTT    = getValue<double>           (EXTMPM_KA_TTT, comp);
        const double              kd        = getValue<double>           (EXTMPM_KD,     comp);
        const double              kd_T      = getValue<double>           (EXTMPM_KD_T,   comp);
        const double              kd_TT     = getValue<double>           (EXTMPM_KD_TT,  comp);
        const double              kd_TTT    = getValue<double>           (EXTMPM_KD_TTT, comp);
        const std::vector<double> qmax      = getValueForAllComp<double> (EXTMPM_QMAX);
        const std::vector<double> qmax_T    = getValueForAllComp<double> (EXTMPM_QMAX_T);
        const std::vector<double> qmax_TT   = getValueForAllComp<double> (EXTMPM_QMAX_TT);
        const std::vector<double> qmax_TTT  = getValueForAllComp<double> (EXTMPM_QMAX_TTT);
        const double              beta      = getValue<double>           (EXTMPM_BETA,      comp);
        const double              beta_T    = getValue<double>           (EXTMPM_BETA_T,    comp);
        const double              beta_TT   = getValue<double>           (EXTMPM_BETA_TT,   comp);
        const double              beta_TTT  = getValue<double>           (EXTMPM_BETA_TTT,  comp);
        const double              gamma     = getValue<double>           (EXTMPM_GAMMA,     comp);
        const double              gamma_T   = getValue<double>           (EXTMPM_GAMMA_T,   comp);
        const double              gamma_TT  = getValue<double>           (EXTMPM_GAMMA_TT,  comp);
        const double              gamma_TTT = getValue<double>           (EXTMPM_GAMMA_TTT, comp);

        // Temperature
        double temp;
        _externalBase->externalProfile(z, t, &temp);

        // Only protein components
        if (comp > 0)
        {
            // Liquid phase concentration
            const double* c = q - _cc.ncomp();

            // Liquid phase salt concentration
            const double c0 = c[-comp];

            const double finalKa = ka + temp * (ka_T + temp * (ka_TT + temp * ka_TTT));
            const double finalKd = kd + temp * (kd_T + temp * (kd_TT + temp * kd_TTT));
            const double finalGamma = gamma + temp * (gamma_T + temp * (gamma_TT + temp * gamma_TTT));
            const double finalBeta = beta + temp * (beta_T + temp * (beta_TT + temp * beta_TTT));

            const double ka_mpm = finalKa * exp(finalGamma * c0);
            const double kd_mpm = finalKd * pow(c0, finalBeta);

            double qsum = 1.0;
            for (int j = 1; j < _cc.ncomp(); ++j)
            {
                const double finalQmax = qmax.at(j) + temp * (qmax_T.at(j) + temp * (qmax_TT.at(j) + temp * qmax_TTT.at(j)));
                qsum -= q[-comp + j] / finalQmax;
            }

            // Jacobian
            const double finalQmax = qmax.at(comp) + temp * (qmax_T.at(comp) + temp * (qmax_TT.at(comp) + temp * qmax_TTT.at(comp)));
            for (int j = 1; j < _cc.ncomp(); ++j)
            {
                const double finalQmaxJ = qmax.at(j) + temp * (qmax_T.at(j) + temp * (qmax_TT.at(j) + temp * qmax_TTT.at(j)));
                jac[-comp + j] = ka_mpm * (*c) * finalQmax / finalQmaxJ;              // dresi/dqj
            }

            jac[0]            += kd_mpm;                                            // dresi/dqi
            jac[-_cc.ncomp()]  = -ka_mpm * finalQmax * qsum;                        // dresi/dci

            jac[-_cc.ncomp() - comp] = -ka_mpm * (*c) * finalQmax * qsum * finalGamma
                    + finalKd * finalBeta * (*q) * pow(c0, finalBeta - 1);            // dresi/dc0
        }
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

private:
    template<typename StateType, typename ResidType, typename ParamType>
    void evaluateResidual(const double t, const double z, const int comp, const StateType *q, ResidType *res) const
    {
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        ParamType              ka        = getValue<ParamType>           (EXTMPM_KA,     comp);
        ParamType              ka_T      = getValue<ParamType>           (EXTMPM_KA_T,   comp);
        ParamType              ka_TT     = getValue<ParamType>           (EXTMPM_KA_TT,  comp);
        ParamType              ka_TTT    = getValue<ParamType>           (EXTMPM_KA_TTT, comp);
        ParamType              kd        = getValue<ParamType>           (EXTMPM_KD,     comp);
        ParamType              kd_T      = getValue<ParamType>           (EXTMPM_KD_T,   comp);
        ParamType              kd_TT     = getValue<ParamType>           (EXTMPM_KD_TT,  comp);
        ParamType              kd_TTT    = getValue<ParamType>           (EXTMPM_KD_TTT, comp);
        std::vector<ParamType> qmax      = getValueForAllComp<ParamType> (EXTMPM_QMAX);
        std::vector<ParamType> qmax_T    = getValueForAllComp<ParamType> (EXTMPM_QMAX_T);
        std::vector<ParamType> qmax_TT   = getValueForAllComp<ParamType> (EXTMPM_QMAX_TT);
        std::vector<ParamType> qmax_TTT  = getValueForAllComp<ParamType> (EXTMPM_QMAX_TTT);
        ParamType              beta      = getValue<ParamType>           (EXTMPM_BETA,      comp);
        ParamType              beta_T    = getValue<ParamType>           (EXTMPM_BETA_T,    comp);
        ParamType              beta_TT   = getValue<ParamType>           (EXTMPM_BETA_TT,   comp);
        ParamType              beta_TTT  = getValue<ParamType>           (EXTMPM_BETA_TTT,  comp);
        ParamType              gamma     = getValue<ParamType>           (EXTMPM_GAMMA,     comp);
        ParamType              gamma_T   = getValue<ParamType>           (EXTMPM_GAMMA_T,   comp);
        ParamType              gamma_TT  = getValue<ParamType>           (EXTMPM_GAMMA_TT,  comp);
        ParamType              gamma_TTT = getValue<ParamType>           (EXTMPM_GAMMA_TTT, comp);

        // Liquid phase concentration
        const StateType* c = q - _cc.ncomp();

        // Temperature
        double temp;
        _externalBase->externalProfile(z, t, &temp);

        if (comp == 0)
        {
            // Salt component
            *res = (ResidType) 0.0;
        }
        else
        {
            // Protein components
            
            // Liquid phase salt concentration
            const StateType c0 = c[-comp];

            ParamType finalKa = ka + temp * (ka_T + temp * (ka_TT + temp * ka_TTT));
            ParamType finalKd = kd + temp * (kd_T + temp * (kd_TT + temp * kd_TTT));
            ParamType finalGamma = gamma + temp * (gamma_T + temp * (gamma_TT + temp * gamma_TTT));
            ParamType finalBeta = beta + temp * (beta_T + temp * (beta_TT + temp * beta_TTT));
            ParamType finalQmax = qmax.at(comp) + temp * (qmax_T.at(comp) + temp * (qmax_TT.at(comp) + temp * qmax_TTT.at(comp)));
            
            ResidType ka_mpm = finalKa * exp(finalGamma * c0);
            ResidType kd_mpm = finalKd * pow(c0, finalBeta);

            ResidType qsum = 1.0;
            for (int j = 1; j < _cc.ncomp(); ++j)
            {
                ParamType finalQmaxJ = qmax.at(j) + temp * (qmax_T.at(j) + temp * (qmax_TT.at(j) + temp * qmax_TTT.at(j)));
                qsum -= q[-comp + j] / finalQmaxJ;
            }

            // Residual
            *res = - (ka_mpm * (*c) * finalQmax * qsum - kd_mpm * (*q));
        }

        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }
};

} // namespace cadet

#endif // ADSORPTIONMODEL_EXTMPM_HPP_
