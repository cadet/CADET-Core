// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright (c) 2008-2012: Eric von Lieres¹, Joel Andersson,
//                           Andreas Puettmann¹, Sebastian Schnittert¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0
//  which accompanies this distribution, and is available at
//  http://www.gnu.org/licenses/gpl.html
//  ---------------------------------------------------------------------------
//  Author    Sebastian Schnittert <schnittert@gmail.com>
//  Version   $Id: AdsorptionModel_THM.hpp 519 2012-10-08 12:55:41Z schnittert $
// =============================================================================

#ifndef ADSORPTIONMODEL_EXT_SMA_HPP_
#define ADSORPTIONMODEL_EXT_SMA_HPP_

#include "AdsorptionModel.hpp"

namespace cadet
{

/// \brief Implementation of the Steric Mass Action adsorption model with externally dependent parameters
/// All parameter-related functions are inherited from the ParameterContainer class
class AdsorptionModel_EXTSMA : public AdsorptionModel
{
public:

    // Constructor
    AdsorptionModel_EXTSMA(const SimulatorPImpl& sim) :
        AdsorptionModel(sim, EXTERNAL_STERIC_MASS_ACTION)
    {
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        double inf = std::numeric_limits<double>::infinity();

        this->configure();
        log::emit<Debug1>() << CURRENT_FUNCTION << ": Configured" << log::endl;

        addParam(Parameter<active> (EXTSMA_LAMBDA,     e2s(EXTSMA_LAMBDA),     -1, -1, 0.0, 0.0, -inf, true, inf, true));
        addParam(Parameter<active> (EXTSMA_LAMBDA_T,   e2s(EXTSMA_LAMBDA_T),   -1, -1, 0.0, 0.0, -inf, true, inf, true));
        addParam(Parameter<active> (EXTSMA_LAMBDA_TT,  e2s(EXTSMA_LAMBDA_TT),  -1, -1, 0.0, 0.0, -inf, true, inf, true));
        addParam(Parameter<active> (EXTSMA_LAMBDA_TTT, e2s(EXTSMA_LAMBDA_TTT), -1, -1, 0.0, 0.0, -inf, true, inf, true));

        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            addParam(Parameter<active> (EXTSMA_KA,         e2s(EXTSMA_KA),         comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTSMA_KA_T,       e2s(EXTSMA_KA_T),       comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTSMA_KA_TT,      e2s(EXTSMA_KA_TT),      comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTSMA_KA_TTT,     e2s(EXTSMA_KA_TTT),     comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTSMA_KD,         e2s(EXTSMA_KD),         comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTSMA_KD_T,       e2s(EXTSMA_KD_T),       comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTSMA_KD_TT,      e2s(EXTSMA_KD_TT),      comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTSMA_KD_TTT,     e2s(EXTSMA_KD_TTT),     comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTSMA_NU,         e2s(EXTSMA_NU),         comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTSMA_NU_T,       e2s(EXTSMA_NU_T),       comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTSMA_NU_TT,      e2s(EXTSMA_NU_TT),      comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTSMA_NU_TTT,     e2s(EXTSMA_NU_TTT),     comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTSMA_SIGMA,      e2s(EXTSMA_SIGMA),      comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTSMA_SIGMA_T,    e2s(EXTSMA_SIGMA_T),    comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTSMA_SIGMA_TT,   e2s(EXTSMA_SIGMA_TT),   comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTSMA_SIGMA_TTT,  e2s(EXTSMA_SIGMA_TTT),  comp, -1, 0.0, 0.0, -inf, true, inf, true));
        }

        log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    // Destructor
    virtual ~AdsorptionModel_EXTSMA()
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

        const double              ka        = getValue<double>           (EXTSMA_KA,     comp);
        const double              ka_T      = getValue<double>           (EXTSMA_KA_T,   comp);
        const double              ka_TT     = getValue<double>           (EXTSMA_KA_TT,  comp);
        const double              ka_TTT    = getValue<double>           (EXTSMA_KA_TTT, comp);
        const double              kd        = getValue<double>           (EXTSMA_KD,     comp);
        const double              kd_T      = getValue<double>           (EXTSMA_KD_T,   comp);
        const double              kd_TT     = getValue<double>           (EXTSMA_KD_TT,  comp);
        const double              kd_TTT    = getValue<double>           (EXTSMA_KD_TTT, comp);
        const std::vector<double> nu        = getValueForAllComp<double> (EXTSMA_NU);
        const std::vector<double> nu_T      = getValueForAllComp<double> (EXTSMA_NU_T);
        const std::vector<double> nu_TT     = getValueForAllComp<double> (EXTSMA_NU_TT);
        const std::vector<double> nu_TTT    = getValueForAllComp<double> (EXTSMA_NU_TTT);
        const std::vector<double> sigma     = getValueForAllComp<double> (EXTSMA_SIGMA);
        const std::vector<double> sigma_T   = getValueForAllComp<double> (EXTSMA_SIGMA_T);
        const std::vector<double> sigma_TT  = getValueForAllComp<double> (EXTSMA_SIGMA_TT);
        const std::vector<double> sigma_TTT = getValueForAllComp<double> (EXTSMA_SIGMA_TTT);

        // Temperature
        double temp;
        _externalBase->externalProfile(z, t, &temp);

        // Liquid phase concentration
        const double* c = q -_cc.ncomp();

        if (comp == 0)  // Salt component
        {
            jac[0] = 1.0;
            for (int j = 1; j < _cc.ncomp(); ++j)
            {
                jac[j] = nu.at(j) + temp * (nu_T.at(j) + temp * (nu_TT.at(j) + temp * nu_TTT.at(j)));
            }
        }
        else  // Protein component
        {
            // Salt concentrations in liquid and solid phase
            const double c0 = c[-comp];
            const double q0 = q[-comp];

            double q0_bar = q0;

            for (int j = 1; j < _cc.ncomp(); ++j)
            {
                const double finalSigma = sigma.at(j) + temp * (sigma_T.at(j) + temp * (sigma_TT.at(j) + temp * sigma_TTT.at(j)));
                q0_bar -= finalSigma * q[-comp + j];
            }

            const double finalNu = nu.at(comp) + temp * (nu_T.at(comp) + temp * (nu_TT.at(comp) + temp * nu_TTT.at(comp)));
            const double c0_pow_nu     = pow(c0, finalNu);
            const double q0_bar_pow_nu = pow(q0_bar, finalNu);

            const double finalKa = ka + temp * (ka_T + temp * (ka_TT + temp * ka_TTT));
            const double finalKd = kd + temp * (kd_T + temp * (kd_TT + temp * kd_TTT));

            // Jacobian
            jac[-_cc.ncomp() - comp] = finalKd * (*q) * finalNu * c0_pow_nu / c0;                      // dres_i / dc0
            jac[-_cc.ncomp()] = -finalKa * q0_bar_pow_nu;                                              // dres_i / dci
            jac[-comp] = -finalKa * (*c) * finalNu * q0_bar_pow_nu / q0_bar;                           // dres_i / dq0

            for (int j = 1; j < _cc.ncomp(); ++j)
            {
                const double finalSigma = sigma.at(j) + temp * (sigma_T.at(j) + temp * (sigma_TT.at(j) + temp * sigma_TTT.at(j)));
                jac[-comp + j] = -finalKa * (*c) * finalNu * q0_bar_pow_nu / q0_bar * (-finalSigma);  // dres_i / dqj
            }

            jac[0] += finalKd * c0_pow_nu;                                                               // dres_i / dqi
        }

        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

private:

    template<typename StateType, typename ResidType, typename ParamType>
    void evaluateResidual(const double t, const double z, const int comp, const StateType *q, ResidType *res) const
    {
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        ParamType              ka         = getValue<ParamType>           (EXTSMA_KA,     comp);
        ParamType              ka_T       = getValue<ParamType>           (EXTSMA_KA_T,   comp);
        ParamType              ka_TT      = getValue<ParamType>           (EXTSMA_KA_TT,  comp);
        ParamType              ka_TTT     = getValue<ParamType>           (EXTSMA_KA_TTT, comp);
        ParamType              kd         = getValue<ParamType>           (EXTSMA_KD,     comp);
        ParamType              kd_T       = getValue<ParamType>           (EXTSMA_KD_T,   comp);
        ParamType              kd_TT      = getValue<ParamType>           (EXTSMA_KD_TT,  comp);
        ParamType              kd_TTT     = getValue<ParamType>           (EXTSMA_KD_TTT, comp);
        ParamType              lambda     = getValue<ParamType>           (EXTSMA_LAMBDA);
        ParamType              lambda_T   = getValue<ParamType>           (EXTSMA_LAMBDA_T);
        ParamType              lambda_TT  = getValue<ParamType>           (EXTSMA_LAMBDA_TT);
        ParamType              lambda_TTT = getValue<ParamType>           (EXTSMA_LAMBDA_TTT);
        std::vector<ParamType> nu         = getValueForAllComp<ParamType> (EXTSMA_NU);
        std::vector<ParamType> nu_T       = getValueForAllComp<ParamType> (EXTSMA_NU_T);
        std::vector<ParamType> nu_TT      = getValueForAllComp<ParamType> (EXTSMA_NU_TT);
        std::vector<ParamType> nu_TTT     = getValueForAllComp<ParamType> (EXTSMA_NU_TTT);
        std::vector<ParamType> sigma      = getValueForAllComp<ParamType> (EXTSMA_SIGMA);
        std::vector<ParamType> sigma_T    = getValueForAllComp<ParamType> (EXTSMA_SIGMA_T);
        std::vector<ParamType> sigma_TT   = getValueForAllComp<ParamType> (EXTSMA_SIGMA_TT);
        std::vector<ParamType> sigma_TTT  = getValueForAllComp<ParamType> (EXTSMA_SIGMA_TTT);

        // Temperature
        double temp;
        _externalBase->externalProfile(z, t, &temp);

        // Liquid phase concentration
        const StateType* c = q -_cc.ncomp();

        if (comp == 0)
        {
            // Salt component
            *res = *q - ( lambda + temp * (lambda_T + temp * (lambda_TT + temp * lambda_TTT)) );

            for (int j = 1; j < _cc.ncomp(); ++j)
            {
                ResidType finalNu = nu.at(j) + temp * (nu_T.at(j) + temp * (nu_TT.at(j) + temp * nu_TTT.at(j)));
                *res += finalNu * q[j];
            }
        }
        else
        {
            // Protein 
            // Salt concentrations in liquid and solid phase
            const StateType c0 = c[-comp];
            const StateType q0 = q[-comp];

            ResidType q0_bar = q0;

            for (int j = 1; j < _cc.ncomp(); ++j)
            {
                ResidType finalSigma = sigma.at(j) + temp * (sigma_T.at(j) + temp * (sigma_TT.at(j) + temp * sigma_TTT.at(j)));
                q0_bar -= finalSigma * q[-comp + j];
            }

            ResidType finalNu = nu.at(comp) + temp * (nu_T.at(comp) + temp * (nu_TT.at(comp) + temp * nu_TTT.at(comp)));
            ResidType c0_pow_nu = pow(c0, finalNu);
            ResidType q0_bar_pow_nu = pow(q0_bar, finalNu);

            ResidType finalKa = ka + temp * (ka_T + temp * (ka_TT + temp * ka_TTT));
            ResidType finalKd = kd + temp * (kd_T + temp * (kd_TT + temp * kd_TTT));

            // Residual
            *res = finalKd * (*q) * c0_pow_nu - finalKa * (*c) * q0_bar_pow_nu;
        }

        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }
};

} // namespace cadet

#endif // ADSORPTIONMODEL_THM_SMA_HPP_
