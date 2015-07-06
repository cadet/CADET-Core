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

        _kA.reserve(_cc.ncomp());
        _kAT.reserve(_cc.ncomp());
        _kATT.reserve(_cc.ncomp());
        _kATTT.reserve(_cc.ncomp());

        _kD.reserve(_cc.ncomp());
        _kDT.reserve(_cc.ncomp());
        _kDTT.reserve(_cc.ncomp());
        _kDTTT.reserve(_cc.ncomp());

        _qMax.reserve(_cc.ncomp());
        _qMaxT.reserve(_cc.ncomp());
        _qMaxTT.reserve(_cc.ncomp());
        _qMaxTTT.reserve(_cc.ncomp());

        _beta.reserve(_cc.ncomp());
        _betaT.reserve(_cc.ncomp());
        _betaTT.reserve(_cc.ncomp());
        _betaTTT.reserve(_cc.ncomp());

        _gamma.reserve(_cc.ncomp());
        _gammaT.reserve(_cc.ncomp());
        _gammaTT.reserve(_cc.ncomp());
        _gammaTTT.reserve(_cc.ncomp());

        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            _kA.push_back(   Parameter<active> (EXTMPM_KA,        e2s(EXTMPM_KA),        comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kA[comp]);
            _kAT.push_back(  Parameter<active> (EXTMPM_KA_T,      e2s(EXTMPM_KA_T),      comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kAT[comp]);
            _kATT.push_back( Parameter<active> (EXTMPM_KA_TT,     e2s(EXTMPM_KA_TT),     comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kATT[comp]);
            _kATTT.push_back(Parameter<active> (EXTMPM_KA_TTT,    e2s(EXTMPM_KA_TTT),    comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kATTT[comp]);
            _kD.push_back(   Parameter<active> (EXTMPM_KD,        e2s(EXTMPM_KD),        comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kD[comp]);
            _kDT.push_back(  Parameter<active> (EXTMPM_KD_T,      e2s(EXTMPM_KD_T),      comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kDT[comp]);
            _kDTT.push_back( Parameter<active> (EXTMPM_KD_TT,     e2s(EXTMPM_KD_TT),     comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kDTT[comp]);
            _kDTTT.push_back(Parameter<active> (EXTMPM_KD_TTT,    e2s(EXTMPM_KD_TTT),    comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kDTTT[comp]);
            _qMax.push_back(   Parameter<active> (EXTMPM_QMAX,      e2s(EXTMPM_QMAX),      comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_qMax[comp]);
            _qMaxT.push_back(  Parameter<active> (EXTMPM_QMAX_T,    e2s(EXTMPM_QMAX_T),    comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_qMaxT[comp]);
            _qMaxTT.push_back( Parameter<active> (EXTMPM_QMAX_TT,   e2s(EXTMPM_QMAX_TT),   comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_qMaxTT[comp]);
            _qMaxTTT.push_back(Parameter<active> (EXTMPM_QMAX_TTT,  e2s(EXTMPM_QMAX_TTT),  comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_qMaxTTT[comp]);
            _beta.push_back(   Parameter<active> (EXTMPM_BETA,      e2s(EXTMPM_BETA),      comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_beta[comp]);
            _betaT.push_back(  Parameter<active> (EXTMPM_BETA_T,    e2s(EXTMPM_BETA_T),    comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_betaT[comp]);
            _betaTT.push_back( Parameter<active> (EXTMPM_BETA_TT,   e2s(EXTMPM_BETA_TT),   comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_betaTT[comp]);
            _betaTTT.push_back(Parameter<active> (EXTMPM_BETA_TTT,  e2s(EXTMPM_BETA_TTT),  comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_betaTTT[comp]);
            _gamma.push_back(   Parameter<active> (EXTMPM_GAMMA,     e2s(EXTMPM_GAMMA),     comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_gamma[comp]);
            _gammaT.push_back(  Parameter<active> (EXTMPM_GAMMA_T,   e2s(EXTMPM_GAMMA_T),   comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_gammaT[comp]);
            _gammaTT.push_back( Parameter<active> (EXTMPM_GAMMA_TT,  e2s(EXTMPM_GAMMA_TT),  comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_gammaTT[comp]);
            _gammaTTT.push_back(Parameter<active> (EXTMPM_GAMMA_TTT, e2s(EXTMPM_GAMMA_TTT), comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_gammaTTT[comp]);
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

        const double              ka        = _kA[comp].getValue<double>();
        const double              ka_T      = _kAT[comp].getValue<double>();
        const double              ka_TT     = _kATT[comp].getValue<double>();
        const double              ka_TTT    = _kATTT[comp].getValue<double>();
        const double              kd        = _kD[comp].getValue<double>();
        const double              kd_T      = _kDT[comp].getValue<double>();
        const double              kd_TT     = _kDTT[comp].getValue<double>();
        const double              kd_TTT    = _kDTTT[comp].getValue<double>();
        const double              beta      = _beta[comp].getValue<double>();
        const double              beta_T    = _betaT[comp].getValue<double>();
        const double              beta_TT   = _betaTT[comp].getValue<double>();
        const double              beta_TTT  = _betaTTT[comp].getValue<double>();
        const double              gamma     = _gamma[comp].getValue<double>();
        const double              gamma_T   = _gammaT[comp].getValue<double>();
        const double              gamma_TT  = _gammaTT[comp].getValue<double>();
        const double              gamma_TTT = _gammaTTT[comp].getValue<double>();

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
                const double finalQmax = _qMax[j].getValue<double>() + temp * (_qMaxT[j].getValue<double>() + temp * (_qMaxTT[j].getValue<double>() + temp * _qMaxTTT[j].getValue<double>()));
                qsum -= q[-comp + j] / finalQmax;
            }

            // Jacobian
            const double finalQmax = _qMax[comp].getValue<double>() + temp * (_qMaxT[comp].getValue<double>() + temp * (_qMaxTT[comp].getValue<double>() + temp * _qMaxTTT[comp].getValue<double>()));
            for (int j = 1; j < _cc.ncomp(); ++j)
            {
                const double finalQmaxJ = _qMax[j].getValue<double>() + temp * (_qMaxT[j].getValue<double>() + temp * (_qMaxTT[j].getValue<double>() + temp * _qMaxTTT[j].getValue<double>()));
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

        const ParamType ka        = _kA[comp].getValue<ParamType>();
        const ParamType ka_T      = _kAT[comp].getValue<ParamType>();
        const ParamType ka_TT     = _kATT[comp].getValue<ParamType>();
        const ParamType ka_TTT    = _kATTT[comp].getValue<ParamType>();
        const ParamType kd        = _kD[comp].getValue<ParamType>();
        const ParamType kd_T      = _kDT[comp].getValue<ParamType>();
        const ParamType kd_TT     = _kDTT[comp].getValue<ParamType>();
        const ParamType kd_TTT    = _kDTTT[comp].getValue<ParamType>();
        const ParamType beta      = _beta[comp].getValue<ParamType>();
        const ParamType beta_T    = _betaT[comp].getValue<ParamType>();
        const ParamType beta_TT   = _betaTT[comp].getValue<ParamType>();
        const ParamType beta_TTT  = _betaTTT[comp].getValue<ParamType>();
        const ParamType gamma     = _gamma[comp].getValue<ParamType>();
        const ParamType gamma_T   = _gammaT[comp].getValue<ParamType>();
        const ParamType gamma_TT  = _gammaTT[comp].getValue<ParamType>();
        const ParamType gamma_TTT = _gammaTTT[comp].getValue<ParamType>();

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

            const ParamType finalKa = ka + temp * (ka_T + temp * (ka_TT + temp * ka_TTT));
            const ParamType finalKd = kd + temp * (kd_T + temp * (kd_TT + temp * kd_TTT));
            const ParamType finalGamma = gamma + temp * (gamma_T + temp * (gamma_TT + temp * gamma_TTT));
            const ParamType finalBeta = beta + temp * (beta_T + temp * (beta_TT + temp * beta_TTT));
            const ParamType finalQmax = _qMax[comp].getValue<ParamType>() + temp * (_qMaxT[comp].getValue<ParamType>() + temp * (_qMaxTT[comp].getValue<ParamType>() + temp * _qMaxTTT[comp].getValue<ParamType>()));
            
            ResidType ka_mpm = finalKa * exp(finalGamma * c0);
            ResidType kd_mpm = finalKd * pow(c0, finalBeta);

            ResidType qsum = 1.0;
            for (int j = 1; j < _cc.ncomp(); ++j)
            {
                const ParamType finalQmaxJ = _qMax[j].getValue<ParamType>() + temp * (_qMaxT[j].getValue<ParamType>() + temp * (_qMaxTT[j].getValue<ParamType>() + temp * _qMaxTTT[j].getValue<ParamType>()));
                qsum -= q[-comp + j] / finalQmaxJ;
            }

            // Residual
            *res = - (ka_mpm * (*c) * finalQmax * qsum - kd_mpm * (*q));
        }

        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    std::vector<Parameter<active>>  _kA;
    std::vector<Parameter<active>>  _kAT;
    std::vector<Parameter<active>>  _kATT;
    std::vector<Parameter<active>>  _kATTT;

    std::vector<Parameter<active>>  _kD;
    std::vector<Parameter<active>>  _kDT;
    std::vector<Parameter<active>>  _kDTT;
    std::vector<Parameter<active>>  _kDTTT;

    std::vector<Parameter<active>>  _qMax;
    std::vector<Parameter<active>>  _qMaxT;
    std::vector<Parameter<active>>  _qMaxTT;
    std::vector<Parameter<active>>  _qMaxTTT;

    std::vector<Parameter<active>>  _beta;
    std::vector<Parameter<active>>  _betaT;
    std::vector<Parameter<active>>  _betaTT;
    std::vector<Parameter<active>>  _betaTTT;

    std::vector<Parameter<active>>  _gamma;
    std::vector<Parameter<active>>  _gammaT;
    std::vector<Parameter<active>>  _gammaTT;
    std::vector<Parameter<active>>  _gammaTTT;
};

} // namespace cadet

#endif // ADSORPTIONMODEL_EXTMPM_HPP_
