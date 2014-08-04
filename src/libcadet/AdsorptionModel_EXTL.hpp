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

#ifndef ADSORPTIONMODEL_EXTL_HPP_
#define ADSORPTIONMODEL_EXTL_HPP_

#include "AdsorptionModel.hpp"

namespace cadet
{

/// \brief Implementation of the Multi Component Langmuir adsorption model with externally dependent parameters
/// All parameter-related functions are inherited from the ParameterContainer class
class AdsorptionModel_EXTL : public AdsorptionModel
{
public:

    // Constructor
    AdsorptionModel_EXTL(const SimulatorPImpl& sim) :
        AdsorptionModel(sim, EXTERNAL_LANGMUIR)
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

        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            _kA.push_back(Parameter<active> (EXTL_KA,       e2s(EXTL_KA),       comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kA[comp]);
            _kAT.push_back(Parameter<active> (EXTL_KA_T,     e2s(EXTL_KA_T),     comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kAT[comp]);
            _kATT.push_back(Parameter<active> (EXTL_KA_TT,    e2s(EXTL_KA_TT),    comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kATT[comp]);
            _kATTT.push_back(Parameter<active> (EXTL_KA_TTT,   e2s(EXTL_KA_TTT),   comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kATTT[comp]);
            _kD.push_back(Parameter<active> (EXTL_KD,       e2s(EXTL_KD),       comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kD[comp]);
            _kDT.push_back(Parameter<active> (EXTL_KD_T,     e2s(EXTL_KD_T),     comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kDT[comp]);
            _kDTT.push_back(Parameter<active> (EXTL_KD_TT,    e2s(EXTL_KD_TT),    comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kDTT[comp]);
            _kDTTT.push_back(Parameter<active> (EXTL_KD_TTT,   e2s(EXTL_KD_TTT),   comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_kDTTT[comp]);
            _qMax.push_back(Parameter<active> (EXTL_QMAX,     e2s(EXTL_QMAX),     comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_qMax[comp]);
            _qMaxT.push_back(Parameter<active> (EXTL_QMAX_T,   e2s(EXTL_QMAX_T),   comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_qMaxT[comp]);
            _qMaxTT.push_back(Parameter<active> (EXTL_QMAX_TT,  e2s(EXTL_QMAX_TT),  comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_qMaxTT[comp]);
            _qMaxTTT.push_back(Parameter<active> (EXTL_QMAX_TTT, e2s(EXTL_QMAX_TTT), comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(_qMaxTTT[comp]);
        }

        log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    // Destructor
    virtual ~AdsorptionModel_EXTL()
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
        const double              ka       = _kA[comp].getValue<double>();
        const double              ka_T     = _kAT[comp].getValue<double>();
        const double              ka_TT    = _kATT[comp].getValue<double>();
        const double              ka_TTT   = _kATTT[comp].getValue<double>();
        const double              kd       = _kD[comp].getValue<double>();
        const double              kd_T     = _kDT[comp].getValue<double>();
        const double              kd_TT    = _kDTT[comp].getValue<double>();
        const double              kd_TTT   = _kDTTT[comp].getValue<double>();

        // Liquid phase concentration
        const double* c = q -_cc.ncomp();

        double temp;
        _externalBase->externalProfile(z, t, &temp);

        const double temp2 = pow(temp, 2); // T^2: we use it again
        const double temp3 = pow(temp, 3); // T^3: we use it again

        // Compute overall parameters from their coefficients, depending on the externally given value
        const double ka_all = ka_TTT * temp3 + ka_TT * temp2 + ka_T * temp + ka;
        const double kd_all = kd_TTT * temp3 + kd_TT * temp2 + kd_T * temp + kd;
        std::vector<double> qmax_all(_cc.ncomp(), 0.0);
        for (int j = 0; j < _cc.ncomp(); ++j)
            qmax_all.at(j) = _qMaxTTT[j].getValue<double>() * temp3 + _qMaxTT[j].getValue<double>() * temp2 + _qMaxT[j].getValue<double>() * temp + _qMax[j].getValue<double>();

        double qsum = 1.0;
        for (int j = 0; j < _cc.ncomp(); ++j)
            qsum -= q[-comp + j] / qmax_all.at(j);

        // Jacobian
        for (int j = 0; j < _cc.ncomp(); ++j)
            jac[-comp + j] = ka_all * (*c) * qmax_all.at(comp) / qmax_all.at(j);  // dres/dq_j

        jac[0]            += kd_all;                                            // dres/dq_i
        jac[-_cc.ncomp()]  = -ka_all * qmax_all.at(comp) * qsum;                // dres/dc
    }

private:

    template<typename StateType, typename ResidType, typename ParamType>
    void evaluateResidual(const double t, const double z, const int comp, const StateType *q, ResidType *res) const
    {
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        const ParamType ka     = _kA[comp].getValue<ParamType>();
        const ParamType ka_T   = _kAT[comp].getValue<ParamType>();
        const ParamType ka_TT  = _kATT[comp].getValue<ParamType>();
        const ParamType ka_TTT = _kATTT[comp].getValue<ParamType>();
        const ParamType kd     = _kD[comp].getValue<ParamType>();
        const ParamType kd_T   = _kDT[comp].getValue<ParamType>();
        const ParamType kd_TT  = _kDTT[comp].getValue<ParamType>();
        const ParamType kd_TTT = _kDTTT[comp].getValue<ParamType>();

        // Liquid phase concentration
        const StateType* c = q -_cc.ncomp();

        double temp;
        _externalBase->externalProfile(z, t, &temp);

        ResidType temp2 = pow(temp, 2); // T^2: we use it again
        ResidType temp3 = pow(temp, 3); // T^3: we use it again

        // Compute overall parameters from their coefficients, depending on the externally given value
        ResidType ka_all = ka_TTT * temp3 + ka_TT * temp2 + ka_T * temp + ka;
        ResidType kd_all = kd_TTT * temp3 + kd_TT * temp2 + kd_T * temp + kd;
        std::vector<ResidType> qmax_all(_cc.ncomp(), 0.0);
        for (int j = 0; j < _cc.ncomp(); ++j)
            qmax_all.at(j) = _qMaxTTT[j].getValue<ParamType>() * temp3 + _qMaxTT[j].getValue<ParamType>() * temp2 + _qMaxT[j].getValue<ParamType>() * temp + _qMax[j].getValue<ParamType>();

        ResidType qsum = 1.0;
        for (int j = 0; j < _cc.ncomp(); ++j)
            qsum -= q[-comp + j] / qmax_all.at(j);

        // Residual
        *res = - (ka_all * (*c) * qmax_all.at(comp) * qsum - kd_all * *q);

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
};

} // namespace cadet

#endif // ADSORPTIONMODEL_EXTL_HPP_
