// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2015: Eric von Lieres¹, Joel Andersson,
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

#ifndef ADSORPTIONMODEL_MCBL_HPP_
#define ADSORPTIONMODEL_MCBL_HPP_

#include "AdsorptionModel.hpp"

namespace cadet
{

/// \brief Implementation of the Multi Component Bi-Langmuir adsorption model
/// All parameter-related functions are inherited from the ParameterContainer class
class AdsorptionModel_MCBL : public AdsorptionModel
{
public:

    // Constructor
    AdsorptionModel_MCBL(const SimulatorPImpl& sim) :
        AdsorptionModel(sim, MULTI_COMPONENT_BILANGMUIR)
    {
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        _nRealComp = _cc.ncomp() / 2;
        const double inf = std::numeric_limits<double>::infinity();

        this->configure();
        log::emit<Debug1>() << CURRENT_FUNCTION << ": Configured" << log::endl;

        _kA1.reserve(_cc.ncomp());
        _kD1.reserve(_cc.ncomp());
        _qMax1.reserve(_cc.ncomp());

        _kA2.reserve(_cc.ncomp());
        _kD2.reserve(_cc.ncomp());
        _qMax2.reserve(_cc.ncomp());

        for (int comp = 0; comp < _nRealComp; ++comp)
        {
            _kA1.push_back(Parameter<active> (MCBL_KA1,   e2s(MCBL_KA1),   comp, -1, 0.0, 0.0, 0.0, false, inf, true));
            addParam(_kA1[comp]);
            _kD1.push_back(Parameter<active> (MCBL_KD1,   e2s(MCBL_KD1),   comp, -1, 0.0, 0.0, 0.0, false, inf, true));
            addParam(_kD1[comp]);
            _qMax1.push_back(Parameter<active> (MCBL_QMAX1, e2s(MCBL_QMAX1), comp, -1, 0.0, 0.0, 0.0, true,  inf, true));
            addParam(_qMax1[comp]);
            _kA2.push_back(Parameter<active> (MCBL_KA2,   e2s(MCBL_KA2),   comp, -1, 0.0, 0.0, 0.0, false, inf, true));
            addParam(_kA2[comp]);
            _kD2.push_back(Parameter<active> (MCBL_KD2,   e2s(MCBL_KD2),   comp, -1, 0.0, 0.0, 0.0, false, inf, true));
            addParam(_kD2[comp]);
            _qMax2.push_back(Parameter<active> (MCBL_QMAX2, e2s(MCBL_QMAX2), comp, -1, 0.0, 0.0, 0.0, true,  inf, true));
            addParam(_qMax2[comp]);
        }

        log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    // Destructor
    virtual ~AdsorptionModel_MCBL()
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
        const int clampedComp = (comp >= _nRealComp) ? comp - _nRealComp : comp;

        const double              ka1   = _kA1[clampedComp].getValue<double>();
        const double              kd1   = _kD1[clampedComp].getValue<double>();
        const double              ka2   = _kA2[clampedComp].getValue<double>();
        const double              kd2   = _kD2[clampedComp].getValue<double>();

        // Residual
        if (comp < _nRealComp)
        {
            double qsum = 1.0;
            for (int j = 0; j < _nRealComp; ++j)
                qsum -= q[-comp + j] / _qMax1[j].getValue<double>();

            // Liquid phase concentration
            const double* c = q - _cc.ncomp();

            const double factor = ka1 * (*c) * _qMax1[comp].getValue<double>();

            // dres / dq1
            for (int j = 0; j < _nRealComp; ++j)
                jac[-comp + j] = factor / _qMax1[j].getValue<double>();  // dres/dq1_j
            jac[0] += kd1; // last summand by product rule

            // dres / dc
            jac[-_cc.ncomp()] = -ka1 * _qMax1[comp].getValue<double>() * qsum;
        }
        else
        {
            double qsum = 1.0;
            for (int j = 0; j < _nRealComp; ++j)
                qsum -= q[-comp + j + _nRealComp] / _qMax2[j].getValue<double>();

            // Liquid phase concentration
            const double* c = q - _cc.ncomp() - _nRealComp;

            const double factor = ka2 * (*c) * _qMax2[clampedComp].getValue<double>();

            // dres / dq2
            for (int j = 0; j < _nRealComp; ++j)
                jac[-comp + j + _nRealComp] = factor / _qMax2[j].getValue<double>();  // dres/dq2_j
            jac[0] += kd2;

            // dres / dc
            jac[-_cc.ncomp() - _nRealComp] = -ka2 * _qMax2[clampedComp].getValue<double>() * qsum;
        }
    }

private:

    template<typename StateType, typename ResidType, typename ParamType>
    void evaluateResidual(const int comp, const StateType *q, ResidType *res) const
    {
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        const int clampedComp = (comp >= _nRealComp) ? comp - _nRealComp : comp;

        const ParamType              ka1   = _kA1[clampedComp].getValue<ParamType>();
        const ParamType              kd1   = _kD1[clampedComp].getValue<ParamType>();
        const ParamType              ka2   = _kA2[clampedComp].getValue<ParamType>();
        const ParamType              kd2   = _kD2[clampedComp].getValue<ParamType>();

        // Residual
        if (comp < _nRealComp)
        {
            ResidType qsum = 1.0;
            for (int j = 0; j < _nRealComp; ++j)
                qsum -= q[-comp + j] / _qMax1[j].getValue<ParamType>();

            // Liquid phase concentration
            const StateType* c = q -_cc.ncomp();

            *res = - (ka1 * (*c) * _qMax1[comp].getValue<ParamType>() * qsum - kd1 * (*q));
        }
        else
        {
            ResidType qsum = 1.0;
            for (int j = 0; j < _nRealComp; ++j)
                qsum -= q[-comp + j + _nRealComp] / _qMax2[j].getValue<ParamType>();

            // Liquid phase concentration
            const StateType* c = q -_cc.ncomp() - _nRealComp;

            *res = - (ka2 * (*c) * _qMax2[clampedComp].getValue<ParamType>() * qsum - kd2 * (*q));
        }

        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    int _nRealComp;

    std::vector<Parameter<active>>  _kA1;
    std::vector<Parameter<active>>  _kD1;
    std::vector<Parameter<active>>  _qMax1;
    std::vector<Parameter<active>>  _kA2;
    std::vector<Parameter<active>>  _kD2;
    std::vector<Parameter<active>>  _qMax2;    
};

} // namespace cadet

#endif // ADSORPTIONMODEL_MCBL_HPP_
