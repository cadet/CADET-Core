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

        for (int comp = 0; comp < _nRealComp; ++comp)
        {
            addParam(Parameter<active> (MCBL_KA1,   e2s(MCBL_KA1),   comp, -1, 0.0, 0.0, 0.0, false, inf, true));
            addParam(Parameter<active> (MCBL_KD1,   e2s(MCBL_KD1),   comp, -1, 0.0, 0.0, 0.0, false, inf, true));
            addParam(Parameter<active> (MCBL_QMAX1, e2s(MCBL_QMAX1), comp, -1, 0.0, 0.0, 0.0, true,  inf, true));
            addParam(Parameter<active> (MCBL_KA2,   e2s(MCBL_KA2),   comp, -1, 0.0, 0.0, 0.0, false, inf, true));
            addParam(Parameter<active> (MCBL_KD2,   e2s(MCBL_KD2),   comp, -1, 0.0, 0.0, 0.0, false, inf, true));
            addParam(Parameter<active> (MCBL_QMAX2, e2s(MCBL_QMAX2), comp, -1, 0.0, 0.0, 0.0, true,  inf, true));
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

        const double              ka1   = getValue<double>              (MCBL_KA1, clampedComp);
        const double              kd1   = getValue<double>              (MCBL_KD1, clampedComp);
        const std::vector<double> qmax1 = getValueForFirstNComp<double> (MCBL_QMAX1, _nRealComp);
        const double              ka2   = getValue<double>              (MCBL_KA2, clampedComp);
        const double              kd2   = getValue<double>              (MCBL_KD2, clampedComp);
        const std::vector<double> qmax2 = getValueForFirstNComp<double> (MCBL_QMAX2, _nRealComp);

        // Residual
        if (comp < _nRealComp)
        {
            double qsum = 1.0;
            for (int j = 0; j < _nRealComp; ++j)
                qsum -= q[-comp + j] / qmax1.at(j);

            // Liquid phase concentration
            const double* c = q - _cc.ncomp();

            const double factor = ka1 * (*c) * qmax1.at(comp);

            // dres / dq1
            for (int j = 0; j < _nRealComp; ++j)
                jac[-comp + j] = factor / qmax1.at(j);  // dres/dq1_j
            jac[0] += kd1; // last summand by product rule

            // dres / dc
            jac[-_cc.ncomp()] = -ka1 * qmax1.at(comp) * qsum;
        }
        else
        {
            double qsum = 1.0;
            for (int j = 0; j < _nRealComp; ++j)
                qsum -= q[-comp + j + _nRealComp] / qmax2.at(j);

            // Liquid phase concentration
            const double* c = q - _cc.ncomp() - _nRealComp;

            const double factor = ka2 * (*c) * qmax2.at(clampedComp);

            // dres / dq2
            for (int j = 0; j < _nRealComp; ++j)
                jac[-comp + j + _nRealComp] = factor / qmax2.at(j);  // dres/dq2_j
            jac[0] += kd2;

            // dres / dc
            jac[-_cc.ncomp() - _nRealComp] = -ka2 * qmax2.at(clampedComp) * qsum;
        }
    }

private:

    template<typename StateType, typename ResidType, typename ParamType>
    void evaluateResidual(const int comp, const StateType *q, ResidType *res) const
    {
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        const int clampedComp = (comp >= _nRealComp) ? comp - _nRealComp : comp;

        ParamType              ka1   = getValue<ParamType>              (MCBL_KA1, clampedComp);
        ParamType              kd1   = getValue<ParamType>              (MCBL_KD1, clampedComp);
        std::vector<ParamType> qmax1 = getValueForFirstNComp<ParamType> (MCBL_QMAX1, _nRealComp);
        ParamType              ka2   = getValue<ParamType>              (MCBL_KA2, clampedComp);
        ParamType              kd2   = getValue<ParamType>              (MCBL_KD2, clampedComp);
        std::vector<ParamType> qmax2 = getValueForFirstNComp<ParamType> (MCBL_QMAX2, _nRealComp);

        // Residual
        if (comp < _nRealComp)
        {
            ResidType qsum = 1.0;
            for (int j = 0; j < _nRealComp; ++j)
                qsum -= q[-comp + j] / qmax1.at(j);

            // Liquid phase concentration
            const StateType* c = q -_cc.ncomp();

            *res = - (ka1 * (*c) * qmax1.at(comp) * qsum - kd1 * (*q));
        }
        else
        {
            ResidType qsum = 1.0;
            for (int j = 0; j < _nRealComp; ++j)
                qsum -= q[-comp + j + _nRealComp] / qmax2.at(j);

            // Liquid phase concentration
            const StateType* c = q -_cc.ncomp() - _nRealComp;

            *res = - (ka2 * (*c) * qmax2.at(clampedComp) * qsum - kd2 * (*q));
        }

        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    int _nRealComp;
};

} // namespace cadet

#endif // ADSORPTIONMODEL_MCBL_HPP_
