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

#ifndef ADSORPTIONMODEL_MPM_HPP_
#define ADSORPTIONMODEL_MPM_HPP_

#include "AdsorptionModel.hpp"

namespace cadet
{

/// \brief Implementation of the Mobile Phyase Modulators adsorption model
/// All parameter-related functions are inherited from the ParameterContainer class
class AdsorptionModel_MPM : public AdsorptionModel
{
public:

    // Constructor
    AdsorptionModel_MPM(const SimulatorPImpl& sim) :
        AdsorptionModel(sim, MOBILE_PHASE_MODULATORS)
    {
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        double inf = std::numeric_limits<double>::infinity();

        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            addParam(Parameter<active> (MPM_KA,    e2s(MPM_KA),    comp, -1, 0.0, 0.0, 0.0, false, inf, true));
            addParam(Parameter<active> (MPM_KD,    e2s(MPM_KD),    comp, -1, 0.0, 0.0, 0.0, false, inf, true));
            addParam(Parameter<active> (MPM_QMAX,  e2s(MPM_QMAX),  comp, -1, 0.0, 0.0, 0.0, false, inf, true));
            addParam(Parameter<active> (MPM_BETA,  e2s(MPM_BETA),  comp, -1, 0.0, 0.0, 0.0, false, inf, true));
            addParam(Parameter<active> (MPM_GAMMA, e2s(MPM_GAMMA), comp, -1, 0.0, 0.0, 0.0, false, inf, true));
        }

        log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    // Destructor
    virtual ~AdsorptionModel_MPM()
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
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        double              ka    = getValue<double>           (MPM_KA, comp);
        double              kd    = getValue<double>           (MPM_KD, comp);
        std::vector<double> qmax  = getValueForAllComp<double> (MPM_QMAX);
        double              beta  = getValue<double>           (MPM_BETA, comp);
        double              gamma = getValue<double>           (MPM_GAMMA, comp);

        // Only protein components
        if (comp > 0)
        {
            // Liquid phase concentration
            const double* c = q - _cc.ncomp();

            // Liquid phase salt concentration
            double c0 = c[-comp];

            double ka_mpm = ka * exp(gamma * c0);
            double kd_mpm = kd * pow(c0, beta);

            double qsum = 1.0;
            for (int j = 1; j < _cc.ncomp(); ++j)
            {
                qsum -= q[-comp + j] / qmax.at(j);
            }

            // Jacobian
            for (int j = 1; j < _cc.ncomp(); ++j)
                jac[-comp + j] = ka_mpm * *c * qmax.at(comp) / qmax.at(j);              // dresi/dqj

            jac[0]            += kd_mpm;                                                // dresi/dqi
            jac[-_cc.ncomp()]  = -ka_mpm * qmax.at(comp) * qsum;                        // dresi/dci

            jac[-_cc.ncomp() - comp] = -ka_mpm * *c * qmax.at(comp) * qsum * gamma
                    + kd * beta * *q * pow(c0, beta - 1);                               // dresi/dc0
        }
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

private:
    template<typename StateType, typename ResidType, typename ParamType>
    void evaluateResidual(const int comp, const StateType *q, ResidType *res) const
    {
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        ParamType              ka    = getValue<ParamType>           (MPM_KA, comp);
        ParamType              kd    = getValue<ParamType>           (MPM_KD, comp);
        std::vector<ParamType> qmax  = getValueForAllComp<ParamType> (MPM_QMAX);
        ParamType              beta  = getValue<ParamType>           (MPM_BETA, comp);
        ParamType              gamma = getValue<ParamType>           (MPM_GAMMA, comp);

        // Liquid phase concentration
        const StateType* c = q - _cc.ncomp();

        if (comp == 0)
        {
            // Salt component
            *res = (ResidType) 0.0;
        }
        else
        {
            // Protein components
            
            // Liquid phase salt concentration
            StateType c0 = c[-comp];
            
            ResidType ka_mpm = ka * exp(gamma * c0);
            ResidType kd_mpm = kd * pow(c0, beta);

            ResidType qsum = 1.0;
            for (int j = 1; j < _cc.ncomp(); ++j)
            {
                qsum -= q[-comp + j] / qmax.at(j);
            }

            // Residual
            *res = - (ka_mpm * *c * qmax.at(comp) * qsum - kd_mpm * *q);
        }

        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }
};

} // namespace cadet

#endif // ADSORPTIONMODEL_MPM_HPP_
