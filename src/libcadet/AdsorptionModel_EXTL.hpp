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

        double inf = std::numeric_limits<double>::infinity();

        this->configure();
        log::emit<Debug1>() << CURRENT_FUNCTION << ": Configured" << log::endl;

        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            addParam(Parameter<active> (EXTL_KA,       e2s(EXTL_KA),       comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTL_KA_T,     e2s(EXTL_KA_T),     comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTL_KA_TT,    e2s(EXTL_KA_TT),    comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTL_KA_TTT,   e2s(EXTL_KA_TTT),   comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTL_KD,       e2s(EXTL_KD),       comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTL_KD_T,     e2s(EXTL_KD_T),     comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTL_KD_TT,    e2s(EXTL_KD_TT),    comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTL_KD_TTT,   e2s(EXTL_KD_TTT),   comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTL_QMAX,     e2s(EXTL_QMAX),     comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTL_QMAX_T,   e2s(EXTL_QMAX_T),   comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTL_QMAX_TT,  e2s(EXTL_QMAX_TT),  comp, -1, 0.0, 0.0, -inf, true, inf, true));
            addParam(Parameter<active> (EXTL_QMAX_TTT, e2s(EXTL_QMAX_TTT), comp, -1, 0.0, 0.0, -inf, true, inf, true));
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
        double              ka       = getValue<double>           (EXTL_KA,     comp);
        double              ka_T     = getValue<double>           (EXTL_KA_T,   comp);
        double              ka_TT    = getValue<double>           (EXTL_KA_TT,  comp);
        double              ka_TTT   = getValue<double>           (EXTL_KA_TTT, comp);
        double              kd       = getValue<double>           (EXTL_KD,     comp);
        double              kd_T     = getValue<double>           (EXTL_KD_T,   comp);
        double              kd_TT    = getValue<double>           (EXTL_KD_TT,  comp);
        double              kd_TTT   = getValue<double>           (EXTL_KD_TTT, comp);
        std::vector<double> qmax     = getValueForAllComp<double> (EXTL_QMAX);
        std::vector<double> qmax_T   = getValueForAllComp<double> (EXTL_QMAX_T);
        std::vector<double> qmax_TT  = getValueForAllComp<double> (EXTL_QMAX_TT);
        std::vector<double> qmax_TTT = getValueForAllComp<double> (EXTL_QMAX_TTT);

        // Liquid phase concentration
        const double* c = q -_cc.ncomp();

        double temp;
        _externalBase->externalProfile(z, t, &temp);

        double temp2 = pow(temp, 2); // T^2: we use it again
        double temp3 = pow(temp, 3); // T^3: we use it again

        // Compute overall parameters from their coefficients, depending on the externally given value
        double ka_all = ka_TTT * temp3 + ka_TT * temp2 + ka_T * temp + ka;
        double kd_all = kd_TTT * temp3 + kd_TT * temp2 + kd_T * temp + kd;
        std::vector<double> qmax_all(_cc.ncomp(), 0.0);
        for (int j = 0; j < _cc.ncomp(); ++j)
            qmax_all.at(j) = qmax_TTT.at(j) * temp3 + qmax_TT.at(j) * temp2 + qmax_T.at(j) * temp + qmax.at(j);

        double qsum = 1.0;
        for (int j = 0; j < _cc.ncomp(); ++j)
            qsum -= q[-comp + j] / qmax_all.at(j);

        // Jacobian
        for (int j = 0; j < _cc.ncomp(); ++j)
            jac[-comp + j] = ka_all * *c * qmax_all.at(comp) / qmax_all.at(j);  // dres/dq_j

        jac[0]            += kd_all;                                            // dres/dq_i
        jac[-_cc.ncomp()]  = -ka_all * qmax_all.at(comp) * qsum;                // dres/dc
    }

private:

    template<typename StateType, typename ResidType, typename ParamType>
    void evaluateResidual(const double t, const double z, const int comp, const StateType *q, ResidType *res) const
    {
        log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

        ParamType              ka       = getValue<ParamType>           (EXTL_KA,     comp);
        ParamType              ka_T     = getValue<ParamType>           (EXTL_KA_T,   comp);
        ParamType              ka_TT    = getValue<ParamType>           (EXTL_KA_TT,  comp);
        ParamType              ka_TTT   = getValue<ParamType>           (EXTL_KA_TTT, comp);
        ParamType              kd       = getValue<ParamType>           (EXTL_KD,     comp);
        ParamType              kd_T     = getValue<ParamType>           (EXTL_KD_T,   comp);
        ParamType              kd_TT    = getValue<ParamType>           (EXTL_KD_TT,  comp);
        ParamType              kd_TTT   = getValue<ParamType>           (EXTL_KD_TTT, comp);
        std::vector<ParamType> qmax     = getValueForAllComp<ParamType> (EXTL_QMAX);
        std::vector<ParamType> qmax_T   = getValueForAllComp<ParamType> (EXTL_QMAX_T);
        std::vector<ParamType> qmax_TT  = getValueForAllComp<ParamType> (EXTL_QMAX_TT);
        std::vector<ParamType> qmax_TTT = getValueForAllComp<ParamType> (EXTL_QMAX_TTT);

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
            qmax_all.at(j) = qmax_TTT.at(j) * temp3 + qmax_TT.at(j) * temp2 + qmax_T.at(j) * temp + qmax.at(j);

        ResidType qsum = 1.0;
        for (int j = 0; j < _cc.ncomp(); ++j)
            qsum -= q[-comp + j] / qmax_all.at(j);

        // Residual
        *res = - (ka_all * *c * qmax_all.at(comp) * qsum - kd_all * *q);

        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

};

} // namespace cadet

#endif // ADSORPTIONMODEL_EXTL_HPP_
