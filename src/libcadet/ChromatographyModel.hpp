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

#ifndef CHROMATOGRAPHYMODEL_HPP_
#define CHROMATOGRAPHYMODEL_HPP_

#include "ParameterContainer.hpp"
#include "Timer.hpp"

namespace cadet {

/// \brief Abstract base class for a chromatography model
class ChromatographyModel : public ParameterContainer<active>
{
public:
    ChromatographyModel(const SimulatorPImpl& sim, ChromatographyType chromType) :
        ParameterContainer<active>(sim),
        _chromType(chromType),
        _am(sim.getAdsorptionModel()),
        _ti(sim.getTimeIntegrator()),
        _pd(sim.getParticleDiscretization())
    {
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    virtual ~ChromatographyModel()
    {
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    virtual int residualDae(double t, N_Vector NV_y, N_Vector NV_yDot, N_Vector NV_res, void* userData) = 0;

    virtual int residualSens(int ns, double t, N_Vector NV_y, N_Vector NV_yDot, N_Vector NV_res,
            N_Vector* NV_yS, N_Vector* NV_ySDot, N_Vector* NV_resS,
            void* user_data, N_Vector NV_tmp1, N_Vector NV_tmp2, N_Vector NV_tmp3) = 0;


    //
    //     \brief Computes consistent initial values for the original DAE system.
    //
    //     ydot is adjusted so as to satisfy the DAE residual \f$F(t, y, ydot)=0\f$ for given t and y.
    //
    //     This affects only column values because the boundary conditions
    //     can only affect the first (and last) column cell.
    //
    //     Example usage:
    //
    //          //=========================================================================
    //          // Compute consistent initial values
    //          //=========================================================================
    //          calc_ic(data, data->t);
    //          calc_ic_sens(data, data->t);
    //          //=========================================================================
    //
    //
    //     \param[in,out]   data    A pointer to the master data structure of the simulator.
    //     \param[in]       t       The time at which the DAE system should be consistent.
    //     \see             calc_ic_sens
    //
    //     \return 0 if no errors occure
    //
    virtual void calcIC(const double t) = 0;

    //  calc_ic_sens
    //  computes consistent initial conditions for the sensitivity systems.
    //
    //  At the beginning of the integration, s and sdot must satisfy the system:
    //       dF/dy * s  +  dF/dydot * sdot  +  dF/dp  =  0
    //
    //  As s is fixed, we must find the respective values for sdot:
    //       sdot =  - (dF/dydot)^-1  *  (dF/dy * s  +  dF/dp)
    //
    //  As for the original DAE system, only the column part needs to be accounted for.
    //
    virtual void calcICSens(const double t) = 0;

    virtual void specialSetup() { }
    virtual void sectionSetup(int section, double t) { }

    // Timer read functions
    virtual double timerResDae()     const { return _timerResDae.getTime(); }
    virtual double timerResSens()    const { return _timerResSens.getTime(); }
    virtual double timerResDaePar()  const { return _timerResDaePar.getTime(); }
    virtual double timerResSensPar() const { return _timerResSensPar.getTime(); }

    // sets the section dependence of a parameter (group)
    virtual void setParameterSectionDependent(const ParameterName id, bool depends) { }

    virtual void setMultipleBoundStatesMode(int mode) { }
    virtual int getMultipleBoundStatesMode() const { return 0; }

protected:

    ChromatographyType _chromType;

    const AdsorptionModel&          _am;
    const TimeIntegrator&           _ti;
    const ParticleDiscretization&   _pd;

    OmpTimer                _timerResDae;           //!< OpenMP timer for evaluation of DAE residual
    OmpTimer                _timerResSens;          //!< OpenMP timer for evaluation of Sensitivity residual
    // Timer for parallel regions
    OmpTimer                _timerResDaePar;        //!< OpenMP timer for evaluation of DAE residual, parallel part
    OmpTimer                _timerResSensPar;       //!< OpenMP timer for evaluation of Sensitivity residual, parallel part

};

} // namespace cadet

#endif /* CHROMATOGRAPHYMODEL_HPP_ */
