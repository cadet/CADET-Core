// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2016: Eric von Lieres¹, Joel Andersson¹,
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

#ifndef GENERALRATEMODEL_HPP_
#define GENERALRATEMODEL_HPP_

#include "ChromatographyModel.hpp"
#include "JacobianData.hpp"

namespace cadet {

#define CADET_STRICT true
#define CADET_LOOSE false

class GeneralRateModel : public ChromatographyModel
{
public:

    GeneralRateModel(SimulatorPImpl& sim);
    virtual ~GeneralRateModel();

    ///todo Check destructors and other stuff... rule of three, for all classes!

    int residualDae(double t, N_Vector NV_y, N_Vector NV_yDot, N_Vector NV_res, void* userData);
    int residualSens(int ns, double t, N_Vector NV_y, N_Vector NV_yDot, N_Vector NV_res,
            N_Vector* NV_yS, N_Vector* NV_ySDot, N_Vector* NV_resS,
            void* userData, N_Vector NV_tmp1, N_Vector NV_tmp2, N_Vector NV_tmp3);

    void calcIC(const double t);
    void calcICSens(const double t);

    void specialSetup();
    void sectionSetup(int section, double t);

    // sets the section dependence of a parameter (group)
    virtual void setParameterSectionDependent(const ParameterName id, bool depends);

    virtual void setMultipleBoundStatesMode(int mode) { _multiBoundMode = mode; }
    virtual int getMultipleBoundStatesMode() const { return _multiBoundMode; }

private:

    // this residual function now only handles column and particles
    // boundaries are treated differently elsewhere
    template <typename StateType, typename ResidType, typename ParamType, bool wantJac>
    int residualColumnParticle(const double t, const StateType* y, const double* ydot, ResidType* res) throw (CadetException);


    //       t     input       current time
    //       y     input       [array] pointer to first column element of state vector
    //      yp     input       [array] pointer to first column element of derivative state vector
    //       p     input       [scalar] pointer to parameter data structure containing all model parameters (sensisitve as well as non-sensitive)
    //     res    output       [array] calculated residual values
    //  csdata                 ChromsimData structure
    //
    template <typename StateType, typename ResidType, typename ParamType, bool wantJac>
    int residualColumn(const double t, const StateType* y, const double* ydot, ResidType* res) throw (CadetException);

    template <typename StateType, typename ResidType, typename ParamType, bool wantJac>
    int residualParticle(const double t, const int pblk, const StateType* y, const double* ydot, ResidType* res) throw (CadetException);

    template <typename ResidType, typename ParamType>
    int residualBoundaries(const double* y, const double* ydot, ResidType* res) throw (CadetException);

    void assembleOffdiagJac(int section, double t) throw (CadetException);

    template <typename ParamType>
    void setInletParamDerivatives(std::vector<ParamType>& concInlet);

    template <typename ParamType>
    ParamType* allocateArray(size_t n);

    template <typename ParamType>
    void destroyArray(ParamType* const ptr, size_t n);

    SimulatorPImpl&         _psim;
    const JacobianData&     _jac;
    const WenoScheme&       _ws;
    std::vector<double>     _c_in;
    std::vector<std::vector<double> > _dc_indp;

    bool        _secDepVelocity;
    bool        _secDepColDispersion;
    bool        _secDepParDiffusion;
    bool        _secDepParSurfDiffusion;
    bool        _secDepFilmDiffusion;
    int         _multiBoundMode;

    // Scalar parameters
    Parameter<active> _colLength;
    Parameter<active> _colPorosity;
    Parameter<active> _parRadius;
    Parameter<active> _parPorosity;

    // Section dependent parameters
    std::vector< Parameter<active> > _colDispersion;
    std::vector< Parameter<active> > _velocity;

    // Vectorial parameters
    std::vector< Parameter<active> > _filmDiffusion;
    std::vector< Parameter<active> > _parDiffusion;
    std::vector< Parameter<active> > _parSurfDiffusion;

    void dFdy_times_s(N_Vector NV_s, N_Vector NV_ret);
    void dFdyDot_times_sDot(N_Vector NV_sDot, N_Vector NV_ret);

    template <typename T>
    T getVelocity() const;

    template <typename T>
    T getDispersion() const;
};

} // namespace cadet


#endif /* GENERALRATEMODEL_HPP_ */
