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

#include <iomanip>
#include <limits>
#include <algorithm>

#include "GeneralRateModel.hpp"
#include "TimeIntegrator.hpp"
#include "JacobianData.hpp"
#include "SchurSolver.hpp"
#include "WenoScheme.hpp"
#include "ParticleDiscretization.hpp"
#include "AdsorptionModel.hpp"

namespace cadet {


GeneralRateModel::GeneralRateModel(SimulatorPImpl& sim) :
    ChromatographyModel(sim, GENERAL_RATE_MODEL),
    _psim(sim),
    _jac (sim.getSchurSolver().getJacobianData()),
    _ws  (sim.getWenoScheme()),
    _c_in(_cc.ncomp(), 0.0),
    _secDepVelocity(false),
    _secDepColDispersion(false),
    _secDepParDiffusion(false),
    _secDepParSurfDiffusion(false),
    _secDepFilmDiffusion(false),
    _multiBoundMode(0),
    _colLength(COL_LENGTH,     e2s(COL_LENGTH),     -1, -1, 0.0, 0.0, 0.0, CADET_STRICT, std::numeric_limits<double>::infinity(), CADET_STRICT),
    _colPorosity(COL_POROSITY,   e2s(COL_POROSITY),   -1, -1, 0.0, 0.0, 0.0, CADET_LOOSE,  1.0, CADET_LOOSE),
    _parRadius(PAR_RADIUS,     e2s(PAR_RADIUS),     -1, -1, 0.0, 0.0, 0.0, CADET_STRICT, std::numeric_limits<double>::infinity(), CADET_STRICT),
    _parPorosity(PAR_POROSITY,   e2s(PAR_POROSITY),   -1, -1, 0.0, 0.0, 0.0, CADET_STRICT, 1.0, CADET_STRICT)
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    double inf = std::numeric_limits<double>::infinity();

    // Section dependent parameters
    _colDispersion.reserve(_cc.nsec()+1);
    _velocity.reserve(_cc.nsec()+1);

    // Vectorial parameters
    _filmDiffusion.reserve((_cc.nsec()+1) * _cc.ncomp());
    _parDiffusion.reserve((_cc.nsec()+1) * _cc.ncomp());
    _parSurfDiffusion.reserve((_cc.nsec()+1) * _cc.ncomp());

    // Scalar parameters
    addParam(_colLength);
    addParam(_colPorosity);
    addParam(_parRadius);
    addParam(_parPorosity);

    // Section dependent parameters
    for (int sec = -1; sec < _cc.nsec(); ++sec)
    {
        _colDispersion.push_back(Parameter<active> (COL_DISPERSION, e2s(COL_DISPERSION), -1, sec, 0.0, 0.0, 0.0, CADET_LOOSE,  inf, CADET_STRICT));
        addParam(_colDispersion[_colDispersion.size() - 1]);
        _velocity.push_back(Parameter<active> (VELOCITY,       e2s(VELOCITY),       -1, sec, 0.0, 0.0, 0.0, CADET_LOOSE, inf, CADET_STRICT));
        addParam(_velocity[_velocity.size()-1]);

        // Vectorial parameters
        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            _filmDiffusion.push_back(Parameter<active> (FILM_DIFFUSION,    e2s(FILM_DIFFUSION),    comp, sec, 0.0, 0.0, 0.0, CADET_LOOSE, inf, CADET_STRICT));
            addParam(_filmDiffusion[_filmDiffusion.size()-1]);
            _parDiffusion.push_back(Parameter<active> (PAR_DIFFUSION,     e2s(PAR_DIFFUSION),     comp, sec, 0.0, 0.0, 0.0, CADET_LOOSE, inf, CADET_STRICT));
            addParam(_parDiffusion[_parDiffusion.size() - 1]);
            _parSurfDiffusion.push_back(Parameter<active> (PAR_SURFDIFFUSION, e2s(PAR_SURFDIFFUSION), comp, sec, 0.0, 0.0, 0.0, CADET_LOOSE, inf, CADET_STRICT));
            addParam(_parSurfDiffusion[_parSurfDiffusion.size() - 1]);
        }
    }

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


GeneralRateModel::~GeneralRateModel()
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}

///todo Checkl destructors and other stuff... rule of three, for all classes!

///todo why are all parameters active? shouldn't it be possible to use double params for computations?

int GeneralRateModel::residualDae(double t, N_Vector NV_y, N_Vector NV_yDot, N_Vector NV_res, void* userData)
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;
    _timerResDae.start();

    //============================================
    double* y    = NV_DATA_S(NV_y);
    double* yDot = NV_DATA_S(NV_yDot);
    double* res  = NV_DATA_S(NV_res);
    //============================================

    _timerResDaePar.start();

#ifndef VERIFY_ANALYTICAL_JAC
    // Residual evaluation including analytical jacobian computation
    if (_psim.getTimeIntegrator().useAnalyticJacobian())
#endif
        residualColumnParticle<double, double, double, true> (t, y, yDot, res);


#ifndef VERIFY_ANALYTICAL_JAC
    // Residual evaluation including AD jacobian computation
    else
#endif
    {
        // reinitialize actives
        for (int i = 0; i < _cc.neq(); ++i)
        {
            // Copy content of state vector (NV) to AD state vector, init directional derivatives with zero
            _ti.getYAd(i).setValue(y[i]);
            _ti.getResAd(i).setValue(0.0);

            for (int dir = 0; dir < _ti.getJacAdDirs(); ++dir)
                _ti.getResAd(i).setADValue(dir, 0.0);
        }

        residualColumnParticle<active, active, double, false> (t, _ti.getYAd(), yDot, _ti.getResAd());

        // copy residual values from AD residual vector to residual vector (NV)
        for (int i = 0; i < _cc.neq(); ++i)
            res[i] = _ti.getResAd(i).getValue();

#ifdef VERIFY_ANALYTICAL_JAC
        // Comparison of analytical an AD jacobian implementation
        _psim.getSchurSolver().getJacobianData().compareWithAd(_ti.getResAd(), _ti.getDiagDir());
#endif

        // Copy Jacobian entries from resAd to band Jacobian data structures...
        _psim.getSchurSolver().getJacobianData().setFromAd(_ti.getResAd(), _ti.getDiagDir());
    }

    _timerResDaePar.stop();

    // now take care of the boundaries (without any sensitivity computation)...
    residualBoundaries<double, double> (y, yDot, res);

    // mark jacobian factorization for next call of schurSolve
    _psim.getTimeIntegrator().setFactorizeJac(true);


    _timerResDae.stop();

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    return 0;
}




int GeneralRateModel::residualSens(int ns, double t, N_Vector NV_y, N_Vector NV_yDot, N_Vector NV_res,
        N_Vector* NV_yS, N_Vector* NV_ySDot, N_Vector* NV_resS,
        void* userData, N_Vector NV_tmp1, N_Vector NV_tmp2, N_Vector NV_tmp3)
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;
    _timerResSens.start();

    NV_tmp1 = _ti.getNvTemp1(); // stores result of (dF / dy) * s
    NV_tmp2 = _ti.getNvTemp2(); // stores result of (dF / dyDot) * sDot

    double* tmp1 = NV_DATA_S(NV_tmp1);
    double* tmp2 = NV_DATA_S(NV_tmp2);

    double* resS;

    double* y    = NV_DATA_S(NV_y);
    double* yDot = NV_DATA_S(NV_yDot);

    _timerResSensPar.start();

    residualColumnParticle<double, active, active, false> (t, y, yDot, _ti.getResAd());

    _timerResSensPar.stop();

    //=============================================================

    residualBoundaries<active, active> (y, yDot, _ti.getResAd());

    for (int param = 0; param < _ti.getNSensParams(); param++)
    {
        //=============================================================================
        // Directional derivative (dF / dy) * s
        dFdy_times_s(NV_yS[param], NV_tmp1);

        // Directional derivative (dF / dyDot) * sDot
        dFdyDot_times_sDot(NV_ySDot[param], NV_tmp2);

        resS = NV_DATA_S(NV_resS[param]);

        _timerResSensPar.start();

        // Complete sens residual is the sum:
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < _cc.neq(); i++)
            resS[i] = tmp1[i] + tmp2[i] + _ti.getResAd(i).getADValue(_ti.getJacAdDirs() + param);

        _timerResSensPar.stop();
        //=============================================================================
    }

    // do not factorize Jacobians at next call of schurSolve
    _psim.getTimeIntegrator().setFactorizeJac(false);

    _timerResSens.stop();
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    return 0;
}



void GeneralRateModel::calcIC(const double t)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    // Compute column residual without time derivative (->nullptr) and store this in ydot.
    // This means ydot will contain the sum of the convective and the dispersive flux term.
    residualColumn<double, double, double, false>(t, _ti.getY(), nullptr, _ti.getYDot());

    const double invBetaC              = 1.0 / _colPorosity.getValue<double>() - 1.0;
    const double surfaceToVolumeRatio  = 3.0 / _parRadius.getValue<double>();

    // ydot shall fulfill ydot = -column_residual - film_flux
    // ==> negate and subtract boundary contribution
    for (int eqc = 0; eqc < _cc.neq_col(); ++eqc)
    {
        _ti.getYDot()[eqc] *= -1.0;
        _ti.getYDot()[eqc] -= invBetaC * surfaceToVolumeRatio * _cc.offsetBnd(_ti.getNvY())[eqc]; // <- boundary contribution!
    }

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


void GeneralRateModel::calcICSens(const double t)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    // Quick return if we have no sensitivities activated
    if (_ti.getNSensParams() < 1) return;

    // call DAE residual to compute Jacobian dF/dy
    residualDae(t, _ti.getNvY(), _ti.getNvYDot(), _ti.getNvTemp1(), nullptr);

    // call residuals for sensitivities
    residualColumnParticle<double, active, active, false>(t, _ti.getY(), _ti.getYDot(), _ti.getResAd());
    residualBoundaries<active, active>(_ti.getY(), _ti.getYDot(), _ti.getResAd());

    for (int param = 0; param < _ti.getNSensParams(); ++param)
    {
        // compute Jacobian times sensitivity -> Js
        dFdy_times_s(_ti.getNvYS(param), _ti.getNvTemp1());

        // For the column, the matrix in front of sdot, dF/dydot, is simply the identity.
        // ==>   sdot = -J*s - dF/dp    if parameter has direct influence
        // ==>   sdot = -J*s            if parameter has no direct influence
        for (int eqc = 0; eqc < _cc.neq_col(); eqc++)
            _ti.getYSDot(param)[eqc] = - _ti.getTemp1()[eqc] - _ti.getResAd(eqc).getADValue(_ti.getJacAdDirs() + param);
    }

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


void GeneralRateModel::specialSetup()
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    _dc_indp = std::vector<std::vector<double> >(_cc.ncomp(), std::vector<double>(_ti.getMaxSensInletParams(), 0.0));

    assembleOffdiagJac(0, 0.0);

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


void GeneralRateModel::sectionSetup(int section, double t)
{
    // Boundary Jacobian blocks only change for section dependent film or particle diffusion coefficients
    if (_secDepFilmDiffusion || _secDepParDiffusion)
        assembleOffdiagJac(section, t);
}


// this residual function now only handles column and particles
// boundaries are treated differently elsewhere
template <typename StateType, typename ResidType, typename ParamType, bool wantJac>
int GeneralRateModel::residualColumnParticle(const double t, const StateType* y, const double* yDot, ResidType* res) throw (CadetException)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    #pragma omp parallel for schedule(static)
    for (int pblk = -1; pblk < _cc.npblk(); ++pblk)
    {
        if (pblk == -1)
            residualColumn<StateType, ResidType, ParamType, wantJac> (t, y, yDot, res);
        else
        {
            const StateType* par_y    = y    + _cc.neq_col() + pblk * _cc.neq_par(); ///todo use offset inline functions
            const double*    par_yDot = yDot + _cc.neq_col() + pblk * _cc.neq_par();
            ResidType*       par_res  = res  + _cc.neq_col() + pblk * _cc.neq_par();

            residualParticle<StateType, ResidType, ParamType, wantJac> (t, pblk, par_y, par_yDot, par_res);
        }
    }

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    return 0;
}


template <>
void GeneralRateModel::setInletParamDerivatives<double>(std::vector<double>& concInlet)
{
    // If ParamType is a double, no derivatives are needed
}

template <>
void GeneralRateModel::setInletParamDerivatives<active>(std::vector<active>& concInlet)
{
    if (_ti.getNSensInletParams() > 0)
    {
        int localParamIndex = 0;
        for (int param = 0; param < _ti.getMaxSensInletParams(); ++param)
        {
            if (_ti.getInletParamIsSensitive(param)) // Found a sensitive inlet parameter
            {
                for (int comp = 0; comp < _cc.ncomp(); ++comp)
                {
                    concInlet.at(comp).setADValue(_ti.getJacAdDirs() + _ti.getNSensModelParams() + localParamIndex,
                            _dc_indp.at(comp).at(param));
                }
                localParamIndex++;
            }
        }
    }
}

//       t     input       current time
//       y     input       [array] pointer to first column element of state vector
//      yp     input       [array] pointer to first column element of derivative state vector
//       p     input       [scalar] pointer to parameter data structure containing all model parameters (sensisitve as well as non-sensitive)
//     res    output       [array] calculated residual values
//  csdata                 ChromsimData structure
//
template <typename StateType, typename ResidType, typename ParamType, bool wantJac>
int GeneralRateModel::residualColumn(const double t, const StateType* y, const double* yDot, ResidType* res) throw (CadetException)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    const ParamType u   = getVelocity<ParamType>();
    const ParamType d_c = getDispersion<ParamType>();
    const ParamType h   = _colLength.getValue<ParamType>() / double(_cc.ncol());
    const ParamType h2  = h * h;

    // WENO help variables
    const int& WK   = _cc.max_wk();           // Max. WENO-order
    StateType* work = new StateType[3 * WK];  // required by "wenoReconstruct"
    StateType* v    = new StateType[2 * WK];  // Stencil space
    double*    Dvm  = new double[2 * WK - 1]; // Derivatives of vm

    StateType* w = v + WK;  // Stencil pointer bend for readability - shortcut for v[WK + x]

    double* jac;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    // Clear content of the inlet concentrations derivative matrix
    for (int i = 0; i < _cc.ncomp(); ++i)
        _dc_indp.at(i).assign(_ti.getMaxSensInletParams(), 0.0);

    // Evaluate the user-specified function for the inlet concentration
    _ti.inletConcentration(t, _ti.getSection(t), _c_in, _ti.getInletParamIsSensitive(), _dc_indp);

    // Copy content of double vector "_c_in" into a vector of ParamTypes
    std::vector<ParamType> concInlet(_c_in.begin(),_c_in.end());

    // In case ParamType is active, copy derivatives to the respective active type, else do nothing
    setInletParamDerivatives<ParamType>(concInlet);
    /////////////////////////////////////////////////////////////////////////////////////////////////


    for (int comp = 0; comp < _cc.ncomp(); ++comp)
    {
        // Add time derivative
        if (yDot != nullptr)
            for (int col = 0; col < _cc.ncol(); ++col)
                _cc.colC<ResidType> (res, col, comp) = _cc.colC<double> (yDot, col, comp);
        else
            for (int col = 0; col < _cc.ncol(); ++col)
                _cc.colC<ResidType> (res, col, comp) = 0.0;

        // Stencil:
        w[2] = _cc.colC<StateType> (y, 2, comp);
        w[1] = _cc.colC<StateType> (y, 1, comp);
        w[0] = _cc.colC<StateType> (y, 0, comp);

        StateType vm = 0.0; // reconstructed value

        for (int i = 0; i < 2 * WK - 1; ++i)
            Dvm[i] = 0.0;

        int wk = 0;  // Current WENO-order

        for (int col = 0; col < _cc.ncol(); ++col)
        {
            // Jacobian entries
            if (wantJac)
            {
                jac = _jac.getJacC(comp,col);
                // Initialize Jacobian to zero
                for (int i = -WK; i < WK; ++i)  // -3 <= i <= 2
                    jac[WK + i] = 0.0;
            }

            // Add dispersion
            if (col < _cc.ncol() - 1) // right side
            {
                _cc.colC<ResidType> (res, col, comp) -= d_c / h2 * (w[ 1] - w[0]);
                // Jacobian entries
                if (wantJac)
                {
                    jac[WK]     += ParamType(d_c / h2);
                    jac[WK + 1] -= ParamType(d_c / h2);
                }
            }

            if (col > 0) // left side
            {
                _cc.colC<ResidType> (res, col, comp) -= d_c / h2 * (w[-1] - w[0]);
                // Jacobian entries
                if (wantJac)
                {
                    jac[WK]     += ParamType(d_c / h2);
                    jac[WK - 1] -= ParamType(d_c / h2);
                }
            }

            //-----------------------------------------------------------------
            // Add convection through this cell's left face
            if (col == 0)
                // for the first cell, the concentration at its left face
                // is determined by the inflow concentration
                _cc.colC<ResidType> (res, col, comp) -= u / h * concInlet.at(comp);
            else
                // ... for all other cells use the reconstructed value ...
                // Remember that vm still contains the reconstructed value
                // of the previous cell's *right* face,
                // which is identical to this cell's left face!
                _cc.colC<ResidType> (res, col, comp) -= u / h * vm;
            // Jacobian entries
            if (wantJac)
            {
                for (int i = 0; i < 2 * wk - 1; ++i)
                    jac[WK - wk + i] -= ParamType(u / h * Dvm[i]);
            }
            //-----------------------------------------------------------------


            // Boundaries
            int bnd = 0;
            switch (_ws.getBoundaryModel())
            {
            case 0: // Lower WENO order
                // This very statement selects the max. weno order for the current column cell
                // wk = min(maxWKleft, maxWKright)
                wk = std::min(std::min(col + 1, _ws.getWenoOrder()), std::min(_cc.ncol() - col, _ws.getWenoOrder()));
                break;

            case 1: // Zero weights
                wk = _ws.getWenoOrder();
                if (col < wk - 1)
                    bnd = -(wk - 1 - col);
                else if (col > _cc.ncol() - wk)
                    bnd = _cc.ncol() - col;
                break;

            case 2: // Zero weights for p != 0
                if (col == 0)
                    wk = 1;
                else
                {
                    wk = _ws.getWenoOrder();
                    if (col < wk - 1)
                        bnd = -(wk - 1 - col);
                    else if (col > _cc.ncol() - wk)
                        bnd = _cc.ncol() - col;
                }
                break;

            case 3: // Large ghost points
                wk = _ws.getWenoOrder();
                if (col == 0)
                {
                    w[-1] = 1e20;
                    w[-2] = 1e50;
                }
                else if (col == _cc.ncol() - 2)
                    w[2] = 1e20;
                else if (col == _cc.ncol() - 1)
                    w[2] = 1e50;
                break;

            default:
                std::ostringstream ss;
                ss << "GeneralRateModel::residualColumn(): Wrong boundary model specified - accepted values [0-3]: " << wk;
                throw CadetException(ss.str());
                break;
            }

            // Reconstruct concentration on this cell's right face
            _ws.wenoReconstruct<StateType, wantJac>(wk, comp, bnd, w - wk + 1, &vm, Dvm, work);

            // Right side
            _cc.colC<ResidType> (res, col, comp) += u / h * vm;
            // Jacobian entries
            if (wantJac)
                for (int i = 0; i < 2 * wk - 1; ++i)
                    jac[WK - wk + i + 1] += ParamType(u / h * Dvm[i]);

            // Update stencil
            w[-3] = w[-2];
            w[-2] = w[-1];
            w[-1] = w[ 0];
            w[ 0] = w[ 1];
            w[ 1] = w[ 2];
            w[ 2] = _cc.colC<StateType> (y, col + 3, comp);
        }
    }

    delete [] work;
    delete [] Dvm;
    delete [] v;

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    return 0;
}


template <>
int GeneralRateModel::residualParticle<active, double, double, true>(const double t, const int pblk, const active* y, const double* ydot, double* res) throw (CadetException)
{
    throw CadetException("You cannot compute sensitivities and analytical jacobian!");
}

template <>
int GeneralRateModel::residualParticle<active, double, active, true>(const double t, const int pblk, const active* y, const double* ydot, double* res) throw (CadetException)
{
    throw CadetException("You cannot compute sensitivities and analytical jacobian!");
}

template <>
int GeneralRateModel::residualParticle<active, active, double, true>(const double t, const int pblk, const active* y, const double* ydot, active* res) throw (CadetException)
{
    throw CadetException("You cannot compute sensitivities and analytical jacobian!");
}

template <>
int GeneralRateModel::residualParticle<active, active, active, true>(const double t, const int pblk, const active* y, const double* ydot, active* res) throw (CadetException)
{
    throw CadetException("You cannot compute sensitivities and analytical jacobian!");
}


template <typename StateType, typename ResidType, typename ParamType, bool wantJac>
int GeneralRateModel::residualParticle(const double t, const int pblk, const StateType* y, const double* ydot, ResidType* res) throw (CadetException)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    // Simulation parameters
    const ParamType radius            = _parRadius.getValue<ParamType>();
    const ParamType inv_beta_p        = 1.0 / _parPorosity.getValue<ParamType>() - 1.0;

    int parDiffOffset = 0;
    if (_secDepParDiffusion)
        parDiffOffset = (_ti.getCurrentSection()+1) * _cc.ncomp();

    int parSurfDiffOffset = 0;
    if (_secDepParSurfDiffusion)
        parSurfDiffOffset = (_ti.getCurrentSection()+1) * _cc.ncomp();

    double z = 1.0 / double(_cc.npblk()) * (0.5 + pblk) ; // Current z coordinate in the column - needed in externally dependent adsorption kinetic
    ParamType dr;

    ParamType pTypeInfo; // This is only to determine the residual function that is called

    // Add ku_par such that jac[0] is always the main diagonal element and negative indices are well-defined
    double* jac = _jac.getJacP(pblk) + _jac.ku_par();

    // Loop over particle cells
    for (int par = 0; par < _cc.npar(); ++par)
    {
        // Geometry
        const ParamType outer_area_per_volume = _pd.getOuterSurfAreaPerVolume(par) / radius;
        const ParamType inner_area_per_volume = _pd.getInnerSurfAreaPerVolume(par) / radius;

        // Mobile phase
        for (int comp = 0; comp < _cc.ncomp(); ++comp, ++res, ++y, ++ydot, jac += _jac.ld_jp())
        {
            // Add time derivatives

            if (_multiBoundMode == 1)
            {
                // Add sum of dq1 / dt + dq2 / dt to the residual of the first half of components
                // Only transport second (virtual) half of components
                if (2 * comp < _cc.ncomp())
                    *res = *ydot + inv_beta_p * (ydot[_cc.ncomp()] + ydot[_cc.ncomp() + _cc.ncomp() / 2]);
                else
                    *res = *ydot;
            }
            else if (_multiBoundMode == 2)
            {
                // Add sum of dq1 / dt + dq2 / dt to the residual of the first half of components, but watch out for salt
                // Only transport second (virtual) half of components
                if (comp == 0)
                    *res = *ydot + inv_beta_p * ydot[_cc.ncomp()]; // Salt has no second bound state
                else if (2 * comp < _cc.ncomp())
                    *res = *ydot + inv_beta_p * (ydot[_cc.ncomp()] + ydot[_cc.ncomp() + (_cc.ncomp() - 1) / 2]);
                else
                    *res = *ydot;
            }
            else
            {
                *res = *ydot + inv_beta_p * ydot[_cc.ncomp()];
            }

            // set dres / dc_i and dres / dq_i = 0
            if (wantJac)
            {
                jac[0]           = 0.0;
                jac[_cc.ncomp()] = 0.0;
            }

            const ParamType dp = _parDiffusion.at(comp + parDiffOffset).getValue<ParamType>();
            const ParamType ds = _parSurfDiffusion.at(comp + parSurfDiffOffset).getValue<ParamType>();

            // Add flow through outer surface
            if (par != 0)
            {
                // difference between two cell-centers
                dr = (_pd.getParCellCoords().at(par - 1) - _pd.getParCellCoords().at(par)) * radius;

                // Gradient approximation
                ResidType grad_c = (y[-2 * _cc.ncomp()] - y[0]) / dr;
                ResidType grad_q = (y[-_cc.ncomp()] - y[_cc.ncomp()]) / dr;

                // Molecular diffusion contribution
                *res -= outer_area_per_volume * dp * grad_c;

                // Surface diffusion contribution
                *res -= outer_area_per_volume * ds * inv_beta_p * grad_q;

                if (wantJac)
                {
                        // Jacobian entrys w.r.t. this cell's concentrations
                        jac[0]           += ParamType(outer_area_per_volume * dp / dr);                      // dres / dc_p,i^(p,j)
                        jac[_cc.ncomp()] += ParamType(outer_area_per_volume * inv_beta_p * ds / dr);         // dres / dq_i^(p,j)

                        // Jacobian entrys w.r.t. neighboring cell's concentrations
                        jac[-2 * _cc.ncomp()] = ParamType(-outer_area_per_volume * dp / dr);                 // dres / dc_p,i^(p,j-1)
                        jac[-_cc.ncomp()]     = ParamType(-outer_area_per_volume * inv_beta_p * ds / dr);    // dres / dq_i^(p,j-1)
                }
//                else
//                    throw CadetException("You cannot compute sensitivities and analytical jacobian!");
            }

            // Add flow through inner surface
            if (par != _cc.npar() - 1)
            {
                // difference between two cell-centers
                dr = (_pd.getParCellCoords().at(par) - _pd.getParCellCoords().at(par + 1)) * radius;

                // Gradient approximation
                ResidType grad_c = (y[0] - y[2 * _cc.ncomp()]) / dr;
                ResidType grad_q = (y[_cc.ncomp()] - y[3 * _cc.ncomp()]) / dr;

                // Molecular diffusion contribution
                *res += inner_area_per_volume * dp * grad_c;

                // Surface diffusion contribution
                *res += inner_area_per_volume * ds * inv_beta_p * grad_q;

                if (wantJac)
                {
                        // Jacobian entrys w.r.t. this cell's concentrations
                        jac[0]           += ParamType(inner_area_per_volume * dp / dr);                    // dres / dc_p,i^(p,j)
                        jac[_cc.ncomp()] += ParamType(inner_area_per_volume * inv_beta_p * ds / dr);       // dres / dq_i^(p,j)

                        // Jacobian entrys w.r.t. neighboring cell's concentrations
                        jac[2 * _cc.ncomp()] = ParamType(-inner_area_per_volume * dp / dr);                // dres / dc_p,i^(p,j+1)
                        jac[3 * _cc.ncomp()] = ParamType(-inner_area_per_volume * inv_beta_p * ds / dr);   // dres / dq_i^(p,j+1)
                }
//                still a bad hack with the operator double() !!! ///todo remove this bad hack
//                else
//                    throw CadetException("You cannot compute sensitivities and analytical jacobian!");
            }
        }

        // Bound phase
        for (int comp = 0; comp < _cc.ncomp(); ++comp, ++res, ++y, ++ydot, jac += _jac.ld_jp())
        {
            // === call the generic isotherm residual function =============
            _am.evaluateResidual(t, z, comp, y, res, &pTypeInfo);
            if (wantJac)
            {
                // static_cast should be sufficient here, but this statement is compiled even when wantJac = false
                _am.setJacobian(t, z, comp, reinterpret_cast<const double*>(y), jac);
            }

            // === Add time derivative for adsorption models ===============
            if (_am.isDifferential(comp)) *res += *ydot;
        }
    }

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    return 0;
}



template <typename ResidType, typename ParamType>
int GeneralRateModel::residualBoundaries(const double* y, const double* yDot, ResidType* res) throw (CadetException)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    const ParamType invBetaC          = 1.0 / _colPorosity.getValue<ParamType>() - 1.0;
    const ParamType epsP              = _parPorosity.getValue<ParamType>();
    const ParamType radius            = _parRadius.getValue<ParamType>();

    int filmDiffOffset = 0;
    if (_secDepFilmDiffusion)
        filmDiffOffset = (_ti.getCurrentSection()+1) * _cc.ncomp();

    int parDiffOffset = 0;
    if (_secDepParDiffusion)
        parDiffOffset = (_ti.getCurrentSection()+1) * _cc.ncomp();

    const ParamType surfaceToVolumeRatio = 3.0 / radius;
    const ParamType outerAreaPerVolume   = _pd.getOuterSurfAreaPerVolume(0) / radius;

    const ParamType jacCB_val = invBetaC * surfaceToVolumeRatio;
    const ParamType jacPB_val = -outerAreaPerVolume / epsP;

    ParamType* const kf_FV = new ParamType[_cc.ncomp()];  // kf for finite volumes

    const double relOuterShellHalfRadius = 0.5 * _pd.getCellSize(0);
    for (int comp = 0; comp < _cc.ncomp(); ++comp)
        kf_FV[comp] = 1.0 / (radius * relOuterShellHalfRadius / epsP / _parDiffusion[comp + parDiffOffset].getValue<ParamType>() + 1.0 / _filmDiffusion[comp + filmDiffOffset].getValue<ParamType>());

    int eq;  // Index for the current equation we work on

    ResidType* res_col   = res;
    ResidType* res_par   = res + _cc.neq_col();
    ResidType* res_bound = res + _cc.neq_col() + _cc.npblk() * _cc.neq_par();

    const double* y_col   = y;
    const double* y_par   = y + _cc.neq_col();
    const double* y_bound = y + _cc.neq_col() + _cc.npblk() * _cc.neq_par();

    //========================================================================
    // J_b part
    //========================================================================
    for (int comp = 0; comp < _cc.ncomp(); ++comp)
        for (int bnd = 0; bnd < _cc.nbnd(); ++bnd)
        {
            eq = bnd + comp * _cc.nbnd();
            res_bound[eq] = y_bound[eq];
        }
    //========================================================================


    //========================================================================
    // J_c,b part
    //========================================================================
    for (eq = 0; eq < _cc.neq_col(); ++eq)
        res_col[eq] += jacCB_val * y_bound[eq];
    //========================================================================


    //========================================================================
    // J_b,c part
    //========================================================================
    for (int bnd = 0; bnd < _cc.nbnd(); ++bnd)
        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            eq = bnd + comp * _cc.nbnd();
            res_bound[eq] += -kf_FV[comp] * y_col[eq];
        }
    //========================================================================


    //========================================================================
    // J_p,b part
    //========================================================================
    for (int pblk = 0; pblk < _cc.npblk(); ++pblk)
        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            eq = pblk + comp * _cc.npblk();
            res_par[pblk * _cc.neq_par() + comp] += jacPB_val * y_bound[eq];
        }
    //========================================================================


    //========================================================================
    // J_b,p part
    //========================================================================
    for (int pblk = 0; pblk < _cc.npblk(); pblk++)
        for (int comp = 0; comp < _cc.ncomp(); comp++)
        {
            eq = pblk + comp * _cc.npblk();
            res_bound[eq] += kf_FV[comp] * y_par[comp + pblk * _cc.neq_par()];
        }
    //========================================================================

    delete [] kf_FV;

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    return 0;
}




void GeneralRateModel::dFdy_times_s(N_Vector NV_s, N_Vector NV_ret)
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    N_VScale(1.0, NV_s, NV_ret);

    lapackInt_t     n;
    lapackInt_t     kl;
    lapackInt_t     ku;
    lapackInt_t     lda;

    double      alpha = 1.0;
    double      beta = 0.0;
    lapackInt_t inc = 1;

    char trans[] = "Trans";

    _timerResSensPar.start();

    #pragma omp parallel for private(n, kl, ku, lda) schedule(static)
    for (int pblk = -1; pblk < _cc.npblk(); ++pblk)
    {
        //==========================================================================
        // Column data
        //==========================================================================
        if (pblk == -1)
        {
            n   = _cc.neq_col();
            kl  = _jac.kl_col();
            ku  = _jac.ku_col();
            lda = _jac.ld_jc();

            DGBMV(trans, &n, &n, &kl, &ku, &alpha, _jac.getJacC(),
                    &lda, _cc.offsetCol(NV_s), &inc, &beta,
                    _cc.offsetCol(NV_ret), &inc);

            _jac.sparseMV(_jac.getJacCB(), _jac.numel_jcb(), 1.0, _cc.offsetBnd(NV_s), _cc.offsetCol(NV_ret));
        }
        //==========================================================================
        // Particle data
        //==========================================================================
        else
        {
            n   = _cc.neq_par();
            kl  = _jac.kl_par();
            ku  = _jac.ku_par();
            lda = _jac.ld_jp();

            DGBMV(trans, &n, &n, &kl, &ku, &alpha, _jac.getJacP(pblk),
                    &lda, _cc.offsetPar(NV_s, pblk), &inc, &beta,
                    _cc.offsetPar(NV_ret, pblk), &inc);

            _jac.sparseMV(_jac.getJacPB(pblk), _jac.numel_jpb(), 1.0, _cc.offsetBnd(NV_s), _cc.offsetPar(NV_ret, pblk));
        }
        //==========================================================================
    }
    _timerResSensPar.stop();

    //==========================================================================
    // Boundary data
    //==========================================================================
    _jac.sparseMV(_jac.getJacBC(), _jac.numel_jbc(), 1.0, _cc.offsetCol(NV_s), _cc.offsetBnd(NV_ret));

    for (int pblk = 0; pblk < _cc.npblk(); ++pblk)
        _jac.sparseMV(_jac.getJacBP(pblk), _jac.numel_jbp(), 1.0, _cc.offsetPar(NV_s, pblk), _cc.offsetBnd(NV_ret));
    //==========================================================================

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}



void GeneralRateModel::dFdyDot_times_sDot(N_Vector NV_sDot, N_Vector NV_ret)
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    const double invBetaP = 1.0 / _parPorosity.getValue<double>() - 1.0;

    _timerResSensPar.start();

    #pragma omp parallel for schedule(static)
    for (int pblk = -1; pblk < _cc.npblk(); ++pblk)
    {
        if (pblk == -1) // column
        {
            double* sDot = _cc.offsetCol(NV_sDot);
            double* ret  = _cc.offsetCol(NV_ret);

            for (int i = 0; i < _cc.neq_col(); ++i) // loop over column equations
                ret[i] = sDot[i];
        }
        else // particle
        {
            double* sDot = _cc.offsetPar(NV_sDot, pblk);
            double* ret  = _cc.offsetPar(NV_ret, pblk);

            for (int par = 0; par < _cc.npar(); ++par)
            {
                for (int comp = 0; comp < _cc.ncomp(); ++comp)
                {
                    // Add time derivatives
                    if (_multiBoundMode == 1)
                    {
                        // Add sum of dq1 / dt + dq2 / dt to the residual of the first half of components
                        // Only transport second (virtual) half of components
                        if (2 * comp < _cc.ncomp())
                            _cc.parC<double> (ret, par, comp) = _cc.parC<double> (sDot, par, comp)
                            + invBetaP * (_cc.parQ<double> (sDot, par, comp) + _cc.parQ<double> (sDot, par, comp + _cc.ncomp() / 2));
                        else
                            _cc.parC<double> (ret, par, comp) = _cc.parC<double> (sDot, par, comp);
                    }
                    else if (_multiBoundMode == 2)
                    {
                        // Add sum of dq1 / dt + dq2 / dt to the residual of the first half of components, but watch out for salt
                        // Only transport second (virtual) half of components
                        if (comp == 0)
                            // Salt has no second bound state
                            _cc.parC<double> (ret, par, comp) = _cc.parC<double> (sDot, par, comp)
                            + invBetaP * _cc.parQ<double> (sDot, par, comp);
                        else if (2 * comp < _cc.ncomp())
                            _cc.parC<double> (ret, par, comp) = _cc.parC<double> (sDot, par, comp)
                            + invBetaP * (_cc.parQ<double> (sDot, par, comp) + _cc.parQ<double> (sDot, par, comp + (_cc.ncomp() - 1) / 2));
                        else
                            _cc.parC<double> (ret, par, comp) = _cc.parC<double> (sDot, par, comp);
                    }
                    else
                    {
                        _cc.parC<double> (ret, par, comp) = _cc.parC<double> (sDot, par, comp)
                        + invBetaP * _cc.parQ<double> (sDot, par, comp);
                    }
                }
            }

            for (int par = 0; par < _cc.npar(); ++par)
                for (int comp = 0; comp < _cc.ncomp(); ++comp)
                    if (_am.isDifferential(comp)) ///todo might be faster like this: parQ = _am.isDifferential(comp) * parQ ...
                        _cc.parQ<double> (ret, par, comp) = _cc.parQ<double> (sDot, par, comp);
                    else
                        _cc.parQ<double> (ret, par, comp) = 0.0;
        }
    }
    _timerResSensPar.stop();

    double* dFdyDot = _cc.offsetBnd(NV_ret);
    for (int eqb = 0; eqb < _cc.neq_bnd(); ++eqb) // loop over boundary equations
        dFdyDot[eqb] = 0.0;

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


void GeneralRateModel::assembleOffdiagJac(int section, double t) throw (CadetException)
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    SparseMatrixElement* jcb;
    SparseMatrixElement* jbc;
    SparseMatrixElement* jpb;
    SparseMatrixElement* jbp;

    const double invBetaC          = 1.0 / _colPorosity.getValue<double>() - 1.0;
    const double epsP              = _parPorosity.getValue<double>();
    const double radius            = _parRadius.getValue<double>();

    int filmDiffOffset = 0;
    if (_secDepFilmDiffusion)
        filmDiffOffset = (_ti.getCurrentSection()+1) * _cc.ncomp();

    int parDiffOffset = 0;
    if (_secDepParDiffusion)
        parDiffOffset = (_ti.getCurrentSection()+1) * _cc.ncomp();

    log::emit<Debug2>() << CURRENT_FUNCTION << ": Cur sec " << _ti.getCurrentSection() << " FilmDiff offset " << filmDiffOffset << " ParDiff offset " << parDiffOffset << log::endl;

    const double surfaceToVolumeRatio = 3.0 / radius;
    const double outerAreaPerVolume   = _pd.getOuterSurfAreaPerVolume(0) / radius;

    log::emit<Debug2>() << CURRENT_FUNCTION << ": radius = " << radius << log::endl;
    log::emit<Debug2>() << CURRENT_FUNCTION << ": _pd.getOuterSurfAreaPerVolume(0) = " << _pd.getOuterSurfAreaPerVolume(0) << log::endl;

    log::emit<Debug2>() << CURRENT_FUNCTION << ": surfaceToVolumeRatio = " << surfaceToVolumeRatio << log::endl;
    log::emit<Debug2>() << CURRENT_FUNCTION << ": outerAreaPerVolume = " << outerAreaPerVolume << log::endl;

    const double jacCB_val = invBetaC * surfaceToVolumeRatio;
    const double jacPB_val = -outerAreaPerVolume / epsP;
    const double relOuterShellHalfRadius = 0.5 * _pd.getCellSize(0);

    log::emit<Debug2>() << CURRENT_FUNCTION << ": jacCB_val = " << jacCB_val << log::endl;
    log::emit<Debug2>() << CURRENT_FUNCTION << ": jacPB_val = " << jacPB_val << log::endl;

    double* const kf_FV = new double[_cc.ncomp()];  // kf for finite volumes ///todo compute and store kf_FV in the specialSetup routine!!!

    for (int comp = 0; comp < _cc.ncomp(); ++comp)
    {
        kf_FV[comp] = 1.0 / (relOuterShellHalfRadius * radius / epsP / _parDiffusion[comp + parDiffOffset].getValue<double>() + 1.0 / _filmDiffusion[comp + filmDiffOffset].getValue<double>());
        log::emit<Debug2>() << CURRENT_FUNCTION << ": kf_FV[" << comp << "] = " << kf_FV[comp] << log::endl;
    }

    int eq;  // Index for the current equation we work on

    //========================================================================
    // J_c,b part
    //========================================================================
    jcb = _jac.getJacCB();
    for (eq = 0; eq < _cc.neq_col(); ++eq)
    {
        jcb->setElement(eq, eq, jacCB_val);
        jcb++;
    }
    //========================================================================


    //========================================================================
    // J_b,c part
    //========================================================================
    jbc = _jac.getJacBC();
    for (int col = 0; col < _cc.ncol(); ++col)
    {
        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            eq = col + comp * _cc.ncol();
            jbc->setElement(eq, eq, -kf_FV[comp]);
            jbc++;
        }
    }
    //========================================================================



    //========================================================================
    // J_p,b part
    //========================================================================
    for (int pblk = 0; pblk < _cc.npblk(); ++pblk)
    {
        jpb = _jac.getJacPB(pblk);
        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            eq = pblk + comp * _cc.ncol();
            jpb->setElement(comp, eq, jacPB_val);
            jpb++;
        }
    }
    //========================================================================


    //========================================================================
    // J_b,p part
    //========================================================================
    for (int pblk = 0; pblk < _cc.npblk(); ++pblk)
    {
        jbp = _jac.getJacBP(pblk);
        for (int comp = 0; comp < _cc.ncomp(); ++comp)
        {
            eq = pblk + comp * _cc.ncol();
            jbp->setElement(eq, comp, kf_FV[comp]);
            jbp++;
        }
    }
    //========================================================================

    delete [] kf_FV;

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}


void GeneralRateModel::setParameterSectionDependent(const ParameterName id, bool depends)
{
    int startSec = -1;
    int endSec = 0;
    int startComp = 0;
    int endComp = _cc.ncomp();

    if (!depends)
    {
        startSec = 0;
        endSec = _cc.nsec();
    }

    switch(id)
    {
        case COL_DISPERSION:
            _secDepColDispersion = depends;
            startComp = -1;
            endComp = 0;
            log::emit<Debug2>() << CURRENT_FUNCTION << ": Setting section dependency of COL_DISPERSION to " << _secDepColDispersion << log::endl;
            break;
        case VELOCITY:
            _secDepVelocity = depends;
            startComp = -1;
            endComp = 0;
            log::emit<Debug2>() << CURRENT_FUNCTION << ": Setting section dependency of VELOCITY to " << _secDepVelocity << log::endl;
            break;
        case PAR_DIFFUSION:
            _secDepParDiffusion = depends;
            log::emit<Debug2>() << CURRENT_FUNCTION << ": Setting section dependency of PAR_DIFFUSION to " << _secDepParDiffusion << log::endl;
            break;    
        case PAR_SURFDIFFUSION:
            _secDepParSurfDiffusion = depends;
            log::emit<Debug2>() << CURRENT_FUNCTION << ": Setting section dependency of PAR_SURFDIFFUSION to " << _secDepParSurfDiffusion << log::endl;
            break;
        case FILM_DIFFUSION:
            _secDepFilmDiffusion = depends;
            log::emit<Debug2>() << CURRENT_FUNCTION << ": Setting section dependency of FILM_DIFFUSION to " << _secDepFilmDiffusion << log::endl;
            break;
        default:
            return;
    }

    // Remove superfluous parameters
    for (int sec = startSec; sec < endSec; ++sec)
    {
        for (int comp = startComp; comp < endComp; ++comp)
        {
            removeParam(id, comp, sec);
        }
    }
}

template <typename T>
T GeneralRateModel::getVelocity() const
{
    if (_secDepVelocity)
        return _velocity[_ti.getCurrentSection()+1].getValue<T>();
    else
        return _velocity[0].getValue<T>();
}

template <typename T>
T GeneralRateModel::getDispersion() const
{
    if (_secDepColDispersion)
        return _colDispersion[_ti.getCurrentSection()+1].getValue<T>();
    else
        return _colDispersion[0].getValue<T>();
}

} // namespace cadet
