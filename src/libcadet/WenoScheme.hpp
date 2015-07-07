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

#ifndef WENOSCHEME_HPP_
#define WENOSCHEME_HPP_

#include <vector>
#include <sstream>

#include "SimulatorPImpl.hpp"

#include "CadetException.hpp"
#include "CadetConvenience.hpp"
#include "CadetLogger.hpp"
#include "MathUtil.hpp"
#include "active.hpp"

namespace cadet {

//! \brief This class implements parts of the WENO scheme,
//!        which is a higher order approximation scheme used for the first order derivative in the convective term.
class WenoScheme
{
public:

    //! Constructor
    WenoScheme(const SimulatorPImpl& sim);

    //!
    //! \param  [in]    wk      WENO order
    //! \param  [in]    comp    Chemical component
    //! \param  [in]    bnd     Boundary treatment
    //! \param  [in]    v       Pointer to first cell centered value used by stencil (array of 5 values)
    //! \param  [out]   vm      Reconstructed value (1 scalar value)
    //! \param  [out]   Dvm     Derivative of vm
    //! \param  [in]    work    Pointer to work space required for intermediate results
    //!
    //! \brief Reconstruction of the cell boundary concentrations and their derivatives
    //!        from multiple cell centered values using the WENO-scheme.
    //!
    //! \return void

    template <typename StateType, bool wantJac>
    void wenoReconstruct(int wk, int comp, int bnd, StateType *v, StateType *vm, double* Dvm, StateType *work) const throw (CadetException);

    void setWenoEpsilons(N_Vector NV_y, const std::vector<double>& sectionTimes, int sec);
    inline double getWenoEpsilon(int comp) const { return _weps.at(comp); }


    // Setter for parameters set by user
    inline void configure(int wenoOrder = 3, double wenoEps = 1e-6, int boundaryModel = 0)
    {
        setWenoOrder(wenoOrder);
        setWenoEps(wenoEps);
        setBoundaryModel(boundaryModel);
        _weps.assign(_weps.size(), _wenoEps);
    }

    inline void   setWenoEps(double wenoEps)        { _wenoEps = wenoEps; }
    inline double getWenoEps()                const { return _wenoEps; }

    inline void setWenoOrder(int wenoOrder)         { _wenoOrder = wenoOrder; }
    inline int  getWenoOrder()                const { return _wenoOrder; }

    inline void setBoundaryModel(int boundaryModel) { _boundaryModel = boundaryModel; }
    inline int  getBoundaryModel()            const { return _boundaryModel; }

private:


    // Parameters set by user
    double _wenoEps;
    int _wenoOrder;
    int _boundaryModel;

    static const double _wenoD2[2];
    static const double _wenoC2[2*2];
    static const double _wenoJbvv2[2*3*3];

    static const double _wenoD3[3];
    static const double _wenoC3[3*3];
    static const double _wenoJbvv3[3*5*5];

    const SimulatorPImpl& _sim;
    const CadetConstants& _cc;

    std::vector<double> _weps;

    // Constants
    const int _samplePointsInletConc; //!< Number of sample points, that are used to estimate the maximum inlet concentration in the setWenoEpsilons function

//    inline StateType& V(int r) { return v[k-1+r]; }
};



// Implementation part for template member functions

template <typename StateType, bool wantJac>
void WenoScheme::wenoReconstruct(int wk, int comp, int bnd, StateType* v, StateType* vm, double* Dvm, StateType* work) const throw (CadetException)
{
    log::emit<Trace2>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

#if defined(ACTIVE_SETFAD) || defined(ACTIVE_SFAD)
    using cadet::sqr;
    using sfad::sqr;
#elif defined(ACTIVE_ADOLC)
    using cadet::sqr;
#endif

    // quick return if we have order 1 (i.e. simple upwind)
    if (wk == 1)
    {
        *vm = *v;
        if (wantJac)
            *Dvm = 1.0;

        log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
        return;
    }

    int i;
    int j;
    int r;

    int sl = 2 * wk - 1; // stencil length

    StateType* w = v + wk - 1; // Stencil pointer bend for readability - shortcut for v[k-1+x] in switch (k)

    StateType* beta  = work;
    StateType* alpha = work + wk;
    StateType* omega = work + wk;
    StateType* vr    = work + 2*wk;  // Reconstructed values

    const double* d;
    const double* c;
    const double* Jbvv;

    // Calculate smoothness measures
    switch (wk)
    {
        case 2:
            beta[0] = sqr(w[1] - w[0]);
            beta[1] = sqr(w[0] - w[-1]);
            d = _wenoD2;
            c = _wenoC2;
            Jbvv = _wenoJbvv2;
            break;
        case 3:
            beta[0] = 13.0/12.0 * sqr(w[ 0] - 2.0 * w[ 1] + w[2]) + 0.25 * sqr(3.0 * w[ 0] - 4.0 * w[ 1] +       w[2]);
            beta[1] = 13.0/12.0 * sqr(w[-1] - 2.0 * w[ 0] + w[1]) + 0.25 * sqr(      w[-1] -       w[ 1]             );
            beta[2] = 13.0/12.0 * sqr(w[-2] - 2.0 * w[-1] + w[0]) + 0.25 * sqr(      w[-2] - 4.0 * w[-1] + 3.0 * w[0]);
            d = _wenoD3;
            c = _wenoC3;
            Jbvv = _wenoJbvv3;
            break;
        default:
            std::ostringstream ss;
            ss << "WenoScheme::wenoReconstruct(): Wrong weno order specified - accepted values [1-3]: " << wk;
            throw CadetException(ss.str());
            break;
    }

    // Add eps to avoid divide-by-zeros
    for (r = 0; r < wk; ++r)
        beta[r] += getWenoEpsilon(comp);

    // Calculate weights
    for (r = 0; r < wk; ++r)
        alpha[r] = d[r] / sqr(beta[r]);

    // Avoid boundaries
    if (bnd != 0)
    {
        if (bnd < 0) // beginning of interval
            for (r = 0; r < -bnd; ++r)
                alpha[wk - 1 - r] = 0;
        else
            // end of interval
            for (r = 0; r < bnd; ++r)
                alpha[r] = 0;
    }

    // Norm weights
    StateType alpha_sum = alpha[0];
    for (r = 1; r < wk; ++r)
        alpha_sum += alpha[r];
    for (r = 0; r < wk; ++r)
        omega[r] /= alpha_sum;

    // Calculate reconstructed values
    for (r = 0; r < wk; ++r)
    {
        vr[r] = 0;
        for (j = 0; j < wk; ++j)
            vr[r] += c[r + wk * j] * w[-r+j];
    }

    // Weighted sum
    *vm = 0;
    for (r = 0; r < wk; ++r)
        *vm += vr[r] * omega[r];

    // Jacobian
    if (wantJac)
    {
        // Dependencies
        // 1. Constant vr in (*)

        // Start with "d(vm)/d(omega)" = vr and
        // multiply with "d(omega)/d(alpha)" to get "d(vm)/d(alpha)"
        double dot = 0;
        for (r = 0; r < wk; ++r)
            dot += static_cast<double>(vr[r]) * static_cast<double>(omega[r]); //StateType(vr[r] * omega[r]);
        for (r = 0; r < wk; ++r)
            vr[r] = (vr[r] - dot) / alpha_sum;

        // Multiply with "d(alpha)/d(beta)" to get "d(vm)/d(beta)"
        for (r = 0; r < wk; ++r)
            vr[r] *= -2 * d[r] / pow(beta[r], 3.0);

        // Multiply with "d(beta)/d(v)" to get Dvm = "d(vm)/d(v)"
        for (j = 0; j < sl; ++j)
        {
            Dvm[j] = 0;
            for (r = 0; r < wk; ++r)
            {
                dot = 0;
                for (i = 0; i < sl; ++i)
                    dot += static_cast<double>(Jbvv[r + wk * j + wk * sl * i]) * static_cast<double>(v[i]); // StateType(Jbvv[r + wk * j + wk * sl * i] * v[i]);
                // To do: re-arange Jbvv to reduce cache misses !
                Dvm[j] += static_cast<double>(vr[r]) * dot; // StateType(vr[r] * dot);
            }
        }

        // 2. Constant omega[r] in (*)
        for (r = 0; r < wk; ++r)
            for (j = 0; j < wk; ++j)
                Dvm[wk - 1 + j - r] += static_cast<double>(omega[r]) * c[r + wk * j]; // StateType(omega[r] * c[r + wk * j]);

    }

    log::emit<Trace2>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}




}

#endif /* WENOSCHEME_HPP_ */
