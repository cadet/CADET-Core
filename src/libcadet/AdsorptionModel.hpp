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

#ifndef ADSORPTIONMODEL_HPP_
#define ADSORPTIONMODEL_HPP_

#include <string>
#include <vector>
#include <limits>

#include "CadetLogger.hpp"
#include "ParameterContainer.hpp"

namespace cadet
{

/// \brief Abstract base class for an adsorption model implementation
class AdsorptionModel : public ParameterContainer<active>
{
public:

    // Constructor
    AdsorptionModel(const SimulatorPImpl& sim, AdsorptionType adsType) :
        ParameterContainer<active>(sim),
        _adsType(adsType),
        _isDifferential(std::vector<bool>(_cc.ncomp(), false))
    {
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    // Destructor
    virtual ~AdsorptionModel()
    {
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;
        log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
    }

    /// \brief Pure virtual function for residual evaluation
    // Has to be implemented in the inheriting children classes
    virtual void evaluateResidual(const double t, const double z, const int comp, const active *q, active *res, const active *p) const = 0;
    virtual void evaluateResidual(const double t, const double z, const int comp, const active *q, active *res, const double *p) const = 0;
    virtual void evaluateResidual(const double t, const double z, const int comp, const double *q, active *res, const active *p) const = 0;
    virtual void evaluateResidual(const double t, const double z, const int comp, const double *q, double *res, const double *p) const = 0;

    virtual void setJacobian(const double t, const double z, const int comp, const double* q, double* jac) const throw (CadetException) {
        throw CadetException("No analytic Jacobian available for the chosen adsorption model!"); }
//    virtual void setJacobian(const int comp, const active* q, double* jac) const throw (CadetException) {
//        throw CadetException("You cannot compute sensitivities with analytical jacobian!"); }

    virtual void setIsKinetic(bool isKinetic = false) = 0;

    inline bool isDifferential(int comp) const { return _isDifferential.at(comp); }
    inline const AdsorptionType& getAdsorptionType() const { return _adsType; }

    inline void setExternalProfile(ExternalBase* externalProfile) { _externalBase = externalProfile; }

    // Setter for parameters set by user
    inline void configure() { setIsKinetic(); }

protected:

    // Parameters set by user
    bool _isKinetic;

    AdsorptionType      _adsType;
    ExternalBase*       _externalBase;    //!< Pointer to a class interface for registering external profile with adsorption kinetics

    std::vector<bool>    _isDifferential;

};

} // namespace cadet

#endif // ADSORPTIONMODEL_HPP_
