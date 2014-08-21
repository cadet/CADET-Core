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

#ifndef LIBCADET_SIMULATOR_HPP_
#define LIBCADET_SIMULATOR_HPP_

//!
//! \mainpage The Chromatography Analysis and Design Toolkit (CADET)
//!
//! The library provides a fast and accurate forward simulator for the general rate model of liquid chromatography
//!
//! \authors    Eric von Lieres, Joel Andersson, Andreas Puettmann, Sebastian Schnittert, Samuel Leweke
//! \version    2.2.1
//!
//! \date       2008-2014
//! \copyright  GNU General Public License v3.0 (or, at your option, any later version).
//!

#include <string>
#include <vector>
#include <tuple>

#include "CadetException.hpp"
#include "CadetEnumeration.hpp"

// Export and import classes when using MS Visual Studio compiler
#ifndef CADET_API
    #ifdef _MSC_VER
        #if defined(libcadet_shared_EXPORTS) || defined(libcadet_mex_EXPORTS)
            #define CADET_API _declspec(dllexport)
        #else
            #define CADET_API _declspec(dllimport)
        #endif
    #else
        #define CADET_API
    #endif
#endif

namespace cadet
{

class SimulatorPImpl;  // Forward declaration of private implementation class

class CADET_API InletBase
{
public:
    virtual ~InletBase() {}
    virtual void inletConcentration(double t, int sec, std::vector<double>& inletConc,
                                    const std::vector<bool>& inletParamIsSensitive,
                                    std::vector<std::vector<double> >& dInletConc_dp) = 0;
};


class CADET_API ExternalBase
{
public:
    virtual ~ExternalBase() {}
    virtual void externalProfile(double z, double t, double* value) = 0;
};


class CADET_API QuadBase
{
    virtual ~QuadBase() {}
    virtual void quadFunction(const std::vector<double>& columnOutletSolution,
                              const std::vector<std::vector<double> >& columnOutletSens,
                              std::vector<double> integrands) = 0;
};

//! \brief This class represents the user interface to the simulator library.
//!        The user has only access to the public functions in this class.
class CADET_API Simulator
{
public:

    //! \brief Constructor for a Simulator object
    //!
    //! The constructor allocates a memory block for the Simulator and stores a pointer to that block in #_sim.
    //! For the construction of a Simulator object all of the parameters are mandatory.
    //! None of them can be changed later. If a change is desired, a new Simulator instance has to be created.
    //!
    //! \param  [in]    ncomp       The number of chemical species (components) taking part in the separation
    //! \param  [in]    ncol        The number of axial discretization cells in the column
    //! \param  [in]    npar        The number of radial discretization cells in the beads
    //! \param  [in]    nsec        The number of sections
    //! \param  [in]    adsType     The type of adsorption model used
    //! \param  [in]    chromType   The type of chromatography model used
    //!
    Simulator(int ncomp, int ncol, int npar, int nsec, AdsorptionType adsType, ChromatographyType chromType);


    //! \brief Destructor for a Simulator object
    //!
    //! The destructor deletes the memory block pointed to by #_sim
    //!
    ~Simulator();


    //! \brief Sets the value of a chromatography- or adsorption-model parameter
    //!
    //! The value of a parameter can be modified through this setter method.
    //! All potentially sensitive parameters (exept for inlet-parameters) are handled by this method.
    //! The parameter is identified by its enum name through \p id, its respective component \p comp,
    //! and its section \p sec.
    //!
    //!     double  newValue = 2e-3;
    //!     int     comp     = 2;
    //!     int     sec      = -1;
    //!     setParameterValue(newValue, MCL_KA, comp, sec);
    //!
    //! \param  [in]    value   The new value assigned to the parameter
    //! \param  [in]    id      The enum name of the parameter
    //! \param  [in]    comp    The chemical species (component) adressed by the parameter
    //! \param  [in]    sec     The section adressed by the parameter
    //!
    void setParameterValue(double value, const ParameterName id, int comp = -1, int sec = -1);


    //! \brief Returns the value of a chromatography- or adsorption-model parameter
    //!
    //!	The value of a parameter can be obtained through this getter method.
    //! All potentially sensitive parameters (exept for inlet-parameters) are handled by this method.
    //! The parameter is identified by its enum name through \p id and its respective component \p comp.
    //!
    //! \param  [in]    id      The enum name of the parameter
    //! \param  [in]    comp    The chemical species (component) adressed by the parameter
    //! \param  [in]    sec     The section adressed by the parameter
    //!
    //! \return The current value of the parameter
    //!
    double getParameterValue(const ParameterName id, int comp = -1, int sec = -1) const;


    //! \brief Sets a parameter sensitive
    //!
    //!	The method enables sensitivity computation for the given parameter.
    //! For each parameter marked as sensitive, another system in size of the original DAE system
    //! is to be solved.
    //!
    //! \param	[in]	id	    The enum name of the parameter
    //! \param	[in]	absTolS Absolute tolerance used in the sensitivity computation for this parameter
    //! \param  [in]    comp    The chemical species (component) adressed by the parameter
    //! \param  [in]    sec     The section adressed by the parameter
    //!
    void setSensitiveParameter(const ParameterName id, double absTolS = 1e-5, int comp = -1, int sec = -1);


    //! \brief Sets a parameter group section dependent
    //!
    //! Some parameters are not section dependent (e.g., the porosities in a General Rate Model).
    //! Other parameters may be section dependent. This function sets whether they are section dependent.
    //!
    //! \param  [in]    id      The enum name of the parameter
    //! \param  [in]    depends Decides whether this parameter (group) is section dependent
    //!
    void setParameterSectionDependent(const ParameterName id, bool depends);


    //! \brief Sets the maximum number of potentially sensitive parameters in the inlet function
    //!
    //!	The boundary condition at the column inlet is controlled by the user via a function that
    //! returns the concentration values for a given time. This inlet function may contain parameters
    //! for which the user wants to compute sensitivities. The exact number of potentially sensitive
    //! parameters in this function must be set before calling #initialize.
    //!
    //! \param	[in]	nInPar	The number of parameters in inlet function
    //!
    void setMaxSensInletParams(const int nInPar);


    //! \brief Reset the simulator to compute no sensitivity at all
    //!
    //!	All parameters are reset so as not to take part in any sensitivity computation.
    //!
    void resetSensParams();


    //! \brief Returns the names of all parameters set sensitive
    //!
    //! \return A vector containing the names of sensitive parameters as strings.
    //!
    std::vector<std::string> getSensModelParamNames() const;


    //! \brief Returns the components adressed by all parameters set sensitive
    //!
    //! \return A vector containing the components adressed by the sensitive parameters.
    //!
    std::vector<int> getSensModelParamComps() const;


    //! \brief Returns the sections adressed by all parameters set sensitive
    //!
    //! \return A vector containing the sections adressed by the sensitive parameters.
    //!
    std::vector<int> getSensModelParamSecs() const;


    //! \brief Returns the number of sensitive parameters
    //!
    //!	All parameters that are set sensitive are counted. This includes model parameters as well as inlet parameters.
    //!
    //! \return The number of sensitive parameters
    //!
    int getNSensParams() const;

    //! \brief Returns the number of sensitive model parameters
    //!
    //!	Only model parameters that are set sensitive are counted.
    //!
    //! \return The number of sensitive model parameters
    //!
    int getNSensModelParams() const;


    //! \brief Returns the number of sensitive inlet parameters
    //!
    //! Only inlet parameters that are set sensitive are counted.
    //!
    //! \return The number of sensitive inlet parameters
    //!
    int getNSensInletParams() const;


    //! \brief Sets the equidistant particle discretization scheme
    //!
    //!	The distribution of the radial nodes inside the beads will be equally spaced.
    //! This is the default behaviour when not calling any of the methods
    //! #setPartDiscSchemeEqVol or #setPartDiscSchemeUser.
    //!
    //! \image html equidistant_par.png "Distribution of radial nodes for equidistant discretization scheme"
    //!
    void setPartDiscSchemeEqDist();


    //! \brief Sets the equivolume particle discretization scheme
    //!
    //! The distribution of the radial nodes inside the beads will lead to a set of sphere shells
    //! which all feature the same volume.
    //!
    //! \image html equivolume_par.png "Distribution of radial nodes for equivolume discretization scheme"
    //!
    void setPartDiscSchemeEqVol();


    //! \brief Sets the user defined particle discretization scheme
    //!
    //! The user can specify a user defined particle discretization scheme through this method.
    //! A vector containing the bead-internal cell interfaces in dimensionless form must be provided.
    //! The coordinates may be in arbitrary order but should not exceed the interval \f$]0.0, 1.0[\f$.
    //!
    //! \image html userdefined_par.png "Distribution of radial nodes for user defined discretization scheme"
    //!
    //! \param  [in]    cellInterfaces  A vector containing all internal cell interface coordinates
    //!                                 in dimensionless form, i.e. in the interval \f$]0.0, 1.0[\f$
    //!                                 (length is <tt>npar - 1</tt>)
    //!
    void setPartDiscSchemeUser(const std::vector<double>& cellInterfaces);

    //! \brief Returns the dimensionless cell centers of the particle discretization cells
    //!
    //!	The method returns the cell-centered dimensionless coordinates of the particle discretization nodes
    //! in the range \f$[0.0, 1.0]\f$. Note: dimensionless values need to be multiplied by the bead radius,
    //! to obtain absolute values.
    //!
    //! \return A vector containing the dimensionless radial cell coordinates
    //!
    std::vector<double> getParCellCoords() const;

//    //! \brief Returns the dimensionless cell centers of the column discretization cells
//    //!
//    //! The method returns the cell-centered dimensionless coordinates of the column discretization nodes
//    //! in the range \f$[0.0, 1.0]\f$. Note: dimensionless values need to be multiplied by the column length,
//    //! to obtain absolute values.
//    //!
//    //! \return A vector containing the dimensionless axial cell coordinates
//    //!
//    std::vector<double> getColDiscCoords() const;

    //! \brief Registers a pointer to a user-defined class which implements the inlet profile
    //!
//     \deprecated
//    	This method registers a function pointer with the simulator. This funtion is called back
//     whenever the simulator needs to evaluate concentrations at the column inlet.
//     The inlet function takes the time as an input and returns the concentrations of
//     all \c ncomp chemical species.
//
//     Type: inletFunction
//     The function that is responsible for returning the concentration values at the column inlet
//     and their derivatives w.r.t. any inlet parameter must have type inletFunction.
//     It takes as input the value of the time \c t and the section \c sec at which
//     the concentrations are sought. \c inletParamIsSensitive provides information about which
//     of the inlet parameters was set sensitive. The vectors \c inletConc and
//     \c dInletConc_dp must be provided as refrence and should be filled within the function.
//     \c inletConc should contain the concentrations of all \c ncomp chemical species.
//     \c dInletConc_dp should contain the derivatives of the inlet concentration function w.r.t.
//     all parameters provided in \c inletParamIsSensitive.
//
    //! \param	[in]	inlet	    A pointer to a class implemnted by the user
    //!
    void setInletProfile(InletBase* inlet);

    //! \brief Sets the user supplied external function profile
    //!
    //! \param [in]     externalProfile     A pointer to a class implemented by the user

    void setExternalProfile(ExternalBase* externalProfile);

    //! \brief Sets the adsortion model to be a kinetic or quasi-stationary model
    //!
    //!	The model describing the adsorption behaviour at the inner bead surfaces can be either
    //! a dynamic model or a quasi-stationary model. The former includes time derivatives and
    //! can be activated by calling this method with a \c true boolean value. The latter leads to
    //! an algebraic equation that mimics instantaneous local equilibrium of concentrations in
    //! the mobile and bound phase. This can be activated by calling with \c false boolean value
    //! and is also set by default.
    //!
    //! \param	[in]	isKinetic	\c true to set a kinetic, \c false to set quasi-stationary adsorption model
    //!
    void setKineticAdsorptionModel(bool isKinetic);


    //! \brief Sets the time points at which a solution is written
    //!
    //!	The simulator returns solutions by default only at it's internal integration timepoints.
    //! This timepoints are usually unevenly spaced and are more dense in regions where
    //! gradients are steep. If the user is interested in solutions at specific timepoints
    //! he can set them using this method. The points may lie in the range \f$[t_{start}, t_{end}]\f$.
    //!
    //! \param	[in]	solutionsTimes	A vector containing the timepoints at which a solution shoudt be computed
    //!
    void setSolutionTimes(const std::vector<double>& solutionsTimes);


    //! \brief Sets the timepoints of (potential) discontinuities
    //!
    //!	Sets the timepoints of (potential) discontinuities where an integrator
    //! restart with adjusted initial conditions might be required.
    //! This applies to all timepoints at which the inlet concentration function is not
    //! continuously differentiable. Hereby the simulation is split into sections which are solved
    //! independently. The \p sectionTimes vector must contain at least \f$t_{start}\f$
    //! and \f$t_{end}\f$ and is of length <tt>nsec + 1</tt>.
    //!
    //! \param	[in]	sectionTimes	A vector containing the timepoints of all sections
    //!
    void setSectionTimes(const std::vector<double>& sectionTimes);

    //! \brief Sets the timepoints of (potential) discontinuities
    //!
    //! Sets the timepoints of (potential) discontinuities where an integrator
    //! restart with adjusted initial conditions might be required.
    //! This applies to all timepoints at which the inlet concentration function is not
    //! continuously differentiable. Hereby the simulation is split into sections which are solved
    //! independently. The \p sectionTimes vector must contain at least \f$t_{start}\f$
    //! and \f$t_{end}\f$ and is of length <tt>nsec + 1</tt>.
    //!
    //! Also sets the continuity of section transitions. A continuous section transition from section i
    //! to section i+1 is indicated by setting the ith entry of \p sectionContinuity to <tt>true</tt>.
    //! The integrator will not be reset on continuous section transitions, thus increasing performance.
    //!
    //! \param  [in]    sectionTimes        A vector containing the timepoints of all sections
    //! \param  [in]    sectionContinuity   A vector determining the continuity of section transitions
    //!
    void setSectionTimes(const std::vector<double>& sectionTimes, const std::vector<bool>& sectionContinuity);


    //! \brief Sets the initial conditions, initializes the time integrator and the sensitivity setup
    //!
    //!	The method sets the initial conditions for all components in the mobile phase as well as in
    //! the bound phase using the concentrations supplied in \p initC and \p initQ. All initialization
    //! for the time integration by IDAS, e.g. specifying tolerances, setting up the linear solver, etc.
    //! is then done. Furthermore sensitivity related setup routines are processed, e.g. activating the
    //! sensitivities, setting seed vectors for AD computations, specifying integration tolerances, etc.
    //!
    //! \param	[in]	initC	A vector containing \c ncomp initial concentration values for the mobile phase
    //! \param	[in]	initQ	A vector containing \c ncomp initial concentration values for the bound phase
    //!
    void initialize(const std::vector<double>& initC, const std::vector<double>& initQ);


    //! \brief Starts the solution of the system specified for this simulator object
    //!
    //!	Checks all model parameters to lie inside their possible bounds and then runs the time integration
    //! from \f$t_{start}\f$ to \f$t_{end}\f$, restarting the time integrator at every section time specified
    //! by #setSectionTimes. Solutions are stored by default at internal timesteps of the integration routine,
    //! or if a call to #setSolutionTimes was made, at the specified user defined solution timepoints.
    //!
    void integrate();


    //! \brief Returns the timepoints at which a solution was stored
    //!
    //!	When there was no call to #setSolutionTimes made, the method returns the internal timepoints at which
    //! the integration routine has computed a solution. Otherwise the returned timepoints are identical to
    //! the ones set through #setSolutionTimes.
    //!
    //! \param	[out]	userVector	A reference to a vector which will be filled with the solution times
    //!
    void getSolutionTimes(std::vector<double>& userVector) const;


    //! \brief Returns the bare state vectors concatenated for all timepoints
    //!
    //!	The method returns the solution as it was written to the memory during solution procedure.
    //! This means, for each timepoint the whole state vetor (including column, particle and boundary parts)
    //! is concatenated and assigned to the supplied vector. The ordering is rather unhandy, since the
    //! column, the particle and the boundary parts are of different size.
    //!
    //! \param	[out]   userVector  A reference to a vector which will be filled with the whole solution
    //!
    void getAllSolutions(std::vector<double>& userVector) const;


    //! \brief Returns the solution in the column at all timepoints
    //!
    //!	The method returns only the column parts of the state vector. The vector can be
    //! viewed as a three-dimensional tensor defined by the dimensions \f$ N_t \cdot N_c \cdot N_z \f$
    //! with
    //! - \f$ N_t \f$ = Number of timepoints
    //! - \f$ N_c \f$ = Number of components (chemical species) (\c ncomp)
    //! - \f$ N_z \f$ = Number of axial discretization nodes (\c ncol)
    //!
    //! Here's how the indexing works:
    //!
    //!     for t in 0 .. ntpts - 1  // Outermost index
    //!         for c in 0 .. ncomp - 1
    //!             for z in 0 .. ncol - 1  // Innermost index
    //!                 // work with userVector.at(t * ncomp * ncol + c * ncol + z);
    //!             end
    //!         end
    //!     end
    //!
    //! \image html column_solution_vector_structure.png "Structure of the column solution vector"
    //!
    //! \param	[out]   userVector  A reference to a vector which will be filled with the column solution
    //!
    void getColSolutions(std::vector<double>& userVector) const;


    //! \brief Returns the solution in the particles at all timepoints
    //!
    //! The method returns only the particle parts of the state vector. The vector can be
    //! viewed as a five-dimensional tensor defined by the dimensions \f$ N_t \cdot N_z \cdot N_r \cdot N_p \cdot N_c \f$
    //! with
    //! - \f$ N_t \f$ = Number of timepoints
    //! - \f$ N_z \f$ = Number of axial discretization nodes (\c ncol)
    //! - \f$ N_r \f$ = Number of radial discretization nodes (\c npar)
    //! - \f$ N_p \f$ = Number of phases a species can be assigned to (\c nphase)
    //! - \f$ N_c \f$ = Number of components (chemical species) (\c ncomp)
    //!
    //! Here's how the indexing works:
    //!
    //!     for t in 0 .. ntpts - 1  // Outermost index
    //!         for z in 0 .. ncol - 1
    //!             for r in 0 .. npar - 1
    //!                 for p in 0 .. nphase - 1
    //!                     for c in 0 .. ncomp - 1  // Innermost index
    //!                         // work with userVector.at(t * ncol * npar * nphase * ncomp +
    //!                                                    z * npar * nphase * ncomp +
    //!                                                    r * nphase * ncomp +
    //!                                                    p * ncomp +
    //!                                                    c );
    //!                     end
    //!                 end
    //!             end
    //!         end
    //!     end
    //!
    //! \image html particle_solution_vector_structure.png "Structure of the particle solution vector"
    //!
    //! \param  [out]   userVector  A reference to a vector which will be filled with the particle solution
    //!
    void getParSolutions(std::vector<double>& userVector) const;
    void getBndSolutions(std::vector<double>& userVector) const;

    void getAllSensitivities(std::vector<double>& userVector) const;
    void getColSensitivities(std::vector<double>& userVector) const;
    void getParSensitivities(std::vector<double>& userVector) const;
    void getBndSensitivities(std::vector<double>& userVector) const;

    void getSolutionColumnOutlet(std::vector<double>& userVector, const int comp) const;
    void getSolutionColumnInlet(std::vector<double>& userVector, const int comp) const;
    void getSensitivityColumnOutlet(std::vector<double>& userVector, const int comp, const std::tuple<ParameterName, int, int>& param) const;
    void getSensitivityColumnOutlet(std::vector<double>& userVector, const int comp, const ParameterName id, const int paramComp = -1, const int paramSec = -1) const;


    // use analytical jacobian implementation for computations instead of automatical derivative information
    void useAnalyticJacobian(const bool analyticJac = false);

    // configure routines with default values for convenient setup of simulator parts
    void configureSchurSolver(double schurSafety = 1e-8, int maxRestarts = 0, int maxKrylov = 0, int gramSchmidtType = 1);
    void configureWenoScheme(int wenoOrder = 3, double wenoEps = 1e-6, int boundaryModel = 0);
    void configureTimeIntegrator(double relTol = 1e-12, double absTol = 1e-9,
            double initStepSize = 1e-6, int maxSteps = 10000);
    void configurePrinting(bool printProgress = false, bool printStatistics = false,
            bool printTiming = false, bool printParamList = false, bool printConfig = false);

private:
    SimulatorPImpl* _sim;
};

} // namespace cadet

#endif  // LIBCADET_SIMULATOR_HPP_
