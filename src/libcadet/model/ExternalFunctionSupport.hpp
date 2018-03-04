// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2018: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

/**
 * @file 
 * Provides support functions for defining externally dependent models.
 */

#ifndef LIBCADET_EXTFUNSUPPORT_HPP_
#define LIBCADET_EXTFUNSUPPORT_HPP_

#include "cadet/ExternalFunction.hpp"
#include "cadet/Exceptions.hpp"
 
#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <order/interpreter.h>

#include <vector>
#include <algorithm>
#include <string>


#define CADET_DEFINE_EXTDEP_VARIABLE(TYPE, VAR) \
	mutable TYPE VAR;                           \
	TYPE VAR##T0;                               \
	TYPE VAR##T1;                               \
	TYPE VAR##T2;                               \
	TYPE VAR##T3;

#define CADET_UPDATE_EXTDEP_VARIABLE_BRACES(VAR, IDXEXPR, EXTVAL) \
	VAR[IDXEXPR] = VAR##T0[IDXEXPR] + EXTVAL * (VAR##T1[IDXEXPR] + EXTVAL * (VAR##T2[IDXEXPR] + EXTVAL * VAR##T3[IDXEXPR]));

#define CADET_UPDATE_EXTDEP_VARIABLE_NATIVE(VAR, IDXEXPR, EXTVAL) \
	VAR.native(IDXEXPR) = VAR##T0.native(IDXEXPR) + EXTVAL * (VAR##T1.native(IDXEXPR) + EXTVAL * (VAR##T2.native(IDXEXPR) + EXTVAL * VAR##T3.native(IDXEXPR)));

#define CADET_UPDATE_EXTDEP_VARIABLE(VAR, EXTVAL) \
	VAR = VAR##T0 + EXTVAL * (VAR##T1 + EXTVAL * (VAR##T2 + EXTVAL * VAR##T3));


#define CADET_UPDATE_EXTDEP_VARIABLE_TDIFF_BRACES(VAR, IDXEXPR, EXTVAL, TDIFFVAL) \
	VAR[IDXEXPR] = TDIFFVAL * (static_cast<double>(VAR##T1[IDXEXPR]) + EXTVAL * (2.0 *  static_cast<double>(VAR##T2[IDXEXPR]) + 3.0 * EXTVAL * static_cast<double>(VAR##T3[IDXEXPR])));

#define CADET_UPDATE_EXTDEP_VARIABLE_TDIFF_NATIVE(VAR, IDXEXPR, EXTVAL, TDIFFVAL) \
	VAR.native(IDXEXPR) = TDIFFVAL * (static_cast<double>(VAR##T1.native(IDXEXPR)) + EXTVAL * (2.0 *  static_cast<double>(VAR##T2.native(IDXEXPR)) + 3.0 * EXTVAL * static_cast<double>(VAR##T3.native(IDXEXPR))));

#define CADET_UPDATE_EXTDEP_VARIABLE_TDIFF(VAR, EXTVAL, TDIFFVAL) \
	VAR = TDIFFVAL * (static_cast<double>(VAR##T1) + EXTVAL * (2.0 * static_cast<double>(VAR##T2) + 3.0 * EXTVAL * static_cast<double>(VAR##T3)));


#define CADET_RESERVE_SPACE(VAR, NUMELEM) \
	VAR##T0.reserve(NUMELEM);             \
	VAR##T1.reserve(NUMELEM);             \
	VAR##T2.reserve(NUMELEM);             \
	VAR##T3.reserve(NUMELEM);

#define CADET_RESERVE_SPACE2(VAR, NUMELEM, NUMSLICES) \
	VAR##T0.reserve(NUMELEM, NUMSLICES);              \
	VAR##T1.reserve(NUMELEM, NUMSLICES);              \
	VAR##T2.reserve(NUMELEM, NUMSLICES);              \
	VAR##T3.reserve(NUMELEM, NUMSLICES);


#define CADET_REGPAR_SCALAR(NAME, PARAMETERS, VAR, UNITOPIDX)                                                                                 \
	PARAMETERS[makeParamId(hashString("EXT_" NAME), UNITOPIDX, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &VAR##T0;          \
	PARAMETERS[makeParamId(hashString("EXT_" NAME "_T"), UNITOPIDX, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &VAR##T1;     \
	PARAMETERS[makeParamId(hashString("EXT_" NAME "_TT"), UNITOPIDX, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &VAR##T2;    \
	PARAMETERS[makeParamId(hashString("EXT_" NAME "_TTT"), UNITOPIDX, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &VAR##T3;

#define CADET_REGPAR_COMPSEC(NAME, PARAMETERS, VAR, UNITOPIDX, NCOMP)                                              \
	VAR = VAR##T0;                                                                                                 \
	registerComponentSectionDependentParam(hashString("EXT_" NAME), PARAMETERS, VAR##T0, UNITOPIDX, NCOMP);        \
	registerComponentSectionDependentParam(hashString("EXT_" NAME "_T"), PARAMETERS, VAR##T1, UNITOPIDX, NCOMP);   \
	registerComponentSectionDependentParam(hashString("EXT_" NAME "_TT"), PARAMETERS, VAR##T2, UNITOPIDX, NCOMP);  \
	registerComponentSectionDependentParam(hashString("EXT_" NAME "_TTT"), PARAMETERS, VAR##T3, UNITOPIDX, NCOMP);

#define CADET_REGPAR_COMPBND_VEC(NAME, PARAMETERS, VAR, UNITOPIDX)                                             \
	VAR = VAR##T0;                                                                                             \
	registerComponentBoundStateDependentParam(hashString("EXT_" NAME), PARAMETERS, VAR##T0, UNITOPIDX);        \
	registerComponentBoundStateDependentParam(hashString("EXT_" NAME "_T"), PARAMETERS, VAR##T1, UNITOPIDX);   \
	registerComponentBoundStateDependentParam(hashString("EXT_" NAME "_TT"), PARAMETERS, VAR##T2, UNITOPIDX);  \
	registerComponentBoundStateDependentParam(hashString("EXT_" NAME "_TTT"), PARAMETERS, VAR##T3, UNITOPIDX);

#define CADET_REGPAR_COMPBND(NAME, PARAMETERS, VAR, UNITOPIDX)                                                 \
	VAR = VAR##T0;                                                                                             \
	registerComponentBoundStateDependentParam(hashString("EXT_" NAME), PARAMETERS, VAR##T0, UNITOPIDX);        \
	registerComponentBoundStateDependentParam(hashString("EXT_" NAME "_T"), PARAMETERS, VAR##T1, UNITOPIDX);   \
	registerComponentBoundStateDependentParam(hashString("EXT_" NAME "_TT"), PARAMETERS, VAR##T2, UNITOPIDX);  \
	registerComponentBoundStateDependentParam(hashString("EXT_" NAME "_TTT"), PARAMETERS, VAR##T3, UNITOPIDX);

#define CADET_REGPAR_COMPBND_COMPMAJOR(NAME, PARAMETERS, VAR, UNITOPIDX)                                                \
	VAR = VAR##T0;                                                                                                      \
	registerComponentBoundStateDependentParamCompMajor(hashString("EXT_" NAME), PARAMETERS, VAR##T0, UNITOPIDX);        \
	registerComponentBoundStateDependentParamCompMajor(hashString("EXT_" NAME "_T"), PARAMETERS, VAR##T1, UNITOPIDX);   \
	registerComponentBoundStateDependentParamCompMajor(hashString("EXT_" NAME "_TT"), PARAMETERS, VAR##T2, UNITOPIDX);  \
	registerComponentBoundStateDependentParamCompMajor(hashString("EXT_" NAME "_TTT"), PARAMETERS, VAR##T3, UNITOPIDX);

#define CADET_REGPAR_SCALARBND(NAME, PARAMETERS, VAR, UNITOPIDX)                                            \
	VAR = VAR##T0;                                                                                          \
	registerScalarBoundStateDependentParam(hashString("EXT_" NAME), PARAMETERS, VAR##T0, UNITOPIDX);        \
	registerScalarBoundStateDependentParam(hashString("EXT_" NAME "_T"), PARAMETERS, VAR##T1, UNITOPIDX);   \
	registerScalarBoundStateDependentParam(hashString("EXT_" NAME "_TT"), PARAMETERS, VAR##T2, UNITOPIDX);  \
	registerScalarBoundStateDependentParam(hashString("EXT_" NAME "_TTT"), PARAMETERS, VAR##T3, UNITOPIDX);


#define CADET_READPAR_MATRIX(VAR, PARAMPROVIDER, NAME, NCOMP, NUM)               \
	readParameterMatrix(VAR##T0, PARAMPROVIDER, "EXT_" NAME, NCOMP, NUM);        \
	readParameterMatrix(VAR##T1, PARAMPROVIDER, "EXT_" NAME "_T", NCOMP, NUM);   \
	readParameterMatrix(VAR##T2, PARAMPROVIDER, "EXT_" NAME "_TT", NCOMP, NUM);  \
	readParameterMatrix(VAR##T3, PARAMPROVIDER, "EXT_" NAME "_TTT", NCOMP, NUM);

#define CADET_READPAR_BOUNDSTATEDEP(TYPE1, TYPE2, VAR, PARAMPROVIDER, NAME, NCOMP, NUM)                     \
	readBoundStateDependentParameter<TYPE1, TYPE2>(VAR##T0, PARAMPROVIDER, "EXT_" NAME, NCOMP, NUM);        \
	readBoundStateDependentParameter<TYPE1, TYPE2>(VAR##T1, PARAMPROVIDER, "EXT_" NAME "_T", NCOMP, NUM);   \
	readBoundStateDependentParameter<TYPE1, TYPE2>(VAR##T2, PARAMPROVIDER, "EXT_" NAME "_TT", NCOMP, NUM);  \
	readBoundStateDependentParameter<TYPE1, TYPE2>(VAR##T3, PARAMPROVIDER, "EXT_" NAME "_TTT", NCOMP, NUM);

#define CADET_READPAR_BOUNDSTATEDEP_MATRIX(TYPE1, TYPE2, VAR, PARAMPROVIDER, NAME, NCOMP, NUM)                          \
	readMatrixValuedBoundStateDependentParameter<TYPE1, TYPE2>(VAR##T0, PARAMPROVIDER, "EXT_" NAME, NCOMP, NUM);        \
	readMatrixValuedBoundStateDependentParameter<TYPE1, TYPE2>(VAR##T1, PARAMPROVIDER, "EXT_" NAME "_T", NCOMP, NUM);   \
	readMatrixValuedBoundStateDependentParameter<TYPE1, TYPE2>(VAR##T2, PARAMPROVIDER, "EXT_" NAME "_TT", NCOMP, NUM);  \
	readMatrixValuedBoundStateDependentParameter<TYPE1, TYPE2>(VAR##T3, PARAMPROVIDER, "EXT_" NAME "_TTT", NCOMP, NUM);

#define CADET_READPAR_SCALAR(VAR, PARAMPROVIDER, NAME)     \
	VAR##T0 = PARAMPROVIDER.getDouble("EXT_" NAME);        \
	VAR##T1 = PARAMPROVIDER.getDouble("EXT_" NAME "_T");   \
	VAR##T2 = PARAMPROVIDER.getDouble("EXT_" NAME "_TT");  \
	VAR##T3 = PARAMPROVIDER.getDouble("EXT_" NAME "_TTT");

#define CADET_READPAR_SCALARBND(VAR, PARAMPROVIDER, NAME, NUM)                  \
	readScalarParameterOrArray(VAR##T0, PARAMPROVIDER, "EXT_" NAME, NUM);       \
	readScalarParameterOrArray(VAR##T1, PARAMPROVIDER, "EXT_" NAME "_T", NUM);  \
	readScalarParameterOrArray(VAR##T2, PARAMPROVIDER, "EXT_" NAME "_TT", NUM); \
	readScalarParameterOrArray(VAR##T3, PARAMPROVIDER, "EXT_" NAME "_TTT", NUM);


/**
 * @brief Generates constant and externally dependent parameter handler classes
 * @param localName Name of the constant parameter handler
 * @param extName Name of the externally dependent parameter handler
 * @param params List of tuples with externally dependent parameters
 * @param constparams List of tuples with constant parameters
 */
#define CADET_BINDINGPARAMS(localName, extName, params, constparams) \
	class localName : public cadet::model::ConstParamHandlerBase \
	{ \
	public: \
		struct ConstParams \
		{ \
			ORDER_PP(8seq_for_each( \
			        	8fn(8X, \
			        		8print( (typename) 8space 8tuple_at_0(8X) (::storage_t) 8space 8tuple_at_1(8X) (;)) \
			        	), \
			        	8vseq_to_seq_of_tuples(8(params)) \
			        ) \
			) \
			ORDER_PP(8seq_for_each( \
			        	8fn(8X, \
			        		8if(8is_tuple(8tuple_at_1(8X)), \
			        			8seq_for_each( \
			        				8fn(8Y, \
			        					8print( (typename) 8space 8tuple_at_0(8X) (::storage_t) 8space 8Y (;) ) \
			        				), \
			        				8tuple_to_seq(8tuple_at_1(8X)) \
			        			), \
			        			8print( (typename) 8space 8tuple_at_0(8X) (::storage_t) 8space 8tuple_at_1(8X) (;)), \
			        			) \
			        	), \
			        	8vseq_to_seq_of_tuples(8(constparams)) \
			        ) \
			) \
		}; \
		typedef ConstParams params_t; \
		static inline const char* identifier() CADET_NOEXCEPT; \
		localName() CADET_NOEXCEPT : \
			ORDER_PP(8seq_for_each_with_delimiter( \
			        	8fn(8X, \
			        		8print( 8cat(8quote(_), 8tuple_at_1(8X)) 8lparen (&_localParams.) 8tuple_at_1(8X) 8rparen ) \
			        	), \
			        	8emit_comma, \
			        	8vseq_to_seq_of_tuples(8(params)) \
			        ) \
			) \
			ORDER_PP(8if( 8isnt_nil( 8(constparams) ), 8print(8comma), 8nil)) \
			ORDER_PP(8seq_for_each( \
			        	8fn(8X, \
			        		8if(8is_tuple(8tuple_at_1(8X)), \
			        			8chain( 8print( 8tuple_at_2(8X) 8lparen ), \
			        			8seq_for_each_with_delimiter( \
			        				8fn(8Y, \
			        					8print( (&_localParams.) 8Y ) \
			        				), \
			        				8emit_comma, \
			        				8tuple_to_seq(8tuple_at_1(8X)) \
			        			), 8print (8rparen) ), \
			        			8print( 8cat(8quote(_), 8tuple_at_1(8X)) 8lparen (&_localParams.) 8tuple_at_1(8X) 8rparen ), \
			        			) \
			        	), \
			        	8vseq_to_seq_of_tuples(8(constparams)) \
			        ) \
			) \
		{ } \
		inline bool configure(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBoundStates) \
		{ \
			ORDER_PP(8seq_for_each( \
			        	8fn(8X, \
			        		8print(8cat(8quote(_), 8tuple_at_1(8X)) (.configure) 8lparen 8tuple_at_2(8X) 8comma (paramProvider) 8comma (nComp) 8comma (nBoundStates) 8rparen (;) ) \
			        	), \
			        	8vseq_to_seq_of_tuples(8(params)) \
			        ) \
			) \
			ORDER_PP(8seq_for_each( \
			        	8fn(8X, \
			        		8if(8is_tuple(8tuple_at_1(8X)), \
			        			8print(8tuple_at_2(8X) (.configure) 8lparen 8tuple_at_3(8X) 8comma (paramProvider) 8comma (nComp) 8comma (nBoundStates) 8rparen (;) ), \
			        			8print(8cat(8quote(_), 8tuple_at_1(8X)) (.configure) 8lparen 8tuple_at_2(8X) 8comma (paramProvider) 8comma (nComp) 8comma (nBoundStates) 8rparen (;) ) \
			        			) \
			        	), \
			        	8vseq_to_seq_of_tuples(8(constparams)) \
			        ) \
			) \
			return validateConfig(nComp, nBoundStates); \
		} \
		inline void registerParameters(std::unordered_map<ParameterId, active*>& parameters, unsigned int unitOpIdx, unsigned int nComp, unsigned int const* nBoundStates) \
		{ \
			ORDER_PP(8seq_for_each( \
			        	8fn(8X, \
			        		8print(8cat(8quote(_), 8tuple_at_1(8X)) (.registerParam) 8lparen 8tuple_at_2(8X) 8comma (parameters) 8comma (unitOpIdx) 8comma (nComp) 8comma (nBoundStates) 8rparen (;) ) \
			        	), \
			        	8vseq_to_seq_of_tuples(8(params)) \
			        ) \
			) \
			ORDER_PP(8seq_for_each( \
			        	8fn(8X, \
			        		8if(8is_tuple(8tuple_at_1(8X)), \
			        			8print(8tuple_at_2(8X) (.registerParam) 8lparen 8tuple_at_3(8X) 8comma (parameters) 8comma (unitOpIdx) 8comma (nComp) 8comma (nBoundStates) 8rparen (;) ), \
			        			8print(8cat(8quote(_), 8tuple_at_1(8X)) (.registerParam) 8lparen 8tuple_at_2(8X) 8comma (parameters) 8comma (unitOpIdx) 8comma (nComp) 8comma (nBoundStates) 8rparen (;) ) \
			        			) \
			        	), \
			        	8vseq_to_seq_of_tuples(8(constparams)) \
			        ) \
			) \
		} \
		inline void reserve(unsigned int numElem, unsigned int numSlices, unsigned int nComp, unsigned int const* nBoundStates) \
		{ \
			ORDER_PP(8seq_for_each( \
			        	8fn(8X, \
			        		8print(8cat(8quote(_), 8tuple_at_1(8X)) (.reserve) 8lparen (numElem) 8comma (numSlices) 8comma (nComp) 8comma (nBoundStates) 8rparen (;) ) \
			        	), \
			        	8vseq_to_seq_of_tuples(8(params)) \
			        ) \
			) \
			ORDER_PP(8seq_for_each( \
			        	8fn(8X, \
			        		8if(8is_tuple(8tuple_at_1(8X)), \
			        			8print(8tuple_at_2(8X) (.reserve) 8lparen (numElem) 8comma (numSlices) 8comma (nComp) 8comma (nBoundStates) 8rparen (;) ), \
			        			8print(8cat(8quote(_), 8tuple_at_1(8X)) (.reserve) 8lparen (numElem) 8comma (numSlices) 8comma (nComp) 8comma (nBoundStates) 8rparen (;) ) \
			        			) \
			        	), \
			        	8vseq_to_seq_of_tuples(8(constparams)) \
			        ) \
			) \
		} \
		inline const params_t& update(double t, double z, double r, unsigned int secIdx, unsigned int nComp, unsigned int const* nBoundStates, void* workSpace) const \
		{ \
			return _localParams; \
		} \
		inline const params_t updateTimeDerivative(double t, double z, double r, unsigned int secIdx, unsigned int nComp, unsigned int const* nBoundStates, void* workSpace) const \
		{ \
			return _localParams; \
		} \
	protected: \
		inline bool validateConfig(unsigned int nComp, unsigned int const* nBoundStates); \
		ConstParams _localParams; \
		ORDER_PP(8seq_for_each( \
		        	8fn(8X, \
		        		8print( 8tuple_at_0(8X) 8space 8cat(8quote(_), 8tuple_at_1(8X)) (;)) \
		        	), \
		        	8vseq_to_seq_of_tuples(8(params)) \
		        ) \
		) \
		ORDER_PP(8seq_for_each( \
		        	8fn(8X, \
		        		8print( 8tuple_at_0(8X) 8space 8tuple_at_2(8X) (;)) \
		        	), \
		        	8vseq_to_seq_of_tuples(8(constparams)) \
		        ) \
		) \
	}; \
	class extName : public cadet::model::ExternalParamHandlerBase \
	{ \
	public: \
		struct ConstParams \
		{ \
			ORDER_PP(8seq_for_each( \
			        	8fn(8X, \
			        		8if(8is_tuple(8tuple_at_1(8X)), \
			        			8seq_for_each( \
			        				8fn(8Y, \
			        					8print( (typename) 8space 8tuple_at_0(8X) (::storage_t) 8space 8Y (;) ) \
			        				), \
			        				8tuple_to_seq(8tuple_at_1(8X)) \
			        			), \
			        			8print( (typename) 8space 8tuple_at_0(8X) (::storage_t) 8space 8tuple_at_1(8X) (;)), \
			        			) \
			        	), \
			        	8vseq_to_seq_of_tuples(8(constparams)) \
			        ) \
			) \
		}; \
		struct VariableParams \
		{ \
			ORDER_PP(8seq_for_each( \
			        	8fn(8X, \
			        		8print( (typename) (util::localVersionOf) (<) (typename) 8space 8tuple_at_0(8X) (::storage_t) (>) (::type) 8space 8tuple_at_1(8X) (;)) \
			        	), \
			        	8vseq_to_seq_of_tuples(8(params)) \
			        ) \
			) \
			ORDER_PP(8seq_for_each( \
			        	8fn(8X, \
			        		8if(8is_tuple(8tuple_at_1(8X)), \
			        			8seq_for_each( \
			        				8fn(8Y, \
			        					8print( (typename) 8tuple_at_0(8X) (::storage_t) 8space 8Y (;) ) \
			        				), \
			        				8tuple_to_seq(8tuple_at_1(8X)) \
			        			), \
			        			8print( (typename) 8tuple_at_0(8X) (::storage_t) 8space 8tuple_at_1(8X) (;)), \
			        			) \
			        	), \
			        	8vseq_to_seq_of_tuples(8(constparams)) \
			        ) \
			) \
		}; \
		typedef VariableParams params_t; \
		static inline const char* identifier() CADET_NOEXCEPT; \
		extName() CADET_NOEXCEPT ORDER_PP(8if( 8isnt_nil( 8(constparams) ), 8print((:)), 8nil)) \
			ORDER_PP(8seq_for_each( \
			        	8fn(8X, \
			        		8if(8is_tuple(8tuple_at_1(8X)), \
			        			8chain( 8print( 8tuple_at_2(8X) 8lparen ), \
			        			8seq_for_each_with_delimiter( \
			        				8fn(8Y, \
			        					8print( (&_constParams.) 8Y ) \
			        				), \
			        				8emit_comma, \
			        				8tuple_to_seq(8tuple_at_1(8X)) \
			        			), 8print (8rparen) ), \
			        			8print( 8cat(8quote(_), 8tuple_at_1(8X)) 8lparen (&_constParams.) 8tuple_at_1(8X) 8rparen ), \
			        			) \
			        	), \
			        	8vseq_to_seq_of_tuples(8(constparams)) \
			        ) \
			) \
		{ } \
		inline bool configure(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBoundStates) \
		{ \
			ORDER_PP(8seq_for_each( \
			        	8fn(8X, \
			        		8print(8cat(8quote(_), 8tuple_at_1(8X)) (.configure) 8lparen 8tuple_at_2(8X) 8comma (paramProvider) 8comma (nComp) 8comma (nBoundStates) 8rparen (;)) \
			        	), \
			        	8vseq_to_seq_of_tuples(8(params)) \
			        ) \
			) \
			ORDER_PP(8seq_for_each( \
			        	8fn(8X, \
			        		8if(8is_tuple(8tuple_at_1(8X)), \
			        			8print(8tuple_at_2(8X) (.configure) 8lparen 8tuple_at_3(8X) 8comma (paramProvider) 8comma (nComp) 8comma (nBoundStates) 8rparen (;) ), \
			        			8print(8cat(8quote(_), 8tuple_at_1(8X)) (.configure) 8lparen 8tuple_at_2(8X) 8comma (paramProvider) 8comma (nComp) 8comma (nBoundStates) 8rparen (;) ) \
			        			) \
			        	), \
			        	8vseq_to_seq_of_tuples(8(constparams)) \
			        ) \
			) \
			ExternalParamHandlerBase::configure(paramProvider, ORDER_PP(8to_lit(8seq_size(8(params))))); \
			return validateConfig(nComp, nBoundStates); \
		} \
		inline void registerParameters(std::unordered_map<ParameterId, active*>& parameters, unsigned int unitOpIdx, unsigned int nComp, unsigned int const* nBoundStates) \
		{ \
			ORDER_PP(8seq_for_each( \
			        	8fn(8X, \
			        		8print(8cat(8quote(_), 8tuple_at_1(8X)) (.registerParam) 8lparen 8tuple_at_2(8X) 8comma (parameters) 8comma (unitOpIdx) 8comma (nComp) 8comma (nBoundStates) 8rparen (;)) \
			        	), \
			        	8vseq_to_seq_of_tuples(8(params)) \
			        ) \
			) \
			ORDER_PP(8seq_for_each( \
			        	8fn(8X, \
			        		8if(8is_tuple(8tuple_at_1(8X)), \
			        			8print(8tuple_at_2(8X) (.registerParam) 8lparen 8tuple_at_3(8X) 8comma (parameters) 8comma (unitOpIdx) 8comma (nComp) 8comma (nBoundStates) 8rparen (;) ), \
			        			8print(8cat(8quote(_), 8tuple_at_1(8X)) (.registerParam) 8lparen 8tuple_at_2(8X) 8comma (parameters) 8comma (unitOpIdx) 8comma (nComp) 8comma (nBoundStates) 8rparen (;) ) \
			        			) \
			        	), \
			        	8vseq_to_seq_of_tuples(8(constparams)) \
			        ) \
			) \
		} \
		inline void reserve(unsigned int numElem, unsigned int numSlices, unsigned int nComp, unsigned int const* nBoundStates) \
		{ \
			ORDER_PP(8seq_for_each( \
			        	8fn(8X, \
			        		8print(8cat(8quote(_), 8tuple_at_1(8X)) (.reserve) 8lparen (numElem) 8comma (numSlices) 8comma (nComp) 8comma (nBoundStates) 8rparen (;)) \
			        	), \
			        	8vseq_to_seq_of_tuples(8(params)) \
			        ) \
			) \
			ORDER_PP(8seq_for_each( \
			        	8fn(8X, \
			        		8if(8is_tuple(8tuple_at_1(8X)), \
			        			8print(8tuple_at_2(8X) (.reserve) 8lparen (numElem) 8comma (numSlices) 8comma (nComp) 8comma (nBoundStates) 8rparen (;) ), \
			        			8print(8cat(8quote(_), 8tuple_at_1(8X)) (.reserve) 8lparen (numElem) 8comma (numSlices) 8comma (nComp) 8comma (nBoundStates) 8rparen (;) ) \
			        			) \
			        	), \
			        	8vseq_to_seq_of_tuples(8(constparams)) \
			        ) \
			) \
		} \
		inline const params_t& update(double t, double z, double r, unsigned int secIdx, unsigned int nComp, unsigned int const* nBoundStates, void* workSpace) const \
		{ \
			params_t* const localParams = reinterpret_cast<params_t*>(workSpace); \
			new (localParams) params_t; \
			double* const extFunBuffer = cadet::util::advancePointer<double>(workSpace, sizeof(params_t)); \
			void* buffer = cadet::util::advancePointer<void>(workSpace, sizeof(params_t) + 2 * ORDER_PP(8to_lit(8seq_size(8(params)))) * sizeof(double)); \
			evaluateExternalFunctions(t, z, r, secIdx, ORDER_PP(8to_lit(8seq_size(8(params)))), extFunBuffer); \
			ORDER_PP(8seq_for_each( \
			        	8fn(8X, \
			        		8if(8is_tuple(8tuple_at_1(8X)), \
			        			8seq_for_each( \
			        				8fn(8Y, \
			        					8print( (localParams->) 8Y (=) (_constParams.) 8Y (;) ) \
			        				), \
			        				8tuple_to_seq(8tuple_at_1(8X)) \
			        			), \
			        			8print( (localParams->) 8tuple_at_1(8X) (=) (_constParams.) 8tuple_at_1(8X) (;) ) \
			        			) \
			        	), \
			        	8vseq_to_seq_of_tuples(8(constparams)) \
			        ) \
			) \
			ORDER_PP(8seq_for_each_with_idx( \
			        	8fn(8I, 8X, \
			        		8print( \
			        			8cat(8quote(_), 8tuple_at_1(8X)) (.prepareCache) 8lparen (localParams->)8tuple_at_1(8X) 8comma (buffer) 8rparen (;) \
			        			8cat(8quote(_), 8tuple_at_1(8X)) (.update) 8lparen (cadet::util::dataOfLocalVersion) 8lparen (localParams->)8tuple_at_1(8X) 8rparen 8comma (extFunBuffer[) 8I (]) 8comma (nComp) 8comma (nBoundStates) 8rparen (;) \
			        			(buffer=) (cadet::util::advancePointer) 8lparen (buffer) 8comma (cadet::util::memoryForDataOf) 8lparen (localParams->)8tuple_at_1(8X) 8rparen 8rparen (;) \
			        		) \
			        	), \
			        	0, \
			        	8vseq_to_seq_of_tuples(8(params)) \
			        ) \
			) \
			return *localParams; \
		} \
		inline params_t updateTimeDerivative(double t, double z, double r, unsigned int secIdx, unsigned int nComp, unsigned int const* nBoundStates, void* workSpace) const \
		{ \
			VariableParams p; \
			params_t* const localParams = reinterpret_cast<params_t*>(workSpace); \
			double* const extFunBuffer = cadet::util::advancePointer<double>(workSpace, sizeof(params_t)); \
			double* const extDerivBuffer = extFunBuffer + ORDER_PP(8to_lit(8seq_size(8(params)))); \
			void* buffer = util::ptrToEndOfData( localParams->ORDER_PP(8tuple_at_1(8seq_last( 8vseq_to_seq_of_tuples(8(params)) ))) ); \
			evaluateTimeDerivativeExternalFunctions(t, z, r, secIdx, ORDER_PP(8to_lit(8seq_size(8(params)))), extDerivBuffer); \
			ORDER_PP(8seq_for_each( \
			        	8fn(8X, \
			        		8if(8is_tuple(8tuple_at_1(8X)), \
			        			8seq_for_each( \
			        				8fn(8Y, \
			        					8print( (p.) 8Y (=) (_constParams.) 8Y (;) ) \
			        				), \
			        				8tuple_to_seq(8tuple_at_1(8X)) \
			        			), \
			        			8print( (p.) 8tuple_at_1(8X) (=) (_constParams.) 8tuple_at_1(8X) (;) ) \
			        			) \
			        	), \
			        	8vseq_to_seq_of_tuples(8(constparams)) \
			        ) \
			) \
			ORDER_PP(8seq_for_each_with_idx( \
			        	8fn(8I, 8X, \
			        		8print( \
			        			8cat(8quote(_), 8tuple_at_1(8X)) (.prepareCache) 8lparen (p.)8tuple_at_1(8X) 8comma (buffer) 8rparen (;) \
			        			8cat(8quote(_), 8tuple_at_1(8X)) (.updateTimeDerivative) 8lparen (cadet::util::dataOfLocalVersion) 8lparen (p.)8tuple_at_1(8X) 8rparen 8comma (extFunBuffer[) 8I (]) 8comma (extDerivBuffer[) 8I (]) 8comma (nComp) 8comma (nBoundStates) 8rparen (;) \
			        			(buffer=) (cadet::util::advancePointer) 8lparen (buffer) 8comma (cadet::util::memoryForDataOf) 8lparen (p.)8tuple_at_1(8X) 8rparen 8rparen (;) \
			        		) \
			        	), \
			        	0, \
			        	8vseq_to_seq_of_tuples(8(params)) \
			        ) \
			) \
			return p; \
		} \
		inline std::size_t cacheSize() const CADET_NOEXCEPT \
		{ \
			return sizeof(params_t) + 2 * ORDER_PP(8to_lit(8seq_size(8(params)))) * sizeof(double) + 2 * ( \
				ORDER_PP(8seq_for_each_with_delimiter( \
		        	8fn(8X, 8print( \
		        		(cadet::util::memoryForDataOf) 8lparen 8cat(8quote(_), 8tuple_at_1(8X)) (.base) 8lparen 8rparen 8rparen \
		        	)), \
		        	8emit(8quote(+)), \
		        	8vseq_to_seq_of_tuples(8(params)) \
		        )) ); \
		} \
	protected: \
		inline bool validateConfig(unsigned int nComp, unsigned int const* nBoundStates); \
		ConstParams _constParams; \
		ORDER_PP(8seq_for_each( \
		        	8fn(8X, \
		        		8print( 8cat(8quote(External), 8tuple_at_0(8X)) 8space 8cat(8quote(_), 8tuple_at_1(8X)) (;)) \
		        	), \
		        	8vseq_to_seq_of_tuples(8(params)) \
		        ) \
		) \
		ORDER_PP(8seq_for_each( \
		        	8fn(8X, \
		        		8print( 8tuple_at_0(8X) 8space 8tuple_at_2(8X) (;)) \
		        	), \
		        	8vseq_to_seq_of_tuples(8(constparams)) \
		        ) \
		) \
	};


namespace cadet
{

template <typename ValType>
inline void readScalarParameterOrArray(std::vector<ValType>& dest, IParameterProvider& paramProvider, const std::string& dataSet, unsigned int nExpand);

namespace model
{

	/**
	 * @brief Base class of model parameter storage classes that do not depend on external functions
	 */
	struct ConstParamHandlerBase
	{
		/**
		 * @brief Sets external functions for this model
		 * @param [in] extFuns Pointer to array of IExternalFunction objects of size @p size
		 * @param [in] size Number of elements in the IExternalFunction array @p extFuns
		 */
		inline void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { }

		/**
		 * @brief Returns whether the model parameters depend on time
		 * @details Model parameters that do not use external functions do not depend on time.
		 * @return @c true if the model parameters depends on time, otherwise @c false
		 */
		static bool dependsOnTime() CADET_NOEXCEPT { return false; }

		/**
		 * @brief Returns how much memory is required for caching in bytes
		 * @details Memory size in bytes.
		 * @return Memory size in bytes
		 */
		inline std::size_t cacheSize() const CADET_NOEXCEPT { return 0; }
	};

	/**
	 * @brief Base class for externally dependent model parameter classes
	 * @details Configures and stores the external function used for the model parameters.
	 */
	struct ExternalParamHandlerBase
	{	
	public:

		/**
		 * @brief Sets external functions for this model
		 * @param [in] extFuns Pointer to array of IExternalFunction objects of size @p size
		 * @param [in] size Number of elements in the IExternalFunction array @p extFuns
		 */
		inline void setExternalFunctions(IExternalFunction** extFuns, int size)
		{
			_extFun.clear();
			_extFun.resize(_extFunIndex.size(), nullptr);
			for (unsigned int i = 0; i < _extFunIndex.size(); ++i)
			{
				if ((_extFunIndex[i] >= 0) && (_extFunIndex[i] < size))
					_extFun[i] = extFuns[_extFunIndex[i]];
				else
				{
					_extFun[i] = nullptr;
					LOG(Warning) << "Index " << _extFunIndex[i] << " exceeds number of passed external functions (" << size << "), external dependence is ignored";
				}
			}
		}

		/**
		 * @brief Returns whether the model parameters depend on time
		 * @details Model parameters that do not use external functions do not depend on time.
		 * @return @c true if the model parameters depends on time, otherwise @c false
		 */
		static bool dependsOnTime() CADET_NOEXCEPT { return true; }

	protected:

		std::vector<IExternalFunction*> _extFun; //!< Pointer to the external function
		std::vector<int> _extFunIndex; //!< Index to the external function

		ExternalParamHandlerBase() : _extFun(), _extFunIndex() { }
		
		/**
		 * @brief Configures the external data source of this externally dependent parameter set
		 * @param [in] paramProvider Parameter provider
		 * @param [in] nParams Number of externally dependent parameters (also size of buffer)
		 */
		inline void configure(IParameterProvider& paramProvider, unsigned int nParams)
		{			
			std::vector<int> idx;
			if (paramProvider.exists("EXTFUN"))
				idx = paramProvider.getIntArray("EXTFUN");

			if (idx.size() >= nParams)
				_extFunIndex = idx;
			else
			{
				_extFunIndex.resize(nParams);
				if (!idx.empty())
				{
					// Use one external function for all parameters
					std::fill(_extFunIndex.begin(), _extFunIndex.end(), idx[0]);
				}
				else
				{
					// There is no external dependence configured
					std::fill(_extFunIndex.begin(), _extFunIndex.end(), -1);
				}
			}
		}

		/**
		 * @brief Evaluates the external functions for the different parameters
		 * @param [in] t Current time
		 * @param [in] z Axial coordinate in the column
		 * @param [in] r Radial coordinate in the bead
		 * @param [in] secIdx Index of the current section
		 * @param [in] nParams Number of externally dependent parameters (also size of buffer)
		 * @param [out] buffer Buffer that holds function evaluations
		 */
		inline void evaluateExternalFunctions(double t, double z, double r, unsigned int secIdx, unsigned int nParams, double* buffer) const
		{
			for (unsigned int i = 0; i < nParams; ++i)
			{
				IExternalFunction* const fun = _extFun[i];
				if (fun)
					buffer[i] = fun->externalProfile(t, z, r, secIdx);
				else
					buffer[i] = 0.0;
			}
		}

		/**
		 * @brief Evaluates the time derivative of the external functions for the different parameters
		 * @param [in] t Current time
		 * @param [in] z Axial coordinate in the column
		 * @param [in] r Radial coordinate in the bead
		 * @param [in] secIdx Index of the current section
		 * @param [in] nParams Number of externally dependent parameters (also size of buffer)
		 * @param [out] buffer Buffer that holds time derivatives of each external function
		 */
		inline void evaluateTimeDerivativeExternalFunctions(double t, double z, double r, unsigned int secIdx, unsigned int nParams, double* buffer) const
		{
			for (unsigned int i = 0; i < nParams; ++i)
			{
				IExternalFunction* const fun = _extFun[i];
				if (fun)
					buffer[i] = fun->timeDerivative(t, z, r, secIdx);
				else
					buffer[i] = 0.0;
			}
		}
	};

}  // namespace model

}  // namespace cadet

#endif  // LIBCADET_EXTFUNSUPPORT_HPP_
