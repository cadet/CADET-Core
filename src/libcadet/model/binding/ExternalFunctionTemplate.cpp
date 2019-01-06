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
 * Template for external function support in binding models.
 * This file serves as a template for generating parameter handler classes that
 * hold, configure, and update parameters, which may also depend on external
 * functions.
 */

/* <codegentemplate> */
class {{ name }} : public cadet::model::ConstParamHandlerBase
{
public:
	struct ConstParams
	{
{% for p in parameters %}
	typename {{ p/type }}::storage_t {{ p/varName }};
{% endfor %}
{% if exists("constantParameters") %}
	{% for p in constantParameters %}

		{% if length(p/varName) > 1 %}

			{% for v in p/varName %}

	typename {{ p/type }}::storage_t {{ v }};

			{% endfor %}

		{% else %}

	typename {{ p/type }}::storage_t {{ p/varName }};

		{% endif %}

	{% endfor %}
{% endif %}
	};

	typedef ConstParams params_t;

	static inline const char* identifier() CADET_NOEXCEPT;

	{{ name }}() CADET_NOEXCEPT :
{% for p in parameters %}
		_{{ p/varName }}(&_localParams.{{ p/varName }})
	{% if not is_last %}
		,
	{% endif %}
{% endfor %}
{% if exists("constantParameters") %}
		,
	{% for p in constantParameters %}
		{% if length(p/varName) > 1 %}
			{{ p/objName }}({% for v in p/varName %} &_localParams.{{ v }} {% if not is_last %},{% endif %} {% endfor %})
		{% else %}
			_{{ p/varName }}(&_localParams.{{ p/varName }})
		{% endif %}
		{% if not is_last %},{% endif %}
	{% endfor %}
{% endif %}
	{ }

	inline bool configure(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBoundStates)
	{
{% for p in parameters %}
		_{{ p/varName }}.configure("{{ p/confName }}", paramProvider, nComp, nBoundStates);
{% endfor %}
{% if exists("constantParameters") %}
	{% for p in constantParameters %}
		{% if length(p/varName) > 1 %}
			{{ p/objName }}.configure("{{ p/confPrefix }}", paramProvider, nComp, nBoundStates);
		{% else %}
			_{{ p/varName }}.configure("{{ p/confName }}", paramProvider, nComp, nBoundStates);
		{% endif %}
	{% endfor %}
{% endif %}
		return validateConfig(nComp, nBoundStates);
	}

	inline void registerParameters(std::unordered_map<ParameterId, active*>& parameters, unsigned int unitOpIdx, unsigned int parTypeIdx, unsigned int nComp, unsigned int const* nBoundStates)
	{
{% for p in parameters %}
		_{{ p/varName }}.registerParam("{{ p/confName }}", parameters, unitOpIdx, parTypeIdx, nComp, nBoundStates);
{% endfor %}
{% if exists("constantParameters") %}
	{% for p in constantParameters %}
		{% if length(p/varName) > 1 %}
			{{ p/objName }}.registerParam("{{ p/confPrefix }}", parameters, unitOpIdx, parTypeIdx, nComp, nBoundStates);
		{% else %}
			_{{ p/varName }}.registerParam("{{ p/confName }}", parameters, unitOpIdx, parTypeIdx, nComp, nBoundStates);
		{% endif %}
	{% endfor %}
{% endif %}
	}

	inline void reserve(unsigned int numElem, unsigned int numSlices, unsigned int nComp, unsigned int const* nBoundStates)
	{
{% for p in parameters %}
		_{{ p/varName }}.reserve(numElem, numSlices, nComp, nBoundStates);
{% endfor %}
{% if exists("constantParameters") %}
	{% for p in constantParameters %}
		{% if length(p/varName) > 1 %}
			{{ p/objName }}.reserve(numElem, numSlices, nComp, nBoundStates);
		{% else %}
			_{{ p/varName }}.reserve(numElem, numSlices, nComp, nBoundStates);
		{% endif %}
	{% endfor %}
{% endif %}
	}

	inline const params_t& update(double t, double z, double r, unsigned int secIdx, unsigned int nComp, unsigned int const* nBoundStates, void* workSpace) const
	{
		return _localParams;
	}

	inline const params_t updateTimeDerivative(double t, double z, double r, unsigned int secIdx, unsigned int nComp, unsigned int const* nBoundStates, void* workSpace) const
	{
		return _localParams;
	}

protected:
	inline bool validateConfig(unsigned int nComp, unsigned int const* nBoundStates);
	
	ConstParams _localParams;

{% for p in parameters %}
		{{ p/type }} _{{ p/varName }};
{% endfor %}
{% if exists("constantParameters") %}
	{% for p in constantParameters %}
		{% if length(p/varName) > 1 %}
			{{ p/type }} {{ p/objName }};
		{% else %}
			{{ p/type }} _{{ p/varName }};
		{% endif %}
	{% endfor %}
{% endif %}
};

class {{ externalName }} : public cadet::model::ExternalParamHandlerBase
{
public:
	struct ConstParams
	{
{% if exists("constantParameters") %}
	{% for p in constantParameters %}

		{% if length(p/varName) > 1 %}

			{% for v in p/varName %}

	typename {{ p/type }}::storage_t {{ v }};

			{% endfor %}

		{% else %}

	typename {{ p/type }}::storage_t {{ p/varName }};

		{% endif %}

	{% endfor %}
{% endif %}
	};

	struct VariableParams
	{
{% for p in parameters %}
	typename util::localVersionOf<{{ p/type }}::storage_t>::type {{ p/varName }};
{% endfor %}
{% if exists("constantParameters") %}
	{% for p in constantParameters %}

		{% if length(p/varName) > 1 %}

			{% for v in p/varName %}

	typename {{ p/type }}::storage_t {{ v }};

			{% endfor %}

		{% else %}

	typename {{ p/type }}::storage_t {{ p/varName }};

		{% endif %}

	{% endfor %}
{% endif %}
	};
	
	typedef VariableParams params_t;
	static inline const char* identifier() CADET_NOEXCEPT;
	
	{{ externalName }}() CADET_NOEXCEPT
{% if exists("constantParameters") %}
		:
	{% for p in constantParameters %}
		{% if length(p/varName) > 1 %}
			{{ p/objName }}({% for v in p/varName %} &_constParams.{{ v }} {% if not is_last %},{% endif %} {% endfor %})
		{% else %}
			_{{ p/varName }}(&_constParams.{{ p/varName }})
		{% endif %}
		{% if not is_last %},{% endif %}
	{% endfor %}
{% endif %}
	{ }

	inline bool configure(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBoundStates)
	{
{% for p in parameters %}
		_{{ p/varName }}.configure("{{ p/confName }}", paramProvider, nComp, nBoundStates);
{% endfor %}
{% if exists("constantParameters") %}
	{% for p in constantParameters %}
		{% if length(p/varName) > 1 %}
			{{ p/objName }}.configure("{{ p/confPrefix }}", paramProvider, nComp, nBoundStates);
		{% else %}
			_{{ p/varName }}.configure("{{ p/confName }}", paramProvider, nComp, nBoundStates);
		{% endif %}
	{% endfor %}
{% endif %}
		ExternalParamHandlerBase::configure(paramProvider, {{ length(parameters) }});
		return validateConfig(nComp, nBoundStates);
	}

	inline void registerParameters(std::unordered_map<ParameterId, active*>& parameters, unsigned int unitOpIdx, unsigned int parTypeIdx, unsigned int nComp, unsigned int const* nBoundStates)
	{
{% for p in parameters %}
		_{{ p/varName }}.registerParam("{{ p/confName }}", parameters, unitOpIdx, parTypeIdx, nComp, nBoundStates);
{% endfor %}
{% if exists("constantParameters") %}
	{% for p in constantParameters %}
		{% if length(p/varName) > 1 %}
			{{ p/objName }}.registerParam("{{ p/confPrefix }}", parameters, unitOpIdx, parTypeIdx, nComp, nBoundStates);
		{% else %}
			_{{ p/varName }}.registerParam("{{ p/confName }}", parameters, unitOpIdx, parTypeIdx, nComp, nBoundStates);
		{% endif %}
	{% endfor %}
{% endif %}
	}

	inline void reserve(unsigned int numElem, unsigned int numSlices, unsigned int nComp, unsigned int const* nBoundStates)
	{
{% for p in parameters %}
		_{{ p/varName }}.reserve(numElem, numSlices, nComp, nBoundStates);
{% endfor %}
{% if exists("constantParameters") %}
	{% for p in constantParameters %}
		{% if length(p/varName) > 1 %}
			{{ p/objName }}.reserve(numElem, numSlices, nComp, nBoundStates);
		{% else %}
			_{{ p/varName }}.reserve(numElem, numSlices, nComp, nBoundStates);
		{% endif %}
	{% endfor %}
{% endif %}
	}

	inline const params_t& update(double t, double z, double r, unsigned int secIdx, unsigned int nComp, unsigned int const* nBoundStates, void* workSpace) const
	{
		params_t* const localParams = reinterpret_cast<params_t*>(workSpace);
		new (localParams) params_t;
		double* const extFunBuffer = cadet::util::advancePointer<double>(workSpace, sizeof(params_t));
		void* buffer = cadet::util::advancePointer<void>(workSpace, sizeof(params_t) + 2 * {{ length(parameters) }} * sizeof(double));
		evaluateExternalFunctions(t, z, r, secIdx, {{ length(parameters) }}, extFunBuffer);
{% if exists("constantParameters") %}
	{% for p in constantParameters %}
		{% if length(p/varName) > 1 %}
			{% for v in p/varName %}
				localParams->{{ v }} = _constParams.{{ v }};
			{% endfor %}
		{% else %}
			localParams->{{ p/varName }} = _constParams.{{ p/varName }};
		{% endif %}
	{% endfor %}
{% endif %}
{% for p in parameters %}
		_{{ p/varName }}.prepareCache(localParams->{{ p/varName }}, buffer);
		_{{ p/varName }}.update(cadet::util::dataOfLocalVersion(localParams->{{ p/varName }}), extFunBuffer[{{ index }}], nComp, nBoundStates);
		{% if not is_last %}buffer = cadet::util::advancePointer(buffer, cadet::util::memoryForDataOf(localParams->{{ p/varName }}));{% endif %}
{% endfor %}
		return *localParams;
	}

	inline params_t updateTimeDerivative(double t, double z, double r, unsigned int secIdx, unsigned int nComp, unsigned int const* nBoundStates, void* workSpace) const
	{
		VariableParams p;
		params_t* const localParams = reinterpret_cast<params_t*>(workSpace);
		double* const extFunBuffer = cadet::util::advancePointer<double>(workSpace, sizeof(params_t));
		double* const extDerivBuffer = extFunBuffer + {{ length(parameters) }};
		void* buffer = util::ptrToEndOfData(localParams->{% for p in parameters %} {% if index1 == length(parameters) %} {{ p/varName }} {% endif %} {% endfor %});
		evaluateTimeDerivativeExternalFunctions(t, z, r, secIdx, {{ length(parameters) }}, extDerivBuffer);

{% if exists("constantParameters") %}
	{% for p in constantParameters %}
		{% if length(p/varName) > 1 %}
			{% for v in p/varName %}
				p.{{ v }} = _constParams.{{ v }};
			{% endfor %}
		{% else %}
			p.{{ p/varName }} = _constParams.{{ p/varName }};
		{% endif %}
	{% endfor %}
{% endif %}
{% for p in parameters %}
		_{{ p/varName }}.prepareCache(p.{{ p/varName }}, buffer);
		_{{ p/varName }}.updateTimeDerivative(cadet::util::dataOfLocalVersion(p.{{ p/varName }}), extFunBuffer[{{ index }}], extDerivBuffer[{{ index }}], nComp, nBoundStates);
		{% if not is_last %}buffer = cadet::util::advancePointer(buffer, cadet::util::memoryForDataOf(p.{{ p/varName }}));{% endif %}
{% endfor %}
		return p;
	}

	inline std::size_t cacheSize(unsigned int nComp, unsigned int totalNumBoundStates, unsigned int const* nBoundStates) const CADET_NOEXCEPT
	{
		return sizeof(params_t) + 2 * {{ length(parameters) }} * sizeof(double) + 2 * (
{% for p in parameters %}
		_{{ p/varName }}.additionalDynamicMemory(nComp, totalNumBoundStates, nBoundStates) {% if not is_last %} + {% endif %}
{% endfor %}
		);
	}

protected:
	inline bool validateConfig(unsigned int nComp, unsigned int const* nBoundStates);
	
	ConstParams _constParams;

{% for p in parameters %}
		External{{ p/type }} _{{ p/varName }};
{% endfor %}
{% if exists("constantParameters") %}
	{% for p in constantParameters %}
		{% if length(p/varName) > 1 %}
			{{ p/type }} {{ p/objName }};
		{% else %}
			{{ p/type }} _{{ p/varName }};
		{% endif %}
	{% endfor %}
{% endif %}
};
