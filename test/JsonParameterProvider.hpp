// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2017: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

/**
 * @file 
 * Defines a ParameterProvider that uses JSON.
 */

#ifndef CADETTEST_JSONPARAMETERPROVIDER_HPP_
#define CADETTEST_JSONPARAMETERPROVIDER_HPP_

#include "cadet/ParameterProvider.hpp"

#include <string>
#include <stack>
#include <ostream>

#ifndef CADETTEST_JSONPARAMETERPROVIDER_NOFORWARD
	namespace nlohmann
	{
		class json;
	}
#endif

namespace cadet
{

class JsonParameterProvider : public cadet::IParameterProvider
{
public:

	JsonParameterProvider(const char* data);
	JsonParameterProvider(const std::string& data);
	JsonParameterProvider(const nlohmann::json& data);
	JsonParameterProvider(const JsonParameterProvider& cpy);
	JsonParameterProvider(JsonParameterProvider&& cpy) CADET_NOEXCEPT;

	virtual ~JsonParameterProvider() CADET_NOEXCEPT;

	JsonParameterProvider& operator=(const JsonParameterProvider& cpy);

#ifdef COMPILER_SUPPORT_NOEXCEPT_DEFAULTED_MOVE
	JsonParameterProvider& operator=(JsonParameterProvider&& cpy) CADET_NOEXCEPT;
#else
	JsonParameterProvider& operator=(JsonParameterProvider&& cpy);
#endif

	virtual double getDouble(const std::string& paramName);
	virtual int getInt(const std::string& paramName);
	virtual uint64_t getUint64(const std::string& paramName);
	virtual bool getBool(const std::string& paramName);
	virtual std::string getString(const std::string& paramName);
	virtual std::vector<double> getDoubleArray(const std::string& paramName);
	virtual std::vector<int> getIntArray(const std::string& paramName);
	virtual std::vector<uint64_t> getUint64Array(const std::string& paramName);
	virtual std::vector<bool> getBoolArray(const std::string& paramName);
	virtual std::vector<std::string> getStringArray(const std::string& paramName);
	virtual bool exists(const std::string& paramName);
	virtual bool isArray(const std::string& paramName);
	virtual void pushScope(const std::string& scope);
	virtual void popScope();

	virtual void addScope(const std::string& scope);

	void set(const std::string& paramName, double val);
	void set(const std::string& paramName, int val);
	void set(const std::string& paramName, uint64_t val);
	void set(const std::string& paramName, bool val);
	void set(const std::string& paramName, char const* val);
	void set(const std::string& paramName, const std::string& val);
	void set(const std::string& paramName, const std::vector<double>& val);
	void set(const std::string& paramName, const std::vector<int>& val);
	void set(const std::string& paramName, const std::vector<uint64_t>& val);
	void set(const std::string& paramName, const std::vector<std::string>& val);

	inline nlohmann::json* data() { return _root; }
	inline nlohmann::json const* data() const { return _root; }
private:
	nlohmann::json* _root;
	std::stack<nlohmann::json*> _opened;
};

std::ostream& operator<<(std::ostream& out, const JsonParameterProvider& jpp);
} // namespace cadet

cadet::JsonParameterProvider createGRMwithSMA();
cadet::JsonParameterProvider createGRMwithLinear();
cadet::JsonParameterProvider createLWE();
cadet::JsonParameterProvider createLinearBenchmark(bool dynamicBinding, bool nonBinding);
cadet::JsonParameterProvider createCSTR(unsigned int nComp);
cadet::JsonParameterProvider createCSTRBenchmark(unsigned int nSec, double endTime, double interval);

#endif  // CADETTEST_JSONPARAMETERPROVIDER_HPP_
