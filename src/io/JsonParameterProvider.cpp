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

#include <json.hpp>

#include <sstream>
#include <fstream>
#include <iomanip>

#define CADET_JSONPARAMETERPROVIDER_NOFORWARD
#include "common/JsonParameterProvider.hpp"

using json = nlohmann::json;

namespace cadet
{

JsonParameterProvider::JsonParameterProvider() : _root(nullptr)
{
}

JsonParameterProvider::JsonParameterProvider(const char* data) : _root(new json(json::parse(data)))
{
	_opened.push(_root);
}

JsonParameterProvider::JsonParameterProvider(const std::string& data) : _root(new json(json::parse(data)))
{
	_opened.push(_root);
}

JsonParameterProvider::JsonParameterProvider(const json& data) : _root(new json(data))
{
	_opened.push(_root);
}

JsonParameterProvider::JsonParameterProvider(json* data) : _root(data)
{
	_opened.push(_root);
}

JsonParameterProvider::JsonParameterProvider(const JsonParameterProvider& cpy)
{
	_root = new json(*cpy._root);
	_opened = cpy._opened;
}

JsonParameterProvider::JsonParameterProvider(JsonParameterProvider&& cpy) CADET_NOEXCEPT : _root(cpy._root), _opened(std::move(cpy._opened))
{
	cpy._root = nullptr;
	cpy._opened = std::stack<json*>();
}

JsonParameterProvider::~JsonParameterProvider() CADET_NOEXCEPT
{
	delete _root;
}

JsonParameterProvider& JsonParameterProvider::operator=(const JsonParameterProvider& cpy)
{
	delete _root;

	_root = new json(*cpy._root);
	_opened = cpy._opened;
	return *this;
}

#ifdef COMPILER_SUPPORT_NOEXCEPT_DEFAULTED_MOVE
	JsonParameterProvider& JsonParameterProvider::operator=(JsonParameterProvider&& cpy) CADET_NOEXCEPT
#else
	JsonParameterProvider& JsonParameterProvider::operator=(JsonParameterProvider&& cpy)
#endif
{
	delete _root;
	_opened = std::stack<nlohmann::json*>();
	_root = cpy._root;
	cpy._root = nullptr;
	cpy._opened = std::stack<nlohmann::json*>();

	_opened.push(_root);
	return *this;
}

double JsonParameterProvider::getDouble(const std::string& paramName)
{
	return _opened.top()->at(paramName).get<double>();
}

int JsonParameterProvider::getInt(const std::string& paramName)
{
	const json p = _opened.top()->at(paramName);
	if (p.is_boolean())
		return p.get<bool>();
	else
		return p.get<int>();
}

uint64_t JsonParameterProvider::getUint64(const std::string& paramName)
{
	return _opened.top()->at(paramName).get<uint64_t>();
}

bool JsonParameterProvider::getBool(const std::string& paramName)
{
	const json p = _opened.top()->at(paramName);
	if (p.is_number_integer())
		return p.get<int>();
	else
		return p.get<bool>();
}

std::string JsonParameterProvider::getString(const std::string& paramName)
{
	return _opened.top()->at(paramName).get<std::string>();
}

std::vector<double> JsonParameterProvider::getDoubleArray(const std::string& paramName)
{
	return _opened.top()->at(paramName).get<std::vector<double>>();
}

std::vector<int> JsonParameterProvider::getIntArray(const std::string& paramName)
{
	return _opened.top()->at(paramName).get<std::vector<int>>();
}

std::vector<uint64_t> JsonParameterProvider::getUint64Array(const std::string& paramName)
{
	return _opened.top()->at(paramName).get<std::vector<uint64_t>>();
}

std::vector<bool> JsonParameterProvider::getBoolArray(const std::string& paramName)
{
	return _opened.top()->at(paramName);
}

std::vector<std::string> JsonParameterProvider::getStringArray(const std::string& paramName)
{
	return _opened.top()->at(paramName).get<std::vector<std::string>>();
}

bool JsonParameterProvider::exists(const std::string& paramName)
{
	return _opened.top()->find(paramName) != _opened.top()->end();
}

bool JsonParameterProvider::isArray(const std::string& paramName)
{
	return _opened.top()->at(paramName).is_array();
}

void JsonParameterProvider::pushScope(const std::string& scope)
{
	_opened.push(&_opened.top()->at(scope));
}

void JsonParameterProvider::popScope()
{
	_opened.pop();
}

void JsonParameterProvider::addScope(const std::string& scope)
{
	if (!exists(scope))
	{
		json j;
		j["blubber"] = 0.0;
		(*_opened.top())[scope] = j;
	}
}

void JsonParameterProvider::set(const std::string& paramName, double val)
{
	(*_opened.top())[paramName] = val;
}

void JsonParameterProvider::set(const std::string& paramName, int val)
{
	(*_opened.top())[paramName] = val;
}

void JsonParameterProvider::set(const std::string& paramName, uint64_t val)
{
	(*_opened.top())[paramName] = val;
}

void JsonParameterProvider::set(const std::string& paramName, bool val)
{
	(*_opened.top())[paramName] = val;
}

void JsonParameterProvider::set(const std::string& paramName, char const* val)
{
	(*_opened.top())[paramName] = std::string(val);
}

void JsonParameterProvider::set(const std::string& paramName, const std::string& val)
{
	(*_opened.top())[paramName] = val;
}

void JsonParameterProvider::set(const std::string& paramName, const std::vector<double>& val)
{
	(*_opened.top())[paramName] = val;
}

void JsonParameterProvider::set(const std::string& paramName, const std::vector<int>& val)
{
	(*_opened.top())[paramName] = val;
}

void JsonParameterProvider::set(const std::string& paramName, const std::vector<uint64_t>& val)
{
	(*_opened.top())[paramName] = val;
}

void JsonParameterProvider::set(const std::string& paramName, const std::vector<std::string>& val)
{
	(*_opened.top())[paramName] = val;
}

JsonParameterProvider JsonParameterProvider::fromFile(const std::string& fileName)
{
	std::ifstream ifs(fileName);
	json* root = new json();
	ifs >> (*root);

	return JsonParameterProvider(root);
}

std::ostream& operator<<(std::ostream& out, const JsonParameterProvider& jpp)
{
	out << (*jpp.data());
	return out;
}

} // namespace cadet
