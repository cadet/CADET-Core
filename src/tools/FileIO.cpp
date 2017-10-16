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

#include <sstream>
#include <iomanip>
#include <stack>
#include <iterator>

#include "io/FileIO.hpp"
#include "io/hdf5/HDF5Reader.hpp"
#include "io/hdf5/HDF5Writer.hpp"
#include "io/xml/XMLReader.hpp"
#include "io/xml/XMLWriter.hpp"
#include "cadet/StringUtil.hpp"

#include <json.hpp>

using json = nlohmann::json;


namespace
{

template <typename reader_t>
class FileReader : public cadet::io::IFileReader
{
public:
	FileReader() { }
	virtual ~FileReader() CADET_NOEXCEPT { }

	virtual double getDouble(const std::string& paramName)
	{
		return _io.template scalar<double>(paramName);
	}

	virtual int getInt(const std::string& paramName)
	{
		return _io.template scalar<int>(paramName);
	}

	virtual bool getBool(const std::string& paramName)
	{
		return _io.template scalar<int>(paramName);
	}

	virtual std::string getString(const std::string& paramName)
	{
		return _io.template scalar<std::string>(paramName);
	}

	virtual std::vector<double> getDoubleArray(const std::string& paramName)
	{
		return _io.template vector<double>(paramName);
	}

	virtual std::vector<int> getIntArray(const std::string& paramName)
	{
		return _io.template vector<int>(paramName);
	}

	virtual std::vector<bool> getBoolArray(const std::string& paramName)
	{
		const std::vector<int> data = _io.template vector<int>(paramName);
		std::vector<bool> bd(data.size());
		for (unsigned int i = 0; i < data.size(); ++i)
			bd[i] = data[i];

		return bd;
	}

	virtual std::vector<std::string> getStringArray(const std::string& paramName)
	{
		return _io.template vector<std::string>(paramName);
	}

	virtual bool isArray(const std::string& paramName)
	{
		return _io.isVector(paramName);
	}

	virtual bool isString(const std::string& elementName)
	{
		return _io.isString(elementName);
	}

	virtual bool isInt(const std::string& elementName)
	{
		return _io.isInt(elementName);
	}

	virtual bool isDouble(const std::string& elementName)
	{
		return _io.isDouble(elementName);
	}

	virtual std::vector<std::size_t> tensorDimensions(const std::string& elementName)
	{
		return _io.tensorDimensions(elementName);
	}

protected:
	reader_t _io;
};


template <typename writer_t>
class FileWriter : public cadet::io::IFileWriter
{
public:
	FileWriter() { }
	virtual ~FileWriter() CADET_NOEXCEPT { }
	
	virtual void writeTensorDouble(const std::string& dataSetName, const std::size_t rank, const std::size_t* dims, const double* buffer, const std::size_t stride)
	{
		_io.template tensor<double>(dataSetName, rank, dims, buffer, stride);
	}

	virtual void writeTensorInt(const std::string& dataSetName, const std::size_t rank, const std::size_t* dims, const int* buffer, const std::size_t stride)
	{
		_io.template tensor<int>(dataSetName, rank, dims, buffer, stride);
	}

	virtual void writeTensorBool(const std::string& dataSetName, const std::size_t rank, const std::size_t* dims, const bool* buffer, const std::size_t stride)
	{
		// Compute buffer size
		std::size_t bufSize = 1;
		for (std::size_t i = 0; i < rank; ++i)
		{
			bufSize *= dims[i];
		}

		// Convert buffer to int
		std::vector<int> bd(bufSize);
		for (std::size_t i = 0; i < bufSize; ++i)
			bd[i] = (buffer[i] ? 1 : 0);

		_io.template tensor<int>(dataSetName, rank, dims, bd.data(), stride);
	}

	virtual void writeTensorString(const std::string& dataSetName, const std::size_t rank, const std::size_t* dims, const std::string* buffer, const std::size_t stride)
	{
		_io.template tensor<std::string>(dataSetName, rank, dims, buffer, stride);
	}

	virtual void writeMatrixDouble(const std::string& dataSetName, const std::size_t rows, const std::size_t cols, const double* buffer, const std::size_t stride)
	{
		_io.template matrix<double>(dataSetName, rows, cols, buffer, stride);
	}

	virtual void writeMatrixInt(const std::string& dataSetName, const std::size_t rows, const std::size_t cols, const int* buffer, const std::size_t stride)
	{
		_io.template matrix<int>(dataSetName, rows, cols, buffer, stride);
	}

	virtual void writeMatrixBool(const std::string& dataSetName, const std::size_t rows, const std::size_t cols, const bool* buffer, const std::size_t stride)
	{
		// Compute buffer size
		const std::size_t bufSize = rows * cols;

		// Convert buffer to int
		std::vector<int> bd(bufSize);
		for (std::size_t i = 0; i < bufSize; ++i)
			bd[i] = (buffer[i] ? 1 : 0);

		_io.template matrix<int>(dataSetName, rows, cols, bd.data(), stride);
	}

	virtual void writeMatrixString(const std::string& dataSetName, const std::size_t rows, const std::size_t cols, const std::string* buffer, const std::size_t stride)
	{
		_io.template matrix<std::string>(dataSetName, rows, cols, buffer, stride);
	}

	virtual void writeVectorDouble(const std::string& dataSetName, const std::size_t length, const double* buffer, const std::size_t stride)
	{
		_io.template vector<double>(dataSetName, length, buffer, stride);
	}

	virtual void writeVectorInt(const std::string& dataSetName, const std::size_t length, const int* buffer, const std::size_t stride)
	{
		_io.template vector<int>(dataSetName, length, buffer, stride);
	}

	virtual void writeVectorBool(const std::string& dataSetName, const std::size_t length, const bool* buffer, const std::size_t stride)
	{
		// Convert buffer to int
		std::vector<int> bd(length);
		for (std::size_t i = 0; i < length; ++i)
			bd[i] = (buffer[i] ? 1 : 0);

		_io.template vector<int>(dataSetName, length, bd.data(), stride);
	}

	virtual void writeVectorString(const std::string& dataSetName, const std::size_t length, const std::string* buffer, const std::size_t stride)
	{
		_io.template vector<std::string>(dataSetName, length, buffer, stride);
	}

	virtual void writeDouble(const std::string& dataSetName, const double buffer)
	{
		_io.template scalar<double>(dataSetName, buffer);
	}

	virtual void writeInt(const std::string& dataSetName, const int buffer)
	{
		_io.template scalar<int>(dataSetName, buffer);
	}

	virtual void writeBool(const std::string& dataSetName, const bool buffer)
	{
		_io.template scalar<int>(dataSetName, buffer);
	}

	virtual void writeString(const std::string& dataSetName, const std::string& buffer)
	{
		_io.template scalar<std::string>(dataSetName, buffer);
	}

	virtual void deleteGroup(const std::string& groupName)
	{
		_io.unlinkGroup(groupName);
	}

	virtual void deleteDataset(const std::string& dsName)
	{
		_io.unlinkDataset(dsName);
	}

protected:
	writer_t _io;
};


template <typename base_t>
class BaseIOWrapper : public base_t
{
public:
	BaseIOWrapper() { }

	virtual void openFile(const std::string& fileName, const char* mode)
	{
		base_t::_io.openFile(fileName, mode);
	}

	virtual void closeFile()
	{
		base_t::_io.closeFile();
	}

	virtual void pushGroup(const std::string& scope)
	{
		base_t::_io.pushGroup(scope);
	}

	virtual void popGroup()
	{
		base_t::_io.popGroup();
	}

	virtual void setGroup(const std::string& scope)
	{
		base_t::_io.setGroup(scope);
	}

	virtual bool exists(const std::string& paramName)
	{
		return base_t::_io.exists(paramName);
	}

	virtual bool isGroup(const std::string& elementName)
	{
		return base_t::_io.isGroup(elementName);
	}

	virtual int numItems()
	{
		return base_t::_io.numItems();
	}

	virtual std::string itemName(int n)
	{
		return base_t::_io.itemName(n);
	}

	virtual std::vector<std::string> itemNames()
	{
		return base_t::_io.itemNames();
	}
};


class JSONMemoryReaderWriter : public cadet::io::IMemoryIO
{
public:

	JSONMemoryReaderWriter() : _root(new json())
	{
		_opened.push(_root);
	}

	virtual ~JSONMemoryReaderWriter() CADET_NOEXCEPT
	{
		delete _root;
	}

	virtual void pushGroup(const std::string& scope)
	{
		if (exists(scope))
			_opened.push(&_opened.top()->at(scope));
		else
		{
			json j;
			j["blubber"] = 0.0;
			(*_opened.top())[scope] = j;
			_opened.push(&_opened.top()->at(scope));
		}
	}

	virtual void popGroup()
	{
		_opened.pop();
	}

	virtual void setGroup(const std::string& scope)
	{
		std::size_t start = 0;
		std::size_t end = 0;
		const std::string delimiter("/");

		// Quick return when called with empty group name
		if (scope.empty() || (scope == "/"))
		{
			_opened = std::stack<nlohmann::json*>();
			_opened.push(_root);
			return;
		}

		// Don't care for a preceding delimiter
		if (scope[0] == delimiter[0]) 
			++start;

		while (end != std::string::npos)
		{
			end = scope.find(delimiter, start);

			// If at end, use length = maxLength.  Else use length = end - start.
			pushGroup(scope.substr(start, (end == std::string::npos) ? std::string::npos : end - start));

			// If at end, use start = maxSize.  Else use start = end + delimiter.
			start = ((end > (std::string::npos - delimiter.size())) ? std::string::npos : end + delimiter.size());
		}		
	}

	virtual bool exists(const std::string& paramName)
	{
		return _opened.top()->find(paramName) != _opened.top()->end();
	}

	virtual double getDouble(const std::string& paramName)
	{
		return _opened.top()->at(paramName).get<double>();
	}

	virtual int getInt(const std::string& paramName)
	{
		const json p = _opened.top()->at(paramName);
		if (p.is_boolean())
			return p.get<bool>();
		else
			return p.get<int>();
	}

	virtual bool getBool(const std::string& paramName)
	{
		const json p = _opened.top()->at(paramName);
		if (p.is_number_integer())
			return p.get<int>();
		else
			return p.get<bool>();
	}

	virtual std::string getString(const std::string& paramName)
	{
		return _opened.top()->at(paramName).get<std::string>();
	}

	virtual std::vector<double> getDoubleArray(const std::string& paramName)
	{
		return _opened.top()->at(paramName).get<std::vector<double>>();
	}

	virtual std::vector<int> getIntArray(const std::string& paramName)
	{
		return _opened.top()->at(paramName).get<std::vector<int>>();
	}

	virtual std::vector<bool> getBoolArray(const std::string& paramName)
	{
		const json p = _opened.top()->at(paramName);
		if (p.is_number_integer())
		{
			const std::vector<int> d = p.get<std::vector<int>>();
			std::vector<bool> bd(d.size());
			for (std::size_t i = 0; i < d.size(); ++i)
				bd[i] = d[i];

			return bd;
		}
		else
			return p.get<std::vector<bool>>();
	}

	virtual std::vector<std::string> getStringArray(const std::string& paramName)
	{
		return _opened.top()->at(paramName).get<std::vector<std::string>>();
	}

	virtual bool isArray(const std::string& elementName)
	{
		return _opened.top()->at(elementName).is_array();
	}

	virtual bool isString(const std::string& elementName)
	{
		return _opened.top()->at(elementName).is_string();
	}

	virtual bool isInt(const std::string& elementName)
	{
		return _opened.top()->at(elementName).is_number_integer();
	}

	virtual bool isDouble(const std::string& elementName)
	{
		return _opened.top()->at(elementName).is_number_float();
	}

	virtual bool isGroup(const std::string& elementName)
	{
		return _opened.top()->at(elementName).is_object();
	}

	virtual int numItems()
	{
		return _opened.top()->size();
	}

	virtual std::string itemName(int n)
	{
		json::iterator it = _opened.top()->begin();
		std::advance(it, n);
		return it.key();
	}

	virtual std::vector<std::string> itemNames()
	{
		std::vector<std::string> names;
		names.reserve(_opened.top()->size());
		for (json::iterator it = _opened.top()->begin(); it != _opened.top()->end(); ++it)
		names.push_back(it.key());
		return names;
	}

	virtual std::vector<std::size_t> tensorDimensions(const std::string& elementName)
	{
		if (_opened.top()->find(elementName) == _opened.top()->end())
			return std::vector<std::size_t>();

		if (exists(elementName + "_rank") && exists(elementName + "_dims"))
		{
			return _opened.top()->at(elementName + "_dims").get<std::vector<std::size_t>>();
		}
		else if (exists(elementName + "_cols") && exists(elementName + "_rows"))
		{
			const std::size_t cols = _opened.top()->at(elementName + "_cols").get<std::size_t>();
			const std::size_t rows = _opened.top()->at(elementName + "_rows").get<std::size_t>();
			return std::vector<std::size_t>({cols, rows});
		}

		const json p = _opened.top()->at(elementName);
		if (p.is_array())
		{
			if (p.is_number_float())
				return std::vector<std::size_t>({ p.get<std::vector<double>>().size() });
			else if (p.is_number_integer())
				return std::vector<std::size_t>({ p.get<std::vector<int>>().size() });
			else if (p.is_string())
				return std::vector<std::size_t>({ p.get<std::string>().size() });
			else if (p.is_boolean())
				return std::vector<std::size_t>({ p.get<std::vector<bool>>().size() });

			return std::vector<std::size_t>({ 0 });
		}
		else
			return std::vector<std::size_t>({ 1 });
	}

	virtual void writeTensorDouble(const std::string& dataSetName, const std::size_t rank, const std::size_t* dims, const double* buffer, const std::size_t stride)
	{
		writeTensor<double, double>(dataSetName, rank, dims, buffer, stride);
	}

	virtual void writeTensorInt(const std::string& dataSetName, const std::size_t rank, const std::size_t* dims, const int* buffer, const std::size_t stride)
	{
		writeTensor<int, int>(dataSetName, rank, dims, buffer, stride);
	}

	virtual void writeTensorBool(const std::string& dataSetName, const std::size_t rank, const std::size_t* dims, const bool* buffer, const std::size_t stride)
	{
		writeTensor<bool, int>(dataSetName, rank, dims, buffer, stride);
	}

	virtual void writeTensorString(const std::string& dataSetName, const std::size_t rank, const std::size_t* dims, const std::string* buffer, const std::size_t stride)
	{
		writeTensor<std::string, std::string>(dataSetName, rank, dims, buffer, stride);
	}

	virtual void writeMatrixDouble(const std::string& dataSetName, const std::size_t rows, const std::size_t cols, const double* buffer, const std::size_t stride)
	{
		writeMatrix<double, double>(dataSetName, rows, cols, buffer, stride);
	}

	virtual void writeMatrixInt(const std::string& dataSetName, const std::size_t rows, const std::size_t cols, const int* buffer, const std::size_t stride)
	{
		writeMatrix<int, int>(dataSetName, rows, cols, buffer, stride);
	}

	virtual void writeMatrixBool(const std::string& dataSetName, const std::size_t rows, const std::size_t cols, const bool* buffer, const std::size_t stride)
	{
		writeMatrix<bool, int>(dataSetName, rows, cols, buffer, stride);
	}

	virtual void writeMatrixString(const std::string& dataSetName, const std::size_t rows, const std::size_t cols, const std::string* buffer, const std::size_t stride)
	{
		writeMatrix<std::string, std::string>(dataSetName, rows, cols, buffer, stride);
	}

	virtual void writeVectorDouble(const std::string& dataSetName, const std::size_t length, const double* buffer, const std::size_t stride)
	{
		writeVector<double, double>(dataSetName, length, buffer, stride);
	}

	virtual void writeVectorInt(const std::string& dataSetName, const std::size_t length, const int* buffer, const std::size_t stride)
	{
		writeVector<int, int>(dataSetName, length, buffer, stride);
	}

	virtual void writeVectorBool(const std::string& dataSetName, const std::size_t length, const bool* buffer, const std::size_t stride)
	{
		writeVector<bool, int>(dataSetName, length, buffer, stride);
	}

	virtual void writeVectorString(const std::string& dataSetName, const std::size_t length, const std::string* buffer, const std::size_t stride)
	{
		writeVector<std::string, std::string>(dataSetName, length, buffer, stride);
	}

	virtual void writeDouble(const std::string& dataSetName, const double buffer)
	{
		(*_opened.top())[dataSetName] = buffer;
	}

	virtual void writeInt(const std::string& dataSetName, const int buffer)
	{
		(*_opened.top())[dataSetName] = buffer;
	}

	virtual void writeBool(const std::string& dataSetName, const bool buffer)
	{
		(*_opened.top())[dataSetName] = buffer;
	}

	virtual void writeString(const std::string& dataSetName, const std::string& buffer)
	{
		(*_opened.top())[dataSetName] = buffer;
	}

	virtual void deleteGroup(const std::string& groupName)
	{
		_opened.top()->erase(groupName);
	}

	virtual void deleteDataset(const std::string& dsName)
	{
		_opened.top()->erase(dsName);
	}

	inline nlohmann::json* data() { return _root; }
	inline nlohmann::json const* data() const { return _root; }

protected:

	template <typename source_t, typename target_t>
	void writeVector(const std::string& dataSetName, const std::size_t length, const source_t* buffer, const std::size_t stride)
	{
		std::vector<target_t> d(length);
		for (std::size_t i = 0; i < length; ++i)
			d[i] = static_cast<target_t>(buffer[i * stride]);

		(*_opened.top())[dataSetName] = d;
	}

	template <typename source_t, typename target_t>
	void writeMatrix(const std::string& dataSetName, const std::size_t rows, const std::size_t cols, const source_t* buffer, const std::size_t stride)
	{
		std::vector<target_t> d(rows * cols);
		for (std::size_t i = 0; i < rows * cols; ++i)
			d[i] = static_cast<target_t>(buffer[i * stride]);

		(*_opened.top())[dataSetName] = d;
		(*_opened.top())[dataSetName + "_rows"] = rows;
		(*_opened.top())[dataSetName + "_cols"] = cols;
	}

	template <typename source_t, typename target_t>
	void writeTensor(const std::string& dataSetName, const std::size_t rank, const std::size_t* dims, const source_t* buffer, const std::size_t stride)
	{
		std::vector<std::size_t> dimsV(rank);
		std::size_t bufSize = 1;
		for (std::size_t i = 0; i < rank; ++i)
		{
			dimsV[i] = dims[i];
			bufSize *= dims[i];
		}

		std::vector<target_t> d(bufSize);
		for (std::size_t i = 0; i < bufSize; ++i)
			d[i] = static_cast<target_t>(buffer[i * stride]);

		(*_opened.top())[dataSetName] = d;
		(*_opened.top())[dataSetName + "_rank"] = rank;
		(*_opened.top())[dataSetName + "_dims"] = dimsV;
	}


private:
	nlohmann::json* _root;
	std::stack<nlohmann::json*> _opened;
};


} // namespace

namespace cadet
{
namespace io
{

template <typename io_t> using FileReaderImpl = BaseIOWrapper<FileReader<io_t>>;
template <typename io_t> using FileWriterImpl = BaseIOWrapper<FileWriter<io_t>>;

IMemoryIO* createMemoryReaderWriter()
{
	return new JSONMemoryReaderWriter();
}

IFileReader* createReader(const std::string& fileExt)
{
	if (cadet::util::caseInsensitiveEquals(fileExt, "h5"))
	{
		return new FileReaderImpl<cadet::io::HDF5Reader>();
	}
	else if (cadet::util::caseInsensitiveEquals(fileExt, "xml"))
	{
		return new FileReaderImpl<cadet::io::XMLReader>();
	}
	return nullptr;
}

IFileWriter* createWriter(const std::string& fileExt)
{
	if (cadet::util::caseInsensitiveEquals(fileExt, "h5"))
	{
		return new FileWriterImpl<cadet::io::HDF5Writer>();
	}
	else if (cadet::util::caseInsensitiveEquals(fileExt, "xml"))
	{
		return new FileWriterImpl<cadet::io::XMLWriter>();
	}
	return nullptr;
}

} // namespace io
} // namespace cadet

