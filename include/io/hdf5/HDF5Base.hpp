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

#ifndef HDF5BASE_HPP_
#define HDF5BASE_HPP_

#include <vector>
#include <string>
#include <sstream>
#include <stack>

#include "cadet/cadetCompilerInfo.hpp"

#include "io/IOException.hpp"
#include <hdf5.h>

namespace cadet
{

namespace io
{

class HDF5Base
{
public:
	/// \brief Constructor
	HDF5Base();

	/// \brief Destructor
	~HDF5Base() CADET_NOEXCEPT;

	/// \brief Open an HDF5 file
	inline void openFile(const std::string& fileName, const std::string& mode = "r");
	inline void openFile(const char* fileName, const std::string& mode = "r") { openFile(std::string(fileName), mode); }

	/// \brief Close the currently opened file
	inline void closeFile();

	/// \brief Set a group to be [read from/written to] in all subsequent calls to [read/write] methods
	inline void setGroup(const std::string& groupName);

	/// \brief Open the subgroup with the given name
	inline void pushGroup(const std::string& groupName);
	/// \brief Close the currently open subgroup 
	inline void popGroup();

	/// \brief Checks if the given dataset or group exists in the file
	inline bool exists(const std::string& elementName);
	inline bool exists(const char* elementName) { return exists(std::string(elementName)); }

	/// \brief Checks if the given dataset is a vector (i.e., has more than one value)
	inline bool isVector(const std::string& elementName);
	inline bool isVector(const char* elementName) { return isVector(std::string(elementName)); }

protected:
	hid_t _file;

	std::stack<hid_t> _groupsOpened;
	std::vector<std::string> _groupNames;

	void openGroup(bool forceCreation = false);
	void closeGroup();
	std::string getFullGroupName();
};


HDF5Base::HDF5Base()
{
	H5Eset_auto(H5E_DEFAULT, NULL, NULL);
	_groupNames.push_back("/");
}


HDF5Base::~HDF5Base() CADET_NOEXCEPT { }


void HDF5Base::openFile(const std::string& fileName, const std::string& mode)
{
	if      (mode == "r" ) _file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT); // open in read mode
	else if (mode == "rw") _file = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);   // open in read / write mode
	else if (mode == "c" ) _file = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); // create new file
	else if (mode == "co") _file = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); // create / overwrite new file
	else throw IOException("Wrong file open mode");

	if (_file < 0)
		throw IOException("Failed to open or create HDF5 file \"" + fileName + "\" in mode " + mode);
}


void HDF5Base::closeFile()
{
	closeGroup();
	H5Fclose(_file);
}


bool HDF5Base::exists(const std::string& elementName)
{
	openGroup();

	// Try to open elementName as group
	hid_t grp = H5Gopen2(_groupsOpened.top(), elementName.c_str(), H5P_DEFAULT);
	if (grp >= 0)
	{
		// Found the group
		H5Gclose(grp);
		closeGroup();
		return true;        
	}

	// Try to open elementName as dataset
	hid_t dset = H5Dopen2(_groupsOpened.top(), elementName.c_str(), H5P_DEFAULT);
	if (dset >= 0)
	{
		// Found the dataset
		H5Dclose(dset);
		closeGroup();
		return true;
	}

	// Not found
	closeGroup();
	return false;
}


bool HDF5Base::isVector(const std::string& elementName)
{
	openGroup();

	hid_t dset = H5Dopen2(_groupsOpened.top(), elementName.c_str(), H5P_DEFAULT);
	if (dset < 0)
	{
		closeGroup();
		return false;
	}

	hid_t dspace = H5Dget_space(dset);
	if (dspace < 0)
	{
		H5Dclose(dset);
		closeGroup();
		return false;
	}

	const bool isVector = H5Sget_simple_extent_npoints(dspace) > 1;

	H5Sclose(dspace);
	H5Dclose(dset);
	closeGroup();

	return isVector;
}


void HDF5Base::setGroup(const std::string& groupName)
{
	_groupNames.clear();

	size_t start   = 0;
	size_t end     = 0;
	std::string delimiter("/");

	// Quick return when called with empty group name
	if (groupName.empty())
	{
		_groupNames.push_back("/");
		return;
	}

	// Don't care for a preceding delimiter
	if (groupName[0] == delimiter[0]) ++start;

	while (end != std::string::npos)
	{
		end = groupName.find(delimiter, start);

		// If at end, use length = maxLength.  Else use length = end - start.
		_groupNames.push_back(delimiter + groupName.substr(start, (end == std::string::npos) ? std::string::npos : end - start));

		// If at end, use start = maxSize.  Else use start = end + delimiter.
		start = ((end > (std::string::npos - delimiter.size())) ? std::string::npos : end + delimiter.size());
	}
}


void HDF5Base::pushGroup(const std::string& groupName)
{
	_groupNames.push_back("/" + groupName);
}


void HDF5Base::popGroup()
{
	_groupNames.pop_back();
}


void HDF5Base::openGroup(bool forceCreation)
{
	std::string dynName;
	for (std::vector<std::string>::const_iterator it = _groupNames.begin(); it < _groupNames.end(); ++it)
	{
		dynName += *it;
		hid_t grp = H5Gopen2(_file, dynName.c_str(), H5P_DEFAULT);
		if (grp >= 0)
			_groupsOpened.push(grp);
		else if (forceCreation)
			_groupsOpened.push(H5Gcreate2(_file, dynName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
		else
			throw IOException("Group '" + dynName + "' doesn't exist in file");
	}
}


void HDF5Base::closeGroup()
{
	while (_groupsOpened.size() >= 1)
	{
		H5Gclose(_groupsOpened.top());
		_groupsOpened.pop();
	}
}


std::string HDF5Base::getFullGroupName()
{
	std::ostringstream oss;
	for (std::vector<std::string>::const_iterator it = _groupNames.begin(); it < _groupNames.end(); ++it)
	{
		oss << *it;
	}
	return oss.str();
}

}  // namespace io

}  // namespace cadet


#endif /* HDF5BASE_HPP_ */
