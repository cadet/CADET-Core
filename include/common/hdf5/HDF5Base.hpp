// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2015: Eric von Lieres¹, Joel Andersson,
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

#ifndef HDF5BASE_HPP_
#define HDF5BASE_HPP_

// Set appropriate compiler macros for windows or linux otherwise
#ifndef CURRENT_FUNCTION
  #ifdef _WIN32
    #define CURRENT_FUNCTION (__FUNCTION__)
  #else
    #define CURRENT_FUNCTION (__FUNCTION__)
  //#define CURRENT_FUNCTION (__PRETTY_FUNCTION__) // GCC: Function names including their type signature of the function as well as its bare name
  #endif
#endif

#include <vector>
#include <string>
#include <stack>

#include <H5Cpp.h>

namespace cadet
{

class HDF5Base
{
public:
    /// \brief Constructor
    HDF5Base();

    /// \brief Destructor
    ~HDF5Base();

    /// \brief Open an HDF5 file
    inline void openFile(const std::string& fileName, const std::string& mode = "r");
    inline void openFile(const char* fileName, const std::string& mode = "r") { openFile(std::string(fileName), mode); }

    /// \brief Close the currently opened file
    inline void closeFile();

    /// \brief Set a group to be [read from/written to] in all subsequent calls to [read/write] methods
    inline void setGroup(const std::string& groupName);

    /// \brief Checks if the given dataset or group exists in the file
    inline bool exists(const std::string& elementName);
    inline bool exists(const char* elementName) { return exists(std::string(elementName)); }

    /// \brief Checks if the given dataset is a vector (i.e., has more than one value)
    inline bool isVector(const std::string& elementName);
    inline bool isVector(const char* elementName) { return isVector(std::string(elementName)); }

protected:
    H5::H5File      _file;
    H5::DataSpace   _dataSpace;
    H5::DataSet     _dataSet;
    H5::DataType    _dataType;

    std::stack<H5::Group>    _groupsOpened;
    std::vector<std::string> _groupNames;

    void openGroup(bool forceCreation = false);
    void closeGroup();
};




// ====================================================================================================================
//    IMPLEMENTATION PART
// ====================================================================================================================
HDF5Base::HDF5Base()
{
    H5::Exception::dontPrint();
    _groupNames.push_back("/");
}


HDF5Base::~HDF5Base(){}


void HDF5Base::openFile(const std::string& fileName, const std::string& mode)
{
    if      (mode == "r" ) _file = H5::H5File(fileName, H5F_ACC_RDONLY);    // open in read mode
    else if (mode == "rw") _file = H5::H5File(fileName, H5F_ACC_RDWR);      // open in read / write mode
    else if (mode == "c" ) _file = H5::H5File(fileName, H5F_ACC_CREAT);     // create new file
    else if (mode == "co") _file = H5::H5File(fileName, H5F_ACC_TRUNC);     // create / overwrite new file
    else if (mode == "cf") _file = H5::H5File(fileName, H5F_ACC_EXCL);      // create / fail on existent file
    else if (mode == "d" ) _file = H5::H5File(fileName, H5F_ACC_DEBUG);     // print debug infos
    else throw H5::Exception(CURRENT_FUNCTION, "Wrong file open mode");
}


void HDF5Base::closeFile()
{
    closeGroup();
    _file.close();
}


bool HDF5Base::exists(const std::string& elementName)
{
    openGroup();

    // Try to open elementName as group
    try
    {
        H5::Group grp = _groupsOpened.top().openGroup(elementName);
        grp.close();

        // Found the group
        closeGroup();
        return true;
    }
    catch ( const H5::Exception& ) { }

    // Try to open elementName as dataset
    try
    {
        H5::DataSet ds = _groupsOpened.top().openDataSet(elementName);
        ds.close();

        // Found the dataset
        closeGroup();
        return true;
    }
    catch ( const H5::Exception& ) { }

    // Not found
    closeGroup();
    return false;
}


bool HDF5Base::isVector(const std::string& elementName)
{
    openGroup();
    bool isVector = false;
    try
    {
        H5::DataSet ds = _groupsOpened.top().openDataSet(elementName);
        isVector = ds.getSpace().getSimpleExtentNpoints() > 1;
        ds.close();
    }
    catch ( const H5::Exception& ) { }

    // Not found
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
    if (groupName.at(0) == delimiter.at(0)) ++start;

    while (end != std::string::npos)
    {
        end = groupName.find(delimiter, start);

        // If at end, use length = maxLength.  Else use length = end - start.
        _groupNames.push_back(delimiter + groupName.substr(start, (end == std::string::npos) ? std::string::npos : end - start));

        // If at end, use start = maxSize.  Else use start = end + delimiter.
        start = ((end > (std::string::npos - delimiter.size())) ? std::string::npos : end + delimiter.size());
    }
}


void HDF5Base::openGroup(bool forceCreation)
{
    std::string dynName;
    for (std::vector<std::string>::const_iterator it = _groupNames.begin(); it < _groupNames.end(); ++it)
    {
        dynName += *it;
        try {
            _groupsOpened.push(_file.openGroup(dynName));
        } catch (H5::Exception& e) {
            if (forceCreation)
                _groupsOpened.push(_file.createGroup(dynName));
            else
                throw (H5::GroupIException(CURRENT_FUNCTION, "Group '" + dynName + "' doesn't exist in file"));
        }
    }
}


void HDF5Base::closeGroup()
{
    while (_groupsOpened.size() >= 1)
    {
        _groupsOpened.top().close();
        _groupsOpened.pop();
    }
}


}  // namespace cadet


#endif /* HDF5BASE_HPP_ */
