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

#ifndef HDF5READER_HPP_
#define HDF5READER_HPP_

#include <vector>
#include <string>

#include <H5Cpp.h>

#include "HDF5Base.hpp"

namespace cadet
{

class HDF5Reader : public HDF5Base
{
public:
    /// \brief Constructor
    HDF5Reader();

    /// \brief Destructor
    ~HDF5Reader();

//    /// \brief Read data as a tensor from a dataset
//    template <typename T>
//    T tensor(const std::string& dataSetName) throw (H5::Exception);
//
//    /// \brief Convenience wrapper for reading matrices
//    template <typename T>
//    T matrix(const std::string& dataSetName);

    /// \brief Convenience wrapper for reading vectors
    template <typename T>
    std::vector<T> vector(const std::string& dataSetName);

    /// \brief Convenience wrapper for reading scalars
    template <typename T>
    T scalar(const std::string& dataSetName, size_t position = 0);

private:
    template <typename T>
    std::vector<T> read(const std::string& dataSetName);
};




// ====================================================================================================================
//   IMPLEMENTATION PART
// ====================================================================================================================


HDF5Reader::HDF5Reader() {}

HDF5Reader::~HDF5Reader() {}


// ============================================================================================================
//   Template specializations of member functions for diffenet data types
// ============================================================================================================
// Double specialization of vector()
template <>
std::vector<double> HDF5Reader::vector<double>(const std::string& dataSetName)
{
    return read<double>(dataSetName);
}

// Integer specialization of vector()
template <>
std::vector<int> HDF5Reader::vector<int>(const std::string& dataSetName)
{
    return read<int>(dataSetName);
}

// std::string specialization of vector()
template <>
std::vector<std::string> HDF5Reader::vector<std::string>(const std::string& dataSetName)
{
    // Get the dataset we want to read from
    openGroup();
    _dataSet = _groupsOpened.top().openDataSet(dataSetName);
    closeGroup();

    // Determine the datatype
    _dataType = _dataSet.getDataType();
    size_t bufSize = _dataSet.getSpace().getSimpleExtentNpoints();

    std::vector<std::string> stringVector;
    if (_dataSet.getDataType().isVariableStr())
    {
        char** buffer  = new char*[bufSize];

        // Read data from file and write it to buffer
        _dataSet.read(buffer, _dataType);

        // Copy read c-strings to a vector of std::strings
        for (size_t i = 0; i < bufSize; ++i)
            stringVector.push_back(std::string(buffer[i]));

        // Free memory alloc'd by the variable length read mechanism
        _dataSet.vlenReclaim((void*)buffer, _dataType, _dataSet.getSpace());
        delete [] buffer;
    }
    else
    {
        // dont know how to get the length of each string, if multiple strings were stored in the dataset...
        if (bufSize > 1) throw (H5::Exception(CURRENT_FUNCTION, "You can only read fixed length strings through scalar()"));

        char* buffer = new char[_dataType.getSize()];

        // Read data from file and write it to buffer
        _dataSet.read(buffer, _dataType);

        stringVector.push_back(std::string(buffer));
        delete [] buffer;
    }
    
    _dataType.close();
    _dataSet.close();

    return stringVector;
}

// Template that matches on every unsupported type and throws an exception
template <typename T>
std::vector<T> HDF5Reader::vector(const std::string& dataSetName)
{
    throw (H5::Exception(CURRENT_FUNCTION, "You may not try to read an unsupported type"));
}
// ============================================================================================================


template <typename T>
T HDF5Reader::scalar(const std::string& dataSetName, size_t position)
{
    return vector<T>(dataSetName).at(position);
}


template <typename T>
std::vector<T> HDF5Reader::read(const std::string& dataSetName)
{
    // Get the dataset we want to read from
    openGroup();
    _dataSet = _groupsOpened.top().openDataSet(dataSetName);
    closeGroup();

    size_t bufSize = _dataSet.getSpace().getSimpleExtentNpoints();
    T* buffer = new T[bufSize];

    // Determine the datatype
    _dataType = _dataSet.getDataType();

    // Read data from file and write it to buffer
    _dataSet.read(buffer, _dataType);
    std::vector<T> bufferVector(buffer, buffer + bufSize);
    delete [] buffer;

    _dataType.close();
    _dataSet.close();

    return bufferVector;
}


}  // namespace cadet


#endif /* HDF5READER_HPP_ */
