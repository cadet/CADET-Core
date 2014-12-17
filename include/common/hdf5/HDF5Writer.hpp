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

#ifndef HDF5WRITER_HPP_
#define HDF5WRITER_HPP_

#include <vector>
#include <string>
#include <cstring>

#include <H5Cpp.h>

#include "HDF5Base.hpp"

namespace cadet
{

class HDF5Writer : public HDF5Base
{
public:
    /// \brief Constructor
    HDF5Writer();

    /// \brief Destructor
    ~HDF5Writer();

    /// \brief Write data from C-array to a dataset
    template <typename T>
    void write(const std::string& dataSetName, const size_t rank, const size_t* dims, const T* buffer) throw (H5::Exception);

    /// \brief Convenience wrapper for writing tensors from C-array
    template <typename T>
    void tensor(const std::string& dataSetName, const size_t rank, const size_t* dims, const T* buffer);

    /// \brief Convenience wrapper for writing tensors from std::vector
    template <typename T>
    void tensor(const std::string& dataSetName, const size_t rank, const size_t* dims, const std::vector<T>& buffer);

    /// \brief Convenience wrapper for writing matrices from C-array
    template <typename T>
    void matrix(const std::string& dataSetName, const size_t rows, const size_t cols, const T* buffer);

    /// \brief Convenience wrapper for writing matrices from std::vector
    template <typename T>
    void matrix(const std::string& dataSetName, const size_t rows, const size_t cols, const std::vector<T>& buffer);

    /// \brief Convenience wrapper for writing vectors from C-array
    template <typename T>
    void vector(const std::string& dataSetName, const size_t length, const T* buffer);

    /// \brief Convenience wrapper for writing vectors from std::vector
    template <typename T>
    void vector(const std::string& dataSetName, const std::vector<T>& buffer);

    /// \brief Convenience wrapper for writing scalars
    template <typename T>
    void scalar(const std::string& dataSetName, const T buffer);

    /// \brief Removes an existing group from the file
    inline void unlinkGroup(const std::string& groupName);

    /// \brief Removes an existing dataset from the current group
    inline void unlinkDataset(const std::string& dsName);

    /// \brief Enable/disable compression for tensors of 2nd order and above
    inline void compressFields(bool setCompression) {_writeCompressed = setCompression;}

    /// \brief Tensors of 2nd order (vectors) and above are written as extendible fields
    ///        (maxsize = unlimited, chunked layout), when set to true.
    inline void extendibleFields(bool setExtendible) {_writeExtendible = setExtendible;}

private:

    void writeWork(const std::string& dataSetName, const size_t rank, const size_t* dims, const void* buffer);

    bool                    _writeScalar;
    bool                    _writeExtendible;
    bool                    _writeCompressed;
    H5::DSetCreatPropList   _dsetCPL;
    hsize_t*                _maxDims;
    hsize_t*                _chunks;
    double                  _chunkFactor;
};




// ====================================================================================================================
//   IMPLEMENTATION PART
// ====================================================================================================================

HDF5Writer::HDF5Writer() :
        _writeScalar(false),
        _writeExtendible(true),
        _writeCompressed(false),
        _dsetCPL(H5::DSetCreatPropList::DEFAULT),
        _maxDims(NULL),
        _chunks(NULL),
        _chunkFactor(1.5)
{}

HDF5Writer::~HDF5Writer() {}


// ============================================================================================================
//   Template specializations of member function write() for diffenet data types
// ============================================================================================================
template <>
void HDF5Writer::write<double>(const std::string& dataSetName, const size_t rank, const size_t* dims, const double* buffer) throw (H5::Exception)
{
    _dataType = H5::PredType::IEEE_F64LE;
    writeWork(dataSetName, rank, dims, buffer);
}

template <>
void HDF5Writer::write<int>(const std::string& dataSetName, const size_t rank, const size_t* dims, const int* buffer) throw (H5::Exception)
{
    _dataType = H5::PredType::STD_I32LE;
    writeWork(dataSetName, rank, dims, buffer);
}

template <>
void HDF5Writer::write<std::string>(const std::string& dataSetName, const size_t rank, const size_t* dims, const std::string* buffer) throw (H5::Exception)
{
    _dataType = H5::StrType(0, H5T_VARIABLE);

    size_t bufSize = 1;
    for (size_t i = 0; i < rank; ++i)
    {
        bufSize *= dims[i];
    }

    // Copy strings to char*
    char** strBuffer = new char*[bufSize];
    for (size_t i = 0; i < bufSize; ++i)
    {
        strBuffer[i] = new char[buffer[i].size() + 1]; // +1 for null-terminating character
        strcpy(strBuffer[i], buffer[i].c_str());
    }

    writeWork(dataSetName, rank, dims, strBuffer);

    // Memory cleanup
    for (size_t i = 0; i < bufSize; ++i)
        delete [] strBuffer[i];
    delete [] strBuffer;
}

// Template that matches on every unsupported type and throws an exception
template <typename T>
void HDF5Writer::write(const std::string& dataSetName, const size_t rank, const size_t* dims, const T* buffer) throw (H5::Exception)
{
    throw (H5::Exception(CURRENT_FUNCTION, "You may not try to write an unsupported type"));
}
// ============================================================================================================


// ============================================================================================================
//   Convenience wrappers
// ============================================================================================================
template <typename T>
void HDF5Writer::tensor(const std::string& dataSetName, const size_t rank, const size_t* dims, const std::vector<T>& buffer)
{
    size_t bufSize = 1;
    for (size_t i = 0; i < rank; ++i)
        bufSize *= dims[i];
    if (bufSize != buffer.size())
        throw (H5::Exception(CURRENT_FUNCTION, "Tensor dimensions must suite size of given vector"));
    write<T>(dataSetName, rank, dims, &buffer[0]);
}

template <typename T>
void HDF5Writer::tensor(const std::string& dataSetName, const size_t rank, const size_t* dims, const T* buffer)
{
    write<T>(dataSetName, rank, dims, buffer);
}

template <typename T>
void HDF5Writer::matrix(const std::string& dataSetName, const size_t rows, const size_t cols, const T* buffer)
{
    size_t dims[2] = {rows, cols};
    write<T>(dataSetName, 2, dims, buffer);
}

template <typename T>
void HDF5Writer::matrix(const std::string& dataSetName, const size_t rows, const size_t cols, const std::vector<T>& buffer)
{
    if (rows*cols not_eq buffer.size())
        throw (H5::Exception(CURRENT_FUNCTION, "Matrix dimensions must suite size of given vector"));
    size_t dims[2] = {rows, cols};
    write<T>(dataSetName, 2, dims, &buffer[0]);
}

template <typename T>
void HDF5Writer::vector(const std::string& dataSetName, const size_t length, const T* buffer)
{
    write<T>(dataSetName, 1, &length, buffer);
}

template <typename T>
void HDF5Writer::vector(const std::string& dataSetName, const std::vector<T>& buffer)
{
    size_t length = buffer.size();
    write<T>(dataSetName, 1, &length, &buffer[0]);
}

template <typename T>
void HDF5Writer::scalar(const std::string& dataSetName, const T buffer)
{
    _writeScalar = true;
    vector<T>(dataSetName, 1, &buffer);
}
// ============================================================================================================


void HDF5Writer::unlinkGroup(const std::string& groupName)
{
    try {
        _file.unlink(groupName);
    } catch (H5::Exception& e)
    {
        return;
    }
}


void HDF5Writer::unlinkDataset(const std::string& dsName)
{
    bool wasOpen = !_groupsOpened.empty();

    if (!wasOpen)
        openGroup(false);

    try {
        _groupsOpened.top().unlink(dsName.c_str());
    } catch (H5::Exception& e)
    {
    }

    if (!wasOpen)
        closeGroup();
}


void HDF5Writer::writeWork(const std::string& dataSetName, const size_t rank, const size_t* dims, const void* buffer)
{
    // init CreatePropertyList with no chunking and no compression
    _dsetCPL = H5::DSetCreatPropList::DEFAULT;

    if (!_writeScalar)
    {
        if (_writeExtendible || _writeCompressed) // we need chunking
        {
            _chunks  = new hsize_t[rank];
            for (size_t i = 0; i < rank; ++i)
                _chunks[i] = (_writeExtendible) ? static_cast<hsize_t>(dims[i] * _chunkFactor) : dims[i]; // leave some space in all dims, if extendible

            _dsetCPL.setChunk(rank, _chunks);
            delete [] _chunks;
        }

        _maxDims = new hsize_t[rank];
        if (_writeExtendible) // we set maxdims unlimited
        {
            for (size_t i = 0; i < rank; ++i)
                _maxDims[i] = H5S_UNLIMITED;
        }
        else // we set maxdims to dims of the buffer
        {
            for (size_t i = 0; i < rank; ++i)
                _maxDims[i] = dims[i];
        }

        hsize_t* convDims = new hsize_t[rank];
        for (size_t i = 0; i < rank; ++i)
        {
            convDims[i] = dims[i];
        }

        _dataSpace = H5::DataSpace(rank, convDims, _maxDims);
        delete[] convDims;
        delete[] _maxDims;

        if (_writeCompressed) // enable compression
            _dsetCPL.setDeflate(9);
    }
    else // reset _writeScalar
    {
        // init dataspace as scalar
        _dataSpace = H5::DataSpace();
        _writeScalar = false;
    }

    bool forceCreation = true;
    openGroup(forceCreation);
    _dataSet = _groupsOpened.top().createDataSet(dataSetName, _dataType, _dataSpace, _dsetCPL);
    closeGroup();

    _dataSet.write(buffer, _dataType);
}


}  // namespace cadet

#endif /* HDF5WRITER_HPP_ */
