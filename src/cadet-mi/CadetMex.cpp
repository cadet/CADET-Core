// =============================================================================
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

#ifndef MATLAB_MEX_FILE
    #define MATLAB_MEX_FILE
#endif

#ifdef _MSC_VER
	#define DLL_EXPORT_SYM _declspec(dllexport)
#else
    #define DLL_EXPORT_SYM
#endif

#include <mex.h>

#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <iostream>

// Comment out to enable logging
#define LOGGING_DISABLE
#include "CadetLogger.hpp"

// Take care of namespace pollution / macros
#ifdef min
    #undef min
#endif
#ifdef max
    #undef max
#endif

#include "driver/CadetCSDriver.hpp"
#include "Cadet.hpp"

class MatlabCadetTranslator
{
public:
    /// \brief Constructor
    MatlabCadetTranslator(mxArray** data) : _root(data), _groupName(""), _group(0), _dataSet(0)
    {
    }
    MatlabCadetTranslator() : _root(0), _groupName(""), _group(0), _dataSet(0)
    {
    }

    MatlabCadetTranslator(const MatlabCadetTranslator& cpy) : _root(cpy._root), _groupName(""), _group(cpy._group), _dataSet(cpy._dataSet)
    {
    }

    /// \brief Destructor
    ~MatlabCadetTranslator() { }

    /// \brief Open an HDF5 file
    inline void openFile(const std::string& fileName, const std::string& mode = "r") { }
    inline void openFile(const char* fileName, const std::string& mode = "r") { }

    /// \brief Close the currently opened file
    inline void closeFile() { }

    /// \brief Set a group to be [read from/written to] in all subsequent calls to [read/write] methods
    inline void setGroup(const std::string& groupName)
    {
        _groupName = groupName;
    }

    /// \brief Checks if the given dataset or group exists in the file
    inline bool exists(const std::string& elementName) { return exists(elementName.c_str()); }
    inline bool exists(const char* elementName);

    /// \brief Checks if the given dataset is a vector (i.e., has more than one value)
    inline bool isVector(const std::string& elementName) { return isVector(elementName.c_str()); }
    inline bool isVector(const char* elementName);

    /// \brief Convenience wrapper for reading vectors
    template <typename T>
    std::vector<T> vector(const std::string& dataSetName);

    /// \brief Convenience wrapper for reading scalars
    template <typename T>
    T scalar(const std::string& dataSetName, size_t position = 0);

    /// \brief Write data from C-array to a dataset
    template <typename T>
    void write(const std::string& dataSetName, const size_t rank, const size_t* dims, const T* buffer);

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

    /// \brief Enable/disable compression for tensors of 2nd order and above
    inline void compressFields(bool setCompression) { }

    /// \brief Tensors of 2nd order (vectors) and above are written as extendible fields
    ///        (maxsize = unlimited, chunked layout), when set to true.
    inline void extendibleFields(bool setExtendible) { }
protected:
    std::string _groupName;
    mxArray**   _root;
    mxArray*    _group;
    mxArray*    _dataSet;

    template <typename T>
    std::vector<T> read();

    mxArray* writeWork(const std::string& dataSetName, const size_t rank, const size_t* dims, const void* buffer, const mxClassID type);

    void openGroup(bool create = false);
    void openAndCreateGroup();

    bool checkForVector() const
    {
        // Validate the current item (expect column or row vector)
        const mwSize numDims = mxGetNumberOfDimensions(_dataSet);
        if (numDims > 2)
        {
            // Error
            std::stringstream str;
            str << "CadetMex: Trying to read vector from field '" << _groupName << "', but got tensor rank " << numDims << ".\n";
            mexErrMsgTxt(str.str().c_str());
            return false;
        }

        const mwSize* dims = mxGetDimensions(_dataSet);
        if ((dims[0] > 1) && (dims[1] > 1))
        {
            // Error
            std::stringstream str;
            str << "CadetMex: Trying to read vector from field '" << _groupName << "', but got " << dims[0] << " x " << dims[1] << " matrix.\n";
            mexErrMsgTxt(str.str().c_str());
            return false;
        }
        return true;        
    }

};






DLL_EXPORT_SYM void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 0)
    {
        // Print version
        if (nlhs == 0)
            mexPrintf("This is CADET version %s built from commit %s\n", cadet::getLibraryVersion(), cadet::getLibraryCommitHash());

        // Return version
        if (nlhs >= 1)
        {
            plhs[0] = mxCreateString(cadet::getLibraryVersion());
        }
        if (nlhs == 2)
        {
            plhs[1] = mxCreateString(cadet::getLibraryCommitHash());
        }
        return;
    }

    if ((nrhs != 1) || (nlhs != 1))
        mexErrMsgIdAndTxt("CADET:mexError", "Expecting exactly one input and one output.\n");

    if (!mxIsStruct(prhs[0]))
        mexErrMsgTxt("CadetMex: Expecting input argument of type 'struct'.\n");

    // IMPORTANT: Reset NumDir variable of ADOL-C
    cadetResetGlobals();

    mxArray* temp = const_cast<mxArray*>(prhs[0]);

    MatlabCadetTranslator input(&temp);
    MatlabCadetTranslator output(&plhs[0]);

    try
    {
        cadet::CadetCS<MatlabCadetTranslator, MatlabCadetTranslator> cs("", input, output);
        cs.run();            
    }
    catch(const cadet::CadetException& e)
    {
        mexErrMsgIdAndTxt("CADET:mexError", "Error in simulation: %s\n", e.msg().c_str());
    }
    catch(...)
    {
        mexErrMsgIdAndTxt("CADET:mexError", "Error in simulation: Unkown error.\n");
    }
}






bool MatlabCadetTranslator::exists(const char* elementName)
{
    openGroup();
    return mxGetField(_group, 0, elementName);
}

bool MatlabCadetTranslator::isVector(const char* elementName)
{
    openGroup();
    return mxGetNumberOfElements(mxGetField(_group, 0, elementName)) > 1;
}

void MatlabCadetTranslator::openGroup(bool create)
{
    _group = *_root;
    if (create)
        return openAndCreateGroup();

    // Early out if empty group
    if (_groupName.empty())
        return;

    // Ignore root slash if present
    std::size_t start = 0;
    if (_groupName[0] == '/')
        start = 1;

    std::size_t end = 0;
    std::string currentPath;
    while (end != std::string::npos)
    {
        end = _groupName.find("/", start);

        // If at end, use length = maxLength.  Else use length = end - start.
        std::string itemName = _groupName.substr(start, (end == std::string::npos) ? std::string::npos : end - start);

        // If at end, use start = maxSize.  Else use start = end + delimiter.
        start = ((end > (std::string::npos - 1)) ? std::string::npos : end + 1);

        // Check if current item is a struct
        if (!mxIsStruct(_group))
        {
            // Error
            std::stringstream str;
            str << "CadetMex: The element '" << currentPath << "' is not a struct (requested field '" << _groupName << "').\n";
            mexErrMsgTxt(str.str().c_str());
        }

        _group = mxGetField(_group, 0, itemName.c_str());

        // Check if field exists
        if (!_group)
        {
            // Error
            std::stringstream str;
            str << "CadetMex: The struct '" << currentPath << "' does not contain the field '" << itemName << "' (requested field '" << _groupName << "').\n";
            mexErrMsgTxt(str.str().c_str());
        }

        if (!currentPath.empty())
            currentPath += "." + itemName;
        else
            currentPath += itemName;
    }
}

void MatlabCadetTranslator::openAndCreateGroup()
{
    _group = *_root;

    // Check if root is a struct
    if (!_group || !mxIsStruct(_group))
    {
        _group = mxCreateStructMatrix(1, 1, 0, NULL);
        *_root = _group;
    }

    // Early out if empty group
    if (_groupName.empty())
        return;

    // Ignore root slash if present
    std::size_t start = 0;
    if (_groupName[0] == '/')
        start = 1;

    std::size_t end = 0;
    std::string currentPath;
    while (end != std::string::npos)
    {
        end = _groupName.find("/", start);

        // If at end, use length = maxLength.  Else use length = end - start.
        std::string itemName = _groupName.substr(start, (end == std::string::npos) ? std::string::npos : end - start);

        // If at end, use start = maxSize.  Else use start = end + delimiter.
        start = ((end > (std::string::npos - 1)) ? std::string::npos : end + 1);

        mxArray* nextGroup = mxGetField(_group, 0, itemName.c_str());

        // Check if field exists
        if (!nextGroup)
        {
            // Create the field and insert a struct
            mwIndex index = mxAddField(_group, itemName.c_str());
            mxSetFieldByNumber(_group, 0, index, mxCreateStructMatrix(1, 1, 0, NULL));
            nextGroup = mxGetFieldByNumber(_group, 0, index);
        }

        _group = nextGroup;

        if (!currentPath.empty())
            currentPath += "." + itemName;
        else
            currentPath += itemName;
    }
}
// ============================================================================================================
// Reading
// ============================================================================================================


// Double specialization of vector()
template <>
std::vector<double> MatlabCadetTranslator::vector<double>(const std::string& dataSetName)
{
    openGroup();

    // Validate data type
    _dataSet = mxGetField(_group, 0, dataSetName.c_str());
    if (!_dataSet)
    {
        // Error
        std::stringstream str;
        str << "CadetMex: Trying to read vector from non-existent field '" << _groupName << "." << dataSetName << "'.\n";
        mexErrMsgTxt(str.str().c_str());
        return std::vector<double>();
    }
    if (!mxIsDouble(_dataSet))
    {
        // Error
        std::stringstream str;
        str << "CadetMex: Trying to read vector of doubles from field '" << _groupName << "." << dataSetName << "', but got non-double values.\n";
        mexErrMsgTxt(str.str().c_str());
        return std::vector<double>();
    }

    return read<double>();
}

// Integer specialization of vector()
template <>
std::vector<int> MatlabCadetTranslator::vector<int>(const std::string& dataSetName)
{
    openGroup();

    // Validate data type
    _dataSet = mxGetField(_group, 0, dataSetName.c_str());
    if (!_dataSet)
    {
        // Error
        std::stringstream str;
        str << "CadetMex: Trying to read vector from non-existent field '" << _groupName << "." << dataSetName << "'.\n";
        mexErrMsgTxt(str.str().c_str());
        return std::vector<int>();
    }
    if (!mxIsInt32(_dataSet))
    {
        // Error
        std::stringstream str;
        str << "CadetMex: Trying to read vector of int32 from field '" << _groupName << "." << dataSetName << "', but got non-int32 values.\n";
        mexErrMsgTxt(str.str().c_str());
        return std::vector<int>();
    }

    return read<int>();
}

// std::string specialization of vector()
template <>
std::vector<std::string> MatlabCadetTranslator::vector<std::string>(const std::string& dataSetName)
{
    openGroup();

    // Validate
    _dataSet = mxGetField(_group, 0, dataSetName.c_str());
    if (!_dataSet)
    {
        // Error
        std::stringstream str;
        str << "CadetMex: Trying to read cell array from non-existent field '" << _groupName << "." << dataSetName << "'.\n";
        mexErrMsgTxt(str.str().c_str());
        return std::vector<std::string>();
    }
    if (!mxIsCell(_dataSet) && !mxIsChar(_dataSet))
    {
        // Error
        std::stringstream str;
        str << "CadetMex: Trying to read cell array of strings from field '" << _groupName << "." << dataSetName << "', but got non-cell array.\n";
        mexErrMsgTxt(str.str().c_str());
        return std::vector<std::string>();
    }

    // Check for single string
    if (mxIsChar(_dataSet))
    {
        std::vector<std::string> data(1);
        
        const char* const strData = mxArrayToString(_dataSet);
        data[0] = strData;
        mxFree(const_cast<char*>(strData));

        return data;
    }

    // We have a cell array
    const mwSize numElems = mxGetNumberOfElements(_dataSet);
    std::vector<std::string> stringVector(numElems);

    for (mwIndex i = 0; i < numElems; ++i)
    {
        mxArray* curCell = mxGetCell(_dataSet, i);
        if (!mxIsChar(_dataSet))
        {
            // Error
            std::stringstream str;
            str << "CadetMex: Expected string in element " << (i+1) << " of cell array " << _groupName << "." << dataSetName << "', but got non-string data.\n";
            mexErrMsgTxt(str.str().c_str());
            return std::vector<std::string>();
        }

        const char* const strData = mxArrayToString(curCell);
        stringVector[i] = strData;
        mxFree(const_cast<char*>(strData));
    }

    return stringVector;
}

// Template that matches on every unsupported type
template <typename T>
std::vector<T> MatlabCadetTranslator::vector(const std::string& dataSetName)
{
    std::stringstream str;
    str << "CadetMex: Trying to read vector of unsupported type from field '" << _groupName << "." << dataSetName << "'.\n";
    mexErrMsgTxt(str.str().c_str());
    return std::vector<T>();
}
// ============================================================================================================



template <typename T>
T MatlabCadetTranslator::scalar(const std::string& dataSetName, size_t position)
{
    return this->template vector<T>(dataSetName).at(position);
}


template <typename T>
std::vector<T> MatlabCadetTranslator::read()
{
    checkForVector();

    const mwSize numElems = mxGetNumberOfElements(_dataSet);

    // Copy data
    if (mxIsDouble(_dataSet))
    {
        double const* const src = mxGetPr(_dataSet);
        return std::vector<T>(src, src + numElems);
    }
    else if (mxIsInt32(_dataSet))
    {
        int const* const src = static_cast<int const* const>(mxGetData(_dataSet));
        return std::vector<T>(src, src + numElems);
    }
    else
    {
        // Unsupported data type
        std::stringstream str;
        str << "CadetMex: Trying to read vector of unsupported type from struct '" << _groupName << "'.\n";
        mexErrMsgTxt(str.str().c_str());
        return std::vector<T>();
    }
}


// ============================================================================================================
// Writing
// ============================================================================================================

mwIndex rowMajorToColMajor(mwIndex idx, const size_t rank, const size_t* fwdProd)
{
    std::vector<size_t> matIdx(rank, 0);

    // Phase 1: Calculate index in every rank from row-major index
    for (size_t i = 0; i < rank; ++i)
    {
        matIdx[i] = idx % fwdProd[rank - i - 1];
        idx -= matIdx[i] * fwdProd[rank - i - 1];
    }

    // Phase 2: Calculate col-major index from linear index
    mwIndex out = 0;
    for (size_t i = 0; i < rank; ++i)
    {
        out += matIdx[0] * fwdProd[i];
    }
    return out;
}

void rowMajorToIndex(mwIndex idx, const size_t rank, const size_t* fwdProd, size_t* matIdx)
{
    for (size_t i = 0; i < rank; ++i)
    {
        matIdx[i] = idx % fwdProd[rank - i - 1];
        idx -= matIdx[i] * fwdProd[rank - i - 1];
    }
}

mxArray* MatlabCadetTranslator::writeWork(const std::string& dataSetName, const size_t rank, const size_t* dims, const void* buffer, const mxClassID type)
{
    openAndCreateGroup();

    // Create Matlab array
    mwSize* dimsMatlab = static_cast<mwSize*>(mxMalloc(rank * sizeof(mwIndex)));
    for (size_t i = 0; i < rank; ++i)
    {
        dimsMatlab[i] = dims[i];
    }

    mxArray* data = mxCreateNumericArray(rank, dimsMatlab, type, mxREAL);

    mxFree(dimsMatlab);

    // Assign it to the field
    mxSetFieldByNumber(_group, 0, mxAddField(_group, dataSetName.c_str()), data);
    
    return data;
}

template <>
void MatlabCadetTranslator::write<double>(const std::string& dataSetName, const size_t rank, const size_t* dims, const double* buffer)
{
    mxArray* matData = writeWork(dataSetName, rank, dims, buffer, mxDOUBLE_CLASS);

    // Allocate memory for array
    size_t numEl = mxGetNumberOfElements(matData);
    mxSetData(matData, mxMalloc(sizeof(double) * numEl));
    double* const data = mxGetPr(matData);

    std::vector<size_t> fwdProd(rank, 1);
    for (size_t i = 1; i < rank; ++i)
    {
        fwdProd[i] = fwdProd[i-1] * dims[i];
    }

    // Write data to matlab
    for (size_t i = 0; i < numEl; ++i)
    {
//        data[rowMajorToColMajor(i, rank, &fwdProd[0])] = buffer[i];
        data[i] = buffer[i];
    }
}

template <>
void MatlabCadetTranslator::write<int>(const std::string& dataSetName, const size_t rank, const size_t* dims, const int* buffer)
{
    mxArray* matData = writeWork(dataSetName, rank, dims, buffer, mxINT32_CLASS);

    // Allocate memory for array
    size_t numEl = mxGetNumberOfElements(matData);
    mxSetData(matData, mxMalloc(sizeof(int) * numEl));
    int* const data = static_cast<int*>(mxGetData(matData));

    std::vector<size_t> fwdProd(rank, 1);
    for (size_t i = 1; i < rank; ++i)
    {
        fwdProd[i] = fwdProd[i-1] * dims[i];
    }

    // Write data to matlab
    for (size_t i = 0; i < numEl; ++i)
    {
//        data[rowMajorToColMajor(i, rank, &fwdProd[0])] = buffer[i];
        data[i] = buffer[i];
    }
}

template <>
void MatlabCadetTranslator::write<std::string>(const std::string& dataSetName, const size_t rank, const size_t* dims, const std::string* buffer)
{
    openAndCreateGroup();

    // Create Matlab array
    mwSize* dimsMatlab = static_cast<mwSize*>(mxMalloc(rank * sizeof(mwIndex)));
    for (size_t i = 0; i < rank; ++i)
        dimsMatlab[i] = dims[i];

    mxArray* cell_array_ptr = mxCreateCellArray(rank, dimsMatlab);
    mxFree(dimsMatlab);

    size_t numEl = mxGetNumberOfElements(cell_array_ptr);
    for (mwIndex i = 0; i < numEl; ++i)
    {
        mxSetCell(cell_array_ptr, i, mxCreateString(buffer[i].c_str()));
    }
    mxSetFieldByNumber(_group, 0, mxAddField(_group, dataSetName.c_str()), cell_array_ptr);
}

// Template that matches on every unsupported type
template <typename T>
void MatlabCadetTranslator::write(const std::string& dataSetName, const size_t rank, const size_t* dims, const T* buffer)
{
    std::stringstream str;
    str << "CadetMex: Trying to write vector of unsupported type to struct '" << _groupName << "." << dataSetName << "'.\n";
    mexErrMsgTxt(str.str().c_str());
}
// ============================================================================================================


// ============================================================================================================
//   Convenience wrappers
// ============================================================================================================

// HDF5: Row Major
// Matlab: Column Major 

template <typename T>
void MatlabCadetTranslator::tensor(const std::string& dataSetName, const size_t rank, const size_t* dims, const std::vector<T>& buffer)
{
    // Check size
    size_t bufSize = 1;
    for (size_t i = 0; i < rank; ++i)
        bufSize *= dims[i];
    
    if (bufSize != buffer.size())
    {
        std::stringstream str;
        str << "CadetMex: Trying to write tensor of size ";
        for (size_t i = 0; i < rank - 1; ++i)
            str << dims[i] << " x ";
        str << dims[rank-1] << " but got " << buffer.size() << " elements.\n";

        mexErrMsgTxt(str.str().c_str());
        return;
    }
    write<T>(dataSetName, rank, dims, &buffer[0]);
}

template <typename T>
void MatlabCadetTranslator::tensor(const std::string& dataSetName, const size_t rank, const size_t* dims, const T* buffer)
{
    write<T>(dataSetName, rank, dims, buffer);
}

template <typename T>
void MatlabCadetTranslator::matrix(const std::string& dataSetName, const size_t rows, const size_t cols, const T* buffer)
{
    size_t dims[2] = {rows, cols};
    write<T>(dataSetName, 2, dims, buffer);
}

template <typename T>
void MatlabCadetTranslator::matrix(const std::string& dataSetName, const size_t rows, const size_t cols, const std::vector<T>& buffer)
{
    if (rows*cols != buffer.size())
    {
        std::stringstream str;
        str << "CadetMex: Trying to write matrix of size " << rows << "x" << cols << " but got " << buffer.size() << " elements.\n";
        mexErrMsgTxt(str.str().c_str());
        return;
    }
    size_t dims[2] = {rows, cols};
    write<T>(dataSetName, 2, dims, &buffer[0]);
}

template <typename T>
void MatlabCadetTranslator::vector(const std::string& dataSetName, const size_t length, const T* buffer)
{
    write<T>(dataSetName, 1, &length, buffer);
}

template <typename T>
void MatlabCadetTranslator::vector(const std::string& dataSetName, const std::vector<T>& buffer)
{
    size_t length = buffer.size();
    write<T>(dataSetName, 1, &length, &buffer[0]);
}

template <typename T>
void MatlabCadetTranslator::scalar(const std::string& dataSetName, const T buffer)
{
    vector<T>(dataSetName, 1, &buffer);
}
// ============================================================================================================


void MatlabCadetTranslator::unlinkGroup(const std::string& groupName)
{
}
