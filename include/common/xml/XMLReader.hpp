// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2016: Eric von Lieres¹, Joel Andersson¹,
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

#ifndef XMLREADER_HPP_
#define XMLREADER_HPP_

#include "XMLBase.hpp"

namespace cadet {

using namespace pugi;

class XMLReader : public XMLBase
{
public:
    /// \brief Constructor
    XMLReader();

    /// \brief Destructor
    ~XMLReader();

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

XMLReader::XMLReader() {}

XMLReader::~XMLReader() {}



// ============================================================================================================
//   Template specializations of member functions for diffenet data types
// ============================================================================================================
// Double specialization of vector()
template <>
std::vector<double> XMLReader::vector<double>(const std::string& dataSetName)
{
    return read<double>(dataSetName);
}

// Integer specialization of vector()
template <>
std::vector<int> XMLReader::vector<int>(const std::string& dataSetName)
{
    return read<int>(dataSetName);
}

// std::string specialization of vector()
template <>
std::vector<std::string> XMLReader::vector<std::string>(const std::string& dataSetName)
{
    return read<std::string>(dataSetName);
}

// Template that matches on every unsupported type and throws an exception
template <typename T>
std::vector<T> XMLReader::vector(const std::string& dataSetName)
{
    throw CadetException("You may not try to read an unsupported type");
}
// ============================================================================================================



template <typename T>
T XMLReader::scalar(const std::string& dataSetName, size_t position)
{
    return vector<T>(dataSetName).at(position);
}


template <typename T>
std::vector<T> XMLReader::read(const std::string& dataSetName)
{
    openGroup();

    // Open dataset and throw if it does not exist
    xml_node dataset = _groupOpened.node().find_child_by_attribute(_nodeDset.c_str(), _attrName.c_str(), dataSetName.c_str());
    if (!dataset)
    {
        std::ostringstream oss;
        oss << "Dataset '" << dataSetName << "' does not exist!";
        throw CadetException(oss.str());
    }

    // Read text and attributes
    std::string data_str = dataset.text().get();
    size_t rank          = dataset.attribute(_attrRank.c_str()).as_int();
    std::string dims_str = dataset.attribute(_attrDims.c_str()).value();

    // Get dims and compute buffer size
    std::vector<std::string> dims_vec = split(dims_str, _dimsSeparator.c_str());
    size_t* dims = new size_t[rank];
    size_t bufSize = 1;
    for (size_t i = 0; i < rank; ++i)
    {
        std::stringstream ss(dims_vec.at(i));
        ss >> dims[i];
        bufSize *= dims[i];
    }
    delete [] dims;

    // Split data and convert to right type
    std::vector<std::string> data_vec = split(data_str, _textSeparator.c_str());
    if (bufSize != data_vec.size())
    {
        std::ostringstream oss;
        oss << "XML file is inconsistent: Possibly wrong no. of entrys in dataset '" << dataset.attribute(_attrName.c_str()).value() << "'";
        throw CadetException(oss.str());
    }

    std::vector<T> data(bufSize);
    for (size_t i = 0; i < bufSize; ++i)
    {
        std::stringstream ss(data_vec.at(i));
        ss >> data.at(i);
    }

    closeGroup();
    return data;
}


} // namespace cadet


#endif /* XMLREADER_HPP_ */
