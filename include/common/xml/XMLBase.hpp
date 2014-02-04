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

#ifndef XMLBASE_HPP_
#define XMLBASE_HPP_

#include <string>
#include <sstream>
#include <vector>
#include <stack>

#include "pugixml.hpp"

namespace cadet {

using namespace pugi;

class XMLBase
{
public:
    /// \brief Constructor
    XMLBase();

    /// \brief Destructor
    ~XMLBase();

    /// \brief Open an XML file
    inline void openFile(const std::string& fileName, const std::string& mode = "r");
    inline void openFile(const char* fileName, const std::string& mode = "r") { openFile(std::string(fileName), mode); }

    /// \brief Close the currently opened file
    inline void closeFile();

    /// \brief Set a group to be [read from/written to] in all subsequent calls to [read/write] methods
    inline void setGroup(const std::string& groupName);

    /// \brief Checks if the given dataset or group exists in the file
    inline bool exists(const std::string& elementName) { return exists(elementName.c_str()); }
    inline bool exists(const char* elementName);

    /// \brief Checks if the given dataset is a vector (i.e., has more than one value)
    inline bool isVector(const std::string& elementName) { return isVector(elementName.c_str()); }
    inline bool isVector(const char* elementName);

protected:

    xml_document _doc;                          //!< Root of the XML document
    bool _enable_write;                         //!< Specifies write permission for currently opened XML file
    std::string _fileName;                      //!< Name of the currently opened XML file

    static const std::string _textSeparator;    //!< Character sequence for separation of text entries
    static const std::string _dimsSeparator;    //!< Character sequence for separation of dimensions

    static const std::string _nodeGrp;          //!< Name used for creation of 'group' nodes
    static const std::string _nodeDset;         //!< Name used for creation of 'dataset' nodes

    static const std::string _attrName;         //!< Name used for creation of 'name' attributes
    static const std::string _attrRank;         //!< Name used for creation of 'rank' attributes
    static const std::string _attrDims;         //!< Name used for creation of 'dims' attributes
    static const std::string _attrType;         //!< Name used for creation of 'type' attributes
    static const std::string _attrValue;        //!< Name used for creation of 'value' attributes

    static const std::string _typeChar;         //!< Name used for 'type' attributes of type 'char'
    static const std::string _typeInt;          //!< Name used for 'type' attributes of type 'int'
    static const std::string _typeDouble;       //!< Name used for 'type' attributes of type 'double'
    static const std::string _typeBool;         //!< Name used for 'type' attributes of type 'bool'

    xpath_node               _groupOpened;      //!< Holds the group that is currently opened
    std::stack<std::string>  _groupExists;      //!< Stack holding the last successfully opened group
    std::vector<std::string> _groupNames;       //!< Vector of group names for the currently selected group

    void openGroup(bool forceCreation = false);
    void closeGroup();

    // Some helper methods
    std::vector<std::string>& split(const std::string& s, const char* delim, std::vector<std::string>& elems);
    std::vector<std::string> split(const std::string& s, const char* delim);
};


// Setting constant values
const std::string XMLBase::_textSeparator = ", ";
const std::string XMLBase::_dimsSeparator = "x";

const std::string XMLBase::_nodeGrp   = "group";
const std::string XMLBase::_nodeDset  = "dataset";

const std::string XMLBase::_attrName  = "name";
const std::string XMLBase::_attrRank  = "rank";
const std::string XMLBase::_attrDims  = "dims";
const std::string XMLBase::_attrType  = "type";
const std::string XMLBase::_attrValue = "value";

const std::string XMLBase::_typeChar   = "char";
const std::string XMLBase::_typeInt    = "int";
const std::string XMLBase::_typeDouble = "double";
const std::string XMLBase::_typeBool   = "bool";




// ====================================================================================================================
//    IMPLEMENTATION PART
// ====================================================================================================================
XMLBase::XMLBase() :
    _enable_write(false)
{}

XMLBase::~XMLBase()
{
    closeFile();
}



void XMLBase::openFile(const std::string& fileName, const std::string& mode)
{
    xml_parse_result flag;

    if      (mode == "r" ) // open in read mode
    {
        if (!_doc.load_file(fileName.c_str())) throw CadetException("XML file does not exist!");
        _enable_write = false;
    }
    else if (mode == "rw") // open in read / write mode
    {
        if (!_doc.load_file(fileName.c_str())) throw CadetException("XML file does not exist!");
        _enable_write = true;
    }
    else if (mode == "c" ) // create new file
    {
        if (_doc.load_file(fileName.c_str())) throw CadetException("XML file already exists");
        _enable_write = true;
    }
    else if (mode == "co") // create / overwrite new file
        _enable_write = true;
    else
        throw CadetException("Wrong file open mode");

    _fileName = fileName;
}




void XMLBase::closeFile()
{
    if (_enable_write) _doc.save_file(_fileName.c_str());
}


bool XMLBase::exists(const char* elementName)
{
    openGroup();
    bool exists = _groupOpened.node().child(elementName);
    closeGroup();
    return exists;
}


bool XMLBase::isVector(const char* elementName)
{
    openGroup();

    // Open dataset and throw if it does not exist
    xml_node dataset = _groupOpened.node().find_child_by_attribute(_nodeDset.c_str(), _attrName.c_str(), elementName);
    if (!dataset)
    {
        std::ostringstream oss;
        oss << "Dataset '" << elementName << "' does not exist!";
        throw CadetException(oss.str());
    }

    // Read text and attributes
    size_t rank          = dataset.attribute(_attrRank.c_str()).as_int();
    std::string dims_str = dataset.attribute(_attrDims.c_str()).value();

    // Get dims and compute buffer size
    std::vector<std::string> dims_vec = split(dims_str, _dimsSeparator.c_str());
    int items = 1;
    for (size_t i = 0; i < rank; ++i)
    {
        int j;
        std::stringstream ss(dims_vec.at(i));
        ss >> j;
        items *= j;
    }

    closeGroup();
    return items > 1;
}


void XMLBase::setGroup(const std::string& groupName)
{
    _groupNames.clear();

    size_t start   = 0;
    size_t end     = 0;
    std::string delimiter("/");

    // Quick return when called with empty group name
    if (groupName.empty())
        return;

    // Don't care for a preceding delimiter
    if (groupName.at(0) == delimiter.at(0)) ++start;

    while (end != std::string::npos)
    {
        end = groupName.find(delimiter, start);

        // If at end, use length = maxLength.  Else use length = end - start.
        _groupNames.push_back(groupName.substr(start, (end == std::string::npos) ? std::string::npos : end - start));

        // If at end, use start = maxSize.  Else use start = end + delimiter.
        start = ((end > (std::string::npos - delimiter.size())) ? std::string::npos : end + delimiter.size());
    }
}




void XMLBase::openGroup(bool forceCreation)
{
    // Clear existent group cache
    while (_groupExists.size() > 0)  _groupExists.pop();

    // Generate xpath query string
    std::ostringstream query;
    for (std::vector<std::string>::const_iterator it = _groupNames.begin(); it < _groupNames.end(); ++it)
    {
        query << "/" << _nodeGrp << "[@" << _attrName << "='" << *it << "']";
        // Store query string of parent groups for potential later usage
        if (_doc.select_single_node(query.str().c_str())) _groupExists.push(query.str());
    }

    // Try to open the selected group
    _groupOpened = _doc.select_single_node(query.str().c_str());

    // Create new group if not existent, creation is forced and write is permitted
    if (!_groupOpened)
    {
        if (forceCreation && _enable_write)
        {
            xml_node newNode;
            size_t ngrpexist = _groupExists.size();
            // Check for any parent group to be existent
            for (std::vector<std::string>::const_iterator it = _groupNames.begin() + ngrpexist; it < _groupNames.end(); ++it)
            {
                query.str("");
                // Select the least group
                if (_groupExists.size() > 0)
                {
                    _groupOpened = _doc.select_single_node(_groupExists.top().c_str());
                    // Create a new group
                    newNode = _groupOpened.node().append_child(_nodeGrp.c_str());
                    query << _groupExists.top();
                }
                else newNode = _doc.append_child(_nodeGrp.c_str());

                // Set attribute name
                newNode.append_attribute(_attrName.c_str()) = it->c_str();

                // Store the new group in existent groups
                query << "/" << _nodeGrp << "[@" << _attrName << "='" << *it << "']";
                _groupExists.push(query.str());
            }
            // Open newly created group
            _groupOpened = _doc.select_single_node(query.str().c_str());
        }
        else throw CadetException("Group was not opened/created! Either not existent, creation not forced or file not opened in write mode");
    }
}



void XMLBase::closeGroup()
{
    while (_groupExists.size() > 0)  _groupExists.pop();
}



std::vector<std::string>& XMLBase::split(const std::string& s, const char* delim, std::vector<std::string>& elems)
{
    std::stringstream ss(s);
    std::string item;

    while(std::getline(ss, item, *delim))
        elems.push_back(item);

    return elems;
}

std::vector<std::string> XMLBase::split(const std::string& s, const char* delim)
{
    std::vector<std::string> elems;
    return split(s, delim, elems);
}


} // namespace cadet


#endif /* XMLBASE_HPP_ */
