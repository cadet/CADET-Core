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

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cstring>

#include <pugixml.hpp>

#include "CadetLogger.hpp"
#include "Simulator.hpp"
#include "hdf5/HDF5Writer.hpp"
#include "xml/XMLBase.hpp"

namespace cadet
{

using namespace pugi;

class XMLtoH5 : public xml_tree_walker , public XMLBase
{
public:
    XMLtoH5(const std::string& h5FileName);
    virtual ~XMLtoH5();

    virtual bool for_each(xml_node& node);

private:
    HDF5Writer _h5w;

    template <typename T>
    void writeToHDF5(const std::string& dataSetName, const size_t rank, const size_t* dims, const size_t bufSize, const std::vector<std::string>& data);

    std::ostringstream& insert_front(std::ostringstream& oss, const std::string& s) const;
};


XMLtoH5::XMLtoH5(const std::string& h5FileName)
{
    H5::Exception::dontPrint();

    // Create/Overwrite HDF5 file
    _h5w.openFile(h5FileName, "co");
    _h5w.compressFields(true);
    _h5w.extendibleFields(false);
}


XMLtoH5::~XMLtoH5()
{
    // Close HDF5 file
    _h5w.closeFile();
}


template <typename T>
void XMLtoH5::writeToHDF5(const std::string& dataSetName, const size_t rank, const size_t* dims, const size_t bufSize, const std::vector<std::string>& data)
{
    T* buffer = new T[bufSize];
    for (size_t i = 0; i < bufSize; ++i)
    {
        std::stringstream ss(data.at(i));
        ss >> buffer[i];
    }

    // Detect scalars
    if (rank == 0)
    {
        _h5w.scalar<T>(dataSetName, buffer[0]);
    }
    else
    {
        _h5w.write<T>(dataSetName, rank, dims, buffer);
    }

    delete [] buffer;
}


// Helper method to insert a string to the front of a output stringstream
std::ostringstream& XMLtoH5::insert_front(std::ostringstream& oss, const std::string& s) const
{
    std::streamsize pos = oss.tellp();
    oss.str(s + oss.str());
    oss.seekp(pos + s.length());
    return oss;
}



bool XMLtoH5::for_each(xml_node& node)
{
    if (strcmp(node.name(), _nodeDset.c_str()) == 0)
    {
        // Read dataset from xml and write to hdf5
        std::string data_str = node.text().get();
        size_t rank          = node.attribute(_attrRank.c_str()).as_int();
        std::string dims_str = node.attribute(_attrDims.c_str()).value();
        std::string type_str = node.attribute(_attrType.c_str()).value();

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

        // Split data and convert to right type
        std::vector<std::string> data_vec = split(data_str, _textSeparator.c_str());
        if (bufSize != data_vec.size())
        {
            std::ostringstream oss;
            oss << "XML file is inconsistent: Possibly wrong no. of entrys in dataset '" << node.attribute(_attrName.c_str()).value() << "'";
            throw CadetException(oss.str());
        }

        if (type_str.compare(_typeInt) == 0)
            writeToHDF5<int>(node.attribute(_attrName.c_str()).value(), rank, dims, bufSize, data_vec);
        else if (type_str.compare(_typeDouble) == 0)
            writeToHDF5<double>(node.attribute(_attrName.c_str()).value(), rank, dims, bufSize, data_vec);
        else if (type_str.compare(_typeChar) == 0)
            writeToHDF5<std::string>(node.attribute(_attrName.c_str()).value(), rank, dims, bufSize, data_vec);
        // Bool is not yet supported!
//        else if (type_str.compare(_typeBool) == 0)
//            writeToHDF5<bool>(node.attribute(_attrName.c_str()).value(), rank, dims, bufSize, data_vec);
        else
            throw CadetException("You may not try to convert an unsupported type");

        delete [] dims;
    }
    else if (strcmp(node.name(), _nodeGrp.c_str()) == 0)
    {
        // Build up the group path
        xml_node mynode(node);
        std::ostringstream groupPath;
        while (mynode != mynode.root())
        {
            insert_front(groupPath, "/" + std::string(mynode.attribute(_attrName.c_str()).value()));
            mynode = mynode.parent();
        }

        // Set the current group
        _h5w.setGroup(groupPath.str());
    }

    return true;
}

} // namespace cadet



void printHelp()
{
    using cadet::log;
    log::emit() << "Usage: xmltoh5 FILE.XML <FILE.H5>" << log::endl;
    log::emit() << "Convert an XML cadet input file to an HDF5 cadet input file." << log::endl;
    log::emit() << "To be used with cadet-cs from the Cromatography Analysis and DEsign Toolkit (CADET)." << log::endl;
    log::emit() << "When called with only input file specified, the output is named accordingly." << log::endl;
    log::emit() << "Examples: xmltoh5 intput.xml             <-- creates input.h5" << log::endl;
    log::emit() << "          xmltoh5 intput.xml output.h5   <-- creates output.h5" << log::endl;
    log::emit() << log::endl;
    log::emit() << "Report bugs to: cadet@fz-juelich.de" << log::endl;
    log::emit() << "CADET homepage: <http://www.cadet-web.de>" << log::endl;
    log::emit() << "Fork CADET on GitHub: <https://github.com/modsim/CADET>" << log::endl;
}


int main(int argc, char** argv)
{
    using namespace cadet;

    std::string outname;
    if (argc == 2)
    {
        outname.assign(argv[1]);
        std::size_t pos = outname.rfind(".xml");
        if (pos != std::string::npos) outname.erase(pos);
        outname.append(".h5");
    }
    if (argc == 3) outname.assign(argv[2]);
    if ((argc < 2) || (argc > 3))
    {
        printHelp();
        return 1;
    }

    using namespace pugi;

    try {
        // Open XML file
        xml_document doc;
        if (!doc.load_file(argv[1])) throw CadetException("XML file could not be opened!");

        XMLtoH5 converter(outname);
        doc.traverse(converter);
    }
    catch (const CadetException& e) {
        log::emit<Except>() << e.msg() << log::endl;
    }
    catch (const H5::Exception& e) {
        log::emit<Except>() << e.getDetailMsg() << log::endl;
    }

    return 0;
}
