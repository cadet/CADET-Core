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

#include <sstream>
#include <string>
#include <vector>
#include <cstring>

#include <H5Cpp.h>

#include "CadetLogger.hpp"
#include "Simulator.hpp"
#include "xml/XMLWriter.hpp"

namespace cadet
{

class H5toXML
{
public:
    H5toXML(const std::string& h5FileName, const std::string& xmlFileName);
    ~H5toXML();

    void convert();
    inline XMLWriter& getWriter() {return _xmlw;}

    static herr_t operate(hid_t obj, const char* name, const H5O_info_t* info, void* op_data);

private:

    hid_t _file;
    herr_t _error;
    XMLWriter _xmlw;
};

H5toXML::H5toXML(const std::string& h5FileName, const std::string& xmlFileName)
{
//    H5::Exception::dontPrint();
    _error = 0;

    // Open HDF5 file
    _file = H5Fopen(h5FileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (_file < 0) throw CadetException("HDF5 file could not be opened!");

    // Create/Overwrite XML file
    _xmlw.openFile(xmlFileName, "co");
}


H5toXML::~H5toXML()
{
    // Close XML file
    _xmlw.closeFile();

    // Close HDF5 file
    _error = H5Fclose(_file);
    if (_error) throw CadetException("HDF5 file could not be closed!");
}


void H5toXML::convert()
{
    _error = H5Ovisit(_file, H5_INDEX_NAME, H5_ITER_NATIVE, &cadet::H5toXML::operate, this);
}


herr_t H5toXML::operate(hid_t obj, const char* name, const H5O_info_t* info, void* op_data)
{
    // Quick return on root group
    if (strcmp(name, ".") == 0) return 0;

    H5toXML* h5toxml = static_cast<H5toXML*>(op_data);

    // Check what type of object we are dealing with
    if (info->type == H5O_TYPE_GROUP)
    {
        // Set the group in XMLWriter
        h5toxml->getWriter().setGroup(std::string(name));
    }
    else if (info->type == H5O_TYPE_DATASET)
    {
        // Extract name of the dataset
        const char* lastSlash = strrchr(name, '/');
        std::string datasetName(lastSlash+1);

        hid_t dataset = H5Dopen(obj, name, H5P_DEFAULT);
        hid_t dataspace = H5Dget_space(dataset);

        // Determine the class of the data
        hid_t datatype = H5Dget_type(dataset);
        H5T_class_t dataclass = H5Tget_class(datatype);

        // Get rank and number of elements
        size_t rank = H5Sget_simple_extent_ndims(dataspace);
        hssize_t bufSize = H5Sget_simple_extent_npoints(dataspace);

        hsize_t* dims = new hsize_t[rank];
        H5Sget_simple_extent_dims(dataspace, dims, NULL);

        if (dataclass == H5T_INTEGER)
        {
            int* buffer = new int[bufSize];
            H5Dread(dataset, H5T_STD_I32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
            h5toxml->getWriter().write<int>(datasetName, rank, (const size_t*) dims, buffer);
            delete [] buffer;
        }
        else if (dataclass == H5T_STRING)
        {
            hid_t memtype = H5Tcopy(H5T_C_S1);
            H5Tset_size(memtype, H5T_VARIABLE);

            char** buffer = new char*[bufSize];
            H5Dread(dataset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);

            // Copy read c-strings to a vector of std::strings
            std::vector<std::string> stringVector;
            for (ssize_t i = 0; i < bufSize; ++i)
                stringVector.push_back(std::string(buffer[i]));

            h5toxml->getWriter().write<std::string>(datasetName, rank, (const size_t*) dims, &stringVector[0]);

            // Free memory alloc'd by the variable length read mechanism
            H5Dvlen_reclaim(memtype, H5Dget_space(dataset), H5P_DEFAULT, buffer);
            delete [] buffer;

            H5Tclose(memtype);
        }
        else if (dataclass == H5T_FLOAT)
        {
            double* buffer = new double[bufSize];
            H5Dread(dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
            h5toxml->getWriter().write<double>(datasetName, rank, (const size_t*) dims, buffer);
            delete [] buffer;
        }
        else
            throw CadetException("Unsupported type in input HDF5 file!");

        delete [] dims;

        H5Tclose(datatype);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }
    else
        throw CadetException("HDF5 file includes unknown node types. Only groups/datasets can be converted!");

    return 0;
}


} // namespace cadet

void printHelp()
{
    using cadet::log;
    log::emit() << "Usage: h5toxml FILE.H5 <FILE.XML>" << log::endl;
    log::emit() << "Convert an HDF5 cadet input file to an XML cadet input file." << log::endl;
    log::emit() << "To be used with cadet-cs from the Cromatography Analysis and DEsign Toolkit (CADET)." << log::endl;
    log::emit() << "When called with only input file specified, the output is named accordingly." << log::endl;
    log::emit() << "Examples: h5toxml intput.h5              <-- creates input.xml" << log::endl;
    log::emit() << "          h5toxml intput.h5 output.xml   <-- creates output.xml" << log::endl;
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
        std::size_t pos = outname.rfind(".h5");
        if (pos != std::string::npos) outname.erase(pos);
        outname.append(".xml");
    }
    if (argc == 3) outname.assign(argv[2]);
    if ((argc < 2) || (argc > 3))
    {
        printHelp();
        return 1;
    }

    try {
        H5toXML h5toxml(argv[1], outname);
        h5toxml.convert();
    }
    catch (const CadetException& e){
        log::emit<Except>() << e.msg() << log::endl;
    }
    catch (const H5::Exception& e){
        log::emit<Except>() << e.getDetailMsg() << log::endl;
    }

    return 0;
}
