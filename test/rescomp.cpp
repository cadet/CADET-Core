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
#include <limits>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <algorithm>
#include <functional>

#include "CadetEnumeration.hpp"

#include <H5Cpp.h>
#include <tclap/CmdLine.h>

struct DatasetDiff
{
	std::string name;
	
	double diffL1;
	double diffL2;
	double diffLinf;

	double refL1;
	double refL2;
	double refLinf;
};

struct Statistics
{
	int elemMissing;
	int datasets;
	int numDiffs;
	int structMismatches;
	int typeMismatches;
	int checked;

	std::vector<DatasetDiff> diffSummary;
};

struct Element
{
	std::string name;
	bool isGroup;
};

struct Config
{
	double absL1;
	double absL2;
	double absLinf;

	double relL1;
	double relL2;
	double relLinf;

	bool verboseSuccess;
	bool summary;
};

const std::vector<Element> getElementsInGroup(const H5::Group& group)
{
    hsize_t items = group.getNumObjs();
    std::vector<Element> elms;

    for (hsize_t i = 0; i < items; ++i)
    {
        elms.push_back(Element { group.getObjnameByIdx(i), group.getObjTypeByIdx(i) == H5G_GROUP });
    }

    return elms;
}

void checkDataSets(const H5::DataSet& ds1, const H5::DataSet& ds2, const std::string& name, const std::string& dsName, Statistics& stat, const Config& cfg)
{
	// Check rank
	const int nDim1 = ds1.getSpace().getSimpleExtentNdims();
	const int nDim2 = ds2.getSpace().getSimpleExtentNdims();

	if (nDim1 != nDim2)
	{
		std::cout << "[DIFF] Dataset \"" << name << "\" has different rank (reference = " << nDim1 << ", trial = " << nDim2 << ")" << std::endl;
		++stat.structMismatches;
		return;
	}

	// Check dimensions
	hsize_t* dims1 = new hsize_t[nDim1];
	hsize_t* dims2 = new hsize_t[nDim1];

	ds1.getSpace().getSimpleExtentDims(dims1);
	ds2.getSpace().getSimpleExtentDims(dims2);

	for (int i = 0; i < nDim1; ++i)
	{
		if (dims1[i] != dims2[i])
		{
			std::cout << "[DIFF] Dataset \"" << name << "\" has different dimensionality (dim = " << i << ", reference = " << dims1[i] << ", trial = " << dims2[i] << ")" << std::endl;
			++stat.structMismatches;

			delete[] dims2;
			delete[] dims1;
			return;
		}
	}

	delete[] dims2;
	delete[] dims1;

	// Check type
	H5::DataType dt1 = ds1.getDataType();
	H5::DataType dt2 = ds2.getDataType();

	if (dt1.getClass() != dt2.getClass())
	{
		std::cout << "[DIFF] Type of dataset \"" << name << "\" has different classes (reference = " << dt1.getClass() << ", trial = " << dt2.getClass() << ")" << std::endl;
		++stat.typeMismatches;
	
		dt2.close();
		dt1.close();
		return;
	}

	if (dt1.getSize() != dt2.getSize())
	{
		std::cout << "[DIFF] Type of dataset \"" << name << "\" has different sizes (reference = " << dt1.getSize() << ", trial = " << dt2.getSize() << ")" << std::endl;
		++stat.typeMismatches;
	
		dt2.close();
		dt1.close();
		return;
	}

	if (!(dt1 == dt2))
	{
		std::cout << "[DIFF] Type of dataset \"" << name << "\" differs" << std::endl;
		++stat.typeMismatches;
	
		dt2.close();
		dt1.close();
		return;
	}

	size_t dataSize = dt1.getSize();
	bool isDouble = (dt1 == H5::PredType::NATIVE_DOUBLE);

	dt2.close();
	dt1.close();

	// Check elements

	if (isDouble)
	{
		size_t bufSize = ds1.getSpace().getSimpleExtentNpoints();
		double* buffer1 = new double[bufSize];
		double* buffer2 = new double[bufSize];

		// Read data from file and write it to buffer
		ds1.read(buffer1, H5::PredType::NATIVE_DOUBLE);
		ds2.read(buffer2, H5::PredType::NATIVE_DOUBLE);
		
		// Compute norms
		double diffNormL1 = 0.0;
		double diffNormL2 = 0.0;
		double diffNormLinf = 0.0;
		double normL1 = 0.0;
		double normL2 = 0.0;
		double normLinf = 0.0;

		for (size_t i = 0; i < bufSize; ++i)
		{
			const double diff = buffer1[i] - buffer2[i];
			diffNormL1 += std::abs(diff);
			diffNormL2 += diff * diff;
			diffNormLinf = std::max(std::abs(diff), diffNormLinf);

			normL1 += std::abs(buffer1[i]);
			normL2 += buffer1[i] * buffer1[i]; 
			normLinf = std::max(std::abs(buffer1[i]), normLinf);
		}

		delete [] buffer2;
		delete [] buffer1;

		normL2 = std::sqrt(normL2);
		diffNormL2 = std::sqrt(diffNormL2);

		stat.diffSummary.push_back({ name, diffNormL1, diffNormL2, diffNormLinf, normL1, normL2, normLinf });

		// Compare to thresholds
		bool diff = false;
		if ( (diffNormL1 > cfg.absL1) || (diffNormL1 > cfg.relL1 * normL1) )
		{
			diff = true;

			std::cout << "[NUM ] Dataset \"" << name << "\" differs ";
			if (diffNormL1 > cfg.absL1)
			{
				std::cout << "L1 " << diffNormL1 << " > Tol " << cfg.absL1 << " ";
			}
			if (diffNormL1 > cfg.relL1 * normL1)
			{
				std::cout << "relL1 " << (diffNormL1 / normL1) << " > Tol " << cfg.relL1 << " ";
			}
		}

		if ( (diffNormL2 > cfg.absL2) || (diffNormL2 > cfg.relL2 * normL2) )
		{
			if (!diff)
				std::cout << "[NUM ] Dataset \"" << name << "\" differs ";
			else
				std::cout << "; ";

			diff = true;
			if (diffNormL2 > cfg.absL2)
			{
				std::cout << "L2 " << diffNormL2 << " > Tol " << cfg.absL2 << " ";
			}
			if (diffNormL2 > cfg.relL2 * normL2)
			{
				std::cout << "relL2 " << (diffNormL2 / normL2) << " > Tol " << cfg.relL2 << " ";
			}
		}

		if ( (diffNormLinf > cfg.absLinf) || (diffNormLinf > cfg.relLinf * normLinf) )
		{
			if (!diff)
				std::cout << "[NUM ] Dataset \"" << name << "\" differs ";
			else
				std::cout << "; ";

			diff = true;
			if (diffNormLinf > cfg.absLinf)
			{
				std::cout << "Linf " << diffNormLinf << " > Tol " << cfg.absLinf << " ";
			}
			if (diffNormLinf > cfg.relLinf * normLinf)
			{
				std::cout << "relLinf " << (diffNormLinf / normLinf) << " > Tol " << cfg.relLinf << " ";
			}
		}

		if (!diff && cfg.verboseSuccess)
		{
			std::cout << "[ OK ] Dataset \"" << name << "\" matches" << std::endl;
		}
		else if (diff)
		{
			++stat.numDiffs;
			std::cout << std::endl;
		}
	}
	else
	{
		size_t bufSize = ds1.getSpace().getSimpleExtentNpoints() * dataSize;
		char* buffer1 = new char[bufSize];
		char* buffer2 = new char[bufSize];

		// Read data from file and write it to buffer
		ds1.read(buffer1, H5::PredType::NATIVE_CHAR);
		ds2.read(buffer2, H5::PredType::NATIVE_CHAR);

		bool diff = false;
		for (size_t i = 0; i < bufSize; ++i)
		{
			if (buffer1[i] != buffer2[i])
			{
				diff = true;
				break;
			}
		}

		delete [] buffer2;
		delete [] buffer1;

		if (!diff && cfg.verboseSuccess)
		{
			std::cout << "[ OK ] Dataset \"" << name << "\" matches" << std::endl;
		}
		else if (diff)
		{
			std::cout << "[NUM ] Dataset \"" << name << "\" differs" << std::endl;
			++stat.numDiffs;
		}
	}

	++stat.checked;
}

void recursiveCheck(const H5::Group& group1, const H5::Group& group2, const std::string& prefix, Statistics& stat, const Config& cfg)
{
    const std::vector<Element> elms1 = getElementsInGroup(group1);
    const std::vector<Element> elms2 = getElementsInGroup(group2);
	
	for (std::vector<Element>::const_iterator it = elms1.begin(); it != elms1.end(); ++it)
	{

		// Check existence in reader2
		std::vector<Element>::const_iterator it2 = std::find_if(elms2.begin(), elms2.end(), [&it](const Element& v) { return it->name == v.name; });
		if (it2 == elms2.end())
		{
			std::cout << "[MISS] Field \"" << prefix << "/" << it->name << "\" in reference but not in trial file" << std::endl;
			++stat.elemMissing;
			continue;
		}

		// Check types
		if (it2->isGroup != it->isGroup)
		{
			std::cout << "[MISS] Field \"" << prefix << "/" << it->name << "\" in reference and trial file are not of the same element type (group, dataset, etc.)" << std::endl;
			++stat.structMismatches;
			continue;
		}

		if (it->isGroup)
		{
			// Dive into this group

			H5::Group gLocal1 = group1.openGroup(it->name);
			H5::Group gLocal2 = group2.openGroup(it->name);
			
			recursiveCheck(gLocal1, gLocal2, prefix + "/" + it->name, stat, cfg);

			gLocal2.close();
			gLocal1.close();
		}
		else
		{
			// It is a dataset
			H5::DataSet ds1 = group1.openDataSet(it->name);
			H5::DataSet ds2 = group2.openDataSet(it->name);

			++stat.datasets;
			checkDataSets(ds1, ds2, prefix + "/" + it->name, it->name, stat, cfg);

			ds2.close();
			ds1.close();
		}
	}
}

int main(int argc, char** argv)
{
    Statistics stat = { 0, 0, 0, 0, 0, 0, std::vector<DatasetDiff>() };
//    Config cfg = { 1e-8, 1e-8, 1e-8, 1e-6, 1e-6, 1e-6, true };
    Config cfg;

    std::string refFile;
    std::string trialFile;

	try
	{
		TCLAP::CmdLine cmd("Compares two HDF5 files (reference and trial) and checks for differing vectors", ' ', "v0.1");

		cmd >> (new TCLAP::ValueArg<double>("", "absL1", "Absolute L1 tolerance", false, 1e-8, "TOL"))->storeIn(&cfg.absL1);
		cmd >> (new TCLAP::ValueArg<double>("", "absL2", "Absolute L2 tolerance", false, 1e-8, "TOL"))->storeIn(&cfg.absL2);
		cmd >> (new TCLAP::ValueArg<double>("", "absLinf", "Absolute Linf tolerance", false, 1e-8, "TOL"))->storeIn(&cfg.absLinf);
		cmd >> (new TCLAP::ValueArg<double>("", "relL1", "Relative L1 tolerance", false, 1e-6, "TOL"))->storeIn(&cfg.relL1);
		cmd >> (new TCLAP::ValueArg<double>("", "relL2", "Relative L2 tolerance", false, 1e-6, "TOL"))->storeIn(&cfg.relL2);
		cmd >> (new TCLAP::ValueArg<double>("", "relLinf", "Relative Linf tolerance", false, 1e-6, "TOL"))->storeIn(&cfg.relLinf);
		cmd >> (new TCLAP::SwitchArg("", "vs", "Outputs a message upon successful matching of two datasets"))->storeIn(&cfg.verboseSuccess);
		cmd >> (new TCLAP::SwitchArg("s", "summary", "Outputs a table with results at the end"))->storeIn(&cfg.summary);
		cmd >> (new TCLAP::UnlabeledValueArg<std::string>("ref", "Reference HDF5 file", true, "", "REF-HDF5"))->storeIn(&refFile);
		cmd >> (new TCLAP::UnlabeledValueArg<std::string>("trial", "Trial HDF5 file", true, "", "TRIAL-HDF5"))->storeIn(&trialFile);

		cmd.parse( argc, argv );
	}
	catch (const TCLAP::ArgException &e)
	{
		std::cerr << "Error: " << e.error() << " for argument " << e.argId() << std::endl;
		return 1;
	}

    std::cout << "===============================================" << std::endl;
    std::cout << "Ref: " << refFile << std::endl;
    std::cout << "Trial: " << trialFile << std::endl;
    std::cout << "Tol L1: Abs " << cfg.absL1 << " Rel " << cfg.relL1 << std::endl;
    std::cout << "Tol L2: Abs " << cfg.absL2 << " Rel " << cfg.relL2 << std::endl;
    std::cout << "Tol Linf: Abs " << cfg.absLinf << " Rel " << cfg.relLinf << std::endl;
    std::cout << "===============================================" << std::endl;

	H5::H5File file1(refFile, H5F_ACC_RDONLY);
	H5::H5File file2(trialFile, H5F_ACC_RDONLY);

	H5::Group g1 = file1.openGroup(cadet::e2s(cadet::GRP_OUT));
	H5::Group g2 = file2.openGroup(cadet::e2s(cadet::GRP_OUT));

    recursiveCheck(g1, g2, cadet::e2s(cadet::GRP_OUT), stat, cfg);

    g2.close();
    g1.close();

    file2.close();
    file1.close();

    std::cout << "===============================================" << std::endl;
    std::cout << "Missing elements: " << stat.elemMissing << std::endl;
    std::cout << "Datasets: " << stat.datasets << std::endl;
    std::cout << "   -> Structural mismatches: " << stat.structMismatches << std::endl;
    std::cout << "   -> Type mismatches: " << stat.typeMismatches << std::endl;
    std::cout << "   -> Checked: " << stat.checked << std::endl;
    std::cout << "   -> Matches: " << (stat.checked - stat.numDiffs) << std::endl;
    std::cout << "   -> Diffs: " << stat.numDiffs << std::endl;

  	std::cout << "===============================================" << std::endl;

	if (cfg.summary)
	{
		const char separator    = ' ';
		unsigned int nameWidth  = 6;
		const int numWidth      = 12;

		// Compute longest length of name
		for (std::vector<DatasetDiff>::iterator it = stat.diffSummary.begin(); it != stat.diffSummary.end(); ++it)
		{
			nameWidth = std::max(nameWidth, static_cast<unsigned int>(it->name.size()));
		}

		// Table headings
		std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) << "Dataset" << " | ";
		std::cout << std::setw(numWidth) << std::setfill(separator) << "L1" << " | ";
		std::cout << std::setw(numWidth) << std::setfill(separator) << "L2" << " | ";
		std::cout << std::setw(numWidth) << std::setfill(separator) << "Linf" << " | ";
		std::cout << std::setw(numWidth) << std::setfill(separator) << "Rel L1" << " | ";
		std::cout << std::setw(numWidth) << std::setfill(separator) << "Rel L2" << " | ";
		std::cout << std::setw(numWidth) << std::setfill(separator) << "Rel Linf" << " | ";
		std::cout << std::setw(numWidth) << std::setfill(separator) << "Ref L1" << " | ";
		std::cout << std::setw(numWidth) << std::setfill(separator) << "Ref L2" << " | ";
		std::cout << std::setw(numWidth) << std::setfill(separator) << "Ref Linf" << std::endl;

		for (unsigned int i = 0; i < nameWidth + (numWidth + 3) * 9; ++i)
		{
			std::cout << "-";
		}
		std::cout << std::endl;

		// Table body
		for (std::vector<DatasetDiff>::iterator it = stat.diffSummary.begin(); it != stat.diffSummary.end(); ++it)
		{
			std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) << it->name << " | ";
			std::cout << std::setw(numWidth) << std::setfill(separator) << it->diffL1 << " | ";
			std::cout << std::setw(numWidth) << std::setfill(separator) << it->diffL2 << " | ";
			std::cout << std::setw(numWidth) << std::setfill(separator) << it->diffLinf << " | ";
			std::cout << std::setw(numWidth) << std::setfill(separator) << it->diffL1 / it->refL1 << " | ";
			std::cout << std::setw(numWidth) << std::setfill(separator) << it->diffL2 / it->refL2 << " | ";
			std::cout << std::setw(numWidth) << std::setfill(separator) << it->diffLinf / it->refLinf << " | ";
			std::cout << std::setw(numWidth) << std::setfill(separator) << it->refL1 << " | ";
			std::cout << std::setw(numWidth) << std::setfill(separator) << it->refL2 << " | ";
			std::cout << std::setw(numWidth) << std::setfill(separator) << it->refLinf << std::endl;
		}
	}

    if ((stat.elemMissing == 0) && (stat.checked == stat.datasets) && (stat.numDiffs == 0))
    {
    	std::cout << "[ OK ] Files match" << std::endl;
    	return 0;
    }
    else
    {
    	std::cout << "[FAIL] Files differ" << std::endl;
    	return 1;
    }
}

