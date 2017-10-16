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

#include <sstream>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>
#include <random>
#include <chrono>
#include <limits>
#include <cmath>

#include "cadet/cadet.hpp"
#include "io/FileIO.hpp"
#include "io/IOException.hpp"

#include <tclap/CmdLine.h>
#include "common/TclapUtils.hpp"
#include "ToolsHelper.hpp"
#include "FormatConverter.hpp"


const int formatVersionLowest = 20000;
const int formatVersionRecent = 30100;


template <class Writer_t>
class FileCloser
{
public:
	FileCloser(Writer_t& writer) : _writer(writer) { }
	~FileCloser() { _writer.closeFile(); }
protected:
	Writer_t& _writer;
};

template <class class_t>
class DeleterScope
{
public:
	DeleterScope(class_t* elem) : _elem(elem) { }
	~DeleterScope() { delete _elem; }
protected:
	class_t* const _elem;
};

class VersionFormatter
{
public:
	VersionFormatter(int version) : _version(version) { }

	friend inline std::ostream& operator<<(std::ostream& os, const VersionFormatter& v);
protected:
	int _version;
};

inline std::ostream& operator<<(std::ostream& os, const VersionFormatter& v)
{
	os << ((v._version / 10000) % 10) << "." << ((v._version / 100) % 10) << "." << (v._version % 10);
	return os;
}

template <class Reader_t>
int inferFileFormatVersion(Reader_t& rd)
{
	const bool hasInput = rd.exists("input");
	if (hasInput)
		rd.pushGroup("input");

	int version = 0;

	if (rd.exists("CHROMATOGRAPHY_TYPE"))
	{
		version = formatVersionLowest;
	}
	else
	{
		if (rd.exists("model"))
		{
			Scope<Reader_t> s1(rd, "model");
			Scope<Reader_t> s2(rd, "connections");
			Scope<Reader_t> s3(rd, "switch_000");

			if (rd.isInt("CONNECTIONS"))
				version = 30000;
			else if (rd.isDouble("CONNECTIONS"))
				version = formatVersionRecent;
		}
	}

	if (hasInput)
		rd.popGroup();

	return version;
}

template <class Reader_t>
int getFileFormatVersion(Reader_t& rd)
{
	if (!rd.exists("meta"))
	{
		// Infer format from contents
		return inferFileFormatVersion(rd);
	}

	rd.pushGroup("meta");
	if (rd.exists("FILE_FORMAT"))
	{
		Scope<Reader_t> s1(rd);
		if (rd.isString("FILE_FORMAT"))
		{
			const std::string strVer = rd.getString("FILE_FORMAT");
			std::vector<std::string> vecElem;

			split(strVer, '.', vecElem);

			int version = 0;
			if (vecElem.size() >= 1)
				version += std::stoi(vecElem[0]) * 10000;
			if (vecElem.size() >= 2)
				version += std::stoi(vecElem[1]) * 100;
			if (vecElem.size() >= 3)
				version += std::stoi(vecElem[2]);

			return version;
		}
		else
		{
			return rd.getInt("FILE_FORMAT");
		}
	}
	else
	{
		// Infer format from contents
		rd.popGroup();
		return inferFileFormatVersion(rd);
	}
}




bool run(const std::string& inFileName, const std::string& inFileExt, const std::string& outFileName, const std::string& outFileExt, bool targetOldFormat)
{
	cadet::io::IFileReader* rd = cadet::io::createReader(inFileExt);
	if (!rd)
	{
		std::cerr << "Input file format ('." << inFileExt << "') not supported" << std::endl;
		return false;
	}
	DeleterScope<cadet::io::IFileReader> dsR(rd);

	cadet::io::IFileWriter* wr = cadet::io::createWriter(outFileExt);
	if (!wr)
	{
		std::cerr << "Output file format ('." << outFileExt << "') not supported" << std::endl;
		return false;
	}
	DeleterScope<cadet::io::IFileWriter> dsW(wr);

	rd->openFile(inFileName, "r");
	FileCloser<cadet::io::IFileReader> fcR(*rd);

	const int curVersion = getFileFormatVersion(*rd);
	if (curVersion <= 0)
	{
		std::cerr << "Cannot detect file format version. Exiting" << std::endl;
		return false;
	}

	std::cout << "Detected file format version " << VersionFormatter(curVersion) << std::endl;

	// Early outs
	if ((curVersion == formatVersionLowest) && targetOldFormat && cadet::util::caseInsensitiveEquals(inFileExt, outFileExt))
	{
		std::cout << "Targeting version " << VersionFormatter(formatVersionLowest) << " so we are done" << std::endl;
		return true;
	}

	if ((curVersion == formatVersionRecent) && !targetOldFormat && cadet::util::caseInsensitiveEquals(inFileExt, outFileExt))
	{
		std::cout << "Targeting version " << VersionFormatter(formatVersionRecent) << " so we are done" << std::endl;
		return true;
	}

	wr->openFile(outFileName, "co");
	FileCloser<cadet::io::IFileWriter> fcW(*wr);

	if (targetOldFormat)
	{
		IFormatConverter* const conv = createConverter(curVersion);
		DeleterScope<IFormatConverter> dsC(conv);

//		dsC.downgrade(rd, wr);
	}
	else
	{

	}

	return true;
}

int main(int argc, char** argv)
{
	std::string inFile;
	std::string outFile;
	bool targetOldFormat = false;
	try
	{
		TCLAP::CustomOutputWithoutVersion customOut("convertFile");
		TCLAP::CmdLine cmd("Converts between CADET file format versions", ' ', "1.0");
		cmd.setOutput(&customOut);

		cmd >> (new TCLAP::SwitchArg("o", "old", "Target file format version 2.0.0"))->storeIn(&targetOldFormat);
		cmd >> (new TCLAP::UnlabeledValueArg<std::string>("input", "Input file", true, "", "File"))->storeIn(&inFile);
		cmd >> (new TCLAP::UnlabeledValueArg<std::string>("output", "Output file", false, "", "File"))->storeIn(&outFile);

		cmd.parse(argc, argv);
	}
	catch (const TCLAP::ArgException &e)
	{
		std::cerr << "ERROR: " << e.error() << " for argument " << e.argId() << std::endl;
		return 1;
	}

	if (outFile == inFile)
	{
		std::cerr << "Input and output file must be distinct" << std::endl;
		return 2;
	}

	// Obtain file extensions for selecting corresponding reader and writer
	const std::size_t dotPosIn = inFile.find_last_of('.');
	if (dotPosIn == std::string::npos)
	{
		std::cerr << "Could not deduce input filetype due to missing extension: " << inFile << std::endl;
		return 2;
	}
	const std::string fileExtIn = inFile.substr(dotPosIn+1);

	// Set output file name
	if (outFile.empty())
	{
		outFile = inFile.substr(0, dotPosIn) + "-converted." + fileExtIn;		
	}

	const std::size_t dotPosOut = outFile.find_last_of('.');
	if (dotPosOut == std::string::npos)
	{
		std::cerr << "Could not deduce output filetype due to missing extension: " << outFile << std::endl;
		return 2;
	}
	const std::string fileExtOut = outFile.substr(dotPosOut+1);

	try
	{
		if (run(inFile, fileExtIn, outFile, fileExtOut, targetOldFormat))
			return 2;
	}
	catch (const cadet::io::IOException& e)
	{
		std::cerr << "IO ERROR: " << e.what() << std::endl;
		return 2;
	}
	catch (const std::exception& e)
	{
		std::cerr << "ERROR: " << e.what() << std::endl;
		return 1;
	}

	return 0;
}
