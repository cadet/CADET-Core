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

#ifndef CADETTOOLS_FORMATCONVERTER_HPP_
#define CADETTOOLS_FORMATCONVERTER_HPP_

namespace cadet
{
namespace io
{
	class IReader;
	class IWriter;
	class IFileReader;
	class IFileWriter;
	class IMemoryIO;

} // namespace io
} // namespace cadet

class IFormatConverter
{
public:
	IFormatConverter() { }
	virtual ~IFormatConverter() { }

	virtual void downgrade(cadet::io::IFileReader& rd, cadet::io::IFileWriter& wr) = 0;
	virtual void upgrade(cadet::io::IReader& rd, cadet::io::IMemoryIO& mem) = 0;
	virtual void upgrade(cadet::io::IMemoryIO& rd, cadet::io::IMemoryIO& mem) = 0;
	virtual void write(cadet::io::IMemoryIO& rd, cadet::io::IFileWriter& wr) = 0;

	virtual int sourceVersion() const = 0;
	virtual int targetVersion() const = 0;
};


IFormatConverter* createConverter(int sourceVersion);

#endif /* CADETTOOLS_FORMATCONVERTER_HPP_ */
