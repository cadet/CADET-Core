// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2018: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#define CATCH_CONFIG_RUNNER
#include <catch.hpp>


// Uncomment the next line to enable logging output of CADET in unit tests
//#define CADETTEST_ENABLE_LOG


#ifdef CADETTEST_ENABLE_LOG
	#include "cadet/Logging.hpp"
	#include <iostream>

	class LogReceiver : public cadet::ILogReceiver
	{
	public:
		LogReceiver() { }

		virtual void message(const char* file, const char* func, const unsigned int line, cadet::LogLevel lvl, const char* lvlStr, const char* message)
		{
			std::cout << '[' << lvlStr << ": " << func << "::" << line << "] " << message << std::flush;
		}
	};
#endif

int main(int argc, char* argv[])
{
#ifdef CADETTEST_ENABLE_LOG
	// Set LogLevel in CADET library
	cadet::LogLevel logLevel = cadet::LogLevel::Trace;
	LogReceiver lr;
	cadetSetLogReceiver(&lr);
	cadetSetLogLevel(static_cast<typename std::underlying_type<cadet::LogLevel>::type>(logLevel));
#endif

	const int result = Catch::Session().run(argc, argv);
	return (result < 0xff ? result : 0xff);
}
