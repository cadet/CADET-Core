// =============================================================================
//  CADET
//  
//  Copyright © 2008-2020: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

#ifdef CADET_PARALLELIZE
	#include <tbb/tbb.h>

	#define TBB_PREVIEW_GLOBAL_CONTROL 1
	#include <tbb/global_control.h>
#endif


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

#ifdef CADET_PARALLELIZE
	int nThreads = tbb::task_scheduler_init::automatic;
#else
	int nThreads = 0;
#endif

	Catch::Session session;

	// Add command line option for threads to CATCH's argument parser
	session.cli(session.cli() | Catch::clara::Opt( nThreads, "number" )["--tbbthreads"]("number of TBB threads"));;

	// Parse the command line and check for error
	const int returnCode = session.applyCommandLine(argc, argv);
	if (returnCode != 0)
		return returnCode;

#ifdef CADET_PARALLELIZE
	if (nThreads <= 0)
		nThreads = tbb::task_scheduler_init::automatic;

	tbb::task_scheduler_init taskSchedulerInit(nThreads);
#endif

	// Run tests
	const int result = session.run();
	return result;
}
