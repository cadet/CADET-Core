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

#include "Logging.hpp"

#ifndef CADET_LOGGING_DISABLE
	namespace
	{
		/**
		 * @brief Receiver of all log messages created in the libcadet library
		 */
		cadet::ILogReceiver* logReceiver = nullptr;
	}

	template <>
	cadet::LogLevel cadet::log::RuntimeFilteringLogger<cadet::log::GlobalLogger>::_minLvl = cadet::LogLevel::Trace;

	#ifdef __clang__
		// Silence -Wundefined-var-template warning
		template class cadet::log::RuntimeFilteringLogger<cadet::log::GlobalLogger>;
	#endif
#endif

namespace cadet
{

#ifdef CADET_LOGGING_DISABLE

	void setLogReceiver(ILogReceiver* const recv) { }
	void setLogLevel(LogLevel lvl) { }
	LogLevel getLogLevel() { return LogLevel::None; }

#else

	namespace log
	{
		void emitLog(const char* file, const char* func, const unsigned int line, LogLevel lvl, const char* message)
		{
			if (logReceiver)
				logReceiver->message(file, func, line, lvl, to_string(lvl), message);
		}
	}

	void setLogReceiver(ILogReceiver* const recv)
	{
		logReceiver = recv;
	}

	void setLogLevel(LogLevel lvl)
	{
		cadet::log::RuntimeFilteringLogger<cadet::log::GlobalLogger>::level(lvl);
	}

	LogLevel getLogLevel()
	{
		return cadet::log::RuntimeFilteringLogger<cadet::log::GlobalLogger>::level();
	}

#endif

} // namespace cadet

extern "C"
{
	void cadetSetLogReceiver(cadet::ILogReceiver* const recv)
	{
		cadet::setLogReceiver(recv);
	}

	void cadetSetLogLevel(unsigned int lvl)
	{
		cadet::setLogLevel(static_cast<cadet::LogLevel>(lvl));
	}

	unsigned int cadetGetLogLevel()
	{
		return static_cast<typename std::underlying_type<cadet::LogLevel>::type>(cadet::getLogLevel());
	}
}
