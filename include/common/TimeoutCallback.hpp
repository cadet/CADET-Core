// =============================================================================
//  CADET
//
//  Copyright © The CADET Authors
//            Please see the CONTRIBUTORS.md file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

/**
 * @file
 * Provides a callback that terminates time integration when a given amount of
 * time has been exceeded.
 */

#ifndef CADET_TIMEOUT_CALLBACK_HPP_
#define CADET_TIMEOUT_CALLBACK_HPP_
 
#include "common/Timer.hpp"
#include "cadet/Notification.hpp"

namespace cadet
{

	class TimeoutCallback : public cadet::INotificationCallback
	{
	public:
		TimeoutCallback() { }
		virtual ~TimeoutCallback() CADET_NOEXCEPT { }

		void setTimeout(double timeout) { _timeout = timeout; }

		virtual void timeIntegrationStart()
		{
			_timer.start();
		}
		virtual void timeIntegrationEnd()
		{
			_timer.reset();
		}

		virtual void timeIntegrationError(char const* message, unsigned int section, double time, double progress) { }

		virtual bool timeIntegrationSection(unsigned int section, double time, double const* state, double const* stateDot, double progress)
		{
			return shouldStop();
		}

		virtual bool timeIntegrationStep(unsigned int section, double time, double const* state, double const* stateDot, double progress)
		{
			return shouldStop();
		}

		virtual bool timeIntegrationLinearSolve(unsigned int section, double time, double const* state, double const* stateDot)
		{
			return shouldStop();
		}
	protected:
		bool shouldStop()
		{
			_timer.stop();
			const bool cont_sim = _timer.totalElapsedTime() <= _timeout;
			_timer.start();
			return cont_sim;
		}

		Timer _timer;
		double _timeout;
	};
}

#endif  // CADET_TIMEOUT_CALLBACK_HPP_
