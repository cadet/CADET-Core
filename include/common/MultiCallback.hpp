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
 * Provides a callback that forwards to multiple other callbacks.
 */

#ifndef CADET_MULTI_CALLBACK_HPP_
#define CADET_MULTI_CALLBACK_HPP_
 
#include <vector>
#include "cadet/Notification.hpp"

namespace cadet
{

	class MultiCallback : public cadet::INotificationCallback
	{
	public:
		MultiCallback() { }
		virtual ~MultiCallback() CADET_NOEXCEPT { }

		virtual void timeIntegrationStart()
		{
			for (auto callback : _callbacks)
				callback->timeIntegrationStart();
		}

		virtual void timeIntegrationEnd()
		{
			for (auto callback : _callbacks)
				callback->timeIntegrationEnd();
		}

		virtual void timeIntegrationError(char const* message, unsigned int section, double time, double progress)
		{
			for (auto callback : _callbacks)
				callback->timeIntegrationError(message, section, time, progress);
		}

		virtual bool timeIntegrationSection(unsigned int section, double time, double const* state, double const* stateDot, double progress)
		{
			bool shouldContinue = true;
			for (auto callback : _callbacks)
				shouldContinue &= callback->timeIntegrationSection(section, time, state, stateDot, progress);

			return shouldContinue;
		}

		virtual bool timeIntegrationStep(unsigned int section, double time, double const* state, double const* stateDot, double progress)
		{
			bool shouldContinue = true;
			for (auto callback : _callbacks)
				shouldContinue &= callback->timeIntegrationStep(section, time, state, stateDot, progress);

			return shouldContinue;
		}

		virtual bool timeIntegrationLinearSolve(unsigned int section, double time, double const* state, double const* stateDot)
		{
			bool shouldContinue = true;
			for (auto callback : _callbacks)
				shouldContinue &= callback->timeIntegrationLinearSolve(section, time, state, stateDot);

			return shouldContinue;
		}

		std::vector<cadet::INotificationCallback*> const& callbacks() const CADET_NOEXCEPT { return _callbacks; }
		std::vector<cadet::INotificationCallback*>& callbacks() CADET_NOEXCEPT { return _callbacks; }

		void addCallback(cadet::INotificationCallback* callback) { _callbacks.push_back(callback); }
		void clear() { _callbacks.clear(); }
	protected:
		std::vector<cadet::INotificationCallback*> _callbacks;
	};
}

#endif  // CADET_MULTI_CALLBACK_HPP_
