/*******************************************************************************
 *
 * Copyright (c) 2008, 2009 Michael Schulze <mschulze@ivs.cs.uni-magdeburg.de>
 * All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without
 *    modification, are permitted provided that the following conditions
 *    are met:
 *
 *    * Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *
 *    * Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *
 *    * Neither the name of the copyright holders nor the names of
 *      contributors may be used to endorse or promote products derived
 *      from this software without specific prior written permission.
 *
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 *    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 * $Id: OutputLevelRunTimeSwitch.h 1 2010-01-11 14:24:28Z perch $
 *
 ******************************************************************************/

#ifndef __OutputLevelRunTimeSwitch_h__
#define __OutputLevelRunTimeSwitch_h__

#include "logging/LoggerLevel.h"

namespace logging {

    /*! \brief Treatment of logging levels can be switched at runtime
     */
    template < typename Base >
    class OutputLevelRunTimeSwitch  : public Base {
            ::logging::Level _level;
            ::logging::Level _current;
        public:
            OutputLevelRunTimeSwitch() {
                _level.l   = ::logging::Level::error |
                             ::logging::Level::warning |
                             ::logging::Level::normal;
                _current.l = ::logging::Level::normal;
            }

            /*! \brief Output is only allowed if the current %level is switched
             *         on in the general %level.
             */
            bool allowed() {
                return !!(_level&_current);
            }

            /*! \brief Matches only on correct type and set the
             *         current %level for the output.
             */
            OutputLevelRunTimeSwitch& operator<<(const ::logging::Level::levels& l) {
                _current = l;
                return *this;
            }

            /*! \brief The operator matches on every type, and delegates further
             *         work to the base class if the output is currently
             *         allowed.
             */
            template< typename T>
            OutputLevelRunTimeSwitch& operator<<(T t) {
                if ( allowed() )
                    Base::operator<<(t);
                return *this;
            }

            OutputLevelRunTimeSwitch& operator-=(const ::logging::Level::levels& l) {
                if (!!(_level.l & l))
                    _level.l ^= l;
                return *this;
            }

            OutputLevelRunTimeSwitch& operator+=(const ::logging::Level::levels& l) {
                if (!(_level.l & l))
                    _level.l |= l;
                return *this;
            }
    };

} /* logging */

#endif
