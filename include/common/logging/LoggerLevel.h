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
 * $Id: LoggerLevel.h 1 2010-01-11 14:24:28Z perch $
 *
 ******************************************************************************/

#ifndef __LoggerLevel_h__
#define __LoggerLevel_h__

/// \brief The macro genrate %logging levels
///
/// \param LEVELNAME is the typename of the generated level.
/// \param LEVEL is the correcsponding level (info, warning, etc.)
/// \param DESC is the string that is outputed if the level is used
///
/// For example, to generate an own level for algorithm XYZ with the info
/// level you have to define it like these:
/// \code
/// LOGGING_GENERATE_LEVEL( AlgoXYZ, ::logging::Level::info, "[ AlgorithmXYZ ]");
/// ...
/// log::emit< AlgoXYZ >() << "Hello World!" << log::endl;
/// // prints -- "[ AlgorithmXYZ ] Hello World!" with a linefeed
/// \endcode
///
#define LOGGING_GENERATE_LEVEL(LEVELNAME, LEVEL, DESC)                        \
struct LEVELNAME {                                                            \
    /*! \brief delivers the current %level of %logging */                     \
    static ::logging::Level::levels level () {                                \
        return LEVEL;                                                         \
    }                                                                         \
    /*! \brief delivers the string reporting the current %level of %logging */\
    static const char * desc() {                                              \
        return DESC;                                                          \
    }                                                                         \
}

namespace logging {

    /*! \brief %Level allows for describing the current %level of %logging.
     */
    struct Level {
        /*! \brief The enum encapsulates the different %levels of the
         *         logging system in order to provide a type used for
         *         matching in some operators like operator<<.
         */
        enum levels {
            disable =   0,     // Disabled
            error   =   1,     // Errors
            warning =   2,     // Warnings
            normal  =   4,     // Standard output of the program
            info    =   8,     // Additional information on program execution
            debug1  =  16,     // Only debug messages that are needed for the current bug search
            debug2  =  32,     // All debug messages
            trace1  = 128,     // Only trace messages that are outside of loops
            trace2  = 256      // All trace messages
        };

        /*! \brief the current set level */
        unsigned char l;

        /*! \brief operator to test if a certain %level is set */
        bool operator & (Level& r) {
            return !!(l & r.l);
        }

        /*! \brief operator to assign a %level */
        Level& operator = (levels b) {
            l = b;
            return *this;
        }
    };

    /*! \brief This class is intended to be used as a template argument for
     *         the logging::log::emit() function.
     *
     *          Prefixes the output with "[ ERROR ] " and enables reporting
     *          of the current logging level.
     */
    struct Error {
        /*! \brief delivers the current %level of %logging */
        static ::logging::Level::levels level () {
            return ::logging::Level::error;
        }
        /*! \brief delivers the string reporting the current %level of %logging */
        static const char * desc() {
            return "\n[ ERROR     ] ";
        }
    };

    /*! \brief This class is intended to be used as a template argument for
     *         the logging::log::emit() function.
     *
     *         Prefixes the output with "[WARNING] " and enables reporting
     *         of the current logging level.
     */
    struct Warning {
        /*! \brief delivers the current %level of %logging */
        static ::logging::Level::levels level () {
            return ::logging::Level::warning;
        }
        /*! \brief delivers the string reporting the current %level of %logging */
        static const char * desc() {
            return "\n[ WARNING   ] ";
        }
    };

    /*! \brief This class is intended to be used as a template argument for
     *         the logging::log::emit() function.
     *
     *         Does not prefix the output, but enables reporting of the
     *         current logging level, too.
     */
    struct Void {
        /*! \brief delivers the current %level of %logging */
        static ::logging::Level::levels level () {
            return ::logging::Level::normal;
        }
        /*! \brief delivers the string reporting the current %level of %logging */
        static const char * desc() {
            return "";
        }
    };

    /*! \brief only for convenience */
    typedef Void Normal;

    /*! \brief This class is intended to be used as a template argument for
     *         the logging::log::emit() function.
     *
     *         Prefixes the output with "[ INFO  ] " and enables reporting
     *         of the current logging level.
     */
    struct Info {
        /*! \brief delivers the current %level of %logging */
        static ::logging::Level::levels level () {
            return ::logging::Level::info;
        }
        /*! \brief delivers the string reporting the current %level of %logging */
        static const char * desc() {
            return "[ INFO      ] ";
        }
    };

    /*! \brief This class is intended to be used as a template argument for
     *         the logging::log::emit() function.
     *
     *         Prefixes the output with "[ TRACE ] " and enables reporting
     *         of the current logging level.
     */
    struct Trace1 {
        /*! \brief delivers the current %level of %logging */
        static ::logging::Level::levels level () {
            return ::logging::Level::trace1;
        }
        /*! \brief delivers the string reporting the current %level of %logging */
        static const char * desc() {
            return "[ TRACE 1   ] ";
        }
    };

    /*! \brief This class is intended to be used as a template argument for
     *         the logging::log::emit() function.
     *
     *         Prefixes the output with "[ TRACE ] " and enables reporting
     *         of the current logging level.
     */
    struct Trace2 {
        /*! \brief delivers the current %level of %logging */
        static ::logging::Level::levels level () {
            return ::logging::Level::trace2;
        }
        /*! \brief delivers the string reporting the current %level of %logging */
        static const char * desc() {
            return "[ TRACE 2   ] ";
        }
    };

    /*! \brief This class is intended to be used as a template argument for
     *         the logging::log::emit() function.
     *
     *         Prefixes the output with "[ DEBUG ] " and enables reporting
     *         of the current logging level.
     */
    struct Debug1 {
        /*! \brief delivers the current %level of %logging */
        static ::logging::Level::levels level () {
            return ::logging::Level::debug1;
        }
        /*! \brief delivers the string reporting the current %level of %logging */
        static const char * desc() {
            return "[ DEBUG 1   ] ";
        }
    };

    /*! \brief This class is intended to be used as a template argument for
     *         the logging::log::emit() function.
     *
     *         Prefixes the output with "[ DEBUG ] " and enables reporting
     *         of the current logging level.
     */
    struct Debug2 {
        /*! \brief delivers the current %level of %logging */
        static ::logging::Level::levels level () {
            return ::logging::Level::debug2;
        }
        /*! \brief delivers the string reporting the current %level of %logging */
        static const char * desc() {
            return "[ DEBUG 2   ] ";
        }
    };

} /* logging */

#endif

