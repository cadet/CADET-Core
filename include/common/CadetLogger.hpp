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

#ifndef CADETLOGGER_HPP_
#define CADETLOGGER_HPP_


#ifdef BENCHMARK_MODE
    #define LOGGING_DISABLE
#else
    // uncomment here to completely switch off logging
    //#define LOGGING_DISABLE
#endif

// uncomment to extend the logging to output own types
#define LOGGING_DEFINE_EXTENDED_OUTPUT_TYPE



#include <sstream>
#include <iomanip>

#include "logging/NullOutput.h"
#include "logging/OutputStream.h"
#include "logging/OutputLevelSwitchDisabled.h"
#include "logging/OutputLevelRunTimeSwitch.h"
#include "logging/Logger.h"

namespace cadet {
    // Embed logging namespace in cadet
    using namespace logging;
}


// Set appropriate compiler macros for windows or linux otherwise
#ifdef _WIN32
    #define CURRENT_FUNCTION (__FUNCTION__)
#else
    #define CURRENT_FUNCTION (__FUNCTION__)
//    #define CURRENT_FUNCTION (__PRETTY_FUNCTION__) // GCC: Function names including their type signature of the function as well as its bare name
#endif



#if !defined(LOGGING_DEFINE_OWN_OUTPUT_TYPE) && !defined(MATLAB_MEX_FILE)

    #ifdef LOGGING_DISABLE

        /*! \brief define NullOutput as output %device */
        LOGGING_DEFINE_OUTPUT( ::logging::NullOutput )

        // undefine the macro and redefine it as empty to switch off all
        // user created output types and ensure that only the NullOutput
        // is used as logging type
        #undef LOGGING_DEFINE_OUTPUT
        #define LOGGING_DEFINE_OUTPUT(NAME)

    #else // LOGGING_DISABLE

        #include "logging/loggingConfigGeneralPurposeOS.h"


        #ifndef LOGGING_DEFINE_EXTENDED_OUTPUT_TYPE
            LOGGING_DEFINE_OUTPUT( ::logging::LoggingType )
        #endif // LOGGING_DEFINE_EXTENDED_OUTPUT_TYPE

    #endif // LOGGING_DISABLE

#elif defined(MATLAB_MEX_FILE)
    
    #include "logging/MexOutput.h"

    /*!\brief StdLogRunTimeSwitchType is a %logging type
     *        supporting %logging to MATLAB by means of mexPrintf
     *        function with the additional feature of switching
     *        the verbosity of the %logging framework at runtime.
     */
    typedef ::logging::OutputLevelRunTimeSwitch <
                ::logging::OutputStream <
                    ::logging::MexOutput
                >
            > StdLogRunTimeSwitchType;

    typedef StdLogRunTimeSwitchType LoggingType;

#endif // LOGGING_DEFINE_OWN_OUTPUT_TYPE



#if !defined(MATLAB_MEX_FILE) && !defined(_WIN32)

    // extensions start
    /*! \brief define a simple ansi console colors struct
     *
     * The contained values correspond to ansi color definitions
     * However, the colors are only working on platforms that are
     * able to interpret these ansi escape sequences, like
     * [L|U]nix
     */
    struct Color {
        enum Colors {
            reset   = 0,
            black   = 30,
            red     = 31,
            green   = 32,
            yellow  = 33,
            blue    = 34,
            magenta = 35,
            cyan    = 36,
            white   = 37
        };
    };

    #define GENERATE_OPERATOR_CONTENT   \
            unsigned char tmp=Base::getBase();  \
            *this << ::logging::log::dec << "\033[" << static_cast<unsigned short>(l) << 'm';  \
            Base::setBase(tmp);


#elif !defined(MATLAB_MEX_FILE)

    #include <windows.h> // WinApi header

    /*! \brief define a simple console colors struct for windows terminals
     */
    struct Color {
        enum Colors {
            reset   = 7,
            blue    = 9,
            green,//= 10
            cyan, //= 11
            red,  //= ...
            magenta,
            yellow,
            white,
            black   = 0
        };
    };

    #define GENERATE_OPERATOR_CONTENT   \
        HANDLE hConsole; \
        hConsole = GetStdHandle(STD_OUTPUT_HANDLE); \
        SetConsoleTextAttribute(hConsole, l);

#else

    struct Color {
        enum Colors {
            reset   = 0,
            black,
            red,
            green,
            yellow,
            blue,
            magenta,
            cyan,
            white
        };
    };

    #define GENERATE_OPERATOR_CONTENT

#endif


/// \brief Defines the extended output type that is capable to interpret
///        the colors and produce the correct escape sequences.

template<typename Base>
struct CadetLogger: public Base
{
    /// \brief catch color type and produce the correct escape sequence

    CadetLogger& operator << (const Color::Colors l)
    {
        GENERATE_OPERATOR_CONTENT
        return *this;
    }

    /// \brief forward all unknown types to the base output type for
    ///        further processing.

    template<typename T>
    CadetLogger& operator<<(const T &t) {
        Base::operator<<(t);
        return *this;
    }

    /// \brief catch std::string type and call the correct method
    ///        for the output generation

    CadetLogger& operator<<(const std::string& s)
    {
        Base::operator<<(s.c_str());
        return *this;
    }

    /// \brief catch double type and call the correct method
    ///        for the output generation

    CadetLogger& operator << (const double& d)
    {
        std::ostringstream oss; oss << d;
        Base::operator << (oss.str().c_str());
        return *this;
    }

};
// extensions end


// Specifying the logging back-end
LOGGING_DEFINE_OUTPUT(CadetLogger<LoggingType>)


// Generate non-standard log levels
LOGGING_GENERATE_LEVEL(Except, ::logging::Level::error, "\n[ EXCEPTION ] ");

// Does not work so far
//LOGGING_GENERATE_LEVEL(Trace1, ::logging::Level::trace, "[ TRACE 1   ] ");
//LOGGING_GENERATE_LEVEL(Trace2, ::logging::Level::trace, "[ TRACE 2   ] ");
//
//LOGGING_GENERATE_LEVEL(Debug1, ::logging::Level::debug, "[ DEBUG 1   ] ");
//LOGGING_GENERATE_LEVEL(Debug2, ::logging::Level::debug, "[ DEBUG 2   ] ");


// Disable certain log levels at compile time
//LOGGING_DISABLE_LEVEL(::logging::Error);
//LOGGING_DISABLE_LEVEL(::logging::Warning);
//LOGGING_DISABLE_LEVEL(::logging::Normal);
//LOGGING_DISABLE_LEVEL(::logging::Info);

// Enabling any of the trace/debug levels may lead to ENORMOUS amounts of outputs
// and can SIGNIFICANTLY decrease the speed of the simulator !!!
LOGGING_DISABLE_LEVEL(::logging::Trace1);
LOGGING_DISABLE_LEVEL(::logging::Trace2);
LOGGING_DISABLE_LEVEL(::logging::Debug1);
LOGGING_DISABLE_LEVEL(::logging::Debug2);


#endif /* CADETLOGGER_HPP_ */
