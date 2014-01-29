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
 * $Id: Logger.h 2 2010-01-12 10:13:58Z perch $
 *
 ******************************************************************************/

#ifndef __Logger_h__
#define __Logger_h__

#include "logging/singleton.h"

#include "logging/NullOutput.h"
#include "logging/LoggerLevel.h"

/*! \brief The namespace %logging contains everything that is related to the
 *         %logging framework (\ref loggingframework) and encloses it to avoid
 *         name interference.
 *
 *         All user-accessible classes and functions (the interface) are
 *         located there.
 */
namespace logging {

    /*! \brief defines the return type of the logging framework that is able to process
     *         different type/value by using the typical operator<< behaviour.
     */
    struct loggingReturnType;

    /*! \brief The %detail namespace shields all details from the upper %logging namespace
     *         in order to separate implementation details from the interface namespace
     *         named %logging.
     */
    namespace detail {

        /*! \brief Handles the type and an instance of the class that performs the %logging.
         *
         *  \tparam Level describes the %logging %Level. It could be logging::Trace,
         *          logging::Info, ... or a user defined type, that fulfils the
         *          correct interface.
         *  \tparam R represents the return_type aka the output type, that is responsible
         *          for the %logging.
         */
        template<typename Level = ::logging::Void, typename R = loggingReturnType>
        class Logger {
            public:
                /*! \brief The typedef is made for function from outer scope,
                 *         determining the return type of the Logger::logging()
                 *         method.
                 */
                typedef R return_type;

                /*! \brief The method gives the held instance of the output %object.
                 *
                 *  \return A reference to the held %object.
                 */
                static return_type& logging () {
                    return Obj<typename return_type::output_base_type, int>::obj();
                }

            private:

                /*! \brief provides a singleton of the output %object
                 *
                 *
                 *  \tparam log_t represents the type of the held output %object
                 *  \tparam T is not used, but is needed because full
                 *          specialisation is not allowed in class scope,
                 *          however partial specialisation is.
                 */
                template <typename log_t, typename T>
                struct Obj {
                    static return_type& obj () {
                        typedef singleton<return_type> log_output;
                        return log_output::instance();
                    }
                };

                /*! \brief It is a specialisation of Logger::Obj. It is used to disable
                 *         the output by delivering always an instance to a NullOutput.
                 *
                 *  \copydoc Obj
                 */
                template <typename T>
                struct Obj<NullOutput, T> {
                    static return_type& obj () {
                        return *reinterpret_cast<return_type*>(0x0);
                    }
                };

        };

        /*! \brief output the given data according to output_base_type's policy
         */
        template<typename TT, typename M>
        static inline void output(TT& t, const M& m){
            static_cast<typename TT::output_base_type&>(t)<<m;
        }

    } /* detail */

    /*! \brief free operator function to dive into the framework
     */
    template<typename T>
    static inline ::logging::loggingReturnType &
    operator<<( ::logging::loggingReturnType& ret, const T& t) {
        ::logging::detail::output(ret,t);
        return ret;
    }


    /*! \brief free operator function overload for disabled parts
     *         of the framework
     */
    template<typename T>
    static inline ::logging::NullOutput&
    operator << (::logging::NullOutput& no, T) {
        return no;
    }

    /*! \brief The struct %log acts as interface to the %logging framework.
     *
     *  Some definitions are provided for convenience, because
     *  it allows the same use case like with the std::streams.
     *
     *  Here is an example and you see it is very intuitive to use
     *  as it is with the std::streams.
     *  \code
     *  log::emit() << log::hex << 15 << log::endl;
     *  // prints 0xf with a line feed
     *  \endcode
     */
    struct log {
        /*! \brief Numerative manipulator definitions */
        enum Numerative {
            bin  =  2,   ///< switch to binary output
            oct  =  8,   ///< switch to octal output
            dec  = 10,   ///< switch to decimal output
            hex  = 16    ///< switch to hexadecimal output
        };

        /*! \brief Manipulator definitions */
        enum Manipulator {
            tab  =  '\t',   ///< prints a tabulator to the output
            endl =  '\n'    ///< adds a line feed output
        };

        /*! \brief The emit method gets the to logged information and emits it to the
         *         output (sink) of the %logging framework
         *
         * \tparam  Level describes the %logging %Level. It could be logging::Trace,
         *          logging::Info, ... or a user defined type, that fulfils the the
         *          correct interface.
         *
         * \returns a type, to that can be logged, but it behaviour depends on the
         *          configuration
         *
         * You can use the log::emit method as follows:
         * \code
         * log::emit() << "Hello World!" << log::endl;
         * // prints "Hello World!" with a line feed.
         *
         * log::emit< ::logging::Info>() << "Hello World!" << log::endl;
         * //prints "[ INFO  ] Hello World!" with a line feed.
         * \endcode
         */
        template<typename Level>
        static inline typename ::logging::detail::Logger<Level>::return_type& emit () {
            return ::logging::detail::Logger<Level>::logging()
                   << Level::level() << Level::desc();
        }

        static inline ::logging::detail::Logger<>::return_type& emit () {
            return ::logging::detail::Logger<>::logging() << ::logging::Void::level();
        }
    };

} /* logging */


/*! \brief The macro generates the correct return type of the logging framework
 *
 *  \param BASE describes the base type of the return type (loggingReturnType)
 *
 *  This macro is used if someone extend the logging framework with further
 *  functionality like customized prefixes or colorization. For an example
 *  and howto read \ref extending
 */
#define LOGGING_DEFINE_OUTPUT(BASE)                                           \
namespace logging {                                                           \
    struct loggingReturnType : public BASE {                                  \
            /*! \brief The provided typedef is used for compile time          \
             *         selection of different implementation of the           \
             *         %logging framework. Thus, it is necessary              \
             *         that any output type supports this type                \
             *         definition, why it is defined here.                    \
             */                                                               \
            typedef BASE    output_base_type;                                 \
        };                                                                    \
}

/*! \brief The macro enables partial disabling of certain %logging levels
 *
 *  \param NAME describes the disabled output level.
 *
 *  For example, to disable the info level you have to define it like these:
 *  \code
 *  LOGGING_DISABLE_NAME( ::logging::INFO);
 *  ...
 *  log::emit< ::logging::Info>() << "Hello World!" << log::endl;
 *  // prints and evaluates to nothing, thus no run-time or other
 *  // ressources have to be paid
 *  \endcode
 */
#define LOGGING_DISABLE_LEVEL(NAME)                                           \
namespace logging {                                                           \
    namespace detail {                                                        \
        template<typename Level, typename R>                                  \
        class Logger;                                                         \
                                                                              \
        template <>                                                           \
        class Logger<NAME, ::logging::loggingReturnType>;                     \
    } /* detail */                                                            \
} /* logging */                                                               \
                                                                              \
template<>                                                                    \
class logging::detail::Logger<NAME, ::logging::loggingReturnType> {           \
    public:                                                                   \
        typedef ::logging::NullOutput return_type;                            \
        static return_type& logging (){                                       \
            return *reinterpret_cast<return_type*>(0x0);                      \
        }                                                                     \
}

#endif

