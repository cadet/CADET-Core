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
 * $Id: OutputStream.h 1 2010-01-11 14:24:28Z perch $
 *
 ******************************************************************************/

#ifndef __OutputStream_h__
#define __OutputStream_h__

#include "logging/Logger.h"

namespace logging {

    /*! \brief The %OutputStream provides different overloading of operator,
     *         enabling the output of different types in an uniform way.
     *
     *         The behaviour of the %OutputStream is the same as known from the
     *         standard streams. The %OutputStream provides manipulators too,
     *         and is very easy to extend, enabling the adaptation by the user.
     *
     * \tparam Base is a type used for the output and usually it is an output
     *         %device.
     */
    template< typename Base>
    class OutputStream : public Base {
            /*!\brief basis for display of digits eg. 2, 8, 10 or 16
             */
            unsigned char base;

            /*!\brief forwards the outputed character to the base class,
             *        that usually performs the output to a %device.
             */
            void put (char c) {
                Base::operator<<(c);
            }
        protected:
            /*!\brief Returns the current active base for number conversions
             *
             * \return the current base of conversions
             */
            unsigned char getBase() {return base;}

            /*!\brief Set the base for number conversions
             *
             * \param tmp the new base for conversions
             */
            void setBase(unsigned char tmp){base=tmp;}

        public:
            /*! \brief Default constructor initialising with dezimal system
             */
            OutputStream () {
                base = 10;
            }

            /*!\brief outputs a characater to the output stream
             *
             * Operator << overloading: Is used to convert the given datatype
             * into a string that can be printed on an output %device.  This
             * operator has to be implemented for every standard data type
             * (char, unsigned char, short, unsigned short, int, unsigned int,
             * long, unsigned long, void*, char*).
             *
             * \param c the character, that is output
             * \return %OutputStream& allows for chaining of operators
             */
            OutputStream& operator << (char c) {
                put(c);
                return *this;
            }

            /*!\brief outputs an unsigned char
             *
             * \param c the unsigned character, that is output
             * \return %OutputStream& allows for chaining of operators
             */
             OutputStream& operator << (unsigned char c) {
                return *this << (char) c;
            }

            /*!\brief Operator that set the numerative
             *
             * \param n is the numerative to be actual unitl a new is set
             * \return %OutputStream& allows for chaining of operators
             */
            OutputStream& operator << (const ::logging::log::Numerative n) {
                base = n;
                return *this;
            }

            /*!\brief Operator for catching manipulators
             *
             * \param m is the manipulator that is feed to the output stream
             * \return %OutputStream& allows for chaining of operators
             */
            OutputStream& operator << (const ::logging::log::Manipulator m) {
                *this << static_cast<char>(m);
                return *this;
            }

            /*!\brief outputs a string
             *
             * \param string the character string, that is output
             * \return %OutputStream& allows for chaining of operators
             */
            OutputStream& operator << (char* string) {
                char* pos = string;
                while (*pos) {
                    put (*pos);
                    pos++;
                }
                return *this;
            }

            /*!\brief outputs a constant string
             *
             * \param string the constant character string, that is output
             * \return %OutputStream& allows for chaining of operators
             */
            OutputStream& operator << (const char* string) {
                return *this << const_cast<char *>(string);
            }

            /*!\brief outputs a string of unsigned character
             *
             * \param string the unsigned character string, that is output
             * \return %OutputStream& allows for chaining of operators
             */
            OutputStream& operator << (unsigned char* string) {
                return *this << reinterpret_cast<char*>(string);
            }

            /*!\brief outputs a constant string of unsigned character
             *
             * \param string the unsigned character string, that is output
             * \return %OutputStream& allows for chaining of operators
             */
            OutputStream& operator << (const unsigned char* string) {
                return *this << reinterpret_cast<char*>(const_cast<unsigned char*>(string));
            }

            /*!\brief displays a short using the set numerative
             *
             * \param ival the short, that is output
             * \return %OutputStream& allows for chaining of operators
             */
            OutputStream& operator << (short ival) {
                return *this << (long) ival;
            }

            /*!\brief displays a short using the set numerative
             *
             * \param ival the unsigned short, that is output
             * \return %OutputStream& allows for chaining of operators
             */
            OutputStream& operator << (unsigned short ival) {
                return *this << (unsigned long) ival;
            }

            /*!\brief displays an integer using the set numerative
             *
             * \param ival the integer, that is output
             * \return %OutputStream& allows for chaining of operators
             */
            OutputStream& operator << (int ival) {
                return *this << (long) ival;
            }

            /*!\brief displays an unsigned short using the set numerative
             *
             * \param ival the unsigned integer, that is output
             * \return %OutputStream& allows for chaining of operators
             */
            OutputStream& operator << (unsigned int ival) {
                return *this << (unsigned long) ival;
            }

            /*!\brief displays an integer using the set numerative
             *
             * \param ival the integer, that is output
             * \return %OutputStream& allows for chaining of operators
             */
            OutputStream& operator << (long ival) {
                return *this << (long long) ival;
            }

            /*!\brief displays an unsigned short using the set numerative
             *
             * \param ival the unsigned integer, that is output
             * \return %OutputStream& allows for chaining of operators
             */
            OutputStream& operator << (unsigned long ival) {
                return *this << (unsigned long long) ival;
            }

            /*!\brief displays a long using the set numerative
             *
             * \param ival the long, that is output
             * \return %OutputStream& allows for chaining of operators
             */
            OutputStream& operator << (long long ival) {
                // if value is negative a minus is outputed first
                if (ival < 0) {
                    put ('-');
                    ival = -ival;
                }
                // than the absolute value of the digit is outputed
                return *this << (unsigned long) ival;
            }

            /*!\brief displays an unsigned long using the set numerative
             *
             * \param ival the unsigned long, that is output
             * \return %OutputStream& allows for chaining of operators
             */
            OutputStream& operator << (unsigned long long ival) {
                unsigned long long div;
                char digit;

                switch (base) {
                    case 2: put('b');break; // binary digits start with a b
                    case 8: put('0');break; // oktal digits start with a NULL
                    case 16:{
                            put ('0');              // hexadezimal digits start with 0x
                            put ('x');
                            break;
                    }
                }

                // computes the max power of the choosen basis, that is smaler than the value
                // of the digit
                for (div = 1; ival / div >= (unsigned long long) base; div *= base);

                // prints the digit character after character
                for (; div > 0; div /= (unsigned long long) base) {
                    digit = ival / div;
                    if (digit < 10)
                        put ('0' + digit);
                    else
                        put ('a' + digit - 10);
                    ival %= div;
                }
                return *this;
            }

            /*!\brief displays a boolean type as 1 (true) or 0 (false)
             *
             * \param b the boolean, that is output
             * \return %OutputStream& allows for chaining of operators
             */
            OutputStream& operator << (bool b) {
                b ? put('1') : put('0');
                return *this;
            }

            /*!\brief displays pointer as hexadezimal digit
             *
             * \param ptr the pointer, that is output
             * \return %OutputStream& allows for chaining of operators
             */
            OutputStream& operator << (void* ptr) {
                int oldbase = base;
                base = 16;
                *this << (unsigned long) ptr;
                base = oldbase;
                return *this;
            }

            /*! \brief The operator matches on every type, and provides an
             *         empty implementation.
             *         The compiler see the empty method or chain of empty
             *         methods, and throwing these away if compiling with
             *         optimizations. However, this operator is only chosen if
             *         no other operator match and thus it catch only unknown
             *         types.
             *
             * \return %OutputStream& allows for chaining of operators
             */
            template<typename T>
            OutputStream& operator << (T) {
                char Logging_Framework_OutputStream;
                char swallowing_an_unsupported_type_leading_to_none_output_of_these_information;
                return *this;
            }

            /*! \brief enable calling of manipulator functions
             *
             * \param f is a function that is called with the %OutputStream&
             *        as parameter.
             * \return %OutputStream& allows for chaining of operators
             */
            OutputStream& operator << (OutputStream& (*f) (OutputStream&)) {
                return f(*this);
            }

    };

} /* logging */

#endif
