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
 * $Id: FileOutput.h 1 2010-01-11 14:24:28Z perch $
 *
 ******************************************************************************/

#ifndef __FileOutput_h__
#define __FileOutput_h__

#include <iostream>
#include <fstream>

namespace logging {

    /*! \brief Provides an interface to the standard
     *         file output facilities of e.g. Linux
     *         or Windows.
     */
    class FileOutput {
            /*! \brief a plain filebufer */
            std::filebuf fb;
            /*! \brief and also a plain ostream */
            std::ostream os;
        public:
            /*! \brief The constructor opens the file, that should
             *         contain the output for writing.
             *
             *  \todo The name of the file should be configurable
             *        via a command line parameter.
             */
            FileOutput() : os(&fb) {
                fb.open ("log.txt", std::ios::out);
            }

            /*! \brief The destructor closes the file and therewith
             *         it should be persitent.
             */
            ~FileOutput() {
                fb.close();
            }

            /*! \brief operator that can output a simple character.
             *
             * \param c the character that needs to be outputed
             * \return a reference to itself allowing chaning of
             *         opertor<< calls.
             */
            FileOutput & operator<<(const char c) {
                os << c;
                os.flush();
                fb.pubsync();
                return *this;
            }
    };

} /* logging */

#endif
