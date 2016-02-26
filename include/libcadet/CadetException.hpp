// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2016: Eric von Lieres¹, Joel Andersson¹,
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

#ifndef LIBCADET_CADETEXCEPTION_HPP_
#define LIBCADET_CADETEXCEPTION_HPP_

#include <string>

namespace cadet
{

class CadetException
{
public:
    CadetException() { }

    CadetException(const std::string & msg) : _msg(msg) { }

    inline const std::string & msg() const { return _msg; }

private:
    std::string _msg;
};

//class CadetException: public std::exception
//{
//public:
//    my_exception(const char* file, const char* line)
//    {
//        log(file, line);
//    }
//};
//... throw my_exception(__FILE__, __LINE__);

} // namespace cadet

#endif // LIBCADET_CADETEXCEPTION_HPP_
