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

#ifndef _ACTIVE_HPP_
#define _ACTIVE_HPP_

// wrapper classes for active data types

#if defined ACTIVE_ADOLC

#define ADOLC_TAPELESS
#define NUMBER_DIRECTIONS 80
#include <adolc/adouble.h>
#include "CadetException.hpp"

#define ACTIVE_INIT ADOLC_TAPELESS_UNIQUE_INTERNALS

namespace cadet
{


class active: public adtl::adouble
{
public:

    // constructors
    active() :
        adtl::adouble()
    {
    }

    active(const double& d) :
        adtl::adouble(d)
    {
    }

    active(const adtl::adouble& ad) :
        adtl::adouble(ad)
    {
    }

    // operators
    const active& operator=(const adtl::adouble& rhs)
    {
        *(dynamic_cast<adtl::adouble*> (this)) = rhs;
        return *this;
    }

    // operators
    const active& operator=(double rhs)
    {
        *(dynamic_cast<adtl::adouble*> (this)) = rhs;
        return *this;
    }


    ///todo remove this cast operator - dirty hack to cast actives to doubles in some residual functions
    // that are actually not called. Should be no problem with C++11 (using constexpr) ???
    operator double() const
    {
        throw CadetException("Cast from active to double is not allowed!");
    }

    static int getMaxDirections() { return adtl::ADOLC_numDir; }

};

}

#elif defined ACTIVE_DCO

#include dco.h

#define ACTIVE_INIT

namespace cadet
{

    class active : public dco::active
    {

    };

}

#else

#error No active data type defined!

#endif // #if defined ACTIVE_
#endif
