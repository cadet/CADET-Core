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

#ifndef _ACTIVE_HPP_
#define _ACTIVE_HPP_

// wrapper classes for active data types

#if defined(ACTIVE_ADOLC)

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

        typedef size_t idx_t;

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
            *(static_cast<adtl::adouble*> (this)) = rhs;
            return *this;
        }

        // operators
        const active& operator=(double rhs)
        {
            *(static_cast<adtl::adouble*> (this)) = rhs;
            return *this;
        }

        inline size_t gradientSize() const { return adtl::ADOLC_numDir; }

        ///todo remove this cast operator - dirty hack to cast actives to doubles in some residual functions
        // that are actually not called. Should be no problem with C++11 (using constexpr) ???
        explicit operator double() const
        {
            return getValue();
        }
    };

    namespace AD
    {
        size_t getMaxDirections();
        void setMaxDirections(size_t numDir);
    }
}

#elif defined(ACTIVE_SFAD)

#define SFAD_DEFAULT_DIR 80

#include "sfad.hpp"
#include "CadetException.hpp"

#define ACTIVE_INIT SFAD_GLOBAL_GRAD_SIZE

namespace cadet
{
    typedef sfad::Fwd<double, sfad::StackStorage> active;

    namespace AD
    {
        size_t getMaxDirections();
        void setMaxDirections(size_t n);
    }
}

#elif defined(ACTIVE_SETFAD)

#define SFAD_DEFAULT_DIR 80

#include "setfad.hpp"
#include "CadetException.hpp"

#define ACTIVE_INIT SFAD_GLOBAL_GRAD_SIZE

namespace cadet
{
    typedef sfad::FwdET<double, sfad::StackStorage> active;

    namespace AD
    {
        size_t getMaxDirections();
        void setMaxDirections(size_t n);
    }
}

#elif defined(ACTIVE_DCO)

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
