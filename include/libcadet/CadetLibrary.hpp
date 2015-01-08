// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2015: Eric von Lieres¹, Joel Andersson,
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

#ifndef LIBCADET_CADETLIBRARY_HPP_
#define LIBCADET_CADETLIBRARY_HPP_

// Export and import functions when using MS Visual Studio compiler
#ifndef CADET_API
    #ifdef _MSC_VER
        #if defined(libcadet_shared_EXPORTS) || defined(libcadet_mex_EXPORTS)
            #define CADET_API _declspec(dllexport)
        #else
            #define CADET_API _declspec(dllimport)
        #endif
    #else
        #define CADET_API
    #endif
#endif

extern "C"
{
    //! \brief Returns the version string of the libcadet library
    CADET_API const char* cadetGetLibraryVersion();

    //! \brief Returns the git commit hash of the source which was used to build the binaries
    CADET_API const char* cadetGetLibraryCommitHash();

    //! \brief Returns the git refspec of the source which was used to build the binaries
    CADET_API const char* cadetGetLibraryBranchRefspec();

    //! \brief Resets global variables
    CADET_API void cadetResetGlobals();
}

namespace cadet
{

    //! \brief Returns the version string of the libcadet library (same as cadetGetLibraryVersion())
    CADET_API const char* getLibraryVersion();

    //! \brief Returns the git commit hash of the source which was used to build the binaries (same as cadetGetLibraryCommitHash())
    CADET_API const char* getLibraryCommitHash();

    //! \brief Returns the git refspec of the source which was used to build the binaries
    CADET_API const char* getLibraryBranchRefspec();

    //! \brief Resets global variables (same as cadetResetGlobals())
    CADET_API void resetGlobals();

} // namespace cadet

#endif  // LIBCADET_CADETLIBRARY_HPP_
