# =============================================================================
#  CADET - The Chromatography Analysis and Design Toolkit
#  
#  Copyright © 2008-2018: The CADET Authors
#            Please see the AUTHORS and CONTRIBUTORS file.
#  
#  All rights reserved. This program and the accompanying materials
#  are made available under the terms of the GNU Public License v3.0 (or, at
#  your option, any later version) which accompanies this distribution, and
#  is available at http://www.gnu.org/licenses/gpl.html
# =============================================================================

function(cadet_target_compile_features TARGET)
	target_compile_features(${TARGET} PUBLIC cxx_alias_templates cxx_defaulted_functions cxx_delegating_constructors
	    cxx_explicit_conversions cxx_inheriting_constructors cxx_rvalue_references cxx_static_assert
	    cxx_lambdas cxx_nullptr cxx_auto_type cxx_range_for cxx_right_angle_brackets cxx_deleted_functions cxx_nullptr 
	    cxx_strong_enums cxx_uniform_initialization cxx_template_template_parameters cxx_defaulted_move_initializers)
endfunction()

function(cadet_apply_compile_options TARGET)
	if (WIN32)
		# Disable Windows headers defining min and max symbol in global namespace
		target_compile_definitions(${TARGET} PRIVATE -DNOMINMAX)
	endif()

	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
		if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
			if (WIN32)
				target_compile_options(${TARGET} PRIVATE /W4 /Wall /Qstd=c++11 /D_DEBUG)
			else()
#				set_target_properties (${TARGET} PROPERTIES LINK_FLAGS "-mkl")
				target_compile_options(${TARGET} PRIVATE -w3 -Wall -std=c++11 -D_DEBUG) # -fno-rtti
			endif()
		elseif ((CMAKE_CXX_COMPILER_ID STREQUAL "Clang") OR (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"))
			target_compile_options(${TARGET} PRIVATE -Wall -Wvla -pedantic) # -fno-rtti
		elseif ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR CMAKE_COMPILER_IS_GNUCXX)
			target_compile_options(${TARGET} PRIVATE -Wall -Wvla -pedantic) # -fno-rtti
		elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
			# GR-= Disable RTTI
			target_compile_options(${TARGET} PRIVATE "/Gw /GR-")
		endif ()
	else()
		if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
			if (WIN32)
				set_target_properties (${TARGET} PROPERTIES LINK_FLAGS "")
				target_compile_options(${TARGET} PRIVATE /W4 /Qstd=c++11 /O3 /Qipo- /Qprec-div- /fp:fast=2 /Qansi-alias /Qalias-args-) # /GR-
			else()
#				set_target_properties (${TARGET} PROPERTIES LINK_FLAGS "-mkl")
				target_compile_options(${TARGET} PRIVATE -w3 -std=c++11 -O3 -no-ipo -no-prec-div -static -fp-model fast=2 -ansi-alias -fargument-noalias) # -fno-rtti
			endif()
		elseif ((CMAKE_CXX_COMPILER_ID STREQUAL "Clang") OR (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"))
			target_compile_options(${TARGET} PRIVATE -Wall -Wvla -pedantic) # -fno-rtti
		elseif ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR CMAKE_COMPILER_IS_GNUCXX)
			target_compile_options(${TARGET} PRIVATE -Wall -Wvla -pedantic) # -fno-rtti
		elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
			# GR-= Disable RTTI, LTCG = Link time code generation, /GL = Whole program optimization
			set_target_properties (${TARGET} PROPERTIES LINK_FLAGS "/LTCG")
			target_compile_options(${TARGET} PRIVATE "/GL /Gw") # "/GR-"
		endif ()
	endif()
endfunction()

function(cadet_choose_ad_lib TARGET)
	if (ADLIB STREQUAL "adolc")
		target_compile_definitions(${TARGET} PRIVATE -DACTIVE_ADOLC)
		target_include_directories(${TARGET} PRIVATE "${CMAKE_SOURCE_DIR}/ThirdParty/ADOL-C/include")
	elseif (ADLIB STREQUAL "sfad")
		target_compile_definitions(${TARGET} PRIVATE -DACTIVE_SFAD)
		target_include_directories(${TARGET} PRIVATE "${CMAKE_SOURCE_DIR}/include/ad")
	elseif (ADLIB STREQUAL "setfad")
		target_compile_definitions(${TARGET} PRIVATE -DACTIVE_SETFAD)
		target_include_directories(${TARGET} PRIVATE "${CMAKE_SOURCE_DIR}/include/ad")
	endif ()
endfunction()
