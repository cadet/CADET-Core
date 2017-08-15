// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2017: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

/**
 * @file 
 * Provides helper functions for ParameterId containers.
 */

#ifndef LIBCADET_PARAMIDUTIL_HPP_
#define LIBCADET_PARAMIDUTIL_HPP_

#include <ostream>
#include <sstream>
#include "cadet/ParameterId.hpp"

namespace std
{
	template<>
	struct hash<cadet::ParameterId>
	{
		inline size_t operator()(const cadet::ParameterId& pId) const CADET_NOEXCEPT
		{
			return cadet::hashParameter(pId);
		}
	};
} // namespace std

namespace cadet
{

	inline const bool operator==(const ParameterId& a, const ParameterId& b) CADET_NOEXCEPT
	{
		return (a.name == b.name) && (a.unitOperation == b.unitOperation) && (a.component == b.component)
			&& (a.boundPhase == b.boundPhase) && (a.reaction == b.reaction) && (a.section == b.section);
	}

	inline const bool operator!=(const ParameterId& a, const ParameterId& b) CADET_NOEXCEPT { return !(a == b); }

	inline const bool operator<(const ParameterId& a, const ParameterId& b) CADET_NOEXCEPT
	{
		return std::tie(a.name, a.unitOperation, a.component, a.boundPhase, a.reaction, a.section) < std::tie(b.name, b.unitOperation, b.component, b.boundPhase, b.reaction, b.section);
	}

	inline std::ostream& operator<<(std::ostream& out, const ParameterId& pId)
	{
		out << "{" << hashParameter(pId) << " = " << pId.name << ", Unit " << static_cast<unsigned int>(pId.unitOperation) << " Comp " << static_cast<unsigned int>(pId.component)
		     << " BoundPhase " << static_cast<unsigned int>(pId.boundPhase) << " Reaction " << static_cast<unsigned int>(pId.reaction) << " Section " << static_cast<unsigned int>(pId.section) << "}";
		return out;
	}

	inline std::string to_string(const ParameterId& pId)
	{
		std::stringstream out;
		out << "{" << hashParameter(pId) << " = " << pId.name << ", Unit " << static_cast<unsigned int>(pId.unitOperation) << " Comp " << static_cast<unsigned int>(pId.component)
		     << " BoundPhase " << static_cast<unsigned int>(pId.boundPhase) << " Reaction " << static_cast<unsigned int>(pId.reaction) << " Section " << static_cast<unsigned int>(pId.section) << "}";
		return out.str();
	}

} // namespace cadet

#endif  // LIBCADET_PARAMIDUTIL_HPP_
