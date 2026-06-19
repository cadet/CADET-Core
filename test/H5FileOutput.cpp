// SPDX-License-Identifier: AGPL-3.0-or-later
// =================================================================================
//  CADET
//
//  Copyright © 2008-present: The CADET-Core Authors
//            Please see the AUTHORS.md file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Affero General Public
//  License v3.0 (or, at your option, any later version).
// =================================================================================

#include <catch.hpp>

//#include "Approx.hpp"
//#include "cadet/cadet.hpp"
//#define CADET_LOGGING_DISABLE
//#include "Logging.hpp"
#include "ColumnTests.hpp"


TEST_CASE("GRM SPLIT_COMPONENTS_DATA output fields", "[GRM],[FV],[Simulation],[Output],[CI]")
{
	cadet::test::column::testOutputSplitComponentsData("GENERAL_RATE_MODEL", "FV");
}

TEST_CASE("GRM IDAS META output fields", "[GRM],[FV],[Simulation],[Output],[CI]")
{
	cadet::test::column::testOutputIDASMetaData("GENERAL_RATE_MODEL", "FV");
}