// =============================================================================
//  CADET
//  
//  Copyright © 2008-2023: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include <catch.hpp>

#include "ColumnTests.hpp"
#include "ReactionModelTests.hpp"
#include "Utils.hpp"
#include "JsonTestModels.hpp"

TEST_CASE("LRMP2D inlet DOF Jacobian", "[LRMP2D],[DG],[DG2D],[UnitOp],[Jacobian],[Inlet],[CI]") // @todo seems to not work anymore
{
	cadet::test::column::testInletDofJacobian("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG");
}

TEST_CASE("LRMP2D transport Jacobian", "[LRMP2DtestHere],[DG],[DG2D],[UnitOp],[Jacobian],[CI]")
{
	const std::string relModelFilePath = std::string("/data/lrmp2d_debug.json");
	cadet::JsonParameterProvider jpp = cadet::test::column::getReferenceFile(relModelFilePath);

	// get the number of radial ports
	jpp.pushScope("model");
	const int nUnits = jpp.getInt("NUNITS"); // there is one column and (nUnits-1)/2 inlet and outlet units
	const int columnIdx = (nUnits - 1) / 2; // assumes the inlet units come first, followed by the column unit
	const int nRad = columnIdx; // number of radial points or ports

	const std::string unitID = std::format("{:03}", columnIdx);

	// we need to set flowRates for 2D models for the JacobianAD test, since velocity is only set with a call to the setFlowRate function,
	// which in turn is only called if we specify flow rates. We get the flow rates from the connections matrix
	jpp.pushScope("connections");
	jpp.pushScope("switch_000");

	const std::vector<double> connections = jpp.getDoubleArray("CONNECTIONS");
	std::vector<cadet::active> flowRate;
	for (int i = 1; i <= nRad; i++)
		flowRate.push_back(connections[i * 7 - 1]);

	jpp.popScope();
	jpp.popScope();


	jpp.pushScope("unit_" + unitID);


	// Success | zElem | rElem
	// yes     |  1-7  | 1-4 
	// no      |  1-5  | 1-7 
	// yes     |  1-4  | 1-7 
	// yes     |  4-5  | 1-4 
	// yes     |  4-5  | 4-5 
	// no      |  4-5  | 5-7 
	// yes     |  4-5  | 5-6 
	// no      |  4-5  | 6-7 
	// no: *1  |  1-7  | 1-6
	// no      |  1-7  | 1-5

	// *1 : error in prepareAD vectors in setADvalue


	for (int zElem = 1; zElem < 5; zElem++) 
	{
		for (int rElem = 1; rElem < 5; rElem++)
		{

			jpp.pushScope("discretization");
			jpp.set("AX_NELEM", zElem);
			jpp.set("RAD_NELEM", rElem);
			jpp.popScope();

			cadet::test::column::testJacobianAD(jpp, 1e10, &flowRate[0]); // @todo figure out why FD Jacobian pattern comparison doesnt work but AD Jacobian comparison does

		}
	}
}
