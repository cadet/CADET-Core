// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2020: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include <catch.hpp>

#include "graph/GraphAlgos.hpp"

TEST_CASE("Linear graph all ports all comps no cycles", "[Graph]")
{
	const std::vector<int> connections = {
		0, 1, -1, -1, -1, -1,
		1, 2, -1, -1, -1, -1,
		2, 3, -1, -1, -1, -1
	};
	cadet::util::SlicedVector<int> adjList = cadet::graph::adjacencyListFromConnectionList(connections.data(), 4, connections.size() / 6);
	std::vector<int> topoOrder;
	const bool cycle = cadet::graph::topologicalSort(adjList, topoOrder);

	REQUIRE(!cycle);
}

TEST_CASE("Linear graph specific ports all comps no cycles", "[Graph]")
{
	const std::vector<int> connections = {
		0, 1, 0, 0, -1, -1,
		1, 2, 0, 0, -1, -1,
		2, 3, 0, 0, -1, -1
	};
	cadet::util::SlicedVector<int> adjList = cadet::graph::adjacencyListFromConnectionList(connections.data(), 4, connections.size() / 6);
	std::vector<int> topoOrder;
	const bool cycle = cadet::graph::topologicalSort(adjList, topoOrder);

	REQUIRE(!cycle);
}

TEST_CASE("Linear graph all ports specific comps no cycles", "[Graph]")
{
	const std::vector<int> connections = {
		0, 1, -1, -1, 0, 0,
		0, 1, -1, -1, 1, 1,
		0, 1, -1, -1, 2, 2,
		1, 2, -1, -1, 0, 0,
		1, 2, -1, -1, 1, 1,
		1, 2, -1, -1, 2, 2,
		2, 3, -1, -1, 0, 0,
		2, 3, -1, -1, 1, 1,
		2, 3, -1, -1, 2, 2
	};
	cadet::util::SlicedVector<int> adjList = cadet::graph::adjacencyListFromConnectionList(connections.data(), 4, connections.size() / 6);
	std::vector<int> topoOrder;
	const bool cycle = cadet::graph::topologicalSort(adjList, topoOrder);

	REQUIRE(!cycle);
}

TEST_CASE("Linear graph specific ports and comps no cycles", "[Graph]")
{
	const std::vector<int> connections = {
		0, 1, 0, 0, 0, 0,
		0, 1, 0, 0, 1, 1,
		0, 1, 0, 0, 2, 2,
		1, 2, 0, 0, 0, 0,
		1, 2, 0, 0, 1, 1,
		1, 2, 0, 0, 2, 2,
		2, 3, 0, 0, 0, 0,
		2, 3, 0, 0, 1, 1,
		2, 3, 0, 0, 2, 2
	};
	cadet::util::SlicedVector<int> adjList = cadet::graph::adjacencyListFromConnectionList(connections.data(), 4, connections.size() / 6);
	std::vector<int> topoOrder;
	const bool cycle = cadet::graph::topologicalSort(adjList, topoOrder);

	REQUIRE(!cycle);
}

TEST_CASE("Two ports linear graph specific ports all comps no cycles", "[Graph]")
{
	/*
	         1     1     1
	      /-->--O-->--O-->--\
	2 -->O                   O--> 2
	      \-->--O-->--O-->--/
	         1     1     1
	*/

	const std::vector<int> connections = {
		0, 1, 0, 0, -1, -1,
		0, 1, 0, 1, -1, -1,
		1, 2, 0, 1, -1, -1,
		1, 2, 1, 0, -1, -1,
		2, 3, 0, 0, -1, -1,
		2, 3, 1, 0, -1, -1
	};
	cadet::util::SlicedVector<int> adjList = cadet::graph::adjacencyListFromConnectionList(connections.data(), 4, connections.size() / 6);
	std::vector<int> topoOrder;
	const bool cycle = cadet::graph::topologicalSort(adjList, topoOrder);

	REQUIRE(!cycle);
}

TEST_CASE("Two ports linear graph specific ports and comps no cycles", "[Graph]")
{
	/*
	         1     1     1
	      /-->--O-->--O-->--\
	2 -->O                   O--> 2
	      \-->--O-->--O-->--/
	         1     1     1
	*/

	const std::vector<int> connections = {
		0, 1, 0, 0, 0, 0,
		0, 1, 0, 0, 1, 1,
		0, 1, 0, 0, 2, 2,
		0, 1, 0, 1, 0, 0,
		0, 1, 0, 1, 1, 1,
		0, 1, 0, 1, 2, 2,
		1, 2, 0, 1, 0, 0,
		1, 2, 0, 1, 1, 1,
		1, 2, 0, 1, 2, 2,
		1, 2, 1, 0, 0, 0,
		1, 2, 1, 0, 1, 1,
		1, 2, 1, 0, 2, 2,
		2, 3, 0, 0, 0, 0,
		2, 3, 0, 0, 1, 1,
		2, 3, 0, 0, 2, 2,
		2, 3, 1, 0, 0, 0,
		2, 3, 1, 0, 1, 1,
		2, 3, 1, 0, 2, 2
	};
	cadet::util::SlicedVector<int> adjList = cadet::graph::adjacencyListFromConnectionList(connections.data(), 4, connections.size() / 6);
	std::vector<int> topoOrder;
	const bool cycle = cadet::graph::topologicalSort(adjList, topoOrder);

	REQUIRE(!cycle);
}

TEST_CASE("N to 1 to N graph specific ports all comps no cycles", "[Graph]")
{
	/*
	    O--\     /--O
	        --O--
	    O--/     \--O
	*/

	const std::vector<int> connections = {
		0, 2, 0, 0, -1, -1,
		1, 2, 0, 0, -1, -1,
		2, 3, 0, 0, -1, -1,
		2, 4, 0, 0, -1, -1
	};
	cadet::util::SlicedVector<int> adjList = cadet::graph::adjacencyListFromConnectionList(connections.data(), 5, connections.size() / 6);
	std::vector<int> topoOrder;
	const bool cycle = cadet::graph::topologicalSort(adjList, topoOrder);

	REQUIRE(!cycle);
}

TEST_CASE("N to 1 to N graph all ports specific comps no cycles", "[Graph]")
{
	/*
	    O--\     /--O
	        --O--
	    O--/     \--O
	*/

	const std::vector<int> connections = {
		0, 2, 0, 0, -1, -1,
		1, 2, 0, 1, -1, -1,
		2, 3, 1, 0, -1, -1,
		2, 4, 0, 0, -1, -1
	};
	cadet::util::SlicedVector<int> adjList = cadet::graph::adjacencyListFromConnectionList(connections.data(), 5, connections.size() / 6);
	std::vector<int> topoOrder;
	const bool cycle = cadet::graph::topologicalSort(adjList, topoOrder);

	REQUIRE(!cycle);
}

TEST_CASE("Linear graph specific ports twisted comps no cycles", "[Graph]")
{
	const std::vector<int> connections = {
		0, 1, 0, 0, 0, 0,
		0, 1, 0, 0, 1, 1,
		0, 1, 0, 0, 2, 2,
		1, 2, 0, 0, 0, 2,
		1, 2, 0, 0, 1, 1,
		1, 2, 0, 0, 2, 0,
		2, 3, 0, 0, 0, 0,
		2, 3, 0, 0, 1, 1,
		2, 3, 0, 0, 2, 2
	};
	cadet::util::SlicedVector<int> adjList = cadet::graph::adjacencyListFromConnectionList(connections.data(), 4, connections.size() / 6);
	std::vector<int> topoOrder;
	const bool cycle = cadet::graph::topologicalSort(adjList, topoOrder);

	REQUIRE(!cycle);
}

TEST_CASE("Two cycles graph mixed ports all comps", "[Graph]")
{
	/*
	    ___________________
	    |                 |
	0---2--\     /--5---  |
	        --4--      |  |
	1---3--/     \--6-----
	    |______________|
	*/

	const std::vector<int> connections = {
		0, 2,  0,  0, -1, -1,
		1, 3,  0,  0, -1, -1,
		2, 4, -1, -1, -1, -1,
		3, 4, -1, -1, -1, -1,
		4, 5,  0,  0, -1, -1,
		4, 6, -1, -1, -1, -1,
		5, 3,  0,  1, -1, -1,
		6, 2,  0,  1, -1, -1
	};
	cadet::util::SlicedVector<int> adjList = cadet::graph::adjacencyListFromConnectionList(connections.data(), 7, connections.size() / 6);
	std::vector<int> topoOrder;
	const bool cycle = cadet::graph::topologicalSort(adjList, topoOrder);

	REQUIRE(cycle);
}

TEST_CASE("N to 1 graph specific ports mixed comps no cycles", "[Graph]")
{
	/*
	    O--\     
	        --O--O
	    O--/     
	*/

	const std::vector<int> connections = {
		0, 2, 0, 0,  0,  0,
		0, 2, 0, 0,  1,  1,
		1, 2, 0, 0,  0,  2,
		1, 2, 0, 0,  0,  3,
		2, 3, 0, 0, -1, -1
	};
	cadet::util::SlicedVector<int> adjList = cadet::graph::adjacencyListFromConnectionList(connections.data(), 4, connections.size() / 6);
	std::vector<int> topoOrder;
	const bool cycle = cadet::graph::topologicalSort(adjList, topoOrder);

	REQUIRE(!cycle);
}

