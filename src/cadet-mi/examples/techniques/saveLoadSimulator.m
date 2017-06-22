function saveLoadSimulator()
%SAVELOADSIMULATOR Demonstrates how to save and load a model to a MAT-file
%
%   A simulator is constructed with a configured model and then saved to a MAT-file.
%   It is then loaded again and used for simulation.
%
%   See also SAVELOADMODEL.

% Copyright: (C) 2008-2017 The CADET Authors
%            See the license note at the end of the file.

	% Create the simulator with the model
	sim = createSimulator();
	
	% Save the simulator to a file
	fileName = [tempname, '.mat'];
	save(fileName, 'sim');
	
	% Load the simulator from the file
	loadedSim = load(fileName, 'sim');
	loadedSim = loadedSim.sim;

	% Run both simulators
	res = sim.run();
	res2 = loadedSim.run();

	% Plot the difference of the two solutions
	solution = [res.solution.time, res.solution.outlet{1}];
	solution2 = [res2.solution.time, res2.solution.outlet{1}];
	plot(solution(:, 1), abs(solution(:, 3:end) - solution2(:, 3:end)));
	legend('Lysozyme', 'Cytochrome', 'Ribonuclease');
	grid on;
	xlabel('Time [s]');
	ylabel('Concentration [mM]');
	title('abs(Simulator - Loaded Simulator)');
end

function sim = createSimulator()
	% General rate model
	mGrm = SingleGRM();

	mGrm.nComponents = 4;
	mGrm.nCellsColumn = 16; % Attention: This is very low and only used for illustration (short runtime)
	mGrm.nCellsParticle = 4; % Attention: This is very low and only used for illustration (short runtime)
	mGrm.nBoundStates = ones(mGrm.nComponents, 1);

	% Initial conditions
	mGrm.initialBulk = [50.0 0.0 0.0 0.0];
	mGrm.initialSolid = [1.2e3 0.0 0.0 0.0];
		
	% Transport
	mGrm.dispersionColumn          = 5.75e-8;
	mGrm.filmDiffusion             = [6.9e-6 6.9e-6 6.9e-6 6.9e-6];
	mGrm.diffusionParticle         = [7e-10 6.07e-11 6.07e-11 6.07e-11];
	mGrm.diffusionParticleSurface  = [0.0 0.0 0.0 0.0];
	mGrm.interstitialVelocity      = 5.75e-4;

	% Geometry
	mGrm.columnLength        = 0.014;
	mGrm.particleRadius      = 4.5e-5;
	mGrm.porosityColumn      = 0.37;
	mGrm.porosityParticle    = 0.75;
	
	% Adsorption
	mSma = StericMassActionBinding();
	mSma.kineticBinding = false;
	mSma.lambda     = 1.2e3;
	mSma.kA         = [0.0 35.5 1.59 7.7];
	mSma.kD         = [0.0 1000 1000 1000];
	mSma.nu         = [0.0 4.7 5.29 3.7];
	mSma.sigma      = [0.0 11.83 10.6 10.0];
	mGrm.bindingModel = mSma;

	mGrm.constant       = zeros(3, mGrm.nComponents);
	mGrm.linear         = zeros(3, mGrm.nComponents);
	mGrm.quadratic      = zeros(3, mGrm.nComponents);
	mGrm.cubic          = zeros(3, mGrm.nComponents);

	% Sec 1
	mGrm.constant(1,1)  = 50.0;  % component 1
	mGrm.constant(1,2)  = 1.0;   % component 2
	mGrm.constant(1,3)  = 1.0;   % component 3
	mGrm.constant(1,4)  = 1.0;   % component 4

	% Sec 2
	mGrm.constant(2,1)  = 50.0;  % component 1

	% Sec 3
	mGrm.constant(3,1)  = 100;  % component 1
	mGrm.linear(3,1)  = 0.2;  % component 1

	% Create simulator
	sim = Simulator.create();
	sim.sectionTimes = [0.0 10.0 90.0 1500.0];
	sim.sectionContinuity = false(2,1);
	sim.solutionTimes = linspace(0, 1500, 1001);

	% Assign model
	sim.model = mGrm;
end

% =============================================================================
%  CADET - The Chromatography Analysis and Design Toolkit
%  
%  Copyright (C) 2008-2017: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
