function parameterizedSimulationWithoutSensitivities()
%PARAMETERIZEDSIMULATIONWITHOUTSENSITIVITIES Shows how to parameterize a simulation without computing sensitivities
%
%   A simulation can be parameterized by specifying parameters to be variables.
%   Using MEXSIMULATOR.RUNWITHPARAMETERS, a simulation can easily be run again with
%   different parameter values, which is helpful for parameter studies or 
%   characterizing operation windows. Since no sensitivities are computed, the
%   simulations are fast.
%
%   See also MAKESENSITIVITY, MEXSIMULATOR.RUNWITHPARAMETERS, MEXSIMULATOR.SETPARAMETERS,
%      PARAMETERIZEDSIMULATIONWITHSENSITIVITIES, PARAMETERIZEDSIMULATIONMIXED

% Copyright: (C) 2008-2018 The CADET Authors
%            See the license note at the end of the file.

	% Create simulator
	sim = createSimulator();

	% Add three parameters (given as cell array)
	params = cell(3, 1);

	% Each parameter is identified by 
	%   - its unit operation id (0-based index or -1 if independent),
	%   - name (according to CADET file format specification),
	%   - component index (0-based index or -1 if independent),
	%   - reaction index (0-based index or -1 if independent),
	%   - bound phase index (0-based index or -1 if independent), and
	%   - time integrator section index (0-based index or -1 if independent).

	% Parameter 1: SMA_KA of component 2 (Lysozyme, component 1 if read as 0-based) in GRM (unit operation 0)
	params{1} = makeSensitivity([0], {'SMA_KA'}, [1], [-1], [0], [-1]);

	% Parameter 2: COL_LENGTH (independent of component) in GRM (unit operation 0)
	params{2} = makeSensitivity([0], {'COL_LENGTH'}, [-1], [-1], [-1], [-1]);

	% Parameter 3: CONST_COEFF of component 2 (Lysozyme, component 1 if read as 0-based)
	%              in inlet unit operation (id 1) and time section 0. This is the loading
	%              concentration of Lysozyme.
	params{3} = makeSensitivity([1], {'CONST_COEFF'}, [1], [-1], [-1], [0]);

	% Set parameters in the Simulator and do not compute sensitivities
	sim.setParameters(params, false(3, 1));

	% Run first simulation with current parameters set in the model
	res1 = sim.runWithParameters([]);
	sol1 = [res1.solution.time, res1.solution.outlet{1}];

	% Run second simulation with the parameters given in the vector (i.e., lower adsoprtion rate, 
	% longer column, and higher loading), skip validation of input data on second run
	res2 = sim.runWithParameters([30, 0.016, 1.1], true);
	sol2 = [res2.solution.time, res2.solution.outlet{1}];

	% Note that the model parameters have changed to the given values:
	assert(sim.model(1).bindingModel.kA(2) == 30);
	assert(sim.model(1).columnLength == 0.016);
	assert(sim.model(1).constant(1, 2) == 1.1);
	
	% Plot the two chromatograms
	figure;
	subplot(2, 1, 1);
	plot(sol1(:, 1), sol1(:,3:end));
	legend('Lysozyme', 'Cytochrome', 'Ribonuclease');
	xlabel('Time [s]');
	ylabel('Concentration [mM]');
	title('KA = 35.5, COL_LENGTH = 0.014, LOAD = 1');
	grid on;

	subplot(2, 1, 2);
	plot(sol2(:, 1), sol2(:,3:end));
	legend('Lysozyme', 'Cytochrome', 'Ribonuclease');
	xlabel('Time [s]');
	ylabel('Concentration [mM]');
	title('KA = 30, COL_LENGTH = 0.016, LOAD = 1.1');
	grid on;
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

	% Create and configure simulator
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
%  Copyright (C) 2008-2018: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
