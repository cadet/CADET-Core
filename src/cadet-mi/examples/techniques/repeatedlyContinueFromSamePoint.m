function repeatedlyContinueFromSamePoint()
%REPEATEDLYCONTINUEFROMSAMEPOINT Shows how to simulate up to a certain point and resume the simulation from there repeatedly with different parameters
%
%   In this example, the load and wash steps of the load-wash-elution example
%   are performed and the simulation is paused at this point. The full system
%   state is saved and the simulation is continued from this point repeatedly
%   with different gradient slopes in the elution step.
%
%   This setup is especially useful in process analysis and process optimization,
%   for instance, where the start phase (loading, washing) is always the same
%   and variations happen afterwards (varying salt gradients).
%
%   See also PARAMETERIZEDSIMULATIONWITHOUTSENSITIVITIES, RESUMECONTINUOUSTIME,
%      RESUMERESETTIME, LOADWASHELUTIONSMASINGLE

% Copyright: (C) 2008-2018 The CADET Authors
%            See the license note at the end of the file.

	% Create simulator
	sim = createSimulator();

	% Compute until elution phase and return the full system state at the end point
	sim.solutionTimes = linspace(0, 90, 91);
	sim.returnLastState = true;
	res1 = sim.run();

	% Extract solution
	sol1 = [res1.solution.time, res1.solution.outlet{1}];

	% Set initial state to last state and disable output of last system state
	sim.model.initStateY = res1.solution.lastState;
	sim.model.initStateYdot = res1.solution.lastStateDot;
	sim.returnLastState = false;

	% Prepare for elution phase:
	%   - Set time from 90 to 1500
	%   - Set gradient as inlet profile
	sim.solutionTimes = linspace(90, 1500, 1500 - 90 +1);
	sim.sectionTimes = [90.0, 1500.0];
	sim.sectionContinuity = [];

	model = sim.model;
	model.constant       = [100, 0, 0, 0]; % Step at elution phase
	model.linear         = [0.2, 0, 0, 0]; % Gradient in elution phase
	model.quadratic      = zeros(1, 4);
	model.cubic          = zeros(1, 4);

	% Add gradient slope as parameter (without sensitivity)
	params = cell(1, 1);
	params{1} = makeSensitivity([1], {'LIN_COEFF'}, [0], [-1], [-1], [0]);
	sim.setParameters(params, [false]);

	% Simulate elution for all slopes
	slopes = [0.1, 0.15, 0.2, 0.25, 0.3];
	shapes = cell(length(slopes), 1);

	for i = 1:length(slopes)
		% Skip validation of input data on subsequent runs
		res2 = sim.runWithParameters(slopes(i), true);

		% Extract solution and glue the pieces together
		shapes{i} = [sol1; ...
		             res2.solution.time, res2.solution.outlet{1}];
	end

	% Plot the results
	figure;
	for i = 1:length(slopes)
		sol = shapes{i};

		subplot(2, 3, i);
		hdAx = plotyy(sol(:, 1), sol(:, 3:end), sol(:, 1), sol(:, 2));
		legend('Lysozyme', 'Cytochrome', 'Ribonuclease', 'Salt');
		xlabel('Time [s]');
		ylabel('Concentration [mM]');
		title(sprintf('Slope %g', slopes(i)));
		grid on;

		ylim(hdAx(1), [0, 0.06]);
		ylim(hdAx(2), [0, 550]);
	end

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

	mGrm.constant       = zeros(2, mGrm.nComponents);
	mGrm.linear         = zeros(2, mGrm.nComponents);
	mGrm.quadratic      = zeros(2, mGrm.nComponents);
	mGrm.cubic          = zeros(2, mGrm.nComponents);

	% Sec 1
	mGrm.constant(1,1)  = 50.0;  % component 1
	mGrm.constant(1,2)  = 1.0;   % component 2
	mGrm.constant(1,3)  = 1.0;   % component 3
	mGrm.constant(1,4)  = 1.0;   % component 4

	% Sec 2
	mGrm.constant(2,1)  = 50.0;  % component 1

	% Create and configure simulator
	sim = Simulator.create();
	sim.sectionTimes = [0.0 10.0 90.0];
	sim.sectionContinuity = false;
	
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
