function makeInletProfileFromResults()
%MAKEINLETPROFILEFROMRESULTS Turns the results of a simulation into an inlet profile and runs the simulation again
%
%   It is shown how the results of a simulation can be used to construct an inlet profile.
%   This inlet profile is then used to run the simulation again, which corresponds to
%   recycling (with some time lag).
%
%   See also LOADWASHELUTIONSMASINGLE.

% Copyright: (C) 2008-2021 The CADET Authors
%            See the license note at the end of the file.

	% Set up a simulation (GRM)
	sim = createSimulator();

	% Run the simulation and extract the solution for plotting it later
	res1 = sim.run();
	sol1 = [res1.solution.time, squeeze(res1.solution.outlet{1})];

	% Convert solution into inlet profile
	sim.model.inlet = PiecewiseCubicPolyProfile.fromResult(res1);
	
	% The constructed inlet profile is a spline (piecewise polynomial)
	% that has its breaks at the time points of the solution. The inlet
	% profile requires that the section times of the simulator match its
	% breaks. Thus, we have to set the sectionTimes to the solution time
	% points. 
	% However, we can assume the profile to be continuous. By handing
	% this information over to the simulator via sectionContinuity, we
	% can exploit the continuity and speed up the simulation.
	sim.sectionTimes = res1.solution.time;
	sim.sectionContinuity = true(length(res1.solution.time) - 2, 1);
	% Note that continuity is specified for each section transition. For
	% N time points, there are N-1 sections and, thus, N-2 section transitions.

	% Run the simulation again with the new inlet profile and extract the solution
	res2 = sim.run();
	sol2 = [res2.solution.time, squeeze(res2.solution.outlet{1})];

	% Plot both solutions side by side
	figure;
	subplot(1, 2, 1);
	plot(sol1(:,1), sol1(:,3:end));
	legend('Lysozyme', 'Cytochrome', 'Ribonuclease');
	xlabel('Time [s]');
	ylabel('Concentration [mM]');
	title('First');
	grid on;
	
	subplot(1, 2, 2);
	plot(sol2(:,1), sol2(:,3:end));
	legend('Lysozyme', 'Cytochrome', 'Ribonuclease');
	xlabel('Time [s]');
	ylabel('Concentration [mM]');
	title('Recycled');
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
	sim.solutionTimes = linspace(0, 1500, 301);
	
	% Assign model
	sim.model = mGrm;
end

% =============================================================================
%  CADET
%  
%  Copyright (C) 2008-2021: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
