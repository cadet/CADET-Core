function [result, solution] = breakthroughLangmuirSingle2D()
%BREAKTHROUGHLANGMUIRSINGLE2D Breakthrough curve with 1 component Langmuir model in 2D GRM
%
%   The model is taken from the publication:
%   Püttmann, A.; Schnittert, S.; Naumann, U.; von Lieres, E. (2013). 
%   Fast and accurate parameter sensitivities for the general rate model 
%   of column liquid chromatography. Computers & Chemical Engineering, 56, 46–57. 
%   doi:10.1016/j.compchemeng.2013.04.021
% 
%   It describes affinity chromatography of lysozyme on Cibacron Blue Sepharose CL-6B
%   at pH 7.2. The column is constantly fed with the same concentration in order to
%   produce a breakthrough curve.
%
%   The column is divided into 4 radial zones that have slightly different
%   flow rates. The output of the column is merged together using an outlet
%   unit operation.
%
%   See also PULSELINEARSINGLE, BREAKTHROUGHLANGMUIRSINGLE.

% Copyright: (C) 2008-2021 The CADET Authors
%            See the license note at the end of the file.


	% Step 1: Construct general rate model unit operation
	% ===================================================

	% General rate model unit operation
	mGrm = GeneralRateModel2D();

	% Discretization
	mGrm.nComponents = 1;
	mGrm.nCellsColumn = 16; % Attention: This is very low and only used for illustration (short runtime)
	mGrm.nCellsRadial = 4; % Number of radial zones
	mGrm.nCellsParticle = 4; % Attention: This is very low and only used for illustration (short runtime)
	mGrm.nBoundStates = ones(mGrm.nComponents, 1); % Number of bound states for each component

	% Initial conditions (empty column)
	mGrm.initialBulk = [0.0]; % [mol / m^3], also used for the particle mobile phase
	mGrm.initialSolid = [0.0]; % [mol / m^3]

	% Transport
	mGrm.dispersionColumn          = 5.75e-8; % [m^2 / s]
	mGrm.dispersionColumnRadial    = 1e-8; % [m^2 / s]
	mGrm.filmDiffusion             = [6.9e-6]; % [m/s]
	mGrm.diffusionParticle         = [6.07e-11]; % [m^2 / s]
	mGrm.diffusionParticleSurface  = [0.0]; % [m^2 / s]
	mGrm.interstitialVelocity      = 5.75e-4; % Forward flow

	% Geometry
	mGrm.columnLength        = 0.014; % [m]
	mGrm.columnRadius        = 0.01; % [m]
	mGrm.particleRadius      = 4.5e-5; % [m]
	mGrm.porosityColumn      = 0.37; % [-]
	mGrm.porosityParticle    = 0.75; % [-]
	
	% Adsorption
	mLangmuir = LangmuirBinding();
	mLangmuir.kineticBinding = true; % Kinetic binding
	mLangmuir.kA         = [1.14]; % Adsorption rate [m^3 / (mol * s)]
	mLangmuir.kD         = [0.002]; % Desorption rate [1 / s]
	mLangmuir.qMax       = [4.88]; % Capacity [mol / m^3]
	mGrm.bindingModel = mLangmuir;


	% Step 2: Construct inlet unit operation
	% ===================================================

	% Inlet unit operation
	mIn = PiecewiseCubicPolyInlet();
	mIn.nComponents = 1;
	
	% Reserve space: nSections x nComponents (a section can be thought of being a 
	% step in the process, see below)
	mIn.constant       = zeros(1, mGrm.nComponents);
	mIn.linear         = zeros(1, mGrm.nComponents);
	mIn.quadratic      = zeros(1, mGrm.nComponents);
	mIn.cubic          = zeros(1, mGrm.nComponents);

	% Section 1: Loading phase
	mIn.constant(1,1)  = 7.14e-3;  % [mol / m^3] component 1


	% Step 3: Construct outlet unit operation
	% ===================================================

	% Inlet unit operation
	mOut = OutletModel();
	mOut.nComponents = 1;
	
	
	% Step 4: Assemble system of unit operations
	% ===================================================

	% Construct ModelSystem and assign unit operations (order determines IDs)
	mSys = ModelSystem();
	mSys.models = [mIn, mGrm, mOut];

	% Calculate cross section areas of GRM
	radialBoundaries = linspace(0, mGrm.columnRadius, mGrm.nCellsRadial + 1);
	crossArea = pi * (radialBoundaries(2:end).^2 - radialBoundaries(1:end-1).^2);

	% Calculate volumetric flow rate required to achieve an interstitial velocity of 5.75e-4 m/s
	flow = 5.75e-4 .* crossArea .* mGrm.porosityColumn;

	% Make the flow rates slightly inhomogeneous (faster at the column walls)
	flow = flow(:) .* [0.95; 0.975; 1.0; 1.025];

	% Define valve configurations / unit operation connections
	% Valve configuration active on entering section 0
	mSys.connectionStartSection = [0];
	% Connect unit 0 port 0 with unit 1 port 0 (all components), flow rate flow(1);
	% Connect unit 0 port 1 with unit 1 port 1 (all components), flow rate flow(2) ...
	mSys.connections = {[0, 1, 0, 0, -1, -1, flow(1); ...
	                     0, 1, 0, 1, -1, -1, flow(2); ...
	                     0, 1, 0, 2, -1, -1, flow(3); ...
	                     0, 1, 0, 3, -1, -1, flow(4); ...
	                     1, 2, 0, 0, -1, -1, flow(1); ...
	                     1, 2, 1, 0, -1, -1, flow(2); ...
	                     1, 2, 2, 0, -1, -1, flow(3); ...
	                     1, 2, 3, 0, -1, -1, flow(4)]};


	% Step 5: Create simulator and configure it
	% ===================================================

	% Construct and configure simulator
	sim = Simulator.create();

	% sectionTimes holds the sections and sectionContinuity indicates whether
	% the transition between two adjacent sections is continuous
	sim.sectionTimes = [0.0 10000.0]; % [s]
	sim.sectionContinuity = [];

	% Hand model over to simulator
	sim.model = mSys;


	% Step 6: Run the model and plot the results
	% ===================================================

	% Run the model
	result = sim.run();

	% Extract solution into a matrix with time being the first column
	solution = [result.solution.time, squeeze(result.solution.outlet{2}), squeeze(result.solution.inlet{3})];

	% Plot the solution
	plot(solution(:, 1), solution(:, end), 'LineWidth', 2, 'LineStyle', '--');
	hold on;
	plot(solution(:, 1), solution(:, 2:end-1));
	hold off;
	legend('Port 0 (Inner)', 'Port 1', 'Port 2', 'Port 3 (Outer)', 'Location', 'northwest');
	grid on;
	xlabel('Time [s]');
	ylabel('Concentration [mM]');
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
