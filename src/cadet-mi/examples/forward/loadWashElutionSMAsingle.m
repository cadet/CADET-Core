function [result, solution] = loadWashElutionSMAsingle()
%LOADWASHELUTIONSMASINGLE Simple Load-Wash-Elution cycle with 4 component SMA model
%
%   Does exactly the same as loadWashElutionSMAsystem(), but accomplishes it in
%   a different way.
%
%   Instead of explicitly building a ModelSystem with general rate model and inlet
%   unit operation, a shortcut is used. The SingleGRM class simplifies the configuration
%   by combining general rate model (main unit operation model), inlet model, and 
%   ModelSystem. 
%
%   This is the preferred way for simple systems (just one main unit operation model).
%
%   See also LOADWASHELUTIONSMASYSTEM.

% Copyright: (C) 2008-2018 The CADET Authors
%            See the license note at the end of the file.


	% Step 1: Construct a system with general rate model 
	%         as main unit operation
	% ===================================================

	% General rate model unit operation
	mGrm = SingleGRM();

	% Discretization
	mGrm.nComponents = 4;
	mGrm.nCellsColumn = 16; % Attention: This is very low and only used for illustration (short runtime)
	mGrm.nCellsParticle = 4; % Attention: This is very low and only used for illustration (short runtime)
	mGrm.nBoundStates = ones(mGrm.nComponents, 1); % Number of bound states for each component

	% Initial conditions (equilibrated empty column)
	mGrm.initialBulk = [50.0 0.0 0.0 0.0]; % [mol / m^3], also used for the particle mobile phase
	mGrm.initialSolid = [1.2e3 0.0 0.0 0.0]; % [mol / m^3]
		
	% Transport
	mGrm.dispersionColumn          = 5.75e-8; % [m^2 / s]
	mGrm.filmDiffusion             = [6.9e-6 6.9e-6 6.9e-6 6.9e-6]; % [m/s]
	mGrm.diffusionParticle         = [7e-10 6.07e-11 6.07e-11 6.07e-11]; % [m^2 / s]
	mGrm.diffusionParticleSurface  = [0.0 0.0 0.0 0.0]; % [m^2 / s]
	mGrm.interstitialVelocity      = 5.75e-4; % [m/s]

	% Geometry
	mGrm.columnLength        = 0.014; % [m]
	mGrm.particleRadius      = 4.5e-5; % [m]
	mGrm.porosityColumn      = 0.37; % [-]
	mGrm.porosityParticle    = 0.75; % [-]
	
	% Adsorption
	mSma = StericMassActionBinding();
	mSma.kineticBinding = false; % Quasi-stationary binding
	mSma.lambda     = 1.2e3; % Ionic capacity [mol / m^3]
	mSma.kA         = [0.0 35.5 1.59 7.7]; % Adsorption rate [(m^3 / mol)^nu / s]
	mSma.kD         = [0.0 1000 1000 1000]; % Desorption rate [(m^3 / mol)^nu / s]
	mSma.nu         = [0.0 4.7 5.29 3.7]; % Characteristic charge [-]
	mSma.sigma      = [0.0 11.83 10.6 10.0]; % Steric factor [-]
	mGrm.bindingModel = mSma;
	
	% Specify inlet profile

	% Reserve space: nSections x nComponents (a section can be thought of being a 
	% step in the process, see below)
	mGrm.constant       = zeros(3, mGrm.nComponents);
	mGrm.linear         = zeros(3, mGrm.nComponents);
	mGrm.quadratic      = zeros(3, mGrm.nComponents);
	mGrm.cubic          = zeros(3, mGrm.nComponents);

	% Section 1: Loading phase
	mGrm.constant(1,1)  = 50.0;  % [mol / m^3] component 1
	mGrm.constant(1,2)  = 1.0;   % [mol / m^3] component 2
	mGrm.constant(1,3)  = 1.0;   % [mol / m^3] component 3
	mGrm.constant(1,4)  = 1.0;   % [mol / m^3] component 4

	% Section 2: Washing phase (no protein feed)
	mGrm.constant(2,1)  = 50.0;  % [mol / m^3] component 1

	% Section 3: Elution phase (linear salt gradient with step at the beginning)
	mGrm.constant(3,1)  = 100;  % [mol / m^3] component 1
	mGrm.linear(3,1)  = 0.2;  % [mol / (m^3 * s)] component 1


	% Step 2: Create simulator and configure it
	% ===================================================

	% Construct and configure simulator
	sim = Simulator.create();
	sim.solutionTimes = linspace(0, 1500, 1001); % [s], time points at which solution is computed

	% sectionTimes holds the sections and sectionContinuity indicates whether
	% the transition between two adjacent sections is continuous
	sim.sectionTimes = [0.0 10.0 90.0 1500.0]; % [s]
	sim.sectionContinuity = false(2,1);

	% Hand model over to simulator	
	sim.model = mGrm;


	% Step 3: Run the model and plot the results
	% ===================================================

	% Run the model
	result = sim.run();

	% Extract solution into a matrix with time being the first column
	% Note that we need to extract the outlet of the first unit operation,
	% which is the general rate model (main unit operation is always first
	% in the SingleXYZ models)
	solution = [result.solution.time, result.solution.outlet{1}];

	% Plot the solution
	plot(solution(:, 1), solution(:, 3:end));
	legend('Lysozyme', 'Cytochrome', 'Ribonuclease');
	grid on;
	xlabel('Time [s]');
	ylabel('Concentration [mM]');
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
