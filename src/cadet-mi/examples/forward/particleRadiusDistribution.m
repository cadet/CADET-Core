function particleRadiusDistribution()
%PARTICLERADIUSDISTRIBUTION Demonstrates multiple particle types with particle radius distribution
%
%   The PULSELINEARSINGLE example is extended with a particle radius distribution.
%   Three particle radii are considered while keeping all other parameters the
%   same.
%
%   See also PULSELINEARSINGLE, MULTIPARTICLETYPES.

% Copyright: (C) 2008-2024 The CADET Authors
%            See the license note at the end of the file.


	% Step 1: Construct a system with general rate model 
	%         as main unit operation
	% ===================================================

	% General rate model unit operation
	mGrm = SingleGRM();

	% Discretization
	mGrm.nComponents = 1;
	mGrm.nParticleTypes = 3;
	mGrm.nCellsColumn = 16; % Attention: This is very low and only used for illustration (short runtime)
	mGrm.nCellsParticle = [4, 6, 5]; % Attention: This is very low and only used for illustration (short runtime)
	mGrm.nBoundStates = ones(mGrm.nComponents * 3, 1); % Number of bound states for each component and particle type

	% Initial conditions
	mGrm.initialBulk = [0.0]; % [mol / m^3], also used for the particle mobile phase
	mGrm.initialSolid = [0.0, 0.0, 0.0]; % For each bound state in each particle type [mol / m^3]
		
	% Transport
	mGrm.dispersionColumn          = 0.002 / (100*100*60); % [m^2 / s]
	mGrm.filmDiffusion             = [0.01 / (100 * 60), ...  % First particle type [m/s]
	                                  0.01 / (100 * 60), ...  % Second particle type [m/s]
	                                  0.01 / (100 * 60)];     % Third particle type [m/s]
	mGrm.diffusionParticle         = [3.003e-6, ...  % First particle type [m^2 / s]
	                                  3.003e-6, ...  % Second particle type [m^2 / s]
	                                  3.003e-6];     % Third particle type [m^2 / s]
	mGrm.diffusionParticleSurface  = [0.0, ...  % First particle type [m^2 / s]
	                                  0.0, ...  % Second particle type [m^2 / s]
	                                  0.0];     % Third particle type [m^2 / s]
	mGrm.interstitialVelocity      = 0.5 / (100*60); % [m/s]

	% Geometry
	mGrm.columnLength        = 0.017; % [m]
	mGrm.particleRadius      = [4e-5, 2e-5, 9e-5]; % For each particle type in [m]
	mGrm.particleCoreRadius  = [0.0, 0.0, 0.0]; % For each particle in [m]
	mGrm.porosityColumn      = 0.4; % [-]
	mGrm.porosityParticle    = [0.333, 0.333, 0.333]; % For each particle type [-]

	% Volume fractions of the particle types (have to sum to 1.0)
	mGrm.particleTypeVolumeFractions = [0.3, 0.5, 0.2];

	% Adsorption
	mLinear = LinearBinding();
	mLinear.kineticBinding = true; % Kinetic binding
	mLinear.kA         = [2.5]; % Adsorption rate [m^3 / (mol * s)]
	mLinear.kD         = [1]; % Desorption rate [1 / s]
	mGrm.bindingModel = mLinear;
	
	% Specify inlet profile

	% Reserve space: nSections x nComponents (a section can be thought of being a 
	% step in the process, see below)
	mGrm.constant       = zeros(2, mGrm.nComponents);
	mGrm.linear         = zeros(2, mGrm.nComponents);
	mGrm.quadratic      = zeros(2, mGrm.nComponents);
	mGrm.cubic          = zeros(2, mGrm.nComponents);

	% Section 1: Pulse
	mGrm.constant(1,1)  = 1;  % [mol / m^3] component 1


	% Step 2: Create simulator and configure it
	% ===================================================

	% Construct and configure simulator
	sim = Simulator.create();
	sim.solutionTimes = linspace(0, 60 * 60, 3001); % [s], time points at which solution is computed

	% sectionTimes holds the sections and sectionContinuity indicates whether
	% the transition between two adjacent sections is continuous
	sim.sectionTimes = [0.0 20*60 60*60]; % [s]
	sim.sectionContinuity = [false];

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
	solution = [result.solution.time, squeeze(result.solution.outlet{1})];

	% Plot the solution
	plot(solution(:, 1), solution(:, 2));
	grid on;
	xlabel('Time [s]');
	ylabel('Concentration [mM]');
end

% =============================================================================
%  CADET
%  
%  Copyright (C) 2008-2024: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
