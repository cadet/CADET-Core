function sectionDependentParameters()
%SECTIONDEPENDENTPARAMETERS Demonstrates how to set up section dependent parameters
%
%   Some (mostly transport related) parameters can be section dependent. This may
%   be helpful if parts of a process are operated under different conditions.
%   In this example, the film diffusion coefficient and the interstitial velocity
%   is made section dependent and varied in different sections of the standard
%   load-wash-elution example.
%
%   See also LOADWASHELUTIONSMASINGLE

% Copyright: (C) 2008-2022 The CADET Authors
%            See the license note at the end of the file.


	% General rate model unit operation
	mGrm = SingleGRM();

	% Discretization
	mGrm.nComponents = 4;
	mGrm.nCellsColumn = 16; % Attention: This is very low and only used for illustration (short runtime)
	mGrm.nCellsParticle = 4; % Attention: This is very low and only used for illustration (short runtime)
	mGrm.nBoundStates = ones(mGrm.nComponents, 1);

	% Initial conditions (equilibrated empty column)
	mGrm.initialBulk = [50.0 0.0 0.0 0.0];
	mGrm.initialSolid = [1.2e3 0.0 0.0 0.0];
		
	% Transport

	% Dispersion does not depend on section here (but it could)
	mGrm.dispersionColumn          = 5.75e-8;

	% Film diffusion depends on section (section-component major order)
	mGrm.filmDiffusion             = [6.9e-6 6.9e-6 6.9e-6 6.9e-6, ... % For components in section 1
	                                  6.9e-7 6.9e-7 6.9e-7 6.9e-7, ... % For components in section 2
	                                  6.9e-5 6.9e-5 6.9e-5 6.9e-5, ... % For components in section 3
	                                  6.9e-5 6.9e-5 6.9e-5 6.9e-5];    % For components in section 4

	% Particle diffusion does not depend on section here (but it could)
	mGrm.diffusionParticle         = [7e-10 6.07e-11 6.07e-11 6.07e-11];

	% Surface diffusion does not depend on section here (but it could)
	mGrm.diffusionParticleSurface  = [0.0 0.0 0.0 0.0];

	% Interstitial velocity depends on the section (increased velocity in last section)
	mGrm.interstitialVelocity      = [5.75e-4 5.75e-4 5.75e-4 9e-4];

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
	
	% Specify inlet profile

	% Reserve space
	mGrm.constant       = zeros(4, mGrm.nComponents);
	mGrm.linear         = zeros(4, mGrm.nComponents);
	mGrm.quadratic      = zeros(4, mGrm.nComponents);
	mGrm.cubic          = zeros(4, mGrm.nComponents);

	% Section 1: Loading phase
	mGrm.constant(1,1)  = 50.0;
	mGrm.constant(1,2)  = 1.0;
	mGrm.constant(1,3)  = 1.0;
	mGrm.constant(1,4)  = 1.0;

	% Section 2: Washing phase (no protein feed)
	mGrm.constant(2,1)  = 50.0;

	% Section 3: Elution phase, first half
	mGrm.constant(3,1)  = 100;
	mGrm.linear(3,1)  = 0.2;

	% Section 4: Elution phase, second half
	mGrm.constant(4,1)  = 100 + 0.2 * (750 - 90);
	mGrm.linear(4,1)  = 0.2;

	% Construct and configure simulator
	sim = Simulator.create();
	sim.solutionTimes = linspace(0, 1500, 1001);
	sim.sectionTimes = [0.0 10.0 90.0 750.0 1500.0];
	sim.sectionContinuity = false(3, 1);

	% Hand model over to simulator	
	sim.model = mGrm;

	% Run the model and extract solution
	result = sim.run();
	solution = [result.solution.time, squeeze(result.solution.outlet{1})];

	% Plot the solution
	subplot(1, 2, 1);
	plot(solution(:, 1), solution(:, 2));
	legend('Salt');
	grid on;
	xlabel('Time [s]');
	ylabel('Concentration [mM]');

	subplot(1, 2, 2);
	plot(solution(:, 1), solution(:, 3:end));
	legend('Lysozyme', 'Cytochrome', 'Ribonuclease');
	grid on;
	xlabel('Time [s]');
	ylabel('Concentration [mM]');
end

% =============================================================================
%  CADET
%  
%  Copyright (C) 2008-2022: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
