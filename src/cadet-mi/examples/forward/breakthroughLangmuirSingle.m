function [result, solution] = breakthroughLangmuirSingle()
%BREAKTHROUGHLANGMUIRSINGLE Breakthrough curve with 1 component Langmuir model
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
%   See also PULSELINEARSINGLE.

% Copyright: (C) 2008-2022 The CADET Authors
%            See the license note at the end of the file.


	% Step 1: Construct a system with general rate model 
	%         as main unit operation
	% ===================================================

	% General rate model unit operation
	mGrm = SingleGRM();

	% Discretization
	mGrm.nComponents = 1;
	mGrm.nCellsColumn = 16; % Attention: This is very low and only used for illustration (short runtime)
	mGrm.nCellsParticle = 4; % Attention: This is very low and only used for illustration (short runtime)
	mGrm.nBoundStates = ones(mGrm.nComponents, 1); % Number of bound states for each component

	% Initial conditions (empty column)
	mGrm.initialBulk = [0.0]; % [mol / m^3], also used for the particle mobile phase
	mGrm.initialSolid = [0.0]; % [mol / m^3]
		
	% Transport
	mGrm.dispersionColumn          = 5.75e-8; % [m^2 / s]
	mGrm.filmDiffusion             = [6.9e-6]; % [m/s]
	mGrm.diffusionParticle         = [6.07e-11]; % [m^2 / s]
	mGrm.diffusionParticleSurface  = [0.0]; % [m^2 / s]
	mGrm.interstitialVelocity      = 5.75e-4; % [m/s]

	% Geometry
	mGrm.columnLength        = 0.014; % [m]
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
	
	% Specify inlet profile

	% Reserve space: nSections x nComponents (a section can be thought of being a 
	% step in the process, see below)
	mGrm.constant       = zeros(1, mGrm.nComponents);
	mGrm.linear         = zeros(1, mGrm.nComponents);
	mGrm.quadratic      = zeros(1, mGrm.nComponents);
	mGrm.cubic          = zeros(1, mGrm.nComponents);

	% Section 1: Loading phase
	mGrm.constant(1,1)  = 7.14e-3;  % [mol / m^3] component 1


	% Step 2: Create simulator and configure it
	% ===================================================

	% Construct and configure simulator
	sim = Simulator.create();
	sim.solutionTimes = linspace(0, 10000, 1001); % [s], time points at which solution is computed

	% sectionTimes holds the sections and sectionContinuity indicates whether
	% the transition between two adjacent sections is continuous
	sim.sectionTimes = [0.0 10000.0]; % [s]
	sim.sectionContinuity = [];

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
	legend('Lysozyme');
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
