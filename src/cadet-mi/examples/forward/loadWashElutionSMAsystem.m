function [result, solution] = loadWashElutionSMAsystem()
%LOADWASHELUTIONSMASYSTEM Simple Load-Wash-Elution cycle with 4 component SMA model
%
%   The model is taken from the publication:
%   Püttmann, A.; Schnittert, S.; Naumann, U.; von Lieres, E. (2013). 
%   Fast and accurate parameter sensitivities for the general rate model 
%   of column liquid chromatography. Computers & Chemical Engineering, 56, 46–57. 
%   doi:10.1016/j.compchemeng.2013.04.021
%   
%   It describes ion-exchange chromatography of lysozyme, cytochrome and ribonuclease 
%   on the strong cation-exchanger SP Sepharose FF. The column is first loaded for 10s 
%   with a constant ionic strength of 50 mol / m^3 and lysozyme, cytochrome and ribonuclease 
%   concentrations of 1 mol / m^3. The column is then washed for 80s at the same ionic 
%   strength without proteins. The bound proteins are then eluted with a salt gradient from 
%   T = 90s to T = 1500s, starting at 100 mol / m^3 and with a slope of 0.2 mol / (m^3 * s). 
%
%   In this function, a ModelSystem consisting of inlet, general rate model, and outlet unit 
%   operations is constructed explicitly. See loadWashElutionSMAsingle() for a shorter
%   way. Also note that for linear cascades an explicit outlet unit operation is not required,
%   since the outlet profile of the last unit operation can be used instead.
%
%   See also LOADWASHELUTIONSMASINGLE.

% Copyright: (C) 2008-2020 The CADET Authors
%            See the license note at the end of the file.


	% Step 1: Construct general rate model unit operation
	% ===================================================

	% General rate model unit operation
	mGrm = GeneralRateModel();

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


	% Step 2: Construct inlet unit operation
	% ===================================================

	% Inlet unit operation
	mIn = PiecewiseCubicPolyInlet();
	mIn.nComponents = 4;
	
	% Reserve space: nSections x nComponents (a section can be thought of being a 
	% step in the process, see below)
	mIn.constant       = zeros(3, mGrm.nComponents);
	mIn.linear         = zeros(3, mGrm.nComponents);
	mIn.quadratic      = zeros(3, mGrm.nComponents);
	mIn.cubic          = zeros(3, mGrm.nComponents);

	% Section 1: Loading phase
	mIn.constant(1,1)  = 50.0;  % [mol / m^3] component 1
	mIn.constant(1,2)  = 1.0;   % [mol / m^3] component 2
	mIn.constant(1,3)  = 1.0;   % [mol / m^3] component 3
	mIn.constant(1,4)  = 1.0;   % [mol / m^3] component 4

	% Section 2: Washing phase (no protein feed)
	mIn.constant(2,1)  = 50.0;  % [mol / m^3] component 1

	% Section 3: Elution phase (linear salt gradient with step at the beginning)
	mIn.constant(3,1)  = 100;  % [mol / m^3] component 1
	mIn.linear(3,1)  = 0.2;  % [mol / (m^3 * s)] component 1


	% Step 3: Construct outlet unit operation
	% ===================================================

	% Inlet unit operation
	mOut = OutletModel();
	mOut.nComponents = 4;
	
	
	% Step 4: Assemble system of unit operations
	% ===================================================

	% Construct ModelSystem and assign unit operations (order determines IDs)
	mSys = ModelSystem();
	mSys.models = [mIn, mGrm, mOut];

	% Define valve configurations / unit operation connections
	% Valve configuration active on entering section 0
	mSys.connectionStartSection = [0];
	% Connect unit 0 with unit 1 (-1 = all ports, all components), flow rate 1.0
	% Connect unit 1 with unit 2 (-1 = all ports, all components), flow rate 1.0
	mSys.connections = {[0, 1, -1, -1, -1, -1, 1.0; ...
	                     1, 2, -1, -1, -1, -1, 1.0]};


	% Step 5: Create simulator and configure it
	% ===================================================

	% Construct and configure simulator
	sim = Simulator.create();
	sim.solutionTimes = linspace(0, 1500, 1001); % [s], time points at which solution is computed

	% sectionTimes holds the sections and sectionContinuity indicates whether
	% the transition between two adjacent sections is continuous
	sim.sectionTimes = [0.0 10.0 90.0 1500.0]; % [s]
	sim.sectionContinuity = false(2,1);

	% Hand model over to simulator	
	sim.model = mSys;


	% Step 6: Run the model and plot the results
	% ===================================================

	% Run the model
	result = sim.run();

	% Extract solution into a matrix with time being the first column
	% Note that we need to extract the outlet of the third unit operation, 
	% which is the OutletModel (order in the ModelSystem matters)
	solution = [result.solution.time, result.solution.outlet{3}];

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
%  Copyright (C) 2008-2020: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
