function parallelSimulations()
%PARALLELSIMULATIONS Shows how to run several simulations in parallel using Matlab's parallel computing toolbox
%
%   CADET itself only offers parallelization on shared memory machines using OpenMP.
%   The number of threads employed is controlled by the MEXSIMULATOR.NTHREADS
%   property.
%
%   Another possibility for parallelization is given by Matlab's parallel computing
%   toolbox. This toolbox spawns several worker processes on which the total work
%   is distributed. The workers communicate via network or proprietary interconnects
%   and can run on remote machines or clusters (nodes). On each machine or node,
%   CADET's intrinsic shared-memory parallelization can be exploited to obtain a
%   hybrid parallelization.
%
%   In this example, a Monte Carlo study of the influence of interstitial velocity
%   and loading in the BREAKTHROUGHLANGMUIRSINGLE example is performed.
%
%   See also BREAKTHROUGHLANGMUIRSINGLE, PARAMETERIZEDSIMULATIONWITHOUTSENSITIVITIES

% Copyright: (C) 2008-2016 The CADET Authors
%            See the license note at the end of the file.

	% Create simulator
	sim = createSimulator();

	% Add parameters without computing sensitivities
	params = cell(2, 1);
	params{1} = makeSensitivity([0], {'VELOCITY'}, [-1], [-1], [-1], [-1]);
	params{2} = makeSensitivity([1], {'CONST_COEFF'}, [0], [-1], [-1], [0]);
	sim.setParameters(params, [false, false]);

	% Allocate space for results and samples (100 samples)
	result = cell(100, 1);
	samples = zeros(100, 2);

	% Run Monte Carlo main loop in parallel
	parfor i = 1:length(result)

		% Draw random values for velocity and loading concentration:
		% Since rand() produces uniform random values in [0, 1], the
		% expression -1 + 2*rand() transforms it to [-1, 1].
		samples(i, :) = [5.75e-4 + (-1 + 2*rand())*1e-4, 7.14e-3 + (-1 + 2*rand())*1e-3];

		% Run simulation and extract results
		localRes = sim.runWithParameters(samples(i, :), false);
		result{i} = [localRes.solution.time, localRes.solution.outlet{1}];
	end
	
	% Plot all curves in one figure
	hold on;
	for i = 1:length(result)
		curRes = result{i};
		plot(curRes(:, 1), curRes(:, 2:end));
	end
	hold off;
	legend('Lysozyme');
	xlabel('Time [s]');
	ylabel('Concentration [mM]');
	grid on;

end

function sim = createSimulator()

	% General rate model unit operation
	mGrm = SingleGRM();

	% Discretization
	mGrm.nComponents = 1;
	mGrm.nCellsColumn = 16; % Attention: This is very low and only used for illustration (short runtime)
	mGrm.nCellsParticle = 4; % Attention: This is very low and only used for illustration (short runtime)
	mGrm.nBoundStates = ones(mGrm.nComponents, 1);

	% Initial conditions (empty column)
	mGrm.initialBulk = [0.0];
	mGrm.initialSolid = [0.0];
		
	% Transport
	mGrm.dispersionColumn          = 5.75e-8;
	mGrm.filmDiffusion             = [6.9e-6];
	mGrm.diffusionParticle         = [6.07e-11];
	mGrm.diffusionParticleSurface  = [0.0];
	mGrm.interstitialVelocity      = 5.75e-4;

	% Geometry
	mGrm.columnLength        = 0.014;
	mGrm.particleRadius      = 4.5e-5;
	mGrm.porosityColumn      = 0.37;
	mGrm.porosityParticle    = 0.75;
	
	% Adsorption
	mLangmuir = LangmuirBinding();
	mLangmuir.kineticBinding = true;
	mLangmuir.kA         = [1.14];
	mLangmuir.kD         = [0.002];
	mLangmuir.qMax       = [4.88];
	mGrm.bindingModel = mLangmuir;
	
	% Specify inlet profile

	mGrm.constant       = zeros(1, mGrm.nComponents);
	mGrm.linear         = zeros(1, mGrm.nComponents);
	mGrm.quadratic      = zeros(1, mGrm.nComponents);
	mGrm.cubic          = zeros(1, mGrm.nComponents);

	% Section 1: Loading phase
	mGrm.constant(1,1)  = 7.14e-3;  % [mol / m^3] component 1

	% Construct and configure simulator
	sim = Simulator.create();
	sim.solutionTimes = linspace(0, 10000, 1001);
	sim.sectionTimes = [0.0 10000.0];
	sim.sectionContinuity = [];

	% Assign model
	sim.model = mGrm;
end

% =============================================================================
%  CADET - The Chromatography Analysis and Design Toolkit
%  
%  Copyright (C) 2008-2016: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
