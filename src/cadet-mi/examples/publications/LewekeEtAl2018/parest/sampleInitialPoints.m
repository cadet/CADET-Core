function [samples, objVals] = sampleInitialPoints(nSamples, transformedSampling)
%SAMPLEINITIALPOINTS Samples initial points for parameter estimation of a three component SMA process.
%
%   The SMA parameters of three components are simultaneously estimated from
%   four gradient and one breakthrough experiments in which only the sum
%   signal of all components is observed. Artificial measurements without
%   noise are used for demonstration purposes.
%
%   The model describes ion-exchange chromatography of lysozyme, cytochrome,
%   and ribonuclease on the strong cation-exchanger SP Sepharose FF. Model
%   parameters are taken from benchmark 2 of the following publication:
%   A. Püttmann, S. Schnittert, U. Naumann & E. von Lieres (2013).
%   Fast and accurate parameter sensitivities for the general rate model of
%   column liquid chromatography.
%   Computers & Chemical Engineering, 56, 46–57.
%   doi:10.1016/j.compchemeng.2013.04.021
%
%   Since there is no constraint on the order in which the peaks appear,
%   multiple optima exist, which are given by permuting the component's
%   parameters. There are two methods implemented in this function to handle
%   this:
%    1. The parameters of the samples are sorted by characteristic charge in
%       ascending order, such that nu_1 <= nu_2 <= nu_3. Each nu_i is sampled
%       from the interval [0, 30].
%    2. Sampling is conducted using a parameter transformation which samples
%       nu_1, an offset to nu_1 (nu_2 = nu_1 + offset_1), and an offset to
%       nu_2 (nu_3 = nu_2 + offset_2 = nu_1 + offset_1 + offset_2). Note that
%       nu_1 and the offsets are sampled from the interval [0, 10],
%       such that nu_i <= 30 for all i.
%
%   SAMPLEINITIALPOINTS(NSAMPLES) generates NSAMPLES parameter sets for
%   parameter estimation. It is guaranteed that each parameter set causes
%   (almost) full elution of the product.
%
%   SAMPLEINITIALPOINTS(..., TRANSFORMEDSAMPLING) determines whether a
%   parameter transformation is used that fixes the ordering of the peaks.
%   Defaults to false (i.e., method 1 described above).
%
%   [SAMPLES, OBJVALS] = SAMPLEINITIALPOINTS(...) returns a matrix SAMPLES
%   in which each row contains the parameters in the following order: 
%      ribonuclease k_a, ribonuclease nu, ribonuclease sigma,
%      lysozyme k_a, lysozyme nu, lysozyme sigma,
%      cytochrome k_a, cytochrome nu, cytochrome sigma.
%   Also returns a vector OBJVALS which contains residuals for each sample:
%      res = sum_i weight_i^2 * sum_j [c_{sim,i}(t_j) - c_{meas,i}(t_j)]^2,
%   where i denotes the index of the experiment and j the index of the
%   time point.

% Copyright: (C) 2008-2016 The CADET Authors
%            See the license note at the end of the file.

	if (nargin <= 1) || isempty(transformedSampling)
		transformedSampling = false;
	end

	% Comment out the following line to avoid resetting the RNG
	rng(0, 'twister');

	% Parameters to estimate
	params = cell(9, 1);
	params{1} = makeSensitivity([0], {'SMA_KA'},    [1], [-1], [-1], [0], [-1]);
	params{2} = makeSensitivity([0], {'SMA_NU'},    [1], [-1], [-1], [0], [-1]);
	params{3} = makeSensitivity([0], {'SMA_SIGMA'}, [1], [-1], [-1], [0], [-1]);
	params{4} = makeSensitivity([0], {'SMA_KA'},    [2], [-1], [-1], [0], [-1]);
	params{5} = makeSensitivity([0], {'SMA_NU'},    [2], [-1], [-1], [0], [-1]);
	params{6} = makeSensitivity([0], {'SMA_SIGMA'}, [2], [-1], [-1], [0], [-1]);
	params{7} = makeSensitivity([0], {'SMA_KA'},    [3], [-1], [-1], [0], [-1]);
	params{8} = makeSensitivity([0], {'SMA_NU'},    [3], [-1], [-1], [0], [-1]);
	params{9} = makeSensitivity([0], {'SMA_SIGMA'}, [3], [-1], [-1], [0], [-1]);
	
	% Gradient lengths in sec
	gradLength = [1000, 1400, 1800, 2200];
	
	% Create ParameterFit object
	parFit = ParameterFit();
	
	% Specify which components are observed in which factor (extinction coefficient) for each wavelength
	idxComp = cell(1, 1);
	idxComp{1} = [0 1 1 1]; % 280nm wavelength: Sum of all components (all factor 1)

	% Reserve space for injected masses
	massesIn = zeros(length(parFit)-1, 3);

	% Assemble experiment descriptions
	for i = 0:length(gradLength)
		% Create simulator and set parameters
		% Breakthrough for i = 0, gradient elution for i = 1..4
		sim = createSimulator(i, gradLength);
		sim.setParameters(params, false(9, 1)); % No parameter sensitivities required

		% Collect artificial data of experiment in cell array (one cell per wavelength)
		% Note that time points of the measurements are given in sim.solutionTimes
		res = sim.run();
		data = cell(1, 1);
		data{1} = sum(res.solution.outlet{1}(:, 2:end), 2);
		
		% Record injected masses for gradient experiments
		if i >= 1
			massesIn(i, :) = sim.model.constant(1, 2:end) .* (sim.sectionTimes(2) - sim.sectionTimes(1));
		end

		% Weight of this particular experiment among all experiments
		weight = 1 / max(data{1});

		% Add the experiment to the collection
		name = 'Breakthrough';
		if i > 0
			name = sprintf('Gradient %d', i);
		end
		parFit.addExperiment(data, sim, [0], idxComp, [], weight, [], [], name, {'280nm'});
	end

	% Link parameters between experiments
	nExperiments = length(gradLength) + 1;
	for i = 1:length(params)
		parFit.addLink(1:nExperiments, i * ones(1, nExperiments));
	end
	
	% Upper and lower bounds
	loBound = repmat([1e-2, 0.1, 0.01], 1, 3);
	if transformedSampling
		upBound = repmat([1e2, 10, 30], 1, 3);
	else
		upBound = repmat([1e2, 30, 30], 1, 3);
	end

	logScale = [true true true];

	% Generate samples
	[samples, objVals] = doSampling(nSamples, loBound, upBound, logScale, parFit, transformedSampling, massesIn);
end

function out = paramTransform(in)
	% Transform Opt (nu_1, offsets) -> Sim (nu_i)
	out = in;
	out(5) = in(2) + in(5);
	out(8) = in(2) + in(5) + in(8);
end

function [samples, objVals] = doSampling(nSamples, loBound, upBound, logScale, parFit, transformedSampling, massesIn)
%DOSAMPLING Runs the parameter sampling using method 1 or 2 for ambiguitiy elimination.

	% Preallocate arrays
	samples = zeros(nSamples, length(loBound));
	objVals = zeros(nSamples, 1);
	i = 1;
	
	% Progress monitoring
	tStart = tic;
	lastProg = -0.1;

	% Main loop
	while i <= nSamples
		% Obtain a (log-)uniform sample from the box specified by loBound and
		% upBound
		r = rand(size(loBound));
		curParams = loBound + (upBound - loBound) .* r;
		curParams(logScale) = exp(log(loBound(logScale)) + log(upBound(logScale) ./ loBound(logScale)) .* r(logScale));
		
		if transformedSampling
			% Method 2: The sample is in transformed space (loBound and 
			% upBound are already transformed) and we need to transform it
			% back to Simulator space (nu_i instead of nu_1, offsets)
			curParams = paramTransform(curParams);
		else
			% Method 1: Sort components in ascending characteristic charges
			idxBefore = [2 5 8];
			[~, idxSorted] = sort(curParams(idxBefore));
			idxAfter = idxBefore(idxSorted);
			% Permute in-place (KA, NU, and SIGMA)
			curParams(idxBefore-1) = curParams(idxAfter-1);
			curParams(idxBefore) = curParams(idxAfter);
			curParams(idxBefore+1) = curParams(idxAfter+1);
		end
		
		try
			% Evaluate objective function
			objVals(i) = parFit.residualSumOfSquares(curParams);

			% Check if protein has eluted completely for gradient experiments
			massBalanceFail = false;
			for j = 2:parFit.nExperiments
				res = parFit.lastSimulationResults{j};

				massesOut = trapz(res.solution.time, res.solution.outlet{1}(:, 2:end), 1);
				if any(abs(massesIn(j-1, :) - massesOut) > 0.05 .* massesIn(j-1, :))
					% More than 5% mass deviation of at least one component
					% Mass balance error (product still on column), try another
					% parameter set
					massBalanceFail = true;
					break;
				end
			end
			if massBalanceFail
				continue;
			end

			samples(i, :) = curParams;
		catch
			% Something went wrong, skip this sample and try again
			continue;
		end

		% Progress monitoring
		curProg = floor(i / nSamples * 100 / 5);
		if (curProg > lastProg)
			tElapsed = toc(tStart);
			fprintf('Completed %d of %d (%g %%) - Elapsed %g sec (%g sec remaining)\n', [i, nSamples, curProg * 5, tElapsed, tElapsed * (nSamples / i - 1)]);
			lastProg = curProg;
		end
		
		i = i + 1;
	end
end

function [sim] = createSimulator(gradIdx, gradLength)
%CREATESIMULATOR Creates the simulator of the process and returns it
% The model parameters are taken from benchmark 2 of
% A. Püttmann, S. Schnittert, U. Naumann & E. von Lieres (2013).
% Fast and accurate parameter sensitivities for the general rate model of
% column liquid chromatography.
% Computers & Chemical Engineering, 56, 46–57.
% doi:10.1016/j.compchemeng.2013.04.021

	% General rate model unit operation
	mGrm = SingleGRM();

	% Discretization
	mGrm.nComponents = 4;
	mGrm.nCellsColumn = 16; % Attention: This is low and only used for illustration (shorter runtime)
	mGrm.nCellsParticle = 4; % Attention: This is low and only used for illustration (shorter runtime)
	mGrm.nBoundStates = ones(mGrm.nComponents, 1); % Number of bound states for each component

	% Components are (in order): Salt, lysozyme, cytochrome, ribonuclease

	% Initial conditions, equilibrated empty column (note that solid phase salt
	% concentration has to match ionic capacity to satisfy equilibrium assumption)
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
	mSma.kineticBinding = false; % Quasi-stationary binding (rapid-equilibrium)
	mSma.lambda     = 1.2e3; % Ionic capacity [mol / m^3]
	mSma.kA         = [0.0 7.7 35.5 1.59]; % Adsorption rate [(m^3 / mol)^nu / s]
	mSma.kD         = [0.0 1000 1000 1000]; % Desorption rate [(m^3 / mol)^nu / s]
	mSma.nu         = [0.0 3.7 4.7 5.29]; % Characteristic charge [-]
	mSma.sigma      = [0.0 10.0 11.83 10.6]; % Steric factor [-]
	mGrm.bindingModel = mSma;
	% The first value in the vectors above is ignored since it corresponds
	% to salt, which is component 0.
	% Note that due to the rapid-equilibrium assumption the equilibrium
	% constant is given by k_a / k_d.
	
	% Construct and configure simulator
	sim = Simulator.create();
	sim.nThreads = 2; % Use 2 CPU cores for computation
	sim.maxSteps = 100000; % Maximum number of (internal) time integrator steps
	sim.consistentInitMode = 5; % Consistently initialize only at beginning of simulation for more speed (but less accurate solutions)

	% Specify inlet profile
	runningSalt = 50;  % Salt concentration in running buffer in mM
	if gradIdx > 0
		% Gradient elution

		highSalt = 350;  % High salt buffer in mM
		stepHeight = 50; % Salt step height in mM
		stripSalt = 500; % Salt concentration in strip step in mM
		
		% Steps: Load, wash, elution, strip

		% Start and end times of the sections in seconds
		sim.sectionTimes = [0.0, 10.0, 90.0, 90 + gradLength(gradIdx), 90 + gradLength(gradIdx) + 120];
		% Sets continuity of transitions between two subsequent sections
		sim.sectionContinuity = false(3, 1);

		% Spline coefficients are initialized with zeros. First index
		% determines the component, second one gives the section index

		% Reserve space: nSections x nComponents (a section can be thought of being a 
		% step in the process, see below)
		mGrm.constant  = zeros(4, mGrm.nComponents);
		mGrm.linear    = zeros(4, mGrm.nComponents);
		mGrm.quadratic = zeros(4, mGrm.nComponents);
		mGrm.cubic     = zeros(4, mGrm.nComponents);

		% Section 1: Loading phase
		mGrm.constant(1, 1)  = runningSalt;  % [mol / m^3] component 1 (salt)
		mGrm.constant(1, 2)  = 1.0;          % [mol / m^3] component 2 (ribonuclease)
		mGrm.constant(1, 3)  = 1.0;          % [mol / m^3] component 3 (lysozyme)
		mGrm.constant(1, 4)  = 1.0;          % [mol / m^3] component 4 (cytochrome)

		% Section 2: Washing phase (no protein feed)
		mGrm.constant(2, 1)  = runningSalt;  % [mol / m^3] component 1 (salt)

		% Sec 3: Elute
		mGrm.constant(3, 1)  = runningSalt + stepHeight;  % [mol / m^3] component 1 (salt)
		mGrm.linear(3, 1)    = (highSalt - stepHeight) / gradLength(gradIdx);  % [mol / (m^3 * s)] component 1 (salt)

		% Sec 4: Strip
		mGrm.constant(4, 1)  = stripSalt;  % [mol / m^3] component 1 (salt)
	else
		% Breakthrough - Only load phase

		% Start and end times of the sections in seconds
		sim.sectionTimes = [0.0 1500.0];
		% Sets continuity of transitions between two subsequent sections
		sim.sectionContinuity = [];

		% Spline coefficients are initialized with zeros. First index
		% determines the component, second one gives the section index

		% Reserve space: nSections x nComponents (a section can be thought of being a 
		% step in the process, see below)
		mGrm.constant  = zeros(1, mGrm.nComponents);
		mGrm.linear    = zeros(1, mGrm.nComponents);
		mGrm.quadratic = zeros(1, mGrm.nComponents);
		mGrm.cubic     = zeros(1, mGrm.nComponents);

		% Section 1: Loading phase
		mGrm.constant(1, 1)  = runningSalt;  % [mol / m^3] component 1 (salt)
		mGrm.constant(1, 2)  = 1.0;          % [mol / m^3] component 2 (ribonuclease)
		mGrm.constant(1, 3)  = 1.0;          % [mol / m^3] component 3 (lysozyme)
		mGrm.constant(1, 4)  = 1.0;          % [mol / m^3] component 4 (cytochrome)
	end

	% Measurement at every second
	sim.solutionTimes = linspace(0, sim.sectionTimes(end), sim.sectionTimes(end)+1); % [s], time points at which solution is computed

	% Hand model over to simulator  
	sim.model = mGrm;
end

% =============================================================================
%  CADET
%  
%  Copyright (C) 2008-2016: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
