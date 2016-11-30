function [samples, purities, yields] = riskAnalysisSampling(nSamples, stdDevRel, cutTimes, gradientShape, colLength)
%RISKANALYSISSAMPLING Samples yield and purity by varying loading concentrations
%
%   The probability of failing a purity requirement is analyzed by a Monte
%   Carlo simulation. It is assumed that the loading concentrations are
%   subject to noise following a multivariate Gaussian distribution. By first
%   sampling loading concentrations from this distribution and computing
%   purity (and yield) for each sample, purity samples are created. In a
%   second step the failure probability is evaluated by computing a mean and
%   the variation of the probability is assessed using a Monte Carlo
%   bootstrap. The cut points are fixed in this study and the the noise on
%   each component is assumed to be uncorrelated with given standard
%   deviation.
%
%   This function performs the sampling only. See RISKANALYSISPOSTPROCESS()
%   for the probability calculation.
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
%   RISKANALYSISSAMPLING(NSAMPLES, STDDEVREL) draws NSAMPLES samples with
%   a relative standard deviation of STDDEVREL for each component (relative
%   to its nominal loading concentration, STDDEVREL = 0.1 means 10 %).
%
%   RISKANALYSISSAMPLING(..., CUTTIMES, GRADIENTSHAPE) uses the cut points
%   given in CUTTIMES and the bilinear gradient shape described by the
%   vector GRADIENTSHAPE that contains start concentration, slope, and length
%   of the first gradient.
%
%   RISKANALYSISSAMPLING(..., CUTTIMES, GRADIENTSHAPE, COLLENGTH) additionally
%   sets the column length which defaults to 0.014m.
%
%   [SAMPLES, PURITIES, YIELDS] = RISKANALYSISSAMPLING(...) returns a matrix
%   SAMPLES in which each row contains the sampled loading concentrations of
%   lysozyme, cytochrome, and ribonuclease (in this order). Also returns a
%   vector PURITIES which contains the purity of each sample and a vector
%   YIELDS with the corresponding yields.
%
%   See also RISKANALYSISPOSTPROCESS

% Copyright: (C) 2008-2016 The CADET Authors
%            See the license note at the end of the file.

	if (nargin <= 2) || isempty(cutTimes)
		cutTimes = [3349.54310141146, 3759.18310763293];
	end

	if (nargin <= 3) || isempty(gradientShape)
		gradientShape = [84.5398519767975, 0.00230614840064746, 3153.77736306280];
	end

	% Set default column length
	if (nargin <= 4) || isempty(colLength)
		colLength = 0.014;
	end

	% Target component is 2 (index would be 3 when including salt)
	idxTarget = 2;

	% Comment out the following line to avoid resetting the RNG
	rng(0, 'twister');

	% Create the model
	sim = createSimulator(gradientShape(1), gradientShape(2), gradientShape(3), colLength, cutTimes);

	% Use loading concentrations of the different components as parameters.
	params = cell(3, 1);

	% Each parameter is identified by 
	%   - its unit operation id (0-based index or -1 if independent),
	%   - name (according to CADET file format specification),
	%   - component index (0-based index or -1 if independent),
	%   - reaction index (0-based index or -1 if independent),
	%   - bound phase index (0-based index or -1 if independent), and
	%   - time integrator section index (0-based index or -1 if independent).

	% Parameter 1: CONST_COEFF of component 2 (lysozyme, component 1 if read as 0-based) in 
	%              inlet unit operation (id 1) and time section 1 (0 if read as 0-based)
	params{1} = makeSensitivity([1], {'CONST_COEFF'}, [1], [-1], [-1], [-1], [0]);

	% Parameter 1: CONST_COEFF of component 3 (cytochrome, component 2 if read as 0-based) in 
	%              inlet unit operation (id 1) and time section 1 (0 if read as 0-based)
	params{2} = makeSensitivity([1], {'CONST_COEFF'}, [2], [-1], [-1], [-1], [0]);

	% Parameter 1: CONST_COEFF of component 4 (ribonuclease, component 3 if read as 0-based) in 
	%              inlet unit operation (id 1) and time section 1 (0 if read as 0-based)
	params{3} = makeSensitivity([1], {'CONST_COEFF'}, [3], [-1], [-1], [-1], [0]);

	% Set parameters and enable sensitivities
	sim.setParameters(params, false(size(params)));

	% Draw samples
	samples = (1 + randn(nSamples, 3) .* stdDevRel) .* repmat(sim.model.constant(1, 2:end), nSamples, 1);
	
	% Initialize storage variables
	purities = zeros(nSamples, 1);
	yields = zeros(nSamples, 1);
	
	% Progress monitoring
	tStart = tic;
	lastProg = -0.1;

	% Sampling main loop
	for i = 1:nSamples

		% Run simulation
		result = sim.runWithParameters(samples(i, :));

		% Total injected mass of each component
		injMass = sim.sectionTimes(2) .* sim.model.constant(1, 2:end);

		% Calculate definite integrals
		time = result.solution.time;
		outlet = result.solution.outlet{1};
		masses = [simpsonRule(time, outlet(:, 2)); ...
				  simpsonRule(time, outlet(:, 3)); ...
				  simpsonRule(time, outlet(:, 4))];
		
		% Extract purity and yield
		purities(i) = purity(masses, idxTarget);
		yields(i) = yield(masses, idxTarget) / injMass(idxTarget);

		% Progress monitoring
		curProg = floor(i / nSamples * 100 / 5);
		if (curProg > lastProg)
			tElapsed = toc(tStart);
			fprintf('Completed %d of %d (%g %%) - Elapsed %g sec (%g sec remaining)\n', [i, nSamples, curProg * 5, tElapsed, tElapsed * (nSamples / i - 1)]);
			lastProg = curProg;
		end
	end
end

function y = yield(masses, idxTarget)
%YIELD Calculates the yield
	y = masses(idxTarget);
end

function p = purity(masses, idxTarget)
%PURITY Calculates the purity
	p = masses(idxTarget) / sum(masses);
end

function r = simpsonRule(t, y)
%SIMPSONRULE Uses Simpson's rule to calculate a definite integral with uniform spacing
	r = (t(2) - t(1)) / 3 * ( y(1) + y(end) + 2 * sum(y(3:2:end-1)) + 4 * sum(y(2:2:end-1)));
end

function [sim] = createSimulator(gradStart, gradSlope, gradLen, colLength, cutTimes)
%CREATESIMULATOR Creates the simulator of the process and returns it
% The model parameters are taken from benchmark 2 of
% A. Püttmann, S. Schnittert, U. Naumann & E. von Lieres (2013).
% Fast and accurate parameter sensitivities for the general rate model of
% column liquid chromatography.
% Computers & Chemical Engineering, 56, 46–57.
% doi:10.1016/j.compchemeng.2013.04.021

	bp = getBasicParams();

	% General rate model unit operation
	mGrm = SingleGRM();

	% Discretization
	mGrm.nComponents = 4;
	mGrm.nCellsColumn = 64; % Attention: This is low and only used for illustration (shorter runtime)
	mGrm.nCellsParticle = 16; % Attention: This is low and only used for illustration (shorter runtime)
	mGrm.nBoundStates = ones(mGrm.nComponents, 1); % Number of bound states for each component

	% Components are (in order): Salt, lysozyme, cytochrome, ribonuclease

	% Initial conditions, equilibrated empty column (note that solid phase salt
	% concentration has to match ionic capacity to satisfy equilibrium assumption)
	mGrm.initialBulk = [bp.initialSalt 0.0 0.0 0.0]; % [mol / m^3], also used for the particle mobile phase
	mGrm.initialSolid = [1.2e3 0.0 0.0 0.0]; % [mol / m^3]
		
	% Transport
	mGrm.dispersionColumn          = 5.75e-8; % [m^2 / s]
	mGrm.filmDiffusion             = [6.9e-6 6.9e-6 6.9e-6 6.9e-6]; % [m/s]
	mGrm.diffusionParticle         = [7e-10 6.07e-11 6.07e-11 6.07e-11]; % [m^2 / s]
	mGrm.diffusionParticleSurface  = [0.0 0.0 0.0 0.0]; % [m^2 / s]
	mGrm.interstitialVelocity      = 5.75e-4; % [m/s]

	% Geometry
	mGrm.columnLength        = colLength; % [m]
	mGrm.particleRadius      = 4.5e-5; % [m]
	mGrm.porosityColumn      = 0.37; % [-]
	mGrm.porosityParticle    = 0.75; % [-]
	
	% Adsorption
	mSma = StericMassActionBinding();
	mSma.kineticBinding = false; % Quasi-stationary binding (rapid-equilibrium)
	mSma.lambda     = 1.2e3; % Ionic capacity [mol / m^3]
	mSma.kA         = [0.0 35.5 1.59 7.7]; % Adsorption rate [(m^3 / mol)^nu / s]
	mSma.kD         = [0.0 1000 1000 1000]; % Desorption rate [(m^3 / mol)^nu / s]
	mSma.nu         = [0.0 4.7 5.29 3.7]; % Characteristic charge [-]
	mSma.sigma      = [0.0 11.83 10.6 10.0]; % Steric factor [-]
	mGrm.bindingModel = mSma;
	% The first value in the vectors above is ignored since it corresponds
	% to salt, which is component 0.
	% Note that due to the rapid-equilibrium assumption the equilibrium
	% constant is given by k_a / k_d.
	
	% Specify inlet profile

	% Reserve space: nSections x nComponents (a section can be thought of being a 
	% step in the process, see below)
	mGrm.constant       = zeros(4, mGrm.nComponents);
	mGrm.linear         = zeros(4, mGrm.nComponents);
	mGrm.quadratic      = zeros(4, mGrm.nComponents);
	mGrm.cubic          = zeros(4, mGrm.nComponents);

	% Section 1: Loading phase
	mGrm.constant(1, 1)  = bp.initialSalt;  % [mol / m^3] component 1 (salt)
	mGrm.constant(1, 2)  = 1.0;   % [mol / m^3] component 2
	mGrm.constant(1, 3)  = 1.0;   % [mol / m^3] component 3
	mGrm.constant(1, 4)  = 1.0;   % [mol / m^3] component 4

	% Section 2: Washing phase (no protein feed)
	mGrm.constant(2, 1)  = bp.initialSalt;  % [mol / m^3] component 1 (salt)

	% Section 3: Gradient 1 (linear salt gradient with step at the beginning)
	mGrm.constant(3, 1)  = gradStart;  % [mol / m^3] component 1 (salt)
	mGrm.linear(3, 1)    = gradSlope;  % [mol / (m^3 * s)] component 1 (salt)

	% Section 4: Gradient 2 (linear salt gradient continuous with respect to first gradient)
	mGrm.constant(4, 1)  = gradStart + gradSlope * gradLen;  % [mol / m^3] component 1 (salt)
	mGrm.linear(4, 1)    = (bp.maxSalt - gradStart - gradSlope * gradLen) / (bp.endTime - gradLen - bp.startTime);  % [mol / (m^3 * s)] component 1 (salt)

	% Construct and configure simulator
	sim = Simulator.create();
	sim.nThreads = 2; % Use 2 CPU cores for computation
	sim.initStepSize = 1e-9; % Initial time step size when beginning a new section

	% Make sure to get an odd number of points because of Simpson's rule
	nPoints = ceil(diff(cutTimes)) * 2;
	if mod(nPoints,2) == 0
		nPoints = nPoints + 1;
	end
	sim.solutionTimes = linspace(cutTimes(1), cutTimes(2), nPoints); % [s], time points at which solution is computed

	% sectionTimes holds the sections and sectionContinuity indicates whether
	% the transition between two adjacent sections is continuous

	% Load, Wash, Gradient1, Gradient2
	sim.sectionTimes = [0.0, 10.0, bp.startTime, bp.startTime + gradLen, bp.endTime]; % [s]
	sim.sectionContinuity = false(3, 1);

	% Hand model over to simulator  
	sim.model = mGrm;
end

function p = getBasicParams()
%GETBASICPARAMETERS Returns a struct with basic process parameters

	% Total process duration in s
	p.endTime = 6000;
	
	% Start time of first gradient in s
	p.startTime = 90;
	
	% Salt buffer concentration in mM for loading and washing
	p.initialSalt = 50;
	
	% Maximum salt buffer concentration in mM
	p.maxSalt = 1000;
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
