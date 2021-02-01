function [params, residual, success] = runFittingTransformed(initParams)
%RUNFITTINGTRANSFORMED Performs a parameter estimation for a three component SMA system
%
%   The SMA parameters of three components are simultaneously estimated from
%   four gradient and one breakthrough experiments in which only the sum
%   signal of all components is observed. Artificial measurements without
%   noise are used for demonstration purposes.
%
%   This function uses a parameter transform on the characteristic charges to
%   prescribe the order in which the peaks appear. This should avoid the
%   problem of multiple global optima (permutation of components). However,
%   the characteristic charges become highly correlated.
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
%   RUNFITTINGTRANSFORMED(INITPARAMS) runs the parameter fit from the initial
%   parameters given in INITPARAMS. The order of the parameters is
%      ribonuclease k_a, ribonuclease nu, ribonuclease sigma,
%      lysozyme k_a, lysozyme nu, lysozyme sigma,
%      cytochrome k_a, cytochrome nu, cytochrome sigma.
%
%   [PARAMS, RESIDUAL, SUCCESS] = RUNFITTINGTRANSFORMED(...) returns a vector
%   PARAMS with the parameters at the end of the optimization, the residual
%   RESIDUAL at this point, and the logical SUCCESS indicating whether the
%   optimizer claims a solution. Returns -1 in RESIDUAL if something went wrong.

% Copyright: (C) 2008-2016 The CADET Authors
%            See the license note at the end of the file.

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
		
	% Assemble experiment descriptions
	for i = 0:length(gradLength)
		% Create simulator and set parameters
		% Breakthrough for i = 0, gradient elution for i = 1..4
		sim = createSimulator(i, gradLength);

		% Collect artificial data of experiment in cell array (one cell per wavelength)
		% Note that time points of the measurements are given in sim.solutionTimes
		res = sim.run();
		data = cell(1, 1);
		data{1} = sum(res.solution.outlet{1}(:, 2:end), 2);

		% Set fit parameters and enable sensitivities
		sim.setParameters(params, true(9, 1));

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

	% Setup parameter transform for nu
	% We know that nu_ribonuclease < nu_lysozyme < nu_cytochrome. We can
	% use this information to fix the order in which the peaks appear.
	% Reparameterize nu = f(NU), where NU are the new characteristic
	% charges.
	%                        / nu_1 \    / NU_1               \
	%      nu = f(NU)   <=>  | nu_2  | = | NU_1 + NU_2         |
	%                        \ nu_3 /    \ NU_1 + NU_2 + NU_3 /
	% The optimizer now uses increments on the previous component's
	% characteristic charge as variable.
	% This makes the characteristic charges highly correlated, since
	% changing one NU parameter also changes the values of all subsequent
	% characteristic charges nu.
	nuTransform = CustomParameterTransformation(...
		@transform, ... % Forward transformation simulator -> optimizer
		@invTransform, ... % Inverse transformation optimizer -> simulator
		@chainRuleInvTransform); % Apply Jacobian of transform to system Jacobian

	% Apply parameter transformation
	parFit.parameterTransform = nuTransform;

	% Upper and lower bounds of parameters (already transformed to optimizer space)
	loBound = repmat([1e-2, 0.01, 0.01], 1, 3);
	upBound = repmat([1e2, 10, 30], 1, 3);
	
	% Transform parameters and bounds to optimizer space
	initParams = parFit.applyTransform(initParams);
	loBound = parFit.applyTransform(loBound);
	upBound = parFit.applyTransform(upBound);

	% Storage for plot handles
	plotHd = [];

	% Set some options (tolerances, enable Jacobian) and request some output (iteration statistics, plots)
	opts = optimset('TolFun', 1e-8, 'TolX', 1e-8, 'MaxIter', 1000, 'MaxFunEvals', 200 * length(initParams), ...
		'Jacobian', 'on', 'Diagnostics', 'off', 'Display', 'iter', 'OutputFcn', @progressMonitor);

	try
		% Invoke optimizer lsqnonlin
		[params, residual, ~, exitflag] = lsqnonlin(@residualFunc, initParams, loBound, upBound, opts);
		success = (exitflag > 0);
	catch exc
		disp('ERROR in lsqnonlin. Probably due to failed CADET simulation.');
		disp(exc.message);

		params = parFit.lastParameters;
		if isempty(params)
			params = initParams;
		end
		residual = -1;
		success = false;
	end

	% Transform parameters back to simulator space
	params = parFit.applyInverseTransform(params);

	function [varargout] = residualFunc(x)
		%RESIDUALFUNC Residual function that is passed to lsqnonlin, just forwards to ParameterFit object
		varargout = cell(nargout,1);
		[varargout{:}] = parFit.residualVector(x);
	end

	function stop = progressMonitor(x, optimValues, state)
		%PROGRESSMONITOR Callback function for progress report invoked by lsqnonlin, just calls ParameterFit's plot function
		stop = false;
		if strcmp(state, 'iter')
			% Call plot function and reuse plot handles from previous call (storage outside this function in captured variable plotHd)
			plotHd = parFit.plot(plotHd, optimValues.iteration, [optimValues.resnorm, optimValues.stepsize]);
		end
	end
end

function out = transform(in)
	% Transform Sim -> Opt
	% nu_i in in([2, 5, 8])
	out = in;
	out(5) = in(5) - in(2);
	out(8) = in(8) - in(5);
end

function out = invTransform(in)
	% Transform Opt -> Sim
	out = in;
	out(5) = in(2) + in(5);
	out(8) = in(2) + in(5) + in(8);
end

function jac = chainRuleInvTransform(jac, transParam, origParam)
	% Apply chain rule, i.e., multiply given Jacobian with Jacobian of the 
	% transform: J_out = J_in * J_transform
	% The Jacobian of the transform is given (with respect to nu's) by
	%                 / 1 0 0 \
	%  J_transform =  | 1 1 0  |
	%                 \ 1 1 1 /
	jac(:,2) = jac(:,8) + jac(:,5) + jac(:,2);
	jac(:,5) = jac(:,8) + jac(:,5);
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
	mGrm.nCellsColumn = 64; % Attention: This is low and only used for illustration (shorter runtime)
	mGrm.nCellsParticle = 16; % Attention: This is low and only used for illustration (shorter runtime)
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
