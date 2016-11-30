function [samples, sampledShape, inTime, optResult] = sampleParameters(idxPulseLength, opMode, configOpts)
%SAMPLEPARAMETERS Samples parameter estimates for perturbed pulse shapes
%
%   The function collects samples of parameter estimates for perturbed pulse
%   shapes in a one component system using a Langmuir isotherm. The model
%   considered here describes affinity chromatography of lysozyme on Cibacron
%   Blue Sepharose CL-6B at pH 2.7. Model parameters are based on
%   benchmark 1 of the following publication:
%   A. Püttmann, S. Schnittert, U. Naumann & E. von Lieres (2013).
%   Fast and accurate parameter sensitivities for the general rate model of
%   column liquid chromatography.
%   Computers & Chemical Engineering, 56, 46–57.
%   doi:10.1016/j.compchemeng.2013.04.021
%
%   The motivation of this study is given by the question whether, given a
%   specific amount of protein, it is better to use long pulses with low
%   concentration or short pulses with high concentrations to estimate
%   isotherm parameters.
%
%   A Monte-Carlo approach is employed to tackle this question. Ten different
%   base pulse shapes (each injecting the same amount of protein) are
%   considered and the effect of measurement uncertainty in injected mass and
%   pulse length are studied. Three different estimation scenarios are
%   considered: Estimation of adsorption constant (k_a) and capacity (q_max),
%   estimation of k_a only, and estimation of q_max only.
%   This function performs a Monte-Carlo run for one case using one base
%   pulse profile.
%
%   Given a base pulse profile, the samples are obtained by perturbing pulse 
%   length and injected mass independently. The artificial measurements 
%   obtained by performing a forward simulation with the perturbed shape are 
%   used to estimate the isotherm parameters (k_a, q_max) starting from the 
%   references values used for measurement generation. The reference isotherm
%   parameters remain constant for all base profiles and samples.
%
%   [SAMPLES, SAMPLEDSHAPE, INTIME, OPTRESULT] = SAMPLEPARAMETERS(IDXPULSELENGTH, OPMODE, CONFIGOPTS)
%   samples optimal parameters via optimization for the base pulse shape
%   denoted by IDXPULSELENGTH (index between 1 and 10), OPMODE determines
%   which parameters are to be sampled (1 = estimate k_a and q_max,
%   2 = estimate k_a, 3 = estimate q_max. Defaults to 1 (k_a and q_max)).
%   CONFIGOPTS is a struct with the following fields:
%      o nSamples: Number of samples to draw
%      o stdInMassAbs: Standard deviation of absolute Gaussian noise on
%           the injected mass
%      o stdInMassRel: Standard deviation of relative Gaussian noise on
%           the injected mass
%      o stdInTimeAbs: Standard deviation of absolute Gaussian noise on
%           the pulse length
%      o stdInTimeRel: Standard deviation of relative Gaussian noise on
%           the pulse length
%      o measNoiseAbs: Standard deviation of absolute Gaussian noise added
%           to the artificial measurements
%      o measNoiseRel: Standard deviation of relative Gaussian noise added
%           to the artificial measurements
%      o saveToFile: Flag which determines whether the results of this
%           function are saved to file.
%   Returns a matrix SAMPLES in which each row contains the result of
%   estimating the parameters for one sample. The first columns give the
%   estimated isotherm parameters, the last but one column holds the obtained
%   residual of the optimizer, and the last column contains a flag whether
%   the optimization was successful. SAMPLEDSHAPE is a matrix in which each
%   row contains the randomly sampled pulse length and height for each sample,
%   and INTIME contains the length of the base pulse in seconds. OPTRESULT
%   is a cell array with optimizer output (e.g., number of iterations) for
%   each sample.
%
%   See also RUNWORKFLOW

% Copyright: (C) 2008-2016 The CADET Authors
%            See the license note at the end of the file.

	% Isotherm parameters to estimate
	params = cell(2, 1);
	params{1} = makeSensitivity([0], {'MCL_KA'}, [0], [-1], [-1], [0], [-1]);
	params{2} = makeSensitivity([0], {'MCL_QMAX'}, [0], [-1], [-1], [0], [-1]);

	% Reference / true isotherm parameters
	paramVals = [1.14, 4.88];

	if opMode == 2
		% Estimate only k_a and fix q_max
		paramVals = paramVals(1);
		params = params(1);
	elseif opMode == 3
		% Estimate only q_max and fix k_a
		paramVals = paramVals(2);
		params = params(2);
	end

	% Create the model
	sim = createSimulator();

	% Set parameters
	sim.setParameters(params, true(numel(params), 1));

	% Length of pulse injection in seconds
	inTime = [60; 120; 180; 240; 300; 360; 420; 450; 480; 510; 540; 570; 600; 660; 720; 780; 840; 900]; % [s]
	inTime = inTime(idxPulseLength);

	% Nominal mass
	nomMass = 1e-6; % [kg]
	nomAmount = nomMass / 14.3; % [mol] (Lysozyme has molar mass of 14.3 kg/mol)

	% Preparative calculations for total injected volume
	columnRadius = 0.005; % [m]
	columnCrossSection = columnRadius^2 * pi; % [m^2]
	volFlowRate = columnCrossSection * sim.model.porosityColumn * sim.model.interstitialVelocity; % [m^3 / s]

	% Parameter bounds
	loBound = [1e-10 1e-8];
	upBound = [];

	if opMode == 2
		loBound = loBound(1);
	elseif opMode == 3
		loBound = loBound(2);
	end

	% Change some options of lsqnonlin
	opts = optimset('TolFun', 1e-8, 'TolX', 1e-8, 'MaxIter', 1000, 'MaxFunEvals', 200 * length(paramVals), 'Jacobian', 'on');

	% Preallocate space for sampled parameters
	samples = zeros(configOpts.nSamples, length(params)+2);
	sampledShape = zeros(configOpts.nSamples, 2);
	optResult = cell(configOpts.nSamples, 1);

	% Generate samples
	for i = 1:configOpts.nSamples

		% Step 1: Sample chromatogram

		% First, draw random numbers for pulse length and concentration
		sampledShape(i, 1) = configOpts.stdInTimeAbs * randn(1) + inTime .* (1 + configOpts.stdInTimeRel * randn(1));
		sampledShape(i, 2) = configOpts.stdInMassAbs * randn(1) + nomMass .* (1 + configOpts.stdInMassRel * randn(1));

		% Apply values to model
		injectedVolume = sampledShape(i, 1) * volFlowRate; % [m^3]
		amount = sampledShape(i, 2) / 14.3; % [mol] (Lysozyme has molar mass of 14.3 kg/mol)

		sim.sectionTimes(2) = sampledShape(i, 1);
		sim.model.constant(1, 1) = amount / injectedVolume;

		% Temporarily disable sensitivities for forward simulation
		disableSensitivities(sim);

		% Create artificial measurements
		measData = sim.runWithParameters(paramVals);
		measData = measData.solution.outlet{1};

		% Turn on sensitivities again for parameter fitting
		enableSensitivities(sim);

		% Create fit object
		pf = ParameterFit();
		idxComp = cell(1, 1);
		idxComp{1} = [1];

		% Add noise and update measurements in experiment
		data = cell(1, 1);
		data{1} = measData .* (1 + configOpts.measNoiseRel * randn(size(measData))) + configOpts.measNoiseAbs .* randn(size(measData));

		% Add the fit
		pf.addExperiment(data, sim, [0], idxComp, [], [], [1e3], [], 'Experiment', {'Lysozyme'});
		% Note that the weight has no theoretical effect since only one
		% experiment is used. However, the scaling applied by the weight
		% helps the optimizer (otherwise, the tolerance for lsqnonlin set
		% above would have to be lowered drastically).

		% Reset model to "true" pulse shape
		sim.sectionTimes(2) = inTime;
		sim.model.constant(1, 1) = nomAmount / (inTime * volFlowRate);

		% Step 2: Estimate parameters starting from true / reference values

		% Invoke optimizer lsqnonlin
		[vals, residual, ~, exitflag, optRes] = lsqnonlin(@residualFunc, paramVals, loBound, upBound, opts);
		samples(i, :) = [vals, residual, exitflag];
		optResult{i} = optRes;

		% Print progress
		fprintf('%g %%: %d of %d\n', i / configOpts.nSamples * 100, i, configOpts.nSamples);
	end

	% Save data to file
	if isfield(configOpts, 'saveToFile') && configOpts.saveToFile
		save(['ExpDes-est' num2str(opMode) '-shape' num2str(idxPulseLength) '.mat'], 'samples', 'configOpts', 'inTime', 'paramVals', 'loBound', 'upBound', 'sampledShape', 'opts', 'params', 'optResult');
	end

	function [varargout] = residualFunc(x)
		%RESIDUALFUNC Residual function that is passed to lsqnonlin, just forwards to ParameterFit object
		varargout = cell(nargout,1);
		[varargout{:}] = pf.residualVector(x);
	end

	function disableSensitivities(sim)
		%DISABLESENSITIVITIES Disables sensitivity computation

		if sim.nSensitiveParameters > 0
			% Remove sensitivities but keep as variable parameters
			sim.clearSensitivities(true);
		end
	end

	function enableSensitivities(sim)
		%ENABLESENSITIVITIES Enables sensitivity computation

		if sim.nSensitiveParameters > 0
			return;
		end

		sim.clearParameters();
		sim.setParameters(params, true(numel(params), 1));
	end
end

function [sim] = createSimulator()
%CREATESIMULATOR Creates the model of the process and returns it
% The model parameters are based on benchmark 1 of
% A. Püttmann, S. Schnittert, U. Naumann & E. von Lieres (2013).
% Fast and accurate parameter sensitivities for the general rate model of
% column liquid chromatography.
% Computers & Chemical Engineering, 56, 46–57.
% doi:10.1016/j.compchemeng.2013.04.021
% The binding mode has been set to quasi-stationary and the equilibrium
% constant has been decreased by two orders of magnitude.

	% General rate model unit operation
	mGrm = SingleGRM();

	% Discretization
	mGrm.nComponents = 1;
	mGrm.nCellsColumn = 64; % Attention: This is very low and only used for illustration (short runtime)
	mGrm.nCellsParticle = 16; % Attention: This is very low and only used for illustration (short runtime)
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
	mLangmuir.kineticBinding = false; % Kinetic binding
	mLangmuir.kA         = [1.14]; % Adsorption rate [m^3 / (mol * s)]
	mLangmuir.kD         = [0.2]; % Desorption rate [1 / s]
	mLangmuir.qMax       = [4.88]; % Capacity [mol / m^3]
	mGrm.bindingModel = mLangmuir;

	% Specify inlet profile

	% Reserve space: nSections x nComponents (a section can be thought of being a 
	% step in the process, see below)
	mGrm.constant       = zeros(2, mGrm.nComponents);
	mGrm.linear         = zeros(2, mGrm.nComponents);
	mGrm.quadratic      = zeros(2, mGrm.nComponents);
	mGrm.cubic          = zeros(2, mGrm.nComponents);

	% Section 1: Loading phase
	mGrm.constant(1, 1)  = 1;  % [mol / m^3] component 1 - will be overwritten later

	% Construct and configure simulator
	sim = Simulator.create();
	sim.solutionTimes = linspace(0, 1600, 1601); % [s], time points at which solution is computed

	% sectionTimes holds the sections and sectionContinuity indicates whether
	% the transition between two adjacent sections is continuous
	sim.sectionTimes = [0, 60, 1600]; % [s]
	sim.sectionContinuity = [false];

	sim.nThreads = 2; % Use 2 CPU cores for computation

	% Hand model over to simulator	
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
