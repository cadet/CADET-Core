function singleExpSeparateComponents()
%SINGLEEXPSEPARATECOMPONENTS Fits adsorption rates in a single experiment in which each component is observed separately
%
%   The goal is to fit the adsorption rates KA (which ultimately correspond to the equilibrium
%   constants due to quasi-stationary binding mode) in the load-wash-elution example using
%   artificial data. One experiment is performed and each component is observed separately
%   (i.e., fractionation data is used).
%
%   See also LOADWASHELUTIONSMASINGLE, PARAMETERIZEDSIMULATIONWITHSENSITIVITIES,
%      SINGLEEXPJOINTCOMPONENTS

% Copyright: (C) 2008-2018 The CADET Authors
%            See the license note at the end of the file.

	% Create simulation and obtain artificial data
	sim = createSimulation();
	res = sim.run();

	% Set model parameters and enable sensitivities
	params = cell(3, 1);
	params{1} = makeSensitivity([0], {'SMA_KA'}, [1], [0], [-1], [0], [-1]);
	params{2} = makeSensitivity([0], {'SMA_KA'}, [2], [0], [-1], [0], [-1]);
	params{3} = makeSensitivity([0], {'SMA_KA'}, [3], [0], [-1], [0], [-1]);
	sim.setParameters(params, true(3, 1));

	% Create fit object
	pf = ParameterFit();

	% Collect data of first experiment in cell array (one cell per observed component / wavelength)
	% Note that time points of the measurements are given in sim.solutionTimes
	data = cell(3, 1);
	data{1} = res.solution.outlet{1}(:, 2);
	data{2} = res.solution.outlet{1}(:, 3);
	data{3} = res.solution.outlet{1}(:, 4);

	% Specify which components are observed in which factor for each observation / wavelength
	idxComp = cell(3, 1);
	idxComp{1} = [0 1 0 0]; % First observation is component 2 (or 1 if read as 0-based) with factor 1
	idxComp{2} = [0 0 1 0]; % Second observation is component 3 (or 2 if read as 0-based) with factor 1
	idxComp{3} = [0 0 0 1]; % Third observation is component 4 (or 3 if read as 0-based) with factor 1

	% Add the experiment to the fit (unit operation 0 is observed, which is the GRM) including name
	% of the experiment and the different observations / wavelengths
	pf.addExperiment(data, sim, [0 0 0], idxComp, [], [], [], [], 'Experiment1', {'Cyt', 'RNase', 'Lys'});

	% Specify initial parameters and their lower and upper bounds
	initParams = [30, 0.4, 10];
	loBound = 1e-3 * ones(1, 3);
	upBound = 1e2 * ones(1, 3);

	% This variable serves as storage for the plot handles returned by the plot function of ParameterFit
	% in the OutputFcn below.
	plotHd = [];

	% Set some options (tolerances, enable Jacobian) and request some output (iteration statistics, plots)
	opts = optimset('TolFun', 1e-8, 'TolX', 1e-8, 'MaxIter', 100, 'Jacobian', 'on', 'Diagnostics', 'off', 'Display', 'iter', 'OutputFcn', @progressMonitor);
	% Invoke optimizer lsqnonlin
	[optimalParams, optimalRes, ~, exitflag] = lsqnonlin(@residual, initParams, loBound, upBound, opts);
	success = (exitflag > 0);

	function [varargout] = residual(x)
		%RESIDUAL Residual function that is passed to lsqnonlin, just forwards to ParameterFit object
		varargout = cell(nargout,1);
		[varargout{:}] = pf.residualVector(x);
	end

	function stop = progressMonitor(x, optimValues, state)
		%PROGRESSMONITOR Callback function for progress report invoked by lsqnonlin, just calls ParameterFit's plot function
		stop = false;
		if strcmp(state, 'iter')
			% Call plot function and reuse plot handles from previous call (storage outside this function in captured variable plotHd)
			plotHd = pf.plot(plotHd, optimValues.iteration, [optimValues.resnorm, optimValues.stepsize]);
		end
	end

end

function sim = createSimulation()
	% General rate model
	mGrm = GeneralRateModel();

	mGrm.nComponents = 4; % Ordering: Salt, Cytochrome, RNase, Lysozyme
	mGrm.nCellsColumn = 16; % Attention: This is very low and only used for illustration (short runtime)
	mGrm.nCellsParticle = 4; % Attention: This is very low and only used for illustration (short runtime)
	mGrm.nBoundStates = ones(mGrm.nComponents, 1);

	% Initial conditions
	mGrm.initialBulk = [50.0 0.0 0.0 0.0];
	mGrm.initialSolid = [1.2e3 0.0 0.0 0.0];
		
	% Transport
	mGrm.dispersionColumn          = 5.75e-8;
	mGrm.filmDiffusion             = [6.9e-6 6.9e-6 6.9e-6 6.9e-6];
	mGrm.diffusionParticle         = [7e-10 6.07e-11 6.07e-11 6.07e-11];
	mGrm.diffusionParticleSurface  = [0.0 0.0 0.0 0.0];
	mGrm.interstitialVelocity      = 5.75e-4;

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

	% Inlet
	mIn = PiecewiseCubicPolyInlet();
	mIn.nComponents = 4;
	
	mIn.constant       = zeros(3, mGrm.nComponents);
	mIn.linear         = zeros(3, mGrm.nComponents);
	mIn.quadratic      = zeros(3, mGrm.nComponents);
	mIn.cubic          = zeros(3, mGrm.nComponents);

	% Sec 1
	mIn.constant(1,1)  = 50.0;  % component 1
	mIn.constant(1,2)  = 1.0;   % component 2
	mIn.constant(1,3)  = 1.0;   % component 3
	mIn.constant(1,4)  = 1.0;   % component 4

	% Sec 2
	mIn.constant(2,1)  = 50.0;  % component 1

	% Sec 3
	mIn.constant(3,1)  = 100;  % component 1
	mIn.linear(3,1)  = 0.2;  % component 1

	% Model system
	mSys = ModelSystem();
	mSys.models = [mGrm, mIn];
	mSys.connectionStartSection = [0];
	mSys.connections = {[1, 0, -1, -1, 1.0]};

	% Configure simulator
	sim = Simulator.create();
	sim.sectionTimes = [0.0 10.0 90.0 1500.0];
	sim.sectionContinuity = false(2, 1);
	sim.solutionTimes = linspace(0, 1500, 1501);

	% Assign model
	sim.model = mSys;
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
