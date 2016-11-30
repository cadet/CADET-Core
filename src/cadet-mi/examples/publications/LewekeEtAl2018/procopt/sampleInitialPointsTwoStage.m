function [gradientShapes, objVals] = sampleInitialPointsTwoStage(nSamples, colLength, doFwd)
%SAMPLEINITIALPOINTSTWOSTAGE Samples initial points for separation optimization of a three component SMA process.
%
%   The optimizer can vary the gradient shape of the elution phase (bilinear
%   gradient).
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
%   A bilinear gradient is assumed. The design parameters are (in order):
%     Start concentration of gradient 1 in mM, slope of gradient 1 in mM / s,
%     and length of gradient 1 in s.
%   The total process time is fixed and the gradients are assumed to be
%   continuous, i.e., there is no step between gradient 1 and 2. Furthermore,
%   the second gradient has to reach the final (pre-defined) high salt
%   concentration.
%
%   If simulations are enabled, the sampled gradient shape is checked for full
%   elution of the componentes (by checking for a closed mass balance). For
%   later use, overlaps are calculated for each accepted sample.
%
%   SAMPLEINITIALPOINTSTWOSTAGE(NSAMPLES) samples until NSAMPLES samples are
%   obtained.
%
%   SAMPLEINITIALPOINTSTWOSTAGE(..., COLLENGTH) additionally sets the column
%   length which defaults to 0.014m.
%
%   SAMPLEINITIALPOINTSTWOSTAGE(..., COLLENGTH, DOFWD) determines via DOFWD
%   whether simulations are performed to calculate overlaps of the samples.
%   Defaults to true.
%
%   [GRADIENTSHAPES, OBJVALS] = SAMPLEINITIALPOINTSTWOSTAGE(...) returns a
%   matrix GRADIENTSHAPES in which each row contains the gradient parameters
%   of a sample and a matrix OBJVALS in which each row contains the overlaps
%   of the corresponding gradient shape. The columns in OBJVALS represent
%      sum_ij sum_k abs( c_i(t_k) * c_j(t_k) ) where i < k,
%      max_ij {sum_k abs( c_i(t_k) * c_j(t_k) )} where i < k.

% Copyright: © 2015 Samuel Leweke, Eric von Lieres
%            See the license note at the end of the file.

	% Set default column length
	if (nargin <= 1) || isempty(colLength)
		colLength = 0.014;
	end

	if (nargin <= 2) || isempty(doFwd)
		doFwd = true;
	end

	% Comment out the following line to avoid resetting the RNG
	rng(0, 'twister');

	if doFwd
		% Parameters for the CADET solver
		params = cell(5, 1);

		% Each parameter is identified by 
		%   - its unit operation id (0-based index or -1 if independent),
		%   - name (according to CADET file format specification),
		%   - component index (0-based index or -1 if independent),
		%   - reaction index (0-based index or -1 if independent),
		%   - bound phase index (0-based index or -1 if independent), and
		%   - time integrator section index (0-based index or -1 if independent).

		% Parameter 1: CONST_COEFF of component 1 (salt, component 0 if read as 0-based) in 
		%              inlet unit operation (id 1) and time section 3 (2 if read as 0-based)
		params{1} = makeSensitivity([1], {'CONST_COEFF'}, [0], [-1], [-1], [-1], [2]);

		% Parameter 1: LIN_COEFF of component 1 (salt, component 0 if read as 0-based) in 
		%              inlet unit operation (id 1) and time section 3 (2 if read as 0-based)
		params{2} = makeSensitivity([1], {'LIN_COEFF'}, [0], [-1], [-1], [-1], [2]);

		% Parameter 1: CONST_COEFF of component 1 (salt, component 0 if read as 0-based) in 
		%              inlet unit operation (id 1) and time section 4 (3 if read as 0-based)
		params{3} = makeSensitivity([1], {'CONST_COEFF'}, [0], [-1], [-1], [-1], [3]);

		% Parameter 1: LIN_COEFF of component 1 (salt, component 0 if read as 0-based) in 
		%              inlet unit operation (id 1) and time section 4 (3 if read as 0-based)
		params{4} = makeSensitivity([1], {'LIN_COEFF'}, [0], [-1], [-1], [-1], [3]);

		% Parameter 5: SECTION_TIMES item 4 (beginning of second gradient, 3 if read as 0-based)
		params{5} = makeSensitivity([-1], {'SECTION_TIMES'}, [-1], [-1], [-1], [-1], [3]);
		
		% Create the model
		sim = createSimulator(colLength);

		% Set parameters (we don't need sensitivities here)
		sim.setParameters(params, false(5, 1));

		% Calculate injected mass
		massesIn = sim.sectionTimes(2) .* sim.model.constant(1, 2:end);

		% Generate rule for numerical integration such that
		% \int f(t) dt = numIntRule' * f(timePoints)
		numIntRule = simps(sim.solutionTimes);
	end
	
	elutionLimit = 1e-6;   % Threshold for complete elution in mM
	bp = getBasicParams();
	
	% Layout: start, slope1, len
	lb = [bp.initialSalt, 0.001, 120];
	ub = [500, 10, bp.endTime - bp.startTime - 120];
	logScale = [false, true, false];

	% Initialize storage variables
	nGenerated = 0;
	gradientShapes = zeros(nSamples, length(lb));
	objVals = zeros(nSamples, 2);
	
	% Sampling main loop
	while nGenerated < nSamples
		
		% Draw uniform random sample in box defined by lower and upper
		% bounds
		r = rand(1, length(lb));
		x = lb + (ub - lb) .* r;
		x(logScale) = exp(log(lb(logScale)) + log(ub(logScale) ./ lb(logScale)) .* r(logScale));
 
		% Check salt monotonicity constraint
		%    start + slope1 * len <= maxSalt
		if (x(1) + x(2) * x(3)) >= bp.maxSalt
			% Draw a new sample
			continue;
		end
		
		if doFwd
			% Set process parameters
			modifyModelParams(sim, x);
			% Perform simulation
			try
				result = sim.run(true);
			catch
				% Something went wrong, try again
				continue;
			end
			outlet = result.solution.outlet{1};
		
			% Check for complete elution: Simple check on end concentration
			if any(outlet(end, 2:end) > elutionLimit)
				% Draw a new sample
				continue;
			end

			% Check for complete elution: Check mass balance
			massesOut = numIntRule.' * outlet(:,2:end);
			if any(abs(massesIn(:) - massesOut(:)) > 0.05 .* massesIn(:))
				% More than 5% mass deviation of at least one component
				continue;
			end
		
			% Valid sample collected
			nGenerated = nGenerated + 1;
			gradientShapes(nGenerated, :) = x;
			
			% Compute objective function values
			
			% Overlap by multiplication
			ov12 = numIntRule.' * abs(outlet(:, 2) .* outlet(:, 3));
			ov13 = numIntRule.' * abs(outlet(:, 2) .* outlet(:, 4));
			ov23 = numIntRule.' * abs(outlet(:, 3) .* outlet(:, 4));

			objVals(nGenerated, 1) = ov12 + ov13 + ov23;
			objVals(nGenerated, 2) = max([ov12, ov13, ov23]);
						
		else
			% Valid sample collected
			nGenerated = nGenerated + 1;
			gradientShapes(nGenerated, :) = x;
		end
	end
	
end

function modifyModelParams(sim, start, slope1, len)
%MODIFYMODELPARAMETERS Modifies the model parameters subject to given process parameters.

	% Check for vector input
	if nargin == 2
		len = start(3);
		slope1 = start(2);
		start = start(1);
	end

	bp = getBasicParams();

	% Set start and slope for second gradient (continuity requires that
	% start(grad2) = end(grad1).)
	start2 = start + slope1 * len;
	slope2 = (bp.maxSalt - start2) / (bp.endTime - len - bp.startTime);

	sim.setVariableParameterValues([start, slope1, start2, slope2, sim.sectionTimes(3) + len]);
end

function [sim] = createSimulator(colLength)
%CREATESIMULATOR Creates the model and simulator of the process and returns it
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
	mGrm.nCellsColumn = 16; % Attention: This is low and only used for illustration (shorter runtime)
	mGrm.nCellsParticle = 4; % Attention: This is low and only used for illustration (shorter runtime)
	mGrm.nBoundStates = ones(mGrm.nComponents, 1); % Number of bound states for each component
	% Discretization is very coarse for demonstration purposes
	% However, a fine discretization / accurate solution is not required at
	% this point, since we are only looking for good initial values.

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
	mGrm.constant(3, 1)  = 100;  % [mol / m^3] component 1 (salt)
	mGrm.linear(3, 1)    = 0.2;  % [mol / (m^3 * s)] component 1 (salt)

	% Section 4: Gradient 2 (linear salt gradient continuous with respect to first gradient)
	mGrm.constant(4, 1)  = mGrm.constant(3, 1) + ((bp.endTime - bp.startTime) * 0.5) * mGrm.linear(3, 1);  % [mol / m^3] component 1 (salt)
	mGrm.linear(4, 1)    = 0.5;  % [mol / (m^3 * s)] component 1 (salt)

	% Construct and configure simulator
	sim = Simulator.create();
	sim.solutionTimes = linspace(0, bp.endTime, bp.endTime+1); % [s], time points at which solution is computed
	sim.nThreads = 2; % Use 2 CPU cores for computation
	sim.initStepSize = 1e-9; % Initial time step size when beginning a new section
	sim.maxSteps = 100000; % Maximum number of (internal) time integrator steps

	% sectionTimes holds the sections and sectionContinuity indicates whether
	% the transition between two adjacent sections is continuous

	% Load, Wash, Gradient1, Gradient2
	sim.sectionTimes = [0.0, 10.0, bp.startTime, (bp.startTime + bp.endTime) * 0.5, bp.endTime]; % [s]
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
%  
%  Copyright © 2015: Samuel Leweke¹, Eric von Lieres¹
%                                      
%    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
