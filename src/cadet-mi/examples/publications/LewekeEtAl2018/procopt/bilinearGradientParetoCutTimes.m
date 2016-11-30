function [optPoints, optYield, optPurity] = bilinearGradientParetoCutTimes(gradientShape, mode, colLength)
%BILINEARGRADIENTPARETOCUTTIMES Calculates yield-purity pareto curve for cut-point optimization
%
%   Calculates a pareto front visualizing the trade-off between yield and
%   purity in cut time optimization. See optimalCutTimes.m for details on the
%   underlying optimization task, the employed model, and explanation of the 
%   different scenarios / modes.
%  
%   The pareto curve is computed by varying the purity (mode 1) or yield
%   (mode 2) constraint from 0 to 1 with a fixed step size. For each
%   constraint the optimal cut time problem is solved as described in
%   optimalCutTimes() function. The cut times, purities, and yields are
%   collected and cleaned (removal of optimization failures). Finally, the
%   points are checked for dominance / pareto-optimality and violating points
%   are removed.
%
%   BILINEARGRADIENTPARETOCUTTIMES(GRADIENTSHAPE, MODE) calculates a yield-purity
%   Pareto curve with respect to the cut points for the given gradient shape
%   GRADIENTSHAPE which consists of height, slope, and length of the first gradient.
%   The objective of the optimization is governed by MODE (1 for maximum yield,
%   2 for maximum purity).
%
%   BILINEARGRADIENTPARETOCUTTIMES(..., COLLENGTH) additionally
%   sets the length of the column which defaults to 0.014m.
%
%   [OPTPOINTS, OPTYIELD, OPTPURITY] = BILINEARGRADIENTPARETOCUTTIMES(...)
%   returns a matrix OPTPOINTS in which each row contains the two optimal
%   cut times t_1, t_2, a vector OPTYIELD with optimal yields, and a vector
%   OPTPURITY with optimal purities.

% Copyright: © 2015 Samuel Leweke, Eric von Lieres
%            See the license note at the end of the file.

	% Set some defaults
	if (nargin <= 0) || isempty(gradientShape)
		gradientShape = [84.5398519767975, 0.00230614840064746, 3153.77736306280];
	end

	if (nargin <= 1) || isempty(mode)
		mode = 1;
	end

	if (nargin <= 2) || isempty(colLength)
		colLength = 0.014;
	end	
	
	% Number of steps taken when varying purity / yield constraint
	nSteps = 1001;
	
	% Target component is 2 (index would be 3 when including salt)
	idxTarget = 2;
	
	% Create model
	sim = createSimulator(gradientShape(1), gradientShape(2), gradientShape(3), colLength);

	% Simulate chromatograms and remove salt (component 1)
	chrom = sim.run();
	chrom = [chrom.solution.time, chrom.solution.outlet{1}(:, 2:end)];
	
	% Create splines and anti-derivatives from chromatograms
	chromSpl = cell(size(chrom, 2) - 1, 1);
	chromInts = cell(size(chrom, 2) - 1, 1);
	for i = 1:length(chromSpl)
		chromSpl{i} = pchip(chrom(:,1), chrom(:,1+i));
		chromInts{i} = ppint(chromSpl{i});
	end

	% Total injected mass of each component
	injMass = sim.sectionTimes(2) .* sim.model.constant(1, 2:end);
	
	% Find time of target component's peak
	[~, initPos] = max(chrom(:, idxTarget + 1));
	initPos = chrom(initPos, 1);
	
	% Bounds and initial cut points
	lb = [0, 0];
	ub = sim.sectionTimes(end) .* [1, 1];
	initPos = initPos .* [0.99, 1.01];
	
	% Set options and perform optimization of the different scenarios
	steps = linspace(0, 1, nSteps);
	optYield = zeros(nSteps, 1);
	optPurity = zeros(nSteps, 1);
	optPoints = zeros(nSteps, 2);
	options = optimoptions('fmincon', 'MaxIter', 500, 'TolCon', 1e-10, 'TolFun', 1e-10, 'TolX', 1e-12, 'Algorithm', 'interior-point', ...
		'Diagnostics', 'off', 'GradObj', 'on', 'GradConstr', 'on', 'Display', 'none');
	if mode == 1
		% Iterate over different constraint magnitudes
		for i = 1:nSteps
			minPurity = steps(i);
			[optPoints(i, :), ~, flag] = fmincon(@yield, initPos, [1, -1], [0], [], [], lb, ub, @purity, options);

			if flag > 0
				% Calculate optimal purity and yield
				optYield(i) = -yield(optPoints(i, :));
				optPurity(i) = -purity(optPoints(i, :));
			else
				% Optimization failed
				optYield(i) = NaN();
				optPurity(i) = NaN();
			end
		end
	else
		% Iterate over different constraint magnitudes
		for i = 1:nSteps
			minYield = steps(i);
			[optPoints(i, :), ~, flag] = fmincon(@purity, initPos, [1, -1], [0], [], [], lb, ub, @yield, options);

			if flag > 0
				% Calculate optimal purity and yield
				optYield(i) = -yield(optPoints(i, :));
				optPurity(i) = -purity(optPoints(i, :));
			else
				% Optimization failed
				optYield(i) = NaN();
				optPurity(i) = NaN();
			end
		end
	end
	
	% Remove invalid points
	idxRemove = isnan(optYield);
	optYield(idxRemove) = [];
	optPurity(idxRemove) = [];
	optPoints(idxRemove, :) = [];

	% Remove non-pareto-optimal points
	idxRemove = [];
	for i = 1:length(optYield)
		idxCand = find(((optYield >= optYield(i)) & (optPurity > optPurity(i))) | ((optYield > optYield(i)) & (optPurity >= optPurity(i))), 1, 'first');
		if ~isempty(idxCand)
			idxRemove = [idxRemove, i];
		end
	end
	optYield(idxRemove) = [];
	optPurity(idxRemove) = [];
	optPoints(idxRemove, :) = [];
	
	% Sort points
	[~, idx] = sort(optPurity);
	optPurity = optPurity(idx);
	optYield = optYield(idx);
	optPoints = optPoints(idx,:);
	
	% Create plots
	subplot(1,2,1);
	plot(optYield, optPurity, '-x');
	xlabel('Yield');
	ylabel('Purity');
	xlim([min(optYield), 1]);
	grid on;
	title('Pareto front');
	
	subplot(1,2,2);
	plot(optYield, optPoints);
	xlabel('Yield');
	ylabel('Cut point [s]');
	xlim([min(optYield), 1]);
	legend('Start', 'Stop');
	grid on;
	title('Optimal cut points');

	function varargout = yield(x)
		% Calculates the yield
		
		% Evaluate yield formula
		y = mass(x, chromInts, idxTarget) / injMass(idxTarget);
		% Compute gradient
		grady = [ppval(chromSpl{idxTarget}, x)] ./ injMass(idxTarget);
		% First element needs to be negated because it is t_1 (lower
		% integral boundary)
		grady(1) = -grady(1);
		
		if nargout == 1
			% Call as objective function => max yield
			varargout{1} = -y;
		elseif nargout == 2
			% Call as objective function => max yield, return [y, grady]
			varargout{1} = -y;
			varargout{2} = -grady;
		elseif nargout == 4
			% Call as constraint => return [y,eq,grady,gradeq]
			varargout{1} = minYield - y;
			varargout{2} = [];
			varargout{3} = -grady';
			varargout{4} = [];
		end
	end

	function varargout = purity(x)
		% Calculates the purity
		
		% Precompute masses
		mSum = mass(x, chromInts, 1) + mass(x, chromInts, 2) + mass(x, chromInts, 3);
		mTarget = mass(x, chromInts, idxTarget);

		% Evaluate purity formula
		y = mTarget / mSum;
		% First element needs to be negated because it is t_1 (lower
		% integral boundary). Note that the product / quotient rule and the
		% chain rule have to be applied.
		grady = [ppval(chromSpl{idxTarget}, x) ./ mSum - mTarget * ( ppval(chromSpl{1}, x) + ppval(chromSpl{2}, x) + ppval(chromSpl{3}, x)) / mSum^2];
		grady(1) = -grady(1);
		
		if nargout == 1
			% Call as objective function => max purity
			varargout{1} = -y;
		elseif nargout == 2
			% Call as objective function => max purity, return [y, grady]
			varargout{1} = -y;
			varargout{2} = -grady;
		elseif nargout == 4
			% Call as constraint => return [y,eq,grady,gradeq]
			varargout{1} = minPurity - y;
			varargout{2} = [];
			varargout{3} = -grady';
			varargout{4} = [];
		end
	end

end

function m = mass(t, chromInts, idxTarget)
%MASS Evaluates the mass of the target component between the given cut points
% Evaluates the anti-derivative of the target component's chromatogram in
% order to calculate the definite integral.
	m = diff(ppval(chromInts{idxTarget}, t));
end

function [sim] = createSimulator(gradStart, gradSlope, gradLen, colLength)
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
	sim.solutionTimes = linspace(0, bp.endTime, bp.endTime+1); % [s], time points at which solution is computed
	sim.nThreads = 2; % Use 2 CPU cores for computation
	sim.initStepSize = 1e-9; % Initial time step size when beginning a new section

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
