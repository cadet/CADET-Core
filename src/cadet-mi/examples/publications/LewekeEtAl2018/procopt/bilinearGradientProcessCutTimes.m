function [cutPoints, optimalYield, optimalPurity] = bilinearGradientProcessCutTimes(gradientShape, mode, colLength, quiet)
%BILINEARGRADIENTPROCESSCUTTIMES Calculates optimal cut times with respect to purity and yield constraints
%  
%   A result of the optimal separation study is used to obtain a chromatogram
%   given the optimal bilinear gradient shape. Having the chromatogram at
%   hand, an optimization problem is solved which determines optimal cut
%   points to collect a target component subject to purity and yield
%   constraints.
%  
%   There are two modes:
%     1. Maximize yield subject to a 95 % purity constraint
%     2. Maximize purity subject to a 70 % yield constraint
%   The target component is chosen as Cytochrome (component 2, neglecting
%   salt).
%  
%   Introducing the mass m_i of component i between the cut points t_1 and 
%   t_2 as
%             / t_2
%      m_i = |      c_i(t) dt,
%            / t_1
%  
%   yield y_i is defined by
%      y_i = m_i / ( c_inj * t_inj )
%   and purity is given by
%      p_i = m_i / (m_1 + m_2 + m_3).
%   For the two cases the optimization problems are formulated as follows:
%     1:  max  y_2
%         s.t. p_2 >= 0.95
%              t_1 <= t_2
%     2:  max  p_2
%         s.t. y_2 >= 0.7
%              t_1 <= t_2
%   
%   The employed model describes ion-exchange chromatography of lysozyme,
%   cytochrome, and ribonuclease on the strong cation-exchanger 
%   SP Sepharose FF. Model parameters are taken from benchmark 2 of the 
%   following publication:
%   A. Püttmann, S. Schnittert, U. Naumann & E. von Lieres (2013).
%   Fast and accurate parameter sensitivities for the general rate model of
%   column liquid chromatography.
%   Computers & Chemical Engineering, 56, 46–57.
%   doi:10.1016/j.compchemeng.2013.04.021
%  
%   The chromatogram is simulated and converted to a cubic spline. In order
%   to evaluate m_i, the spline, which is a piecewise polynomial, is
%   integrated analytically by calculating its anti-derivative.
%
%   BILINEARGRADIENTPROCESSCUTTIMES(GRADIENTSHAPE, MODE) performs
%   cut time optimization for the given gradient shape GRADIENTSHAPE which
%   consists of height, slope, and length of the first gradient. The
%   objective of the optimization is governed by MODE (1 for maximum yield,
%   2 for maximum purity).
%
%   BILINEARGRADIENTPROCESSCUTTIMES(..., COLLENGTH) additionally
%   sets the length of the column which defaults to 0.014m.
%
%   BILINEARGRADIENTPROCESSCUTTIMES(..., COLLENGTH, QUIET) additionally
%   determines whether the optimizer's iterations are visualized and
%   statistics printed (FALSE, default) or no output is generated during 
%   optimization (TRUE).
%
%   [CUTPOINTS, OPTIMALYIELD, OPTIMALPURITY] = BILINEARGRADIENTPROCESSCUTTIMES(...)
%   returns the cut times in CUTPOINTS, the achieved yield OPTIMALYIELD
%   at the optimum and the achieved purity OPTIMALPURITY at the optimum.
%
% See also LOADWASHELUTIONSMASINGLE

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
	
	if (nargin <= 3) || isempty(quiet)
		quiet = false;
	end
	
	% Target component is 2 (index would be 3 when including salt)
	idxTarget = 2;
		
	% Create model
	sim = createSimulator(gradientShape(1), gradientShape(2), gradientShape(3), colLength);

	% Simulate chromatograms and remove salt (component 1)
	chrom = sim.run();
	chrom = [chrom.solution.time, chrom.solution.outlet{1}(:, 2:end)];
	
	% Total injected mass of each component
	injMass = sim.sectionTimes(2) .* sim.model.constant(1, 2:end);

	% Perform optimization
	[cutPoints, optimalYield, optimalPurity] = optimalCutTimesChromatogram(chrom, idxTarget, injMass, mode, false, quiet);
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
