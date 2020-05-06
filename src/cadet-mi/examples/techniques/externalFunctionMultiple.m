function externalFunctionMultiple()
%EXTERNALFUNCTIONMULTIPLE Shows how to make the binding model depend on a single external function
%
%   In this example, the BREAKTHROUGHLANGMUIRSINGLE example is modified for a series
%   of pulse injections. The capacity of the Langmuir binding model is assumed to
%   depend on an external quantity that is measured at the column outlet and
%   transported inside the column with a certain velocity. Additionally, the
%   equilibrium constant is depending on another external quantity.
%
%   See also BREAKTHROUGHLANGMUIRSINGLE, EXTERNALFUNCTIONSINGLE

% Copyright: (C) 2008-2020 The CADET Authors
%            See the license note at the end of the file.

	% General rate model
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
	mLangmuir = ExtFunLangmuirBinding();
	mLangmuir.kineticBinding = false;

	% Basic values
	mLangmuir.kA         = [570.0];
	mLangmuir.kD         = [1.0];
	mLangmuir.qMax       = [4.88];

	% Linear dependence on external profile
	mLangmuir.kA_T       = [-500.0];
	mLangmuir.kD_T       = [0.0];
	mLangmuir.qMax_T     = [1.0];

	% Quadratic dependence on external profile
	mLangmuir.kA_TT      = [0.0];
	mLangmuir.kD_TT      = [0.0];
	mLangmuir.qMax_TT    = [0.0];

	% Cubic dependence on external profile
	mLangmuir.kA_TTT     = [0.0];
	mLangmuir.kD_TTT     = [0.0];
	mLangmuir.qMax_TTT   = [0.0];

	% Capacity QMAX depends on external function in a linear way:
	%   QMAX(T) = 4.88 + 1.0 * T
	% Equilibrium constant (adsorption rate) KA depends on external function in a linear way:
	%   EQ(T) = 570 - 500.0 * T

	% Use external profile 0 for KA and KD, external profile 1 for QMAX
	mLangmuir.externalSource = [0, 0, 1];

	mGrm.bindingModel = mLangmuir;
	
	% Specify inlet profile

	mGrm.constant       = zeros(6, mGrm.nComponents);
	mGrm.linear         = zeros(6, mGrm.nComponents);
	mGrm.quadratic      = zeros(6, mGrm.nComponents);
	mGrm.cubic          = zeros(6, mGrm.nComponents);

	% Sections 1, 3, 5: Pulse
	mGrm.constant([1, 3, 5], 1)  = 1;  % [mol / m^3] component 1

	% Construct external profile 1 with linear ramp
	extFun1 = LinearInterpolationExtFun();
	extFun1.time = [0, 600];
	extFun1.profile = [0, 1];
	
	% Transport of external quantity is as fast as convective transport in
	% the column: interstitialVelocity / columnLength
	extFun1.velocity = 5.75e-4 / 0.014;

	% Construct external profile 2 with sine signal
	extFun2 = LinearInterpolationExtFun();
	extFun2.time = linspace(0, 6 * 100, 5001);
	extFun2.profile = sin(linspace(0, 6 * 100, 5001) .* pi ./ 71);
	extFun2.velocity = 5.75e-4 / 0.014;
	
	% Assign external functions
	mGrm.externalFunctions = [extFun1, extFun2];

	% Construct and configure simulator
	sim = Simulator.create();
	sim.solutionTimes = linspace(0, 7 * 100, 5001);
	sim.sectionTimes = [[0:5] .* 100.0, 700];
	sim.sectionContinuity = [];
	
	% Lower error tolerances to make simulation pass (default values given
	% in parentheses)
	sim.absTol = 1e-5; %1e-10;
	sim.relTol = 1e-3; %1e-8;

	% Assign model
	sim.model = mGrm;

	% Run simulation and extract results
	res = sim.run();
	sol = [res.solution.time, squeeze(res.solution.outlet{1})];
	
	% Plot the chromatogram
	subplot(1, 2, 1);
	plot(sol(:, 1), sol(:, 2:end));
	h = legend('Lysozyme');
	set(h, 'Location', 'NorthWest');
	xlabel('Time [s]');
	ylabel('Concentration [mM]');
	title('Chromatogram with external dependence');
	grid on;
	
	% Disable external functions by setting their values to 0.0
	extFun1.time = [0, sim.sectionTimes(end)];
	extFun1.profile = [0, 0];
	extFun2.time = [0, sim.sectionTimes(end)];
	extFun2.profile = [0, 0];
	
	% Run the simulation again with default error settings
	sim.absTol = 1e-10;
	sim.relTol = 1e-8;
	res2 = sim.run();
	sol2 = [res2.solution.time, squeeze(res2.solution.outlet{1})];
	
	% Plot the chromatogram
	subplot(1, 2, 2);
	plot(sol2(:, 1), sol2(:, 2:end));
	h = legend('Lysozyme');
	set(h, 'Location', 'NorthWest');
	xlabel('Time [s]');
	ylabel('Concentration [mM]');
	title('Chromatogram w/o external dependence');
	grid on;
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