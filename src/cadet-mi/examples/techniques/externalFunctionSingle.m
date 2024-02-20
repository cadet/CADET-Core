function externalFunctionSingle()
%EXTERNALFUNCTIONSINGLE Shows how to make the binding model depend on a single external function
%
%   In this example, the BREAKTHROUGHLANGMUIRSINGLE example is modified for a series
%   of pulse injections. The capacity of the Langmuir binding model is assumed to
%   depend on an external quantity that is measured at the column outlet and
%   transported inside the column with a certain velocity. 
%
%   See also BREAKTHROUGHLANGMUIRSINGLE, EXTERNALFUNCTIONMULTIPLE

% Copyright: (C) 2008-2024 The CADET Authors
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
	mLangmuir.kA_T       = [0.0];
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

	% Only capacity QMAX depends on external function in a linear way:
	%   QMAX(T) = 4.88 + 1.0 * T

	% Use external profile with index 0 for all parameters
	mLangmuir.externalSource = 0;

	mGrm.bindingModel = mLangmuir;
	
	% Specify inlet profile

	mGrm.constant       = zeros(6, mGrm.nComponents);
	mGrm.linear         = zeros(6, mGrm.nComponents);
	mGrm.quadratic      = zeros(6, mGrm.nComponents);
	mGrm.cubic          = zeros(6, mGrm.nComponents);

	% Sections 1, 3, 5: Pulse
	mGrm.constant([1, 3, 5], 1)  = 1;  % [mol / m^3] component 1

	% Construct external profile with sine signal
	extFun = LinearInterpolationExtFun();
	extFun.time = linspace(0, 6 * 100, 5001);
	extFun.profile = sin(linspace(0, 6 * 100, 5001) .* pi ./ 71);
	
	% Transport of external quantity is as fast as convective transport in
	% the column: interstitialVelocity / columnLength
	extFun.velocity = 5.75e-4 / 0.014;
	
	% Assign external function
	mGrm.externalFunctions = [extFun];

	% Construct and configure simulator
	sim = Simulator.create();
	sim.solutionTimes = linspace(0, 7 * 100, 5001);
	sim.sectionTimes = [[0:5] .* 100.0, 700];
	sim.sectionContinuity = [];
	
	% Lower error tolerances to make simulation pass (default values given
	% in parentheses)
	sim.absTol = 1e-6; %1e-10;
	sim.relTol = 1e-4; %1e-8;

	% Assign model
	sim.model = mGrm;

	% Run simulation and extract results
	res = sim.run();
	sol = [res.solution.time, squeeze(res.solution.outlet{1})];
	
	% Plot the chromatogram
	subplot(1, 3, 1);
	plot(sol(:, 1), sol(:, 2:end));
	h = legend('Lysozyme');
	set(h, 'Location', 'NorthWest');
	xlabel('Time [s]');
	ylabel('Concentration [mM]');
	title('Chromatogram QMAX = 4.88 + EXT');
	grid on;
	
	% Plot external function profile
	subplot(1, 3, 2);
	[X, Y] = meshgrid(linspace(0, 7 * 100, 501), linspace(0,1,51));
	Z = extFun.evaluate(X, Y, []);
	surf(X, Y, Z, 'EdgeColor', 'none');
	shading interp;
	grid on;
	xlabel('Time [s]');
	ylabel('Normalized axial position [-]');
	zlabel('External function');
	title('EXT profile');

	% Disable external function by setting its value to 0.0
	extFun.time = [0, sim.sectionTimes(end)];
	extFun.profile = [0, 0];
	
	% Run the simulation again with default error settings
	sim.absTol = 1e-10;
	sim.relTol = 1e-8;
	res2 = sim.run();
	sol2 = [res2.solution.time, squeeze(res2.solution.outlet{1})];
	
	% Plot the chromatogram
	subplot(1, 3, 3);
	plot(sol2(:, 1), sol2(:, 2:end));
	h = legend('Lysozyme');
	set(h, 'Location', 'NorthWest');
	xlabel('Time [s]');
	ylabel('Concentration [mM]');
	title('Chromatogram w/o external dependence');
	grid on;
end

% =============================================================================
%  CADET
%  
%  Copyright (C) 2008-2024: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================