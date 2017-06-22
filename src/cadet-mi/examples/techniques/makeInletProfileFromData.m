function makeInletProfileFromData()
%MAKEINLETPROFILEFROMDATA Constructs an inlet profile from data (measurements) and runs a simulation
%
%   It is shown how a measured profile is turned into an inlet profile suitable for 
%   simulation. A simulation with that inlet profile is then performed. The example
%   is based on BREAKTHROUGHLANGMUIRSINGLE.
%
%   See also BREAKTHROUGHLANGMUIRSINGLE.

% Copyright: (C) 2008-2017 The CADET Authors
%            See the license note at the end of the file.

	% Set up a general rate model
	mGrm = createModel();

	% Sample datapoints from a profile
	t = linspace(0, 1000, 1001).';
	f = @(t) (t <= 100) .* t ./ 100 + ((t > 100) & (t <= 200)) .* 1 + ((t > 200) & (t <= 400)) .* 0.5;
	data = f(t);

	% Convert solution into inlet profile
	inletProfile = PiecewiseCubicPolyProfile.fromUniformData(t, data);
	mGrm.inlet = inletProfile;
	
	% Create and configure simulator
	sim = Simulator.create();

	sim.solutionTimes = t;

	% The constructed inlet profile is a spline (piecewise polynomial)
	% that has its breaks at the time points of the solution. The inlet
	% profile requires that the section times of the simulator match its
	% breaks. Thus, we have to apply the breaks of the inlet profile to
	% the simulator.
	% We can also profit from setting the transition continuity to the
	% one of the profile.
	sim.sectionTimes = inletProfile.breaks;
	sim.sectionContinuity = inletProfile.continuity;
	
	% Assign model
	sim.model = mGrm;

	% Run the simulation and extract the solution for plotting it later
	res = sim.run();
	sol = [res.solution.time, res.solution.outlet{1}];

	% Convert inlet profile (component 1) to Matlab style piecewise polynomial
	pp = inletProfile.makePiecewisePoly([], 1);
	
	% Plot inlet profile
	figure;
	subplot(1, 2, 1);
	plot(linspace(0, 1000, 4001), ppval(pp, linspace(0, 1000, 4001)));
	legend('Lysozyme');
	xlabel('Time [s]');
	ylabel('Concentration [mM]');
	grid on;
	title('Inlet profile');
	
	% Plot chromatogram
	subplot(1, 2, 2);
	plot(sol(:,1), sol(:,2));
	legend('Lysozyme');
	xlabel('Time [s]');
	ylabel('Concentration [mM]');
	grid on;
	title('Chromatogram');
end

function mGrm = createModel()
	% General rate model
	mGrm = SingleGRM();

	mGrm.nComponents = 1;
	mGrm.nCellsColumn = 16; % Attention: This is very low and only used for illustration (short runtime)
	mGrm.nCellsParticle = 4; % Attention: This is very low and only used for illustration (short runtime)
	mGrm.nBoundStates = ones(mGrm.nComponents, 1);

	% Initial conditions
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
	mLangmuir = LangmuirBinding();
	mLangmuir.kineticBinding = true;
	mLangmuir.kA         = [1.14];
	mLangmuir.kD         = [0.002];
	mLangmuir.qMax       = [4.88];
	mGrm.bindingModel = mLangmuir;
end

% =============================================================================
%  CADET - The Chromatography Analysis and Design Toolkit
%  
%  Copyright (C) 2008-2017: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
