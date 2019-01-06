function resumeContinuousTime()
%RESUMECONTINUOUSTIME Shows how to resume a previous simulation with continuous time
%
%   A simulation can be resumed from the point it has ended. This feature is
%   demonstrated with the standard load-wash-elution example as shown in
%   LOADWASHELUTIONSMASINGLE. The simulation is split into two halves and the
%   first one is performed. Then, the second half is computed by resuming the
%   simulation from the first half.
%
%   In this example, the time in the second half continues as well. Another
%   option is to reset the time for the second half and start from 0.0 again.
%   This approach is presented in RESUMERESETTIME.
%
%   Note that also parameter sensitivity computations are continued, which is
%   also demonstrated. See PARAMETERIZEDSIMULATIONWITHSENSITIVITIES on how
%   to configure parameter sensitivities.
%
%   See also PARAMETERIZEDSIMULATIONWITHSENSITIVITIES, RESUMERESETTIME,
%      LOADWASHELUTIONSMASINGLE, MEXSIMULATOR.RESUME

% Copyright: (C) 2008-2018 The CADET Authors
%            See the license note at the end of the file.

	% Create simulator
	sim = createSimulator();

	% Add a joint sensitivity
	params = cell(1, 1);
	params{1} = makeSensitivity([0,0], {'SMA_KA', 'SMA_KA'}, [1, 2], [0, 0], [-1, -1], [0, 0], [-1, -1], [], [1, 1.59 / 35.5]);
	sim.setParameters(params, [true]);

	% Compute first half from 0s to 750s
	sim.solutionTimes = linspace(0, 750, 1001);

	% Set the inlet profile such that it resembles the first half of the load-wash-elution example
	setInlet(sim, 0, 750);

	% Run the first half and extract solution and sensitivity
	res1 = sim.runWithParameters([]);
	sol1 = [res1.solution.time, res1.solution.outlet{1}];
	sens1 = [res1.solution.time, res1.sensitivity.jacobian{1}];

	% Prepare for second half:
	%   - Since time is not resetted, solutionTimes is set from 750s to 1500s
	%   - Set the inlet profile to the second half of the full profile

	sim.solutionTimes = linspace(750, 1500, 1001);
	setInlet(sim, 750, 1500);

	% Resume the simulation for the second half
	res2 = sim.resumeWithParameters([]);

	% Extract solution and sensitivity
	sol2 = [res2.solution.time, res2.solution.outlet{1}];
	sens2 = [res2.solution.time, res2.sensitivity.jacobian{1}];

	% Glue the two pieces together by appending the second half to the first one
	sol = [sol1; sol2];
	sens = [sens1; sens2];

	% Plot the results
	figure;
	subplot(1, 2, 1);
	plot(sol(:, 1), sol(:, 3:end));
	legend('Lysozyme', 'Cytochrome', 'Ribonuclease');
	xlabel('Time [s]');
	ylabel('Concentration [mM]');
	title('Chromatogram');
	grid on;

	subplot(1, 2, 2);
	plot(sens(:, 1, 1), sens(:, 3:end, 1));
	legend('Lysozyme', 'Cytochrome', 'Ribonuclease');
	xlabel('Time [s]');
	title('Sensitivity wrt. SMA_KA');
	grid on;
end

function sim = createSimulator()
	% General rate model
	mGrm = SingleGRM();

	mGrm.nComponents = 4;
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

	mGrm.constant       = zeros(3, mGrm.nComponents);
	mGrm.linear         = zeros(3, mGrm.nComponents);
	mGrm.quadratic      = zeros(3, mGrm.nComponents);
	mGrm.cubic          = zeros(3, mGrm.nComponents);

	% Sec 1
	mGrm.constant(1,1)  = 50.0;  % component 1
	mGrm.constant(1,2)  = 1.0;   % component 2
	mGrm.constant(1,3)  = 1.0;   % component 3
	mGrm.constant(1,4)  = 1.0;   % component 4

	% Sec 2
	mGrm.constant(2,1)  = 50.0;  % component 1

	% Sec 3
	mGrm.constant(3,1)  = 100;  % component 1
	mGrm.linear(3,1)  = 0.2;  % component 1

	% Create and configure simulator
	sim = Simulator.create();
	sim.sectionTimes = [0.0 10.0 90.0 1500.0];
	sim.sectionContinuity = false(2,1);
	
	% Assign model
	sim.model = mGrm;
end

function setInlet(sim, tStart, tEnd)
%SETINLET Sets the current inlet profile of SIM to an excerpt of the full load-wash-elution profile specified by TSTART and TEND

	% Original section times of the full cycle
	nomSecTimes = [0.0 10.0 90.0 1500.0];

	% Find index of section that contains tStart and tEnd
	secStart = find(nomSecTimes <= tStart, 1, 'last');
	secEnd = find(nomSecTimes >= tEnd, 1);

	% Number of sections the excerpt contains
	nSec = secEnd - secStart;

	% Calculate elapsed time from beginning of each section to tStart
	tInSec = tStart - nomSecTimes(1:end-1);
	idx = 1;

	% Begin with empty profile
	mGrm = sim.model;

	mGrm.constant       = zeros(nSec, mGrm.nComponents);
	mGrm.linear         = zeros(nSec, mGrm.nComponents);
	mGrm.quadratic      = zeros(nSec, mGrm.nComponents);
	mGrm.cubic          = zeros(nSec, mGrm.nComponents);

	% Fill in each section in the excerpt. A section
	% may not start at its own beginning and the
	% coefficients of the piecewise polynomials have
	% to be adjusted for the time shift.

	% Sec 1
	if (secStart <= 1) && (secEnd >= 1)
		mGrm.constant(idx,1)  = 50.0;  % component 1
		mGrm.constant(idx,2)  = 1.0;   % component 2
		mGrm.constant(idx,3)  = 1.0;   % component 3
		mGrm.constant(idx,4)  = 1.0;   % component 4
		idx = idx + 1;
	end

	% Sec 2
	if (secStart <= 2) && (secEnd >= 2)
		mGrm.constant(idx,1)  = 50.0;  % component 1
		idx = idx + 1;
	end

	% Sec 3
	if (secStart <= 3) && (secEnd >= 3)
		mGrm.constant(idx,1)  = 100 + 0.2 * max(0, tInSec(3));  % component 1
		mGrm.linear(idx,1)  = 0.2;  % component 1
		idx = idx + 1;
	end

	% Update simulator
	sim.sectionTimes = max(min(nomSecTimes(secStart:secEnd), tEnd) - tStart, 0) + tStart;
	sim.sectionContinuity = false(nSec-1, 1);
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
