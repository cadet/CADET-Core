function jac = jacobianFiniteDifferences(sim, params, fdSteps)
%JACOBIANFINITEDIFFERENCES Calculates the Jacobian by central finite differences
%   JAC = JACOBIANFINITEDIFFERENCES(SIM, PARAMS, FDSTEPS) estimates the Jacobian
%   of a simulation by using central finite differences. This function can replace
%   parameter sensitivities. The current simulation represented by SIM is run
%   several times with slightly changed parameters centered on PARAMS, the current
%   parameters. The vector FDSTEPS determines the relative step size for each
%   parameter. The returned cell array JAC contains the Jacobian for each unit
%   operation with the format [nTime, nComponents, nParams].

% Copyright: (C) 2008-2018 The CADET Authors
%            See the license note at the end of the file.

	for i = 1:length(params)
		% Right point of FD
		curParams = params;
		curParams(i) = curParams(i) * (1 + fdSteps(i));
		right = sim.runWithParameters(curParams, true);
		right = right.solution.outlet; % Format is [nTime, nComp]

		% Left point of FD
		curParams = params;
		curParams(i) = curParams(i) * (1 - fdSteps(i));
		left = sim.runWithParameters(curParams, true);
		left = left.solution.outlet; % Format is [nTime, nComp]

		if i == 1
			% Reserve space
			jac = cell(size(right));
			for j = 1:length(jac)
				jac{j} = zeros([size(right{j}), length(params)]);
			end
		end

		% Calculate FDs
		for j = 1:length(jac)
			jac{j}(:, :, i) = (right{j} - left{j}) ./ (2 * params(i) * fdSteps(i));
		end
	end
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
