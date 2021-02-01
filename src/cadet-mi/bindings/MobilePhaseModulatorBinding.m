
classdef MobilePhaseModulatorBinding < KineticQuasiStationaryBindingModel
	%MobilePhaseModulatorBinding Mobile phase modulator binding model
	%
	% See also BINDINGMODEL, KINETICQUASISTATIONARYBINDINGMODEL

	% Copyright: (C) 2008-2021 The CADET Authors
	%            See the license note at the end of the file.

	properties(Constant)
		name = 'MOBILE_PHASE_MODULATOR'; % Name of the binding model according to CADET file format specs
	end

	properties (Dependent, Transient)
		% Adsorption rate in [m^3_MP / (mol * s)]
		kA;
		MPM_KA;
		% Desorption rate in [1 / s]
		kD;
		MPM_KD;
		% Capacity in [mol / m^3_SP]
		qMax;
		MPM_QMAX;
		% Ion exchange characteristics [-]
		beta;
		MPM_BETA;
		% Hydrophobicity in [m^3_MP / mol]
		gamma;
		MPM_GAMMA;
	end

	methods

		function obj = MobilePhaseModulatorBinding(kA, kD, qMax)
			%MOBILEPHASEMODULATORSBINDING Constructs a MobilePhaseModulatorBinding object
			%   OBJ = MOBILEPHASEMODULATORSBINDING(KA) creates a MobilePhaseModulatorBinding model with
			%   the given adsorption rates KA.
			%
			%   OBJ = MOBILEPHASEMODULATORSBINDING(..., KD) also sets the desorption rates to KD.
			%
			%   OBJ = MOBILEPHASEMODULATORSBINDING(..., KD, QMAX) also sets the desorption rates
			%   to KD and the capacity to QMAX.

			obj = obj@KineticQuasiStationaryBindingModel();

			if nargin >= 1
				obj.kA = kA;
			end
			if nargin >= 2
				obj.kD = kD;
			end
			if nargin >= 3
				obj.qMax = qMax;
			end
		end

		function res = validate(obj, nComponents, nBoundStates)
			%VALIDATE Validates the binding model parameters
			%   RES = VALIDATE(NCOMPONENTS, NBOUNDSTATES) uses the number of
			%   components NCOMPONENTS of the model and the number of bound
			%   states of each of these components as given by the vector
			%   NBOUNDSTATES to validate the parameters of the binding model.
			%   Returns true in RES if everything is fine and false otherwise.

			validateattributes(obj.kA, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'kA');
			validateattributes(obj.kD, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'kD');
			validateattributes(obj.qMax, {'double'}, {'vector', 'nonempty', '>', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'qMax');
			validateattributes(obj.beta, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'beta');
			validateattributes(obj.gamma, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'gamma');
			res = obj.validate@KineticQuasiStationaryBindingModel(nComponents, nBoundStates);
		end

		function val = get.kA(obj)
			val = obj.data.MPM_KA;
		end

		function set.kA(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'kA');
			obj.data.MPM_KA = val;
			obj.hasChanged = true;
		end

		function val = get.kD(obj)
			val = obj.data.MPM_KD;
		end

		function set.kD(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'kD');
			obj.data.MPM_KD = val;
			obj.hasChanged = true;
		end

		function val = get.qMax(obj)
			val = obj.data.MPM_QMAX;
		end

		function set.qMax(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'qMax');
			obj.data.MPM_QMAX = val;
			obj.hasChanged = true;
		end

		function val = get.beta(obj)
			val = obj.data.MPM_BETA;
		end

		function set.beta(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'beta');
			obj.data.MPM_BETA = val;
			obj.hasChanged = true;
		end

		function val = get.gamma(obj)
			val = obj.data.MPM_GAMMA;
		end

		function set.gamma(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'gamma');
			obj.data.MPM_GAMMA = val;
			obj.hasChanged = true;
		end


		function val = get.MPM_KA(obj)
			val = obj.kA;
		end
		function set.MPM_KA(obj, val)
			obj.kA = val;
		end
		function val = get.MPM_KD(obj)
			val = obj.kD;
		end
		function set.MPM_KD(obj, val)
			obj.kD = val;
		end
		function val = get.MPM_QMAX(obj)
			val = obj.qMax;
		end
		function set.MPM_QMAX(obj, val)
			obj.qMax = val;
		end
		function val = get.MPM_BETA(obj)
			val = obj.beta;
		end
		function set.MPM_BETA(obj, val)
			obj.beta = val;
		end
		function val = get.MPM_GAMMA(obj)
			val = obj.gamma;
		end
		function set.MPM_GAMMA(obj, val)
			obj.gamma = val;
		end

	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = MobilePhaseModulatorBinding();
				obj.loadobjInternal(S);
			end
		end
		
	end

end

% =============================================================================
%  CADET
%  
%  Copyright (C) 2008-2021: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. obj program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies obj distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
