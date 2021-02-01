
classdef ExtendedMobilePhaseModulatorBinding < KineticQuasiStationaryBindingModel
	%ExtendedMobilePhaseModulatorBinding Extended mobile phase modulator binding model
	%
	% See also BINDINGMODEL, KINETICQUASISTATIONARYBINDINGMODEL

	% Copyright: (C) 2008-2021 The CADET Authors
	%            See the license note at the end of the file.

	properties(Constant)
		name = 'EXTENDED_MOBILE_PHASE_MODULATOR'; % Name of the binding model according to CADET file format specs
	end

	properties (Dependent, Transient)
		% Component binding mode
		compMode;
		EMPM_COMP_MODE;
		% Adsorption rate in [m^3_MP / (mol * s)]
		kA;
		EMPM_KA;
		% Desorption rate in [1 / s]
		kD;
		EMPM_KD;
		% Capacity in [mol / m^3_SP]
		qMax;
		EMPM_QMAX;
		% Ion exchange characteristics [-]
		beta;
		EMPM_BETA;
		% Hydrophobicity in [m^3_MP / mol]
		gamma;
		EMPM_GAMMA;
	end

	methods

		function obj = ExtendedMobilePhaseModulatorBinding(kA, kD, qMax)
			%EXTENDEDMOBILEPHASEMODULATORSBINDING Constructs an ExtendedMobilePhaseModulatorBinding object
			%   OBJ = EXTENDEDMOBILEPHASEMODULATORSBINDING(KA) creates an ExtendedMobilePhaseModulatorBinding model with
			%   the given adsorption rates KA.
			%
			%   OBJ = EXTENDEDMOBILEPHASEMODULATORSBINDING(..., KD) also sets the desorption rates to KD.
			%
			%   OBJ = EXTENDEDMOBILEPHASEMODULATORSBINDING(..., KD, QMAX) also sets the desorption rates
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

			validateattributes(obj.compMode, {'double'}, {'vector', 'nonempty', '>=', 0.0, '<=', 2.0, 'finite', 'real', 'numel', nComponents}, '', 'compMode');
			validateattributes(obj.kA, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'kA');
			validateattributes(obj.kD, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'kD');
			validateattributes(obj.qMax, {'double'}, {'vector', 'nonempty', '>', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'qMax');
			validateattributes(obj.beta, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'beta');
			validateattributes(obj.gamma, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'gamma');
			res = obj.validate@KineticQuasiStationaryBindingModel(nComponents, nBoundStates);
		end

		function val = get.compMode(obj)
			val = double(obj.data.EMPM_COMP_MODE);
		end

		function set.compMode(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, '<=', 2.0, 'finite', 'real'}, '', 'compMode');
			obj.data.EMPM_COMP_MODE = int32(val);
			obj.hasChanged = true;
		end

		function val = get.kA(obj)
			val = obj.data.EMPM_KA;
		end

		function set.kA(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'kA');
			obj.data.EMPM_KA = val;
			obj.hasChanged = true;
		end

		function val = get.kD(obj)
			val = obj.data.EMPM_KD;
		end

		function set.kD(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'kD');
			obj.data.EMPM_KD = val;
			obj.hasChanged = true;
		end

		function val = get.qMax(obj)
			val = obj.data.EMPM_QMAX;
		end

		function set.qMax(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'qMax');
			obj.data.EMPM_QMAX = val;
			obj.hasChanged = true;
		end

		function val = get.beta(obj)
			val = obj.data.EMPM_BETA;
		end

		function set.beta(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'beta');
			obj.data.EMPM_BETA = val;
			obj.hasChanged = true;
		end

		function val = get.gamma(obj)
			val = obj.data.EMPM_GAMMA;
		end

		function set.gamma(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'gamma');
			obj.data.EMPM_GAMMA = val;
			obj.hasChanged = true;
		end


		function val = get.EMPM_COMP_MODE(obj)
			val = obj.compMode;
		end
		function set.EMPM_COMP_MODE(obj, val)
			obj.compMode = val;
		end
		function val = get.EMPM_KA(obj)
			val = obj.kA;
		end
		function set.EMPM_KA(obj, val)
			obj.kA = val;
		end
		function val = get.EMPM_KD(obj)
			val = obj.kD;
		end
		function set.EMPM_KD(obj, val)
			obj.kD = val;
		end
		function val = get.EMPM_QMAX(obj)
			val = obj.qMax;
		end
		function set.EMPM_QMAX(obj, val)
			obj.qMax = val;
		end
		function val = get.EMPM_BETA(obj)
			val = obj.beta;
		end
		function set.EMPM_BETA(obj, val)
			obj.beta = val;
		end
		function val = get.EMPM_GAMMA(obj)
			val = obj.gamma;
		end
		function set.EMPM_GAMMA(obj, val)
			obj.gamma = val;
		end

	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = ExtendedMobilePhaseModulatorBinding();
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
