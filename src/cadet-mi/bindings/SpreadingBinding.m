
classdef SpreadingBinding < KineticQuasiStationaryBindingModel
	%SpreadingBinding Multi component spreading binding model
	%
	% See also BINDINGMODEL, KINETICQUASISTATIONARYBINDINGMODEL

	% Copyright: (C) 2008-2024 The CADET Authors
	%            See the license note at the end of the file.

	properties(Constant)
		name = 'MULTI_COMPONENT_SPREADING'; % Name of the binding model according to CADET file format specs
	end

	properties (Dependent, Transient)
		% Adsorption rate in [m^3_MP / (mol * s)]
		kA;
		MCSPR_KA;
		% Desorption rate in [1 / s]
		kD;
		MCSPR_KD;
		% Capacity in [mol / m^3_SP]
		qMax;
		MCSPR_QMAX;
		% Exchange rates from first to second bound state in [1 / s]
		k12;
		MCSPR_K12;
		% Exchange rates from second to first bound state in [1 / s]
		k21;
		MCSPR_K21;
	end

	methods

		function obj = SpreadingBinding(kA, kD)
			%SPREADINGBINDING Constructs a SpreadingBinding object
			%   OBJ = SPREADINGBINDING(KA) creates a SpreadingBinding model with the
			%   given adsorption rates KA.
			%
			%   OBJ = SPREADINGBINDING(..., KD) also sets the desorption rates to KD.

			obj = obj@KineticQuasiStationaryBindingModel();

			if nargin >= 1
				obj.kA = kA;
			end
			if nargin >= 2
				obj.kD = kD;
			end
		end

		function res = validate(obj, nComponents, nBoundStates)
			%VALIDATE Validates the binding model parameters
			%   RES = VALIDATE(NCOMPONENTS, NBOUNDSTATES) uses the number of
			%   components NCOMPONENTS of the model and the number of bound
			%   states of each of these components as given by the vector
			%   NBOUNDSTATES to validate the parameters of the binding model.
			%   Returns true in RES if everything is fine and false otherwise.

			if any((nBoundStates ~= 0) & (nBoundStates ~= 2))
				error('CADET:invalidParameter', 'SpreadingBinding requires each component to have either 0 or 2 bound states.');
			end

			validateattributes(obj.kA, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', 2 * nComponents}, '', 'kA');
			validateattributes(obj.kD, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', 2 * nComponents}, '', 'kD');
			validateattributes(obj.qMax, {'double'}, {'vector', 'nonempty', '>', 0.0, 'finite', 'real', 'numel', 2 * nComponents}, '', 'qMax');
			validateattributes(obj.k12, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'k12');
			validateattributes(obj.k21, {'double'}, {'scalar', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'k21');
			res = obj.validate@KineticQuasiStationaryBindingModel(nComponents, nBoundStates);
		end

		function val = get.kA(obj)
			val = obj.data.MCSPR_KA;
		end

		function set.kA(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'kA');
			obj.data.MCSPR_KA = val;
			obj.hasChanged = true;
		end

		function val = get.kD(obj)
			val = obj.data.MCSPR_KD;
		end

		function set.kD(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'kD');
			obj.data.MCSPR_KD = val;
			obj.hasChanged = true;
		end

		function val = get.qMax(obj)
			val = obj.data.MCSPR_QMAX;
		end

		function set.qMax(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'qMax');
			obj.data.MCSPR_QMAX = val;
			obj.hasChanged = true;
		end

		function val = get.k12(obj)
			val = obj.data.MCSPR_K12;
		end

		function set.k12(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'k12');
			obj.data.MCSPR_K12 = val;
			obj.hasChanged = true;
		end

		function val = get.k21(obj)
			val = obj.data.MCSPR_K21;
		end

		function set.k21(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'k21');
			obj.data.MCSPR_K21 = val;
			obj.hasChanged = true;
		end


		function val = get.MCSPR_KA(obj)
			val = obj.kA;
		end
		function set.MCSPR_KA(obj, val)
			obj.kA = val;
		end
		function val = get.MCSPR_KD(obj)
			val = obj.kD;
		end
		function set.MCSPR_KD(obj, val)
			obj.kD = val;
		end
		function val = get.MCSPR_QMAX(obj)
			val = obj.qMax;
		end
		function set.MCSPR_QMAX(obj, val)
			obj.qMax = val;
		end
		function val = get.MCSPR_K12(obj)
			val = obj.k12;
		end
		function set.MCSPR_K12(obj, val)
			obj.k12 = val;
		end
		function val = get.MCSPR_K21(obj)
			val = obj.k21;
		end
		function set.MCSPR_K21(obj, val)
			obj.k21 = val;
		end

	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = SpreadingBinding();
				obj.loadobjInternal(S);
			end
		end
		
	end
end

% =============================================================================
%  CADET
%  
%  Copyright (C) 2008-2024: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. obj program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies obj distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
