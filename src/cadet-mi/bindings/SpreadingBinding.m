
classdef SpreadingBinding < KineticQuasiStationaryBindingModel
	%SpreadingBinding Multi component spreading binding model
	%
	% See also BINDINGMODEL, KINETICQUASISTATIONARYBINDINGMODEL

	% Copyright: (C) 2008-2020 The CADET Authors
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

			validateattributes(obj.kA, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', sum(nBoundStates)}, '', 'kA');
			validateattributes(obj.kD, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', sum(nBoundStates)}, '', 'kD');
			validateattributes(obj.qMax, {'double'}, {'vector', 'nonempty', '>', 0.0, 'finite', 'real', 'numel', sum(nBoundStates)}, '', 'qMax');
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

	methods (Access = 'protected')

		function offset = offsetToParameter(obj, nBoundStates, param)
			%OFFSETTOPARAMETER Computes the (zero-based) offset to the given parameter in a linearized array
			%   OFFSET = OFFSETTOPARAMETER(NBOUNDSTATES, PARAM) uses the number of bound states of each
			%   component given in the vector NBOUNDSTATES and the parameter struct PARAM to compute the
			%   0-based offset of the parameter in a linearized array. The fields SENS_COMP, SENS_BOUNDPHASE,
			%   SENS_REACTION, and SENS_SECTION of PARAM are used to calculate the offset.
			%
			%   This implementation assumes component-state-major ordering for MCSPR_KA, MCSPR_KD, and MCSPR_QMAX
			%   and section-boundphase-component-major ordering for all other parameters.
			%
			% See also BINDINGMODEL.OFFSETTOPARAMETER

			switch param.SENS_NAME
				case {'MCSPR_KA', 'MCSPR_KD', 'MCSPR_QMAX'}
					offset = 0;
					if (param.SENS_BOUNDPHASE ~= -1)
						offset = offset + param.SENS_BOUNDPHASE;
					end
					if (param.SENS_COMP ~= -1)
						offset = offset + sum(nBoundStates(1:param.SENS_COMP));
					end
				otherwise
					offset = obj.offsetToParameter@KineticQuasiStationaryBindingModel(nBoundStates, param);
			end
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
%  CADET - The Chromatography Analysis and Design Toolkit
%  
%  Copyright (C) 2008-2020: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. obj program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies obj distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
