
classdef ExtFunSpreadingBinding < KineticQuasiStationaryBindingModel
	%ExtFunSpreadingBinding Multi component spreading binding model with external function support
	%
	% See also BINDINGMODEL, SPREADINGBINDING, KINETICQUASISTATIONARYBINDINGMODEL

	% Copyright: (C) 2008-2016 The CADET Authors
	%            See the license note at the end of the file.

	properties (Constant, Access = 'protected')
		hasConsistencySolver = true; % Determines whether this binding model has a consistency solver
	end

	properties(Constant)
		name = 'EXT_MULTI_COMPONENT_SPREADING'; % Name of the binding model according to CADET file format specs
	end

	properties (Dependent, Transient)
		% Adsorption rate in [m^3_MP / (mol * s)], constant term of external dependence
		kA;
		EXT_MCSPR_KA;
		% Adsorption rate in [m^3_MP / (mol * s * [T])], linear term of external dependence
		kA_T;
		EXT_MCSPR_KA_T;
		% Adsorption rate in [m^3_MP / (mol * s * [T]^2)], quadratic term of external dependence
		kA_TT;
		EXT_MCSPR_KA_TT;
		% Adsorption rate in [m^3_MP / (mol * s * [T]^3)], cubic term of external dependence
		kA_TTT;
		EXT_MCSPR_KA_TTT;
		% Desorption rate in [1 / s], constant term of external dependence
		kD;
		EXT_MCSPR_KD;
		% Desorption rate in [1 / (s * [T])], linear term of external dependence
		kD_T;
		EXT_MCSPR_KD_T;
		% Desorption rate in [1 / (s * [T]^2)], quadratic term of external dependence
		kD_TT;
		EXT_MCSPR_KD_TT;
		% Desorption rate in [1 / (s * [T]^3)], cubic term of external dependence
		kD_TTT;
		EXT_MCSPR_KD_TTT;
		% Capacity in [mol / m^3_SP], constant term of external dependence
		qMax;
		EXT_MCSPR_QMAX;
		% Capacity in [mol / (m^3_SP * [T])], linear term of external dependence
		qMax_T;
		EXT_MCSPR_QMAX_T;
		% Capacity in [mol / (m^3_SP * [T]^2)], quadratic term of external dependence
		qMax_TT;
		EXT_MCSPR_QMAX_TT;
		% Capacity in [mol / (m^3_SP * [T]^3)], cubic term of external dependence
		qMax_TTT;
		EXT_MCSPR_QMAX_TTT;
		% Exchange rates from first to second bound state in [1 / s], constant term of external dependence
		k12;
		EXT_MCSPR_K12;
		% Exchange rates from first to second bound state in [1 / (s * [T])], linear term of external dependence
		k12_T;
		EXT_MCSPR_K12_T;
		% Exchange rates from first to second bound state in [1 / (s * [T]^2)], quadratic term of external dependence
		k12_TT;
		EXT_MCSPR_K12_TT;
		% Exchange rates from first to second bound state in [1 / (s * [T]^3)], cubic term of external dependence
		k12_TTT;
		EXT_MCSPR_K12_TTT;
		% Exchange rates from second to first bound state in [1 / s], constant term of external dependence
		k21;
		EXT_MCSPR_K21;
		% Exchange rates from second to first bound state in [1 / (s * [T])], linear term of external dependence
		k21_T;
		EXT_MCSPR_K21_T;
		% Exchange rates from second to first bound state in [1 / (s * [T]^2)], quadratic term of external dependence
		k21_TT;
		EXT_MCSPR_K21_TT;
		% Exchange rates from second to first bound state in [1 / (s * [T]^3)], cubic term of external dependence
		k21_TTT;
		EXT_MCSPR_K21_TTT;
		% Indices of external functions (0-based)
		externalSource;
		EXTFUN;
	end

	methods

		function obj = ExtFunSpreadingBinding(kA, kD)
			%EXTFUNSPREADINGBINDING Constructs an ExtFunSpreadingBinding object with external function support
			%   OBJ = EXTFUNSPREADINGBINDING(KA) creates an ExtFunSpreadingBinding model
			%   with the given adsorption rates KA.
			%
			%   OBJ = EXTFUNSPREADINGBINDING(..., KD) also sets the desorption rates to KD.

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

			validateattributes(obj.kA, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', sum(nBoundStates)}, '', 'kA');
			validateattributes(obj.kD, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', sum(nBoundStates)}, '', 'kD');
			validateattributes(obj.qMax, {'double'}, {'vector', 'nonempty', '>', 0.0, 'finite', 'real', 'numel', sum(nBoundStates)}, '', 'qMax');
			validateattributes(obj.k12, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'k12');
			validateattributes(obj.k21, {'double'}, {'scalar', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'k21');

			validateattributes(obj.kA_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', sum(nBoundStates)}, '', 'kA_T');
			validateattributes(obj.kD_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', sum(nBoundStates)}, '', 'kD_T');
			validateattributes(obj.qMax_T, {'double'}, {'vector', 'nonempty', '>', 0.0, 'finite', 'real', 'numel', sum(nBoundStates)}, '', 'qMax_T');
			validateattributes(obj.k12_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'k12_T');
			validateattributes(obj.k21_T, {'double'}, {'scalar', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'k21_T');

			validateattributes(obj.kA_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', sum(nBoundStates)}, '', 'kA_TT');
			validateattributes(obj.kD_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', sum(nBoundStates)}, '', 'kD_TT');
			validateattributes(obj.qMax_TT, {'double'}, {'vector', 'nonempty', '>', 0.0, 'finite', 'real', 'numel', sum(nBoundStates)}, '', 'qMax_TT');
			validateattributes(obj.k12_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'k12_TT');
			validateattributes(obj.k21_TT, {'double'}, {'scalar', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'k21_TT');

			validateattributes(obj.kA_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', sum(nBoundStates)}, '', 'kA_TTT');
			validateattributes(obj.kD_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', sum(nBoundStates)}, '', 'kD_TTT');
			validateattributes(obj.qMax_TTT, {'double'}, {'vector', 'nonempty', '>', 0.0, 'finite', 'real', 'numel', sum(nBoundStates)}, '', 'qMax_TTT');
			validateattributes(obj.k12_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'k12_TTT');
			validateattributes(obj.k21_TTT, {'double'}, {'scalar', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'k21_TTT');

			if ~isempty(obj.externalSource)
				validateattributes(obj.externalSource, {'double', 'int32'}, {'vector', 'nonempty', 'finite', 'real', '>=', 0}, '', 'externalSource');
				if ~any(numel(obj.externalSource) == [1, 5])
					error('CADET:invalidConfig', 'Expected externalSource to have either 1 or 5 entries but got %d entries.', numel(obj.externalSource));
				end
			end

			res = obj.validate@KineticQuasiStationaryBindingModel(nComponents, nBoundStates);
		end

		function val = get.kA(obj)
			val = obj.data.EXT_MCSPR_KA;
		end

		function set.kA(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kA');
			obj.data.EXT_MCSPR_KA = val;
			obj.hasChanged = true;
		end

		function val = get.kD(obj)
			val = obj.data.EXT_MCSPR_KD;
		end

		function set.kD(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kD');
			obj.data.EXT_MCSPR_KD = val;
			obj.hasChanged = true;
		end

		function val = get.qMax(obj)
			val = obj.data.EXT_MCSPR_QMAX;
		end

		function set.qMax(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'qMax');
			obj.data.EXT_MCSPR_QMAX = val;
			obj.hasChanged = true;
		end

		function val = get.k12(obj)
			val = obj.data.EXT_MCSPR_K12;
		end

		function set.k12(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'k12');
			obj.data.EXT_MCSPR_K12 = val;
			obj.hasChanged = true;
		end

		function val = get.k21(obj)
			val = obj.data.EXT_MCSPR_K21;
		end

		function set.k21(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'k21');
			obj.data.EXT_MCSPR_K21 = val;
			obj.hasChanged = true;
		end

		function val = get.kA_T(obj)
			val = obj.data.EXT_MCSPR_KA_T;
		end

		function set.kA_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kA_T');
			obj.data.EXT_MCSPR_KA_T = val;
			obj.hasChanged = true;
		end

		function val = get.kD_T(obj)
			val = obj.data.EXT_MCSPR_KD_T;
		end

		function set.kD_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kD_T');
			obj.data.EXT_MCSPR_KD_T = val;
			obj.hasChanged = true;
		end

		function val = get.qMax_T(obj)
			val = obj.data.EXT_MCSPR_QMAX_T;
		end

		function set.qMax_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'qMax_T');
			obj.data.EXT_MCSPR_QMAX_T = val;
			obj.hasChanged = true;
		end

		function val = get.k12_T(obj)
			val = obj.data.EXT_MCSPR_K12_T;
		end

		function set.k12_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'k12_T');
			obj.data.EXT_MCSPR_K12_T = val;
			obj.hasChanged = true;
		end

		function val = get.k21_T(obj)
			val = obj.data.EXT_MCSPR_K21_T;
		end

		function set.k21_T(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'k21_T');
			obj.data.EXT_MCSPR_K21_T = val;
			obj.hasChanged = true;
		end

		function val = get.kA_TT(obj)
			val = obj.data.EXT_MCSPR_KA_TT;
		end

		function set.kA_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kA_TT');
			obj.data.EXT_MCSPR_KA_TT = val;
			obj.hasChanged = true;
		end

		function val = get.kD_TT(obj)
			val = obj.data.EXT_MCSPR_KD_TT;
		end

		function set.kD_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kD_TT');
			obj.data.EXT_MCSPR_KD_TT = val;
			obj.hasChanged = true;
		end

		function val = get.qMax_TT(obj)
			val = obj.data.EXT_MCSPR_QMAX_TT;
		end

		function set.qMax_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'qMax_TT');
			obj.data.EXT_MCSPR_QMAX_TT = val;
			obj.hasChanged = true;
		end

		function val = get.k12_TT(obj)
			val = obj.data.EXT_MCSPR_K12_TT;
		end

		function set.k12_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'k12_TT');
			obj.data.EXT_MCSPR_K12_TT = val;
			obj.hasChanged = true;
		end

		function val = get.k21_TT(obj)
			val = obj.data.EXT_MCSPR_K21_TT;
		end

		function set.k21_TT(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'k21_TT');
			obj.data.EXT_MCSPR_K21_TT = val;
			obj.hasChanged = true;
		end

		function val = get.kA_TTT(obj)
			val = obj.data.EXT_MCSPR_KA_TTT;
		end

		function set.kA_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kA_TTT');
			obj.data.EXT_MCSPR_KA_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.kD_TTT(obj)
			val = obj.data.EXT_MCSPR_KD_TTT;
		end

		function set.kD_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kD_TTT');
			obj.data.EXT_MCSPR_KD_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.qMax_TTT(obj)
			val = obj.data.EXT_MCSPR_QMAX_TTT;
		end

		function set.qMax_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'qMax_TTT');
			obj.data.EXT_MCSPR_QMAX_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.k12_TTT(obj)
			val = obj.data.EXT_MCSPR_K12_TTT;
		end

		function set.k12_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'k12_TTT');
			obj.data.EXT_MCSPR_K12_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.k21_TTT(obj)
			val = obj.data.EXT_MCSPR_K21_TTT;
		end

		function set.k21_TTT(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'k21_TTT');
			obj.data.EXT_MCSPR_K21_TTT = val;
			obj.hasChanged = true;
		end


		function val = get.EXT_MCSPR_KA(obj)
			val = obj.kA;
		end
		function set.EXT_MCSPR_KA(obj, val)
			obj.kA = val;
		end
		function val = get.EXT_MCSPR_KD(obj)
			val = obj.kD;
		end
		function set.EXT_MCSPR_KD(obj, val)
			obj.kD = val;
		end
		function val = get.EXT_MCSPR_QMAX(obj)
			val = obj.qMax;
		end
		function set.EXT_MCSPR_QMAX(obj, val)
			obj.qMax = val;
		end
		function val = get.EXT_MCSPR_K12(obj)
			val = obj.k12;
		end
		function set.EXT_MCSPR_K12(obj, val)
			obj.k12 = val;
		end
		function val = get.EXT_MCSPR_K21(obj)
			val = obj.k21;
		end
		function set.EXT_MCSPR_K21(obj, val)
			obj.k21 = val;
		end

		function val = get.EXT_MCSPR_KA_T(obj)
			val = obj.kA_T;
		end
		function set.EXT_MCSPR_KA_T(obj, val)
			obj.kA_T = val;
		end
		function val = get.EXT_MCSPR_KD_T(obj)
			val = obj.kD_T;
		end
		function set.EXT_MCSPR_KD_T(obj, val)
			obj.kD_T = val;
		end
		function val = get.EXT_MCSPR_QMAX_T(obj)
			val = obj.qMax_T;
		end
		function set.EXT_MCSPR_QMAX_T(obj, val)
			obj.qMax_T = val;
		end
		function val = get.EXT_MCSPR_K12_T(obj)
			val = obj.k12_T;
		end
		function set.EXT_MCSPR_K12_T(obj, val)
			obj.k12_T = val;
		end
		function val = get.EXT_MCSPR_K21_T(obj)
			val = obj.k21_T;
		end
		function set.EXT_MCSPR_K21_T(obj, val)
			obj.k21_T = val;
		end

		function val = get.EXT_MCSPR_KA_TT(obj)
			val = obj.kA_TT;
		end
		function set.EXT_MCSPR_KA_TT(obj, val)
			obj.kA_TT = val;
		end
		function val = get.EXT_MCSPR_KD_TT(obj)
			val = obj.kD_TT;
		end
		function set.EXT_MCSPR_KD_TT(obj, val)
			obj.kD_TT = val;
		end
		function val = get.EXT_MCSPR_QMAX_TT(obj)
			val = obj.qMax_TT;
		end
		function set.EXT_MCSPR_QMAX_TT(obj, val)
			obj.qMax_TT = val;
		end
		function val = get.EXT_MCSPR_K12_TT(obj)
			val = obj.k12_TT;
		end
		function set.EXT_MCSPR_K12_TT(obj, val)
			obj.k12_TT = val;
		end
		function val = get.EXT_MCSPR_K21_TT(obj)
			val = obj.k21_TT;
		end
		function set.EXT_MCSPR_K21_TT(obj, val)
			obj.k21_TT = val;
		end

		function val = get.EXT_MCSPR_KA_TTT(obj)
			val = obj.kA_TTT;
		end
		function set.EXT_MCSPR_KA_TTT(obj, val)
			obj.kA_TTT = val;
		end
		function val = get.EXT_MCSPR_KD_TTT(obj)
			val = obj.kD_TTT;
		end
		function set.EXT_MCSPR_KD_TTT(obj, val)
			obj.kD_TTT = val;
		end
		function val = get.EXT_MCSPR_QMAX_TTT(obj)
			val = obj.qMax_TTT;
		end
		function set.EXT_MCSPR_QMAX_TTT(obj, val)
			obj.qMax_TTT = val;
		end
		function val = get.EXT_MCSPR_K12_TTT(obj)
			val = obj.k12_TTT;
		end
		function set.EXT_MCSPR_K12_TTT(obj, val)
			obj.k12_TTT = val;
		end
		function val = get.EXT_MCSPR_K21_TTT(obj)
			val = obj.k21_TTT;
		end
		function set.EXT_MCSPR_K21_TTT(obj, val)
			obj.k21_TTT = val;
		end

		function val = get.externalSource(obj)
			if isfield(obj.data, 'EXTFUN')
				val = double(obj.data.EXTFUN);
			else
				val = [];
			end
		end

		function set.externalSource(obj, val)
			obj.hasChanged = true;
			if isempty(val)
				obj.data = rmfield(obj.data, 'EXTFUN');
			else
				validateattributes(val, {'double', 'int32'}, {'vector', 'nonempty', 'finite', 'real', '>=', 0}, '', 'externalSource');
				if ~any(numel(val) == [1, 5])
					error('CADET:invalidConfig', 'Expected externalSource to have either 1 or 5 entries but got %d entries.', numel(val));
				end
				obj.data.EXTFUN = int32(val);
			end
		end

		function val = get.EXTFUN(obj)
			val = obj.externalSource;
		end
		function set.EXTFUN(obj, val)
			obj.externalSource = val;
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
			%   This implementation assumes component-state-major ordering for EXT_MCSPR_KA, EXT_MCSPR_KD, and EXT_MCSPR_QMAX
			%   and section-boundphase-component-major ordering for all other parameters.
			%
			% See also BINDINGMODEL.OFFSETTOPARAMETER

			switch param.SENS_NAME
				case {'EXT_MCSPR_KA', 'EXT_MCSPR_KD', 'EXT_MCSPR_QMAX', ...
					'EXT_MCSPR_KA_T', 'EXT_MCSPR_KD_T', 'EXT_MCSPR_QMAX_T', ...
					'EXT_MCSPR_KA_TT', 'EXT_MCSPR_KD_TT', 'EXT_MCSPR_QMAX_TT', ...
					'EXT_MCSPR_KA_TTT', 'EXT_MCSPR_KD_TTT', 'EXT_MCSPR_QMAX_TTT'}
					offset = 0;
					if (param.SENS_BOUNDPHASE ~= -1)
						offset = offset + param.SENS_BOUNDPHASE;
					end
					if (param.SENS_COMP ~= -1)
						offset = offset + sum(nBoundStates(1:param.SENS_COMP));
					end
				otherwise
					offset = obj.offsetToParameter@KineticQuasiStationaryBindingModel(nBoundStates, param)
			end
		end

	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = ExtFunSpreadingBinding();
				obj.loadobjInternal(S);
			end
		end
		
	end
end

% =============================================================================
%  CADET - The Chromatography Analysis and Design Toolkit
%  
%  Copyright (C) 2008-2016: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. obj program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies obj distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
