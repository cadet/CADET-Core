
classdef ExtFunAntiLangmuirBinding < KineticQuasiStationaryBindingModel
	%ExtFunAntiLangmuirBinding Anti-Langmuir binding model with external function support
	%
	% See also BINDINGMODEL, ANTILANGMUIRBINDING, KINETICQUASISTATIONARYBINDINGMODEL

	% Copyright: (C) 2008-2020 The CADET Authors
	%            See the license note at the end of the file.

	properties(Constant)
		name = 'EXT_MULTI_COMPONENT_ANTILANGMUIR'; % Name of the binding model according to CADET file format specs
	end

	properties (Dependent, Transient)
		% Adsorption rate in [m^3_MP / (mol * s)], constant term of external dependence
		kA;
		EXT_MCAL_KA;
		% Adsorption rate in [m^3_MP / (mol * s * [T])], linear term of external dependence
		kA_T;
		EXT_MCAL_KA_T;
		% Adsorption rate in [m^3_MP / (mol * s * [T]^2)], quadratic term of external dependence
		kA_TT;
		EXT_MCAL_KA_TT;
		% Adsorption rate in [m^3_MP / (mol * s * [T]^3)], cubic term of external dependence
		kA_TTT;
		EXT_MCAL_KA_TTT;
		% Desorption rate in [1 / s], constant term of external dependence
		kD;
		EXT_MCAL_KD;
		% Desorption rate in [1 / (s * [T])], linear term of external dependence
		kD_T;
		EXT_MCAL_KD_T;
		% Desorption rate in [1 / (s * [T]^2)], quadratic term of external dependence
		kD_TT;
		EXT_MCAL_KD_TT;
		% Desorption rate in [1 / (s * [T]^3)], cubic term of external dependence
		kD_TTT;
		EXT_MCAL_KD_TTT;
		% Capacity in [mol / m^3_SP], constant term of external dependence
		qMax;
		EXT_MCAL_QMAX;
		% Capacity in [mol / (m^3_SP * [T])], linear term of external dependence
		qMax_T;
		EXT_MCAL_QMAX_T;
		% Capacity in [mol / (m^3_SP * [T]^2)], quadratic term of external dependence
		qMax_TT;
		EXT_MCAL_QMAX_TT;
		% Capacity in [mol / (m^3_SP * [T]^3)], cubic term of external dependence
		qMax_TTT;
		EXT_MCAL_QMAX_TTT;
		% Anti-Langmuir coefficient [-], constant term of external dependence
		antiFactor;
		EXT_MCAL_ANTILANGMUIR;
		% Anti-Langmuir coefficient [1 / [T]], linear term of external dependence
		antiFactor_T;
		EXT_MCAL_ANTILANGMUIR_T;
		% Anti-Langmuir coefficient [1 / [T]^2], quadratic term of external dependence
		antiFactor_TT;
		EXT_MCAL_ANTILANGMUIR_TT;
		% Anti-Langmuir coefficient [1 / [T]^3], cubic term of external dependence
		antiFactor_TTT;
		EXT_MCAL_ANTILANGMUIR_TTT;
		% Indices of external functions (0-based)
		externalSource;
		EXTFUN;
	end

	methods

		function obj = ExtFunAntiLangmuirBinding(kA, kD, qMax, antiFactor)
			%EXTFUNANTILANGMUIRBINDING Constructs an Anti-Langmuir binding model object with external function support
			%   OBJ = EXTFUNANTILANGMUIRBINDING(KA) creates an Anti-Langmuir binding model with
			%   the given adsorption rates KA.
			%
			%   OBJ = EXTFUNANTILANGMUIRBINDING(..., KD) also sets the desorption rates to KD.
			%
			%   OBJ = EXTFUNANTILANGMUIRBINDING(..., KD, QMAX) also sets the desorption rates to KD
			%   and the capacity to QMAX.
			%
			%   OBJ = EXTFUNANTILANGMUIRBINDING(..., KD, QMAX, ANTIFACTOR) also sets the desorption
			%   rates to KD, the capacity to QMAX, and the anti-Langmuir coefficients to ANTIFACTOR.

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
			if nargin >= 4
				obj.antiFactor = antiFactor;
			end
		end

		function res = validate(obj, nComponents, nBoundStates)
			%VALIDATE Validates the binding model parameters
			%   RES = VALIDATE(NCOMPONENTS, NBOUNDSTATES) uses the number of
			%   components NCOMPONENTS of the model and the number of bound
			%   states of each of these components as given by the vector
			%   NBOUNDSTATES to validate the parameters of the binding model.
			%   Returns true in RES if everything is fine and false otherwise.

			validateattributes(obj.kA, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kA');
			validateattributes(obj.kD, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kD');
			validateattributes(obj.qMax, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'qMax');
			validateattributes(obj.antiFactor, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'antiFactor');

			validateattributes(obj.kA_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kA_T');
			validateattributes(obj.kD_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kD_T');
			validateattributes(obj.qMax_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'qMax_T');
			validateattributes(obj.antiFactor_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'antiFactor_T');

			validateattributes(obj.kA_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kA_TT');
			validateattributes(obj.kD_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kD_TT');
			validateattributes(obj.qMax_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'qMax_TT');
			validateattributes(obj.antiFactor_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'antiFactor_TT');

			validateattributes(obj.kA_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kA_TTT');
			validateattributes(obj.kD_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kD_TTT');
			validateattributes(obj.qMax_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'qMax_TTT');
			validateattributes(obj.antiFactor_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'antiFactor_TTT');

			if ~isempty(obj.externalSource)
				validateattributes(obj.externalSource, {'double', 'int32'}, {'vector', 'nonempty', 'finite', 'real', '>=', 0}, '', 'externalSource');
				if ~any(numel(obj.externalSource) == [1, 4])
					error('CADET:invalidConfig', 'Expected externalSource to have either 1 or 4 entries but got %d entries.', numel(obj.externalSource));
				end
			end

			res = obj.validate@KineticQuasiStationaryBindingModel(nComponents, nBoundStates);
		end

		function val = get.kA(obj)
			val = obj.data.EXT_MCAL_KA;
		end

		function set.kA(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kA');
			obj.data.EXT_MCAL_KA = val;
			obj.hasChanged = true;
		end

		function val = get.kD(obj)
			val = obj.data.EXT_MCAL_KD;
		end

		function set.kD(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kD');
			obj.data.EXT_MCAL_KD = val;
			obj.hasChanged = true;
		end

		function val = get.qMax(obj)
			val = obj.data.EXT_MCAL_QMAX;
		end

		function set.qMax(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'qMax');
			obj.data.EXT_MCAL_QMAX = val;
			obj.hasChanged = true;
		end

		function val = get.antiFactor(obj)
			val = obj.data.EXT_MCAL_ANTILANGMUIR;
		end

		function set.antiFactor(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'antiFactor');
			obj.data.EXT_MCAL_ANTILANGMUIR = val;
			obj.hasChanged = true;
		end

		function val = get.kA_T(obj)
			val = obj.data.EXT_MCAL_KA_T;
		end

		function set.kA_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kA_T');
			obj.data.EXT_MCAL_KA_T = val;
			obj.hasChanged = true;
		end

		function val = get.kD_T(obj)
			val = obj.data.EXT_MCAL_KD_T;
		end

		function set.kD_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kD_T');
			obj.data.EXT_MCAL_KD_T = val;
			obj.hasChanged = true;
		end

		function val = get.qMax_T(obj)
			val = obj.data.EXT_MCAL_QMAX;
		end

		function set.qMax_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'qMax_T');
			obj.data.EXT_MCAL_QMAX_T = val;
			obj.hasChanged = true;
		end

		function val = get.antiFactor_T(obj)
			val = obj.data.EXT_MCAL_ANTILANGMUIR_T;
		end

		function set.antiFactor_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'antiFactor_T');
			obj.data.EXT_MCAL_ANTILANGMUIR_T = val;
			obj.hasChanged = true;
		end

		function val = get.kA_TT(obj)
			val = obj.data.EXT_MCAL_KA_TT;
		end

		function set.kA_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kA_TT');
			obj.data.EXT_MCAL_KA_TT = val;
			obj.hasChanged = true;
		end

		function val = get.kD_TT(obj)
			val = obj.data.EXT_MCAL_KD_TT;
		end

		function set.kD_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kD_TT');
			obj.data.EXT_MCAL_KD_TT = val;
			obj.hasChanged = true;
		end

		function val = get.qMax_TT(obj)
			val = obj.data.EXT_MCAL_QMAX_TT;
		end

		function set.qMax_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'qMax_TT');
			obj.data.EXT_MCAL_QMAX_TT = val;
			obj.hasChanged = true;
		end

		function val = get.antiFactor_TT(obj)
			val = obj.data.EXT_MCAL_ANTILANGMUIR_TT;
		end

		function set.antiFactor_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'antiFactor_TT');
			obj.data.EXT_MCAL_ANTILANGMUIR_TT = val;
			obj.hasChanged = true;
		end

		function val = get.kA_TTT(obj)
			val = obj.data.EXT_MCAL_KA_TTT;
		end

		function set.kA_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kA_TTT');
			obj.data.EXT_MCAL_KA_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.kD_TTT(obj)
			val = obj.data.EXT_MCAL_KD_TTT;
		end

		function set.kD_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kD_TTT');
			obj.data.EXT_MCAL_KD_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.qMax_TTT(obj)
			val = obj.data.EXT_MCAL_QMAX_TTT;
		end

		function set.qMax_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'qMax_TTT');
			obj.data.EXT_MCAL_QMAX_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.antiFactor_TTT(obj)
			val = obj.data.EXT_MCAL_ANTILANGMUIR_TTT;
		end

		function set.antiFactor_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'antiFactor_TTT');
			obj.data.EXT_MCAL_ANTILANGMUIR_TTT = val;
			obj.hasChanged = true;
		end


		function val = get.EXT_MCAL_KA(obj)
			val = obj.kA;
		end
		function set.EXT_MCAL_KA(obj, val)
			obj.kA = val;
		end
		function val = get.EXT_MCAL_KD(obj)
			val = obj.kD;
		end
		function set.EXT_MCAL_KD(obj, val)
			obj.kD = val;
		end
		function val = get.EXT_MCAL_QMAX(obj)
			val = obj.qMax;
		end
		function set.EXT_MCAL_QMAX(obj, val)
			obj.qMax = val;
		end
		function val = get.EXT_MCAL_ANTILANGMUIR(obj)
			val = obj.antiFactor;
		end
		function set.EXT_MCAL_ANTILANGMUIR(obj, val)
			obj.antiFactor = val;
		end

		function val = get.EXT_MCAL_KA_T(obj)
			val = obj.kA_T;
		end
		function set.EXT_MCAL_KA_T(obj, val)
			obj.kA_T = val;
		end
		function val = get.EXT_MCAL_KD_T(obj)
			val = obj.kD_T;
		end
		function set.EXT_MCAL_KD_T(obj, val)
			obj.kD_T = val;
		end
		function val = get.EXT_MCAL_QMAX_T(obj)
			val = obj.qMax_T;
		end
		function set.EXT_MCAL_QMAX_T(obj, val)
			obj.qMax_T = val;
		end
		function val = get.EXT_MCAL_ANTILANGMUIR_T(obj)
			val = obj.antiFactor_T;
		end
		function set.EXT_MCAL_ANTILANGMUIR_T(obj, val)
			obj.antiFactor_T = val;
		end

		function val = get.EXT_MCAL_KA_TT(obj)
			val = obj.kA_TT;
		end
		function set.EXT_MCAL_KA_TT(obj, val)
			obj.kA_TT = val;
		end
		function val = get.EXT_MCAL_KD_TT(obj)
			val = obj.kD_TT;
		end
		function set.EXT_MCAL_KD_TT(obj, val)
			obj.kD_TT = val;
		end
		function val = get.EXT_MCAL_QMAX_TT(obj)
			val = obj.qMax_TT;
		end
		function set.EXT_MCAL_QMAX_TT(obj, val)
			obj.qMax_TT = val;
		end
		function val = get.EXT_MCAL_ANTILANGMUIR_TT(obj)
			val = obj.antiFactor_TT;
		end
		function set.EXT_MCAL_ANTILANGMUIR_TT(obj, val)
			obj.antiFactor_TT = val;
		end

		function val = get.EXT_MCAL_KA_TTT(obj)
			val = obj.kA_TTT;
		end
		function set.EXT_MCAL_KA_TTT(obj, val)
			obj.kA_TTT = val;
		end
		function val = get.EXT_MCAL_KD_TTT(obj)
			val = obj.kD_TTT;
		end
		function set.EXT_MCAL_KD_TTT(obj, val)
			obj.kD_TTT = val;
		end
		function val = get.EXT_MCAL_QMAX_TTT(obj)
			val = obj.qMax_TTT;
		end
		function set.EXT_MCAL_QMAX_TTT(obj, val)
			obj.qMax_TTT = val;
		end
		function val = get.EXT_MCAL_ANTILANGMUIR_TTT(obj)
			val = obj.antiFactor_TTT;
		end
		function set.EXT_MCAL_ANTILANGMUIR_TTT(obj, val)
			obj.antiFactor_TTT = val;
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
				if ~any(numel(val) == [1, 4])
					error('CADET:invalidConfig', 'Expected externalSource to have either 1 or 4 entries but got %d entries.', numel(val));
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

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = ExtFunAntiLangmuirBinding();
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
