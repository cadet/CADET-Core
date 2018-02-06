
classdef ExtFunLinearBinding < KineticQuasiStationaryBindingModel
	%ExtFunLinearBinding Linear binding model with external function support
	%
	% See also BINDINGMODEL, LINEARBINDING, KINETICQUASISTATIONARYBINDINGMODEL

	% Copyright: (C) 2008-2018 The CADET Authors
	%            See the license note at the end of the file.

	properties (Constant, Access = 'protected')
		hasConsistencySolver = true; % Determines whether this binding model has a consistency solver
	end

	properties(Constant)
		name = 'EXT_LINEAR'; % Name of the binding model according to CADET file format specs
	end

	properties (Dependent, Transient)
		% Adsorption rate in [m^3_MP / (m^3_SP * s)], constant term of external dependence
		kA;
		EXT_LIN_KA;
		% Adsorption rate in [m^3_MP / (m^3_SP * s * [T])], linear term of external dependence
		kA_T;
		EXT_LIN_KA_T;
		% Adsorption rate in [m^3_MP / (m^3_SP * s * [T]^2)], quadratic term of external dependence
		kA_TT;
		EXT_LIN_KA_TT;
		% Adsorption rate in [m^3_MP / (m^3_SP * s * [T]^3)], cubic term of external dependence
		kA_TTT;
		EXT_LIN_KA_TTT;
		% Desorption rate in [1 / s], constant term of external dependence
		kD;
		EXT_LIN_KD;
		% Desorption rate in [1 / (s * [T])], linear term of external dependence
		kD_T;
		EXT_LIN_KD_T;
		% Desorption rate in [1 / (s * [T]^2)], quadratic term of external dependence
		kD_TT;
		EXT_LIN_KD_TT;
		% Desorption rate in [1 / (s * [T]^3)], cubic term of external dependence
		kD_TTT;
		EXT_LIN_KD_TTT;
	end

	methods

		function obj = ExtFunLinearBinding(kA, kD)
			%EXTFUNLINEARBINDING Constructs an ExtFunLinearBinding object with external function support
			%   OBJ = EXTFUNLINEARBINDING(KA) creates an ExtFunLinearBinding model with
			%   the given adsorption rates KA.
			%
			%   OBJ = EXTFUNLINEARBINDING(..., KD) also sets the desorption rates to KD.

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

			validateattributes(obj.kA, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kA');
			validateattributes(obj.kD, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kD');

			validateattributes(obj.kA_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kA_T');
			validateattributes(obj.kD_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kD_T');

			validateattributes(obj.kA_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kA_TT');
			validateattributes(obj.kD_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kD_TT');

			validateattributes(obj.kA_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kA_TTT');
			validateattributes(obj.kD_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kD_TTT');

			if ~isempty(obj.externalSource)
				validateattributes(obj.externalSource, {'double', 'int32'}, {'vector', 'nonempty', 'finite', 'real', '>=', 0}, '', 'externalSource');
				if ~any(numel(obj.externalSource) == [1, 2])
					error('CADET:invalidConfig', 'Expected externalSource to have either 1 or 2 entries but got %d entries.', numel(obj.externalSource));
				end
			end

			res = obj.validate@KineticQuasiStationaryBindingModel(nComponents, nBoundStates);
		end

		function val = get.kA(obj)
			val = obj.data.EXT_LIN_KA;
		end

		function set.kA(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kA');
			obj.data.EXT_LIN_KA = val;
			obj.hasChanged = true;
		end

		function val = get.kD(obj)
			val = obj.data.EXT_LIN_KD;
		end

		function set.kD(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kD');
			obj.data.EXT_LIN_KD = val;
			obj.hasChanged = true;
		end

		function val = get.kA_T(obj)
			val = obj.data.EXT_LIN_KA_T;
		end

		function set.kA_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kA_T');
			obj.data.EXT_LIN_KA_T = val;
			obj.hasChanged = true;
		end

		function val = get.kD_T(obj)
			val = obj.data.EXT_LIN_KD_T;
		end

		function set.kD_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kD_T');
			obj.data.EXT_LIN_KD_T = val;
			obj.hasChanged = true;
		end

		function val = get.kA_TT(obj)
			val = obj.data.EXT_LIN_KA_TT;
		end

		function set.kA_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kA_TT');
			obj.data.EXT_LIN_KA_TT = val;
			obj.hasChanged = true;
		end

		function val = get.kD_TT(obj)
			val = obj.data.EXT_LIN_KD_TT;
		end

		function set.kD_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kD_TT');
			obj.data.EXT_LIN_KD_TT = val;
			obj.hasChanged = true;
		end

		function val = get.kA_TTT(obj)
			val = obj.data.EXT_LIN_KA_TTT;
		end

		function set.kA_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kA_TTT');
			obj.data.EXT_LIN_KA_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.kD_TTT(obj)
			val = obj.data.EXT_LIN_KD_TTT;
		end

		function set.kD_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kD_TTT');
			obj.data.EXT_LIN_KD_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.EXT_LIN_KA(obj)
			val = obj.kA;
		end
		function set.EXT_LIN_KA(obj, val)
			obj.kA = val;
		end
		function val = get.EXT_LIN_KD(obj)
			val = obj.kD;
		end
		function set.EXT_LIN_KD(obj, val)
			obj.kD = val;
		end

		function val = get.EXT_LIN_KA_T(obj)
			val = obj.kA_T;
		end
		function set.EXT_LIN_KA_T(obj, val)
			obj.kA_T = val;
		end
		function val = get.EXT_LIN_KD_T(obj)
			val = obj.kD_T;
		end
		function set.EXT_LIN_KD_T(obj, val)
			obj.kD_T = val;
		end

		function val = get.EXT_LIN_KA_TT(obj)
			val = obj.kA_TT;
		end
		function set.EXT_LIN_KA_TT(obj, val)
			obj.kA_TT = val;
		end
		function val = get.EXT_LIN_KD_TT(obj)
			val = obj.kD_TT;
		end
		function set.EXT_LIN_KD_TT(obj, val)
			obj.kD_TT = val;
		end

		function val = get.EXT_LIN_KA_TTT(obj)
			val = obj.kA_TTT;
		end
		function set.EXT_LIN_KA_TTT(obj, val)
			obj.kA_TTT = val;
		end
		function val = get.EXT_LIN_KD_TTT(obj)
			val = obj.kD_TTT;
		end
		function set.EXT_LIN_KD_TTT(obj, val)
			obj.kD_TTT = val;
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
				if ~any(numel(val) == [1, 2])
					error('CADET:invalidConfig', 'Expected externalSource to have either 1 or 2 entries but got %d entries.', numel(val));
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
				obj = ExtFunLinearBinding();
				obj.loadobjInternal(S);
			end
		end
		
	end

end

% =============================================================================
%  CADET - The Chromatography Analysis and Design Toolkit
%  
%  Copyright (C) 2008-2018: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. obj program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies obj distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
