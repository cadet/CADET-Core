
classdef ExtFunSaskaBinding < KineticQuasiStationaryBindingModel
	%ExtFunSaskaBinding Saska binding model with external function support
	%
	% See also BINDINGMODEL, SASKABINDING, KINETICQUASISTATIONARYBINDINGMODEL

	% Copyright: (C) 2008-2019 The CADET Authors
	%            See the license note at the end of the file.

	properties(Constant)
		name = 'EXT_SASKA'; % Name of the binding model according to CADET file format specs
	end

	properties (Dependent, Transient)
		% Henri coefficient in [m^3_MP / (m^3_SP * s)], constant term of external dependence
		h;
		EXT_SASKA_H;
		% Henri coefficient in [m^3_MP / (m^3_SP * s * [T])], linear term of external dependence
		h_T;
		EXT_SASKA_H_T;
		% Henri coefficient in [m^3_MP / (m^3_SP * s * [T]^2)], quadratic term of external dependence
		h_TT;
		EXT_SASKA_H_TT;
		% Henri coefficient in [m^3_MP / (m^3_SP * s * [T]^3)], cubic term of external dependence
		h_TTT;
		EXT_SASKA_H_TTT;
		% Quadratic factors in [m^6_MP / (m^3_SP * mol * s)], constant term of external dependence
		k;
		EXT_SASKA_K;
		% Quadratic factors in [m^6_MP / (m^3_SP * mol * s * [T])], linear term of external dependence
		k_T;
		EXT_SASKA_K_T;
		% Quadratic factors in [m^6_MP / (m^3_SP * mol * s * [T]^2)], quadratic term of external dependence
		k_TT;
		EXT_SASKA_K_TT;
		% Quadratic factors in [m^6_MP / (m^3_SP * mol * s * [T]^3)], cubic term of external dependence
		k_TTT;
		EXT_SASKA_K_TTT;
		% Indices of external functions (0-based)
		externalSource;
		EXTFUN;
	end

	methods

		function obj = ExtFunSaskaBinding(h, k)
			%EXTFUNSASKABINDING Constructs an ExtFunSaskaBinding object with external function support
			%   OBJ = EXTFUNSASKABINDING(H) creates an ExtFunSaskaBinding model with the
			%   given Henri coefficients H.
			%
			%   OBJ = EXTFUNSASKABINDING(..., K) also sets the quadratic factors to K.

			obj = obj@KineticQuasiStationaryBindingModel();

			if nargin >= 1
				obj.h = h;
			end
			if nargin >= 2
				obj.k = k;
			end
		end

		function res = validate(obj, nComponents, nBoundStates)
			%VALIDATE Validates the binding model parameters
			%   RES = VALIDATE(NCOMPONENTS, NBOUNDSTATES) uses the number of
			%   components NCOMPONENTS of the model and the number of bound
			%   states of each of these components as given by the vector
			%   NBOUNDSTATES to validate the parameters of the binding model.
			%   Returns true in RES if everything is fine and false otherwise.

			validateattributes(obj.h, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'h');
			validateattributes(obj.k, {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [nComponents, nComponents]}, '', 'k');

			validateattributes(obj.h_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'h_T');
			validateattributes(obj.k_T, {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [nComponents, nComponents]}, '', 'k_T');

			validateattributes(obj.h_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'h_TT');
			validateattributes(obj.k_TT, {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [nComponents, nComponents]}, '', 'k_TT');

			validateattributes(obj.h_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'h_TTT');
			validateattributes(obj.k_TTT, {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [nComponents, nComponents]}, '', 'k_TTT');

			if ~isempty(obj.externalSource)
				validateattributes(obj.externalSource, {'double', 'int32'}, {'vector', 'nonempty', 'finite', 'real', '>=', 0}, '', 'externalSource');
				if ~any(numel(obj.externalSource) == [1, 2])
					error('CADET:invalidConfig', 'Expected externalSource to have either 1 or 2 entries but got %d entries.', numel(obj.externalSource));
				end
			end

			res = obj.validate@KineticQuasiStationaryBindingModel(nComponents, nBoundStates);
		end

		function val = get.h(obj)
			val = obj.data.EXT_SASKA_H;
		end

		function set.h(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'h');
			obj.data.EXT_SASKA_H = val;
			obj.hasChanged = true;
		end

		function val = get.k(obj)
			val = obj.data.EXT_SASKA_K;
			val = reshape(val, sqrt(numel(val)), sqrt(numel(val))).';
		end

		function set.k(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'k');
			nVals = sqrt(numel(val));
			if nVals ~= ceil(nVals)
				error('CADET:invalidConfig', 'Expected k to have a square number of values.');
			end
			val = val.';
			obj.data.EXT_SASKA_K = val(:);
			obj.hasChanged = true;
		end

		function val = get.k_T(obj)
			val = obj.data.EXT_SASKA_K_T;
			val = reshape(val, sqrt(numel(val)), sqrt(numel(val))).';
		end

		function set.k_T(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'k_T');
			nVals = sqrt(numel(val));
			if nVals ~= ceil(nVals)
				error('CADET:invalidConfig', 'Expected k to have a square number of values.');
			end
			val = val.';
			obj.data.EXT_SASKA_K_T = val(:);
			obj.hasChanged = true;
		end

		function val = get.h_T(obj)
			val = obj.data.EXT_SASKA_H_T;
		end

		function set.h_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'h_T');
			obj.data.EXT_SASKA_H_T = val;
			obj.hasChanged = true;
		end

		function val = get.k_TT(obj)
			val = obj.data.EXT_SASKA_K_TT;
			val = reshape(val, sqrt(numel(val)), sqrt(numel(val))).';
		end

		function set.k_TT(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'k_TT');
			nVals = sqrt(numel(val));
			if nVals ~= ceil(nVals)
				error('CADET:invalidConfig', 'Expected k to have a square number of values.');
			end
			val = val.';
			obj.data.EXT_SASKA_K_TT = val(:);
			obj.hasChanged = true;
		end

		function val = get.h_TTT(obj)
			val = obj.data.EXT_SASKA_H_TTT;
		end

		function set.h_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'h_TTT');
			obj.data.EXT_SASKA_H_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.k_TTT(obj)
			val = obj.data.EXT_SASKA_K_TTT;
			val = reshape(val, sqrt(numel(val)), sqrt(numel(val))).';
		end

		function set.k_TTT(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'k_TTT');
			nVals = sqrt(numel(val));
			if nVals ~= ceil(nVals)
				error('CADET:invalidConfig', 'Expected k to have a square number of values.');
			end
			val = val.';
			obj.data.EXT_SASKA_K_TTT = val(:);
			obj.hasChanged = true;
		end

		function val = get.EXT_SASKA_H(obj)
			val = obj.h;
		end
		function set.EXT_SASKA_H(obj, val)
			obj.h = val;
		end
		function val = get.EXT_SASKA_K(obj)
			val = obj.k;
		end
		function set.EXT_SASKA_K(obj, val)
			obj.k = val;
		end

		function val = get.EXT_SASKA_H_T(obj)
			val = obj.h_T;
		end
		function set.EXT_SASKA_H_T(obj, val)
			obj.h_T = val;
		end
		function val = get.EXT_SASKA_K_T(obj)
			val = obj.k_T;
		end
		function set.EXT_SASKA_K_T(obj, val)
			obj.k_T = val;
		end

		function val = get.EXT_SASKA_H_TT(obj)
			val = obj.h_TT;
		end
		function set.EXT_SASKA_H_TT(obj, val)
			obj.h_TT = val;
		end
		function val = get.EXT_SASKA_K_TT(obj)
			val = obj.k_TT;
		end
		function set.EXT_SASKA_K_TT(obj, val)
			obj.k_TT = val;
		end

		function val = get.EXT_SASKA_H_TTT(obj)
			val = obj.h_TTT;
		end
		function set.EXT_SASKA_H_TTT(obj, val)
			obj.h_TTT = val;
		end
		function val = get.EXT_SASKA_K_TTT(obj)
			val = obj.k_TTT;
		end
		function set.EXT_SASKA_K_TTT(obj, val)
			obj.k_TTT = val;
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
				obj = ExtFunSaskaBinding();
				obj.loadobjInternal(S);
			end
		end
		
	end

end

% =============================================================================
%  CADET - The Chromatography Analysis and Design Toolkit
%  
%  Copyright (C) 2008-2019: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. obj program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies obj distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
