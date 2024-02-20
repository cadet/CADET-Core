
classdef SaskaBinding < KineticQuasiStationaryBindingModel
	%SaskaBinding Saska binding model
	%
	% See also BINDINGMODEL, KINETICQUASISTATIONARYBINDINGMODEL

	% Copyright: (C) 2008-2024 The CADET Authors
	%            See the license note at the end of the file.

	properties(Constant)
		name = 'SASKA'; % Name of the binding model according to CADET file format specs
	end

	properties (Dependent, Transient)
		% Henri coefficient in [m^3_MP / (m^3_SP * s)]
		h;
		SASKA_H;
		% Quadratic factors in [m^6_MP / (m^3_SP * mol * s)]
		k;
		SASKA_K;
	end

	methods

		function obj = SaskaBinding(h, k)
			%SASKABINDING Constructs a SaskaBinding object
			%   OBJ = SASKABINDING(H) creates a SaskaBinding model with the
			%   given Henri coefficients H.
			%
			%   OBJ = SASKABINDING(..., K) also sets the quadratic factors to K.

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
			res = obj.validate@KineticQuasiStationaryBindingModel(nComponents, nBoundStates);
		end

		function val = get.h(obj)
			val = obj.data.SASKA_H;
		end

		function set.h(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'h');
			obj.data.SASKA_H = val;
			obj.hasChanged = true;
		end

		function val = get.k(obj)
			val = obj.data.SASKA_K;
			val = reshape(val, sqrt(numel(val)), sqrt(numel(val))).';
		end

		function set.k(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'k');
			nVals = sqrt(numel(val));
			if nVals ~= ceil(nVals)
				error('CADET:invalidConfig', 'Expected k to have a square number of values.');
			end
			val = val.';
			obj.data.SASKA_K = val(:);
			obj.hasChanged = true;
		end

		function val = get.SASKA_H(obj)
			val = obj.h;
		end
		function set.SASKA_H(obj, val)
			obj.h = val;
		end
		function val = get.SASKA_K(obj)
			val = obj.k;
		end
		function set.SASKA_K(obj, val)
			obj.k = val;
		end
	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = SaskaBinding();
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
