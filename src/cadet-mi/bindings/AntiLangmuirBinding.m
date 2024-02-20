
classdef AntiLangmuirBinding < KineticQuasiStationaryBindingModel
	%AntiLangmuirBinding Anti-Langmuir binding model
	%
	% See also BINDINGMODEL

	% Copyright: (C) 2008-2024 The CADET Authors
	%            See the license note at the end of the file.

	properties(Constant)
		name = 'MULTI_COMPONENT_ANTILANGMUIR'; % Name of the binding model according to CADET file format specs
	end

	properties (Dependent, Transient)
		% Adsorption rate in [m^3_MP / (mol * s)]
		kA;
		MCAL_KA;
		% Desorption rate in [1 / s]
		kD;
		MCAL_KD;
		% Capacity in [mol / m^3_SP]
		qMax;
		MCAL_QMAX;
		% Anti-Langmuir coefficient [-]
		antiFactor;
		MCAL_ANTILANGMUIR;
	end

	methods

		function obj = AntiLangmuirBinding(kA, kD, qMax, antiFactor)
			%ANTILANGMUIRBINDING Constructs an Anti-Langmuir binding model object
			%   OBJ = ANTILANGMUIRBINDING(KA) creates an Anti-Langmuir binding model with
			%   the given adsorption rates KA.
			%
			%   OBJ = ANTILANGMUIRBINDING(..., KD) also sets the desorption rates to KD.
			%
			%   OBJ = ANTILANGMUIRBINDING(..., KD, QMAX) also sets the desorption rates to KD
			%   and the capacity to QMAX.
			%
			%   OBJ = ANTILANGMUIRBINDING(..., KD, QMAX, ANTIFACTOR) also sets the desorption
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

			validateattributes(obj.kA, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'kA');
			validateattributes(obj.kD, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'kD');
			validateattributes(obj.qMax, {'double'}, {'vector', 'nonempty', '>', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'qMax');
			validateattributes(obj.antiFactor, {'double'}, {'vector', 'nonempty', '>=', -1.0, '<=', 1.0, 'finite', 'real', 'numel', nComponents}, '', 'antiFactor');
			if ~all((obj.antiFactor == -1) | (obj.antiFactor == 1))
				error('CADET:invalidConfig', 'Expected antiFactor to be an array consisting of -1 and 1.')
			end
			res = obj.validate@KineticQuasiStationaryBindingModel(nComponents, nBoundStates);
		end

		function val = get.kA(obj)
			val = obj.data.MCAL_KA;
		end

		function set.kA(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'kA');
			obj.data.MCAL_KA = val;
			obj.hasChanged = true;
		end

		function val = get.kD(obj)
			val = obj.data.MCAL_KD;
		end

		function set.kD(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'kD');
			obj.data.MCAL_KD = val;
			obj.hasChanged = true;
		end

		function val = get.qMax(obj)
			val = obj.data.MCAL_QMAX;
		end

		function set.qMax(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'qMax');
			obj.data.MCAL_QMAX = val;
			obj.hasChanged = true;
		end

		function val = get.antiFactor(obj)
			val = obj.data.MCAL_ANTILANGMUIR;
		end

		function set.antiFactor(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'antiFactor');
			if ~all((val == -1) | (val == 1))
				error('CADET:invalidConfig', 'Expected antiFactor to be an array consisting of -1 and 1.')
			end
			obj.data.MCAL_ANTILANGMUIR = val;
			obj.hasChanged = true;
		end

		function val = get.MCAL_KA(obj)
			val = obj.kA;
		end
		function set.MCAL_KA(obj, val)
			obj.kA = val;
		end
		function val = get.MCAL_KD(obj)
			val = obj.kD;
		end
		function set.MCAL_KD(obj, val)
			obj.kD = val;
		end
		function val = get.MCAL_QMAX(obj)
			val = obj.qMax;
		end
		function set.MCAL_QMAX(obj, val)
			obj.qMax = val;
		end
		function val = get.MCAL_ANTILANGMUIR(obj)
			val = obj.antiFactor;
		end
		function set.MCAL_ANTILANGMUIR(obj, val)
			obj.antiFactor = val;
		end

	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = AntiLangmuirBinding();
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
