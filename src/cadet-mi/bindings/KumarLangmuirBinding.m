
classdef KumarLangmuirBinding < KineticQuasiStationaryBindingModel
	%KumarLangmuirBinding Kumar-Langmuir binding model
	%
	% See also BINDINGMODEL, KINETICQUASISTATIONARYBINDINGMODEL
	
	% Copyright: (C) 2008-2017 The CADET Authors
	%            See the license note at the end of the file.

	properties (Constant, Access = 'protected')
		hasConsistencySolver = true; % Determines whether this binding model has a consistency solver
	end

	properties(Constant)
		name = 'KUMAR_MULTI_COMPONENT_LANGMUIR'; % Name of the binding model according to CADET file format specs
	end

	properties (Dependent, Transient)
		% Adsorption rate in [m^3_MP / (mol * s)]
		kA;
		KMCL_KA;
		% Desorption rate in [m^(3 * nu_i)_MP / (mol^(nu_i) * s)]
		kD;
		KMCL_KD;
		% Activation temperature in [K]
		kAct;
		KMCL_KACT;
		% Capacity in [mol / m^3_SP]
		qMax;
		KMCL_QMAX;
		% Characteristic charge [-]
		nu;
		KMCL_NU;
		% Temperature in [K]
		temperature;
		KMCL_TEMP;
	end

	methods

		function obj = KumarLangmuirBinding(kA, kD, qMax)
			%KUMARLANGMUIRBINDING Constructs a KumarLangmuirBinding object
			%   OBJ = KUMARLANGMUIRBINDING(KA) creates a KumarLangmuirBinding model with
			%   the given adsorption rates KA.
			%
			%   OBJ = KUMARLANGMUIRBINDING(..., KD) also sets the desorption rates to KD.
			%
			%   OBJ = KUMARLANGMUIRBINDING(..., KD, QMAX) also sets the desorption rates
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
			validateattributes(obj.kAct, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'kAct');
			validateattributes(obj.qMax, {'double'}, {'vector', 'nonempty', '>', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'qMax');
			validateattributes(obj.nu, {'double'}, {'vector', 'nonempty', '>', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'nu');
			validateattributes(obj.temperature, {'double'}, {'scalar', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'temperature');
			res = obj.validate@KineticQuasiStationaryBindingModel(nComponents, nBoundStates);
		end

		function val = get.kA(obj)
			val = obj.data.KMCL_KA;
		end

		function set.kA(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'kA');
			obj.data.KMCL_KA = val;
			obj.hasChanged = true;
		end

		function val = get.kD(obj)
			val = obj.data.KMCL_KD;
		end

		function set.kD(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'kD');
			obj.data.KMCL_KD = val;
			obj.hasChanged = true;
		end

		function val = get.kAct(obj)
			val = obj.data.KMCL_KACT;
		end

		function set.kAct(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'kAct');
			obj.data.KMCL_KACT = val;
			obj.hasChanged = true;
		end

		function val = get.qMax(obj)
			val = obj.data.KMCL_QMAX;
		end

		function set.qMax(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'qMax');
			obj.data.KMCL_QMAX = val;
			obj.hasChanged = true;
		end

		function val = get.nu(obj)
			val = obj.data.KMCL_NU;
		end

		function set.nu(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'nu');
			obj.data.KMCL_NU = val;
			obj.hasChanged = true;
		end

		function val = get.temperature(obj)
			val = obj.data.KMCL_TEMP;
		end

		function set.temperature(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'temperature');
			obj.data.KMCL_TEMP = val;
			obj.hasChanged = true;
		end

		function val = get.KMCL_KA(obj)
			val = obj.kA;
		end
		function set.KMCL_KA(obj, val)
			obj.kA = KMCL_KA;
		end
		function val = get.KMCL_KD(obj)
			val = obj.kD;
		end
		function set.KMCL_KD(obj, val)
			obj.kD = KMCL_KD;
		end
		function val = get.KMCL_KACT(obj)
			val = obj.kAct;
		end
		function set.KMCL_KACT(obj, val)
			obj.kAct = KMCL_KACT;
		end
		function val = get.KMCL_QMAX(obj)
			val = obj.qMax;
		end
		function set.KMCL_QMAX(obj, val)
			obj.qMax = KMCL_QMAX;
		end
		function val = get.KMCL_NU(obj)
			val = obj.nu;
		end
		function set.KMCL_NU(obj, val)
			obj.nu = KMCL_NU;
		end
		function val = get.KMCL_TEMP(obj)
			val = obj.temperature;
		end
		function set.KMCL_TEMP(obj, val)
			obj.temperature = KMCL_TEMP;
		end

	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = KumarLangmuirBinding();
				obj.loadobjInternal(S);
			end
		end
		
	end

end

% =============================================================================
%  CADET - The Chromatography Analysis and Design Toolkit
%  
%  Copyright (C) 2008-2017: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. obj program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies obj distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
