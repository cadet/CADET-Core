
classdef SelfAssociationBinding < KineticQuasiStationaryBindingModel
	%SelfAssociationBinding Self association binding model
	%
	% See also BINDINGMODEL, KINETICQUASISTATIONARYBINDINGMODEL

	% Copyright: (C) 2008-2024 The CADET Authors
	%            See the license note at the end of the file.

	properties(Constant)
		name = 'SELF_ASSOCIATION'; % Name of the binding model according to CADET file format specs
	end

	properties (Dependent, Transient)
		% Adsorption rate in [m^3_MP / (m^3_SP * s)]
		kA;
		SAI_KA1;
		% Adsorption rate of dimerization in [m^6_MP / (m^6_SP * s)]
		kAdimer;
		SAI_KA2;
		% Desorption rate in [1 / s]
		kD;
		SAI_KD;
		% Characteristic charge [-]
		nu;
		SAI_NU;
		% Steric factor [-]
		sigma;
		SAI_SIGMA;
		% Ionic capacity in [mol / m^3_SP]
		lambda;
		SAI_LAMBDA;
		% Reference solid phase concentration [mol / m^3_SP]
		refSolid;
		SAI_REFQ;
		% Reference liquid phase concentration [mol / m^3_MP]
		refLiquid;
		SAI_REFC0;
	end

	methods

		function obj = SelfAssociationBinding(kA, kD)
			%SELFASSOCIATIONBINDING Constructs a SelfAssociationBinding object
			%   OBJ = SELFASSOCIATIONBINDING(KA) creates a SelfAssociationBinding model with
			%   the given adsorption rates KA.
			%
			%   OBJ = SELFASSOCIATIONBINDING(..., KD) also sets the desorption rates to KD.

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

			validateattributes(obj.kA, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'kA');
			validateattributes(obj.kAdimer, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'kAdimer');
			validateattributes(obj.kD, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'kD');
			validateattributes(obj.nu, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'nu');
			validateattributes(obj.sigma, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'sigma');
			validateattributes(obj.lambda, {'double'}, {'scalar', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'lambda');

			if ~isempty(obj.refSolid)
				validateattributes(obj.refSolid, {'double'}, {'scalar', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refSolid');
			end
			if ~isempty(obj.refLiquid)
				validateattributes(obj.refLiquid, {'double'}, {'scalar', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refLiquid');
			end

			res = obj.validate@KineticQuasiStationaryBindingModel(nComponents, nBoundStates);
		end

		function val = get.kA(obj)
			val = obj.data.SAI_KA1;
		end

		function set.kA(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'kA');
			obj.data.SAI_KA1 = val;
			obj.hasChanged = true;
		end

		function val = get.kAdimer(obj)
			val = obj.data.SAI_KA2;
		end

		function set.kAdimer(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'kAdimer');
			obj.data.SAI_KA2 = val;
			obj.hasChanged = true;
		end

		function val = get.kD(obj)
			val = obj.data.SAI_KD;
		end

		function set.kD(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'kD');
			obj.data.SAI_KD = val;
			obj.hasChanged = true;
		end

		function val = get.nu(obj)
			val = obj.data.SAI_NU;
		end

		function set.nu(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'nu');
			obj.data.SAI_NU = val;
			obj.hasChanged = true;
		end

		function val = get.sigma(obj)
			val = obj.data.SAI_SIGMA;
		end

		function set.sigma(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'sigma');
			obj.data.SAI_SIGMA = val;
			obj.hasChanged = true;
		end

		function val = get.lambda(obj)
			val = obj.data.SAI_LAMBDA;
		end

		function set.lambda(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'lambda');
			obj.data.SAI_LAMBDA = val;
			obj.hasChanged = true;
		end

		function val = get.refSolid(obj)
			if isfield(obj.data, 'SAI_REFQ')
				val = obj.data.SAI_REFQ;
			else
				val = [];
			end
		end

		function set.refSolid(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'SAI_REFQ');
			else
				validateattributes(val, {'double'}, {'scalar', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refSolid');
				obj.data.SAI_REFQ = val;
			end
			obj.hasChanged = true;
		end

		function val = get.refLiquid(obj)
			if isfield(obj.data, 'SAI_REFC0')
				val = obj.data.SAI_REFC0;
			else
				val = [];
			end
		end

		function set.refLiquid(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'SAI_REFC0');
			else
				validateattributes(val, {'double'}, {'scalar', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refLiquid');
				obj.data.SAI_REFC0 = val;
			end
			obj.hasChanged = true;
		end

		
		function val = get.SAI_KA1(obj)
			val = obj.kA;
		end
		function set.SAI_KA1(obj, val)
			obj.kA = val;
		end
		function val = get.SAI_KA2(obj)
			val = obj.kAdimer;
		end
		function set.SAI_KA2(obj, val)
			obj.kAdimer = val;
		end
		function val = get.SAI_KD(obj)
			val = obj.kD;
		end
		function set.SAI_KD(obj, val)
			obj.kD = val;
		end
		function val = get.SAI_LAMBDA(obj)
			val = obj.lambda;
		end
		function set.SAI_LAMBDA(obj, val)
			obj.lambda = val;
		end
		function val = get.SAI_NU(obj)
			val = obj.nu;
		end
		function set.SAI_NU(obj, val)
			obj.nu = val;
		end
		function val = get.SAI_SIGMA(obj)
			val = obj.sigma;
		end
		function set.SAI_SIGMA(obj, val)
			obj.sigma = val;
		end
		function val = get.SAI_REFQ(obj)
			val = obj.refSolid;
		end
		function set.SAI_REFQ(obj, val)
			obj.refSolid = val;
		end
		function val = get.SAI_REFC0(obj)
			val = obj.refLiquid;
		end
		function set.SAI_REFC0(obj, val)
			obj.refLiquid = val;
		end

	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = SelfAssociationBinding();
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
