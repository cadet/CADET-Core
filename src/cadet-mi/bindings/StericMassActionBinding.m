
classdef StericMassActionBinding < KineticQuasiStationaryBindingModel
	%StericMassActionBinding Steric mass action binding model
	%
	% See also BINDINGMODEL, KINETICQUASISTATIONARYBINDINGMODEL

	% Copyright: (C) 2008-2017 The CADET Authors
	%            See the license note at the end of the file.

	properties (Constant, Access = 'protected')
		hasConsistencySolver = true; % Determines whether this binding model has a consistency solver
	end

	properties(Constant)
		name = 'STERIC_MASS_ACTION'; % Name of the binding model according to CADET file format specs
	end

	properties (Dependent, Transient)
		% Adsorption rate in [m^3_MP / (m^3_SP * s)]
		kA;
		SMA_KA;
		% Desorption rate in [1 / s]
		kD;
		SMA_KD;
		% Characteristic charge [-]
		nu;
		SMA_NU;
		% Steric factor [-]
		sigma;
		SMA_SIGMA;
		% Ionic capacity in [mol / m^3_SP]
		lambda;
		SMA_LAMBDA;
		% Reference solid phase concentration [mol / m^3_SP]
		refSolid;
		SMA_REFQ;
		% Reference liquid phase concentration [mol / m^3_MP]
		refLiquid;
		SMA_REFC0;
	end

	methods

		function obj = StericMassActionBinding(kA, kD)
			%STERICMASSACTIONBINDING Constructs a StericMassActionBinding object
			%   OBJ = STERICMASSACTIONBINDING(KA) creates a StericMassActionBinding model with
			%   the given adsorption rates KA.
			%
			%   OBJ = STERICMASSACTIONBINDING(..., KD) also sets the desorption rates to KD.

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
			val = obj.data.SMA_KA;
		end

		function set.kA(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'kA');
			obj.data.SMA_KA = val;
			obj.hasChanged = true;
		end

		function val = get.kD(obj)
			val = obj.data.SMA_KD;
		end

		function set.kD(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'kD');
			obj.data.SMA_KD = val;
			obj.hasChanged = true;
		end

		function val = get.nu(obj)
			val = obj.data.SMA_NU;
		end

		function set.nu(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'nu');
			obj.data.SMA_NU = val;
			obj.hasChanged = true;
		end

		function val = get.sigma(obj)
			val = obj.data.SMA_SIGMA;
		end

		function set.sigma(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'sigma');
			obj.data.SMA_SIGMA = val;
			obj.hasChanged = true;
		end

		function val = get.lambda(obj)
			val = obj.data.SMA_LAMBDA;
		end

		function set.lambda(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'lambda');
			obj.data.SMA_LAMBDA = val;
			obj.hasChanged = true;
		end

		function val = get.refSolid(obj)
			if isfield(obj.data, 'SMA_REFQ')
				val = obj.data.SMA_REFQ;
			else
				val = [];
			end
		end

		function set.refSolid(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'SMA_REFQ');
			else
				validateattributes(val, {'double'}, {'scalar', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refSolid');
				obj.data.SMA_REFQ = val;
			end
			obj.hasChanged = true;
		end

		function val = get.refLiquid(obj)
			if isfield(obj.data, 'SMA_REFC0')
				val = obj.data.SMA_REFC0;
			else
				val = [];
			end
		end

		function set.refLiquid(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'SMA_REFC0');
			else
				validateattributes(val, {'double'}, {'scalar', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refLiquid');
				obj.data.SMA_REFC0 = val;
			end
			obj.hasChanged = true;
		end

		function val = get.SMA_KA(obj)
			val = obj.kA;
		end
		function set.SMA_KA(obj, val)
			obj.kA = val;
		end
		function val = get.SMA_KD(obj)
			val = obj.kD;
		end
		function set.SMA_KD(obj, val)
			obj.kD = val;
		end
		function val = get.SMA_LAMBDA(obj)
			val = obj.lambda;
		end
		function set.SMA_LAMBDA(obj, val)
			obj.lambda = val;
		end
		function val = get.SMA_NU(obj)
			val = obj.nu;
		end
		function set.SMA_NU(obj, val)
			obj.nu = val;
		end
		function val = get.SMA_SIGMA(obj)
			val = obj.sigma;
		end
		function set.SMA_SIGMA(obj, val)
			obj.sigma = val;
		end
		function val = get.SMA_REFQ(obj)
			val = obj.refSolid;
		end
		function set.SMA_REFQ(obj, val)
			obj.refSolid = val;
		end
		function val = get.SMA_REFC0(obj)
			val = obj.refLiquid;
		end
		function set.SMA_REFC0(obj, val)
			obj.refLiquid = val;
		end

	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = StericMassActionBinding();
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
