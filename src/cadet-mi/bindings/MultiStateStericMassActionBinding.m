
classdef MultiStateStericMassActionBinding < KineticQuasiStationaryBindingModel
	%MultiStateStericMassActionBinding Multi-state steric mass action binding model
	%
	% See also BINDINGMODEL, KINETICQUASISTATIONARYBINDINGMODEL

	% Copyright: (C) 2008-2016 The CADET Authors
	%            See the license note at the end of the file.

	properties (Constant, Access = 'protected')
		hasConsistencySolver = true; % Determines whether this binding model has a consistency solver
	end

	properties(Constant)
		name = 'MULTISTATE_STERIC_MASS_ACTION'; % Name of the binding model according to CADET file format specs
	end

	properties (Dependent, Transient)
		% Adsorption rate in [m^3_MP / (m^3_SP * s)] as linearized component-major array
		kA;
		MSSMA_KA;
		% Desorption rate in [1 / s] as linearized component-major array
		kD;
		MSSMA_KD;
		% Characteristic charge [-] as linearized component-major array
		nu;
		MSSMA_NU;
		% Steric factor [-] as linearized component-major array
		sigma;
		MSSMA_SIGMA;
		% Conversion rates between different bound states [1 / s] as linearized component-row-major array
		rates;
		MSSMA_RATES;
		% Ionic capacity in [mol / m^3_SP] as linearized component-major array
		lambda;
		MSSMA_LAMBDA;
		% Reference solid phase concentration [mol / m^3_SP]
		refSolid;
		MSSMA_REFQ;
		% Reference liquid phase concentration [mol / m^3_MP]
		refLiquid;
		MSSMA_REFC0;
	end

	methods

		function obj = MultiStateStericMassActionBinding(kA, kD)
			%MULTISTATESTERICMASSACTIONBINDING Constructs a MultiStateStericMassActionBinding object
			%   OBJ = MULTISTATESTERICMASSACTIONBINDING(KA) creates a MultiStateStericMassActionBinding
			%   model with the given adsorption rates KA.
			%
			%   OBJ = MULTISTATESTERICMASSACTIONBINDING(..., KD) also sets the desorption rates to KD.

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
			validateattributes(obj.nu, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', sum(nBoundStates)}, '', 'nu');
			validateattributes(obj.sigma, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', sum(nBoundStates)}, '', 'sigma');
			validateattributes(obj.rates, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', sum(nBoundStates.^2)}, '', 'rates');
			validateattributes(obj.lambda, {'double'}, {'scalar', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'lambda');
			
			% Check if nu is nondecreasing for each component
			idx = 1;
			for i = 1:nComponents
				validateattributes(obj.nu(idx:idx+nBoundStates(i)-1), {'double'}, {'nondecreasing'}, '', sprintf('nu component %d', i-1));
				idx = idx + nBoundStates(i);
			end
			
			if ~isempty(obj.refSolid)
				validateattributes(obj.refSolid, {'double'}, {'scalar', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refSolid');
			end
			if ~isempty(obj.refLiquid)
				validateattributes(obj.refLiquid, {'double'}, {'scalar', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refLiquid');
			end

			res = obj.validate@KineticQuasiStationaryBindingModel(nComponents, nBoundStates);
		end

		function val = get.kA(obj)
			val = obj.data.MSSMA_KA;
		end

		function set.kA(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'kA');
			obj.data.MSSMA_KA = val;
			obj.hasChanged = true;
		end

		function val = get.kD(obj)
			val = obj.data.MSSMA_KD;
		end

		function set.kD(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'kD');
			obj.data.MSSMA_KD = val;
			obj.hasChanged = true;
		end

		function val = get.nu(obj)
			val = obj.data.MSSMA_NU;
		end

		function set.nu(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'nu');
			obj.data.MSSMA_NU = val;
			obj.hasChanged = true;
		end

		function val = get.sigma(obj)
			val = obj.data.MSSMA_SIGMA;
		end

		function set.sigma(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'sigma');
			obj.data.MSSMA_SIGMA = val;
			obj.hasChanged = true;
		end

		function val = get.lambda(obj)
			val = obj.data.MSSMA_LAMBDA;
		end

		function set.lambda(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'lambda');
			obj.data.MSSMA_LAMBDA = val;
			obj.hasChanged = true;
		end

		function val = get.rates(obj)
			val = obj.data.MSSMA_RATES;
		end

		function set.rates(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'rates');
			obj.data.MSSMA_RATES = val;
			obj.hasChanged = true;
		end

		function val = get.refSolid(obj)
			if isfield(obj.data, 'MSSMA_REFQ')
				val = obj.data.MSSMA_REFQ;
			else
				val = [];
			end
		end

		function set.refSolid(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'MSSMA_REFQ');
			else
				validateattributes(val, {'double'}, {'scalar', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refSolid');
				obj.data.MSSMA_REFQ = val;
			end
			obj.hasChanged = true;
		end

		function val = get.refLiquid(obj)
			if isfield(obj.data, 'MSSMA_REFC0')
				val = obj.data.MSSMA_REFC0;
			else
				val = [];
			end
		end

		function set.refLiquid(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'MSSMA_REFC0');
			else
				validateattributes(val, {'double'}, {'scalar', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refLiquid');
				obj.data.MSSMA_REFC0 = val;
			end
			obj.hasChanged = true;
		end


		function val = get.MSSMA_KA(obj)
			val = obj.kA;
		end
		function set.MSSMA_KA(obj, val)
			obj.kA = val;
		end
		function val = get.MSSMA_KD(obj)
			val = obj.kD;
		end
		function set.MSSMA_KD(obj, val)
			obj.kD = val;
		end
		function val = get.MSSMA_LAMBDA(obj)
			val = obj.lambda;
		end
		function set.MSSMA_LAMBDA(obj, val)
			obj.lambda = val;
		end
		function val = get.MSSMA_NU(obj)
			val = obj.nu;
		end
		function set.MSSMA_NU(obj, val)
			obj.nu = val;
		end
		function val = get.MSSMA_SIGMA(obj)
			val = obj.sigma;
		end
		function set.MSSMA_SIGMA(obj, val)
			obj.sigma = val;
		end
		function val = get.MSSMA_RATES(obj)
			val = obj.rates;
		end
		function set.MSSMA_RATES(obj, val)
			obj.rates = val;
		end
		function val = get.MSSMA_REFQ(obj)
			val = obj.refSolid;
		end
		function set.MSSMA_REFQ(obj, val)
			obj.refSolid = val;
		end
		function val = get.MSSMA_REFC0(obj)
			val = obj.refLiquid;
		end
		function set.MSSMA_REFC0(obj, val)
			obj.refLiquid = val;
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
			%   This implementation assumes component-state-major ordering, which differs from the base
			%   class that uses section-boundphase-component-major ordering.
			%
			% See also BINDINGMODEL.OFFSETTOPARAMETER

			offset = 0;
			if (param.SENS_BOUNDPHASE ~= -1)
				offset = offset + param.SENS_BOUNDPHASE;
			end
			if (param.SENS_COMP ~= -1)
				offset = offset + param.SENS_COMP * sum(nBoundStates(1:param.SENS_COMP));
			end
		end

	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = MultiStateStericMassActionBinding();
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
