
classdef BiStericMassActionBinding < KineticQuasiStationaryBindingModel
	%BiStericMassActionBinding Bi steric mass action binding model with multiple binding sites
	%
	% See also BINDINGMODEL, KINETICQUASISTATIONARYBINDINGMODEL

	% Copyright: (C) 2008-2020 The CADET Authors
	%            See the license note at the end of the file.

	properties(Constant)
		name = 'BI_STERIC_MASS_ACTION'; % Name of the binding model according to CADET file format specs
	end

	properties (Dependent, Transient)
		% Adsorption rate in [m^3_MP / (m^3_SP * s)], rows correspond to binding sites
		kA;
		BISMA_KA;
		% Desorption rate in [1 / s], rows correspond to binding sites
		kD;
		BISMA_KD;
		% Characteristic charge [-], rows correspond to binding sites
		nu;
		BISMA_NU;
		% Steric factor [-], rows correspond to binding sites
		sigma;
		BISMA_SIGMA;
		% Ionic capacity in [mol / m^3_SP] for each binding site
		lambda;
		BISMA_LAMBDA;
		% Reference solid phase concentrations [mol / m^3_SP]
		refSolid;
		BISMA_REFQ;
		% Reference liquid phase concentrations [mol / m^3_MP]
		refLiquid;
		BISMA_REFC0;
	end

	methods

		function obj = BiStericMassActionBinding(kA, kD)
			%BISTERICMASSACTIONBINDING Constructs a BiStericMassActionBinding object
			%   OBJ = BISTERICMASSACTIONBINDING(KA) creates a BiStericMassActionBinding model with
			%   the given adsorption rates KA.
			%
			%   OBJ = BISTERICMASSACTIONBINDING(..., KD) also sets the desorption rates to KD.

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

			if any((nBoundStates ~= 0) & (nBoundStates ~= max(nBoundStates)))
				error('CADET:invalidConfig', 'Expected all components to have the same number of bound states (or 0).');
			end
			numStates = max(nBoundStates);

			% kA, kD, ... are stored as transpose, so each row corresponds to a binding site
			validateattributes(obj.data.BISMA_KA.', {'double'}, {'2d', 'nonempty', '>=', 0.0, 'finite', 'real', 'size', [numStates, nComponents]}, '', 'kA');
			validateattributes(obj.data.BISMA_KD.', {'double'}, {'2d', 'nonempty', '>=', 0.0, 'finite', 'real', 'size', [numStates, nComponents]}, '', 'kD');
			validateattributes(obj.data.BISMA_NU.', {'double'}, {'2d', 'nonempty', '>=', 0.0, 'finite', 'real', 'size', [numStates, nComponents]}, '', 'nu');
			validateattributes(obj.data.BISMA_SIGMA.', {'double'}, {'2d', 'nonempty', '>=', 0.0, 'finite', 'real', 'size', [numStates, nComponents]}, '', 'sigma');
			validateattributes(obj.data.BISMA_LAMBDA, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', numStates}, '', 'lambda');

			if ~isempty(obj.refSolid)
				validateattributes(obj.refSolid, {'double'}, {'vector', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refSolid');
				if ~any(numel(obj.refSolid) == [1, numStates])
					error('CADET:invalidConfig', 'Expected refSolid to have either 1 or %d entries but got %d entries.', numStates, numel(obj.refSolid));
				end
			end
			if ~isempty(obj.refLiquid)
				validateattributes(obj.refLiquid, {'double'}, {'vector', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refLiquid');
				if ~any(numel(obj.refLiquid) == [1, numStates])
					error('CADET:invalidConfig', 'Expected refSolid to have either 1 or %d entries but got %d entries.', numStates, numel(obj.refLiquid));
				end
			end

			res = obj.validate@KineticQuasiStationaryBindingModel(nComponents, nBoundStates);
		end

		function val = get.kA(obj)
			val = obj.data.BISMA_KA.';
		end

		function set.kA(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'kA');
			obj.data.BISMA_KA = val.';
			obj.hasChanged = true;
		end

		function val = get.kD(obj)
			val = obj.data.BISMA_KD.';
		end

		function set.kD(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'kD');
			obj.data.BISMA_KD = val.';
			obj.hasChanged = true;
		end

		function val = get.nu(obj)
			val = obj.data.BISMA_NU.';
		end

		function set.nu(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'nu');
			obj.data.BISMA_NU = val.';
			obj.hasChanged = true;
		end

		function val = get.sigma(obj)
			val = obj.data.BISMA_SIGMA.';
		end

		function set.sigma(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'sigma');
			obj.data.BISMA_SIGMA = val.';
			obj.hasChanged = true;
		end

		function val = get.lambda(obj)
			val = obj.data.BISMA_LAMBDA.';
		end

		function set.lambda(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'lambda');
			obj.data.BISMA_LAMBDA = val.';
			obj.hasChanged = true;
		end

		function val = get.refSolid(obj)
			if isfield(obj.data, 'BISMA_REFQ')
				val = obj.data.BISMA_REFQ;
			else
				val = [];
			end
		end

		function set.refSolid(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'BISMA_REFQ');
			else
				validateattributes(val, {'double'}, {'vector', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refSolid');
				obj.data.BISMA_REFQ = val;
			end
			obj.hasChanged = true;
		end

		function val = get.refLiquid(obj)
			if isfield(obj.data, 'BISMA_REFC0')
				val = obj.data.BISMA_REFC0;
			else
				val = [];
			end
		end

		function set.refLiquid(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'BISMA_REFC0');
			else
				validateattributes(val, {'double'}, {'vector', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refLiquid');
				obj.data.BISMA_REFC0 = val;
			end
			obj.hasChanged = true;
		end


		function val = get.BISMA_KA(obj)
			val = obj.kA;
		end
		function set.BISMA_KA(obj, val)
			obj.kA = val;
		end
		function val = get.BISMA_KD(obj)
			val = obj.kD;
		end
		function set.BISMA_KD(obj, val)
			obj.kD = val;
		end
		function val = get.BISMA_NU(obj)
			val = obj.nu;
		end
		function set.BISMA_NU(obj, val)
			obj.nu = val;
		end
		function val = get.BISMA_SIGMA(obj)
			val = obj.sigma;
		end
		function set.BISMA_SIGMA(obj, val)
			obj.sigma = val;
		end
		function val = get.BISMA_LAMBDA(obj)
			val = obj.lambda;
		end
		function set.BISMA_LAMBDA(obj, val)
			obj.lambda = val;
		end
		function val = get.BISMA_REFQ(obj)
			val = obj.refSolid;
		end
		function set.BISMA_REFQ(obj, val)
			obj.refSolid = val;
		end
		function val = get.BISMA_REFC0(obj)
			val = obj.refLiquid;
		end
		function set.BISMA_REFC0(obj, val)
			obj.refLiquid = val;
		end

		function res = assembleConfig(obj)
			%ASSEMBLECONFIG Assembles the configuration according to the CADET file format spec
			%   RES = ASSEMBLECONFIG() returns a Matlab struct RES that represents the
			%   binding model as detailed in the CADET file format spec.
			%
			% See also BINDINGMODEL.ASSEMBLECONFIG

			res = obj.data;

			% Linearize
			res.BISMA_KA = obj.data.BISMA_KA(:);
			res.BISMA_KD = obj.data.BISMA_KD(:);
			res.BISMA_NU = obj.data.BISMA_NU(:);
			res.BISMA_SIGMA = obj.data.BISMA_SIGMA(:);
		end
	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = BiStericMassActionBinding();
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
