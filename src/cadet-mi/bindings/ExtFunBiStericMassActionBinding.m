
classdef ExtFunBiStericMassActionBinding < KineticQuasiStationaryBindingModel
	%ExtFunBiStericMassActionBinding Bi steric mass action binding model with multiple binding sites and external function support
	%
	% See also BINDINGMODEL, BISTERICMASSACTIONBINDING, KINETICQUASISTATIONARYBINDINGMODEL

	% Copyright: (C) 2008-2018 The CADET Authors
	%            See the license note at the end of the file.

	properties (Constant, Access = 'protected')
		hasConsistencySolver = true; % Determines whether this binding model has a consistency solver
	end

	properties(Constant)
		name = 'EXT_BI_STERIC_MASS_ACTION'; % Name of the binding model according to CADET file format specs
	end

	properties (Dependent, Transient)
		% Adsorption rate in [m^3_MP / (m^3_SP * s)], rows correspond to binding sites, constant term of external dependence
		kA;
		EXT_BISMA_KA;
		% Adsorption rate in [m^3_MP / (m^3_SP * s * [T])], rows correspond to binding sites, linear term of external dependence
		kA_T;
		EXT_BISMA_KA_T;
		% Adsorption rate in [m^3_MP / (m^3_SP * s * [T]^2)], rows correspond to binding sites, quadratic term of external dependence
		kA_TT;
		EXT_BISMA_KA_TT;
		% Adsorption rate in [m^3_MP / (m^3_SP * s * [T]^3)], rows correspond to binding sites, cubic term of external dependence
		kA_TTT;
		EXT_BISMA_KA_TTT;
		% Desorption rate in [1 / s], rows correspond to binding sites, constant term of external dependence
		kD;
		EXT_BISMA_KD;
		% Desorption rate in [1 / (s * [T])], rows correspond to binding sites, linear term of external dependence
		kD_T;
		EXT_BISMA_KD_T;
		% Desorption rate in [1 / (s * [T]^2)], rows correspond to binding sites, quadratic term of external dependence
		kD_TT;
		EXT_BISMA_KD_TT;
		% Desorption rate in [1 / (s * [T]^3)], rows correspond to binding sites, cubic term of external dependence
		kD_TTT;
		EXT_BISMA_KD_TTT;
		% Characteristic charge [-], rows correspond to binding sites, constant term of external dependence
		nu;
		EXT_BISMA_NU;
		% Characteristic charge [1 / [T]], rows correspond to binding sites, linear term of external dependence
		nu_T;
		EXT_BISMA_NU_T;
		% Characteristic charge [1 / [T]^2], rows correspond to binding sites, quadratic term of external dependence
		nu_TT;
		EXT_BISMA_NU_TT;
		% Characteristic charge [1 / [T]^3], rows correspond to binding sites, cubic term of external dependence
		nu_TTT;
		EXT_BISMA_NU_TTT;
		% Steric factor [-], rows correspond to binding sites, constant term of external dependence
		sigma;
		EXT_BISMA_SIGMA;
		% Steric factor [1 / [T]], rows correspond to binding sites, linear term of external dependence
		sigma_T;
		EXT_BISMA_SIGMA_T;
		% Steric factor [1 / [T]^2], rows correspond to binding sites, quadratic term of external dependence
		sigma_TT;
		EXT_BISMA_SIGMA_TT;
		% Steric factor [1 / [T]^3], rows correspond to binding sites, cubic term of external dependence
		sigma_TTT;
		EXT_BISMA_SIGMA_TTT;
		% Ionic capacity in [mol / m^3_SP] for each binding site, constant term of external dependence
		lambda;
		EXT_BISMA_LAMBDA;
		% Ionic capacity in [mol / (m^3_SP * [T])] for each binding site, linear term of external dependence
		lambda_T;
		EXT_BISMA_LAMBDA_T;
		% Ionic capacity in [mol / (m^3_SP * [T]^2)] for each binding site, quadratic term of external dependence
		lambda_TT;
		EXT_BISMA_LAMBDA_TT;
		% Ionic capacity in [mol / (m^3_SP * [T]^3)] for each binding site, cubic term of external dependence
		lambda_TTT;
		EXT_BISMA_LAMBDA_TTT;
		% Reference solid phase concentrations [mol / m^3_SP]
		refSolid;
		EXT_BISMA_REFQ;
		% Reference liquid phase concentrations [mol / m^3_MP]
		refLiquid;
		EXT_BISMA_REFC0;
		% Indices of external functions (0-based)
		externalSource;
		EXTFUN;
	end

	methods

		function obj = ExtFunBiStericMassActionBinding(kA, kD)
			%EXTFUNBISTERICMASSACTIONBINDING Constructs an ExtFunBiStericMassActionBinding object with external function support
			%   OBJ = EXTFUNBISTERICMASSACTIONBINDING(KA) creates an ExtFunBiStericMassActionBinding
			%   model with the given adsorption rates KA.
			%
			%   OBJ = EXTFUNBISTERICMASSACTIONBINDING(..., KD) also sets the desorption rates to KD.

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
			validateattributes(obj.data.EXT_BISMA_KA.', {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'kA');
			validateattributes(obj.data.EXT_BISMA_KD.', {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'kD');
			validateattributes(obj.data.EXT_BISMA_NU.', {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'nu');
			validateattributes(obj.data.EXT_BISMA_SIGMA.', {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'sigma');
			validateattributes(obj.data.EXT_BISMA_LAMBDA, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', numStates}, '', 'lambda');

			validateattributes(obj.data.EXT_BISMA_KA_T.', {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'kA_T');
			validateattributes(obj.data.EXT_BISMA_KD_T.', {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'kD_T');
			validateattributes(obj.data.EXT_BISMA_NU_T.', {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'nu_T');
			validateattributes(obj.data.EXT_BISMA_SIGMA_T.', {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'sigma_T');
			validateattributes(obj.data.EXT_BISMA_LAMBDA_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', numStates}, '', 'lambda_T');

			validateattributes(obj.data.EXT_BISMA_KA_TT.', {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'kA_TT');
			validateattributes(obj.data.EXT_BISMA_KD_TT.', {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'kD_TT');
			validateattributes(obj.data.EXT_BISMA_NU_TT.', {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'nu_TT');
			validateattributes(obj.data.EXT_BISMA_SIGMA_TT.', {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'sigma_TT');
			validateattributes(obj.data.EXT_BISMA_LAMBDA_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', numStates}, '', 'lambda_TT');

			validateattributes(obj.data.EXT_BISMA_KA_TTT.', {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'kA_TTT');
			validateattributes(obj.data.EXT_BISMA_KD_TTT.', {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'kD_TTT');
			validateattributes(obj.data.EXT_BISMA_NU_TTT.', {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'nu_TTT');
			validateattributes(obj.data.EXT_BISMA_SIGMA_TTT.', {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'sigma_TTT');
			validateattributes(obj.data.EXT_BISMA_LAMBDA_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', numStates}, '', 'lambda_TTT');

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

			if ~isempty(obj.externalSource)
				validateattributes(obj.externalSource, {'double', 'int32'}, {'vector', 'nonempty', 'finite', 'real', '>=', 0}, '', 'externalSource');
				if ~any(numel(obj.externalSource) == [1, 5])
					error('CADET:invalidConfig', 'Expected externalSource to have either 1 or 5 entries but got %d entries.', numel(obj.externalSource));
				end
			end

			res = obj.validate@KineticQuasiStationaryBindingModel(nComponents, nBoundStates);
		end

		function val = get.kA(obj)
			val = obj.data.EXT_BISMA_KA.';
		end

		function set.kA(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'kA');
			obj.data.EXT_BISMA_KA = val.';
			obj.hasChanged = true;
		end

		function val = get.kD(obj)
			val = obj.data.EXT_BISMA_KD.';
		end

		function set.kD(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'kD');
			obj.data.EXT_BISMA_KD = val.';
			obj.hasChanged = true;
		end

		function val = get.nu(obj)
			val = obj.data.EXT_BISMA_NU.';
		end

		function set.nu(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'nu');
			obj.data.EXT_BISMA_NU = val.';
			obj.hasChanged = true;
		end

		function val = get.sigma(obj)
			val = obj.data.EXT_BISMA_SIGMA.';
		end

		function set.sigma(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'sigma');
			obj.data.EXT_BISMA_SIGMA = val.';
			obj.hasChanged = true;
		end

		function val = get.lambda(obj)
			val = obj.data.EXT_BISMA_LAMBDA.';
		end

		function set.lambda(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'lambda');
			obj.data.EXT_BISMA_LAMBDA = val.';
			obj.hasChanged = true;
		end

		function val = get.kA_T(obj)
			val = obj.data.EXT_BISMA_KA_T.';
		end

		function set.kA_T(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'kA_T');
			obj.data.EXT_BISMA_KA_T = val.';
			obj.hasChanged = true;
		end

		function val = get.kD_T(obj)
			val = obj.data.EXT_BISMA_KD_T.';
		end

		function set.kD_T(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'kD_T');
			obj.data.EXT_BISMA_KD_T = val.';
			obj.hasChanged = true;
		end

		function val = get.nu_T(obj)
			val = obj.data.EXT_BISMA_NU_T.';
		end

		function set.nu_T(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'nu_T');
			obj.data.EXT_BISMA_NU_T = val.';
			obj.hasChanged = true;
		end

		function val = get.sigma_T(obj)
			val = obj.data.EXT_BISMA_SIGMA_T.';
		end

		function set.sigma_T(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'sigma_T');
			obj.data.EXT_BISMA_SIGMA_T = val.';
			obj.hasChanged = true;
		end

		function val = get.lambda_T(obj)
			val = obj.data.EXT_BISMA_LAMBDA_T.';
		end

		function set.lambda_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'lambda_T');
			obj.data.EXT_BISMA_LAMBDA_T = val.';
			obj.hasChanged = true;
		end

		function val = get.kA_TT(obj)
			val = obj.data.EXT_BISMA_KA_TT.';
		end

		function set.kA_TT(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'kA_TT');
			obj.data.EXT_BISMA_KA_TT = val.';
			obj.hasChanged = true;
		end

		function val = get.kD_TT(obj)
			val = obj.data.EXT_BISMA_KD_TT.';
		end

		function set.kD_TT(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'kD_TT');
			obj.data.EXT_BISMA_KD_TT = val.';
			obj.hasChanged = true;
		end

		function val = get.nu_TT(obj)
			val = obj.data.EXT_BISMA_NU_TT.';
		end

		function set.nu_TT(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'nu_TT');
			obj.data.EXT_BISMA_NU_TT = val.';
			obj.hasChanged = true;
		end

		function val = get.sigma_TT(obj)
			val = obj.data.EXT_BISMA_SIGMA_TT.';
		end

		function set.sigma_TT(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'sigma_TT');
			obj.data.EXT_BISMA_SIGMA_TT = val.';
			obj.hasChanged = true;
		end

		function val = get.lambda_TT(obj)
			val = obj.data.EXT_BISMA_LAMBDA_TT.';
		end

		function set.lambda_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'lambda_TT');
			obj.data.EXT_BISMA_LAMBDA_TT = val.';
			obj.hasChanged = true;
		end

		function val = get.kA_TTT(obj)
			val = obj.data.EXT_BISMA_KA_TTT.';
		end

		function set.kA_TTT(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'kA_TTT');
			obj.data.EXT_BISMA_KA_TTT = val.';
			obj.hasChanged = true;
		end

		function val = get.kD_TTT(obj)
			val = obj.data.EXT_BISMA_KD_TTT.';
		end

		function set.kD_TTT(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'kD_TTT');
			obj.data.EXT_BISMA_KD_TTT = val.';
			obj.hasChanged = true;
		end

		function val = get.nu_TTT(obj)
			val = obj.data.EXT_BISMA_NU_TTT.';
		end

		function set.nu_TTT(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'nu_TTT');
			obj.data.EXT_BISMA_NU_TTT = val.';
			obj.hasChanged = true;
		end

		function val = get.sigma_TTT(obj)
			val = obj.data.EXT_BISMA_SIGMA_TTT.';
		end

		function set.sigma_TTT(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'sigma_TTT');
			obj.data.EXT_BISMA_SIGMA_TTT = val.';
			obj.hasChanged = true;
		end

		function val = get.lambda_TTT(obj)
			val = obj.data.EXT_BISMA_LAMBDA_TTT.';
		end

		function set.lambda_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'lambda_TTT');
			obj.data.EXT_BISMA_LAMBDA_TTT = val.';
			obj.hasChanged = true;
		end

		function val = get.refSolid(obj)
			if isfield(obj.data, 'EXT_BISMA_REFQ')
				val = obj.data.EXT_BISMA_REFQ;
			else
				val = [];
			end
		end

		function set.refSolid(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'EXT_BISMA_REFQ');
			else
				validateattributes(val, {'double'}, {'vector', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refSolid');
				obj.data.EXT_BISMA_REFQ = val;
			end
			obj.hasChanged = true;
		end

		function val = get.refLiquid(obj)
			if isfield(obj.data, 'EXT_BISMA_REFC0')
				val = obj.data.EXT_BISMA_REFC0;
			else
				val = [];
			end
		end

		function set.refLiquid(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'EXT_BISMA_REFC0');
			else
				validateattributes(val, {'double'}, {'vector', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refLiquid');
				obj.data.EXT_BISMA_REFC0 = val;
			end
			obj.hasChanged = true;
		end

		function val = get.EXT_BISMA_KA(obj)
			val = obj.kA;
		end
		function set.EXT_BISMA_KA(obj, val)
			obj.kA = val;
		end
		function val = get.EXT_BISMA_KD(obj)
			val = obj.kD;
		end
		function set.EXT_BISMA_KD(obj, val)
			obj.kD = val;
		end
		function val = get.EXT_BISMA_NU(obj)
			val = obj.nu;
		end
		function set.EXT_BISMA_NU(obj, val)
			obj.nu = val;
		end
		function val = get.EXT_BISMA_SIGMA(obj)
			val = obj.sigma;
		end
		function set.EXT_BISMA_SIGMA(obj, val)
			obj.sigma = val;
		end
		function val = get.EXT_BISMA_LAMBDA(obj)
			val = obj.lambda;
		end
		function set.EXT_BISMA_LAMBDA(obj, val)
			obj.lambda = val;
		end

		function val = get.EXT_BISMA_KA_T(obj)
			val = obj.kA_T;
		end
		function set.EXT_BISMA_KA_T(obj, val)
			obj.kA_T = val;
		end
		function val = get.EXT_BISMA_KD_T(obj)
			val = obj.kD_T;
		end
		function set.EXT_BISMA_KD_T(obj, val)
			obj.kD_T = val;
		end
		function val = get.EXT_BISMA_NU_T(obj)
			val = obj.nu_T;
		end
		function set.EXT_BISMA_NU_T(obj, val)
			obj.nu_T = val;
		end
		function val = get.EXT_BISMA_SIGMA_T(obj)
			val = obj.sigma_T;
		end
		function set.EXT_BISMA_SIGMA_T(obj, val)
			obj.sigma_T = val;
		end
		function val = get.EXT_BISMA_LAMBDA_T(obj)
			val = obj.lambda_T;
		end
		function set.EXT_BISMA_LAMBDA_T(obj, val)
			obj.lambda_T = val;
		end

		function val = get.EXT_BISMA_KA_TT(obj)
			val = obj.kA_TT;
		end
		function set.EXT_BISMA_KA_TT(obj, val)
			obj.kA_TT = val;
		end
		function val = get.EXT_BISMA_KD_TT(obj)
			val = obj.kD_TT;
		end
		function set.EXT_BISMA_KD_TT(obj, val)
			obj.kD_TT = val;
		end
		function val = get.EXT_BISMA_NU_TT(obj)
			val = obj.nu_TT;
		end
		function set.EXT_BISMA_NU_TT(obj, val)
			obj.nu_TT = val;
		end
		function val = get.EXT_BISMA_SIGMA_TT(obj)
			val = obj.sigma_TT;
		end
		function set.EXT_BISMA_SIGMA_TT(obj, val)
			obj.sigma_TT = val;
		end
		function val = get.EXT_BISMA_LAMBDA_TT(obj)
			val = obj.lambda_TT;
		end
		function set.EXT_BISMA_LAMBDA_TT(obj, val)
			obj.lambda_TT = val;
		end

		function val = get.EXT_BISMA_KA_TTT(obj)
			val = obj.kA_TTT;
		end
		function set.EXT_BISMA_KA_TTT(obj, val)
			obj.kA_TTT = val;
		end
		function val = get.EXT_BISMA_KD_TTT(obj)
			val = obj.kD_TTT;
		end
		function set.EXT_BISMA_KD_TTT(obj, val)
			obj.kD_TTT = val;
		end
		function val = get.EXT_BISMA_NU_TTT(obj)
			val = obj.nu_TTT;
		end
		function set.EXT_BISMA_NU_TTT(obj, val)
			obj.nu_TTT = val;
		end
		function val = get.EXT_BISMA_SIGMA_TTT(obj)
			val = obj.sigma_TTT;
		end
		function set.EXT_BISMA_SIGMA_TTT(obj, val)
			obj.sigma_TTT = val;
		end
		function val = get.EXT_BISMA_LAMBDA_TTT(obj)
			val = obj.lambda_TTT;
		end
		function set.EXT_BISMA_LAMBDA_TTT(obj, val)
			obj.lambda_TTT = val;
		end
		function val = get.EXT_BISMA_REFQ(obj)
			val = obj.refSolid;
		end
		function set.EXT_BISMA_REFQ(obj, val)
			obj.refSolid = val;
		end
		function val = get.EXT_BISMA_REFC0(obj)
			val = obj.refLiquid;
		end
		function set.EXT_BISMA_REFC0(obj, val)
			obj.refLiquid = val;
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
				if ~any(numel(val) == [1, 5])
					error('CADET:invalidConfig', 'Expected externalSource to have either 1 or 5 entries but got %d entries.', numel(val));
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

		function res = assembleConfig(obj)
			%ASSEMBLECONFIG Assembles the configuration according to the CADET file format spec
			%   RES = ASSEMBLECONFIG() returns a Matlab struct RES that represents the
			%   binding model as detailed in the CADET file format spec.
			%
			% See also BINDINGMODEL.ASSEMBLECONFIG

			res = obj.assembleConfig@KineticQuasiStationaryBindingModel();

			% Linearize
			res.EXT_BISMA_KA = obj.data.EXT_BISMA_KA(:);
			res.EXT_BISMA_KD = obj.data.EXT_BISMA_KD(:);
			res.EXT_BISMA_NU = obj.data.EXT_BISMA_NU(:);
			res.EXT_BISMA_SIGMA = obj.data.EXT_BISMA_SIGMA(:);
		end
	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = ExtFunBiStericMassActionBinding();
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
