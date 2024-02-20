
classdef ExtFunSelfAssociationBinding < KineticQuasiStationaryBindingModel
	%ExtFunSelfAssociationBinding Self association binding model with external function support
	%
	% See also BINDINGMODEL, SELFASSOCIATIONBINDING, KINETICQUASISTATIONARYBINDINGMODEL

	% Copyright: (C) 2008-2024 The CADET Authors
	%            See the license note at the end of the file.

	properties(Constant)
		name = 'EXT_SELF_ASSOCIATION'; % Name of the binding model according to CADET file format specs
	end

	properties (Dependent, Transient)
		% Adsorption rate in [m^3_MP / (m^3_SP * s)], constant term of external dependence
		kA;
		EXT_SAI_KA1;
		% Adsorption rate in [m^3_MP / (m^3_SP * s * [T])], linear term of external dependence
		kA_T;
		EXT_SAI_KA1_T;
		% Adsorption rate in [m^3_MP / (m^3_SP * s * [T]^2)], quadratic term of external dependence
		kA_TT;
		EXT_SAI_KA1_TT;
		% Adsorption rate in [m^3_MP / (m^3_SP * s * [T]^3)], cubic term of external dependence
		kA_TTT;
		EXT_SAI_KA1_TTT;
		% Adsorption rate of dimerization in [m^6_MP / (m^6_SP * s)], constant term of external dependence
		kAdimer;
		EXT_SAI_KA2;
		% Adsorption rate of dimerization in [m^6_MP / (m^6_SP * s * [T])], linear term of external dependence
		kAdimer_T;
		EXT_SAI_KA2_T;
		% Adsorption rate of dimerization in [m^6_MP / (m^6_SP * s * [T]^2)], quadratic term of external dependence
		kAdimer_TT;
		EXT_SAI_KA2_TT;
		% Adsorption rate of dimerization in [m^6_MP / (m^6_SP * s * [T]^3)], cubic term of external dependence
		kAdimer_TTT;
		EXT_SAI_KA2_TTT;
		% Desorption rate in [1 / s], constant term of external dependence
		kD;
		EXT_SAI_KD;
		% Desorption rate in [1 / (s * [T])], linear term of external dependence
		kD_T;
		EXT_SAI_KD_T;
		% Desorption rate in [1 / (s * [T]^2)], quadratic term of external dependence
		kD_TT;
		EXT_SAI_KD_TT;
		% Desorption rate in [1 / (s * [T]^3)], cubic term of external dependence
		kD_TTT;
		EXT_SAI_KD_TTT;
		% Characteristic charge [-], constant term of external dependence
		nu;
		EXT_SAI_NU;
		% Characteristic charge [1 / [T]], linear term of external dependence
		nu_T;
		EXT_SAI_NU_T;
		% Characteristic charge [1 / [T]^2], quadratic term of external dependence
		nu_TT;
		EXT_SAI_NU_TT;
		% Characteristic charge [1 / [T]^3], cubic term of external dependence
		nu_TTT;
		EXT_SAI_NU_TTT;
		% Steric factor [-], constant term of external dependence
		sigma;
		EXT_SAI_SIGMA;
		% Steric factor [1 / [T]], linear term of external dependence
		sigma_T;
		EXT_SAI_SIGMA_T;
		% Steric factor [1 / [T]^2], quadratic term of external dependence
		sigma_TT;
		EXT_SAI_SIGMA_TT;
		% Steric factor [1 / [T]^3], cubic term of external dependence
		sigma_TTT;
		EXT_SAI_SIGMA_TTT;
		% Ionic capacity in [mol / m^3_SP], constant term of external dependence
		lambda;
		EXT_SAI_LAMBDA;
		% Ionic capacity in [mol / (m^3_SP * [T])], linear term of external dependence
		lambda_T;
		EXT_SAI_LAMBDA_T;
		% Ionic capacity in [mol / (m^3_SP * [T]^2)], quadratic term of external dependence
		lambda_TT;
		EXT_SAI_LAMBDA_TT;
		% Ionic capacity in [mol / (m^3_SP * [T]^3)], cubic term of external dependence
		lambda_TTT;
		EXT_SAI_LAMBDA_TTT;
		% Reference solid phase concentration [mol / m^3_SP]
		refSolid;
		EXT_SAI_REFQ;
		% Reference liquid phase concentration [mol / m^3_MP]
		refLiquid;
		EXT_SAI_REFC0;
		% Indices of external functions (0-based)
		externalSource;
		EXTFUN;
	end

	methods

		function obj = ExtFunSelfAssociationBinding(kA, kD)
			%EXTFUNSELFASSOCIATIONBINDING Constructs an ExtFunSelfAssociationBinding object with external function support
			%   OBJ = EXTFUNSELFASSOCIATIONBINDING(KA) creates an ExtFunSelfAssociationBinding model with
			%   the given adsorption rates KA.
			%
			%   OBJ = EXTFUNSELFASSOCIATIONBINDING(..., KD) also sets the desorption rates to KD.

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
			validateattributes(obj.kAdimer, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kAdimer');
			validateattributes(obj.kD, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kD');
			validateattributes(obj.nu, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'nu');
			validateattributes(obj.sigma, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'sigma');
			validateattributes(obj.lambda, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'lambda');

			validateattributes(obj.kA_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kA_T');
			validateattributes(obj.kAdimer_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kAdimer_T');
			validateattributes(obj.kD_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kD_T');
			validateattributes(obj.nu_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'nu_T');
			validateattributes(obj.sigma_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'sigma_T');
			validateattributes(obj.lambda_T, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'lambda_T');

			validateattributes(obj.kA_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kA_TT');
			validateattributes(obj.kAdimer_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kAdimer_TT');
			validateattributes(obj.kD_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kD_TT');
			validateattributes(obj.nu_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'nu_TT');
			validateattributes(obj.sigma_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'sigma_TT');
			validateattributes(obj.lambda_TT, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'lambda_TT');

			validateattributes(obj.kA_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kA_TTT');
			validateattributes(obj.kAdimer_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kAdimer_TTT');
			validateattributes(obj.kD_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kD_TTT');
			validateattributes(obj.nu_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'nu_TTT');
			validateattributes(obj.sigma_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'sigma_TTT');
			validateattributes(obj.lambda_TTT, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'lambda_TTT');

			if ~isempty(obj.refSolid)
				validateattributes(obj.refSolid, {'double'}, {'scalar', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refSolid');
			end
			if ~isempty(obj.refLiquid)
				validateattributes(obj.refLiquid, {'double'}, {'scalar', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refLiquid');
			end

			if ~isempty(obj.externalSource)
				validateattributes(obj.externalSource, {'double', 'int32'}, {'vector', 'nonempty', 'finite', 'real', '>=', 0}, '', 'externalSource');
				if ~any(numel(obj.externalSource) == [1, 6])
					error('CADET:invalidConfig', 'Expected externalSource to have either 1 or 6 entries but got %d entries.', numel(obj.externalSource));
				end
			end

			res = obj.validate@KineticQuasiStationaryBindingModel(nComponents, nBoundStates);
		end

		function val = get.kA(obj)
			val = obj.data.EXT_SAI_KA1;
		end

		function set.kA(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kA');
			obj.data.EXT_SAI_KA1 = val;
			obj.hasChanged = true;
		end

		function val = get.kAdimer(obj)
			val = obj.data.EXT_SAI_KA2;
		end

		function set.kAdimer(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kAdimer');
			obj.data.EXT_SAI_KA2 = val;
			obj.hasChanged = true;
		end

		function val = get.kD(obj)
			val = obj.data.EXT_SAI_KD;
		end

		function set.kD(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kD');
			obj.data.EXT_SAI_KD = val;
			obj.hasChanged = true;
		end

		function val = get.nu(obj)
			val = obj.data.EXT_SAI_NU;
		end

		function set.nu(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'nu');
			obj.data.EXT_SAI_NU = val;
			obj.hasChanged = true;
		end

		function val = get.sigma(obj)
			val = obj.data.EXT_SAI_SIGMA;
		end

		function set.sigma(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'sigma');
			obj.data.EXT_SAI_SIGMA = val;
			obj.hasChanged = true;
		end

		function val = get.lambda(obj)
			val = obj.data.EXT_SAI_LAMBDA;
		end

		function set.lambda(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'lambda');
			obj.data.EXT_SAI_LAMBDA = val;
			obj.hasChanged = true;
		end

		function val = get.kA_T(obj)
			val = obj.data.EXT_SAI_KA1_T;
		end

		function set.kA_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kA_T');
			obj.data.EXT_SAI_KA1_T = val;
			obj.hasChanged = true;
		end

		function val = get.kAdimer_T(obj)
			val = obj.data.EXT_SAI_KA2_T;
		end

		function set.kAdimer_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kAdimer_T');
			obj.data.EXT_SAI_KA2_T = val;
			obj.hasChanged = true;
		end

		function val = get.kD_T(obj)
			val = obj.data.EXT_SAI_KD_T;
		end

		function set.kD_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kD_T');
			obj.data.EXT_SAI_KD_T = val;
			obj.hasChanged = true;
		end

		function val = get.nu_T(obj)
			val = obj.data.EXT_SAI_NU_T;
		end

		function set.nu_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'nu_T');
			obj.data.EXT_SAI_NU_T = val;
			obj.hasChanged = true;
		end

		function val = get.sigma_T(obj)
			val = obj.data.EXT_SAI_SIGMA_T;
		end

		function set.sigma_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'sigma_T');
			obj.data.EXT_SAI_SIGMA_T = val;
			obj.hasChanged = true;
		end

		function val = get.lambda_T(obj)
			val = obj.data.EXT_SAI_LAMBDA_T;
		end

		function set.lambda_T(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'lambda_T');
			obj.data.EXT_SAI_LAMBDA_T = val;
			obj.hasChanged = true;
		end

		function val = get.kA_TT(obj)
			val = obj.data.EXT_SAI_KA1_TT;
		end

		function set.kA_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kA_TT');
			obj.data.EXT_SAI_KA1_TT = val;
			obj.hasChanged = true;
		end

		function val = get.kAdimer_TT(obj)
			val = obj.data.EXT_SAI_KA2_TT;
		end

		function set.kAdimer_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kAdimer_TT');
			obj.data.EXT_SAI_KA2_TT = val;
			obj.hasChanged = true;
		end

		function val = get.kD_TT(obj)
			val = obj.data.EXT_SAI_KD_TT;
		end

		function set.kD_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kD_TT');
			obj.data.EXT_SAI_KD_TT = val;
			obj.hasChanged = true;
		end

		function val = get.nu_TT(obj)
			val = obj.data.EXT_SAI_NU_TT;
		end

		function set.nu_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'nu_TT');
			obj.data.EXT_SAI_NU_TT = val;
			obj.hasChanged = true;
		end

		function val = get.sigma_TT(obj)
			val = obj.data.EXT_SAI_SIGMA_TT;
		end

		function set.sigma_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'sigma_TT');
			obj.data.EXT_SAI_SIGMA_TT = val;
			obj.hasChanged = true;
		end

		function val = get.lambda_TT(obj)
			val = obj.data.EXT_SAI_LAMBDA_TT;
		end

		function set.lambda_TT(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'lambda_TT');
			obj.data.EXT_SAI_LAMBDA_TT = val;
			obj.hasChanged = true;
		end

		function val = get.kA_TTT(obj)
			val = obj.data.EXT_SAI_KA1_TTT;
		end

		function set.kA_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kA_TTT');
			obj.data.EXT_SAI_KA1_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.kAdimer_TTT(obj)
			val = obj.data.EXT_SAI_KA2_TTT;
		end

		function set.kAdimer_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kAdimer_TTT');
			obj.data.EXT_SAI_KA2_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.kD_TTT(obj)
			val = obj.data.EXT_SAI_KD_TTT;
		end

		function set.kD_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kD_TTT');
			obj.data.EXT_SAI_KD_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.nu_TTT(obj)
			val = obj.data.EXT_SAI_NU_TTT;
		end

		function set.nu_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'nu_TTT');
			obj.data.EXT_SAI_NU_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.sigma_TTT(obj)
			val = obj.data.EXT_SAI_SIGMA_TTT;
		end

		function set.sigma_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'sigma_TTT');
			obj.data.EXT_SAI_SIGMA_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.lambda_TTT(obj)
			val = obj.data.EXT_SAI_LAMBDA_TTT;
		end

		function set.lambda_TTT(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'lambda_TTT');
			obj.data.EXT_SAI_LAMBDA_TTT = val;
			obj.hasChanged = true;
		end
		
		function val = get.refSolid(obj)
			if isfield(obj.data, 'EXT_SAI_REFQ')
				val = obj.data.EXT_SAI_REFQ;
			else
				val = [];
			end
		end

		function set.refSolid(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'EXT_SAI_REFQ');
			else
				validateattributes(val, {'double'}, {'scalar', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refSolid');
				obj.data.EXT_SAI_REFQ = val;
			end
			obj.hasChanged = true;
		end

		function val = get.refLiquid(obj)
			if isfield(obj.data, 'EXT_SAI_REFC0')
				val = obj.data.EXT_SAI_REFC0;
			else
				val = [];
			end
		end

		function set.refLiquid(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'EXT_SAI_REFC0');
			else
				validateattributes(val, {'double'}, {'scalar', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refLiquid');
				obj.data.EXT_SAI_REFC0 = val;
			end
			obj.hasChanged = true;
		end

		function val = get.EXT_SAI_KA1(obj)
			val = obj.kA;
		end
		function set.EXT_SAI_KA1(obj, val)
			obj.kA = val;
		end
		function val = get.EXT_SAI_KA2(obj)
			val = obj.kAdimer;
		end
		function set.EXT_SAI_KA2(obj, val)
			obj.kAdimer = val;
		end
		function val = get.EXT_SAI_KD(obj)
			val = obj.kD;
		end
		function set.EXT_SAI_KD(obj, val)
			obj.kD = val;
		end
		function val = get.EXT_SAI_LAMBDA(obj)
			val = obj.lambda;
		end
		function set.EXT_SAI_LAMBDA(obj, val)
			obj.lambda = val;
		end
		function val = get.EXT_SAI_NU(obj)
			val = obj.nu;
		end
		function set.EXT_SAI_NU(obj, val)
			obj.nu = val;
		end
		function val = get.EXT_SAI_SIGMA(obj)
			val = obj.sigma;
		end
		function set.EXT_SAI_SIGMA(obj, val)
			obj.sigma = val;
		end

		function val = get.EXT_SAI_KA1_T(obj)
			val = obj.kA_T;
		end
		function set.EXT_SAI_KA1_T(obj, val)
			obj.kA_T = val;
		end
		function val = get.EXT_SAI_KA2_T(obj)
			val = obj.kAdimer_T;
		end
		function set.EXT_SAI_KA2_T(obj, val)
			obj.kAdimer_T = val;
		end
		function val = get.EXT_SAI_KD_T(obj)
			val = obj.kD_T;
		end
		function set.EXT_SAI_KD_T(obj, val)
			obj.kD_T = val;
		end
		function val = get.EXT_SAI_LAMBDA_T(obj)
			val = obj.lambda_T;
		end
		function set.EXT_SAI_LAMBDA_T(obj, val)
			obj.lambda_T = val;
		end
		function val = get.EXT_SAI_NU_T(obj)
			val = obj.nu_T;
		end
		function set.EXT_SAI_NU_T(obj, val)
			obj.nu_T = val;
		end
		function val = get.EXT_SAI_SIGMA_T(obj)
			val = obj.sigma_T;
		end
		function set.EXT_SAI_SIGMA_T(obj, val)
			obj.sigma_T = val;
		end

		function val = get.EXT_SAI_KA1_TT(obj)
			val = obj.kA_TT;
		end
		function set.EXT_SAI_KA1_TT(obj, val)
			obj.kA_TT = val;
		end
		function val = get.EXT_SAI_KA2_TT(obj)
			val = obj.kAdimer_TT;
		end
		function set.EXT_SAI_KA2_TT(obj, val)
			obj.kAdimer_TT = val;
		end
		function val = get.EXT_SAI_KD_TT(obj)
			val = obj.kD_TT;
		end
		function set.EXT_SAI_KD_TT(obj, val)
			obj.kD_TT = val;
		end
		function val = get.EXT_SAI_LAMBDA_TT(obj)
			val = obj.lambda_TT;
		end
		function set.EXT_SAI_LAMBDA_TT(obj, val)
			obj.lambda_TT = val;
		end
		function val = get.EXT_SAI_NU_TT(obj)
			val = obj.nu_TT;
		end
		function set.EXT_SAI_NU_TT(obj, val)
			obj.nu_TT = val;
		end
		function val = get.EXT_SAI_SIGMA_TT(obj)
			val = obj.sigma_TT;
		end
		function set.EXT_SAI_SIGMA_TT(obj, val)
			obj.sigma_TT = val;
		end

		function val = get.EXT_SAI_KA1_TTT(obj)
			val = obj.kA_TTT;
		end
		function set.EXT_SAI_KA1_TTT(obj, val)
			obj.kA_TTT = val;
		end
		function val = get.EXT_SAI_KA2_TTT(obj)
			val = obj.kAdimer_TTT;
		end
		function set.EXT_SAI_KA2_TTT(obj, val)
			obj.kAdimer_TTT = val;
		end
		function val = get.EXT_SAI_KD_TTT(obj)
			val = obj.kD_TTT;
		end
		function set.EXT_SAI_KD_TTT(obj, val)
			obj.kD_TTT = val;
		end
		function val = get.EXT_SAI_LAMBDA_TTT(obj)
			val = obj.lambda_TTT;
		end
		function set.EXT_SAI_LAMBDA_TTT(obj, val)
			obj.lambda_TTT = val;
		end
		function val = get.EXT_SAI_NU_TTT(obj)
			val = obj.nu_TTT;
		end
		function set.EXT_SAI_NU_TTT(obj, val)
			obj.nu_TTT = val;
		end
		function val = get.EXT_SAI_SIGMA_TTT(obj)
			val = obj.sigma_TTT;
		end
		function set.EXT_SAI_SIGMA_TTT(obj, val)
			obj.sigma_TTT = val;
		end

		function val = get.EXT_SAI_REFQ(obj)
			val = obj.refSolid;
		end
		function set.EXT_SAI_REFQ(obj, val)
			obj.refSolid = val;
		end
		function val = get.EXT_SAI_REFC0(obj)
			val = obj.refLiquid;
		end
		function set.EXT_SAI_REFC0(obj, val)
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
				if ~any(numel(val) == [1, 6])
					error('CADET:invalidConfig', 'Expected externalSource to have either 1 or 6 entries but got %d entries.', numel(val));
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
				obj = ExtFunSelfAssociationBinding();
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
