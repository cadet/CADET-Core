
classdef ExtFunGeneralizedIonExchangeBinding < KineticQuasiStationaryBindingModel
	%ExtFunGeneralizedIonExchangeBinding Generalized ion exchange binding model with external function support
	%
	% See also BINDINGMODEL, GENERALIZEDIONEXCHANGEBINDING, KINETICQUASISTATIONARYBINDINGMODEL

	% Copyright: (C) 2008-2024 The CADET Authors
	%            See the license note at the end of the file.

	properties(Constant)
		name = 'EXT_GENERALIZED_ION_EXCHANGE'; % Name of the binding model according to CADET file format specs
	end

	properties (Dependent, Transient)
		% Adsorption rate in [m^3_MP / (m^3_SP * s)], constant term of external dependence
		kA;
		EXT_GIEX_KA;
		% Adsorption rate in [m^3_MP / (m^3_SP * s * [T])], linear term of external dependence
		kA_T;
		EXT_GIEX_KA_T;
		% Adsorption rate in [m^3_MP / (m^3_SP * s * [T]^2)], quadratic term of external dependence
		kA_TT;
		EXT_GIEX_KA_TT;
		% Adsorption rate in [m^3_MP / (m^3_SP * s * [T]^3)], cubic term of external dependence
		kA_TTT;
		EXT_GIEX_KA_TTT;
		% Adsorption rate linear modifier coefficient in [1 / Mod], constant term of external dependence
		kALin;
		EXT_GIEX_KA_LIN;
		% Adsorption rate linear modifier coefficient in [1 / ([Mod] * [T])], linear term of external dependence
		kALin_T;
		EXT_GIEX_KA_LIN_T;
		% Adsorption rate linear modifier coefficient in [1 / ([Mod] * [T]^2)], quadratic term of external dependence
		kALin_TT;
		EXT_GIEX_KA_LIN_TT;
		% Adsorption rate linear modifier coefficient in [1 / ([Mod] * [T]^3)], cubic term of external dependence
		kALin_TTT;
		EXT_GIEX_KA_LIN_TTT;
		% Adsorption rate quadratic modifier coefficient in [1 / Mod^2], constant term of external dependence
		kAQuad;
		EXT_GIEX_KA_QUAD;
		% Adsorption rate quadratic modifier coefficient in [1 / ([Mod]^2 * [T])], linear term of external dependence
		kAQuad_T;
		EXT_GIEX_KA_QUAD_T;
		% Adsorption rate quadratic modifier coefficient in [1 / ([Mod]^2 * [T]^2)], quadratic term of external dependence
		kAQuad_TT;
		EXT_GIEX_KA_QUAD_TT;
		% Adsorption rate quadratic modifier coefficient in [1 / ([Mod]^2 * [T]^3)], cubic term of external dependence
		kAQuad_TTT;
		EXT_GIEX_KA_QUAD_TTT;
		% Adsorption rate protein-protein modifier coefficient in [m^3_MP / mol], constant term of external dependence
		kAProt;
		EXT_GIEX_KA_PROT;
		% Adsorption rate protein-protein modifier coefficient in [m^3_MP / (mol * [T])], linear term of external dependence
		kAProt_T;
		EXT_GIEX_KA_PROT_T;
		% Adsorption rate protein-protein modifier coefficient in [m^3_MP / (mol * [T]^2)], quadratic term of external dependence
		kAProt_TT;
		EXT_GIEX_KA_PROT_TT;
		% Adsorption rate protein-protein modifier coefficient in [m^3_MP / (mol * [T]^3)], cubic term of external dependence
		kAProt_TTT;
		EXT_GIEX_KA_PROT_TTT;
		% Adsorption rate salt modifier coefficient [-], constant term of external dependence
		kASalt;
		EXT_GIEX_KA_SALT;
		% Adsorption rate salt modifier coefficient in [1 / ([T])], linear term of external dependence
		kASalt_T;
		EXT_GIEX_KA_SALT_T;
		% Adsorption rate salt modifier coefficient in [1 / ([T]^2)], quadratic term of external dependence
		kASalt_TT;
		EXT_GIEX_KA_SALT_TT;
		% Adsorption rate salt modifier coefficient in [1 / ([T]^3)], cubic term of external dependence
		kASalt_TTT;
		EXT_GIEX_KA_SALT_TTT;
		% Desorption rate in [1 / s], constant term of external dependence
		kD;
		EXT_GIEX_KD;
		% Desorption rate in [1 / (s * [T])], linear term of external dependence
		kD_T;
		EXT_GIEX_KD_T;
		% Desorption rate in [1 / (s * [T]^2)], quadratic term of external dependence
		kD_TT;
		EXT_GIEX_KD_TT;
		% Desorption rate in [1 / (s * [T]^3)], cubic term of external dependence
		kD_TTT;
		EXT_GIEX_KD_TTT;
		% Desorption rate linear modifier coefficient in [1 / Mod], constant term of external dependence
		kDLin;
		EXT_GIEX_KD_LIN;
		% Desorption rate linear modifier coefficient in [1 / ([Mod] * [T])], linear term of external dependence
		kDLin_T;
		EXT_GIEX_KD_LIN_T;
		% Desorption rate linear modifier coefficient in [1 / ([Mod] * [T]^2)], quadratic term of external dependence
		kDLin_TT;
		EXT_GIEX_KD_LIN_TT;
		% Desorption rate linear modifier coefficient in [1 / ([Mod] * [T]^3)], cubic term of external dependence
		kDLin_TTT;
		EXT_GIEX_KD_LIN_TTT;
		% Desorption rate quadratic modifier coefficient in [1 / Mod^2], constant term of external dependence
		kDQuad;
		EXT_GIEX_KD_QUAD;
		% Desorption rate quadratic modifier coefficient in [1 / ([Mod]^2 * [T])], linear term of external dependence
		kDQuad_T;
		EXT_GIEX_KD_QUAD_T;
		% Desorption rate quadratic modifier coefficient in [1 / ([Mod]^2 * [T]^2)], quadratic term of external dependence
		kDQuad_TT;
		EXT_GIEX_KD_QUAD_TT;
		% Desorption rate quadratic modifier coefficient in [1 / ([Mod]^2 * [T]^3)], cubic term of external dependence
		kDQuad_TTT;
		EXT_GIEX_KD_QUAD_TTT;
		% Desorption rate protein-protein modifier coefficient in [m^3_MP / mol], constant term of external dependence
		kDProt;
		EXT_GIEX_KD_PROT;
		% Desorption rate protein-protein modifier coefficient in [m^3_MP / (mol * [T])], linear term of external dependence
		kDProt_T;
		EXT_GIEX_KD_PROT_T;
		% Desorption rate protein-protein modifier coefficient in [m^3_MP / (mol * [T]^2)], quadratic term of external dependence
		kDProt_TT;
		EXT_GIEX_KD_PROT_TT;
		% Desorption rate protein-protein modifier coefficient in [m^3_MP / (mol * [T]^3)], cubic term of external dependence
		kDProt_TTT;
		EXT_GIEX_KD_PROT_TTT;
		% Desorption rate salt modifier coefficient [-], constant term of external dependence
		kDSalt;
		EXT_GIEX_KD_SALT;
		% Desorption rate salt modifier coefficient in [1 / ([T])], linear term of external dependence
		kDSalt_T;
		EXT_GIEX_KD_SALT_T;
		% Desorption rate salt modifier coefficient in [1 / ([T]^2)], quadratic term of external dependence
		kDSalt_TT;
		EXT_GIEX_KD_SALT_TT;
		% Desorption rate salt modifier coefficient in [1 / ([T]^3)], cubic term of external dependence
		kDSalt_TTT;
		EXT_GIEX_KD_SALT_TTT;
		% Characteristic charge [-], constant term of external dependence
		nu;
		EXT_GIEX_NU;
		% Characteristic charge in [1 / ([T])], linear term of external dependence
		nu_T;
		EXT_GIEX_NU_T;
		% Characteristic charge in [1 / ([T]^2)], quadratic term of external dependence
		nu_TT;
		EXT_GIEX_NU_TT;
		% Characteristic charge in [1 / ([T]^3)], cubic term of external dependence
		nu_TTT;
		EXT_GIEX_NU_TTT;
		% Linear dependence on modifier of characteristic charge [1 / Mod], constant term of external dependence
		nuLin;
		EXT_GIEX_NU_LIN;
		% Linear dependence on modifier of characteristic charge [1 / ([Mod] * [T])], linear term of external dependence
		nuLin_T;
		EXT_GIEX_NU_LIN_T;
		% Linear dependence on modifier of characteristic charge [1 / ([Mod] * [T]^2)], quadratic term of external dependence
		nuLin_TT;
		EXT_GIEX_NU_LIN_TT;
		% Linear dependence on modifier of characteristic charge [1 / ([Mod] * [T]^3)], cubic term of external dependence
		nuLin_TTT;
		EXT_GIEX_NU_LIN_TTT;
		% Quadratic dependence on modifier of characteristic charge [1 / Mod^2], constant term of external dependence
		nuQuad;
		EXT_GIEX_NU_QUAD;
		% Quadratic dependence on modifier of characteristic charge [1 / ([Mod]^2 * [T])], linear term of external dependence
		nuQuad_T;
		EXT_GIEX_NU_QUAD_T;
		% Quadratic dependence on modifier of characteristic charge [1 / ([Mod]^2 * [T]^2)], quadratic term of external dependence
		nuQuad_TT;
		EXT_GIEX_NU_QUAD_TT;
		% Quadratic dependence on modifier of characteristic charge [1 / ([Mod]^2 * [T]^3)], cubic term of external dependence
		nuQuad_TTT;
		EXT_GIEX_NU_QUAD_TTT;
		% Steric factor [-], constant term of external dependence
		sigma;
		EXT_GIEX_SIGMA;
		% Steric factor [1 / [T]], linear term of external dependence
		sigma_T;
		EXT_GIEX_SIGMA_T;
		% Steric factor [1 / [T]^2], quadratic term of external dependence
		sigma_TT;
		EXT_GIEX_SIGMA_TT;
		% Steric factor [1 / [T]^3], cubic term of external dependence
		sigma_TTT;
		EXT_GIEX_SIGMA_TTT;
		% Ionic capacity in [mol / m^3_SP], constant term of external dependence
		lambda;
		EXT_GIEX_LAMBDA;
		% Ionic capacity in [mol / (m^3_SP * [T])], linear term of external dependence
		lambda_T;
		EXT_GIEX_LAMBDA_T;
		% Ionic capacity in [mol / (m^3_SP * [T]^2)], quadratic term of external dependence
		lambda_TT;
		EXT_GIEX_LAMBDA_TT;
		% Ionic capacity in [mol / (m^3_SP * [T]^3)], cubic term of external dependence
		lambda_TTT;
		EXT_GIEX_LAMBDA_TTT;
		% Reference solid phase concentration [mol / m^3_SP]
		refSolid;
		EXT_GIEX_REFQ;
		% Reference liquid phase concentration [mol / m^3_MP]
		refLiquid;
		EXT_GIEX_REFC0;
	end

	methods

		function obj = ExtFunGeneralizedIonExchangeBinding(kA, kD)
			%EXTFUNGENERALIZEDIONEXCHANGEBINDING Constructs an ExtFunGeneralizedIonExchangeBinding object
			%   OBJ = EXTFUNGENERALIZEDIONEXCHANGEBINDING(KA) creates an ExtFunGeneralizedIonExchangeBinding model with
			%   the given adsorption rates KA.
			%
			%   OBJ = EXTFUNGENERALIZEDIONEXCHANGEBINDING(..., KD) also sets the desorption rates to KD.

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

			if nBoundStates(1) ~= 1
				error('CADET:invalidConfig', 'ExtFunGeneralizedIonExchangeBinding requires first component (salt) to have one bound state');
			end

			if nBoundStates(2) ~= 0
				error('CADET:invalidConfig', 'ExtFunGeneralizedIonExchangeBinding requires second component (modifier) to be non-binding');
			end

			validateattributes(obj.kA, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kA');
			validateattributes(obj.kALin, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kALin');
			validateattributes(obj.kAQuad, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kAQuad');
			validateattributes(obj.kAProt, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kAProt');
			validateattributes(obj.kASalt, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kASalt');
			validateattributes(obj.kD, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kD');
			validateattributes(obj.kDLin, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kDLin');
			validateattributes(obj.kDQuad, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kDQuad');
			validateattributes(obj.kDProt, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kDProt');
			validateattributes(obj.kDSalt, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kDSalt');
			validateattributes(obj.nu, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'nu');
			validateattributes(obj.nuLin, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'nuLin');
			validateattributes(obj.nuQuad, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'nuQuad');
			validateattributes(obj.sigma, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'sigma');
			validateattributes(obj.lambda, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'lambda');

			validateattributes(obj.kA_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kA_T');
			validateattributes(obj.kALin_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kALin_T');
			validateattributes(obj.kAQuad_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kAQuad_T');
			validateattributes(obj.kAProt_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kAProt_T');
			validateattributes(obj.kASalt_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kASalt_T');
			validateattributes(obj.kD_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kD_T');
			validateattributes(obj.kDLin_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kDLin_T');
			validateattributes(obj.kDQuad_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kDQuad_T');
			validateattributes(obj.kDProt_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kDProt_T');
			validateattributes(obj.kDSalt_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kDSalt_T');
			validateattributes(obj.nu_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'nu_T');
			validateattributes(obj.nuLin_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'nuLin_T');
			validateattributes(obj.nuQuad_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'nuQuad_T');
			validateattributes(obj.sigma_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'sigma_T');
			validateattributes(obj.lambda_T, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'lambda_T');

			validateattributes(obj.kA_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kA_TT');
			validateattributes(obj.kALin_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kALin_TT');
			validateattributes(obj.kAQuad_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kAQuad_TT');
			validateattributes(obj.kAProt_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kAProt_TT');
			validateattributes(obj.kASalt_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kASalt_TT');
			validateattributes(obj.kD_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kD_TT');
			validateattributes(obj.kDLin_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kDLin_TT');
			validateattributes(obj.kDQuad_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kDQuad_TT');
			validateattributes(obj.kDProt_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kDProt_TT');
			validateattributes(obj.kDSalt_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kDSalt_TT');
			validateattributes(obj.nu_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'nu_TT');
			validateattributes(obj.nuLin_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'nuLin_TT');
			validateattributes(obj.nuQuad_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'nuQuad_TT');
			validateattributes(obj.sigma_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'sigma_TT');
			validateattributes(obj.lambda_TT, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'lambda_TT');

			validateattributes(obj.kA_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kA_TTT');
			validateattributes(obj.kALin_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kALin_TTT');
			validateattributes(obj.kAQuad_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kAQuad_TTT');
			validateattributes(obj.kAProt_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kAProt_TTT');
			validateattributes(obj.kASalt_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kASalt_TTT');
			validateattributes(obj.kD_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kD_TTT');
			validateattributes(obj.kDLin_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kDLin_TTT');
			validateattributes(obj.kDQuad_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kDQuad_TTT');
			validateattributes(obj.kDProt_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kDProt_TTT');
			validateattributes(obj.kDSalt_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kDSalt_TTT');
			validateattributes(obj.nu_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'nu_TTT');
			validateattributes(obj.nuLin_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'nuLin_TTT');
			validateattributes(obj.nuQuad_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'nuQuad_TTT');
			validateattributes(obj.sigma_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'sigma_TTT');
			validateattributes(obj.lambda_TTT, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'lambda_TTT');

			if ~isempty(obj.refSolid)
				validateattributes(obj.refSolid, {'double'}, {'scalar', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refSolid');
			end
			if ~isempty(obj.refLiquid)
				validateattributes(obj.refLiquid, {'double'}, {'scalar', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refLiquid');
			end

			res = obj.validate@KineticQuasiStationaryBindingModel(nComponents, nBoundStates);
		end

		function val = get.kA(obj)
			val = obj.data.EXT_GIEX_KA;
		end

		function set.kA(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kA');
			obj.data.EXT_GIEX_KA = val;
			obj.hasChanged = true;
		end

		function val = get.kALin(obj)
			val = obj.data.EXT_GIEX_KA_LIN;
		end

		function set.kALin(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kALin');
			obj.data.EXT_GIEX_KA_LIN = val;
			obj.hasChanged = true;
		end

		function val = get.kAQuad(obj)
			val = obj.data.EXT_GIEX_KA_QUAD;
		end

		function set.kAQuad(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kAQuad');
			obj.data.EXT_GIEX_KA_QUAD = val;
			obj.hasChanged = true;
		end

		function val = get.kASalt(obj)
			val = obj.data.EXT_GIEX_KA_SALT;
		end

		function set.kASalt(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kASalt');
			obj.data.EXT_GIEX_KA_SALT = val;
			obj.hasChanged = true;
		end

		function val = get.kAProt(obj)
			val = obj.data.EXT_GIEX_KA_PROT;
		end

		function set.kAProt(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kAProt');
			obj.data.EXT_GIEX_KA_PROT = val;
			obj.hasChanged = true;
		end

		function val = get.kD(obj)
			val = obj.data.EXT_GIEX_KD;
		end

		function set.kD(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kD');
			obj.data.EXT_GIEX_KD = val;
			obj.hasChanged = true;
		end

		function val = get.kDLin(obj)
			val = obj.data.EXT_GIEX_KD_LIN;
		end

		function set.kDLin(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kDLin');
			obj.data.EXT_GIEX_KD_LIN = val;
			obj.hasChanged = true;
		end

		function val = get.kDQuad(obj)
			val = obj.data.EXT_GIEX_KD_QUAD;
		end

		function set.kDQuad(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kDQuad');
			obj.data.EXT_GIEX_KD_QUAD = val;
			obj.hasChanged = true;
		end

		function val = get.kDSalt(obj)
			val = obj.data.EXT_GIEX_KD_SALT;
		end

		function set.kDSalt(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kDSalt');
			obj.data.EXT_GIEX_KD_SALT = val;
			obj.hasChanged = true;
		end

		function val = get.kDProt(obj)
			val = obj.data.EXT_GIEX_KD_PROT;
		end

		function set.kDProt(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kDProt');
			obj.data.EXT_GIEX_KD_PROT = val;
			obj.hasChanged = true;
		end

		function val = get.nu(obj)
			val = obj.data.EXT_GIEX_NU;
		end

		function set.nu(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'nu');
			obj.data.EXT_GIEX_NU = val;
			obj.hasChanged = true;
		end

		function val = get.nuLin(obj)
			val = obj.data.EXT_GIEX_NU_LIN;
		end

		function set.nuLin(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'nuLin');
			obj.data.EXT_GIEX_NU_LIN = val;
			obj.hasChanged = true;
		end

		function val = get.nuQuad(obj)
			val = obj.data.EXT_GIEX_NU_QUAD;
		end

		function set.nuQuad(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'nuQuad');
			obj.data.EXT_GIEX_NU_QUAD = val;
			obj.hasChanged = true;
		end

		function val = get.sigma(obj)
			val = obj.data.EXT_GIEX_SIGMA;
		end

		function set.sigma(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'sigma');
			obj.data.EXT_GIEX_SIGMA = val;
			obj.hasChanged = true;
		end

		function val = get.lambda(obj)
			val = obj.data.EXT_GIEX_LAMBDA;
		end

		function set.lambda(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'lambda');
			obj.data.EXT_GIEX_LAMBDA = val;
			obj.hasChanged = true;
		end

		function val = get.kA_T(obj)
			val = obj.data.EXT_GIEX_KA_T;
		end

		function set.kA_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kA_T');
			obj.data.EXT_GIEX_KA_T = val;
			obj.hasChanged = true;
		end

		function val = get.kALin_T(obj)
			val = obj.data.EXT_GIEX_KA_LIN_T;
		end

		function set.kALin_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kALin_T');
			obj.data.EXT_GIEX_KA_LIN_T = val;
			obj.hasChanged = true;
		end

		function val = get.kAQuad_T(obj)
			val = obj.data.EXT_GIEX_KA_QUAD_T;
		end

		function set.kAQuad_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kAQuad_T');
			obj.data.EXT_GIEX_KA_QUAD_T = val;
			obj.hasChanged = true;
		end

		function val = get.kASalt_T(obj)
			val = obj.data.EXT_GIEX_KA_SALT_T;
		end

		function set.kASalt_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kASalt_T');
			obj.data.EXT_GIEX_KA_SALT_T = val;
			obj.hasChanged = true;
		end

		function val = get.kAProt_T(obj)
			val = obj.data.EXT_GIEX_KA_PROT_T;
		end

		function set.kAProt_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kAProt_T');
			obj.data.EXT_GIEX_KA_PROT_T = val;
			obj.hasChanged = true;
		end

		function val = get.kD_T(obj)
			val = obj.data.EXT_GIEX_KD_T;
		end

		function set.kD_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kD_T');
			obj.data.EXT_GIEX_KD_T = val;
			obj.hasChanged = true;
		end

		function val = get.kDLin_T(obj)
			val = obj.data.EXT_GIEX_KD_LIN_T;
		end

		function set.kDLin_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kDLin_T');
			obj.data.EXT_GIEX_KD_LIN_T = val;
			obj.hasChanged = true;
		end

		function val = get.kDQuad_T(obj)
			val = obj.data.EXT_GIEX_KD_QUAD_T;
		end

		function set.kDQuad_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kDQuad_T');
			obj.data.EXT_GIEX_KD_QUAD_T = val;
			obj.hasChanged = true;
		end

		function val = get.kDSalt_T(obj)
			val = obj.data.EXT_GIEX_KD_SALT_T;
		end

		function set.kDSalt_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kDSalt_T');
			obj.data.EXT_GIEX_KD_SALT_T = val;
			obj.hasChanged = true;
		end

		function val = get.kDProt_T(obj)
			val = obj.data.EXT_GIEX_KD_PROT_T;
		end

		function set.kDProt_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kDProt_T');
			obj.data.EXT_GIEX_KD_PROT_T = val;
			obj.hasChanged = true;
		end

		function val = get.nu_T(obj)
			val = obj.data.EXT_GIEX_NU_T;
		end

		function set.nu_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'nu_T');
			obj.data.EXT_GIEX_NU_T = val;
			obj.hasChanged = true;
		end

		function val = get.nuLin_T(obj)
			val = obj.data.EXT_GIEX_NU_LIN_T;
		end

		function set.nuLin_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'nuLin_T');
			obj.data.EXT_GIEX_NU_LIN_T = val;
			obj.hasChanged = true;
		end

		function val = get.nuQuad_T(obj)
			val = obj.data.EXT_GIEX_NU_QUAD_T;
		end

		function set.nuQuad_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'nuQuad_T');
			obj.data.EXT_GIEX_NU_QUAD_T = val;
			obj.hasChanged = true;
		end

		function val = get.sigma_T(obj)
			val = obj.data.EXT_GIEX_SIGMA_T;
		end

		function set.sigma_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'sigma_T');
			obj.data.EXT_GIEX_SIGMA_T = val;
			obj.hasChanged = true;
		end

		function val = get.lambda_T(obj)
			val = obj.data.EXT_GIEX_LAMBDA_T;
		end

		function set.lambda_T(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'lambda_T');
			obj.data.EXT_GIEX_LAMBDA_T = val;
			obj.hasChanged = true;
		end

		function val = get.kA_TT(obj)
			val = obj.data.EXT_GIEX_KA_TT;
		end

		function set.kA_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kA_TT');
			obj.data.EXT_GIEX_KA_TT = val;
			obj.hasChanged = true;
		end

		function val = get.kALin_TT(obj)
			val = obj.data.EXT_GIEX_KA_LIN_TT;
		end

		function set.kALin_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kALin_TT');
			obj.data.EXT_GIEX_KA_LIN_TT = val;
			obj.hasChanged = true;
		end

		function val = get.kAQuad_TT(obj)
			val = obj.data.EXT_GIEX_KA_QUAD_TT;
		end

		function set.kAQuad_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kAQuad_TT');
			obj.data.EXT_GIEX_KA_QUAD_TT = val;
			obj.hasChanged = true;
		end

		function val = get.kASalt_TT(obj)
			val = obj.data.EXT_GIEX_KA_SALT_TT;
		end

		function set.kASalt_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kASalt_TT');
			obj.data.EXT_GIEX_KA_SALT_TT = val;
			obj.hasChanged = true;
		end

		function val = get.kAProt_TT(obj)
			val = obj.data.EXT_GIEX_KA_PROT_TT;
		end

		function set.kAProt_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kAProt_TT');
			obj.data.EXT_GIEX_KA_PROT_TT = val;
			obj.hasChanged = true;
		end

		function val = get.kD_TT(obj)
			val = obj.data.EXT_GIEX_KD_TT;
		end

		function set.kD_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kD_TT');
			obj.data.EXT_GIEX_KD_TT = val;
			obj.hasChanged = true;
		end

		function val = get.kDLin_TT(obj)
			val = obj.data.EXT_GIEX_KD_LIN_TT;
		end

		function set.kDLin_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kDLin_TT');
			obj.data.EXT_GIEX_KD_LIN_TT = val;
			obj.hasChanged = true;
		end

		function val = get.kDQuad_TT(obj)
			val = obj.data.EXT_GIEX_KD_QUAD_TT;
		end

		function set.kDQuad_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kDQuad_TT');
			obj.data.EXT_GIEX_KD_QUAD_TT = val;
			obj.hasChanged = true;
		end

		function val = get.kDSalt_TT(obj)
			val = obj.data.EXT_GIEX_KD_SALT_TT;
		end

		function set.kDSalt_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kDSalt_TT');
			obj.data.EXT_GIEX_KD_SALT_TT = val;
			obj.hasChanged = true;
		end

		function val = get.kDProt_TT(obj)
			val = obj.data.EXT_GIEX_KD_PROT_TT;
		end

		function set.kDProt_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kDProt_TT');
			obj.data.EXT_GIEX_KD_PROT_TT = val;
			obj.hasChanged = true;
		end

		function val = get.nu_TT(obj)
			val = obj.data.EXT_GIEX_NU_TT;
		end

		function set.nu_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'nu_TT');
			obj.data.EXT_GIEX_NU_TT = val;
			obj.hasChanged = true;
		end

		function val = get.nuLin_TT(obj)
			val = obj.data.EXT_GIEX_NU_LIN_TT;
		end

		function set.nuLin_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'nuLin_TT');
			obj.data.EXT_GIEX_NU_LIN_TT = val;
			obj.hasChanged = true;
		end

		function val = get.nuQuad_TT(obj)
			val = obj.data.EXT_GIEX_NU_QUAD_TT;
		end

		function set.nuQuad_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'nuQuad_TT');
			obj.data.EXT_GIEX_NU_QUAD_TT = val;
			obj.hasChanged = true;
		end

		function val = get.sigma_TT(obj)
			val = obj.data.EXT_GIEX_SIGMA_TT;
		end

		function set.sigma_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'sigma_TT');
			obj.data.EXT_GIEX_SIGMA_TT = val;
			obj.hasChanged = true;
		end

		function val = get.lambda_TT(obj)
			val = obj.data.EXT_GIEX_LAMBDA_TT;
		end

		function set.lambda_TT(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'lambda_TT');
			obj.data.EXT_GIEX_LAMBDA_TT = val;
			obj.hasChanged = true;
		end

		function val = get.kA_TTT(obj)
			val = obj.data.EXT_GIEX_KA_TTT;
		end

		function set.kA_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kA_TTT');
			obj.data.EXT_GIEX_KA_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.kALin_TTT(obj)
			val = obj.data.EXT_GIEX_KA_LIN_TTT;
		end

		function set.kALin_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kALin_TTT');
			obj.data.EXT_GIEX_KA_LIN_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.kAQuad_TTT(obj)
			val = obj.data.EXT_GIEX_KA_QUAD_TTT;
		end

		function set.kAQuad_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kAQuad_TTT');
			obj.data.EXT_GIEX_KA_QUAD_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.kASalt_TTT(obj)
			val = obj.data.EXT_GIEX_KA_SALT_TTT;
		end

		function set.kASalt_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kASalt_TTT');
			obj.data.EXT_GIEX_KA_SALT_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.kAProt_TTT(obj)
			val = obj.data.EXT_GIEX_KA_PROT_TTT;
		end

		function set.kAProt_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kAProt_TTT');
			obj.data.EXT_GIEX_KA_PROT_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.kD_TTT(obj)
			val = obj.data.EXT_GIEX_KD_TTT;
		end

		function set.kD_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kD_TTT');
			obj.data.EXT_GIEX_KD_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.kDLin_TTT(obj)
			val = obj.data.EXT_GIEX_KD_LIN_TTT;
		end

		function set.kDLin_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kDLin_TTT');
			obj.data.EXT_GIEX_KD_LIN_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.kDQuad_TTT(obj)
			val = obj.data.EXT_GIEX_KD_QUAD_TTT;
		end

		function set.kDQuad_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kDQuad_TTT');
			obj.data.EXT_GIEX_KD_QUAD_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.kDSalt_TTT(obj)
			val = obj.data.EXT_GIEX_KD_SALT_TTT;
		end

		function set.kDSalt_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kDSalt_TTT');
			obj.data.EXT_GIEX_KD_SALT_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.kDProt_TTT(obj)
			val = obj.data.EXT_GIEX_KD_PROT_TTT;
		end

		function set.kDProt_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kDProt_TTT');
			obj.data.EXT_GIEX_KD_PROT_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.nu_TTT(obj)
			val = obj.data.EXT_GIEX_NU_TTT;
		end

		function set.nu_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'nu_TTT');
			obj.data.EXT_GIEX_NU_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.nuLin_TTT(obj)
			val = obj.data.EXT_GIEX_NU_LIN_TTT;
		end

		function set.nuLin_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'nuLin_TTT');
			obj.data.EXT_GIEX_NU_LIN_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.nuQuad_TTT(obj)
			val = obj.data.EXT_GIEX_NU_QUAD_TTT;
		end

		function set.nuQuad_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'nuQuad_TTT');
			obj.data.EXT_GIEX_NU_QUAD_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.sigma_TTT(obj)
			val = obj.data.EXT_GIEX_SIGMA_TTT;
		end

		function set.sigma_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'sigma_TTT');
			obj.data.EXT_GIEX_SIGMA_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.lambda_TTT(obj)
			val = obj.data.EXT_GIEX_LAMBDA_TTT;
		end

		function set.lambda_TTT(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'lambda_TTT');
			obj.data.EXT_GIEX_LAMBDA_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.refSolid(obj)
			if isfield(obj.data, 'EXT_GIEX_REFQ')
				val = obj.data.EXT_GIEX_REFQ;
			else
				val = [];
			end
		end

		function set.refSolid(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'EXT_GIEX_REFQ');
			else
				validateattributes(val, {'double'}, {'scalar', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refSolid');
				obj.data.EXT_GIEX_REFQ = val;
			end
			obj.hasChanged = true;
		end

		function val = get.refLiquid(obj)
			if isfield(obj.data, 'EXT_GIEX_REFC0')
				val = obj.data.EXT_GIEX_REFC0;
			else
				val = [];
			end
		end

		function set.refLiquid(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'EXT_GIEX_REFC0');
			else
				validateattributes(val, {'double'}, {'scalar', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refLiquid');
				obj.data.EXT_GIEX_REFC0 = val;
			end
			obj.hasChanged = true;
		end

		function val = get.EXT_GIEX_KA(obj)
			val = obj.kA;
		end
		function set.EXT_GIEX_KA(obj, val)
			obj.kA = val;
		end
		function val = get.EXT_GIEX_KA_LIN(obj)
			val = obj.kALin;
		end
		function set.EXT_GIEX_KA_LIN(obj, val)
			obj.kALin = val;
		end
		function val = get.EXT_GIEX_KA_QUAD(obj)
			val = obj.kAQuad;
		end
		function set.EXT_GIEX_KA_QUAD(obj, val)
			obj.kAQuad = val;
		end
		function val = get.EXT_GIEX_KA_SALT(obj)
			val = obj.kASalt;
		end
		function set.EXT_GIEX_KA_SALT(obj, val)
			obj.kASalt = val;
		end
		function val = get.EXT_GIEX_KA_PROT(obj)
			val = obj.kAProt;
		end
		function set.EXT_GIEX_KA_PROT(obj, val)
			obj.kAProt = val;
		end
		function val = get.EXT_GIEX_KD(obj)
			val = obj.kD;
		end
		function set.EXT_GIEX_KD(obj, val)
			obj.kD = val;
		end
		function val = get.EXT_GIEX_KD_LIN(obj)
			val = obj.kDLin;
		end
		function set.EXT_GIEX_KD_LIN(obj, val)
			obj.kDLin = val;
		end
		function val = get.EXT_GIEX_KD_QUAD(obj)
			val = obj.kDQuad;
		end
		function set.EXT_GIEX_KD_QUAD(obj, val)
			obj.kDQuad = val;
		end
		function val = get.EXT_GIEX_KD_SALT(obj)
			val = obj.kDSalt;
		end
		function set.EXT_GIEX_KD_SALT(obj, val)
			obj.kDSalt = val;
		end
		function val = get.EXT_GIEX_KD_PROT(obj)
			val = obj.kDProt;
		end
		function set.EXT_GIEX_KD_PROT(obj, val)
			obj.kDProt = val;
		end
		function val = get.EXT_GIEX_LAMBDA(obj)
			val = obj.lambda;
		end
		function set.EXT_GIEX_LAMBDA(obj, val)
			obj.lambda = val;
		end
		function val = get.EXT_GIEX_NU(obj)
			val = obj.nu;
		end
		function set.EXT_GIEX_NU(obj, val)
			obj.nu = val;
		end
		function val = get.EXT_GIEX_NU_LIN(obj)
			val = obj.nuLin;
		end
		function set.EXT_GIEX_NU_LIN(obj, val)
			obj.nuLin = val;
		end
		function val = get.EXT_GIEX_NU_QUAD(obj)
			val = obj.nuQuad;
		end
		function set.EXT_GIEX_NU_QUAD(obj, val)
			obj.nuQuad = val;
		end
		function val = get.EXT_GIEX_SIGMA(obj)
			val = obj.sigma;
		end
		function set.EXT_GIEX_SIGMA(obj, val)
			obj.sigma = val;
		end
		function val = get.EXT_GIEX_KA_T(obj)
			val = obj.kA_T;
		end
		function set.EXT_GIEX_KA_T(obj, val)
			obj.kA_T = val;
		end
		function val = get.EXT_GIEX_KA_LIN_T(obj)
			val = obj.kALin_T;
		end
		function set.EXT_GIEX_KA_LIN_T(obj, val)
			obj.kALin_T = val;
		end
		function val = get.EXT_GIEX_KA_QUAD_T(obj)
			val = obj.kAQuad_T;
		end
		function set.EXT_GIEX_KA_QUAD_T(obj, val)
			obj.kAQuad_T = val;
		end
		function val = get.EXT_GIEX_KA_SALT_T(obj)
			val = obj.kASalt_T;
		end
		function set.EXT_GIEX_KA_SALT_T(obj, val)
			obj.kASalt_T = val;
		end
		function val = get.EXT_GIEX_KA_PROT_T(obj)
			val = obj.kAProt_T;
		end
		function set.EXT_GIEX_KA_PROT_T(obj, val)
			obj.kAProt_T = val;
		end
		function val = get.EXT_GIEX_KD_T(obj)
			val = obj.kD_T;
		end
		function set.EXT_GIEX_KD_T(obj, val)
			obj.kD_T = val;
		end
		function val = get.EXT_GIEX_KD_LIN_T(obj)
			val = obj.kDLin_T;
		end
		function set.EXT_GIEX_KD_LIN_T(obj, val)
			obj.kDLin_T = val;
		end
		function val = get.EXT_GIEX_KD_QUAD_T(obj)
			val = obj.kDQuad_T;
		end
		function set.EXT_GIEX_KD_QUAD_T(obj, val)
			obj.kDQuad_T = val;
		end
		function val = get.EXT_GIEX_KD_SALT_T(obj)
			val = obj.kDSalt_T;
		end
		function set.EXT_GIEX_KD_SALT_T(obj, val)
			obj.kDSalt_T = val;
		end
		function val = get.EXT_GIEX_KD_PROT_T(obj)
			val = obj.kDProt_T;
		end
		function set.EXT_GIEX_KD_PROT_T(obj, val)
			obj.kDProt_T = val;
		end
		function val = get.EXT_GIEX_LAMBDA_T(obj)
			val = obj.lambda_T;
		end
		function set.EXT_GIEX_LAMBDA_T(obj, val)
			obj.lambda_T = val;
		end
		function val = get.EXT_GIEX_NU_T(obj)
			val = obj.nu_T;
		end
		function set.EXT_GIEX_NU_T(obj, val)
			obj.nu_T = val;
		end
		function val = get.EXT_GIEX_NU_LIN_T(obj)
			val = obj.nuLin_T;
		end
		function set.EXT_GIEX_NU_LIN_T(obj, val)
			obj.nuLin_T = val;
		end
		function val = get.EXT_GIEX_NU_QUAD_T(obj)
			val = obj.nuQuad_T;
		end
		function set.EXT_GIEX_NU_QUAD_T(obj, val)
			obj.nuQuad_T = val;
		end
		function val = get.EXT_GIEX_SIGMA_T(obj)
			val = obj.sigma_T;
		end
		function set.EXT_GIEX_SIGMA_T(obj, val)
			obj.sigma_T = val;
		end
		function val = get.EXT_GIEX_KA_TT(obj)
			val = obj.kA_TT;
		end
		function set.EXT_GIEX_KA_TT(obj, val)
			obj.kA_TT = val;
		end
		function val = get.EXT_GIEX_KA_LIN_TT(obj)
			val = obj.kALin_TT;
		end
		function set.EXT_GIEX_KA_LIN_TT(obj, val)
			obj.kALin_TT = val;
		end
		function val = get.EXT_GIEX_KA_QUAD_TT(obj)
			val = obj.kAQuad_TT;
		end
		function set.EXT_GIEX_KA_QUAD_TT(obj, val)
			obj.kAQuad_TT = val;
		end
		function val = get.EXT_GIEX_KA_SALT_TT(obj)
			val = obj.kASalt_TT;
		end
		function set.EXT_GIEX_KA_SALT_TT(obj, val)
			obj.kASalt_TT = val;
		end
		function val = get.EXT_GIEX_KA_PROT_TT(obj)
			val = obj.kAProt_TT;
		end
		function set.EXT_GIEX_KA_PROT_TT(obj, val)
			obj.kAProt_TT = val;
		end
		function val = get.EXT_GIEX_KD_TT(obj)
			val = obj.kD_TT;
		end
		function set.EXT_GIEX_KD_TT(obj, val)
			obj.kD_TT = val;
		end
		function val = get.EXT_GIEX_KD_LIN_TT(obj)
			val = obj.kDLin_TT;
		end
		function set.EXT_GIEX_KD_LIN_TT(obj, val)
			obj.kDLin_TT = val;
		end
		function val = get.EXT_GIEX_KD_QUAD_TT(obj)
			val = obj.kDQuad_TT;
		end
		function set.EXT_GIEX_KD_QUAD_TT(obj, val)
			obj.kDQuad_TT = val;
		end
		function val = get.EXT_GIEX_KD_SALT_TT(obj)
			val = obj.kDSalt_TT;
		end
		function set.EXT_GIEX_KD_SALT_TT(obj, val)
			obj.kDSalt_TT = val;
		end
		function val = get.EXT_GIEX_KD_PROT_TT(obj)
			val = obj.kDProt_TT;
		end
		function set.EXT_GIEX_KD_PROT_TT(obj, val)
			obj.kDProt_TT = val;
		end
		function val = get.EXT_GIEX_LAMBDA_TT(obj)
			val = obj.lambda_TT;
		end
		function set.EXT_GIEX_LAMBDA_TT(obj, val)
			obj.lambda_TT = val;
		end
		function val = get.EXT_GIEX_NU_TT(obj)
			val = obj.nu_TT;
		end
		function set.EXT_GIEX_NU_TT(obj, val)
			obj.nu_TT = val;
		end
		function val = get.EXT_GIEX_NU_LIN_TT(obj)
			val = obj.nuLin_TT;
		end
		function set.EXT_GIEX_NU_LIN_TT(obj, val)
			obj.nuLin_TT = val;
		end
		function val = get.EXT_GIEX_NU_QUAD_TT(obj)
			val = obj.nuQuad_TT;
		end
		function set.EXT_GIEX_NU_QUAD_TT(obj, val)
			obj.nuQuad_TT = val;
		end
		function val = get.EXT_GIEX_SIGMA_TT(obj)
			val = obj.sigma_TT;
		end
		function set.EXT_GIEX_SIGMA_TT(obj, val)
			obj.sigma_TT = val;
		end
		function val = get.EXT_GIEX_KA_TTT(obj)
			val = obj.kA_TTT;
		end
		function set.EXT_GIEX_KA_TTT(obj, val)
			obj.kA_TTT = val;
		end
		function val = get.EXT_GIEX_KA_LIN_TTT(obj)
			val = obj.kALin_TTT;
		end
		function set.EXT_GIEX_KA_LIN_TTT(obj, val)
			obj.kALin_TTT = val;
		end
		function val = get.EXT_GIEX_KA_QUAD_TTT(obj)
			val = obj.kAQuad_TTT;
		end
		function set.EXT_GIEX_KA_QUAD_TTT(obj, val)
			obj.kAQuad_TTT = val;
		end
		function val = get.EXT_GIEX_KA_SALT_TTT(obj)
			val = obj.kASalt_TTT;
		end
		function set.EXT_GIEX_KA_SALT_TTT(obj, val)
			obj.kASalt_TTT = val;
		end
		function val = get.EXT_GIEX_KA_PROT_TTT(obj)
			val = obj.kAProt_TTT;
		end
		function set.EXT_GIEX_KA_PROT_TTT(obj, val)
			obj.kAProt_TTT = val;
		end
		function val = get.EXT_GIEX_KD_TTT(obj)
			val = obj.kD_TTT;
		end
		function set.EXT_GIEX_KD_TTT(obj, val)
			obj.kD_TTT = val;
		end
		function val = get.EXT_GIEX_KD_LIN_TTT(obj)
			val = obj.kDLin_TTT;
		end
		function set.EXT_GIEX_KD_LIN_TTT(obj, val)
			obj.kDLin_TTT = val;
		end
		function val = get.EXT_GIEX_KD_QUAD_TTT(obj)
			val = obj.kDQuad_TTT;
		end
		function set.EXT_GIEX_KD_QUAD_TTT(obj, val)
			obj.kDQuad_TTT = val;
		end
		function val = get.EXT_GIEX_KD_SALT_TTT(obj)
			val = obj.kDSalt_TTT;
		end
		function set.EXT_GIEX_KD_SALT_TTT(obj, val)
			obj.kDSalt_TTT = val;
		end
		function val = get.EXT_GIEX_KD_PROT_TTT(obj)
			val = obj.kDProt_TTT;
		end
		function set.EXT_GIEX_KD_PROT_TTT(obj, val)
			obj.kDProt_TTT = val;
		end
		function val = get.EXT_GIEX_LAMBDA_TTT(obj)
			val = obj.lambda_TTT;
		end
		function set.EXT_GIEX_LAMBDA_TTT(obj, val)
			obj.lambda_TTT = val;
		end
		function val = get.EXT_GIEX_NU_TTT(obj)
			val = obj.nu_TTT;
		end
		function set.EXT_GIEX_NU_TTT(obj, val)
			obj.nu_TTT = val;
		end
		function val = get.EXT_GIEX_NU_LIN_TTT(obj)
			val = obj.nuLin_TTT;
		end
		function set.EXT_GIEX_NU_LIN_TTT(obj, val)
			obj.nuLin_TTT = val;
		end
		function val = get.EXT_GIEX_NU_QUAD_TTT(obj)
			val = obj.nuQuad_TTT;
		end
		function set.EXT_GIEX_NU_QUAD_TTT(obj, val)
			obj.nuQuad_TTT = val;
		end
		function val = get.EXT_GIEX_SIGMA_TTT(obj)
			val = obj.sigma_TTT;
		end
		function set.EXT_GIEX_SIGMA_TTT(obj, val)
			obj.sigma_TTT = val;
		end
		function val = get.EXT_GIEX_REFQ(obj)
			val = obj.refSolid;
		end
		function set.EXT_GIEX_REFQ(obj, val)
			obj.refSolid = val;
		end
		function val = get.EXT_GIEX_REFC0(obj)
			val = obj.refLiquid;
		end
		function set.EXT_GIEX_REFC0(obj, val)
			obj.refLiquid = val;
		end

	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = ExtFunGeneralizedIonExchangeBinding();
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
