
classdef ExtFunKumarLangmuirBinding < KineticQuasiStationaryBindingModel
	%ExtFunKumarLangmuirBinding Kumar-Langmuir binding model with external function support
	%
	% See also BINDINGMODEL, KUMARLANGMUIRBINDING, KINETICQUASISTATIONARYBINDINGMODEL
	
	% Copyright: (C) 2008-2021 The CADET Authors
	%            See the license note at the end of the file.

	properties(Constant)
		name = 'EXT_KUMAR_MULTI_COMPONENT_LANGMUIR'; % Name of the binding model according to CADET file format specs
	end

	properties (Dependent, Transient)
		% Adsorption rate in [m^3_MP / (mol * s)], constant term of external dependence
		kA;
		EXT_KMCL_KA;
		% Adsorption rate in [m^3_MP / (mol * s * [T])], linear term of external dependence
		kA_T;
		EXT_KMCL_KA_T;
		% Adsorption rate in [m^3_MP / (mol * s * [T]^2)], quadratic term of external dependence
		kA_TT;
		EXT_KMCL_KA_TT;
		% Adsorption rate in [m^3_MP / (mol * s * [T]^3)], cubic term of external dependence
		kA_TTT;
		EXT_KMCL_KA_TTT;
		% Desorption rate in [m^(3 * nu_i)_MP / (mol^(nu_i) * s)], constant term of external dependence
		kD;
		EXT_KMCL_KD;
		% Desorption rate in [m^(3 * nu_i)_MP / (mol^(nu_i) * s * [T])], linear term of external dependence
		kD_T;
		EXT_KMCL_KD_T;
		% Desorption rate in [m^(3 * nu_i)_MP / (mol^(nu_i) * s * [T]^2)], quadratic term of external dependence
		kD_TT;
		EXT_KMCL_KD_TT;
		% Desorption rate in [m^(3 * nu_i)_MP / (mol^(nu_i) * s * [T]^3)], cubic term of external dependence
		kD_TTT;
		EXT_KMCL_KD_TTT;
		% Activation temperature in [K], constant term of external dependence
		kAct;
		EXT_KMCL_KACT;
		% Activation temperature in [K / [T]], linear term of external dependence
		kAct_T;
		EXT_KMCL_KACT_T;
		% Activation temperature in [K / [T]^2], quadratic term of external dependence
		kAct_TT;
		EXT_KMCL_KACT_TT;
		% Activation temperature in [K / [T]^3], cubic term of external dependence
		kAct_TTT;
		EXT_KMCL_KACT_TTT;
		% Capacity in [mol / m^3_SP], constant term of external dependence
		qMax;
		EXT_KMCL_QMAX;
		% Capacity in [mol / (m^3_SP * [T])], linear term of external dependence
		qMax_T;
		EXT_KMCL_QMAX_T;
		% Capacity in [mol / (m^3_SP * [T]^2)], quadratic term of external dependence
		qMax_TT;
		EXT_KMCL_QMAX_TT;
		% Capacity in [mol / (m^3_SP * [T]^3)], cubic term of external dependence
		qMax_TTT;
		EXT_KMCL_QMAX_TTT;
		% Characteristic charge [-], constant term of external dependence
		nu;
		EXT_KMCL_NU;
		% Characteristic charge [1 / [T]], linear term of external dependence
		nu_T;
		EXT_KMCL_NU_T;
		% Characteristic charge [1 / [T]^2], quadratic term of external dependence
		nu_TT;
		EXT_KMCL_NU_TT;
		% Characteristic charge [1 / [T]^3], cubic term of external dependence
		nu_TTT;
		EXT_KMCL_NU_TTT;
		% Temperature in [K], constant term of external dependence
		temperature;
		EXT_KMCL_TEMP;
		% Temperature in [K / [T]], linear term of external dependence
		temperature_T;
		EXT_KMCL_TEMP_T;
		% Temperature in [K / [T]^2], quadratic term of external dependence
		temperature_TT;
		EXT_KMCL_TEMP_TT;
		% Temperature in [K / [T]^3], cubic term of external dependence
		temperature_TTT;
		EXT_KMCL_TEMP_TTT;
		% Indices of external functions (0-based)
		externalSource;
		EXTFUN;
	end

	methods

		function obj = ExtFunKumarLangmuirBinding(kA, kD, qMax)
			%EXTFUNKUMARLANGMUIRBINDING Constructs an ExtFunKumarLangmuirBinding object with external function support
			%   OBJ = EXTFUNKUMARLANGMUIRBINDING(KA) creates an ExtFunKumarLangmuirBinding
			%   model with the given adsorption rates KA.
			%
			%   OBJ = EXTFUNKUMARLANGMUIRBINDING(..., KD) also sets the desorption rates to KD.
			%
			%   OBJ = EXTFUNKUMARLANGMUIRBINDING(..., KD, QMAX) also sets the desorption rates
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

			validateattributes(obj.kA, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kA');
			validateattributes(obj.kD, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kD');
			validateattributes(obj.kAct, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kAct');
			validateattributes(obj.qMax, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'qMax');
			validateattributes(obj.nu, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'nu');
			validateattributes(obj.temperature, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'temperature');

			validateattributes(obj.kA_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kA_T');
			validateattributes(obj.kD_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kD_T');
			validateattributes(obj.kAct_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kAct_T');
			validateattributes(obj.qMax_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'qMax_T');
			validateattributes(obj.nu_T, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'nu_T');
			validateattributes(obj.temperature_T, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'temperature_T');

			validateattributes(obj.kA_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kA_TT');
			validateattributes(obj.kD_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kD_TT');
			validateattributes(obj.kAct_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kAct_TT');
			validateattributes(obj.qMax_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'qMax_TT');
			validateattributes(obj.nu_TT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'nu_TT');
			validateattributes(obj.temperature_TT, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'temperature_TT');

			validateattributes(obj.kA_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kA_TTT');
			validateattributes(obj.kD_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kD_TTT');
			validateattributes(obj.kAct_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kAct_TTT');
			validateattributes(obj.qMax_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'qMax_TTT');
			validateattributes(obj.nu_TTT, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'nu_TTT');
			validateattributes(obj.temperature_TTT, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'temperature_TTT');

			if ~isempty(obj.externalSource)
				validateattributes(obj.externalSource, {'double', 'int32'}, {'vector', 'nonempty', 'finite', 'real', '>=', 0}, '', 'externalSource');
				if ~any(numel(obj.externalSource) == [1, 6])
					error('CADET:invalidConfig', 'Expected externalSource to have either 1 or 6 entries but got %d entries.', numel(obj.externalSource));
				end
			end
			res = obj.validate@KineticQuasiStationaryBindingModel(nComponents, nBoundStates);
		end

		function val = get.kA(obj)
			val = obj.data.EXT_KMCL_KA;
		end

		function set.kA(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kA');
			obj.data.EXT_KMCL_KA = val;
			obj.hasChanged = true;
		end

		function val = get.kD(obj)
			val = obj.data.EXT_KMCL_KD;
		end

		function set.kD(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kD');
			obj.data.EXT_KMCL_KD = val;
			obj.hasChanged = true;
		end

		function val = get.kAct(obj)
			val = obj.data.EXT_KMCL_KACT;
		end

		function set.kAct(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kAct');
			obj.data.EXT_KMCL_KACT = val;
			obj.hasChanged = true;
		end

		function val = get.qMax(obj)
			val = obj.data.EXT_KMCL_QMAX;
		end

		function set.qMax(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'qMax');
			obj.data.EXT_KMCL_QMAX = val;
			obj.hasChanged = true;
		end

		function val = get.nu(obj)
			val = obj.data.EXT_KMCL_NU;
		end

		function set.nu(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'nu');
			obj.data.EXT_KMCL_NU = val;
			obj.hasChanged = true;
		end

		function val = get.temperature(obj)
			val = obj.data.EXT_KMCL_TEMP;
		end

		function set.temperature(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'temperature');
			obj.data.EXT_KMCL_TEMP = val;
			obj.hasChanged = true;
		end

		function val = get.kA_T(obj)
			val = obj.data.EXT_KMCL_KA_T;
		end

		function set.kA_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kA_T');
			obj.data.EXT_KMCL_KA_T = val;
			obj.hasChanged = true;
		end

		function val = get.kD_T(obj)
			val = obj.data.EXT_KMCL_KD_T;
		end

		function set.kD_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kD_T');
			obj.data.EXT_KMCL_KD_T = val;
			obj.hasChanged = true;
		end

		function val = get.kAct_T(obj)
			val = obj.data.EXT_KMCL_KACT_T;
		end

		function set.kAct_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kAct_T');
			obj.data.EXT_KMCL_KACT_T = val;
			obj.hasChanged = true;
		end

		function val = get.qMax_T(obj)
			val = obj.data.EXT_KMCL_QMAX_T;
		end

		function set.qMax_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'qMax_T');
			obj.data.EXT_KMCL_QMAX_T = val;
			obj.hasChanged = true;
		end

		function val = get.nu_T(obj)
			val = obj.data.EXT_KMCL_NU_T;
		end

		function set.nu_T(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'nu_T');
			obj.data.EXT_KMCL_NU_T = val;
			obj.hasChanged = true;
		end

		function val = get.temperature_T(obj)
			val = obj.data.EXT_KMCL_TEMP_T;
		end

		function set.temperature_T(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'temperature_T');
			obj.data.EXT_KMCL_TEMP_T = val;
			obj.hasChanged = true;
		end

		function val = get.kA_TT(obj)
			val = obj.data.EXT_KMCL_KA_TT;
		end

		function set.kA_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kA_TT');
			obj.data.EXT_KMCL_KA_TT = val;
			obj.hasChanged = true;
		end

		function val = get.kD_TT(obj)
			val = obj.data.EXT_KMCL_KD_TT;
		end

		function set.kD_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kD_TT');
			obj.data.EXT_KMCL_KD_TT = val;
			obj.hasChanged = true;
		end

		function val = get.kAct_TT(obj)
			val = obj.data.EXT_KMCL_KACT_TT;
		end

		function set.kAct_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kAct_TT');
			obj.data.EXT_KMCL_KACT_TT = val;
			obj.hasChanged = true;
		end

		function val = get.qMax_TT(obj)
			val = obj.data.EXT_KMCL_QMAX_TT;
		end

		function set.qMax_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'qMax_TT');
			obj.data.EXT_KMCL_QMAX_TT = val;
			obj.hasChanged = true;
		end

		function val = get.nu_TT(obj)
			val = obj.data.EXT_KMCL_NU_TT;
		end

		function set.nu_TT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'nu_TT');
			obj.data.EXT_KMCL_NU_TT = val;
			obj.hasChanged = true;
		end

		function val = get.temperature_TT(obj)
			val = obj.data.EXT_KMCL_TEMP_TT;
		end

		function set.temperature_TT(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'temperature_TT');
			obj.data.EXT_KMCL_TEMP_TT = val;
			obj.hasChanged = true;
		end

		function val = get.kA_TTT(obj)
			val = obj.data.EXT_KMCL_KA_TTT;
		end

		function set.kA_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kA_TTT');
			obj.data.EXT_KMCL_KA_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.kD_TTT(obj)
			val = obj.data.EXT_KMCL_KD_TTT;
		end

		function set.kD_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kD_TTT');
			obj.data.EXT_KMCL_KD_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.kAct_TTT(obj)
			val = obj.data.EXT_KMCL_KACT_TTT;
		end

		function set.kAct_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kAct_TTT');
			obj.data.EXT_KMCL_KACT_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.qMax_TTT(obj)
			val = obj.data.EXT_KMCL_QMAX_TTT;
		end

		function set.qMax_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'qMax_TTT');
			obj.data.EXT_KMCL_QMAX_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.nu_TTT(obj)
			val = obj.data.EXT_KMCL_NU_TTT;
		end

		function set.nu_TTT(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'nu_TTT');
			obj.data.EXT_KMCL_NU_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.temperature_TTT(obj)
			val = obj.data.EXT_KMCL_TEMP_TTT;
		end

		function set.temperature_TTT(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'temperature_TTT');
			obj.data.EXT_KMCL_TEMP_TTT = val;
			obj.hasChanged = true;
		end

		function val = get.EXT_KMCL_KA(obj)
			val = obj.kA;
		end
		function set.EXT_KMCL_KA(obj, val)
			obj.kA = val;
		end
		function val = get.EXT_KMCL_KD(obj)
			val = obj.kD;
		end
		function set.EXT_KMCL_KD(obj, val)
			obj.kD = val;
		end
		function val = get.EXT_KMCL_KACT(obj)
			val = obj.kAct;
		end
		function set.EXT_KMCL_KACT(obj, val)
			obj.kAct = val;
		end
		function val = get.EXT_KMCL_QMAX(obj)
			val = obj.qMax;
		end
		function set.EXT_KMCL_QMAX(obj, val)
			obj.qMax = val;
		end
		function val = get.EXT_KMCL_NU(obj)
			val = obj.nu;
		end
		function set.EXT_KMCL_NU(obj, val)
			obj.nu = val;
		end
		function val = get.EXT_KMCL_TEMP(obj)
			val = obj.temperature;
		end
		function set.EXT_KMCL_TEMP(obj, val)
			obj.temperature = val;
		end

		function val = get.EXT_KMCL_KA_T(obj)
			val = obj.kA_T;
		end
		function set.EXT_KMCL_KA_T(obj, val)
			obj.kA_T = val;
		end
		function val = get.EXT_KMCL_KD_T(obj)
			val = obj.kD_T;
		end
		function set.EXT_KMCL_KD_T(obj, val)
			obj.kD_T = val;
		end
		function val = get.EXT_KMCL_KACT_T(obj)
			val = obj.kAct_T;
		end
		function set.EXT_KMCL_KACT_T(obj, val)
			obj.kAct_T = val;
		end
		function val = get.EXT_KMCL_QMAX_T(obj)
			val = obj.qMax_T;
		end
		function set.EXT_KMCL_QMAX_T(obj, val)
			obj.qMax_T = val;
		end
		function val = get.EXT_KMCL_NU_T(obj)
			val = obj.nu_T;
		end
		function set.EXT_KMCL_NU_T(obj, val)
			obj.nu_T = val;
		end
		function val = get.EXT_KMCL_TEMP_T(obj)
			val = obj.temperature_T;
		end
		function set.EXT_KMCL_TEMP_T(obj, val)
			obj.temperature_T = val;
		end

		function val = get.EXT_KMCL_KA_TT(obj)
			val = obj.kA_TT;
		end
		function set.EXT_KMCL_KA_TT(obj, val)
			obj.kA_TT = val;
		end
		function val = get.EXT_KMCL_KD_TT(obj)
			val = obj.kD_TT;
		end
		function set.EXT_KMCL_KD_TT(obj, val)
			obj.kD_TT = val;
		end
		function val = get.EXT_KMCL_KACT_TT(obj)
			val = obj.kAct_TT;
		end
		function set.EXT_KMCL_KACT_TT(obj, val)
			obj.kAct_TT = val;
		end
		function val = get.EXT_KMCL_QMAX_TT(obj)
			val = obj.qMax_TT;
		end
		function set.EXT_KMCL_QMAX_TT(obj, val)
			obj.qMax_TT = val;
		end
		function val = get.EXT_KMCL_NU_TT(obj)
			val = obj.nu_TT;
		end
		function set.EXT_KMCL_NU_TT(obj, val)
			obj.nu_TT = val;
		end
		function val = get.EXT_KMCL_TEMP_TT(obj)
			val = obj.temperature_TT;
		end
		function set.EXT_KMCL_TEMP_TT(obj, val)
			obj.temperature_TT = val;
		end

		function val = get.EXT_KMCL_KA_TTT(obj)
			val = obj.kA_TTT;
		end
		function set.EXT_KMCL_KA_TTT(obj, val)
			obj.kA_TTT = val;
		end
		function val = get.EXT_KMCL_KD_TTT(obj)
			val = obj.kD_TTT;
		end
		function set.EXT_KMCL_KD_TTT(obj, val)
			obj.kD_TTT = val;
		end
		function val = get.EXT_KMCL_KACT_TTT(obj)
			val = obj.kAct_TTT;
		end
		function set.EXT_KMCL_KACT_TTT(obj, val)
			obj.kAct_TTT = val;
		end
		function val = get.EXT_KMCL_QMAX_TTT(obj)
			val = obj.qMax_TTT;
		end
		function set.EXT_KMCL_QMAX_TTT(obj, val)
			obj.qMax_TTT = val;
		end
		function val = get.EXT_KMCL_NU_TTT(obj)
			val = obj.nu_TTT;
		end
		function set.EXT_KMCL_NU_TTT(obj, val)
			obj.nu_TTT = val;
		end
		function val = get.EXT_KMCL_TEMP_TTT(obj)
			val = obj.temperature_TTT;
		end
		function set.EXT_KMCL_TEMP_TTT(obj, val)
			obj.temperature_TTT = val;
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
				obj = ExtFunKumarLangmuirBinding();
				obj.loadobjInternal(S);
			end
		end
		
	end

end

% =============================================================================
%  CADET
%  
%  Copyright (C) 2008-2021: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. obj program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies obj distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
