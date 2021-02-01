
classdef GeneralizedIonExchangeBinding < KineticQuasiStationaryBindingModel
	%GeneralizedIonExchangeBinding Generalized ion exchange binding model
	%
	% See also BINDINGMODEL, KINETICQUASISTATIONARYBINDINGMODEL

	% Copyright: (C) 2008-2020 The CADET Authors
	%            See the license note at the end of the file.

	properties(Constant)
		name = 'GENERALIZED_ION_EXCHANGE'; % Name of the binding model according to CADET file format specs
	end

	properties (Dependent, Transient)
		% Adsorption rate in [m^3_MP / (m^3_SP * s)]
		kA;
		GIEX_KA;
		% Adsorption rate linear modifier coefficient in [1 / Mod]
		kALin;
		GIEX_KA_LIN;
		% Adsorption rate quadratic modifier coefficient in [1 / Mod^2]
		kAQuad;
		GIEX_KA_QUAD;
		% Adsorption rate protein-protein modifier coefficient in [m^3_MP / mol]
		kAProt;
		GIEX_KA_PROT;
		% Adsorption rate salt modifier coefficient [-]
		kASalt;
		GIEX_KA_SALT;
		% Desorption rate in [1 / s]
		kD;
		GIEX_KD;
		% Desorption rate linear modifier coefficient in [1 / Mod]
		kDLin;
		GIEX_KD_LIN;
		% Desorption rate quadratic modifier coefficient in [1 / Mod^2]
		kDQuad;
		GIEX_KD_QUAD;
		% Desorption rate protein-protein modifier coefficient in [m^3_MP / mol]
		kDProt;
		GIEX_KD_PROT;
		% Desorption rate salt modifier coefficient [-]
		kDSalt;
		GIEX_KD_SALT;
		% Characteristic charge [-]
		nu;
		GIEX_NU;
		% Linear dependence on modifier of characteristic charge [1 / Mod]
		nuLin;
		GIEX_NU_LIN;
		% Quadratic dependence on modifier of characteristic charge [1 / Mod^2]
		nuQuad;
		GIEX_NU_QUAD;
		% Steric factor [-]
		sigma;
		GIEX_SIGMA;
		% Ionic capacity in [mol / m^3_SP]
		lambda;
		GIEX_LAMBDA;
		% Reference solid phase concentration [mol / m^3_SP]
		refSolid;
		GIEX_REFQ;
		% Reference liquid phase concentration [mol / m^3_MP]
		refLiquid;
		GIEX_REFC0;
	end

	methods

		function obj = GeneralizedIonExchangeBinding(kA, kD)
			%GENERALIZEDIONEXCHANGEBINDING Constructs a GeneralizedIonExchangeBinding object
			%   OBJ = GENERALIZEDIONEXCHANGEBINDING(KA) creates a GeneralizedIonExchangeBinding model with
			%   the given adsorption rates KA.
			%
			%   OBJ = GENERALIZEDIONEXCHANGEBINDING(..., KD) also sets the desorption rates to KD.

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
				error('CADET:invalidConfig', 'GeneralizedIonExchangeBinding requires first component (salt) to have one bound state');
			end

			if nBoundStates(2) ~= 0
				error('CADET:invalidConfig', 'GeneralizedIonExchangeBinding requires second component (modifier) to be non-binding');
			end

			validateattributes(obj.kA, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'kA');
			validateattributes(obj.kALin, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kALin');
			validateattributes(obj.kAQuad, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kAQuad');
			validateattributes(obj.kAProt, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kAProt');
			validateattributes(obj.kASalt, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kASalt');
			validateattributes(obj.kD, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'kD');
			validateattributes(obj.kDLin, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kDLin');
			validateattributes(obj.kDQuad, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kDQuad');
			validateattributes(obj.kDProt, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kDProt');
			validateattributes(obj.kDSalt, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'kDSalt');
			validateattributes(obj.nu, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'nu');
			validateattributes(obj.nuLin, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'nuLin');
			validateattributes(obj.nuQuad, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'nuQuad');
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
			val = obj.data.GIEX_KA;
		end

		function set.kA(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'kA');
			obj.data.GIEX_KA = val;
			obj.hasChanged = true;
		end

		function val = get.kALin(obj)
			val = obj.data.GIEX_KA_LIN;
		end

		function set.kALin(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kALin');
			obj.data.GIEX_KA_LIN = val;
			obj.hasChanged = true;
		end

		function val = get.kAQuad(obj)
			val = obj.data.GIEX_KA_QUAD;
		end

		function set.kAQuad(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kAQuad');
			obj.data.GIEX_KA_QUAD = val;
			obj.hasChanged = true;
		end

		function val = get.kASalt(obj)
			val = obj.data.GIEX_KA_SALT;
		end

		function set.kASalt(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kASalt');
			obj.data.GIEX_KA_SALT = val;
			obj.hasChanged = true;
		end

		function val = get.kAProt(obj)
			val = obj.data.GIEX_KA_PROT;
		end

		function set.kAProt(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kAProt');
			obj.data.GIEX_KA_PROT = val;
			obj.hasChanged = true;
		end

		function val = get.kD(obj)
			val = obj.data.GIEX_KD;
		end

		function set.kD(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'kD');
			obj.data.GIEX_KD = val;
			obj.hasChanged = true;
		end

		function val = get.kDLin(obj)
			val = obj.data.GIEX_KD_LIN;
		end

		function set.kDLin(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kDLin');
			obj.data.GIEX_KD_LIN = val;
			obj.hasChanged = true;
		end

		function val = get.kDQuad(obj)
			val = obj.data.GIEX_KD_QUAD;
		end

		function set.kDQuad(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kDQuad');
			obj.data.GIEX_KD_QUAD = val;
			obj.hasChanged = true;
		end

		function val = get.kDSalt(obj)
			val = obj.data.GIEX_KD_SALT;
		end

		function set.kDSalt(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kDSalt');
			obj.data.GIEX_KD_SALT = val;
			obj.hasChanged = true;
		end

		function val = get.kDProt(obj)
			val = obj.data.GIEX_KD_PROT;
		end

		function set.kDProt(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'kDProt');
			obj.data.GIEX_KD_PROT = val;
			obj.hasChanged = true;
		end

		function val = get.nu(obj)
			val = obj.data.GIEX_NU;
		end

		function set.nu(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'nu');
			obj.data.GIEX_NU = val;
			obj.hasChanged = true;
		end

		function val = get.nuLin(obj)
			val = obj.data.GIEX_NU_LIN;
		end

		function set.nuLin(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'nuLin');
			obj.data.GIEX_NU_LIN = val;
			obj.hasChanged = true;
		end

		function val = get.nuQuad(obj)
			val = obj.data.GIEX_NU_QUAD;
		end

		function set.nuQuad(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'nuQuad');
			obj.data.GIEX_NU_QUAD = val;
			obj.hasChanged = true;
		end

		function val = get.sigma(obj)
			val = obj.data.GIEX_SIGMA;
		end

		function set.sigma(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'sigma');
			obj.data.GIEX_SIGMA = val;
			obj.hasChanged = true;
		end

		function val = get.lambda(obj)
			val = obj.data.GIEX_LAMBDA;
		end

		function set.lambda(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'lambda');
			obj.data.GIEX_LAMBDA = val;
			obj.hasChanged = true;
		end

		function val = get.refSolid(obj)
			if isfield(obj.data, 'GIEX_REFQ')
				val = obj.data.GIEX_REFQ;
			else
				val = [];
			end
		end

		function set.refSolid(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'GIEX_REFQ');
			else
				validateattributes(val, {'double'}, {'scalar', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refSolid');
				obj.data.GIEX_REFQ = val;
			end
			obj.hasChanged = true;
		end

		function val = get.refLiquid(obj)
			if isfield(obj.data, 'GIEX_REFC0')
				val = obj.data.GIEX_REFC0;
			else
				val = [];
			end
		end

		function set.refLiquid(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'GIEX_REFC0');
			else
				validateattributes(val, {'double'}, {'scalar', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refLiquid');
				obj.data.GIEX_REFC0 = val;
			end
			obj.hasChanged = true;
		end

		function val = get.GIEX_KA(obj)
			val = obj.kA;
		end
		function set.GIEX_KA(obj, val)
			obj.kA = val;
		end
		function val = get.GIEX_KA_LIN(obj)
			val = obj.kALin;
		end
		function set.GIEX_KA_LIN(obj, val)
			obj.kALin = val;
		end
		function val = get.GIEX_KA_QUAD(obj)
			val = obj.kAQuad;
		end
		function set.GIEX_KA_QUAD(obj, val)
			obj.kAQuad = val;
		end
		function val = get.GIEX_KA_SALT(obj)
			val = obj.kASalt;
		end
		function set.GIEX_KA_SALT(obj, val)
			obj.kASalt = val;
		end
		function val = get.GIEX_KA_PROT(obj)
			val = obj.kAProt;
		end
		function set.GIEX_KA_PROT(obj, val)
			obj.kAProt = val;
		end
		function val = get.GIEX_KD(obj)
			val = obj.kD;
		end
		function set.GIEX_KD(obj, val)
			obj.kD = val;
		end
		function val = get.GIEX_KD_LIN(obj)
			val = obj.kDLin;
		end
		function set.GIEX_KD_LIN(obj, val)
			obj.kDLin = val;
		end
		function val = get.GIEX_KD_QUAD(obj)
			val = obj.kDQuad;
		end
		function set.GIEX_KD_QUAD(obj, val)
			obj.kDQuad = val;
		end
		function val = get.GIEX_KD_SALT(obj)
			val = obj.kDSalt;
		end
		function set.GIEX_KD_SALT(obj, val)
			obj.kDSalt = val;
		end
		function val = get.GIEX_KD_PROT(obj)
			val = obj.kDProt;
		end
		function set.GIEX_KD_PROT(obj, val)
			obj.kDProt = val;
		end
		function val = get.GIEX_LAMBDA(obj)
			val = obj.lambda;
		end
		function set.GIEX_LAMBDA(obj, val)
			obj.lambda = val;
		end
		function val = get.GIEX_NU(obj)
			val = obj.nu;
		end
		function set.GIEX_NU(obj, val)
			obj.nu = val;
		end
		function val = get.GIEX_NU_LIN(obj)
			val = obj.nuLin;
		end
		function set.GIEX_NU_LIN(obj, val)
			obj.nuLin = val;
		end
		function val = get.GIEX_NU_QUAD(obj)
			val = obj.nuQuad;
		end
		function set.GIEX_NU_QUAD(obj, val)
			obj.nuQuad = val;
		end
		function val = get.GIEX_SIGMA(obj)
			val = obj.sigma;
		end
		function set.GIEX_SIGMA(obj, val)
			obj.sigma = val;
		end
		function val = get.GIEX_REFQ(obj)
			val = obj.refSolid;
		end
		function set.GIEX_REFQ(obj, val)
			obj.refSolid = val;
		end
		function val = get.GIEX_REFC0(obj)
			val = obj.refLiquid;
		end
		function set.GIEX_REFC0(obj, val)
			obj.refLiquid = val;
		end

	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = GeneralizedIonExchangeBinding();
				obj.loadobjInternal(S);
			end
		end
		
	end
end

% =============================================================================
%  CADET
%  
%  Copyright (C) 2008-2020: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. obj program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies obj distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
