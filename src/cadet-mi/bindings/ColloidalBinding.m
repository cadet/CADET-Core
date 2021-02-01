
classdef ColloidalBinding < KineticQuasiStationaryBindingModel
	%ColloidalBinding Colloidal binding model
	%   A multi component (i.e., competitive) colloidal binding model.
	%
	% See also BINDINGMODEL, KINETICQUASISTATIONARYBINDINGMODEL

	% Copyright: (C) 2008-2021 The CADET Authors
	%            See the license note at the end of the file.

	properties(Constant)
		name = 'MULTI_COMPONENT_COLLOIDAL'; % Name of the binding model according to CADET file format specs
	end

	properties (Dependent, Transient)
		phi;
		COL_PHI;
		kappaExponent;
		COL_KAPPA_EXP;
		kappaFactor;
		COL_KAPPA_FACT;
		kappaConstant;
		COL_KAPPA_CONST;
		cordNum;
		COL_CORDNUM;
		logKeqPhExponent;
		COL_LOGKEQ_PH_EXP;
		logKeqSaltPowerLawExponent;
		COL_LOGKEQ_SALT_POWEXP;
		logKeqSaltPowerLawFactor;
		COL_LOGKEQ_SALT_POWFACT;
		logKeqSaltExpLawFactor;
		COL_LOGKEQ_SALT_EXPFACT;
		logKeqSaltExpLawExponentFactor;
		COL_LOGKEQ_SALT_EXPARGMULT;
		bppPhExponent;
		COL_BPP_PH_EXP;
		bppSaltPowerLawExponent;
		COL_BPP_SALT_POWEXP;
		bppSaltPowerLawFactor;
		COL_BPP_SALT_POWFACT;
		bppSaltExpLawFactor;
		COL_BPP_SALT_EXPFACT;
		bppSaltExpLawExponentFactor;
		COL_BPP_SALT_EXPARGMULT;
		radius;
		COL_RADIUS;
		kKin;
		COL_KKIN;
		linearizationThreshold;
		COL_LINEAR_THRESHOLD;
		usePh;
		COL_USE_PH;
	end

	methods

		function obj = ColloidalBinding()
			%COLLOIDALBINDING Constructs a ColloidalBinding object

			obj = obj@KineticQuasiStationaryBindingModel();
			obj.linearizationThreshold = 1e-8;
		end

		function res = validate(obj, nComponents, nBoundStates)
			%VALIDATE Validates the binding model parameters
			%   RES = VALIDATE(NCOMPONENTS, NBOUNDSTATES) uses the number of
			%   components NCOMPONENTS of the model and the number of bound
			%   states of each of these components as given by the vector
			%   NBOUNDSTATES to validate the parameters of the binding model.
			%   Returns true in RES if everything is fine and false otherwise.

			if (nBoundStates(1) ~= 0)
				error('CADET:invalidConfig', 'Salt component has to be non-binding (component 0).');
			end

			validateattributes(obj.usePh, {'logical'}, {'scalar', 'nonempty'}, '', 'usePh');
			if (obj.usePh)
				if (nComponents <= 2)
					error('CADET:invalidConfig', 'Expected at least 3 components (at least 1 protein).');
				end
				if (nBoundStates(2) ~= 0)
					error('CADET:invalidConfig', 'PH pseudo component has to be non-binding (component 2).');
				end
			end

			validateattributes(obj.linearizationThreshold, {'double'}, {'scalar', 'nonempty', 'finite', 'real', '>', 0.0}, '', 'linearizationThreshold');
			validateattributes(obj.phi, {'double'}, {'scalar', 'nonempty', 'finite', 'real', '>=', 0.0}, '', 'phi');
			validateattributes(obj.kappaExponent, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'kappaExponent');
			validateattributes(obj.kappaFactor, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'kappaFactor');
			validateattributes(obj.kappaConstant, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'kappaConstant');
			validateattributes(obj.cordNum, {'double'}, {'scalar', 'nonempty', 'finite', 'real', '>', 0.0}, '', 'cordNum');

			validateattributes(obj.logKeqPhExponent, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'logKeqPhExponent');
			validateattributes(obj.logKeqSaltPowerLawExponent, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'logKeqSaltPowerLawExponent');
			validateattributes(obj.logKeqSaltPowerLawFactor, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'logKeqSaltPowerLawFactor');
			validateattributes(obj.logKeqSaltExpLawFactor, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'logKeqSaltExpLawFactor');
			validateattributes(obj.logKeqSaltExpLawExponentFactor, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'logKeqSaltExpLawExponentFactor');
			validateattributes(obj.bppPhExponent, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'bppPhExponent');
			validateattributes(obj.bppSaltPowerLawExponent, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'bppSaltPowerLawExponent');
			validateattributes(obj.bppSaltPowerLawFactor, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'bppSaltPowerLawFactor');
			validateattributes(obj.bppSaltExpLawFactor, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'bppSaltExpLawFactor');
			validateattributes(obj.bppSaltExpLawExponentFactor, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'bppSaltExpLawExponentFactor');
			validateattributes(obj.radius, {'double'}, {'vector', 'nonempty', 'finite', 'real', '>=', 0.0, 'numel', nComponents}, '', 'radius');
			validateattributes(obj.kKin, {'double'}, {'vector', 'nonempty', 'finite', 'real', '>=', 0.0, 'numel', nComponents}, '', 'kKin');

			res = obj.validate@KineticQuasiStationaryBindingModel(nComponents, nBoundStates);
		end

		function val = get.usePh(obj)
			val = logical(obj.data.COL_USE_PH);
		end

		function set.usePh(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'usePh');
			obj.data.COL_USE_PH = int32(logical(val));
			obj.hasChanged = true;
		end

		function val = get.phi(obj)
			val = obj.data.COL_PHI;
		end
		function set.phi(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real', '>=', 0.0}, '', 'phi');
			obj.data.COL_PHI = val;
			obj.hasChanged = true;
		end

		function val = get.kappaExponent(obj)
			val = obj.data.COL_KAPPA_EXP;
		end
		function set.kappaExponent(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'kappaExponent');
			obj.data.COL_KAPPA_EXP = val;
			obj.hasChanged = true;
		end

		function val = get.kappaFactor(obj)
			val = obj.data.COL_KAPPA_FACT;
		end
		function set.kappaFactor(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'kappaFactor');
			obj.data.COL_KAPPA_FACT = val;
			obj.hasChanged = true;
		end

		function val = get.kappaConstant(obj)
			val = obj.data.COL_KAPPA_CONST;
		end
		function set.kappaConstant(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real'}, '', 'kappaConstant');
			obj.data.COL_KAPPA_CONST = val;
			obj.hasChanged = true;
		end

		function val = get.cordNum(obj)
			val = obj.data.COL_CORDNUM;
		end
		function set.cordNum(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real', '>', 0.0}, '', 'cordNum');
			obj.data.COL_CORDNUM = val;
			obj.hasChanged = true;
		end

		function val = get.logKeqPhExponent(obj)
			val = obj.data.COL_LOGKEQ_PH_EXP;
		end
		function set.logKeqPhExponent(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'logKeqPhExponent');
			obj.data.COL_LOGKEQ_PH_EXP = val;
			obj.hasChanged = true;
		end

		function val = get.logKeqSaltPowerLawExponent(obj)
			val = obj.data.COL_LOGKEQ_SALT_POWEXP;
		end
		function set.logKeqSaltPowerLawExponent(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'logKeqSaltPowerLawExponent');
			obj.data.COL_LOGKEQ_SALT_POWEXP = val;
			obj.hasChanged = true;
		end

		function val = get.logKeqSaltPowerLawFactor(obj)
			val = obj.data.COL_LOGKEQ_SALT_POWFACT;
		end
		function set.logKeqSaltPowerLawFactor(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'logKeqSaltPowerLawFactor');
			obj.data.COL_LOGKEQ_SALT_POWFACT = val;
			obj.hasChanged = true;
		end

		function val = get.logKeqSaltExpLawFactor(obj)
			val = obj.data.COL_LOGKEQ_SALT_EXPFACT;
		end
		function set.logKeqSaltExpLawFactor(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'logKeqSaltExpLawFactor');
			obj.data.COL_LOGKEQ_SALT_EXPFACT = val;
			obj.hasChanged = true;
		end

		function val = get.logKeqSaltExpLawExponentFactor(obj)
			val = obj.data.COL_LOGKEQ_SALT_EXPARGMULT;
		end
		function set.logKeqSaltExpLawExponentFactor(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'logKeqSaltExpLawExponentFactor');
			obj.data.COL_LOGKEQ_SALT_EXPARGMULT = val;
			obj.hasChanged = true;
		end

		function val = get.bppPhExponent(obj)
			val = obj.data.COL_BPP_PH_EXP;
		end
		function set.bppPhExponent(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'bppPhExponent');
			obj.data.COL_BPP_PH_EXP = val;
			obj.hasChanged = true;
		end

		function val = get.bppSaltPowerLawExponent(obj)
			val = obj.data.COL_BPP_SALT_POWEXP;
		end
		function set.bppSaltPowerLawExponent(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'bppSaltPowerLawExponent');
			obj.data.COL_BPP_SALT_POWEXP = val;
			obj.hasChanged = true;
		end

		function val = get.bppSaltPowerLawFactor(obj)
			val = obj.data.COL_BPP_SALT_POWFACT;
		end
		function set.bppSaltPowerLawFactor(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'bppSaltPowerLawFactor');
			obj.data.COL_BPP_SALT_POWFACT = val;
			obj.hasChanged = true;
		end

		function val = get.bppSaltExpLawFactor(obj)
			val = obj.data.COL_BPP_SALT_EXPFACT;
		end
		function set.bppSaltExpLawFactor(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'bppSaltExpLawFactor');
			obj.data.COL_BPP_SALT_EXPFACT = val;
			obj.hasChanged = true;
		end

		function val = get.bppSaltExpLawExponentFactor(obj)
			val = obj.data.COL_BPP_SALT_EXPARGMULT;
		end
		function set.bppSaltExpLawExponentFactor(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'bppSaltExpLawExponentFactor');
			obj.data.COL_BPP_SALT_EXPARGMULT = val;
			obj.hasChanged = true;
		end

		function val = get.radius(obj)
			val = obj.data.COL_RADIUS;
		end
		function set.radius(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real', '>=', 0.0}, '', 'radius');
			obj.data.COL_RADIUS = val;
			obj.hasChanged = true;
		end

		function val = get.kKin(obj)
			val = obj.data.COL_KKIN;
		end
		function set.kKin(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real', '>=', 0.0}, '', 'kKin');
			obj.data.COL_KKIN = val;
			obj.hasChanged = true;
		end

		function val = get.linearizationThreshold(obj)
			val = obj.data.COL_LINEAR_THRESHOLD;
		end
		function set.linearizationThreshold(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real', '>', 0.0}, '', 'linearizationThreshold');
			obj.data.COL_LINEAR_THRESHOLD = val;
			obj.hasChanged = true;
		end

		function val = get.COL_PHI(obj)
			val = obj.phi;
		end
		function set.COL_PHI(obj, val)
			obj.phi = val;
		end
		function val = get.COL_KAPPA_EXP(obj)
			val = obj.kappaExponent;
		end
		function set.COL_KAPPA_EXP(obj, val)
			obj.kappaExponent = val;
		end
		function val = get.COL_KAPPA_FACT(obj)
			val = obj.kappaFactor;
		end
		function set.COL_KAPPA_FACT(obj, val)
			obj.kappaFactor = val;
		end
		function val = get.COL_KAPPA_CONST(obj)
			val = obj.kappaConstant;
		end
		function set.COL_KAPPA_CONST(obj, val)
			obj.kappaConstant = val;
		end
		function val = get.COL_CORDNUM(obj)
			val = obj.cordNum;
		end
		function set.COL_CORDNUM(obj, val)
			obj.cordNum = val;
		end
		function val = get.COL_LOGKEQ_PH_EXP(obj)
			val = obj.logKeqPhExponent;
		end
		function set.COL_LOGKEQ_PH_EXP(obj, val)
			obj.logKeqPhExponent = val;
		end
		function val = get.COL_LOGKEQ_SALT_POWEXP(obj)
			val = obj.logKeqSaltPowerLawExponent;
		end
		function set.COL_LOGKEQ_SALT_POWEXP(obj, val)
			obj.logKeqSaltPowerLawExponent = val;
		end
		function val = get.COL_LOGKEQ_SALT_POWFACT(obj)
			val = obj.logKeqSaltPowerLawFactor;
		end
		function set.COL_LOGKEQ_SALT_POWFACT(obj, val)
			obj.logKeqSaltPowerLawFactor = val;
		end
		function val = get.COL_LOGKEQ_SALT_EXPFACT(obj)
			val = obj.logKeqSaltExpLawFactor;
		end
		function set.COL_LOGKEQ_SALT_EXPFACT(obj, val)
			obj.logKeqSaltExpLawFactor = val;
		end
		function val = get.COL_LOGKEQ_SALT_EXPARGMULT(obj)
			val = obj.logKeqSaltExpLawExponentFactor;
		end
		function set.COL_LOGKEQ_SALT_EXPARGMULT(obj, val)
			obj.logKeqSaltExpLawExponentFactor = val;
		end
		function val = get.COL_BPP_PH_EXP(obj)
			val = obj.bppPhExponent;
		end
		function set.COL_BPP_PH_EXP(obj, val)
			obj.bppPhExponent = val;
		end
		function val = get.COL_BPP_SALT_POWEXP(obj)
			val = obj.bppSaltPowerLawExponent;
		end
		function set.COL_BPP_SALT_POWEXP(obj, val)
			obj.bppSaltPowerLawExponent = val;
		end
		function val = get.COL_BPP_SALT_POWFACT(obj)
			val = obj.bppSaltPowerLawFactor;
		end
		function set.COL_BPP_SALT_POWFACT(obj, val)
			obj.bppSaltPowerLawFactor = val;
		end
		function val = get.COL_BPP_SALT_EXPFACT(obj)
			val = obj.bppSaltExpLawFactor;
		end
		function set.COL_BPP_SALT_EXPFACT(obj, val)
			obj.bppSaltExpLawFactor = val;
		end
		function val = get.COL_BPP_SALT_EXPARGMULT(obj)
			val = obj.bppSaltExpLawExponentFactor;
		end
		function set.COL_BPP_SALT_EXPARGMULT(obj, val)
			obj.bppSaltExpLawExponentFactor = val;
		end
		function val = get.COL_RADIUS(obj)
			val = obj.radius;
		end
		function set.COL_RADIUS(obj, val)
			obj.radius = val;
		end
		function val = get.COL_KKIN(obj)
			val = obj.kKin;
		end
		function set.COL_KKIN(obj, val)
			obj.kKin = val;
		end
		function val = get.COL_LINEAR_THRESHOLD(obj)
			val = obj.linearizationThreshold;
		end
		function set.COL_LINEAR_THRESHOLD(obj, val)
			obj.linearizationThreshold = val;
		end

	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = ColloidalBinding();
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
