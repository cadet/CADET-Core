
classdef SimpleMultiStateStericMassActionBinding < KineticQuasiStationaryBindingModel
	%SimpleMultiStateStericMassActionBinding Simplified multi-state steric mass action binding model
	%
	% See also BINDINGMODEL, KINETICQUASISTATIONARYBINDINGMODEL

	% Copyright: (C) 2008-2019 The CADET Authors
	%            See the license note at the end of the file.

	properties (Constant, Access = 'protected')
		hasConsistencySolver = true; % Determines whether this binding model has a consistency solver
	end

	properties(Constant)
		name = 'SIMPLE_MULTISTATE_STERIC_MASS_ACTION'; % Name of the binding model according to CADET file format specs
	end

	properties (Dependent, Transient)
		% Adsorption rate in [m^3_MP / (m^3_SP * s)]
		kA;
		SMSSMA_KA;
		% Desorption rate in [1 / s]
		kD;
		SMSSMA_KD;
		% Ionic capacity in [mol / m^3_SP]
		lambda;
		SMSSMA_LAMBDA;
		% Minimal characteristic charge [-] (per component)
		nuMin;
		SMSSMA_NU_MIN;
		% Maximal characteristic charge [-] (per component)
		nuMax;
		SMSSMA_NU_MAX;
		% Quadratic factor of characteristic charge [-] (per component)
		nuQuad;
		SMSSMA_NU_QUAD;
		% Minimal steric factor [-] (per component)
		sigmaMin;
		SMSSMA_SIGMA_MIN;
		% Maximal steric factor [-] (per component)
		sigmaMax;
		SMSSMA_SIGMA_MAX;
		% Quadratic factor steric factor [-] (per component)
		sigmaQuad;
		SMSSMA_SIGMA_QUAD;
		% Conversion rates [1 / s] from weakly bound states to strongly bound states (per component)
		weakToStrong;
		SMSSMA_KWS;
		% Linear factor of conversion rates [1 / s] from weakly bound states to strongly bound states (per component)
		weakToStrongLin;
		SMSSMA_KWS_LIN;
		% Quadratic factor of conversion rates [1 / s] from weakly bound states to strongly bound states
		weakToStrongQuad;
		SMSSMA_KWS_QUAD;
		% Conversion rates [1 / s] from strongly bound states to weakly bound states (per component)
		strongToWeak;
		SMSSMA_KSW;
		% Linear factor of conversion rates [1 / s] from strongly bound states to weakly bound states (per component)
		strongToWeakLin;
		SMSSMA_KSW_LIN;
		% Quadratic factor of conversion rates [1 / s] from strongly bound states to weakly bound states (per component)
		strongToWeakQuad;
		SMSSMA_KSW_QUAD;
		% Reference solid phase concentration [mol / m^3_SP]
		refSolid;
		SMSSMA_REFQ;
		% Reference liquid phase concentration [mol / m^3_MP]
		refLiquid;
		SMSSMA_REFC0;
	end

	methods

		function obj = SimpleMultiStateStericMassActionBinding(kA, kD)
			%SIMPLEMULTISTATESTERICMASSACTIONBINDING Constructs a SimpleMultiStateStericMassActionBinding object
			%   OBJ = SIMPLEMULTISTATESTERICMASSACTIONBINDING(KA) creates a SimpleMultiStateStericMassActionBinding
			%   model with the given adsorption rates KA.
			%
			%   OBJ = SIMPLEMULTISTATESTERICMASSACTIONBINDING(..., KD) also sets the desorption rates to KD.

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
			validateattributes(obj.nuMin, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'nuMin');
			validateattributes(obj.nuMax, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'nuMax');
			validateattributes(obj.nuQuad, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'nuQuad');
			validateattributes(obj.sigmaMin, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'sigmaMin');
			validateattributes(obj.sigmaMax, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'sigmaMax');
			validateattributes(obj.sigmaQuad, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'sigmaQuad');
			validateattributes(obj.weakToStrong, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'weakToStrong');
			validateattributes(obj.weakToStrongLin, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'weakToStrongLin');
			validateattributes(obj.weakToStrongQuad, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'weakToStrongQuad');
			validateattributes(obj.strongToWeak, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real', 'numel', nComponents}, '', 'strongToWeak');
			validateattributes(obj.strongToWeakLin, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'strongToWeakLin');
			validateattributes(obj.strongToWeakQuad, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nComponents}, '', 'strongToWeakQuad');
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
			val = obj.data.SMSSMA_KA;
		end

		function set.kA(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'kA');
			obj.data.SMSSMA_KA = val;
			obj.hasChanged = true;
		end

		function val = get.kD(obj)
			val = obj.data.SMSSMA_KD;
		end

		function set.kD(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'kD');
			obj.data.SMSSMA_KD = val;
			obj.hasChanged = true;
		end

		function val = get.nuMin(obj)
			val = obj.data.SMSSMA_NU_MIN;
		end

		function set.nuMin(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'nuMin');
			obj.data.SMSSMA_NU_MIN = val;
			obj.hasChanged = true;
		end

		function val = get.nuMax(obj)
			val = obj.data.SMSSMA_NU_MAX;
		end

		function set.nuMax(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'nuMax');
			obj.data.SMSSMA_NU_MAX = val;
			obj.hasChanged = true;
		end

		function val = get.nuQuad(obj)
			val = obj.data.SMSSMA_NU_QUAD;
		end

		function set.nuQuad(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'nuQuad');
			obj.data.SMSSMA_NU_QUAD = val;
			obj.hasChanged = true;
		end

		function val = get.sigmaMin(obj)
			val = obj.data.SMSSMA_SIGMA_MIN;
		end

		function set.sigmaMin(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'sigmaMin');
			obj.data.SMSSMA_SIGMA_MIN = val;
			obj.hasChanged = true;
		end

		function val = get.sigmaMax(obj)
			val = obj.data.SMSSMA_SIGMA_MAX;
		end

		function set.sigmaMax(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'sigmaMax');
			obj.data.SMSSMA_SIGMA_MAX = val;
			obj.hasChanged = true;
		end

		function val = get.sigmaQuad(obj)
			val = obj.data.SMSSMA_SIGMA_QUAD;
		end

		function set.sigmaQuad(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'sigmaQuad');
			obj.data.SMSSMA_SIGMA_QUAD = val;
			obj.hasChanged = true;
		end

		function val = get.lambda(obj)
			val = obj.data.SMSSMA_LAMBDA;
		end

		function set.lambda(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'lambda');
			obj.data.SMSSMA_LAMBDA = val;
			obj.hasChanged = true;
		end

		function val = get.weakToStrong(obj)
			val = obj.data.SMSSMA_KWS;
		end

		function set.weakToStrong(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'weakToStrong');
			obj.data.SMSSMA_KWS = val;
			obj.hasChanged = true;
		end

		function val = get.weakToStrongLin(obj)
			val = obj.data.SMSSMA_KWS_LIN;
		end

		function set.weakToStrongLin(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'weakToStrongLin');
			obj.data.SMSSMA_KWS_LIN = val;
			obj.hasChanged = true;
		end

		function val = get.weakToStrongQuad(obj)
			val = obj.data.SMSSMA_KWS_QUAD;
		end

		function set.weakToStrongQuad(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'weakToStrongQuad');
			obj.data.SMSSMA_KWS_QUAD = val;
			obj.hasChanged = true;
		end

		function val = get.strongToWeak(obj)
			val = obj.data.SMSSMA_KSW;
		end

		function set.strongToWeak(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'strongToWeak');
			obj.data.SMSSMA_KSW = val;
			obj.hasChanged = true;
		end

		function val = get.strongToWeakLin(obj)
			val = obj.data.SMSSMA_KSW_LIN;
		end

		function set.strongToWeakLin(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'strongToWeakLin');
			obj.data.SMSSMA_KSW_LIN = val;
			obj.hasChanged = true;
		end

		function val = get.strongToWeakQuad(obj)
			val = obj.data.SMSSMA_KSW_QUAD;
		end

		function set.strongToWeakQuad(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'strongToWeakQuad');
			obj.data.SMSSMA_KSW_QUAD = val;
			obj.hasChanged = true;
		end

		function val = get.refSolid(obj)
			if isfield(obj.data, 'SMSSMA_REFQ')
				val = obj.data.SMSSMA_REFQ;
			else
				val = [];
			end
		end

		function set.refSolid(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'SMSSMA_REFQ');
			else
				validateattributes(val, {'double'}, {'scalar', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refSolid');
				obj.data.SMSSMA_REFQ = val;
			end
			obj.hasChanged = true;
		end

		function val = get.refLiquid(obj)
			if isfield(obj.data, 'SMSSMA_REFC0')
				val = obj.data.SMSSMA_REFC0;
			else
				val = [];
			end
		end

		function set.refLiquid(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'SMSSMA_REFC0');
			else
				validateattributes(val, {'double'}, {'scalar', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'refLiquid');
				obj.data.SMSSMA_REFC0 = val;
			end
			obj.hasChanged = true;
		end


		function val = get.SMSSMA_KA(obj)
			val = obj.kA;
		end
		function set.SMSSMA_KA(obj, val)
			obj.kA = val;
		end
		function val = get.SMSSMA_KD(obj)
			val = obj.kD;
		end
		function set.SMSSMA_KD(obj, val)
			obj.kD = val;
		end
		function val = get.SMSSMA_LAMBDA(obj)
			val = obj.lambda;
		end
		function set.SMSSMA_LAMBDA(obj, val)
			obj.lambda = val;
		end
		function val = get.SMSSMA_NU_MIN(obj)
			val = obj.nuMin;
		end
		function set.SMSSMA_NU_MIN(obj, val)
			obj.nuMin = val;
		end
		function val = get.SMSSMA_NU_MAX(obj)
			val = obj.nuMax;
		end
		function set.SMSSMA_NU_MAX(obj, val)
			obj.nuMax = val;
		end
		function val = get.SMSSMA_NU_QUAD(obj)
			val = obj.nuQuad;
		end
		function set.SMSSMA_NU_QUAD(obj, val)
			obj.nuQuad = val;
		end
		function val = get.SMSSMA_SIGMA_MIN(obj)
			val = obj.sigmaMin;
		end
		function set.SMSSMA_SIGMA_MIN(obj, val)
			obj.sigmaMin = val;
		end
		function val = get.SMSSMA_SIGMA_MAX(obj)
			val = obj.sigmaMax;
		end
		function set.SMSSMA_SIGMA_MAX(obj, val)
			obj.sigmaMax = val;
		end
		function val = get.SMSSMA_SIGMA_QUAD(obj)
			val = obj.sigmaQuad;
		end
		function set.SMSSMA_SIGMA_QUAD(obj, val)
			obj.sigmaQuad = val;
		end
		function val = get.SMSSMA_KWS(obj)
			val = obj.weakToStrong;
		end
		function set.SMSSMA_KWS(obj, val)
			obj.weakToStrong = val;
		end
		function val = get.SMSSMA_KWS_LIN(obj)
			val = obj.weakToStrongLin;
		end
		function set.SMSSMA_KWS_LIN(obj, val)
			obj.weakToStrongLin = val;
		end
		function val = get.SMSSMA_KWS_QUAD(obj)
			val = obj.weakToStrongQuad;
		end
		function set.SMSSMA_KWS_QUAD(obj, val)
			obj.weakToStrongQuad = val;
		end
		function val = get.SMSSMA_KSW(obj)
			val = obj.strongToWeak;
		end
		function set.SMSSMA_KSW(obj, val)
			obj.strongToWeak = val;
		end
		function val = get.SMSSMA_KSW_LIN(obj)
			val = obj.strongToWeakLin;
		end
		function set.SMSSMA_KSW_LIN(obj, val)
			obj.strongToWeakLin = val;
		end
		function val = get.SMSSMA_KSW_QUAD(obj)
			val = obj.strongToWeakQuad;
		end
		function set.SMSSMA_KSW_QUAD(obj, val)
			obj.strongToWeakQuad = val;
		end
		function val = get.SMSSMA_REFQ(obj)
			val = obj.refSolid;
		end
		function set.SMSSMA_REFQ(obj, val)
			obj.refSolid = val;
		end
		function val = get.SMSSMA_REFC0(obj)
			val = obj.refLiquid;
		end
		function set.SMSSMA_REFC0(obj, val)
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
				obj = SimpleMultiStateStericMassActionBinding();
				obj.loadobjInternal(S);
			end
		end
		
	end
end

% =============================================================================
%  CADET - The Chromatography Analysis and Design Toolkit
%  
%  Copyright (C) 2008-2019: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. obj program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies obj distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
