
classdef ExtFunBiLangmuirBinding < KineticQuasiStationaryBindingModel
	%ExtFunBiLangmuirBinding Bi-Langmuir binding model with external function support
	%
	% See also BINDINGMODEL, BILANGMUIRBINDING, KINETICQUASISTATIONARYBINDINGMODEL
	
	% Copyright: (C) 2008-2020 The CADET Authors
	%            See the license note at the end of the file.

	properties(Constant)
		name = 'EXT_MULTI_COMPONENT_BILANGMUIR'; % Name of the binding model according to CADET file format specs
	end

	properties (Dependent, Transient)
		% Adsorption rate in [m^3_MP / (mol * s)], rows correspond to binding sites, constant term of external dependence
		kA;
		EXT_MCBL_KA;
		% Adsorption rate in [m^3_MP / (mol * s * [T])], rows correspond to binding sites, linear term of external dependence
		kA_T;
		EXT_MCBL_KA_T;
		% Adsorption rate in [m^3_MP / (mol * s * [T^2])], rows correspond to binding sites, quadratic term of external dependence
		kA_TT;
		EXT_MCBL_KA_TT;
		% Adsorption rate in [m^3_MP / (mol * s * [T^3])], rows correspond to binding sites, cubic term of external dependence
		kA_TTT;
		EXT_MCBL_KA_TTT;
		% Desorption rate in [1 / s], rows correspond to binding sites, constant term of external dependence
		kD;
		EXT_MCBL_KD;
		% Desorption rate in [1 / (s * [T])], rows correspond to binding sites, linear term of external dependence
		kD_T;
		EXT_MCBL_KD_T;
		% Desorption rate in [1 / (s * [T]^2)], rows correspond to binding sites, quadratic term of external dependence
		kD_TT;
		EXT_MCBL_KD_TT;
		% Desorption rate in [1 / (s * [T]^3)], rows correspond to binding sites, cubic term of external dependence
		kD_TTT;
		EXT_MCBL_KD_TTT;
		% Capacity in [mol / m^3_SP], rows correspond to binding sites, constant term of external dependence
		qMax;
		EXT_MCBL_QMAX;
		% Capacity in [mol / (m^3_SP * [T])], rows correspond to binding sites, linear term of external dependence
		qMax_T;
		EXT_MCBL_QMAX_T;
		% Capacity in [mol / (m^3_SP * [T]^2)], rows correspond to binding sites, quadratic term of external dependence
		qMax_TT;
		EXT_MCBL_QMAX_TT;
		% Capacity in [mol / (m^3_SP * [T]^3)], rows correspond to binding sites, cubic term of external dependence
		qMax_TTT;
		EXT_MCBL_QMAX_TTT;
		% Indices of external functions (0-based)
		externalSource;
		EXTFUN;
	end

	methods

		function obj = ExtFunBiLangmuirBinding(kA, kD, qMax)
			%EXTFUNBILANGMUIRBINDING Constructs a Bi-Langmuir binding model object with external function support
			%   OBJ = EXTFUNBILANGMUIRBINDING(KA) Constructs a Bi-Langmuir binding model
			%   object using the given adsorption rates KA.
			%			
			%   OBJ = EXTFUNBILANGMUIRBINDING(..., KD) Constructs a Bi-Langmuir binding
			%   model object using the given adsorption rates KA and desorption rates
			%   KD.
			%			
			%   OBJ = EXTFUNBILANGMUIRBINDING(..., KD, QMAX) Constructs a Bi-Langmuir binding
			%   model object using the given adsorption rates KA, desorption rates KD,
			%   and capacity QMAX.

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

			if any((nBoundStates ~= 0) & (nBoundStates ~= max(nBoundStates)))
				error('CADET:invalidConfig', 'Expected all components to have the same number of bound states (or 0).');
			end
			numStates = max(nBoundStates);

			% Each row corresponds to a binding site
			validateattributes(obj.kA, {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'kA');
			validateattributes(obj.kD, {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'kD');
			validateattributes(obj.qMax, {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'qMax');

			validateattributes(obj.kA_T, {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'kA_T');
			validateattributes(obj.kD_T, {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'kD_T');
			validateattributes(obj.qMax_T, {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'qMax_T');

			validateattributes(obj.kA_TT, {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'kA_TT');
			validateattributes(obj.kD_TT, {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'kD_TT');
			validateattributes(obj.qMax_TT, {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'qMax_TT');

			validateattributes(obj.kA_TTT, {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'kA_TTT');
			validateattributes(obj.kD_TTT, {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'kD_TTT');
			validateattributes(obj.qMax_TTT, {'double'}, {'2d', 'nonempty', 'finite', 'real', 'size', [numStates, nComponents]}, '', 'qMax_TTT');

			if ~isempty(obj.externalSource)
				validateattributes(obj.externalSource, {'double', 'int32'}, {'vector', 'nonempty', 'finite', 'real', '>=', 0}, '', 'externalSource');
				if ~any(numel(obj.externalSource) == [1, 3])
					error('CADET:invalidConfig', 'Expected externalSource to have either 1 or 3 entries but got %d entries.', numel(obj.externalSource));
				end
			end

			res = obj.validate@KineticQuasiStationaryBindingModel(nComponents, nBoundStates);
		end

		function val = get.kA(obj)
			val = obj.data.EXT_MCBL_KA.';
		end

		function set.kA(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'kA');
			obj.data.EXT_MCBL_KA = val.';
			obj.hasChanged = true;
		end

		function val = get.kD(obj)
			val = obj.data.EXT_MCBL_KD.';
		end

		function set.kD(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'kD');
			obj.data.EXT_MCBL_KD = val.';
			obj.hasChanged = true;
		end

		function val = get.qMax(obj)
			val = obj.data.EXT_MCBL_QMAX.';
		end

		function set.qMax(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'qMax');
			obj.data.EXT_MCBL_QMAX = val.';
			obj.hasChanged = true;
		end

		function val = get.kA_T(obj)
			val = obj.data.EXT_MCBL_KA_T.';
		end

		function set.kA_T(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'kA_T');
			obj.data.EXT_MCBL_KA_T = val.';
			obj.hasChanged = true;
		end

		function val = get.kD_T(obj)
			val = obj.data.EXT_MCBL_KD_T.';
		end

		function set.kD_T(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'kD_T');
			obj.data.EXT_MCBL_KD_T = val.';
			obj.hasChanged = true;
		end

		function val = get.qMax_T(obj)
			val = obj.data.EXT_MCBL_QMAX_T.';
		end

		function set.qMax_T(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'qMax_T');
			obj.data.EXT_MCBL_QMAX_T = val.';
			obj.hasChanged = true;
		end

		function val = get.kA_TT(obj)
			val = obj.data.EXT_MCBL_KA_TT.';
		end

		function set.kA_TT(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'kA_TT');
			obj.data.EXT_MCBL_KA_TT = val.';
			obj.hasChanged = true;
		end

		function val = get.kD_TT(obj)
			val = obj.data.EXT_MCBL_KD_TT.';
		end

		function set.kD_TT(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'kD_TT');
			obj.data.EXT_MCBL_KD_TT = val.';
			obj.hasChanged = true;
		end

		function val = get.qMax_TT(obj)
			val = obj.data.EXT_MCBL_QMAX_TT.';
		end

		function set.qMax_TT(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'qMax_TT');
			obj.data.EXT_MCBL_QMAX_TT = val.';
			obj.hasChanged = true;
		end

		function val = get.kA_TTT(obj)
			val = obj.data.EXT_MCBL_KA_TTT.';
		end

		function set.kA_TTT(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'kA_TTT');
			obj.data.EXT_MCBL_KA_TTT = val.';
			obj.hasChanged = true;
		end

		function val = get.kD_TTT(obj)
			val = obj.data.EXT_MCBL_KD_TTT.';
		end

		function set.kD_TTT(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'kD_TTT');
			obj.data.EXT_MCBL_KD_TTT = val.';
			obj.hasChanged = true;
		end

		function val = get.qMax_TTT(obj)
			val = obj.data.EXT_MCBL_QMAX_TTT.';
		end

		function set.qMax_TTT(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', 'finite', 'real'}, '', 'qMax_TTT');
			obj.data.EXT_MCBL_QMAX_TTT = val.';
			obj.hasChanged = true;
		end

		function val = get.EXT_MCBL_KA(obj)
			val = obj.kA;
		end
		function set.EXT_MCBL_KA(obj, val)
			obj.kA = val;
		end
		function val = get.EXT_MCBL_KD(obj)
			val = obj.kD;
		end
		function set.EXT_MCBL_KD(obj, val)
			obj.kD = val;
		end
		function val = get.EXT_MCBL_QMAX(obj)
			val = obj.qMax;
		end
		function set.EXT_MCBL_QMAX(obj, val)
			obj.qMax = val;
		end

		function val = get.EXT_MCBL_KA_T(obj)
			val = obj.kA_T;
		end
		function set.EXT_MCBL_KA_T(obj, val)
			obj.kA_T = val;
		end
		function val = get.EXT_MCBL_KD_T(obj)
			val = obj.kD_T;
		end
		function set.EXT_MCBL_KD_T(obj, val)
			obj.kD_T = val;
		end
		function val = get.EXT_MCBL_QMAX_T(obj)
			val = obj.qMax_T;
		end
		function set.EXT_MCBL_QMAX_T(obj, val)
			obj.qMax_T = val;
		end

		function val = get.EXT_MCBL_KA_TT(obj)
			val = obj.kA_TT;
		end
		function set.EXT_MCBL_KA_TT(obj, val)
			obj.kA_TT = val;
		end
		function val = get.EXT_MCBL_KD_TT(obj)
			val = obj.kD_TT;
		end
		function set.EXT_MCBL_KD_TT(obj, val)
			obj.kD_TT = val;
		end
		function val = get.EXT_MCBL_QMAX_TT(obj)
			val = obj.qMax_TT;
		end
		function set.EXT_MCBL_QMAX_TT(obj, val)
			obj.qMax_TT = val;
		end

		function val = get.EXT_MCBL_KA_TTT(obj)
			val = obj.kA_TTT;
		end
		function set.EXT_MCBL_KA_TTT(obj, val)
			obj.kA_TTT = val;
		end
		function val = get.EXT_MCBL_KD_TTT(obj)
			val = obj.kD_TTT;
		end
		function set.EXT_MCBL_KD_TTT(obj, val)
			obj.kD_TTT = val;
		end
		function val = get.EXT_MCBL_QMAX_TTT(obj)
			val = obj.qMax_TTT;
		end
		function set.EXT_MCBL_QMAX_TTT(obj, val)
			obj.qMax_TTT = val;
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
				if ~any(numel(val) == [1, 3])
					error('CADET:invalidConfig', 'Expected externalSource to have either 1 or 3 entries but got %d entries.', numel(val));
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
			res.EXT_MCBL_KA = obj.data.EXT_MCBL_KA(:);
			res.EXT_MCBL_KD = obj.data.EXT_MCBL_KD(:);
			res.EXT_MCBL_QMAX = obj.data.EXT_MCBL_QMAX(:);
			res.EXT_MCBL_KA_T = obj.data.EXT_MCBL_KA_T(:);
			res.EXT_MCBL_KD_T = obj.data.EXT_MCBL_KD_T(:);
			res.EXT_MCBL_QMAX_T = obj.data.EXT_MCBL_QMAX_T(:);
			res.EXT_MCBL_KA_TT = obj.data.EXT_MCBL_KA_TT(:);
			res.EXT_MCBL_KD_TT = obj.data.EXT_MCBL_KD_TT(:);
			res.EXT_MCBL_QMAX_TT = obj.data.EXT_MCBL_QMAX_TT(:);
			res.EXT_MCBL_KA_TTT = obj.data.EXT_MCBL_KA_TTT(:);
			res.EXT_MCBL_KD_TTT = obj.data.EXT_MCBL_KD_TTT(:);
			res.EXT_MCBL_QMAX_TTT = obj.data.EXT_MCBL_QMAX_TTT(:);
		end

	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = ExtFunBiLangmuirBinding();
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
