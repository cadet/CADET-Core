
classdef BiLangmuirBinding < KineticQuasiStationaryBindingModel
	%BiLangmuirBinding Bi-Langmuir binding model
	%
	% See also BINDINGMODEL, KINETICQUASISTATIONARYBINDINGMODEL
	
	% Copyright: (C) 2008-2020 The CADET Authors
	%            See the license note at the end of the file.

	properties(Constant)
		name = 'MULTI_COMPONENT_BILANGMUIR'; % Name of the binding model according to CADET file format specs
	end

	properties (Dependent, Transient)
		% Adsorption rate in [m^3_MP / (mol * s)], rows correspond to binding sites
		kA;
		MCBL_KA;
		% Desorption rate in [1 / s], rows correspond to binding sites
		kD;
		MCBL_KD;
		% Capacity in [mol / m^3_SP], rows correspond to binding sites
		qMax;
		MCBL_QMAX;
	end

	methods

		function obj = BiLangmuirBinding(kA, kD, qMax)
			%BILANGMUIRBINDING Constructs a Bi-Langmuir binding model object
			%   OBJ = BILANGMUIRBINDING(KA) Constructs a Bi-Langmuir binding model
			%   object using the given adsorption rates KA.
			%			
			%   OBJ = BILANGMUIRBINDING(..., KD) Constructs a Bi-Langmuir binding
			%   model object using the given adsorption rates KA and desorption rates
			%   KD.
			%			
			%   OBJ = BILANGMUIRBINDING(..., KD, QMAX) Constructs a Bi-Langmuir binding
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

			% In kA, kD, qMax each row corresponds to a binding site
			validateattributes(obj.kA, {'double'}, {'2d', 'nonempty', '>=', 0.0, 'finite', 'real', 'size', [numStates, nComponents]}, '', 'kA');
			validateattributes(obj.kD, {'double'}, {'2d', 'nonempty', '>=', 0.0, 'finite', 'real', 'size', [numStates, nComponents]}, '', 'kD');
			validateattributes(obj.qMax, {'double'}, {'2d', 'nonempty', '>', 0.0, 'finite', 'real', 'size', [numStates, nComponents]}, '', 'qMax');
			res = obj.validate@KineticQuasiStationaryBindingModel(nComponents, nBoundStates);
		end

		function val = get.kA(obj)
			val = obj.data.MCBL_KA.';
		end

		function set.kA(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'kA');
			obj.data.MCBL_KA = val.';
			obj.hasChanged = true;
		end

		function val = get.kD(obj)
			val = obj.data.MCBL_KD.';
		end

		function set.kD(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'kD');
			obj.data.MCBL_KD = val.';
			obj.hasChanged = true;
		end

		function val = get.qMax(obj)
			val = obj.data.MCBL_QMAX.';
		end

		function set.qMax(obj, val)
			validateattributes(val, {'double'}, {'2d', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'qMax');
			obj.data.MCBL_QMAX = val.';
			obj.hasChanged = true;
		end


		function val = get.MCBL_KA(obj)
			val = obj.kA;
		end
		function set.MCBL_KA(obj, val)
			obj.kA = val;
		end
		function val = get.MCBL_KD(obj)
			val = obj.kD;
		end
		function set.MCBL_KD(obj, val)
			obj.kD = val;
		end
		function val = get.MCBL_QMAX(obj)
			val = obj.qMax;
		end
		function set.MCBL_QMAX(obj, val)
			obj.qMax = val;
		end

		function res = assembleConfig(obj)
			%ASSEMBLECONFIG Assembles the configuration according to the CADET file format spec
			%   RES = ASSEMBLECONFIG() returns a Matlab struct RES that represents the
			%   binding model as detailed in the CADET file format spec.
			%
			% See also BINDINGMODEL.ASSEMBLECONFIG

			res = obj.data;

			% Linearize
			res.MCBL_KA = obj.data.MCBL_KA(:);
			res.MCBL_KD = obj.data.MCBL_KD(:);
			res.MCBL_QMAX = obj.data.MCBL_QMAX(:);
		end

	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = BiLangmuirBinding();
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
