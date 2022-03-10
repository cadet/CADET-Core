
classdef KineticQuasiStationaryBindingModel < BindingModel
	%KineticQuasiStationaryBindingModel Base class for binding models with kinetic and quasi-stationary binding
	%   Represents binding models whose binding can either by kinetic or quasi-stationary.
	%
	% See also BINDINGMODEL

	% Copyright: (C) 2008-2022 The CADET Authors
	%            See the license note at the end of the file.
	
	properties (Dependent, Transient)
		% Determines whether kinetic binding mode is used (true) or quasi-stationary (false)
		kineticBinding;
	end

	methods

		function obj = KineticQuasiStationaryBindingModel()
			%KINETICQUASISTATIONARYBINDINGMODEL Creates an object of the KineticQuasiStationaryBindingModel
			obj = obj@BindingModel();
		end

		function val = get.kineticBinding(obj)
			val = logical(obj.data.IS_KINETIC);
		end

		function set.kineticBinding(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'kineticBinding');
			obj.data.IS_KINETIC = int32(logical(val));
			obj.hasChanged = true;
		end

		function res = validate(obj, nComponents, nBoundStates)
			%VALIDATE Validates the binding model parameters
			%   RES = VALIDATE(NCOMPONENTS, NBOUNDSTATES) uses the number of
			%   components NCOMPONENTS of the model and the number of bound
			%   states of each of these components as given by the vector
			%   NBOUNDSTATES to validate the parameters of the binding model.
			%   Returns true in RES if everything is fine and false otherwise.
			%
			%   Derived classes are supposed to extend this method.

			if ~isfield(obj.data, 'IS_KINETIC')
				error('CADET:invalidConfig', 'Expected kineticBinding to be set.');
			end

			res = obj.validate@BindingModel(nComponents, nBoundStates);
		end

	end

end

% =============================================================================
%  CADET
%  
%  Copyright (C) 2008-2022: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. obj program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies obj distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
