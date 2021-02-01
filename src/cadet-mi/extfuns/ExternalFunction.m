
classdef ExternalFunction < handle & matlab.mixin.Heterogeneous
	%ExternalFunction Base class for external functions
	%   Provides properties and functionalities common to all external
	%   functions.
	%
	%   Derived classes are supposed to use the field 'data' for storing
	%   their configuration conforming to the CADET file format spec.
	%
	%   Changes in the parameters are tracked using the field hasChanged.
	%   If a value is changed using the properties, hasChanged is set to
	%   true. It is reset by calling notifySync(). Derived classes have to
	%   take care of monitoring property changes and updating hasChanged.
	%
	% See also MODELSYSTEM
	
	% Copyright: (C) 2008-2021 The CADET Authors
	%            See the license note at the end of the file.
	
	properties (Hidden, Access = 'protected')
		data; % Struct with stored property values
	end

	properties (Transient)
		hasChanged; % Determines whether this object has changed after the last synchronization with CADET (C++ layer)
	end

	properties (Constant, Transient, Abstract)
		name; % Name of the external function according to CADET file format specs
	end

	methods

		function obj = ExternalFunction()
			%EXTERNALFUNCTION Constructs an external function base class object

			obj.hasChanged = true;
			obj.data.EXTFUN_TYPE = obj.name;
		end

		function res = validate(obj, sectionTimes)
			%VALIDATE Validates the configuration of the external function
			%   RES = VALIDATE(SECTIONTIMES) uses the section times of the simulator given in
			%   SECTIONTIMES to validate the external function object. Returns true in RES if
			%   everything is fine and false otherwise.
			%
			%   Derived classes are supposed to overwrite this method.

			res = true;
		end

		function res = assembleConfig(obj)
			%ASSEMBLECONFIG Assembles the configuration according to the CADET file format spec
			%   RES = ASSEMBLECONFIG() returns a Matlab struct RES that represents the
			%   external function as detailed in the CADET file format spec.
			%
			%   External functions are supposed to store their configuration (conforming to
			%   the CADET file format spec) in the data field of this base class. However, by
			%   overwriting this function, derived classes can customize the assembly process.
			%
			% See also MODELSYSTEM.ASSEMBLECONFIG

			res = obj.data;
		end

		function notifySync(obj)
			%NOTIFYSYNC Called after parameters have been sent to CADET
			%   NOTIFYSYNC() resets the HASCHANGED property.

			obj.hasChanged = false;
		end

		function S = saveobj(obj)
			S.data = obj.data;
		end

	end

	methods (Abstract)
		val = evaluate(obj, t, z, r)
			%EVALUATE Evaluates the external function
			%   VAL = EVALUATE(T, Z, R) evaluates the external function at time T,
			%   normalized axial position Z (in [0, 1]), and normalized radial position
			%   R (in [0, 1]).
	end

	methods (Access = 'protected')

		function loadobjInternal(obj, S)
			obj.hasChanged = false;
			obj.data = S.data;
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
