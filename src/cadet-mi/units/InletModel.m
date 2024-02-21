
classdef (Abstract) InletModel < Model
	%InletModel Represents an inlet in a system of unit operations
	%   An inlet is a pseudo unit operation that injects mass into the
	%   system. The INLETMODEL serves as a base class for all types of
	%   inlet models that differ in the injection profile.
	%
	% See also MODEL, MODELSYSTEM

	% Copyright: (C) 2008-2024 The CADET Authors
	%            See the license note at the end of the file.
	
	properties(Constant)
		name = 'INLET'; % Type of the model according to CADET file format specs
		hasInlet = false; % Determines whether the unit operation has an inlet
		hasOutlet = true; % Determines whether the unit operation has an outlet
	end

	properties(Dependent, Transient)
		nInletPorts; % Number of inlet ports
		nOutletPorts; % Number of outlet ports
	end

	properties (Constant, Access = 'protected')
		hasConsistencySolver = false; % Determines whether this unit operation model has a consistency solver
	end

	methods
		
		function obj = InletModel()
			%INLETMODEL Constructs an InletModel object and inserts as much default values as possible

			obj = obj@Model();
		end
		
		function res = validate(obj, sectionTimes)
			%VALIDATE Validates the configuration of the Model
			%   RES = VALIDATE(SECTIONTIMES) uses the section times of the simulator given in
			%   SECTIONTIMES to validate the Model object. Returns true in RES if everything
			%   is fine and false otherwise.
			%
			%   Derived classes are supposed to extend this function.
			%
			% See also MODEL.VALIDATE

			res = obj.validate@Model(sectionTimes);
		end

		function res = assembleConfig(obj)
			%ASSEMBLECONFIG Assembles the configuration according to the CADET file format spec
			%   RES = ASSEMBLECONFIG() returns a nested Matlab struct RES that represents the
			%   model as detailed in the CADET file format spec.
			%
			% See also MODEL.ASSEMBLECONFIG, MEXSIMULATOR.ASSEMBLECONFIG

			res = obj.assembleConfig@Model();
		end

		function notifySync(obj)
			%NOTIFYSYNC Called after parameters have been sent to CADET
			%   NOTIFYSYNC() resets the HASCHANGED property.

			obj.notifySync@Model();
		end

		function val = get.nInletPorts(obj)
			val = 0;
		end

		function val = get.nOutletPorts(obj)
			val = 1;
		end

	end
	
end

% =============================================================================
%  CADET
%  
%  Copyright (C) 2008-2024: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
