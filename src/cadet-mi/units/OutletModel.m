
classdef OutletModel < Model
	%OutletModel Represents an outlet in a system of unit operations
	%   An outlet is a pseudo unit operation that represents a physical
	%   outlet that removes mass from the system.
	%
	%   Note that outlets are not necessary if they are simply a leaf
	%   or terminal node in a network of unit operations (e.g., a simple
	%   chain of unit operations does not require an outlet at its end).
	%
	% See also MODEL, MODELSYSTEM
	
	% Copyright: (C) 2008-2020 The CADET Authors
	%            See the license note at the end of the file.
	
	properties(Constant)
		name = 'OUTLET'; % Type of the model according to CADET file format specs
		hasInlet = true; % Determines whether the unit operation has an inlet
		hasOutlet = false; % Determines whether the unit operation has an outlet
	end

	properties(Dependent, Transient)
		nInletPorts; % Number of inlet ports
		nOutletPorts; % Number of outlet ports
	end

	properties (Constant, Access = 'protected')
		hasConsistencySolver = false; % Determines whether this unit operation model has a consistency solver
	end

	methods
		
		function obj = OutletModel()
			%OUTLETMODEL Constructs an OutletModel object

			obj = obj@Model();
		end

		function val = get.nInletPorts(obj)
			val = 1;
		end

		function val = get.nOutletPorts(obj)
			val = 0;
		end
		
		function S = saveobj(obj)
			S = obj.saveobj@Model();
		end

		function res = validate(obj, sectionTimes)
			%VALIDATE Validates the configuration of the Model
			%   RES = VALIDATE(SECTIONTIMES) uses the section times of the simulator given in
			%   SECTIONTIMES to validate the Model object. Returns true in RES if everything
			%   is fine and false otherwise.
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

		function val = getParameterValue(obj, param)
			%GETPARAMETERVALUE Retrieves a parameter value from the model
			%   VAL = GETPARAMETERVALUE(PARAM) searches for the (single) parameter identified by
			%   the struct PARAM with the fields SENS_NAME, SENS_UNIT, SENS_COMP, SENS_PARTYPE,
			%   SENS_BOUNDPHASE, SENS_REACTION, and SENS_SECTION (as returned by MAKESENSITIVITY()).
			%   The returned value VAL contains the current value of the parameter (on the Matlab side,
			%   not in the current CADET configuration) or NaN if the parameter could not be found.
			%
			% See also MODEL.GETPARAMETERVALUE, MEXSIMULATOR.GETPARAMETERVALUE, OUTLETMODEL.SETPARAMETERVALUE,
			%   MAKESENSITIVITY

			val = nan;
		end

		function oldVal = setParameterValue(obj, param, newVal)
			%SETPARAMETERVALUE Sets a parameter value in the model
			%   OLDVAL = SETPARAMETERVALUE(PARAM, NEWVAL) searches for the parameter
			%   identified by the struct PARAM with the fields SENS_NAME, SENS_UNIT,
			%   SENS_COMP, SENS_PARTYPE, SENS_BOUNDPHASE, SENS_REACTION, and SENS_SECTION (as
			%   returned by MAKESENSITIVITY()). The returned value OLDVAL contains the old value
			%   of the parameter (on the Matlab side, not in the current CADET configuration)
			%   or NaN if the parameter could not be found. The value of the parameter is
			%   set to NEWVAL. The changes are not propagated to the underlying CADET simulator.
			%
			% See also MODEL.SETPARAMETERVALUE, MEXSIMULATOR.SETPARAMETERVALUE, OUTLETMODEL.GETPARAMETERVALUE,
			%   MAKESENSITIVITY

			oldVal = nan;
		end

	end

	methods (Access = 'protected')

		function loadobjInternal(obj, S)
			obj.loadobjInternal@Model(S);
		end

	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = OutletModel();
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
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
