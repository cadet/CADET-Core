
classdef Model < handle & matlab.mixin.Heterogeneous
	%Model Base class for unit operation models
	%   Every unit operation model class is derived from this class.
	%
	%   Derived classes are supposed to use the field 'data' for storing
	%   their configuration conforming to the CADET file format spec.

	% Copyright: (C) 2008-2017 The CADET Authors
	%            See the license note at the end of the file.
	
	properties (Hidden, Access = 'protected')
		data; % Struct with stored property values
		hasChangedInternal; % Actually stores the value of the hasChanged property
	end

	properties (Constant, Transient, Abstract)
		name; % Name of the model according to CADET file format specs
		hasInlet; % Determines whether the unit operation has an inlet
		hasOutlet; % Determines whether the unit operation has an outlet
	end

	properties (Dependent, Transient)
		nComponents; % Number of chemical components
		hasChanged; % Determines whether this object has changed after the last synchronization with CADET (C++ layer)
	end

	properties
		unitOpIdx; % Index of this unit operation in a system (0-based)
		returnSolutionInlet; % Determines whether the solution at the inlet is returned
		returnSolutionOutlet; % Determines whether the solution at the outlet is returned
		returnSolutionColumn; % Determines whether the solution in the bulk volume is returned
		returnSolutionParticle; % Determines whether the solution in the particles is returned
		returnSolutionFlux; % Determines whether the solution of the bulk-particle flux is returned
		returnSolDotInlet; % Determines whether the time derivative of the solution at the inlet is returned
		returnSolDotOutlet; % Determines whether the time derivative of the solution at the outlet is returned
		returnSolDotColumn; % Determines whether the time derivative of the solution in the bulk volume is returned
		returnSolDotParticle; % Determines whether the time derivative of the solution in the particles is returned
		returnSolDotFlux; % Determines whether the time derivative of the solution of the bulk-particle flux is returned
		returnSensInlet; % Determines whether the sensitivities at the inlet are returned
		returnSensOutlet; % Determines whether the sensitivities at the outlet are returned
		returnSensColumn; % Determines whether the sensitivities in the bulk volume are returned
		returnSensParticle; % Determines whether the sensitivities in the particles are returned
		returnSensFlux; % Determines whether the sensitivities of the bulk-particle fluxes are returned
		returnSensDotInlet; % Determines whether the time derivatives of the sensitivities at the inlet are returned
		returnSensDotOutlet; % Determines whether the time derivatives of the sensitivities at the outlet are returned
		returnSensDotColumn; % Determines whether the time derivatives of the sensitivities in the bulk volume are returned
		returnSensDotParticle; % Determines whether the time derivatives of the sensitivities in the particles are returned
		returnSensDotFlux; % Determines whether the time derivatives of the sensitivities of the bulk-particle fluxes are returned
	end

	methods

		function obj = Model()
			%MODEL Constructs a model base class object
			%   By default only returns the solution and sensitivity at the outlet.

			obj.hasChangedInternal = true;

			obj.unitOpIdx = -1;
			obj.data = [];
			obj.data.UNIT_TYPE = obj.name;

			obj.returnSolutionInlet = false;
			obj.returnSolutionOutlet = true;
			obj.returnSolutionColumn = false;
			obj.returnSolutionParticle = false;
			obj.returnSolutionFlux = false;

			obj.returnSolDotInlet = false;
			obj.returnSolDotOutlet = false;
			obj.returnSolDotColumn = false;
			obj.returnSolDotParticle = false;
			obj.returnSolDotFlux = false;

			obj.returnSensInlet = false;
			obj.returnSensOutlet = true;
			obj.returnSensColumn = false;
			obj.returnSensParticle = false;
			obj.returnSensFlux = false;

			obj.returnSensDotInlet = false;
			obj.returnSensDotOutlet = false;
			obj.returnSensDotColumn = false;
			obj.returnSensDotParticle = false;
			obj.returnSensDotFlux = false;
		end

		function res = validate(obj, sectionTimes)
			%VALIDATE Validates the configuration of the Model
			%   RES = VALIDATE(SECTIONTIMES) uses the section times of the simulator given in
			%   SECTIONTIMES to validate the Model object. Returns true in RES if everything
			%   is fine and false otherwise.
			%
			%   Derived classes are supposed to extend this function.
			%
			% See also MODELSYSTEM.VALIDATE

			if ~isfield(obj.data, 'NCOMP')
				error('CADET:invalidConfig', 'Property nComponents must be set.');
			end
			validateattributes(obj.nComponents, {'numeric'}, {'scalar', 'nonempty', '>=', 1, 'finite', 'real'}, '', 'nComponents');
			validateattributes(obj.unitOpIdx, {'numeric'}, {'scalar', 'nonempty', '>=', 0, 'finite', 'real'}, '', 'unitOpIdx');

			validateattributes(obj.returnSolutionInlet, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolutionInlet');
			validateattributes(obj.returnSolutionOutlet, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolutionOutlet');
			validateattributes(obj.returnSolutionColumn, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolutionColumn');
			validateattributes(obj.returnSolutionParticle, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolutionParticle');
			validateattributes(obj.returnSolutionFlux, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolutionFlux');

			validateattributes(obj.returnSolDotInlet, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolDotInlet');
			validateattributes(obj.returnSolDotOutlet, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolDotOutlet');
			validateattributes(obj.returnSolDotColumn, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolDotColumn');
			validateattributes(obj.returnSolDotParticle, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolDotParticle');
			validateattributes(obj.returnSolDotFlux, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolDotFlux');

			validateattributes(obj.returnSensInlet, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensInlet');
			validateattributes(obj.returnSensOutlet, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensOutlet');
			validateattributes(obj.returnSensColumn, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensColumn');
			validateattributes(obj.returnSensParticle, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensParticle');
			validateattributes(obj.returnSensFlux, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensFlux');

			validateattributes(obj.returnSensDotInlet, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensDotInlet');
			validateattributes(obj.returnSensDotOutlet, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensDotOutlet');
			validateattributes(obj.returnSensDotColumn, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensDotColumn');
			validateattributes(obj.returnSensDotParticle, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensDotParticle');
			validateattributes(obj.returnSensDotFlux, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensDotFlux');

			res = true;
		end

		function res = assembleConfig(obj)
			%ASSEMBLECONFIG Assembles the configuration according to the CADET file format spec
			%   RES = ASSEMBLECONFIG() returns a nested Matlab struct RES that represents the
			%   Model as detailed in the CADET file format spec.
			%
			%   Models are supposed to store their configuration (conforming to the CADET file
			%   format spec) in the data field of this base class. However, by overwriting this
			%   function, derived classes can customize the assembly process.
			%
			% See also MODELSYSTEM.ASSEMBLECONFIG, MEXSIMULATOR.ASSEMBLECONFIG

			res = obj.data;
		end

		function res = assembleReturnConfig(obj)
			%ASSEMBLERETURNCONFIG Assembles the return configuration according to the CADET file format spec
			%   RES = ASSEMBLERETURNCONFIG() returns a nested Matlab struct RES that contains the return
			%   configuration (return group in the file format spec) of the model.
			%
			% See also MODELSYSTEM.ASSEMBLERETURNCONFIG, MEXSIMULATOR.ASSEMBLERETURNCONFIG

			res = [];

			res.WRITE_SOLUTION_COLUMN_INLET = int32(logical(obj.returnSolutionInlet));
			res.WRITE_SOLUTION_COLUMN_OUTLET = int32(logical(obj.returnSolutionOutlet));
			res.WRITE_SOLUTION_COLUMN = int32(logical(obj.returnSolutionColumn));
			res.WRITE_SOLUTION_PARTICLE = int32(logical(obj.returnSolutionParticle));
			res.WRITE_SOLUTION_FLUX = int32(logical(obj.returnSolutionFlux));

			res.WRITE_SOLDOT_COLUMN_INLET = int32(logical(obj.returnSolDotInlet));
			res.WRITE_SOLDOT_COLUMN_OUTLET = int32(logical(obj.returnSolDotOutlet));
			res.WRITE_SOLDOT_COLUMN = int32(logical(obj.returnSolDotColumn));
			res.WRITE_SOLDOT_PARTICLE = int32(logical(obj.returnSolDotParticle));
			res.WRITE_SOLDOT_FLUX = int32(logical(obj.returnSolDotFlux));

			res.WRITE_SENS_COLUMN_INLET = int32(logical(obj.returnSensInlet));
			res.WRITE_SENS_COLUMN_OUTLET = int32(logical(obj.returnSensOutlet));
			res.WRITE_SENS_COLUMN = int32(logical(obj.returnSensColumn));
			res.WRITE_SENS_PARTICLE = int32(logical(obj.returnSensParticle));
			res.WRITE_SENS_FLUX = int32(logical(obj.returnSensFlux));

			res.WRITE_SENSDOT_COLUMN_INLET = int32(logical(obj.returnSensDotInlet));
			res.WRITE_SENSDOT_COLUMN_OUTLET = int32(logical(obj.returnSensDotOutlet));
			res.WRITE_SENSDOT_COLUMN = int32(logical(obj.returnSensDotColumn));
			res.WRITE_SENSDOT_PARTICLE = int32(logical(obj.returnSensDotParticle));
			res.WRITE_SENSDOT_FLUX = int32(logical(obj.returnSensDotFlux));
		end

		function res = assembleInitialConditions(obj)
			%ASSEMBLEINITIALCONDITIONS Assembles the initial conditions according to the CADET file format spec
			%   RES = ASSEMBLEINITIALCONDITIONS() returns a nested Matlab struct RES that represents only the
			%   initial conditions part of the model as detailed in the (full configuration) CADET file format
			%   spec.
			%
			%   This function is supposed to be overwritten by derived unit operation classes.
			%
			% See also MODELSYSTEM.ASSEMBLEINITIALCONDITIONS

			res = [];
		end

		function notifySync(obj)
			%NOTIFYSYNC Called after parameters have been sent to CADET
			%   NOTIFYSYNC() resets the HASCHANGED property.

			obj.hasChangedInternal = false;
		end

		function S = saveobj(obj)
			S.data = obj.data;
			S.unitOpIdx = obj.unitOpIdx;

			S.returnSolutionInlet = obj.returnSolutionInlet;
			S.returnSolutionOutlet = obj.returnSolutionOutlet;
			S.returnSolutionColumn = obj.returnSolutionColumn;
			S.returnSolutionParticle = obj.returnSolutionParticle;
			S.returnSolutionFlux = obj.returnSolutionFlux;

			S.returnSolDotInlet = obj.returnSolDotInlet;
			S.returnSolDotOutlet = obj.returnSolDotOutlet;
			S.returnSolDotColumn = obj.returnSolDotColumn;
			S.returnSolDotParticle = obj.returnSolDotParticle;
			S.returnSolDotFlux = obj.returnSolDotFlux;

			S.returnSensInlet = obj.returnSensInlet;
			S.returnSensOutlet = obj.returnSensOutlet;
			S.returnSensColumn = obj.returnSensColumn;
			S.returnSensParticle = obj.returnSensParticle;
			S.returnSensFlux = obj.returnSensFlux;

			S.returnSensDotInlet = obj.returnSensDotInlet;
			S.returnSensDotOutlet = obj.returnSensDotOutlet;
			S.returnSensDotColumn = obj.returnSensDotColumn;
			S.returnSensDotParticle = obj.returnSensDotParticle;
			S.returnSensDotFlux = obj.returnSensDotFlux;
		end

		function val = get.nComponents(obj)
			val = double(obj.data.NCOMP);
		end

		function set.nComponents(obj, val)
			validateattributes(val, {'numeric'}, {'scalar', 'nonempty', '>=', 1, 'finite', 'real'}, '', 'nComponents');
			obj.data.NCOMP = int32(val);
			obj.hasChanged = true;
		end

		function val = get.hasChanged(obj)
			val = obj.getHasChanged();
		end

		function set.hasChanged(obj, val)
			obj.hasChangedInternal = val;
		end
	end

	methods (Abstract)
		val = getParameterValue(obj, param)
			%GETPARAMETERVALUE Retrieves a parameter value from the model
			%   VAL = GETPARAMETERVALUE(PARAM) searches for the (single) parameter identified by
			%   the struct PARAM with the fields SENS_NAME, SENS_UNIT, SENS_COMP, SENS_BOUNDPHASE,
			%   SENS_REACTION, and SENS_SECTION (as returned by MAKESENSITIVITY()). The returned
			%   value VAL contains the current value of the parameter (on the Matlab side, not in
			%   the current CADET configuration) or NaN if the parameter could not be found.
			%
			% See also MEXSIMULATOR.GETPARAMETERVALUE, MODEL.SETPARAMETERVALUE, MAKESENSITIVITY

		oldVal = setParameterValue(obj, param, newVal)
			%SETPARAMETERVALUE Sets a parameter value in the model
			%   OLDVAL = SETPARAMETERVALUE(PARAM, NEWVAL) searches for the parameter
			%   identified by the struct PARAM with the fields SENS_NAME, SENS_UNIT,
			%   SENS_COMP, SENS_BOUNDPHASE, SENS_REACTION, and SENS_SECTION (as returned
			%   by MAKESENSITIVITY()). The returned value OLDVAL contains the old value
			%   of the parameter (on the Matlab side, not in the current CADET configuration)
			%   or NaN if the parameter could not be found. The value of the parameter is
			%   set to NEWVAL. The changes are not propagated to the underlying CADET simulator.
			%
			% See also MEXSIMULATOR.SETPARAMETERVALUE, MODEL.GETPARAMETERVALUE, MAKESENSITIVITY
	end

	methods (Access = 'protected')

		function val = getHasChanged(obj)
			val = obj.hasChangedInternal;
		end

		function loadobjInternal(obj, S)
			obj.hasChanged = false;
			obj.data = S.data;
			obj.unitOpIdx = S.unitOpIdx;

			obj.returnSolutionInlet = S.returnSolutionInlet;
			obj.returnSolutionOutlet = S.returnSolutionOutlet;
			obj.returnSolutionColumn = S.returnSolutionColumn;
			obj.returnSolutionParticle = S.returnSolutionParticle;
			obj.returnSolutionFlux = S.returnSolutionFlux;

			obj.returnSolDotInlet = S.returnSolDotInlet;
			obj.returnSolDotOutlet = S.returnSolDotOutlet;
			obj.returnSolDotColumn = S.returnSolDotColumn;
			obj.returnSolDotParticle = S.returnSolDotParticle;
			obj.returnSolDotFlux = S.returnSolDotFlux;

			obj.returnSensInlet = S.returnSensInlet;
			obj.returnSensOutlet = S.returnSensOutlet;
			obj.returnSensColumn = S.returnSensColumn;
			obj.returnSensParticle = S.returnSensParticle;
			obj.returnSensFlux = S.returnSensFlux;

			obj.returnSensDotInlet = S.returnSensDotInlet;
			obj.returnSensDotOutlet = S.returnSensDotOutlet;
			obj.returnSensDotColumn = S.returnSensDotColumn;
			obj.returnSensDotParticle = S.returnSensDotParticle;
			obj.returnSensDotFlux = S.returnSensDotFlux;
		end

	end
end

% =============================================================================
%  CADET - The Chromatography Analysis and Design Toolkit
%  
%  Copyright (C) 2008-2017: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. obj program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies obj distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
