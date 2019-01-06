
classdef BindingModel < handle & matlab.mixin.Heterogeneous
	%BindingModel Base class for binding models
	%   Provides properties and functionalities common to all binding models.
	%
	%   Derived classes are supposed to use the field 'data' for storing
	%   their configuration conforming to the CADET file format spec.
	%
	%   Changes in the parameters are tracked using the field hasChanged.
	%   If a value is changed using the properties, hasChanged is set to
	%   true. It is reset by calling notifySync(). Derived classes have to
	%   take care of monitoring property changes and updating hasChanged.
	%   Changes made by setParameterValue() are not tracked since they may
	%   have already been propagated to CADET (C++ layer).
	%
	% See also MODEL
	
	% Copyright: (C) 2008-2018 The CADET Authors
	%            See the license note at the end of the file.
	
	properties (Hidden, Access = 'protected')
		data; % Struct with stored property values
	end

	properties (Abstract, Constant, Access = 'protected')
		hasConsistencySolver; % Determines whether this binding model has a consistency solver
	end

	properties (Transient)
		hasChanged; % Determines whether this object has changed after the last synchronization with CADET (C++ layer)
	end

	properties (Constant, Transient, Abstract)
		name; % Name of the binding model according to CADET file format specs
	end

	properties (Dependent, Transient)
		solverName; % Name of the consistency solver
		maxIterations; % Maximum number of solver iterations
		initialDamping; % Initial damping of the consistency solver
		minDamping; % Minimum damping of the consistency solver
	end

	methods

		function obj = BindingModel()
			%BINDINGMODEL Constructs a binding model base class object

			obj.hasChanged = true;
			obj.data.consistency_solver = [];

			% Set some defaults
			obj.solverName = {'ATRN_ERR', 'LEVMAR'};
			obj.maxIterations = 50;
			obj.initialDamping = 1e-2;
			obj.minDamping = 1e-4;
		end

		function res = validate(obj, nComponents, nBoundStates)
			%VALIDATE Validates the binding model parameters
			%   RES = VALIDATE(NCOMPONENTS, NBOUNDSTATES) uses the number of
			%   components NCOMPONENTS of the model and the number of bound
			%   states of each of these components as given by the vector
			%   NBOUNDSTATES to validate the parameters of the binding model.
			%   Returns true in RES if everything is fine and false otherwise.
			%
			%   Derived classes are supposed to overwrite this method.

			res = true;
		end

		function res = assembleConfig(obj)
			%ASSEMBLECONFIG Assembles the configuration according to the CADET file format spec
			%   RES = ASSEMBLECONFIG() returns a Matlab struct RES that represents the
			%   binding model as detailed in the CADET file format spec.
			%
			%   Binding models are supposed to store their configuration (conforming to the
			%   CADET file format spec) in the data field of this base class. However, by
			%   overwriting this function, derived classes can customize the assembly process.
			%
			% See also MODEL.ASSEMBLECONFIG

			res = obj.data;
		end

		function val = getParameterValue(obj, param, nBoundStates)
			%GETPARAMETERVALUE Retrieves a parameter value from the model
			%   VAL = GETPARAMETERVALUE(PARAM, NBOUNDSTATES) searches for the (single) parameter
			%   identified by the struct PARAM with the fields SENS_NAME, SENS_COMP, SENS_PARTYPE,
			%   SENS_BOUNDPHASE, SENS_REACTION, and SENS_SECTION (as returned by MAKESENSITIVITY())
			%   using the number of bound states for each component given in the vector NBOUNDSTATES.
			%   The returned value VAL contains the current value of the parameter (on the Matlab side,
			%   not in the current CADET configuration) or NaN if the parameter could not be found.
			%
			%   Note that the fields SENS_UNIT and SENS_PARTYPE are ignored. Thus, it is assumed that
			%   this object is assigned to the unit operation of the queried parameter.
			%
			%   If the parameter is an element of a vector- or matrix-valued parameter, its offset
			%   in the linearized array is determined by BINDINGMODEL.OFFSETTOPARAMETER, which can
			%   be customized for each binding model.
			%
			% See also BINDINGMODEL.SETPARAMETERVALUE, MAKESENSITIVITY, BINDINGMODEL.OFFSETTOPARAMETER

			if ~isfield(obj.data, param.SENS_NAME)
				% We don't have this parameter
				val = nan;
				return;
			end
			
			offset = obj.offsetToParameter(nBoundStates, param) + 1;
			val = obj.data.(param.SENS_NAME)(offset);
		end

		function oldVal = setParameterValue(obj, param, nBoundStates, newVal)
			%SETPARAMETERVALUE Sets a parameter value in the model
			%   OLDVAL = SETPARAMETERVALUE(PARAM, NBOUNDSTATES, NEWVAL) searches for the
			%   parameter identified by the struct PARAM with the fields SENS_NAME, SENS_COMP,
			%   SENS_PARTYPE, SENS_BOUNDPHASE, SENS_REACTION, and SENS_SECTION (as returned by
			%   MAKESENSITIVITY()) using the number of bound states for each component given
			%   in the vector NBOUNDSTATES. The returned value OLDVAL contains the old value
			%   of the parameter (on the Matlab side, not in the current CADET configuration)
			%   or NaN if the parameter could not be found. The value of the parameter is
			%   set to NEWVAL. The changes are not propagated to the underlying CADET simulator.
			%
			%   Note that the fields SENS_UNIT and SENS_PARTYPE are ignored. Thus, it is assumed
			%   that this object is assigned to the unit operation of the queried parameter.
			%
			%   If the parameter is an element of a vector- or matrix-valued parameter, its offset
			%   in the linearized array is determined by BINDINGMODEL.OFFSETTOPARAMETER, which can
			%   be customized for each binding model.
			%
			% See also BINDINGMODEL.GETPARAMETERVALUE, MAKESENSITIVITY, BINDINGMODEL.OFFSETTOPARAMETER

			if ~isfield(obj.data, param.SENS_NAME)
				% We don't have this parameter
				oldVal = nan;
				return;
			end

			offset = obj.offsetToParameter(nBoundStates, param) + 1;
			oldVal = obj.data.(param.SENS_NAME)(offset);
			obj.data.(param.SENS_NAME)(offset) = newVal;
		end

		function notifySync(obj)
			%NOTIFYSYNC Called after parameters have been sent to CADET
			%   NOTIFYSYNC() resets the HASCHANGED property.

			obj.hasChanged = false;
		end

		function val = get.solverName(obj)
			if isfield(obj.data.consistency_solver, 'SUBSOLVERS')
				val = obj.data.consistency_solver.SUBSOLVERS;
			else
				val = obj.data.consistency_solver.SOLVER_NAME;
			end
		end

		function set.solverName(obj, val)
			% Only accept single string or cell array of strings
			if ~iscell(val) && ~ischar(val)
				error('CADET:invalidConfig', 'Expected solverName to be a string or cell array of strings.');
			end

			% If it is a single cell, treat it as a string
			if iscell(val) && (numel(val) == 1)
				val = val{1};
				validateattributes(val, {'char'}, {}, '', 'solverName');
			end

			if iscell(val)
				for i = 1:length(val)
					validateattributes(val{i}, {'char'}, {}, '', 'solverName');
					val{i} = validatestring(val{i}, {'LEVMAR', 'ATRN_RES', 'ATRN_ERR'}, '', 'solverName');
				end

				obj.data.consistency_solver.SOLVER_NAME = 'COMPOSITE';
				obj.data.consistency_solver.SUBSOLVERS = val;
			else
				val = validatestring(val, {'LEVMAR', 'ATRN_RES', 'ATRN_ERR'}, '', 'solverName');
				obj.data.consistency_solver = rmfield(obj.data.consistency_solver, 'SUBSOLVERS');
				obj.data.consistency_solver.SOLVER_NAME = val;
			end
			obj.hasChanged = true;
		end

		function val = get.maxIterations(obj)
			val = double(obj.data.consistency_solver.MAX_ITERATIONS);
		end

		function set.maxIterations(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real', '>=', 1}, '', 'maxIterations');
			obj.data.consistency_solver.MAX_ITERATIONS = int32(val);
			obj.hasChanged = true;
		end

		function val = get.initialDamping(obj)
			val = obj.data.consistency_solver.INIT_DAMPING;
		end

		function set.initialDamping(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real', '>=', 0.0}, '', 'initialDamping');
			obj.data.consistency_solver.INIT_DAMPING = val;
			obj.hasChanged = true;
		end

		function val = get.minDamping(obj)
			val = obj.data.consistency_solver.MIN_DAMPING;
		end

		function set.minDamping(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real', '>=', 0.0}, '', 'minDamping');
			obj.data.consistency_solver.MIN_DAMPING = val;
			obj.hasChanged = true;
		end

		function S = saveobj(obj)
			S.data = obj.data;
		end

	end

	methods (Access = 'protected')

		function loadobjInternal(obj, S)
			obj.hasChanged = false;
			obj.data = S.data;
		end

		function offset = offsetToParameter(obj, nBoundStates, param)
			%OFFSETTOPARAMETER Computes the (zero-based) offset to the given parameter in a linearized array
			%   OFFSET = OFFSETTOPARAMETER(NBOUNDSTATES, PARAM) uses the number of bound states of each
			%   component given in the vector NBOUNDSTATES and the parameter struct PARAM to compute the
			%   0-based offset of the parameter in a linearized array. The fields SENS_COMP, SENS_BOUNDPHASE,
			%   SENS_REACTION, and SENS_SECTION of PARAM are used to calculate the offset.
			%
			%   This implementation assumes section-boundphase-component-major ordering. Subclasses may change
			%   this default implementation.
			%
			% See also MAKESENSITIVITY

			offset = 0;
			if (param.SENS_COMP ~= -1)
				offset = offset + param.SENS_COMP;
			end
			if (param.SENS_BOUNDPHASE ~= -1)
				offset = offset + param.SENS_BOUNDPHASE * numel(nBoundStates);
			end
			if (param.SENS_SECTION ~= -1)
				offset = offset + param.SENS_SECTION * sum(nBoundStates);
			end
		end

	end
end

% =============================================================================
%  CADET - The Chromatography Analysis and Design Toolkit
%  
%  Copyright (C) 2008-2018: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. obj program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies obj distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
