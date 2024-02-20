
classdef SingleUnitOpSystem < handle
	%SingleUnitOpSystem Represents a unit operation system with a single (main) unit operation and an implicit inlet
	%   In order to simulate just a (main) unit operation, an inlet is required.
	%   This class saves the user from creating and wiring the additional inlet
	%   unit operation by providing it implicitly. The class can be used as if
	%   it were a single (main) unit operation with inlet profile.
	%
	%   The main unit operation model is assigned the id 0 and the implicit (but
	%   mandatory) inlet model is assigned id 1.
	%
	%   SingleUnitOpSystem is close to the behavior and settings of the CADET 2.0
	%   Matlab interface. The return configuration of the inlet is ignored.
	%
	% See also MODELSYSTEM

	% Copyright: (C) 2008-2024 The CADET Authors
	%            See the license note at the end of the file.

	properties
		externalFunctions; % Array with external functions
		initStateY; % Initial state vector of the full system (all DOFs)
		initStateYdot; % Initial time derivative state vector of the full system (all DOFs)
		initSensY; % Initial state of each sensitivity system (columns)
		initSensYdot; % Initial time derivative state of each sensitivity system (columns)
		flowRate; % Volumetric flow rate
	end

	properties (Dependent)
		systemGramSchmidtType; % Type of Gram-Schmidt orthogonalization process
		systemMaxKrylovSize; % Maximum size of the Krylov subspace
		systemMaxRestarts; % Maximum number of restarts in GMRES
		systemSchurSafetyTol; % Schur-complement safety error tolerance
		systemLinearSolutionMode; % Determines whether parallel (1) or sequential (2) solution mode is used for linear systems
	end

	properties (Constant)
		numUnitOperations = 2; % Number of unit operations in the system
	end

	properties (Transient, Abstract)
		hasChanged; % Determines whether the main unit operation has changed after the last synchronization with CADET (C++ layer)
	end

	properties (Access = 'protected')
		solverOptions; % Schur complement solver options
	end
	
	properties (Access = 'protected', Transient)
		hasSystemChanged; % Determines whether this object has changed after the last synchronization with CADET (C++ layer)
		hasInletChanged; % Determines whether the inlet has changed after the last synchronization with CADET (C++ layer)
	end

	methods

		function obj = SingleUnitOpSystem()
			%SINGLEUNITOPSYSTEM Constructs a model system with a single (main) unit operation

			obj.hasSystemChanged = true;
			obj.hasInletChanged = true;
			
			obj.flowRate = 1.0;

			obj.systemLinearSolutionMode = 2; % Sequential solution mode
			obj.systemGramSchmidtType = 1; % Modified Gram-Schmidt (more stable)
			obj.systemMaxKrylovSize = 0; % Use largest possible size
			obj.systemMaxRestarts = 0;
			obj.systemSchurSafetyTol = 1e-8;
		end

		function res = validate(obj, sectionTimes, subModels)
			%VALIDATE Validates the configuration of the main unit operation
			%   RES = VALIDATE(SECTIONTIMES) uses the section times of the simulator given in
			%   SECTIONTIMES to validate the main unit operation in this system and the system
			%   itself. Returns true in RES if everything is fine and false otherwise.
			%
			%   RES = VALIDATE(..., SUBMODELS) determines whether the main unit operation and
			%   inlet model are also validated (SUBMODELS = true, default) or not (SUBMODELS = false).
			%
			%   This function is supposed to be extended by derived classes.
			%
			% See also MODELSYSTEM.VALIDATE, MODEL.VALIDATE

			res = true;
			for i = 1:length(obj.externalFunctions)
				res = obj.externalFunctions(i).validate(sectionTimes) && res;
			end

			validateattributes(obj.flowRate, {'double'}, {'>=', 0.0, 'scalar', 'nonempty', 'finite', 'real'}, '', 'flowRate');

			if ~isempty(obj.initStateY)
				validateattributes(obj.initStateY, {'double'}, {'vector', 'finite', 'real'}, '', 'initStateY');
			end
			if ~isempty(obj.initStateYdot)
				validateattributes(obj.initStateYdot, {'double'}, {'vector', 'finite', 'real'}, '', 'initStateYdot');
			end
			if ~isempty(obj.initSensY)
				validateattributes(obj.initSensY, {'double'}, {'2d', 'finite', 'real'}, '', 'initSensY');
			end
			if ~isempty(obj.initSensYdot)
				validateattributes(obj.initSensYdot, {'double'}, {'2d', 'finite', 'real'}, '', 'initSensYdot');
			end
			res = true;
		end

		function S = saveobj(obj)
			S = [];
			S.initStateY = obj.initStateY;
			S.initStateYdot = obj.initStateYdot;
			S.initSensY = obj.initSensY;
			S.initSensYdot = obj.initSensYdot;
			S.flowRate = obj.flowRate;
			S.solverOptions = obj.solverOptions;

			S.extfuns = cell(numel(obj.externalFunctions), 1);
			for i = 1:length(obj.externalFunctions)
				curExtFun = [];
				curExtFun.extFunClass = class(obj.externalFunctions(i));
				curExtFun.extFun = obj.externalFunctions(i).saveobj();
				S.extfuns{i} = curExtFun;
			end
		end

		function set.initStateY(obj, val)
			obj.initStateY = val;
			obj.hasSystemChanged = true;
		end
		
		function set.initStateYdot(obj, val)
			obj.initStateYdot = val;
			obj.hasSystemChanged = true;
		end
		
		function set.initSensY(obj, val)
			obj.initSensY = val;
			obj.hasSystemChanged = true;
		end
		
		function set.initSensYdot(obj, val)
			obj.initSensYdot = val;
			obj.hasSystemChanged = true;
		end

		function set.externalFunctions(obj, val)
			if isempty(val)
				obj.externalFunctions = ExternalFunction.empty();
				obj.hasSystemChanged = true;
				return;
			end

			if iscell(val)
				obj.externalFunctions = ExternalFunction.empty();
				for i = 1:length(val)
					if ~isa(val{i}, 'ExternalFunction')
						error('CADET:invalidConfig', 'Expected argument %d to be an object derived from ExternalFunction.', i);
					end
					obj.externalFunctions(i) = val{i};
				end
			else
				validateattributes(val, {'ExternalFunction'}, {'vector'}, '', 'externalFunctions');
				obj.externalFunctions = val;
			end

			obj.hasSystemChanged = true;
		end

		function val = get.systemGramSchmidtType(obj)
			val = double(obj.solverOptions.GS_TYPE);
		end

		function set.systemGramSchmidtType(obj, val)
			validateattributes(val, {'numeric'}, {'>=', 0, '<=', 1, 'scalar', 'nonempty', 'finite', 'real'}, '', 'gramSchmidtType');
			obj.solverOptions.GS_TYPE = int32(val);
			obj.hasSystemChanged = true;
		end

		function val = get.systemMaxKrylovSize(obj)
			val = double(obj.solverOptions.MAX_KRYLOV);
		end

		function set.systemMaxKrylovSize(obj, val)
			validateattributes(val, {'numeric'}, {'nonnegative', 'scalar', 'nonempty', 'finite', 'real'}, '', 'maxKrylovSize');
			obj.solverOptions.MAX_KRYLOV = int32(val);
			obj.hasSystemChanged = true;
		end

		function val = get.systemMaxRestarts(obj)
			val = double(obj.solverOptions.MAX_RESTARTS);
		end

		function set.systemMaxRestarts(obj, val)
			validateattributes(val, {'numeric'}, {'nonnegative', 'scalar', 'nonempty', 'finite', 'real'}, '', 'maxRestarts');
			obj.solverOptions.MAX_RESTARTS = int32(val);
			obj.hasSystemChanged = true;
		end

		function val = get.systemSchurSafetyTol(obj)
			val = obj.solverOptions.SCHUR_SAFETY;
		end

		function set.systemSchurSafetyTol(obj, val)
			validateattributes(val, {'double'}, {'positive', 'scalar', 'nonempty', 'finite', 'real'}, '', 'schurSafetyTol');
			obj.solverOptions.SCHUR_SAFETY = val;
			obj.hasSystemChanged = true;
		end
		
		function val = get.systemLinearSolutionMode(obj)
			val = double(obj.solverOptions.LINEAR_SOLUTION_MODE);
		end

		function set.systemLinearSolutionMode(obj, val)
			validateattributes(val, {'numeric'}, {'>=', 0, '<=', 2, 'nonnegative', 'scalar', 'nonempty', 'finite', 'real'}, '', 'linearSolutionMode');
			obj.solverOptions.LINEAR_SOLUTION_MODE = int32(val);
			obj.hasSystemChanged = true;
		end

		function set.flowRate(obj, val)
			validateattributes(val, {'double'}, {'>=', 0.0, 'scalar', 'nonempty', 'finite', 'real'}, '', 'flowRate');
			obj.flowRate = val;
			obj.hasSystemChanged = true;
		end

		function val = get.hasSystemChanged(obj)
			val = obj.hasSystemChanged;
			for i = 1:length(obj.externalFunctions)
				val = val || obj.externalFunctions(i).hasChanged;
			end
		end

		function val = getParameterValue(obj, param)
			%GETPARAMETERVALUE Retrieves a parameter value from the model or inlet
			%   VAL = GETPARAMETERVALUE(PARAM) searches for the (single) parameter identified by
			%   the struct PARAM with the fields SENS_NAME, SENS_UNIT, SENS_COMP, SENS_PARTYPE,
			%   SENS_BOUNDPHASE, SENS_REACTION, and SENS_SECTION (as returned by MAKESENSITIVITY()).
			%   The returned value VAL contains the current value of the parameter (on the Matlab side,
			%   not in the current CADET configuration) or NaN if the parameter could not be found.
			%
			% See also MODEL.GETPARAMETERVALUE, MEXSIMULATOR.GETPARAMETERVALUE, SINGLEGRM.SETPARAMETERVALUE,
			%   MAKESENSITIVITY

			val = nan;
			if (strcmp(param.SENS_NAME, 'CONNECTION'))
				if ((param.SENS_SECTION == 0) && (param.SENS_BOUNDPHASE == 1) && (param.SENS_REACTION == 0))
					val = obj.flowRate;
				end
			end
		end

		function oldVal = setParameterValue(obj, param, newVal)
			%SETPARAMETERVALUE Sets a parameter value in the model or inlet
			%   OLDVAL = SETPARAMETERVALUE(PARAM, NEWVAL) searches for the parameter
			%   identified by the struct PARAM with the fields SENS_NAME, SENS_UNIT,
			%   SENS_COMP, SENS_PARTYPE, SENS_BOUNDPHASE, SENS_REACTION, and SENS_SECTION (as
			%   returned by MAKESENSITIVITY()). The returned value OLDVAL contains the old value
			%   of the parameter (on the Matlab side, not in the current CADET configuration)
			%   or NaN if the parameter could not be found. The value of the parameter is
			%   set to NEWVAL. The changes are not propagated to the underlying CADET simulator.
			%
			% See also MODEL.SETPARAMETERVALUE, MEXSIMULATOR.SETPARAMETERVALUE, SINGLEGRM.GETPARAMETERVALUE,
			%   MAKESENSITIVITY

			oldVal = nan;
			if (strcmp(param.SENS_NAME, 'CONNECTION'))
				if ((param.SENS_SECTION == 0) && (param.SENS_BOUNDPHASE == 1) && (param.SENS_REACTION == 0))
					oldVal = obj.flowRate;
					obj.flowRate = newVal;
				end
			end
		end

		function notifySync(obj, systemOnly)
			%NOTIFYSYNC Called after parameters have been sent to CADET
			%   NOTIFYSYNC() resets the HASCHANGED property for the system and all model objects
			%   registered in the system.
			%
			%   NOTIFYSYNC(SYSTEMONLY) resets the HASCHANGED property for the system and, depending
			%   on SYSTEMONLY, for the registered models too (SYSTEMONLY = false, default).

			if (nargin <= 1) || isempty(systemOnly)
				systemOnly = false;
			end

			obj.hasSystemChanged = false;
			for i = 1:length(obj.externalFunctions)
				obj.externalFunctions(i).notifySync();
			end

			if ~systemOnly
				obj.hasChanged = false;
				obj.hasInletChanged = false;
			end
		end

		function changedUnits = getChangedUnits(obj)
			%GETCHANGEDUNITS Returns unit operation ids of the registered models that have been changed since the last synchronization
			%   CHANGEDUNITS = GETCHANGEDUNITS() returns a vector of (0-based) unit operation ids
			%   of the registered models that have been changed since the last synchronization.
			%   Also returns a -1 id for the ModelSystem itself.

			changedUnits = [];

			if obj.hasChanged
				changedUnits = [0];
			end

			if obj.hasInletChanged
				changedUnits = [changedUnits; 1];
			end

			if obj.hasSystemChanged
				changedUnits = [changedUnits; -1];
			end
		end
	end

	methods (Access = 'protected')

		function loadobjInternal(obj, S)
			obj.hasSystemChanged = false;
			obj.hasInletChanged = false;
			obj.hasChanged = false;
			obj.initStateY = S.initStateY;
			obj.initStateYdot = S.initStateYdot;
			obj.initSensY = S.initSensY;
			obj.initSensYdot = S.initSensYdot;
			obj.flowRate = S.flowRate;
			obj.solverOptions = S.solverOptions;

			tempExtFuns = ExternalFunction.empty();
			for i = 1:length(S.extfuns)
				ctor = str2func([S.extfuns{i}.extFunClass '.loadobj']);
				tempExtFuns(i) = ctor(S.extfuns{i}.extFun);
			end
			obj.externalFunctions = tempExtFuns;
		end

		function res = assembleConfigBase(obj, modelConfig, inletConfig)
			%ASSEMBLECONFIGBASE Assembles the configuration of the surrounding system according to the CADET file format spec
			%   RES = ASSEMBLECONFIGBASE(MODELCONFIG, INLETCONFIG) returns a nested Matlab struct RES
			%   that represents the system of unit operations and its main unit operation together
			%   with the implicit inlet model as detailed in the CADET file format spec. The
			%   configuration of the main unit operation is given by the nested Matlab struct
			%   MODELCONFIG and the configuration of the inlet model is given by INLETCONFIG.
			%
			%   This function is supposed to be called by derived classes in their implementation
			%   of MODEL.ASSEMBLECONFIG. This function will then use the ingredients (configuration
			%   of main unit operation and inlet) and combine them to a configuration of a model
			%   system.
			%
			% See also MODEL.ASSEMBLECONFIG, MODELSYSTEM.ASSEMBLECONFIG, MEXSIMULATOR.ASSEMBLECONFIG

			res.connections = [];
			res.external = [];

			% Assemble unit operation models
			res.NUNITS = int32(2);
			if ~isempty(modelConfig)
				res.unit_000 = modelConfig;
			end

			if ~isempty(inletConfig)			
				res.unit_001 = inletConfig;
				res.unit_001.UNIT_TYPE = 'INLET';
				res.unit_001.NCOMP = int32(obj.nComponents);
			end

			% Assemble switches
			conMat = [1, 0, -1, -1, -1, -1, obj.flowRate].';
			res.connections.NSWITCHES = int32(1);
			res.connections.CONNECTIONS_INCLUDE_PORTS = int32(1);
			res.connections.switch_000 = [];
			res.connections.switch_000.SECTION = int32(0);
			res.connections.switch_000.CONNECTIONS = conMat(:);

			res.solver = obj.solverOptions;
			
			if ~isempty(obj.initStateY)
				res.INIT_STATE_Y = obj.initStateY;
			end
			if ~isempty(obj.initStateYdot)
				res.INIT_STATE_YDOT = obj.initStateYdot;
			end

			if ~isempty(obj.initSensY)
				res = MultiFields.write(res, 'INIT_STATE_SENSY', obj.initSensY);
			end
			if ~isempty(obj.initSensYdot)
				res = MultiFields.write(res, 'INIT_STATE_SENSYDOT', obj.initSensYdot);
			end

			% Assemble external functions
			for i = 1:length(obj.externalFunctions)
				res.external.(sprintf('source_%03d', i-1)) = obj.externalFunctions(i).assembleConfig();
			end

		end

		function res = assembleReturnConfigBase(obj, modelConfig)
			%ASSEMBLERETURNCONFIGBASE Assembles the return configuration of the surrounding system according to the CADET file format spec
			%   RES = ASSEMBLERETURNCONFIGBASE(MODELCONFIG) returns a nested Matlab struct RES that contains
			%   the return configuration (return group in the file format spec) of the system, its main unit
			%   operation model, and its inlet model. The configuration of the main unit operation is given
			%   by the nested Matlab struct MODELCONFIG and augmented to yield a return configuration for the
			%   full system of unit operations.
			%
			%   This function is supposed to be called by derived classes in their implementation
			%   of MODEL.ASSEMBLERETURNCONFIG. This function will then use the ingredients (return
			%   configuration of main unit operation) and combine them to a return configuration of
			%   a model system.
			%
			% See also MODEL.ASSEMBLERETURNCONFIG, MODELSYSTEM.ASSEMBLERETURNCONFIG, MEXSIMULATOR.ASSEMBLERETURNCONFIG

			res = [];
			res.unit_000 = modelConfig;
		end

		function res = assembleInitialConditionsBase(obj, modelConfig)
			%ASSEMBLEINITIALCONDITIONSBASE Assembles the initial conditions of the surrounding system according to the CADET file format spec
			%   RES = ASSEMBLEINITIALCONDITIONSBASE() returns a nested Matlab struct RES that represents only
			%   the initial conditions part of the ModelSystem, its main unit operation model, and its inlet 
			%   as detailed in the (full configuration) CADET file format spec. The initial conditions of the
			%   main unit operation model given in MODELCONFIG are augmented to yield an initial condition
			%   of the full system of unit operations.
			%
			%   This function is supposed to be called by derived classes in their implementation
			%   of MODEL.ASSEMBLEINITIALCONDITIONS. This function will then use the ingredients (initial
			%   conditions of main unit operation) and combine them to initial conditions of a model
			%   system.
			%
			% See also MODEL.ASSEMBLEINITIALCONDITIONS, MODELSYSTEM.ASSEMBLEINITIALCONDITIONS

			res = [];
			if ~isempty(obj.initStateY)
				res.INIT_STATE_Y = obj.initStateY;
			end
			if ~isempty(obj.initStateYdot)
				res.INIT_STATE_YDOT = obj.initStateYdot;
			end

			if ~isempty(obj.initSensY)
				res = MultiFields.write(res, 'INIT_STATE_SENSY', obj.initSensY);
			end
			if ~isempty(obj.initSensYdot)
				res = MultiFields.write(res, 'INIT_STATE_SENSYDOT', obj.initSensYdot);
			end

			res.unit_000 = modelConfig;
			res.unit_001 = [];
		end

	end

end

% =============================================================================
%  CADET
%  
%  Copyright (C) 2008-2024: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. obj program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies obj distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
