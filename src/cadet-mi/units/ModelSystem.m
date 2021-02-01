
classdef ModelSystem < handle
	%ModelSystem Represents systems of unit operation models
	%   A ModelSystem contains one or multiple unit operations and their
	%   connections. The connections of the unit operations can change
	%   at the beginning of each section by issuing a valve switch. This
	%   makes complicated setups like SMB (simulated moving bed) or other
	%   multi column processes (e.g., MCSGP - multicolumn countercurrent 
	%   solvent gradient purification) possible.
	%
	%   A connection setup is activated by a valve switch. A valve switch
	%   is indicated by 0-based section index. The connection setup remains
	%   active until the next valve switch occurs. A connection setup
	%   itself is represented by a matrix in which each row encodes a
	%   connection of two unit operations. The first two columns denote the
	%   source and destination unit operation id (0-based), respectively.
	%   The third and fourth column denote the port index (0-based) of the
	%   source and destination port, respectively. The fifth and sixth column
	%   give the component index (0-based) of the source unit that is wired
	%   to the component index (0-based) of the destination unit. 
	%   If both port indices are set to -1, all ports of the two unit
	%   operations are connected to their respective counterpart (requires
	%   the same number of ports). If both component indices are set to -1,
	%   all components of the two unit operations are connected to their
	%   respective counterpart (requires the same number of components).
	%   The connection setups are packed into a cell array to allow for
	%   different matrix sizes in each setup.
	%
	%   The ModelSystem assigns the unit operation IDs to the Model objects
	%   registered in the system. The numbering is sequential and starts
	%   with 0.
	
	% Copyright: (C) 2008-2020 The CADET Authors
	%            See the license note at the end of the file.
	
	properties (Hidden, Access = 'protected')
		data; % Struct with stored property values
	end

	properties
		models; % Array with unit operation models
		externalFunctions; % Array with external functions
		connectionStartSection; % Vector with section indices (0-based) that trigger a connection setup
		connections; % Cell array of connection matrices (7 columns: UOID from, UOID to, port from, port to, comp from, comp to, flow rate)
	end

	properties (Transient, Access = 'protected')
		hasSystemChanged; % Determines whether this system object has changed after the last synchronization with CADET (C++ layer)
	end

	properties (Dependent)
		systemGramSchmidtType; % Type of Gram-Schmidt orthogonalization process
		systemMaxKrylovSize; % Maximum size of the Krylov subspace
		systemMaxRestarts; % Maximum number of restarts in GMRES
		systemSchurSafetyTol; % Schur-complement safety error tolerance
		systemLinearSolutionMode; % Determines whether parallel (1) or sequential (2) solution mode is used for linear systems
	end
	
	properties (Dependent, Transient)
		hasChanged; % Determines whether an object in this system has changed after the last synchronization with CADET (C++ layer)
		initStateY; % Initial state vector of the full system (all DOFs)
		initStateYdot; % Initial time derivative state vector of the full system (all DOFs)
		initSensY; % Initial state of each sensitivity system (columns)
		initSensYdot; % Initial time derivative state of each sensitivity system (columns)
		numUnitOperations; % Number of unit operations in the system
	end

	methods

		function obj = ModelSystem()
			%MODELSYSTEM Constructs an empty model system

			obj.hasSystemChanged = true;

			obj.data = [];
			obj.data.connections = [];
			obj.data.external = [];
			obj.data.solver = [];

			obj.systemLinearSolutionMode = 0; % Automatically detect linear solution mode
			obj.systemGramSchmidtType = 1; % Modified Gram-Schmidt (more stable)
			obj.systemMaxKrylovSize = 0; % Use largest possible size
			obj.systemMaxRestarts = 0;
			obj.systemSchurSafetyTol = 1e-8;
		end

		function res = validate(obj, sectionTimes, subModels)
			%VALIDATE Validates the configuration of the ModelSystem and the assigned Models
			%   RES = VALIDATE(SECTIONTIMES) uses the section times of the simulator given in
			%   SECTIONTIMES to validate the ModelSystem object and all model objects assigned
			%   to this system. Returns true in RES if everything is fine and false otherwise.
			%
			%   RES = VALIDATE(..., SUBMODELS) determines whether the (sub-)models registered
			%   in the system are also validated (SUBMODELS = true, default) or not (SUBMODELS = false).

			if (nargin <= 2) || isempty(subModels)
				subModels = true;
			end

			res = true;

			validateattributes(obj.models, {'Model'}, {'vector', 'nonempty'}, '', 'models');

			if subModels
				for i = 1:length(obj.models)
					res = obj.models(i).validate(sectionTimes) && res;
				end
			end

			for i = 1:length(obj.externalFunctions)
				res = obj.externalFunctions(i).validate(sectionTimes) && res;
			end

			if ~isempty(obj.initStateY)
				validateattributes(obj.initStateY, {'double'}, {'vector', 'real', 'finite'}, '', 'initStateY');
			end
			if ~isempty(obj.initStateYdot)
				validateattributes(obj.initStateYdot, {'double'}, {'vector', 'real', 'finite'}, '', 'initStateYdot');
			end
			if ~isempty(obj.initSensY)
				validateattributes(obj.initSensY, {'double'}, {'2d', 'real', 'finite'}, '', 'initSensY');
			end
			if ~isempty(obj.initSensYdot)
				validateattributes(obj.initSensYdot, {'double'}, {'2d', 'real', 'finite'}, '', 'initSensYdot');
			end

			validateattributes(obj.connectionStartSection, {'numeric'}, {'vector', 'nonempty', 'increasing', '>=', 0, 'finite', 'real'}, '', 'connectionStartSection');
			validateattributes(obj.connections, {'cell'}, {'vector', 'numel', numel(obj.connectionStartSection)}, '', 'connections');
			arrayfun(@(x) validateattributes(obj.connections{x}, {'numeric'}, {'2d', 'finite', 'real', 'size', [NaN, 7], '>=', -1}, '', sprintf('connections{%d}', x)), 1:length(obj.connections));

			for i = 1:length(obj.connections)
				curCon = obj.connections{i};
				for j = 1:size(curCon, 1)
					% Detect invalid flow rate
					if (curCon(j, 7) < 0.0)
						error('CADET:invalidConfig', 'Expected non-negative flow rate in connection %d in valve configuration connections{%d} in the last column.', j, i);
					end

					% Detect invalid mixed 'connect all' operations
					if ((curCon(j, 3) < 0) && (curCon(j, 4) >= 0)) || ((curCon(j, 4) < 0) && (curCon(j, 3) >= 0))
						error('CADET:invalidConfig', 'Expected connection %d in valve configuration connections{%d} to either connect all ports or just one.', j, i);
					end
					if ((curCon(j, 5) < 0) && (curCon(j, 6) >= 0)) || ((curCon(j, 6) < 0) && (curCon(j, 5) >= 0))
						error('CADET:invalidConfig', 'Expected connection %d in valve configuration connections{%d} to either connect all components or just one.', j, i);
					end

					% Find invalid unit operation ids
					if any((curCon(j, 1:2) < 0) | (curCon(j, 1:2) >= length(obj.models)))
						error('CADET:invalidConfig', 'Expected connection %d in valve configuration connections{%d} to use 0-based unit operation ids in the first two columns.', j, i);
					end

					uoFrom = curCon(j, 1) + 1;
					uoTo = curCon(j, 2) + 1;

					% Check if unit operations possess inlets and outlets
					if (~obj.models(uoFrom).hasOutlet)
						error('CADET:invalidConfig', 'Expected source unit operation %d (id) of connection %d of valve configuration connections{%d} to possess an outlet.', curCon(j, 1), j, i);
					end
					if (~obj.models(uoTo).hasInlet)
						error('CADET:invalidConfig', 'Expected destination unit operation %d (id) of connection %d of valve configuration connections{%d} to possess an inlet.', curCon(j, 2), j, i);
					end

					% Find invalid port ids
					if (curCon(j, 3) >= obj.models(uoFrom).nOutletPorts)
						error('CADET:invalidConfig', 'Expected valid 0-based port index in third column of connection %d of valve configuration connections{%d}.', j, i);
					end
					if (curCon(j, 4) >= obj.models(uoTo).nInletPorts)
						error('CADET:invalidConfig', 'Expected valid 0-based port index in fourth column of connection %d of valve configuration connections{%d}.', j, i);
					end

					% Find invalid component ids
					if (curCon(j, 5) >= obj.models(uoFrom).nComponents)
						error('CADET:invalidConfig', 'Expected valid 0-based component index in fifth column of connection %d of valve configuration connections{%d}.', j, i);
					end
					if (curCon(j, 6) >= obj.models(uoTo).nComponents)
						error('CADET:invalidConfig', 'Expected valid 0-based component index in sixth column of connection %d of valve configuration connections{%d}.', j, i);
					end
				end
			end
		end

		function res = assembleConfig(obj, unitOp)
			%ASSEMBLECONFIG Assembles the configuration according to the CADET file format spec
			%   RES = ASSEMBLECONFIG() returns a nested Matlab struct RES that represents the
			%   ModelSystem and its registered Model objects as detailed in the CADET file format
			%   spec.
			%
			%   RES = ASSEMBLECONFIG(UNITOP) returns the configuration of the full system
			%   (UNITOP = 'all' or []), the model system itself (UNITOP = 'system' or -1), or
			%   a specific unit operation identified by its id as detailed in the CADET file
			%   format spec.
			%
			% See also MODEL.ASSEMBLECONFIG, MEXSIMULATOR.ASSEMBLECONFIG

			if (nargin <= 1)
				unitOp = [];
			end

			if ischar(unitOp)
				switch unitOp
					case 'all'
						unitOp = [];
					case 'system'
						unitOp = -1;
					otherwise
						error('CADET:funcParamError', 'Invalid option "%s" for unitOp.', unitOp);
				end
			end

			% Assemble single unit operation only
			if ~isempty(unitOp) && (unitOp >= 0)
				res = obj.models(unitOp + 1).assembleConfig();
				return;
			end

			res = obj.data;

			% Assemble unit operation models
			res.NUNITS = int32(numel(obj.models));
			if isempty(unitOp)
				for i = 1:length(obj.models)
					res.(sprintf('unit_%03d', obj.models(i).unitOpIdx)) = obj.models(i).assembleConfig();
				end
			end

			% Assemble switches
			res.connections.NSWITCHES = int32(length(obj.connectionStartSection));
			res.connections.CONNECTIONS_INCLUDE_PORTS = int32(1);
			for i = 1:length(obj.connectionStartSection)
				curMat = obj.connections{i}.';
				curSwitch = [];
				curSwitch.SECTION = int32(obj.connectionStartSection(i));
				curSwitch.CONNECTIONS = curMat(:);
				res.connections.(sprintf('switch_%03d', i-1)) = curSwitch;
			end

			% Assemble external functions
			for i = 1:length(obj.externalFunctions)
				res.external.(sprintf('source_%03d', i-1)) = obj.externalFunctions(i).assembleConfig();
			end
		end

		function res = assembleReturnConfig(obj)
			%ASSEMBLERETURNCONFIG Assembles the return configuration according to the CADET file format spec
			%   RES = ASSEMBLERETURNCONFIG() returns a nested Matlab struct RES that contains the return
			%   configuration (return group in the file format spec) of the ModelSystem and its registered
			%   Model objects.
			%
			% See also MODEL.ASSEMBLERETURNCONFIG, MEXSIMULATOR.ASSEMBLERETURNCONFIG

			res = [];

			% Assemble unit operation models
			for i = 1:length(obj.models)
				res.(sprintf('unit_%03d', obj.models(i).unitOpIdx)) = obj.models(i).assembleReturnConfig();
			end
		end

		function res = assembleInitialConditions(obj)
			%ASSEMBLEINITIALCONDITIONS Assembles the initial conditions according to the CADET file format spec
			%   RES = ASSEMBLEINITIALCONDITIONS() returns a nested Matlab struct RES that represents only the
			%   initial conditions part of the ModelSystem and each of its registered Model objects as detailed
			%   in the (full configuration) CADET file format spec.
			%
			% See also MODEL.ASSEMBLEINITIALCONDITIONS

			res = obj.data;
			res = rmfield(res, {'connections', 'external'});

			% Assemble unit operation models
			for i = 1:length(obj.models)
				res.(sprintf('unit_%03d', obj.models(i).unitOpIdx)) = obj.models(i).assembleInitialConditions();
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
				for i = 1:length(obj.models)
					obj.models(i).notifySync();
				end
			end
		end

		function S = saveobj(obj)
			S.data = obj.data;
			S.connectionStartSection = obj.connectionStartSection;
			S.connections = obj.connections;
			S.models = cell(numel(obj.models), 1);
			for i = 1:length(obj.models)
				curMod = [];
				curMod.modelClass = class(obj.models(i));
				curMod.model = obj.models(i).saveobj();
				S.models{i} = curMod;
			end

			S.extfuns = cell(numel(obj.externalFunctions), 1);
			for i = 1:length(obj.externalFunctions)
				curExtFun = [];
				curExtFun.extFunClass = class(obj.externalFunctions(i));
				curExtFun.extFun = obj.externalFunctions(i).saveobj();
				S.extfuns{i} = curExtFun;
			end
		end

		function val = get.initStateY(obj)
			if isfield(obj.data, 'INIT_STATE_Y')
				val = obj.data.INIT_STATE_Y;
			else
				val = [];
			end
		end
		
		function set.initStateY(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'INIT_STATE_Y');
			else
				validateattributes(val, {'double'}, {'vector', 'finite', 'real'}, '', 'initStateY');
				obj.data.INIT_STATE_Y = val;
			end

			obj.hasSystemChanged = true;
		end

		function val = get.initStateYdot(obj)
			if isfield(obj.data, 'INIT_STATE_YDOT')
				val = obj.data.INIT_STATE_YDOT;
			else
				val = [];
			end
		end
		
		function set.initStateYdot(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'INIT_STATE_YDOT');
			else
				validateattributes(val, {'double'}, {'vector', 'finite', 'real'}, '', 'initStateYdot');
				obj.data.INIT_STATE_YDOT = val;
			end

			obj.hasSystemChanged = true;
		end

		function val = get.initSensY(obj)
			if isfield(obj.data, 'INIT_STATE_SENSY_000')
				% Extract multiple fields into matrix
				val = MultiFields.read(obj.data, 'INIT_STATE_SENSY');
			else
				val = [];
			end
		end
		
		function set.initSensY(obj, val)
			if ~isempty(val)
				validateattributes(val, {'double'}, {'2d', 'finite', 'real'}, '', 'initSensY');

				% Remove multiple fields
				obj.data = MultiFields.remove(obj.data, 'INIT_STATE_SENSY');
				% Set multiple fields
				obj.data = MultiFields.write(obj.data, 'INIT_STATE_SENSY', val);
			else
				% Remove multiple fields
				obj.data = MultiFields.remove(obj.data, 'INIT_STATE_SENSY');
			end

			obj.hasSystemChanged = true;
		end

		function val = get.initSensYdot(obj)
			if isfield(obj.data, 'INIT_STATE_SENSY_000')
				% Extract multiple fields into matrix
				val = MultiFields.read(obj.data, 'INIT_STATE_SENSYDOT');
			else
				val = [];
			end
		end
		
		function set.initSensYdot(obj, val)
			if ~isempty(val)
				validateattributes(val, {'double'}, {'2d', 'finite', 'real'}, '', 'initSensYdot');

				% Remove multiple fields
				obj.data = MultiFields.remove(obj.data, 'INIT_STATE_SENSYDOT');
				% Set multiple fields
				obj.data = MultiFields.write(obj.data, 'INIT_STATE_SENSYDOT', val);
			else
				% Remove multiple fields
				obj.data = MultiFields.remove(obj.data, 'INIT_STATE_SENSYDOT');
			end

			obj.hasSystemChanged = true;
		end

		function val = get.numUnitOperations(obj)
			val = length(obj.models);
		end
		
		function set.models(obj, val)
			if iscell(val)
				obj.models = Model.empty();
				for i = 1:length(val)
					if ~isa(val{i}, 'Model')
						error('CADET:invalidConfig', 'Expected argument %d to be an object derived from Model.', i);
					end
					obj.models(i) = val{i};
				end
			else
				validateattributes(val, {'Model'}, {'vector'}, '', 'models');
				obj.models = val;
			end

			% Assign unit operation ids
			for i = 1:length(obj.models)
				obj.models(i).unitOpIdx = i - 1;
			end			

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
			val = double(obj.data.solver.GS_TYPE);
		end

		function set.systemGramSchmidtType(obj, val)
			validateattributes(val, {'numeric'}, {'>=', 0, '<=', 1, 'scalar', 'nonempty', 'finite', 'real'}, '', 'gramSchmidtType');
			obj.data.solver.GS_TYPE = int32(val);
			obj.hasSystemChanged = true;
		end

		function val = get.systemMaxKrylovSize(obj)
			val = double(obj.data.solver.MAX_KRYLOV);
		end

		function set.systemMaxKrylovSize(obj, val)
			validateattributes(val, {'numeric'}, {'nonnegative', 'scalar', 'nonempty', 'finite', 'real'}, '', 'maxKrylovSize');
			obj.data.solver.MAX_KRYLOV = int32(val);
			obj.hasSystemChanged = true;
		end

		function val = get.systemMaxRestarts(obj)
			val = double(obj.data.solver.MAX_RESTARTS);
		end

		function set.systemMaxRestarts(obj, val)
			validateattributes(val, {'numeric'}, {'nonnegative', 'scalar', 'nonempty', 'finite', 'real'}, '', 'maxRestarts');
			obj.data.solver.MAX_RESTARTS = int32(val);
			obj.hasSystemChanged = true;
		end

		function val = get.systemSchurSafetyTol(obj)
			val = obj.data.solver.SCHUR_SAFETY;
		end

		function set.systemSchurSafetyTol(obj, val)
			validateattributes(val, {'double'}, {'positive', 'scalar', 'nonempty', 'finite', 'real'}, '', 'schurSafetyTol');
			obj.data.solver.SCHUR_SAFETY = val;
			obj.hasSystemChanged = true;
		end

		function val = get.systemLinearSolutionMode(obj)
			val = double(obj.data.solver.LINEAR_SOLUTION_MODE);
		end

		function set.systemLinearSolutionMode(obj, val)
			validateattributes(val, {'numeric'}, {'>=', 0, '<=', 2, 'nonnegative', 'scalar', 'nonempty', 'finite', 'real'}, '', 'linearSolutionMode');
			obj.data.solver.LINEAR_SOLUTION_MODE = int32(val);
			obj.hasSystemChanged = true;
		end
		
		function val = get.hasSystemChanged(obj)
			val = obj.hasSystemChanged;
			for i = 1:length(obj.externalFunctions)
				val = val || obj.externalFunctions(i).hasChanged;
			end
		end

		function val = get.hasChanged(obj)
			val = obj.hasSystemChanged;
			for i = 1:length(obj.externalFunctions)
				val = val || obj.externalFunctions(i).hasChanged;
			end
			for i = 1:length(obj.models)
				val = val || obj.models(i).hasChanged;
			end
		end		

		function addModel(obj, model)
			%ADDMODEL Adds a Model object to the ModelSystem
			%   ADDMODEL(MODEL) adds the given Model object MODEL to the system.

			validateattributes(model, {'Model'}, {'scalar', 'nonempty'}, '', 'model');

			% Assign unit operation id
			model.unitOpIdx = length(obj.models);
			obj.models(end+1) = model;

			obj.hasSystemChanged = true;
		end

		function val = getParameterValue(obj, param)
			%GETPARAMETERVALUE Retrieves a parameter value from the ModelSystem or one of its registered Model objects
			%   VAL = GETPARAMETERVALUE(PARAM) searches for the (single) parameter identified by
			%   the struct PARAM with the fields SENS_NAME, SENS_UNIT, SENS_COMP, SENS_PARTYPE,
			%   SENS_BOUNDPHASE, SENS_REACTION, and SENS_SECTION (as returned by MAKESENSITIVITY()).
			%   The returned value VAL contains the current value of the parameter (on the Matlab side,
			%   not in the current CADET configuration) or NaN if the parameter could not be found.
			%
			% See also MEXSIMULATOR.GETPARAMETERVALUE, MODELSYSTEM.SETPARAMETERVALUE, MAKESENSITIVITY

			val = nan;
			if (strcmp(param.SENS_NAME, 'CONNECTION') && (param.SENS_UNIT == -1))
				if ((param.SENS_SECTION >= 0) && (param.SENS_SECTION <= max(obj.connectionStartSection)))
					idxSource = param.SENS_BOUNDPHASE;
					idxDest = param.SENS_REACTION;
					portSource = param.SENS_COMP;
					portDest = param.SENS_PARTYPE;
					idxCon = find(param.SENS_SECTION == obj.connectionStartSection, 1, 'first');
					if isempty(idxCon)
						return;
					end

					conn = obj.connections{idxCon};
					if (portSource >= 0)
						idx = find((conn(:, 1) == idxSource) & (conn(:, 2) == idxDest) & (conn(:, 3) == portSource) & (conn(:, 4) == portDest), 1, 'first');
					else
						idx = find((conn(:, 1) == idxSource) & (conn(:, 2) == idxDest), 1, 'first');
					end
					if ~isempty(idx)
						val = conn(idx, 7);
					end
				end
				return;
			end

			for i = 1:length(obj.models)
				if (param.SENS_UNIT == obj.models(i).unitOpIdx)
					val = obj.models(i).getParameterValue(param);
					return;
				end
			end
		end

		function oldVal = setParameterValue(obj, param, newVal)
			%SETPARAMETERVALUE Sets a parameter value in the ModelSystem or one of its registered Model objects
			%   OLDVAL = SETPARAMETERVALUE(PARAM, NEWVAL) searches for the parameter
			%   identified by the struct PARAM with the fields SENS_NAME, SENS_UNIT,
			%   SENS_COMP, SENS_PARTYPE, SENS_BOUNDPHASE, SENS_REACTION, and SENS_SECTION (as
			%   returned by MAKESENSITIVITY()). The returned value OLDVAL contains the old value
			%   of the parameter (on the Matlab side, not in the current CADET configuration)
			%   or NaN if the parameter could not be found. The value of the parameter is
			%   set to NEWVAL. The changes are not propagated to the underlying CADET simulator.
			%
			% See also MEXSIMULATOR.SETPARAMETERVALUE, MODELSYSTEM.GETPARAMETERVALUE, MAKESENSITIVITY

			oldVal = nan;
			if (strcmp(param.SENS_NAME, 'CONNECTION') && (param.SENS_UNIT == -1))
				if ((param.SENS_SECTION >= 0) && (param.SENS_SECTION <= max(obj.connectionStartSection)))
					idxSource = param.SENS_BOUNDPHASE;
					idxDest = param.SENS_REACTION;
					portSource = param.SENS_COMP;
					portDest = param.SENS_PARTYPE;
					idxCon = find(param.SENS_SECTION == obj.connectionStartSection, 1, 'first');
					if isempty(idxCon)
						return;
					end

					conn = obj.connections{idxCon};
					if (portSource >= 0)
						idx = (conn(:, 1) == idxSource) & (conn(:, 2) == idxDest) & (conn(:, 3) == portSource) & (conn(:, 4) == portDest);
					else
						idx = (conn(:, 1) == idxSource) & (conn(:, 2) == idxDest);
					end
					idxFirst = find(idx, 1, 'first');
					if ~isempty(idxFirst)
						oldVal = conn(idx, 7);
						conn(idx, 7) = newVal;
						obj.connections{idxCon} = conn;
					end
				end
				return;
			end

			for i = 1:length(obj.models)
				if (param.SENS_UNIT == obj.models(i).unitOpIdx)
					oldVal = obj.models(i).setParameterValue(param, newVal);
					return;
				end
			end
		end

		function changedUnits = getChangedUnits(obj)
			%GETCHANGEDUNITS Returns unit operation ids of the registered models that have been changed since the last synchronization
			%   CHANGEDUNITS = GETCHANGEDUNITS() returns a vector of (0-based) unit operation ids
			%   of the registered models that have been changed since the last synchronization.
			%   Also returns a -1 id for the ModelSystem itself.

			changedUnits = [];

			for i = 1:length(obj.models)
				if obj.models(i).hasChanged
					changedUnits = [changedUnits; obj.models(i).unitOpIdx];
				end
			end

			if obj.hasSystemChanged
				changedUnits = [changedUnits; -1];
			end
		end

		function retChanged = hasReturnConfigurationChanged(obj)
			%HASRETURNCONFIGURATIONCHANGED Determines whether the return configuration has changed since
			%   the last synchronization.

			for i = 1:length(obj.models)
				if obj.models(i).hasChangedReturnConfig
					retChanged = true;
					return;
				end
			end
			retChanged = false;
		end

	end

	methods (Access = 'protected')

		function loadobjInternal(obj, S)
			obj.hasSystemChanged = false;
			obj.data = S.data;
			obj.connectionStartSection = S.connectionStartSection;
			obj.connections = S.connections;
			tempModels = Model.empty();
			for i = 1:length(S.models)
				ctor = str2func([S.models{i}.modelClass '.loadobj']);
				tempModels(i) = ctor(S.models{i}.model);
			end
			obj.models = tempModels;

			tempExtFuns = ExternalFunction.empty();
			for i = 1:length(S.extfuns)
				ctor = str2func([S.extfuns{i}.extFunClass '.loadobj']);
				tempExtFuns(i) = ctor(S.extfuns{i}.extFun);
			end
			obj.externalFunctions = tempExtFuns;
		end

	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = ModelSystem();
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
