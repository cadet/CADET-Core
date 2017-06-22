
classdef MexSimulator < handle
	%MexSimulator Represents the CADET simulator using a MEX interface
	%   Uses the MEX interface to create and run simulations. A simulator is kept
	%   in memory for the time the instance of this class exists.
	%
	%   A simulation moves through one or multiple sections. A section is a time
	%   span that denotes a certain set of active parameters or a certain phase
	%   in a process simulation (e.g., load, wash, etc.). 
	%
	%   After a simulation has been set up, it can be parameterized. This allows the
	%   simulation to be run multiple times with a (small) changing set of parameters.
	%   This can be useful for optimizations (e.g., inferring binding parameters), or
	%   process analysis (e.g., parameter studies, characterization of operation
	%   window), for instance. These use cases can benefit from a classification 
	%   of changeable model parameters (including SECTION_TIMES) into sensitive and
	%   insensitive (only variable) ones. Whereas for sensitive parameters the 
	%   derivatives of the simulation result with respect to those parameters
	%   (i.e., parameter sensitivities) are computed, variable parameters are only
	%   changed. This allows to speed up use cases that do not require parameter
	%   sensitivties by using only variable parameters.
	%
	%   Model parameters can be joined. If there are two parameters, p1 and p2, 
	%   that have a linear relationship 
	%      a * p1 = b * p2,
	%   with some coefficients a and b, then a joined parameter P is introduced with
	%      P = a * p1 = b * p2.
	%   Assigning a value to P automatically sets p1 and p2 to their respective
	%   values. Sensitivities of the joined parameter P are efficiently computed:
	%   f(p1, p2) = f(P / a, P / b)  =>  df/dP = f_{p1} * 1/a + f_{p2} * 1/b.
	%   This happens automatically in CADET (via automatic differentiation).
	%   Combined parameters can be useful, for instance, when a chain of unit
	%   operations is optimized for flow rate, but the flow rates of subsequent
	%   unit operations are a function of the rate of the first one (e.g., u_1 = 2 * u_0).
	%
	%   Objects of this class can be saved to and loaded from MAT files. The latter
	%   process naturally requires the existence of all model classes that are used.
	%   Note that the restoration of simulator state (e.g., state vectors) may not
	%   be complete.
	%
	%   The MEX interface uses a CADET simulator that is kept in memory. The
	%   classes of the underlying CADET simulator are replicated in Matlab to
	%   some extent. The Matlab layer is used as input for the C++ layer. Since
	%   the parameters of the C++ models never change, this resembles a top-down
	%   architecture without the need of propagating changes back from the C++
	%   to the Matlab layer. However, not all changes in the Matlab layer are
	%   directly communicated to the C++ layer. Consequently, the Matlab layer
	%   can be in a different (newer) state than the C++ layer, which is eventually
	%   synchronized upon starting a simulation. A synchronization can also be
	%   forced (for specific parts of the layer, that is, the simulator, or 
	%   certain models) by calling MEXSIMULATOR.RECONFIGURE.
	
	% Copyright: (C) 2008-2017 The CADET Authors
	%            See the license note at the end of the file.
	
	properties (Hidden, Transient, Access = 'protected')
		mexHandle; % Handle to the underlying CADET simulator in the MEX file
		variableParameters; % Parameters that can be changed but do not have to be sensitive
		sensitivities; % Parameter sensitivities
		sensParamToVariableParam; % Maps indices of sensitivities to indices of corresponding variableParameters
	end

	properties
		model; % Model to be simulated
	end

	properties (Hidden, Access = 'protected')
		data; % Struct with stored property values
	end

	properties (Dependent, Transient)
		% Requested time points at which solution is returned (all time points returned if empty) in [s]
		solutionTimes;
		% Number of threads used by CADET
		nThreads;
		% Absolute error tolerance (scalar or vector with number of DOFs elements)
		absTol;
		% Relative error tolerance
		relTol;
		% Error tolerance of algebraic equations in consistent initialization
		algTol;
		% Relative error tolerance for sensitivity systems
		relTolSens;
		% Initial step size (scalar or vector with value for each section) in [s]
		initStepSize;
		% Maximum number of time integration steps
		maxSteps;
		% Maximum step size of the time integrator
		maxStepSize;
		% Determines how the system is initialized consistently
		consistentInitMode;
		% Determines how the sensitivity systems are initialized consistently
		consistentInitModeSens;
		% Section times [t0, t1, t2, t3, ...] in [s]
		sectionTimes;
		% Determines whether transition from section [t_n, t_{n+1}] to [t_{n+1}, t_{n+2}] is continuous
		sectionContinuity;
		% Determines whether the solution time points are returned from the simulator
		returnSolutionTimes;
		% Determines whether the last state of the system is returned by the simulator
		returnLastState;
		% Determines whether the last states of the sensitivity systems are returned by the simulator
		returnLastStateSens;
		% Returns whether the underlying CADET simulator is configured and ready to run
		isConfigured;
		% Returns the number of variable parameters
		nVariableParameters;
		% Returns the number of sensitive parameters
		nSensitiveParameters;
		% Returns the last duration of a simulation run in seconds
		lastSimulationDuration;
		% Returns the accumulated duration of all simulation runs in seconds
		totalSimulationDuration;
	end

	methods

		function obj = MexSimulator(model)
			%MEXSIMULATOR Creates an object of the simulator class
			%   MEXSIMULATOR() creates a simulator object without a model assigned.
			%   MEXSIMULATOR(MODEL) creates a simulator object with the given MODEL.
			%
			%   In this function, an instance of CADET is created and kept int 
			%   the Simulator object.

			% Initialize CADET
			obj.mexHandle = CadetMex('create');

			% Create data holder
			obj.data = [];
			obj.data.return = [];
			obj.data.solver = [];
			obj.data.solver.time_integrator = [];
			obj.data.solver.sections = [];

			% Set default values
			obj.absTol = 1e-10;
			obj.relTol = 1e-8;
			obj.algTol = 1e-10;
			obj.relTolSens = 1e-8;
			obj.initStepSize = 1e-6;
			obj.maxSteps = 10000;
			obj.maxStepSize = 0.0;
			obj.nThreads = 1;
			obj.consistentInitMode = 1;
			obj.consistentInitModeSens = 1;

			obj.returnSolutionTimes = true;
			obj.returnLastState = false;
			obj.returnLastStateSens = false;

			if nargin >= 1
				obj.model = model;
			end
		end

		function delete(obj)
			%DELETE Destructor
			%   Destroys the CADET instance that is associated with this simulator object.

			% Free the memory of the CADET instance
			CadetMex('destroy', obj.mexHandle);
		end

		function val = get.solutionTimes(obj)
			if isfield(obj.data.solver, 'USER_SOLUTION_TIMES')
				val = obj.data.solver.USER_SOLUTION_TIMES;
			else
				val = [];
			end
		end
		
		function set.solutionTimes(obj, val)			
			if isempty(val)
				if isfield(obj.data.solver.sections, 'USER_SOLUTION_TIMES')
					obj.data.solver = rmfield(obj.data.solver, 'USER_SOLUTION_TIMES');
				end
			else
				validateattributes(val, {'double'}, {'nonnegative', 'vector', 'increasing', 'finite', 'real'}, '', 'solutionTimes');
				obj.data.solver.USER_SOLUTION_TIMES = val;
			end
			if obj.isConfigured
				CadetMex('setsoltimes', obj.mexHandle, val);
			end
		end
		
		function val = get.nThreads(obj)
			val = double(obj.data.solver.NTHREADS);
		end

		function set.nThreads(obj, val)
			validateattributes(val, {'numeric'}, {'>=', 1.0, 'scalar', 'nonempty', 'finite', 'real'}, '', 'nThreads');
			obj.data.solver.NTHREADS = int32(val);
			if obj.isConfigured
				CadetMex('setnumthreads', obj.mexHandle, int32(val));
			end
		end

		function val = get.absTol(obj)
			val = obj.data.solver.time_integrator.ABSTOL;
		end

		function set.absTol(obj, val)
			validateattributes(val, {'double'}, {'nonnegative', 'scalar', 'nonempty', 'finite', 'real'}, '', 'absTol');
			obj.data.solver.time_integrator.ABSTOL = val;
			if obj.isConfigured
				CadetMex('settimeintopts', obj.mexHandle, [], val, [], [], [], [], []);
			end
		end

		function val = get.relTol(obj)
			val = obj.data.solver.time_integrator.RELTOL;
		end

		function set.relTol(obj, val)
			validateattributes(val, {'double'}, {'nonnegative', 'scalar', 'nonempty', 'finite', 'real'}, '', 'relTol');
			obj.data.solver.time_integrator.RELTOL = val;
			if obj.isConfigured
				CadetMex('settimeintopts', obj.mexHandle, val, [], [], [], [], [], []);
			end
		end

		function val = get.relTolSens(obj)
			val = obj.data.solver.time_integrator.RELTOL_SENS;
		end

		function set.relTolSens(obj, val)
			validateattributes(val, {'double'}, {'nonnegative', 'scalar', 'nonempty', 'finite', 'real'}, '', 'relTolSens');
			obj.data.solver.time_integrator.RELTOL_SENS = val;
			if obj.isConfigured
				CadetMex('settimeintopts', obj.mexHandle, [], [], [], [], [], [], val);
			end
		end

		function val = get.algTol(obj)
			val = obj.data.solver.time_integrator.ALGTOL;
		end

		function set.algTol(obj, val)
			validateattributes(val, {'double'}, {'nonnegative', 'scalar', 'nonempty', 'finite', 'real'}, '', 'algTol');
			obj.data.solver.time_integrator.ALGTOL = val;
			if obj.isConfigured
				CadetMex('settimeintopts', obj.mexHandle, [], [], val, [], [], [], []);
			end
		end

		function val = get.initStepSize(obj)
			val = obj.data.solver.time_integrator.INIT_STEP_SIZE;
		end

		function set.initStepSize(obj, val)
			validateattributes(val, {'double'}, {'positive', 'vector', 'nonempty', 'finite', 'real'}, '', 'initStepSize');
			obj.data.solver.time_integrator.INIT_STEP_SIZE = val;
			if obj.isConfigured
				CadetMex('settimeintopts', obj.mexHandle, [], [], [], val, [], [], []);
			end
		end

		function val = get.maxSteps(obj)
			val = double(obj.data.solver.time_integrator.MAX_STEPS);
		end

		function set.maxSteps(obj, val)
			validateattributes(val, {'numeric'}, {'>=', 1.0, 'scalar', 'nonempty', 'finite', 'real'}, '', 'maxSteps');
			obj.data.solver.time_integrator.MAX_STEPS = int32(val);
			if obj.isConfigured
				CadetMex('settimeintopts', obj.mexHandle, [], [], [], [], int32(val), [], []);
			end
		end

		function val = get.maxStepSize(obj)
			val = double(obj.data.solver.time_integrator.MAX_STEP_SIZE);
		end

		function set.maxStepSize(obj, val)
			validateattributes(val, {'double'}, {'>=', 0.0, 'scalar', 'nonempty', 'finite', 'real'}, '', 'maxStepSize');
			obj.data.solver.time_integrator.MAX_STEP_SIZE = val;
			if obj.isConfigured
				CadetMex('settimeintopts', obj.mexHandle, [], [], [], [], [], val, []);
			end
		end

		function val = get.consistentInitMode(obj)
			val = double(obj.data.solver.CONSISTENT_INIT_MODE);
		end

		function set.consistentInitMode(obj, val)
			validateattributes(val, {'numeric'}, {'>=', 0.0, '<=', 7.0 'scalar', 'nonempty', 'finite', 'real'}, '', 'consistentInitMode');
			obj.data.solver.CONSISTENT_INIT_MODE = int32(val);
			if obj.isConfigured
				CadetMex('setconsinitmode', obj.mexHandle, int32(val));
			end
		end

		function val = get.consistentInitModeSens(obj)
			val = double(obj.data.solver.CONSISTENT_INIT_MODE_SENS);
		end

		function set.consistentInitModeSens(obj, val)
			validateattributes(val, {'numeric'}, {'>=', 0.0, '<=', 7.0 'scalar', 'nonempty', 'finite', 'real'}, '', 'consistentInitModeSens');
			obj.data.solver.CONSISTENT_INIT_MODE_SENS = int32(val);
			if obj.isConfigured
				CadetMex('setconsinitmode', obj.mexHandle, [], int32(val));
			end
		end

		function val = get.sectionTimes(obj)
			val = obj.data.solver.sections.SECTION_TIMES;
		end
		
		function set.sectionTimes(obj, val)
			validateattributes(val, {'double'}, {'nonnegative', 'vector', 'increasing', 'nonempty', 'finite', 'real'}, '', 'sectionTimes');
			if numel(val) < 2
				error('CADET:invalidConfig', 'Expected sectionTimes to have at least 2 elements.');
			end
			obj.data.solver.sections.SECTION_TIMES = val;
			if obj.isConfigured
				CadetMex('setsectimes', obj.mexHandle, val);
			end
		end

		function val = get.sectionContinuity(obj)
			if isfield(obj.data.solver.sections, 'SECTION_CONTINUITY')
				val = logical(obj.data.solver.sections.SECTION_CONTINUITY);
			else
				val = [];
			end
		end
		
		function set.sectionContinuity(obj, val)			
			if isempty(val)
				if isfield(obj.data.solver.sections, 'SECTION_CONTINUITY')
					obj.data.solver.sections = rmfield(obj.data.solver.sections, 'SECTION_CONTINUITY');
				end
			else
				validateattributes(val, {'logical'}, {'vector'}, '', 'sectionContinuity');
				obj.data.solver.sections.SECTION_CONTINUITY = int32(logical(val));
			end
			if obj.isConfigured
				CadetMex('setsectimes', obj.mexHandle, obj.sectionTimes, int32(logical(val)));
			end
		end

		function val = get.returnSolutionTimes(obj)
			val = logical(obj.data.return.WRITE_SOLUTION_TIMES);
		end

		function set.returnSolutionTimes(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolutionTimes');
			obj.data.return.WRITE_SOLUTION_TIMES = int32(logical(val));
			CadetMex('setwritetimeandlaststate', obj.mexHandle, int32(logical(val)), [], []);
		end

		function val = get.returnLastState(obj)
			val = logical(obj.data.return.WRITE_SOLUTION_LAST);
		end

		function set.returnLastState(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnLastState');
			obj.data.return.WRITE_SOLUTION_LAST = int32(logical(val));
			CadetMex('setwritetimeandlaststate', obj.mexHandle, [], int32(logical(val)), []);
		end

		function val = get.returnLastStateSens(obj)
			val = logical(obj.data.return.WRITE_SENS_LAST);
		end

		function set.returnLastStateSens(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnLastStateSens');
			obj.data.return.WRITE_SENS_LAST = int32(logical(val));
			CadetMex('setwritetimeandlaststate', obj.mexHandle, [], [], int32(logical(val)));
		end

		function set.model(obj, val)
			if ~isempty(val) && ~(isa(val, 'ModelSystem') || isa(val, 'SingleUnitOpSystem'))
				error('CADET:invalidConfig', 'Expected a valid model.');
			end
			if ~isempty(obj.model)
				CadetMex('clearsim', obj.mexHandle);
			end
			obj.model = val;

			% Configure and validate
			obj.configureCadet(false);
			obj.model.notifySync();
		end

		function val = get.isConfigured(obj)
			val = CadetMex('isconf', obj.mexHandle);
		end

		function val = get.nVariableParameters(obj)
			val = length(obj.variableParameters);
		end

		function val = get.nSensitiveParameters(obj)
			val = length(obj.sensParamToVariableParam);
		end

		function val = get.lastSimulationDuration(obj)
			val = CadetMex('getsimtime', obj.mexHandle);
		end

		function val = get.totalSimulationDuration(obj)
			[~, val] = CadetMex('getsimtime', obj.mexHandle);
		end

		function res = validate(obj)
			%VALIDATE Validates the configuration of the simulator
			%   RES = VALIDATE() returns true if the settings are valid, otherwise false.
			%   Requires a (valid) model.
			%
			% See also MEXSIMULATOR.CONFIGURECADET

			res = false;

			validateattributes(obj.sectionTimes, {'double'}, {'nonnegative', 'vector', 'increasing', 'nonempty', 'finite', 'real'}, '', 'sectionTimes');
			if numel(obj.sectionTimes) < 2
				error('CADET:invalidConfig', 'Expected sectionTimes to have at least 2 elements.');
			end
			if ~isempty(obj.sectionContinuity) && (numel(obj.sectionContinuity) ~= numel(obj.sectionTimes) - 2)
				error('CADET:invalidConfig', 'Expected sectionContinuity to have %d elements.', numel(obj.sectionTimes) - 2);
			end

			if ~isempty(obj.solutionTimes) && ((obj.solutionTimes(1) < obj.sectionTimes(1)) || (obj.solutionTimes(end) > obj.sectionTimes(end)))
				error('CADET:invalidConfig', 'Expected solutionTimes to be in the range of sectionTimes.');
			end

			if isempty(obj.model) || ~(isa(obj.model, 'ModelSystem') || isa(obj.model, 'SingleUnitOpSystem'))
				error('CADET:invalidConfig', 'Expected a valid model.');
			end

			res = obj.model.validate(obj.sectionTimes);
		end

		function S = saveobj(obj)
			S.data = obj.data;
			S.variableParameters = obj.variableParameters;
			S.sensitivities = obj.sensitivities;
			S.sensParamToVariableParam = obj.sensParamToVariableParam;
			S.isConfigured = obj.isConfigured;
			S.model = [];

			if ~isempty(obj.model)
				S.modelClass = class(obj.model);
				S.model = obj.model.saveobj();
			end
		end

		function res = run(obj, skipValidation)
			%RUN Runs a simulation from the beginning
			%   RES = RUN() resets CADET (solver state) and runs the simulation from
			%   the beginning. Returns the results in a nested Matlab struct that 
			%   contains cell arrays (one cell per unit operation) as leaves. The
			%   configuration is validated by default.
			%
			%   RES = RUN(..., SKIPVALIDATION) toggles whether validation of the
			%   simulator configuration is skipped or not by setting SKIPVALIDATION
			%   appropriately.
			%
			% See also MEXSIMULATOR.RUNWITHPARAMETERS, MEXSIMULATOR.RESUME, 
			%   MEXSIMULATOR.RESUMEWITHPARAMETERS

			if (nargin <= 2) || isempty(skipValidation)
				skipValidation = true;
			end
			res = obj.runWithParameters([], skipValidation);
		end

		function [res] = runWithParameters(obj, paramVals, skipValidation)
			%RUNWITHPARAMETERS Runs a simulation with given parameters from the beginning
			%   RES = RUNWITHPARAMETERS(PARAMVALS) resets CADET (solver state), sets the
			%   configured parameters to PARAMVALS, and runs the simulation from the
			%   beginning. Returns the results in a nested Matlab struct that contains
			%   cell arrays (one cell per unit operation) as leaves. The configuration is
			%   validated by default.
			%
			%   RES = RUNWITHPARAMETERS(..., SKIPVALIDATION) toggles whether validation
			%   of the simulator configuration is skipped or not by setting SKIPVALIDATION
			%   appropriately.
			%
			% See also MEXSIMULATOR.RUN, MEXSIMULATOR.RESUME, MEXSIMULATOR.RESUMEWITHPARAMETERS
			
			if (nargin <= 2) || isempty(skipValidation)
				skipValidation = true;
			end

			obj.updateCadetConfig(skipValidation);

			if ~isempty(paramVals)
				obj.setVariableParameterValues(paramVals);
			end
			res = CadetMex('rerun', obj.mexHandle, obj.model.assembleInitialConditions());
			res = ResultsHelper.extract(res, obj.model.numUnitOperations);
		end

		function [res] = resume(obj, skipValidation)
			%RESUME Resumes a simulation from the current state
			%   RES = RESUME() does not reset solver state (e.g., state vectors) and continues
			%   the simulation from the last point on. Returns the results in a nested Matlab
			%   struct that contains cell arrays (one cell per unit operation) as leaves. The
			%   configuration is validated by default.
			%
			%   RES = RESUME(..., SKIPVALIDATION) toggles whether validation of the simulator
			%   configuration is skipped or not by setting SKIPVALIDATION appropriately.
			%
			% See also MEXSIMULATOR.RESUMEWITHPARAMETERS, MEXSIMULATOR.RUN, MEXSIMULATOR.RUNWITHPARAMETERS

			if (nargin <= 2) || isempty(skipValidation)
				skipValidation = true;
			end
			res = obj.resumeWithParameters([], skipValidation);
		end

		function [res] = resumeWithParameters(obj, paramVals, skipValidation)
			%RESUMEWITHPARAMETERS Resumes a simulation from the current state with given parameters
			%   RES = RESUMEWITHPARAMETERS(PARAMVALS) does not reset solver state (e.g.,
			%   state vectors) and continues the simulation from the last point on with the
			%   parameters PARAMVALS. Returns the results in a nested Matlab struct that
			%   contains cell arrays (one cell per unit operation) as leaves.
			%   CADET is configured if it has not been yet, and the configuration is
			%   validated by default.
			%
			%   RES = RESUMEWITHPARAMETERS(..., SKIPVALIDATION) toggles whether validation
			%   of the simulator configuration is skipped or not by setting SKIPVALIDATION
			%   appropriately.
			%
			% See also MEXSIMULATOR.RESUME, MEXSIMULATOR.RUN, MEXSIMULATOR.RUNWITHPARAMETERS

			if (nargin <= 2) || isempty(skipValidation)
				skipValidation = true;
			end

			obj.updateCadetConfig(skipValidation);

			if ~isempty(paramVals)
				obj.setVariableParameterValues(paramVals);
			end
			res = CadetMex('rerun', obj.mexHandle);
			res = ResultsHelper.extract(res, obj.model.numUnitOperations);
		end

		function res = parameterExists(obj, sens)
			%PARAMETEREXISTS Checks whether the given parameters exist in the model / simulator
			%   RES = PARAMETEREXISTS(SENS) checks whether each parameters given in SENS exists
			%   and returns a logical vector indicating the result. SENS is a struct with the 
			%   fields SENS_NAME, SENS_UNIT, SENS_COMP, SENS_BOUNDPHASE, SENS_REACTION, and
			%   SENS_SECTION. The fields are (cell) arrays such that taking one element of every
			%   field identifies a parameter. Such a struct is easily created by MAKESENSITIVITY().
			%   If CADET is configured, the request is relayed to the MEXsimulator. Otherwise, the
			%   existence is checked purely in Matlab using the assigned model.
			%
			% See also MEXSIMULATOR.GETALLPARAMETERVALUES, MEXSIMULATOR.SETPARAMETERS,
			%   MEXSIMULATOR.GETPARAMETERVALUE, MEXSIMULATOR.SETPARAMETERVALUE, MAKESENSITIVITY
			
			if obj.isConfigured
				res = CadetMex('checkpar', obj.mexHandle, sens.SENS_NAME, sens.SENS_UNIT, sens.SENS_COMP, ...
					sens.SENS_BOUNDPHASE, sens.SENS_REACTION, sens.SENS_SECTION);
			else
				res = false(size(sens.SENS_NAME));
				for i = 1:length(sens.SENS_NAME)
					res(i) = ~isnan(obj.getParameterValue(makeSensitivity(sens.SENS_UNIT(i), sens.SENS_NAME{i}, sens.SENS_COMP(i), ...
						sens.SENS_REACTION(i), sens.SENS_BOUNDPHASE(i), sens.SENS_SECTION(i))));
				end
			end
		end

		function [params, vals] = getAllParameters(obj)
			%GETALLPARAMETERS Returns all existing parameters and their current values
			%   PARAMS = GETALLPARAMETERS() returns a struct array with the fields NAMEHASH,
			%   UNIT, COMP, REACTION, BOUNDPHASE, SECTION that consists of all model
			%   parameters in the current CADET instance. The simulator needs to be
			%   configured prior to calling this function.
			%
			%   [PARAMS, VALS] = GETALLPARAMETERS() additionally returns a vector that
			%   contains the current value of each parameters.
			%
			% See also MEXSIMULATOR.PARAMETEREXISTS, MEXSIMULATOR.SETPARAMETERS,
			%   MEXSIMULATOR.GETPARAMETERVALUE, MEXSIMULATOR.SETPARAMETERVALUE

			if nargout == 1
				params = CadetMex('getallpar', obj.mexHandle);
			elseif nargout == 2
				[params, vals] = CadetMex('getallpar', obj.mexHandle);
			end
		end

		function setParameters(obj, sens, isSensitive, autoAbsTol)
			%SETPARAMETERS Sets variable parameters of the simulation that may be sensitive
			%   SETPARAMETERS(SENS) sets the model parameters (including SECTION_TIMES) 
			%   identified by the struct of arrays SENS as variable. SENS is supposed to be
			%   a cell array of structs created by MAKESENSITIVITY(). The structs consist of
			%   the fields SENS_NAME, SENS_COMP, SENS_UNIT, SENS_REACTION, SENS_BOUNDPHASE,
			%   SENS_SECTION, and SENS_FACTOR. Those fields contain (cell) arrays that
			%   together describe a (joined) parameter. The parameters are assumed to be 
			%   sensitive, that is, derivatives of the simulation results with respect to
			%   those parameters are computed. If the field SENS_ABSTOL is present, it is
			%   used for assigning the absolute error tolerance of the sensitivity equations,
			%   which is computed by dividing SENS_ABSTOL by the current parameter value
			%   (if it is non-zero). Otherwise, the absolute tolerance is set to MEXSIMULATOR.ABSTOL.
			%
			%   SETPARAMETERS(..., ISSENSITIVE) uses the logical vector ISSENSITIVE to
			%   determine which variable parameters are sensitive.
			%
			%   SETPARAMETERS(..., ISSENSITIVE, AUTOABSTOL) controls whether the absolute
			%   tolerance of the sensitivity equations is adjusted to the current parameter
			%   value (AUTOABSTOL is true for this sensitivity). If AUTOABSTOL is
			%   disabled (false) for a sensitivity, the raw absolute tolerance (either
			%   SENS_ABSTOL if set, or MEXSIMULATOR.ABSTOL) is passed through directly without
			%   dividing by the current parameter value.
			%
			% See also MEXSIMULATOR.CLEARPARAMETERS, MEXSIMULATOR.SETPARAMETERVALUE, MEXSIMULATOR.GETPARAMETERVALUE,
			%   MEXSIMULATOR.GETVARIABLEPARAMETERVALUE, MEXSIMULATOR.SETVARIABLEPARAMETERVALUE, MAKESENSITIVITY
			
			if (nargin <= 2) || isempty(isSensitive)
				isSensitive = true(length(sens), 1);
			end

			if (nargin <= 3) || isempty(autoAbsTol)
				autoAbsTol = true(length(sens), 1);
			end
			
			validateattributes(sens, {'cell'}, {'vector', 'nonempty'}, '', 'sens');
			validateattributes(isSensitive, {'logical'}, {'vector', 'numel', numel(sens)}, '', 'isSensitive');
			validateattributes(autoAbsTol, {'logical'}, {'vector', 'numel', numel(sens)}, '', 'autoAbsTol');

			if ~isempty(obj.variableParameters)
				CadetMex('clearsens', obj.mexHandle);
			end

			% Check parameters
			obj.variableParameters = cell(length(sens), 1);
			obj.sensitivities = [];
			obj.sensParamToVariableParam = [];
			j = 1;
			for i = 1:length(sens)
				s = sens{i};

				if ~obj.parameterExists(s)
					obj.variableParameters = [];
					obj.sensitivities = [];
					obj.sensParamToVariableParam = [];
					error('CADET:invalidConfig', 'Expected valid parameters, but parameter %d is invalid or does not exist.', i);
				end

				s.autoAbsTol = autoAbsTol(i);
				s.isSensitive = isSensitive(i);

				obj.variableParameters{i} = s;
				if (s.isSensitive)
					obj.sensitivities{j} = sens{i};
					obj.sensParamToVariableParam(j) = i;
					j = j + 1;
				end
			end

			% Enable sensitivity in CADET if the model has already been set up
			if obj.isConfigured
				for i = 1:length(obj.sensParamToVariableParam)
					curPar = obj.variableParameters{obj.sensParamToVariableParam(i)};
					val = obj.getParameterValue(curPar);

					curAbsTol = obj.absTol;

					% Calculate absolute tolerance for sensitivity systems with autoAbsTol enabled
					if isfield(curPar, 'SENS_ABSTOL') && (curPar.SENS_ABSTOL > 0.0)
						curAbsTol = curPar.SENS_ABSTOL;
					end

					if curPar.autoAbsTol && (val ~= 0)
						curAbsTol = curAbsTol / abs(val);
					end

					CadetMex('setsenspar', obj.mexHandle, curPar.SENS_NAME, curPar.SENS_UNIT, curPar.SENS_COMP, curPar.SENS_BOUNDPHASE, ...
						curPar.SENS_REACTION, curPar.SENS_SECTION, curPar.SENS_FACTOR, curAbsTol);
				end
			end
		end
		
		function clearParameters(obj)
			%CLEARPARAMETERS Clears all sensitive and variable parameters
			%   CLEARPARAMETERS() removes all sensitivities and clears the list of
			%   variable parameters.
			%
			% See also MEXSIMULATOR.SETPARAMETERS
			
			obj.variableParameters = [];
			obj.sensitivities = [];
			obj.sensParamToVariableParam = [];

			CadetMex('clearsens', obj.mexHandle);
		end
		
		function clearSensitivities(obj, keepAsVariable)
			%CLEARSENSITIVITIES Disables computation of parameter sensitivities and (optionally) removes sensitive parameters
			%   CLEARSENSITIVITIES() Disables computation of parameter sensitivities but keeps
			%   sensitive parameters as variable parameters.
			%
			%   CLEARSENSITIVITIES(KEEPASVARIABLE) controls whether sensitive parameters are kept
			%   as variable parameters or removed completely.

			if (nargin <= 1) || isempty(keepAsVariable)
				keepAsVariable = true;
			end

			% Remove sensitivities in CADET
			CadetMex('clearsens', obj.mexHandle);

			if ~keepAsVariable
				obj.variableParameters(obj.sensParamToVariableParam) = [];
			end

			obj.sensitivities = [];
			obj.sensParamToVariableParam = [];
		end

		function saveAsHDF5(obj, fileName)
			%SAVEASHDF5 Saves the current setup to a CADET compatible HDF5 file
			%   SAVEASHDF5(FILENAME) saves the current setup (the simulator does not need to
			%   be configured) to the HDF5 file FILENAME. The HDF5 file follows the CADET
			%   file format specifications and can be run by compatiable CADET executables.

			HDF5Tools.struct2hdf(fileName, obj.assembleConfig(), '/', [], []);
		end

		function val = getParameterValue(obj, param)
			%GETPARAMETERVALUE Retrieves parameter values from the simulator or model
			%   VAL = GETPARAMETERVALUE(PARAM) searches for the parameters identified by
			%   the structs in the cell array PARAM with the fields SENS_NAME, SENS_UNIT,
			%   SENS_COMP, SENS_BOUNDPHASE, SENS_REACTION, and SENS_SECTION (as returned
			%   by MAKESENSITIVITY()). The returned value VAL is an array that contains
			%   the current value of the parameter (on the Matlab side, not in the current
			%   CADET configuration) or NaN if a parameter could not be found.
			%
			% See also MEXSIMULATOR.SETPARAMETERVALUE, MAKESENSITIVITY

			validateattributes(param, {'struct', 'cell'}, {'vector', 'nonempty'}, '', 'param');
			if isstruct(param) && numel(param) ~= 1
				error('CADET:funcParamError', 'Expected a single struct or a cell array of structs, but got struct array.');
			end
			
			if isstruct(param)
				param = {param};
			end

			val = nan(numel(param), 1);
			for i = 1:numel(param)
				[curPar, factor] = extractParam(param{i}, 1);
				if iscell(curPar.SENS_NAME)
					curPar.SENS_NAME = curPar.SENS_NAME{1};
				end

				if curPar.SENS_UNIT == -1
					if (strcmp('SECTION_TIMES', curPar.SENS_NAME) && (curPar.SENS_SECTION >= 0) && (curPar.SENS_SECTION < length(obj.sectionTimes)) ...
						&& (curPar.SENS_REACTION == -1) && (curPar.SENS_BOUNDPHASE == -1) && (curPar.SENS_COMP == -1))
						val(i) = obj.sectionTimes(curPar.SENS_SECTION + 1) / factor;
					end
				elseif ~isempty(obj.model)
					val(i) = obj.model.getParameterValue(curPar) / factor;
				end
			end
		end

		function oldVal = setParameterValue(obj, param, newVal, propagateToCadet)
			%SETPARAMETERVALUE Sets parameter values in the simulator or model
			%   OLDVAL = SETPARAMETERVALUE(PARAM, NEWVAL) searches for the parameters
			%   identified by the structs in the cell array PARAM with the fields SENS_NAME,
			%   SENS_UNIT, SENS_COMP, SENS_BOUNDPHASE, SENS_REACTION, and SENS_SECTION (as
			%   returned by MAKESENSITIVITY()). The returned value OLDVAL is an array that
			%   contains the old value of the parameter (on the Matlab side, not in the
			%   current CADET configuration) or NaN if a parameter could not be found.
			%   The values of the parameters are set to the respective elements in the
			%   NEWVAL vector. The changes are propagated to the underlying CADET simulator
			%   if it has already been configured.
			%
			%   OLDVAL = SETPARAMETERVALUE(..., PROPAGATETOCADET) controls whether parameter
			%   changes are sent to the underlying CADET MEXsimulator. A scalar PROPAGATETOCADET
			%   decides for all parameters, a logical vector decides for each single parameter.
			%
			%   Note that parameters are only sent to CADET if it has been configured yet.
			%
			% See also MEXSIMULATOR.GETPARAMETERVALUE, MAKESENSITIVITY

			validateattributes(param, {'struct', 'cell'}, {'vector', 'nonempty'}, '', 'param');
			if isstruct(param) && numel(param) ~= 1
				error('CADET:funcParamError', 'Expected a single struct or a cell array of structs, but got struct array.');
			end
			
			if isstruct(param)
				param = {param};
			end
			validateattributes(newVal, {'double'}, {'vector', 'numel', length(param)}, '', 'newVal');

			if (nargin <= 3) || isempty(propagateToCadet)
				propagateToCadet = true(length(param), 1);
			end
			if numel(propagateToCadet) == 1
				propagateToCadet = repmat(propagateToCadet, length(param), 1);
			end
			if ~obj.isConfigured
				propagateToCadet = false(length(param), 1);
			end
			validateattributes(propagateToCadet, {'logical'}, {'vector', 'numel', length(param)}, '', 'propagateToCadet');

			oldVal = zeros(length(param), 1);
			for i = 1:numel(param)
				for j = 1:numel(param{i}.SENS_UNIT)
					[curPar, factor] = extractParam(param{i}, j);
					if iscell(curPar.SENS_NAME)
						curPar.SENS_NAME = curPar.SENS_NAME{1};
					end

					if curPar.SENS_UNIT == -1
						if (strcmp('SECTION_TIMES', curPar.SENS_NAME) && (curPar.SENS_SECTION >= 0) && (curPar.SENS_SECTION < length(obj.sectionTimes)) ...
								&& (curPar.SENS_REACTION == -1) && (curPar.SENS_BOUNDPHASE == -1) && (curPar.SENS_COMP == -1))
							oldVal(i) = obj.sectionTimes(curPar.SENS_SECTION + 1) / factor;
							obj.sectionTimes(curPar.SENS_SECTION + 1) = newVal(i) * factor;
						end
					elseif ~isempty(obj.model)
						oldVal(i) = obj.model.setParameterValue(curPar, newVal(i) * factor) / factor;
					end

					% Set values in CADET
					if propagateToCadet(i)
						CadetMex('setparval', obj.mexHandle, {curPar.SENS_NAME}, curPar.SENS_UNIT, curPar.SENS_COMP, curPar.SENS_BOUNDPHASE, ...
							curPar.SENS_REACTION, curPar.SENS_SECTION, newVal);
					end
				end
			end
		end

		function val = getVariableParameterValues(obj)
			%GETVARIABLEPARAMETERVALUES Returns the values of the variable parameters
			%   VAL = GETVARIABLEPARAMETERVALUES() returns the values of the variable
			%   parameters as set by MEXSIMULATOR.SETPARAMETERS(). The array VAL contains
			%   the values of the (joined) parameters in the order they were set when
			%   calling MEXSIMULATOR.SETPARAMETERS().
			%
			% See also MEXSIMULATOR.GETPARAMETERVALUE, MEXSIMULATOR.SETPARAMETERS,
			%   MEXSIMULATOR.SETVARIABLEPARAMETERVALUES

			val = obj.getParameterValue(obj.variableParameters);
		end

		function oldVal = setVariableParameterValues(obj, newVals)
			%SETVARIABLEPARAMETERVALUES Sets the values of the variable parameters
			%   OLDVAL = GETVARIABLEPARAMETERVALUES(NEWVALS) assigned the values in the
			%   vector NEWVALS to the variable parameters as set by MEXSIMULATOR.SETPARAMETERS().
			%   The ordering is given by the ordering in which the variable parameters
			%   have been configured in the call to MEXSIMULATOR.SETPARAMETERS(). The array
			%   OLDVAL contains the old values prior to setting NEWVALS.
			%
			%   The changes are made in the Matlab models, but also propagated to CADET.
			%
			% See also MEXSIMULATOR.GETPARAMETERVALUE, MEXSIMULATOR.SETPARAMETERS,
			%   MEXSIMULATOR.GETVARIABLEPARAMETERVALUES

			if isempty(newVals)
				return;
			end

			validateattributes(newVals, {'double'}, {'vector', 'finite', 'real', 'numel', length(obj.variableParameters)}, '', 'newVals');

			oldVal = zeros(length(obj.variableParameters), 1);

			% Allocate memory
			names = cell(length(obj.variableParameters), 1);
			unitOpIds = int32(zeros(length(obj.variableParameters), 1));
			comps = int32(zeros(length(obj.variableParameters), 1));
			boundPhases = int32(zeros(length(obj.variableParameters), 1));
			reactions = int32(zeros(length(obj.variableParameters), 1));
			sections = int32(zeros(length(obj.variableParameters), 1));
			changedVals = zeros(length(obj.variableParameters), 1);

			isSensitive = false(length(obj.variableParameters), 1);
			idx = 1;

			absTols = ones(length(obj.variableParameters), 1) .* obj.absTol;

			for i = 1:length(obj.variableParameters)
				curPar = obj.variableParameters{i};

				isSensitive(i) = any(obj.sensParamToVariableParam == i);
				nJoined = length(curPar.SENS_UNIT);

				% Collect parameters for batch assignment
				if ~isSensitive(i)
					names(idx:idx+nJoined-1) = curPar.SENS_NAME;
					unitOpIds(idx:idx+nJoined-1) = curPar.SENS_UNIT;
					comps(idx:idx+nJoined-1) = curPar.SENS_COMP;
					boundPhases(idx:idx+nJoined-1) = curPar.SENS_BOUNDPHASE;
					reactions(idx:idx+nJoined-1) = curPar.SENS_REACTION;
					sections(idx:idx+nJoined-1) = curPar.SENS_SECTION;
					changedVals(idx:idx+nJoined-1) = curPar.SENS_FACTOR(:) .* newVals(i);

					idx = idx + nJoined;
				end

				% Apply changes in matlab interface, but don't propagate to CADET yet
				oldVal(i) = obj.setParameterValue(curPar, newVals(i), false);

				if isSensitive(i)
					% Calculate absolute tolerance for sensitivity systems with autoAbsTol enabled
					if isfield(curPar, 'SENS_ABSTOL') && (curPar.SENS_ABSTOL > 0.0)
						absTols(i) = curPar.SENS_ABSTOL;
					end

					if curPar.autoAbsTol && (newVals(i) ~= 0)
						absTols(i) = absTols(i) / abs(newVals(i));
					end
				end
			end

			% Apply changed values in CADET (higher efficiency due to batched changes -> only 1 MEX call)
			if any(~isSensitive)
				CadetMex('setparval', obj.mexHandle, names(1:idx-1), unitOpIds(1:idx-1), comps(1:idx-1), boundPhases(1:idx-1), ...
					reactions(1:idx-1), sections(1:idx-1), changedVals(1:idx-1));
			end

			% Apply changed sensitive parameters
			if any(isSensitive)
				CadetMex('setsensparval', obj.mexHandle, newVals(obj.sensParamToVariableParam));

				% Update error tolerances for sensitivities in CADET
				CadetMex('setsenserror', obj.mexHandle, obj.relTolSens, absTols(obj.sensParamToVariableParam));
			end

		end

		function reconfigure(obj, model, conf)
			%RECONFIGURE Forces the simulator to reread all model parameters of the given model or simulator
			%   RECONFIGURE() forces the CADET simulator to reread all parameters of the simulator itself
			%   (e.g., error tolerances) and some global model parameters like SECTION_TIMES.
			%
			%   RECONFIGURE(MODEL) forces the CADET simulator to reread all model parameters of the given
			%   model. MODEL can be an object of a unit operation Model class, an object of a ModelSystem
			%   class, or a unit operation id. If the unit operation id in MODEL is -1, the underlying
			%   ModelSystem is reconfigured. Note that a reconfiguration of a ModelSystem does not recurse
			%   into the Model objects the system contains.
			%
			%   RECONFIGURE(MODEL, CONF) uses the given Matlab struct CONF as source for reconfiguring the
			%   MODEL (ModelSystem, unit operation Model, or unit operation id).
			%
			%   Note that no structural changes can be made here (i.e., unit operations cannot be added or
			%   removed, discretizations cannot be changed, etc.).
			%
			% See also MEXSIMULATOR.CONFIGURECADET, MODEL.ASSEMBLECONFIG

			if (nargin <= 1) || isempty(model)
				if ~obj.validate()
					error('CADET:invalidConfig', 'Expected valid simulator configuration.');					
				end
				CadetMex('reconf', obj.mexHandle, obj.data.solver);

				if isfield(obj.data.solver.sections, 'SECTION_CONTINUITY') && ~isempty(obj.data.solver.sections.SECTION_CONTINUITY)
					CadetMex('setsectimes', obj.mexHandle, obj.data.solver.sections.SECTION_TIMES, obj.data.solver.sections.SECTION_CONTINUITY);
				else
					CadetMex('setsectimes', obj.mexHandle, obj.data.solver.sections.SECTION_TIMES);
				end

				retConfig = [];
				retConfig.input = [];
				retConfig.input.return = obj.assembleReturnConfig();
				CadetMex('setreturnconf', obj.mexHandle, retConfig);

				return;
			end

			if isa(model, 'Model')
				model = int32(model.unitOpIdx);
			elseif isa(model, 'ModelSystem') || isa(model, 'SingleUnitOpSystem')
				model = int32(-1);
			else
				% Assume model is of type double (unit operation id)
				model = int32(model);
			end

			if (nargin <= 2) || isempty(conf)
				if ~model.validate()
					error('CADET:invalidConfig', 'Expected valid model configuration.');					
				end
				conf = model.assembleConfig();
			end

			CadetMex('reconf', obj.mexHandle, conf, model);
		end
	end

	methods (Access = 'protected')

		function res = assembleConfig(obj)
			%ASSEMBLECONFIG Assembles the configuration according to the CADET file format spec
			%   RES = ASSEMBLECONFIG() assembles the configuration for a full simulation. The
			%   configuration is a nested Matlab struct that follows the CADET file format spec
			%   and encompasses the configuration of the simulator, the model system, and the
			%   unit operations. It also configures parameter sensitivities and traverses down
			%   into all class.
			%
			%   Since the resulting Matlab struct is constructed according to the CADET file
			%   format spec, it can be serialized to an HDF5 file and executed with a different
			%   CADET frontend.
			%
			% See also MEXSIMULATOR.ASSEMBLERETURNCONFIG, MEXSIMULATOR.ASSEMBLESENSITIVITIES,
			%   MEXSIMULATOR.CONFIGURECADET

			res = [];
			res.input = obj.data;

			if isempty(obj.model) || ~(isa(obj.model, 'ModelSystem') || isa(obj.model, 'SingleUnitOpSystem'))
				error('CADET:invalidConfig', 'Expected a valid model.');
			end

			res.input.model = obj.model.assembleConfig();
			res.input.return = obj.assembleReturnConfig();

			res.input.sensitivity = obj.assembleSensitivities();
			res.input.sensitivity.NSENS = int32(length(obj.sensitivities));
			res.input.sensitivity.SENS_METHOD = 'ad1';
		end

		function res = assembleReturnConfig(obj)
			%ASSEMBLERETURNCONFIG Assembles the return group of the CADET file format spec
			%   RES = ASSEMBLERETURNCONFIG() assembles the configuration of the returned data
			%   as detailed in the CADET file format spec. The configuration is returned as
			%   nested Matlab struct RES.
			%
			%   This function recurses into the model if there is one. If there is no model,
			%   the return configuration of the simulator is constructed.
			%
			% See also MEXSIMULATOR.ASSEMBLECONFIG

			if isempty(obj.model) || ~(isa(obj.model, 'ModelSystem') || isa(obj.model, 'SingleUnitOpSystem'))
				res = [];
			else
				res = obj.model.assembleReturnConfig();
			end
			res.WRITE_SOLUTION_TIMES = obj.data.return.WRITE_SOLUTION_TIMES;
			res.WRITE_SOLUTION_LAST = obj.data.return.WRITE_SOLUTION_LAST;
			res.WRITE_SENS_LAST = obj.data.return.WRITE_SENS_LAST;
		end

		function res = assembleSensitivities(obj)
			%ASSEMBLESENSITIVITIES Assembles the sensitivities according to the CADET file format spec
			%   RES = ASSEMBLESENSITIVITIES() Assembles the parameter sensitivities as detailed in the
			%   CADET file format spec and returns them as nested Matlab struct in RES.
			%
			% See also MEXSIMULATOR.ASSEMBLECONFIG

			res = [];
			for i = 1:length(obj.sensitivities)
				curSens = obj.sensitivities{i};
				curVarPar = obj.variableParameters{obj.sensParamToVariableParam(i)};

				% Update absolute error tolerance
				if curVarPar.autoAbsTol
					curVal = obj.getParameterValue(curVarPar);

					if curSens.SENS_ABSTOL <= 0.0
						curSens.SENS_ABSTOL = obj.absTol;
					end

					if curVal ~= 0.0
						curSens.SENS_ABSTOL = curSens.SENS_ABSTOL / abs(curVal);
					end
				end

				res.(sprintf('param_%03d', i-1)) = curSens;
			end
		end

		function configureCadet(obj, skipValidation)
			%CONFIGURECADET Configures the underlying CADET simulator by creating all models and submodels
			%   CONFIGURECADET() actually creates and configures the CADET simulator (i.e., time integrator)
			%   and models (including submodels). The full configuration is validated before it is sent
			%   to CADET.
			%
			%   CONFIGURECADET(SKIPVALIDATION) determines whether validation of the full configuration is
			%   performed (false, default) or skipped (true).
			%
			% See also MEXSIMULATOR.ASSEMBLECONFIG, MEXSIMULATOR.UPDATECADETCONFIG

			if (nargin <= 1) || isempty(skipValidation)
				skipValidation = false;
			end

			if ~skipValidation && ~obj.validate()
				error('CADET:invalidConfig', 'Expected valid model and simulator configuration.');
			end

			cadetInput = obj.assembleConfig();
			CadetMex('conf', obj.mexHandle, cadetInput);
		end

		function updateCadetConfig(obj, skipValidation)
			%UPDATECADETCONFIG Updates the underlying CADET simulator by propagating changed configurations in all models and submodels
			%   UPDATECADETCONFIG() checks whether configurations in models or submodels have changed and
			%   updates the modified parts in CADET (C++ layer) after a validation.
			%
			%   UPDATECADETCONFIG(SKIPVALIDATION) determines whether validation of the modified configurations
			%   are performed (false, default) or skipped (true).
			%
			% See also MEXSIMULATOR.CONFIGURECADET

			if (nargin <= 1) || isempty(skipValidation)
				skipValidation = false;
			end

			if isempty(obj.model)
				error('CADET:invalidConfig', 'Expected valid model for simulation.');
			end

			changedUnitOps = obj.model.getChangedUnits();
			for i = 1:length(changedUnitOps)
				unitOpIdx = changedUnitOps(i);

				% Differentiate between ModelSystem and Model objects
				if unitOpIdx == -1
					if ~skipValidation && ~obj.model.validate(obj.sectionTimes)
						error('CADET:invalidConfig', 'Expected valid model system configuration.');
					end

					% Reconfigure only the system itself
					obj.reconfigure(unitOpIdx, obj.model.assembleConfig('system'));
					
					% Reset only the system itself
					obj.model.notifySync(true);
				else
					% SingleUnitOpSystem get a special treatment since they don't have a 'models' property
					if isa(obj.model, 'SingleUnitOpSystem')
						if ~skipValidation && ~obj.model.validate(obj.sectionTimes)
							error('CADET:invalidConfig', 'Expected valid model configuration of unit operation model with id %d.', unitOpIdx);
						end
						obj.reconfigure(unitOpIdx, obj.model.assembleConfig(unitOpIdx));
						obj.model.notifySync();
					else
						m = obj.model.models(unitOpIdx + 1);
						if ~skipValidation && ~m.validate(obj.sectionTimes)
							error('CADET:invalidConfig', 'Expected valid model configuration of unit operation model with id %d.', unitOpIdx);
						end
						obj.reconfigure(unitOpIdx, m.assembleConfig());
						m.notifySync();
					end
				end
			end
		end

		function loadobjInternal(obj, S)
			obj.data = S.data;
			obj.variableParameters = S.variableParameters;
			obj.sensitivities = S.sensitivities;
			obj.sensParamToVariableParam = S.sensParamToVariableParam;

			if ~isempty(S.model)
				ctor = str2func([S.modelClass '.loadobj']);
				obj.model = ctor(S.model);
			end

			% Try to configure the model in CADET
			try
				if ~isempty(obj.model) && S.isConfigured
					obj.configureCadet();
				end
			end
		end
		
	end
	
	methods (Static, Access = 'public')

		function [version, commit, branch] = getVersion()
			%GETVERSION Returns the version of CADET, commit hash and branch it was built from
			%   [VERSION] = getVersion() returns the version string VERSION of the underlying CADET
			%   simulator.
			%
			%   [VERSION, COMMIT] = getVersion() returns version string VERSION and git commit hash
			%   COMMIT.
			%
			%   [VERSION, COMMIT, BRANCH] = getVersion() returns version string VERSION, git commit
			%   hash COMMIT, and branch name BRANCH.
			%
			%   Note that a successful call of this function indicates that the MEX library has been
			%   loaded correctly and CADET is expected to work.
			
			[version, commit, branch] = CadetMex();
		end
		
		function ll = logLevel(newLevel)
			%LOGLEVEL Sets or gets the current log level
			%   LL = LOGLEVEL() returns the current log level.
			%
			%   LL = LOGLEVEL(NEWLEVEL) sets the log level to NEWLEVEL and returns the previous
			%   log level. Valid choices for NEWLEVEL are 'None' (0), 'Fatal' (1), 'Warning' (2),
			%   'Normal' (4), 'Info' (5), 'Debug' (6), and 'Trace' (7). NEWLEVEL can either be
			%   one of the strings or an integer / double (values in the parentheses).
			%
			%   Note that the log level is set for all CADET simulators in the current process,
			%   not just for this particular simulator object.

			if nargin == 0
				ll = CadetMex('loglevel');
			else
				ll = CadetMex('loglevel', newLevel);
			end
		end
		
		function obj = loadobj(S)
			if isstruct(S)
				obj = MexSimulator();
				obj.loadobjInternal(S);
			end
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
