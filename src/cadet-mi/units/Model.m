
classdef Model < handle & matlab.mixin.Heterogeneous
	%Model Base class for unit operation models
	%   Every unit operation model class is derived from this class.
	%
	%   Derived classes are supposed to use the field 'data' for storing
	%   their configuration conforming to the CADET file format spec.

	% Copyright: (C) 2008-2019 The CADET Authors
	%            See the license note at the end of the file.
	
	properties (Hidden, Access = 'protected')
		data; % Struct with stored property values
		hasChangedInternal; % Actually stores the value of the hasChanged property
	end

	properties (Abstract, Constant, Access = 'protected')
		hasConsistencySolver; % Determines whether this unit operation model has a consistency solver
	end

	properties (Constant, Transient, Abstract)
		name; % Name of the model according to CADET file format specs
		hasInlet; % Determines whether the unit operation has an inlet
		hasOutlet; % Determines whether the unit operation has an outlet
	end

	properties (Dependent, Transient, Abstract)
		nInletPorts; % Number of inlet ports
		nOutletPorts; % Number of outlet ports
	end

	properties (Dependent, Transient)
		nComponents; % Number of chemical components
		hasChanged; % Determines whether this object has changed after the last synchronization with CADET (C++ layer)
		solverName; % Name of the consistency solver
		maxIterations; % Maximum number of solver iterations
		initialDamping; % Initial damping of the consistency solver
		minDamping; % Minimum damping of the consistency solver
	end

	properties
		hasChangedReturnConfig; % Stores whether the return configuration has changed
		unitOpIdx; % Index of this unit operation in a system (0-based)
		returnCoordinates; % Determines whether node coordinates are returned
		returnSolutionInlet; % Determines whether the solution at the inlet is returned
		returnSolutionOutlet; % Determines whether the solution at the outlet is returned
		returnSolutionBulk; % Determines whether the solution in the bulk volume is returned
		returnSolutionParticle; % Determines whether the solution in the particle mobile phase is returned
		returnSolutionSolid; % Determines whether the solution in the solid phase is returned
		returnSolutionFlux; % Determines whether the solution of the bulk-particle flux is returned
		returnSolutionVolume; % Determines whether the solution of the volume DOFs is returned
		returnSolDotInlet; % Determines whether the time derivative of the solution at the inlet is returned
		returnSolDotOutlet; % Determines whether the time derivative of the solution at the outlet is returned
		returnSolDotBulk; % Determines whether the time derivative of the solution in the bulk volume is returned
		returnSolDotParticle; % Determines whether the time derivative of the solution in the particle mobile phase is returned
		returnSolDotSolid; % Determines whether the time derivative of the solution in the solid phase is returned
		returnSolDotFlux; % Determines whether the time derivative of the solution of the bulk-particle flux is returned
		returnSolDotVolume; % Determines whether the time derivative of the solution of the volume DOFs is returned
		returnSensInlet; % Determines whether the sensitivities at the inlet are returned
		returnSensOutlet; % Determines whether the sensitivities at the outlet are returned
		returnSensBulk; % Determines whether the sensitivities in the bulk volume are returned
		returnSensParticle; % Determines whether the sensitivities in the particle mobile phase are returned
		returnSensSolid; % Determines whether the sensitivities in the solid phase are returned
		returnSensFlux; % Determines whether the sensitivities of the bulk-particle fluxes are returned
		returnSensVolume; % Determines whether the sensitivities of the volume DOFs are returned
		returnSensDotInlet; % Determines whether the time derivatives of the sensitivities at the inlet are returned
		returnSensDotOutlet; % Determines whether the time derivatives of the sensitivities at the outlet are returned
		returnSensDotBulk; % Determines whether the time derivatives of the sensitivities in the bulk volume are returned
		returnSensDotParticle; % Determines whether the time derivatives of the sensitivities in the particle mobile phase are returned
		returnSensDotSolid; % Determines whether the time derivatives of the sensitivities in the solid phase are returned
		returnSensDotFlux; % Determines whether the time derivatives of the sensitivities of the bulk-particle fluxes are returned
		returnSensDotVolume; % Determines whether the time derivatives of the sensitivities of the volume DOFs are returned
	end

	methods

		function obj = Model()
			%MODEL Constructs a model base class object
			%   By default only returns the solution and sensitivity at the outlet.

			obj.hasChangedInternal = true;
			obj.hasChangedReturnConfig = true;

			obj.unitOpIdx = -1;
			obj.data = [];
			obj.data.UNIT_TYPE = obj.name;
			obj.data.discretization = [];
			obj.data.discretization.consistency_solver = [];

			obj.returnCoordinates = true;

			obj.returnSolutionInlet = false;
			obj.returnSolutionOutlet = true;
			obj.returnSolutionBulk = false;
			obj.returnSolutionParticle = false;
			obj.returnSolutionSolid = false;
			obj.returnSolutionFlux = false;
			obj.returnSolutionVolume = false;

			obj.returnSolDotInlet = false;
			obj.returnSolDotOutlet = false;
			obj.returnSolDotBulk = false;
			obj.returnSolDotParticle = false;
			obj.returnSolDotSolid = false;
			obj.returnSolDotFlux = false;
			obj.returnSolDotVolume = false;

			obj.returnSensInlet = false;
			obj.returnSensOutlet = true;
			obj.returnSensBulk = false;
			obj.returnSensParticle = false;
			obj.returnSensSolid = false;
			obj.returnSensFlux = false;
			obj.returnSensVolume = false;

			obj.returnSensDotInlet = false;
			obj.returnSensDotOutlet = false;
			obj.returnSensDotBulk = false;
			obj.returnSensDotParticle = false;
			obj.returnSensDotSolid = false;
			obj.returnSensDotFlux = false;
			obj.returnSensDotVolume = false;

			obj.solverName = {'ATRN_ERR', 'LEVMAR'};
			obj.maxIterations = 50;
			obj.initialDamping = 1e-2;
			obj.minDamping = 1e-4;
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

			validateattributes(obj.returnCoordinates, {'logical'}, {'scalar', 'nonempty'}, '', 'returnCoordinates');

			validateattributes(obj.returnSolutionInlet, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolutionInlet');
			validateattributes(obj.returnSolutionOutlet, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolutionOutlet');
			validateattributes(obj.returnSolutionBulk, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolutionBulk');
			validateattributes(obj.returnSolutionParticle, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolutionParticle');
			validateattributes(obj.returnSolutionSolid, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolutionSolid');
			validateattributes(obj.returnSolutionFlux, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolutionFlux');
			validateattributes(obj.returnSolutionVolume, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolutionVolume');

			validateattributes(obj.returnSolDotInlet, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolDotInlet');
			validateattributes(obj.returnSolDotOutlet, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolDotOutlet');
			validateattributes(obj.returnSolDotBulk, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolDotBulk');
			validateattributes(obj.returnSolDotParticle, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolDotParticle');
			validateattributes(obj.returnSolDotSolid, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolDotSolid');
			validateattributes(obj.returnSolDotFlux, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolDotFlux');
			validateattributes(obj.returnSolDotVolume, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolDotVolume');

			validateattributes(obj.returnSensInlet, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensInlet');
			validateattributes(obj.returnSensOutlet, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensOutlet');
			validateattributes(obj.returnSensBulk, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensBulk');
			validateattributes(obj.returnSensParticle, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensParticle');
			validateattributes(obj.returnSensSolid, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensSolid');
			validateattributes(obj.returnSensFlux, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensFlux');
			validateattributes(obj.returnSensVolume, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensVolume');

			validateattributes(obj.returnSensDotInlet, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensDotInlet');
			validateattributes(obj.returnSensDotOutlet, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensDotOutlet');
			validateattributes(obj.returnSensDotBulk, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensDotBulk');
			validateattributes(obj.returnSensDotParticle, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensDotParticle');
			validateattributes(obj.returnSensDotSolid, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensDotSolid');
			validateattributes(obj.returnSensDotFlux, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensDotFlux');
			validateattributes(obj.returnSensDotVolume, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensDotVolume');

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

			res.WRITE_COORDINATES = int32(logical(obj.returnCoordinates));

			res.WRITE_SOLUTION_INLET = int32(logical(obj.returnSolutionInlet));
			res.WRITE_SOLUTION_OUTLET = int32(logical(obj.returnSolutionOutlet));
			res.WRITE_SOLUTION_BULK = int32(logical(obj.returnSolutionBulk));
			res.WRITE_SOLUTION_PARTICLE = int32(logical(obj.returnSolutionParticle));
			res.WRITE_SOLUTION_SOLID = int32(logical(obj.returnSolutionSolid));
			res.WRITE_SOLUTION_FLUX = int32(logical(obj.returnSolutionFlux));
			res.WRITE_SOLUTION_VOLUME = int32(logical(obj.returnSolutionVolume));

			res.WRITE_SOLDOT_INLET = int32(logical(obj.returnSolDotInlet));
			res.WRITE_SOLDOT_OUTLET = int32(logical(obj.returnSolDotOutlet));
			res.WRITE_SOLDOT_BULK = int32(logical(obj.returnSolDotBulk));
			res.WRITE_SOLDOT_PARTICLE = int32(logical(obj.returnSolDotParticle));
			res.WRITE_SOLDOT_SOLID = int32(logical(obj.returnSolDotSolid));
			res.WRITE_SOLDOT_FLUX = int32(logical(obj.returnSolDotFlux));
			res.WRITE_SOLDOT_VOLUME = int32(logical(obj.returnSolDotVolume));

			res.WRITE_SENS_INLET = int32(logical(obj.returnSensInlet));
			res.WRITE_SENS_OUTLET = int32(logical(obj.returnSensOutlet));
			res.WRITE_SENS_BULK = int32(logical(obj.returnSensBulk));
			res.WRITE_SENS_PARTICLE = int32(logical(obj.returnSensParticle));
			res.WRITE_SENS_SOLID = int32(logical(obj.returnSensSolid));
			res.WRITE_SENS_FLUX = int32(logical(obj.returnSensFlux));
			res.WRITE_SENS_VOLUME = int32(logical(obj.returnSensVolume));

			res.WRITE_SENSDOT_INLET = int32(logical(obj.returnSensDotInlet));
			res.WRITE_SENSDOT_OUTLET = int32(logical(obj.returnSensDotOutlet));
			res.WRITE_SENSDOT_BULK = int32(logical(obj.returnSensDotBulk));
			res.WRITE_SENSDOT_PARTICLE = int32(logical(obj.returnSensDotParticle));
			res.WRITE_SENSDOT_SOLID = int32(logical(obj.returnSensDotSolid));
			res.WRITE_SENSDOT_FLUX = int32(logical(obj.returnSensDotFlux));
			res.WRITE_SENSDOT_VOLUME = int32(logical(obj.returnSensDotVolume));
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
			obj.hasChangedReturnConfig = false;
		end

		function S = saveobj(obj)
			S.data = obj.data;
			S.unitOpIdx = obj.unitOpIdx;

			S.returnCoordinates = obj.returnCoordinates;

			S.returnSolutionInlet = obj.returnSolutionInlet;
			S.returnSolutionOutlet = obj.returnSolutionOutlet;
			S.returnSolutionBulk = obj.returnSolutionBulk;
			S.returnSolutionParticle = obj.returnSolutionParticle;
			S.returnSolutionSolid = obj.returnSolutionSolid;
			S.returnSolutionFlux = obj.returnSolutionFlux;
			S.returnSolutionVolume = obj.returnSolutionVolume;

			S.returnSolDotInlet = obj.returnSolDotInlet;
			S.returnSolDotOutlet = obj.returnSolDotOutlet;
			S.returnSolDotBulk = obj.returnSolDotBulk;
			S.returnSolDotParticle = obj.returnSolDotParticle;
			S.returnSolDotSolid = obj.returnSolDotSolid;
			S.returnSolDotFlux = obj.returnSolDotFlux;
			S.returnSolDotVolume = obj.returnSolDotVolume;

			S.returnSensInlet = obj.returnSensInlet;
			S.returnSensOutlet = obj.returnSensOutlet;
			S.returnSensBulk = obj.returnSensBulk;
			S.returnSensParticle = obj.returnSensParticle;
			S.returnSensSolid = obj.returnSensSolid;
			S.returnSensFlux = obj.returnSensFlux;
			S.returnSensVolume = obj.returnSensVolume;

			S.returnSensDotInlet = obj.returnSensDotInlet;
			S.returnSensDotOutlet = obj.returnSensDotOutlet;
			S.returnSensDotBulk = obj.returnSensDotBulk;
			S.returnSensDotParticle = obj.returnSensDotParticle;
			S.returnSensDotSolid = obj.returnSensDotSolid;
			S.returnSensDotFlux = obj.returnSensDotFlux;
			S.returnSensDotVolume = obj.returnSensDotVolume;
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

		function set.returnCoordinates(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnCoordinates');
			obj.returnCoordinates = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSolutionInlet(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolutionInlet');
			obj.returnSolutionInlet = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSolutionOutlet(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolutionOutlet');
			obj.returnSolutionOutlet = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSolutionBulk(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolutionBulk');
			obj.returnSolutionBulk = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSolutionParticle(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolutionParticle');
			obj.returnSolutionParticle = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSolutionSolid(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolutionSolid');
			obj.returnSolutionSolid = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSolutionFlux(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolutionFlux');
			obj.returnSolutionFlux = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSolutionVolume(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolutionVolume');
			obj.returnSolutionVolume = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSolDotInlet(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolDotInlet');
			obj.returnSolDotInlet = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSolDotOutlet(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolDotOutlet');
			obj.returnSolDotOutlet = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSolDotBulk(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolDotBulk');
			obj.returnSolDotBulk = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSolDotParticle(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolDotParticle');
			obj.returnSolDotParticle = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSolDotSolid(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolDotSolid');
			obj.returnSolDotSolid = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSolDotFlux(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolDotFlux');
			obj.returnSolDotFlux = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSolDotVolume(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSolDotVolume');
			obj.returnSolDotVolume = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSensInlet(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensInlet');
			obj.returnSensInlet = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSensOutlet(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensOutlet');
			obj.returnSensOutlet = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSensBulk(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensBulk');
			obj.returnSensBulk = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSensParticle(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensParticle');
			obj.returnSensParticle = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSensSolid(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensSolid');
			obj.returnSensSolid = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSensFlux(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensFlux');
			obj.returnSensFlux = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSensVolume(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensVolume');
			obj.returnSensVolume = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSensDotInlet(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensDotInlet');
			obj.returnSensDotInlet = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSensDotOutlet(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensDotOutlet');
			obj.returnSensDotOutlet = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSensDotBulk(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensDotBulk');
			obj.returnSensDotBulk = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSensDotParticle(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensDotParticle');
			obj.returnSensDotParticle = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSensDotSolid(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensDotSolid');
			obj.returnSensDotSolid = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSensDotFlux(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensDotFlux');
			obj.returnSensDotFlux = val;
			obj.hasChangedReturnConfig = true;
		end

		function set.returnSensDotVolume(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'returnSensDotVolume');
			obj.returnSensDotVolume = val;
			obj.hasChangedReturnConfig = true;
		end

		function retChanged = hasReturnConfigurationChanged(obj)
			%HASRETURNCONFIGURATIONCHANGED Determines whether the return configuration has changed since
			%   the last synchronization.
			retChanged = obj.hasChangedReturnConfig;
		end

		function val = get.solverName(obj)
			if isfield(obj.data.discretization.consistency_solver, 'SUBSOLVERS')
				val = obj.data.discretization.consistency_solver.SUBSOLVERS;
			else
				val = obj.data.discretization.consistency_solver.SOLVER_NAME;
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

				obj.data.discretization.consistency_solver.SOLVER_NAME = 'COMPOSITE';
				obj.data.discretization.consistency_solver.SUBSOLVERS = val;
			else
				val = validatestring(val, {'LEVMAR', 'ATRN_RES', 'ATRN_ERR'}, '', 'solverName');
				obj.data.discretization.consistency_solver = rmfield(obj.data.discretization.consistency_solver, 'SUBSOLVERS');
				obj.data.discretization.consistency_solver.SOLVER_NAME = val;
			end
			obj.hasChanged = true;
		end

		function val = get.maxIterations(obj)
			val = double(obj.data.discretization.consistency_solver.MAX_ITERATIONS);
		end

		function set.maxIterations(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real', '>=', 1}, '', 'maxIterations');
			obj.data.discretization.consistency_solver.MAX_ITERATIONS = int32(val);
			obj.hasChanged = true;
		end

		function val = get.initialDamping(obj)
			val = obj.data.discretization.consistency_solver.INIT_DAMPING;
		end

		function set.initialDamping(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real', '>=', 0.0}, '', 'initialDamping');
			obj.data.discretization.consistency_solver.INIT_DAMPING = val;
			obj.hasChanged = true;
		end

		function val = get.minDamping(obj)
			val = obj.data.discretization.consistency_solver.MIN_DAMPING;
		end

		function set.minDamping(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', 'finite', 'real', '>=', 0.0}, '', 'minDamping');
			obj.data.discretization.consistency_solver.MIN_DAMPING = val;
			obj.hasChanged = true;
		end
	end

	methods (Abstract)
		val = getParameterValue(obj, param)
			%GETPARAMETERVALUE Retrieves a parameter value from the model
			%   VAL = GETPARAMETERVALUE(PARAM) searches for the (single) parameter identified by
			%   the struct PARAM with the fields SENS_NAME, SENS_UNIT, SENS_COMP, SENS_PARTYPE,
			%   SENS_BOUNDPHASE, SENS_REACTION, and SENS_SECTION (as returned by MAKESENSITIVITY()).
			%   The returned value VAL contains the current value of the parameter (on the Matlab side,
			%   not in the current CADET configuration) or NaN if the parameter could not be found.
			%
			% See also MEXSIMULATOR.GETPARAMETERVALUE, MODEL.SETPARAMETERVALUE, MAKESENSITIVITY

		oldVal = setParameterValue(obj, param, newVal)
			%SETPARAMETERVALUE Sets a parameter value in the model
			%   OLDVAL = SETPARAMETERVALUE(PARAM, NEWVAL) searches for the parameter
			%   identified by the struct PARAM with the fields SENS_NAME, SENS_UNIT,
			%   SENS_COMP, SENS_PARTYPE, SENS_BOUNDPHASE, SENS_REACTION, and SENS_SECTION (as
			%   returned by MAKESENSITIVITY()). The returned value OLDVAL contains the old value
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

			obj.returnCoordinates = S.returnCoordinates;

			obj.returnSolutionInlet = S.returnSolutionInlet;
			obj.returnSolutionOutlet = S.returnSolutionOutlet;
			obj.returnSolutionBulk = S.returnSolutionBulk;
			obj.returnSolutionParticle = S.returnSolutionParticle;
			obj.returnSolutionSolid = S.returnSolutionSolid;
			obj.returnSolutionFlux = S.returnSolutionFlux;
			obj.returnSolutionVolume = S.returnSolutionVolume;

			obj.returnSolDotInlet = S.returnSolDotInlet;
			obj.returnSolDotOutlet = S.returnSolDotOutlet;
			obj.returnSolDotBulk = S.returnSolDotBulk;
			obj.returnSolDotParticle = S.returnSolDotParticle;
			obj.returnSolDotSolid = S.returnSolDotSolid;
			obj.returnSolDotFlux = S.returnSolDotFlux;
			obj.returnSolDotVolume = S.returnSolDotVolume;

			obj.returnSensInlet = S.returnSensInlet;
			obj.returnSensOutlet = S.returnSensOutlet;
			obj.returnSensBulk = S.returnSensBulk;
			obj.returnSensParticle = S.returnSensParticle;
			obj.returnSensSolid = S.returnSensSolid;
			obj.returnSensFlux = S.returnSensFlux;
			obj.returnSensVolume = S.returnSensVolume;

			obj.returnSensDotInlet = S.returnSensDotInlet;
			obj.returnSensDotOutlet = S.returnSensDotOutlet;
			obj.returnSensDotBulk = S.returnSensDotBulk;
			obj.returnSensDotParticle = S.returnSensDotParticle;
			obj.returnSensDotSolid = S.returnSensDotSolid;
			obj.returnSensDotFlux = S.returnSensDotFlux;
			obj.returnSensDotVolume = S.returnSensDotVolume;
		end

	end
end

% =============================================================================
%  CADET - The Chromatography Analysis and Design Toolkit
%  
%  Copyright (C) 2008-2019: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. obj program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies obj distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
