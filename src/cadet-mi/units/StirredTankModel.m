
classdef StirredTankModel < Model
	%StirredTankModel Represents a continuous stirred (reaction) tank model (short CSTR)
	%   Holds all model specific parameters which are necessary for a CADET simulation.
	%
	% See also MODEL, SINGLECSTR, MODELSYSTEM
	
	% Copyright: (C) 2008-2024 The CADET Authors
	%            See the license note at the end of the file.

	properties
		bindingModel; % Object of the used binding model class
		reactionModelBulk; % Object of the used reaction model class for the bulk volume
		reactionModelParticle; % Object of the used reaction model class for the particle volume
	end

	properties (Dependent)
		% Parameters
		flowRateFilter; % Flow rate of filtered out liquid in [m^3 / s]
		porosity; % Porosity
		nBoundStates; % Number of bound states for each component
		particleTypeVolumeFractions; % Volume fractions of particle types

		% Initial values
		initialConcentration; % Initial concentrations for each component in [mol / m^3]
		initialSolid; % Initial concentrations for each component in the solid phase in [mol / m^3_SP]
		initialVolume; % Initial volume in [m^3]
		initialState; % Initial values for each degree of freedom

		useAnalyticJacobian; % Determines whether Jacobian is calculated analytically or via AD
	end

	properties(Dependent, Transient)
		nInletPorts; % Number of inlet ports
		nOutletPorts; % Number of outlet ports
	end

	properties (Constant)
		name = 'CSTR'; % Type of the model according to CADET file format specs
		hasInlet = true; % Determines whether the unit operation has an inlet
		hasOutlet = true; % Determines whether the unit operation has an outlet
	end

	properties (Constant, Access = 'protected')
		hasConsistencySolver = true; % Determines whether this unit operation model has a consistency solver
	end

	methods
		
		function obj = StirredTankModel()
			%STIRREDTANKMODEL Constructs a StirredTankModel object and inserts as much default values as possible

			obj = obj@Model();

			% Set some default values
			obj.bindingModel = [];
			obj.reactionModelBulk = [];
			obj.reactionModelParticle = [];
			obj.porosity = 1.0;
			obj.flowRateFilter = 0.0;
			obj.useAnalyticJacobian = true;
			obj.particleTypeVolumeFractions = 1;

			% Return volume by default
			obj.returnSolutionVolume = true;
			obj.returnSensVolume = true;
		end
		
		% Parameters

		function val = get.flowRateFilter(obj)
			val = obj.data.FLOWRATE_FILTER;
		end

		function set.flowRateFilter(obj, val)
			validateattributes(val, {'double'}, {'nonnegative', 'vector', 'nonempty', 'finite', 'real'}, '', 'flowRateFilter');
			obj.data.FLOWRATE_FILTER = val(:);
			obj.hasChanged = true;
		end

		function val = get.porosity(obj)
			val = obj.data.POROSITY;
		end

		function set.porosity(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', '>=', 0.0, '<=', 1.0, 'finite', 'real'}, '', 'porosity');
			obj.data.POROSITY = val;
			obj.hasChanged = true;
		end

		function val = get.nBoundStates(obj)
			if isfield(obj.data, 'NBOUND')
				val = double(obj.data.NBOUND);
			else
				val = [];
			end
		end

		function set.nBoundStates(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'NBOUND');
			else
				validateattributes(val, {'numeric'}, {'nonnegative', 'vector', 'nonempty', 'finite', 'real'}, '', 'nBoundStates');
				obj.data.NBOUND = int32(val);
			end
			obj.hasChanged = true;
		end

		function val = get.particleTypeVolumeFractions(obj)
			val = obj.data.PAR_TYPE_VOLFRAC;
		end

		function set.particleTypeVolumeFractions(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, '<=', 1.0, 'finite', 'real'}, '', 'particleTypeVolumeFractions');
			obj.data.PAR_TYPE_VOLFRAC = val;
			obj.hasChanged = true;
		end

		% Initial values

		function val = get.initialConcentration(obj)
			val = obj.data.INIT_C;
		end

		function set.initialConcentration(obj, val)
			validateattributes(val, {'double'}, {'nonnegative', 'vector', 'nonempty', 'finite', 'real'}, '', 'initialConcentration');
			obj.data.INIT_C = val;
			obj.hasChanged = true;
		end

		function val = get.initialSolid(obj)
			if ~isfield(obj.data, 'INIT_Q')
				val = [];
			else
				val = obj.data.INIT_Q;
			end
		end

		function set.initialSolid(obj, val)
			if ~isempty(val)
				validateattributes(val, {'double'}, {'nonnegative', 'vector', 'finite', 'real'}, '', 'initialSolid');
			end
			obj.data.INIT_Q = val;
			obj.hasChanged = true;
		end

		function val = get.initialVolume(obj)
			val = obj.data.INIT_VOLUME;
		end

		function set.initialVolume(obj, val)
			validateattributes(val, {'double'}, {'nonnegative', 'scalar', 'nonempty', 'finite', 'real'}, '', 'initialVolume');
			obj.data.INIT_VOLUME = val;
			obj.hasChanged = true;
		end

		function val = get.initialState(obj)
			if isfield(obj.data, 'INIT_STATE')
				val = obj.data.INIT_STATE;
			else
				val = [];
			end
		end

		function set.initialState(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'INIT_STATE');
			else
				validateattributes(val, {'double'}, {'nonnegative', 'vector', 'nonempty', 'finite', 'real'}, '', 'initialState');
				obj.data.INIT_STATE = val;
			end
			obj.hasChanged = true;
		end

		function val = get.useAnalyticJacobian(obj)
			val = logical(obj.data.USE_ANALYTIC_JACOBIAN);
		end

		function set.useAnalyticJacobian(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'useAnalyticJacobian');
			obj.data.USE_ANALYTIC_JACOBIAN = int32(logical(val));
			obj.hasChanged = true;
		end

		function set.bindingModel(obj, val)
			for i = 1:numel(val)
				if ~isempty(val(i)) && ~isa(val(i), 'BindingModel')
					error('CADET:invalidConfig', 'Expected a valid binding model at index %d.', i);
				end
			end
			obj.bindingModel = val;
			obj.hasChanged = true;
		end

		function set.reactionModelBulk(obj, val)
			if numel(val) > 1
				error('CADET:invalidConfig', 'Expected a single reaction model instead of array.');
			end
			if ~isempty(val) && ~isa(val, 'ReactionModel')
				error('CADET:invalidConfig', 'Expected a valid reaction model.');
			end
			obj.reactionModelBulk = val;
			obj.hasChanged = true;
		end

		function set.reactionModelParticle(obj, val)
			for i = 1:numel(val)
				if ~isempty(val(i)) && ~isa(val(i), 'ReactionModel')
					error('CADET:invalidConfig', 'Expected a valid reaction model at index %d.', i);
				end
			end
			obj.reactionModelParticle = val;
			obj.hasChanged = true;
		end

		function val = get.nInletPorts(obj)
			val = 1;
		end

		function val = get.nOutletPorts(obj)
			val = 1;
		end


		function S = saveobj(obj)
			S = obj.saveobj@Model();

			S.bindingModel = [];
			S.bindingModelClass = [];

			if ~isempty(obj.bindingModel)
				S.bindingModelClass = cell(numel(obj.bindingModel), 1);
				S.bindingModel = cell(numel(obj.bindingModel), 1);
				for i = 1:numel(obj.bindingModel)
					S.bindingModelClass{i} = class(obj.bindingModel(i));
					S.bindingModel{i} = obj.bindingModel(i).saveobj();
				end
			end

			S.reactionModelBulk = [];
			S.reactionModelBulkClass = [];

			if ~isempty(obj.reactionModelBulk)
				S.reactionModelBulkClass = class(obj.reactionModelBulk);
				S.reactionModelBulk = obj.reactionModelBulk.saveobj();
			end

			S.reactionModelParticle = [];
			S.reactionModelParticleClass = [];

			if ~isempty(obj.reactionModelParticle)
				S.reactionModelParticleClass = cell(numel(obj.reactionModelParticle), 1);
				S.reactionModelParticle = cell(numel(obj.reactionModelParticle), 1);
				for i = 1:numel(obj.reactionModelParticle)
					S.reactionModelParticleClass{i} = class(obj.reactionModelParticle(i));
					S.reactionModelParticle{i} = obj.reactionModelParticle(i).saveobj();
				end
			end
		end


		function res = validate(obj, sectionTimes)
			%VALIDATE Validates the configuration of the Model
			%   RES = VALIDATE(SECTIONTIMES) uses the section times of the simulator given in
			%   SECTIONTIMES to validate the Model object. Returns true in RES if everything
			%   is fine and false otherwise.
			%
			% See also MODEL.VALIDATE

			res = obj.validate@Model(sectionTimes);

			if ~isfield(obj.data, 'POROSITY')
				error('CADET:invalidConfig', 'Property porosity must be set.');
			end

			if isfield(obj.data, 'NBOUND')
				validateattributes(obj.nBoundStates, {'numeric'}, {'nonnegative', 'vector', 'numel', obj.nComponents, 'finite', 'real'}, '', 'nBoundStates');
				nTotalBnd = sum(obj.nBoundStates);
			else
				nTotalBnd = 0;
			end
			
			validateattributes(obj.initialConcentration, {'double'}, {'nonnegative', 'vector', 'numel', obj.nComponents, 'finite', 'real'}, '', 'initialConcentration');
			validateattributes(obj.initialVolume, {'double'}, {'nonnegative', 'scalar', 'finite', 'real'}, '', 'initialVolume');
			if ~isempty(obj.initialSolid)
				validateattributes(obj.initialSolid, {'double'}, {'nonnegative', 'vector', 'numel', nTotalBnd, 'finite', 'real'}, '', 'initialSolid');
			end
			if ~isempty(obj.initialState)
				nDof = nComponents + 1 + nTotalBnd;
				if (numel(obj.initialState) ~= nDof) && (numel(obj.initialState) ~= 2 * nDof)
					error('CADET:invalidConfig', 'Expected initialState to be of size %d or %d.', nDof, 2*nDof);
				end
				validateattributes(obj.initialState(1:nDof), {'double'}, {'nonnegative', 'finite', 'real'}, '', 'initialState');
			end

			if (numel(obj.flowRateFilter) ~= 1) && (numel(obj.flowRateFilter) ~= nSections)
				error('CADET:invalidConfig', 'Expected flowRateFilter to be of size %d or %d (number of time sections).', 1, nSections);
			end

			validateattributes(obj.porosity, {'double'}, {'scalar', 'nonempty', '>=', 0.0, '<=', 1.0, 'finite', 'real'}, '', 'porosity');
			if ~isempty(obj.bindingModel)
				validateattributes(obj.particleTypeVolumeFractions, {'double'}, {'vector', 'numel', length(obj.nBoundStates) / obj.nComponents, 'nonempty', '>=', 0.0, '<=', 1.0, 'finite', 'real'}, '', 'particleTypeVolumeFractions');
				if abs(sum(obj.particleTypeVolumeFractions) - 1.0) >= 1e-10
					error('CADET:invalidConfig', 'Expected particleTypeVolumeFractions to sum to 1.0.');
				end
			end

			for i = 1:numel(obj.bindingModel)
				if ~isempty(obj.bindingModel(i)) && ~isa(obj.bindingModel(i), 'BindingModel')
					error('CADET:invalidConfig', 'Expected a valid binding model.');
				end
				if isempty(obj.bindingModel(i)) && (sum(obj.nBoundStates(((i-1) * obj.nComponents + 1):(i * obj.nComponents))) ~= 0)
					error('CADET:invalidConfig', 'Expected no bound states when using no binding model.');
				end

				if ~isempty(obj.bindingModel(i))
					res = obj.bindingModel(i).validate(obj.nComponents, obj.nBoundStates(((i-1) * obj.nComponents + 1):(i * obj.nComponents))) && res;
				end
			end

			if ~isempty(obj.reactionModelBulk)
				if ~isa(obj.reactionModelBulk, 'ReactionModel')
					error('CADET:invalidConfig', 'Expected a valid reaction model in bulk volume.');
				end

				res = obj.reactionModelBulk.validate(obj.nComponents, obj.nBoundStates(((i-1) * obj.nComponents + 1):(i * obj.nComponents))) && res;

				if ~obj.reactionModelBulk.hasBulkPhaseReactions
					warning('CADET:unexptectedConfig', 'Bulk reaction model does not contain bulk reactions.');
				end
			end

			for i = 1:numel(obj.reactionModelParticle)
				if ~isempty(obj.reactionModelParticle(i)) && ~isa(obj.reactionModelParticle(i), 'ReactionModel')
					error('CADET:invalidConfig', 'Expected a valid reaction model in particle volume.');
				end

				if ~isempty(obj.reactionModelParticle(i))
					res = obj.reactionModelParticle(i).validate(obj.nComponents, obj.nBoundStates(((i-1) * obj.nComponents + 1):(i * obj.nComponents))) && res;

					if ~obj.reactionModelParticle(i).hasLiquidPhaseReactions && ~obj.reactionModelParticle(i).hasSolidPhaseReactions
						warning('CADET:unexptectedConfig', 'Particle reaction model %d does neither contain liquid nor solid phase reactions.', i);
					end
				end
			end
		end

		function res = assembleConfig(obj)
			%ASSEMBLECONFIG Assembles the configuration according to the CADET file format spec
			%   RES = ASSEMBLECONFIG() returns a nested Matlab struct RES that represents the
			%   model as detailed in the CADET file format spec.
			%
			% See also MODEL.ASSEMBLECONFIG, MEXSIMULATOR.ASSEMBLECONFIG

			res = obj.assembleConfig@Model();

			if isempty(obj.bindingModel)
				res.ADSORPTION_MODEL = 'NONE';
			else
				res.ADSORPTION_MODEL = cell(numel(obj.bindingModel), 1);
				for i = 1:length(obj.bindingModel)
					if isempty(obj.bindingModel(i))
						res.ADSORPTION_MODEL{i} = 'NONE';
					else
						res.ADSORPTION_MODEL{i} = obj.bindingModel(i).name;
						res.(sprintf('adsorption_%03d', i-1)) = obj.bindingModel(i).assembleConfig();
					end
				end
			end

			if isempty(obj.reactionModelBulk)
				res.REACTION_MODEL = 'NONE';
			else
				res.REACTION_MODEL = obj.reactionModelBulk.name;
				res.reaction_bulk = obj.reactionModelBulk.assembleConfig();
			end

			if isempty(obj.reactionModelParticle)
				res.REACTION_MODEL_PARTICLES = 'NONE';
			else
				res.REACTION_MODEL_PARTICLES = cell(numel(obj.reactionModelParticle), 1);
				for i = 1:length(obj.reactionModelParticle)
					if isempty(obj.reactionModelParticle(i))
						res.REACTION_MODEL_PARTICLES{i} = 'NONE';
					else
						res.REACTION_MODEL_PARTICLES{i} = obj.reactionModelParticle(i).name;
						res.(sprintf('reaction_particle_%03d', i-1)) = obj.reactionModelParticle(i).assembleConfig();
					end
				end
			end
		end

		function res = assembleInitialConditions(obj)
			%ASSEMBLEINITIALCONDITIONS Assembles the initial conditions according to the CADET file format spec
			%   RES = ASSEMBLEINITIALCONDITIONS() returns a nested Matlab struct RES that represents only the
			%   initial conditions part of the model as detailed in the (full configuration) CADET file format
			%   spec.
			%
			% See also MODEL.ASSEMBLEINITIALCONDITIONS, MODELSYSTEM.ASSEMBLEINITIALCONDITIONS

			res = obj.assembleInitialConditions@Model();

			res.INIT_C = obj.data.INIT_C;
			if isfield(obj.data, 'INIT_Q')
				res.INIT_Q = obj.data.INIT_Q;
			else
				res.INIT_Q = [];
			end
			res.INIT_VOLUME = obj.data.INIT_VOLUME;

			if isfield(obj.data, 'INIT_STATE')
				res.INIT_STATE = obj.data.INIT_STATE;
			end
		end

		function val = getParameterValue(obj, param)
			%GETPARAMETERVALUE Retrieves a parameter value from the model
			%   VAL = GETPARAMETERVALUE(PARAM) searches for the (single) parameter identified by
			%   the struct PARAM with the fields SENS_NAME, SENS_UNIT, SENS_COMP, SENS_PARTYPE,
			%   SENS_BOUNDPHASE, SENS_REACTION, and SENS_SECTION (as returned by MAKESENSITIVITY()).
			%   The returned value VAL contains the current value of the parameter (on the Matlab side,
			%   not in the current CADET configuration) or NaN if the parameter could not be found.
			%
			% See also MODEL.GETPARAMETERVALUE, MEXSIMULATOR.GETPARAMETERVALUE,
			%   STIRREDTANKMODEL.SETPARAMETERVALUE, MAKESENSITIVITY

			if param.SENS_UNIT ~= obj.unitOpIdx
				% Wrong unit operation
				val = nan;
				return;
			end
			
			if ~isfield(obj.data, param.SENS_NAME)
				% We don't have this parameter

				if ~isempty(obj.bindingModel)
					% Try binding model
					if (param.SENS_PARTYPE >= 0)
						val = obj.bindingModel(param.SENS_PARTYPE+1).getParameterValue(param, obj.nBoundStates((param.SENS_PARTYPE * obj.nComponents + 1):((param.SENS_PARTYPE + 1) * obj.nComponents)));
					else
						val = obj.bindingModel(1).getParameterValue(param, obj.nBoundStates(1:obj.nComponents));
					end
					return;
				end

				if ~isempty(obj.reactionModelBulk)
					% Try reaction model
					val = obj.reactionModelBulk.getParameterValue(param, obj.nBoundStates(1:obj.nComponents));
					return;
				end

				if ~isempty(obj.reactionModelParticle)
					% Try reaction model
					if (param.SENS_PARTYPE >= 0)
						val = obj.reactionModelParticle(param.SENS_PARTYPE+1).getParameterValue(param, obj.nBoundStates((param.SENS_PARTYPE * obj.nComponents + 1):((param.SENS_PARTYPE + 1) * obj.nComponents)));
					else
						val = obj.reactionModelParticle(1).getParameterValue(param, obj.nBoundStates(1:obj.nComponents));
					end
					return;
				end

				val = nan;
				return;
			end

			% The only parameter is section but not component dependent
			val = obj.data.(param.SENS_NAME);
			offset = 1;
			if (param.SENS_SECTION ~= -1)
				offset = offset + param.SENS_SECTION;
			end
			val = val(offset);
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
			% See also MODEL.SETPARAMETERVALUE, MEXSIMULATOR.SETPARAMETERVALUE,
			%   STIRREDTANKMODEL.GETPARAMETERVALUE, MAKESENSITIVITY

			if param.SENS_UNIT ~= obj.unitOpIdx
				% Wrong unit operation
				oldVal = nan;
				return;
			end

			if ~isfield(obj.data, param.SENS_NAME)
				% We don't have this parameter

				if ~isempty(obj.bindingModel)
					% Try binding model
					if (param.SENS_PARTYPE >= 0)
						oldVal = obj.bindingModel(param.SENS_PARTYPE+1).setParameterValue(param, obj.nBoundStates((param.SENS_PARTYPE * obj.nComponents + 1):((param.SENS_PARTYPE + 1) * obj.nComponents)), newVal);
					else
						oldVal = obj.bindingModel(1).setParameterValue(param, obj.nBoundStates(1:obj.nComponents), newVal);
					end
					return;
				end

				if ~isempty(obj.reactionModelBulk)
					% Try reaction model
					oldVal = obj.reactionModelBulk.setParameterValue(param, obj.nBoundStates(1:obj.nComponents), newVal);
					return;
				end

				if ~isempty(obj.reactionModelParticle)
					% Try reaction model
					if (param.SENS_PARTYPE >= 0)
						oldVal = obj.reactionModelParticle(param.SENS_PARTYPE+1).setParameterValue(param, obj.nBoundStates((param.SENS_PARTYPE * obj.nComponents + 1):((param.SENS_PARTYPE + 1) * obj.nComponents)), newVal);
					else
						oldVal = obj.reactionModelParticle(1).setParameterValue(param, obj.nBoundStates(1:obj.nComponents), newVal);
					end
					return;
				end

				oldVal = nan;
				return;
			end

			offset = 0;
			if (param.SENS_COMP ~= - 1)
				offset = offset + param.SENS_COMP;
			end
			if (param.SENS_SECTION ~= -1)
				offset = offset + param.SENS_SECTION * obj.nComponents;
			end
			oldVal = obj.data.(param.SENS_NAME)(offset + 1);
			obj.data.(param.SENS_NAME)(offset + 1) = newVal;
		end

		function notifySync(obj)
			%NOTIFYSYNC Called after parameters have been sent to CADET
			%   NOTIFYSYNC() resets the HASCHANGED property.

			obj.notifySync@Model();

			if ~isempty(obj.bindingModel)
				for i = 1:length(obj.bindingModel)
					obj.bindingModel(i).notifySync();
				end
			end

			if ~isempty(obj.reactionModelBulk)
				obj.reactionModelBulk.notifySync();
			end

			if ~isempty(obj.reactionModelParticle)
				for i = 1:length(obj.reactionModelParticle)
					obj.reactionModelParticle(i).notifySync();
				end
			end
		end
	end

	methods (Access = 'protected')

		function val = getHasChanged(obj)
			val = obj.getHasChanged@Model();

			% Check binding models
			if ~isempty(obj.bindingModel)
				for i = 1:length(obj.bindingModel)
					if obj.bindingModel(i).hasChanged
						val = true;
						return
					end
				end
			end

			% Check reaction models
			if ~isempty(obj.reactionModelBulk)
				if obj.reactionModelBulk.hasChanged
					val = true;
					return
				end
			end

			if ~isempty(obj.reactionModelParticle)
				for i = 1:length(obj.reactionModelParticle)
					if obj.reactionModelParticle(i).hasChanged
						val = true;
						return
					end
				end
			end
		end

		function loadobjInternal(obj, S)
			obj.loadobjInternal@Model(S);

			obj.bindingModel = BindingModel.empty();
			if ~isempty(S.bindingModel)
				for i = 1:length(S.bindingModel)
					ctor = str2func([S.bindingModelClass{i} '.loadobj']);
					obj.bindingModel(i) = ctor(S.bindingModel{i});
				end
			end

			if ~isempty(S.reactionModelBulk)
				ctor = str2func([S.reactionModelBulkClass '.loadobj']);
				obj.reactionModelBulk = ctor(S.reactionModelBulk);
			else
				obj.reactionModelBulk = ReactionModel.empty();
			end

			obj.reactionModelParticle = ReactionModel.empty();
			if ~isempty(S.reactionModelParticle)
				for i = 1:length(S.reactionModelParticle)
					ctor = str2func([S.reactionModelParticleClass{i} '.loadobj']);
					obj.reactionModelParticle(i) = ctor(S.reactionModelParticle{i});
				end
			end
		end

	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = StirredTankModel();
				obj.loadobjInternal(S);
			end
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
