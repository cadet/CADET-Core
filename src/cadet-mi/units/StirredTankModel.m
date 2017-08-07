
classdef StirredTankModel < Model
	%StirredTankModel Represents a continuously stirred (reaction) tank model (short CSTR)
	%   Holds all model specific parameters which are necessary for a CADET simulation.
	%
	% See also MODEL, SINGLECSTR, MODELSYSTEM
	
	% Copyright: (C) 2008-2017 The CADET Authors
	%            See the license note at the end of the file.

	properties
		bindingModel; % Object of the used binding model class
	end

	properties (Dependent)
		% Parameters
		flowRateFilter; % Flow rate of filtered out liquid in [m^3 / s]
		porosity; % Porosity
		nBoundStates; % Number of bound states for each component

		% Initial values
		initialConcentration; % Initial concentrations for each component in [mol / m^3]
		initialSolid; % Initial concentrations for each component in the solid phase in [mol / m^3_SP]
		initialVolume; % Initial volume in [m^3]
		initialState; % Initial values for each degree of freedom

		useAnalyticJacobian; % Determines whether Jacobian is calculated analytically or via AD
	end
	
	properties (Constant)
		name = 'CSTR'; % Type of the model according to CADET file format specs
		hasInlet = true; % Determines whether the unit operation has an inlet
		hasOutlet = true; % Determines whether the unit operation has an outlet
	end

	methods
		
		function obj = StirredTankModel()
			%STIRREDTANKMODEL Constructs a StirredTankModel object and inserts as much default values as possible

			obj = obj@Model();

			% Set some default values
			obj.bindingModel = [];
			obj.flowRateFilter = 0.0;
			obj.useAnalyticJacobian = true;

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
			val = double(obj.data.NBOUND);
		end

		function set.nBoundStates(obj, val)
			validateattributes(val, {'numeric'}, {'nonnegative', 'vector', 'nonempty', 'finite', 'real'}, '', 'nBoundStates');
			obj.data.NBOUND = int32(val);
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
			val = obj.data.INIT_Q;
		end

		function set.initialSolid(obj, val)
			validateattributes(val, {'double'}, {'nonnegative', 'vector', 'empty', 'finite', 'real'}, '', 'initialSolid');
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


		function S = saveobj(obj)
			S = obj.saveobj@Model();
			S.bindingModel = [];

			if ~isempty(obj.bindingModel)
				S.bindingModelClass = class(obj.bindingModel);
				S.bindingModel = obj.bindingModel.saveobj();
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

			validateattributes(obj.nBoundStates, {'numeric'}, {'nonnegative', 'vector', 'numel', obj.nComponents, 'finite', 'real'}, '', 'nBoundStates');
			validateattributes(obj.initialConcentration, {'double'}, {'nonnegative', 'vector', 'numel', obj.nComponents, 'finite', 'real'}, '', 'initialConcentration');
			validateattributes(obj.initialSolid, {'double'}, {'nonnegative', 'vector', 'numel', sum(obj.nBoundStates), 'finite', 'real'}, '', 'initialSolid');
			validateattributes(obj.initialVolume, {'double'}, {'nonnegative', 'scalar', 'finite', 'real'}, '', 'initialVolume');
			if ~isempty(obj.initialState)
				nDof = nComponents + 1 + sum(obj.nBoundStates);
				if (numel(obj.initialState) ~= nDof) && (numel(obj.initialState) ~= 2 * nDof)
					error('CADET:invalidConfig', 'Expected initialState to be of size %d or %d.', nDof, 2*nDof);
				end
				validateattributes(obj.initialState(1:nDof), {'double'}, {'nonnegative', 'finite', 'real'}, '', 'initialState');
			end

			if (numel(obj.flowRateFilter) ~= 1) && (numel(obj.flowRateFilter) ~= nSections)
				error('CADET:invalidConfig', 'Expected flowRateFilter to be of size %d or %d (number of time sections).', 1, nSections);
			end

			validateattributes(obj.porosity, {'double'}, {'scalar', 'nonempty', '>=', 0.0, '<=', 1.0, 'finite', 'real'}, '', 'porosity');
			if ~isempty(obj.bindingModel) && ~isa(obj.bindingModel, 'BindingModel')
				error('CADET:invalidConfig', 'Expected a valid binding model.');
			end
			if isempty(obj.bindingModel) && (sum(obj.nBoundStates) ~= 0)
				error('CADET:invalidConfig', 'Expected no bound states when using no binding model.');
			end

			if ~isempty(obj.bindingModel)
				res = obj.bindingModel.validate(obj.nComponents, obj.nBoundStates) && res;
			end
		end

		function res = assembleConfig(obj)
			%ASSEMBLECONFIG Assembles the configuration according to the CADET file format spec
			%   RES = ASSEMBLECONFIG() returns a nested Matlab struct RES that represents the
			%   model as detailed in the CADET file format spec.
			%
			% See also MODEL.ASSEMBLECONFIG, MEXSIMULATOR.ASSEMBLECONFIG

			res = obj.assembleConfig@Model();

			if ~isempty(obj.bindingModel) && ~isa(obj.bindingModel, 'BindingModel')
				error('CADET:invalidConfig', 'Expected a valid binding model.');
			end

			if isempty(obj.bindingModel)
				res.ADSORPTION_MODEL = 'NONE';
			else
				res.ADSORPTION_MODEL = obj.bindingModel.name;
				res.adsorption = obj.bindingModel.assembleConfig();
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
			res.INIT_Q = obj.data.INIT_Q;
			res.INIT_VOLUME = obj.data.INIT_VOLUME;

			if isfield(obj.data, 'INIT_STATE')
				res.INIT_STATE = obj.data.INIT_STATE;
			end
		end

		function val = getParameterValue(obj, param)
			%GETPARAMETERVALUE Retrieves a parameter value from the model
			%   VAL = GETPARAMETERVALUE(PARAM) searches for the (single) parameter identified by
			%   the struct PARAM with the fields SENS_NAME, SENS_UNIT, SENS_COMP, SENS_BOUNDPHASE,
			%   SENS_REACTION, and SENS_SECTION (as returned by MAKESENSITIVITY()). The returned
			%   value VAL contains the current value of the parameter (on the Matlab side, not in
			%   the current CADET configuration) or NaN if the parameter could not be found.
			%
			% See also MODEL.GETPARAMETERVALUE, MEXSIMULATOR.GETPARAMETERVALUE,
			%   STIRREDTANKMODEL.SETPARAMETERVALUE, MAKESENSITIVITY

			if param.SENS_UNIT ~= obj.unitOpIdx
				% Wrong unit operation
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
			%   SENS_COMP, SENS_BOUNDPHASE, SENS_REACTION, and SENS_SECTION (as returned
			%   by MAKESENSITIVITY()). The returned value OLDVAL contains the old value
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
		end
	end

	methods (Access = 'protected')

		function val = getHasChanged(obj)
			val = obj.getHasChanged@Model();
		end

		function loadobjInternal(obj, S)
			obj.loadobjInternal@Model(S);

			if ~isempty(S.bindingModel)
				ctor = str2func([S.bindingModelClass '.loadobj']);
				obj.bindingModel = ctor(S.bindingModel);
			else
				obj.bindingModel = [];
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
%  CADET - The Chromatography Analysis and Design Toolkit
%  
%  Copyright (C) 2008-2017: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
