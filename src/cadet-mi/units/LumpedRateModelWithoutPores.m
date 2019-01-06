
classdef LumpedRateModelWithoutPores < Model
	%LumpedRateModelWithoutPores Represents the lumped rate model of liquid column chromatography
	%   Holds all model specific parameters (e.g., interstitial velocity, dispersion, etc.) which are necessary
	%   for a CADET simulation. This class is mainly concerned with the transport model. Binding is handled by a
	%   binding model class.
	%
	% See also MODEL, SINGLELRM, MODELSYSTEM
	
	% Copyright: (C) 2008-2018 The CADET Authors
	%            See the license note at the end of the file.

	properties
		bindingModel; % Object of the used binding model class
	end

	properties (Dependent)
		
		% Discretization

		nCellsColumn; % Number of axial cells in the column
		nBoundStates; % Number of bound states for each component

		useAnalyticJacobian; % Determines whether Jacobian is calculated analytically or via AD

		% Initial values
		initialBulk; % Initial concentrations for each component in the bulk volume in [mol / m^3_IV]
		initialSolid; % Initial concentrations for each component in the bead solid volume in [mol / m^3_SP]
		initialState; % Initial concentrations for each degree of freedom

		% Transport
		
		dispersionColumn; % Dispersion coefficient of the mobile phase transport inside the column in [m^2_IV / s]
		interstitialVelocity; % Interstitial velocity of the mobile phase transport inside the column in [m_IV / s]

		% Geometry
		
		porosity; % Total porosity
		columnLength; % Length of the column in [m]
		crossSectionArea; % Cross section area of the column in [m]

		% Numerical method for advection
		
		reconstructionType; % Type of reconstruction
		wenoBoundaryHandling; % How WENO handles left and right boundary
		wenoEpsilon; % Epsilon in the WENO method (for continuity detection)
		wenoOrder; % Order of the WENO method

	end
	
	properties (Constant)
		name = 'LUMPED_RATE_MODEL_WITHOUT_PORES'; % Type of the model according to CADET file format specs
		hasInlet = true; % Determines whether the unit operation has an inlet
		hasOutlet = true; % Determines whether the unit operation has an outlet
	end

	methods
		
		function obj = LumpedRateModelWithoutPores()
			%LUMPEDRATEMODELWITHOUTPORES Constructs a LumpedRateModelWithoutPores object and inserts as much default values as possible

			obj = obj@Model();
			obj.data.discretization = [];
			obj.data.discretization.weno = [];

			% Set some default values
			obj.bindingModel = [];

			obj.useAnalyticJacobian = true;

			obj.reconstructionType = 'WENO';
			obj.wenoBoundaryHandling = 0; % Decrease order of WENO scheme at boundary
			obj.wenoEpsilon = 1e-12;
			obj.wenoOrder = 3; % Largest possible order
		end
		

		% Discretization

		function val = get.nCellsColumn(obj)
			val = double(obj.data.discretization.NCOL);
		end

		function set.nCellsColumn(obj, val)
			validateattributes(val, {'numeric'}, {'scalar', 'nonempty', '>=', 1, 'finite', 'real'}, '', 'nCellsColumn');
			obj.data.discretization.NCOL = int32(val);
			obj.hasChanged = true;
		end

		function val = get.nBoundStates(obj)
			if isfield(obj.data.discretization, 'NBOUND')
				val = double(obj.data.discretization.NBOUND);
			else
				val = [];
			end
		end

		function set.nBoundStates(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data.discretization, 'NBOUND');
			else
				validateattributes(val, {'numeric'}, {'nonnegative', 'vector', 'nonempty', 'finite', 'real'}, '', 'nBoundStates');
				obj.data.discretization.NBOUND = int32(val);
			end
			obj.hasChanged = true;
		end

		function val = get.useAnalyticJacobian(obj)
			val = logical(obj.data.discretization.USE_ANALYTIC_JACOBIAN);
		end

		function set.useAnalyticJacobian(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'useAnalyticJacobian');
			obj.data.discretization.USE_ANALYTIC_JACOBIAN = int32(logical(val));
			obj.hasChanged = true;
		end

		% Initial values

		function val = get.initialBulk(obj)
			val = obj.data.INIT_C;
		end

		function set.initialBulk(obj, val)
			validateattributes(val, {'double'}, {'nonnegative', 'vector', 'nonempty', 'finite', 'real'}, '', 'initialBulk');
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

		% Transport
		
		function val = get.dispersionColumn(obj)
			val = obj.data.COL_DISPERSION;
		end

		function set.dispersionColumn(obj, val)
			validateattributes(val, {'double'}, {'nonnegative', 'vector', 'nonempty', 'finite', 'real'}, '', 'dispersionColumn');
			obj.data.COL_DISPERSION = val(:);
			obj.hasChanged = true;
		end

		function val = get.interstitialVelocity(obj)
			if isfield(obj.data, 'VELOCITY')
				val = obj.data.VELOCITY;
			else
				val = [];
			end
		end

		function set.interstitialVelocity(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'VELOCITY');
			else
				validateattributes(val, {'double'}, {'positive', 'vector', 'nonempty', 'finite', 'real'}, '', 'interstitialVelocity');
				obj.data.VELOCITY = val(:);
			end
			obj.hasChanged = true;
		end

		% Geometry
		
		function val = get.porosity(obj)
			val = obj.data.TOTAL_POROSITY;
		end

		function set.porosity(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', '>=', 0.0, '<=', 1.0, 'finite', 'real'}, '', 'porosity');
			obj.data.TOTAL_POROSITY = val;
			obj.hasChanged = true;
		end

		function val = get.columnLength(obj)
			val = obj.data.COL_LENGTH;
		end

		function set.columnLength(obj, val)
			validateattributes(val, {'double'}, {'positive', 'scalar', 'nonempty', 'finite', 'real'}, '', 'columnLength');
			obj.data.COL_LENGTH = val;
			obj.hasChanged = true;
		end

		function val = get.crossSectionArea(obj)
			if isfield(obj.data, 'CROSS_SECTION_AREA')
				val = obj.data.CROSS_SECTION_AREA;
			else
				val = [];
			end
		end

		function set.crossSectionArea(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'CROSS_SECTION_AREA');
			else
				validateattributes(val, {'double'}, {'positive', 'scalar', 'nonempty', 'finite', 'real'}, '', 'crossSectionArea');
				obj.data.CROSS_SECTION_AREA = val;
			end
			obj.hasChanged = true;
		end

		% Numerical method for advection
		
		function val = get.reconstructionType(obj)
			val = obj.data.discretization.RECONSTRUCTION;
		end

		function set.reconstructionType(obj, val)
			obj.data.discretization.RECONSTRUCTION = validatestring(val, {'WENO'}, '', 'reconstructionType');
			obj.hasChanged = true;
		end

		function val = get.wenoBoundaryHandling(obj)
			val = double(obj.data.discretization.weno.BOUNDARY_MODEL);
		end

		function set.wenoBoundaryHandling(obj, val)
			validateattributes(val, {'numeric'}, {'>=', 0, '<=', 3, 'scalar', 'nonempty', 'finite', 'real'}, '', 'wenoBoundaryHandling');
			obj.data.discretization.weno.BOUNDARY_MODEL = int32(val);
			obj.hasChanged = true;
		end

		function val = get.wenoEpsilon(obj)
			val = obj.data.discretization.weno.WENO_EPS;
		end

		function set.wenoEpsilon(obj, val)
			validateattributes(val, {'double'}, {'positive', 'scalar', 'nonempty', 'finite', 'real'}, '', 'wenoEpsilon');
			obj.data.discretization.weno.WENO_EPS = val;
			obj.hasChanged = true;
		end

		function val = get.wenoOrder(obj)
			val = double(obj.data.discretization.weno.WENO_ORDER);
		end

		function set.wenoOrder(obj, val)
			validateattributes(val, {'numeric'}, {'>=', 1, 'scalar', 'nonempty', '<=', 3, 'finite', 'real'}, '', 'wenoOrder');
			obj.data.discretization.weno.WENO_ORDER = int32(val);
			obj.hasChanged = true;
		end

		function set.bindingModel(obj, val)
			if ~isempty(val) && ~isa(val, 'BindingModel')
				error('CADET:invalidConfig', 'Expected a valid binding model.');
			end
			obj.bindingModel = val;
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
			nSections = numel(sectionTimes) - 1;
			if ~isfield(obj.data, 'COL_DISPERSION')
				error('CADET:invalidConfig', 'Property dispersionColumn must be set.');
			end
			if (~isfield(obj.data, 'VELOCITY')) && (~isfield(obj.data, 'CROSS_SECTION_AREA'))
				error('CADET:invalidConfig', 'Property interstitialVelocity or crossSectionArea must be set.');
			end
			if ~isfield(obj.data, 'TOTAL_POROSITY')
				error('CADET:invalidConfig', 'Property porosity must be set.');
			end
			if ~isfield(obj.data, 'COL_LENGTH')
				error('CADET:invalidConfig', 'Property columnLength must be set.');
			end

			validateattributes(obj.nCellsColumn, {'numeric'}, {'scalar', 'nonempty', '>=', 1, 'finite', 'real'}, '', 'nCellsColumn');
			if isfield(obj.data.discretization, 'NBOUND')
				validateattributes(obj.nBoundStates, {'numeric'}, {'nonnegative', 'vector', 'numel', obj.nComponents, 'finite', 'real'}, '', 'nBoundStates');
				nTotalBnd = sum(obj.nBoundStates);
			else
				nTotalBnd = 0;
			end
			
			validateattributes(obj.initialBulk, {'double'}, {'nonnegative', 'vector', 'numel', obj.nComponents, 'finite', 'real'}, '', 'initialBulk');
			if ~isempty(obj.initialSolid)
				validateattributes(obj.initialSolid, {'double'}, {'nonnegative', 'vector', 'numel', nTotalBnd, 'finite', 'real'}, '', 'initialSolid');
			end
			if ~isempty(obj.initialState)
				nDof = obj.nCellsColumn * (obj.nComponents + nTotalBnd);
				if (numel(obj.initialState) ~= nDof) && (numel(obj.initialState) ~= 2 * nDof)
					error('CADET:invalidConfig', 'Expected initialState to be of size %d or %d.', nDof, 2*nDof);
				end
			end

			if (numel(obj.dispersionColumn) ~= 1) && (numel(obj.dispersionColumn) ~= nSections)
				error('CADET:invalidConfig', 'Expected dispersionColumn to be of size %d or %d (number of time sections).', 1, nSections);
			end
			if ~isempty(obj.interstitialVelocity)
				if (numel(obj.interstitialVelocity) ~= 1) && (numel(obj.interstitialVelocity) ~= nSections)
					error('CADET:invalidConfig', 'Expected interstitialVelocity to be of size %d or %d (number of time sections).', 1, nSections);
				end
			end

			validateattributes(obj.porosity, {'double'}, {'scalar', 'nonempty', '>=', 0.0, '<=', 1.0, 'finite', 'real'}, '', 'porosity');
			validateattributes(obj.columnLength, {'double'}, {'positive', 'scalar', 'nonempty', 'finite', 'real'}, '', 'columnLength');

			if ~isempty(obj.bindingModel) && ~isa(obj.bindingModel, 'BindingModel')
				error('CADET:invalidConfig', 'Expected a valid binding model.');
			end
			if isempty(obj.bindingModel) && (nTotalBnd ~= 0)
				error('CADET:invalidConfig', 'Expected no bound states when using no binding model.');
			end

			if ~isempty(obj.bindingModel)
				validateattributes(obj.nBoundStates, {'numeric'}, {'nonnegative', 'vector', 'numel', obj.nComponents, 'finite', 'real'}, '', 'nBoundStates');
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

			if ~isfield(obj.data.discretization, 'NBOUND') || isempty(obj.data.discretization.NBOUND)
				res.discretization.NBOUND = int32(zeros(obj.nComponents, 1));
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
			if isfield(obj.data, 'INIT_Q')
				res.INIT_Q = obj.data.INIT_Q;
			else
				res.INIT_Q = [];
			end

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
			%   LUMPEDRATEMODELWITHOUTPORES.SETPARAMETERVALUE, MAKESENSITIVITY

			if param.SENS_UNIT ~= obj.unitOpIdx
				% Wrong unit operation
				val = nan;
				return;
			end

			if ~isfield(obj.data, param.SENS_NAME) && ~isempty(obj.bindingModel)
				% We don't have this parameter, so try binding model
				val = obj.bindingModel.getParameterValue(param, obj.nBoundStates);
				return;
			end
			
			val = obj.data.(param.SENS_NAME);
			offset = 0;
			if (param.SENS_COMP ~= - 1)
				offset = offset + param.SENS_COMP;
			end
			if (param.SENS_SECTION ~= -1)
				offset = offset + param.SENS_SECTION * obj.nComponents;
			end
			val = val(offset + 1);
		end

		function oldVal = setParameterValue(obj, param, newVal)
			%SETPARAMETERVALUE Sets a parameter value in the model
			%   OLDVAL = SETPARAMETERVALUE(PARAM, NEWVAL) searches for the parameter
			%   identified by the struct PARAM with the fields SENS_NAME, SENS_UNIT,
			%   SENS_COMP, SENS_PARTYPE, SENS_BOUNDPHASE, SENS_REACTION, and SENS_SECTION
			%   (as returned by MAKESENSITIVITY()). The returned value OLDVAL contains the old value
			%   of the parameter (on the Matlab side, not in the current CADET configuration)
			%   or NaN if the parameter could not be found. The value of the parameter is
			%   set to NEWVAL. The changes are not propagated to the underlying CADET simulator.
			%
			% See also MODEL.SETPARAMETERVALUE, MEXSIMULATOR.SETPARAMETERVALUE,
			%   LUMPEDRATEMODELWITHOUTPORES.GETPARAMETERVALUE, MAKESENSITIVITY

			if param.SENS_UNIT ~= obj.unitOpIdx
				% Wrong unit operation
				oldVal = nan;
				return;
			end

			if ~isfield(obj.data, param.SENS_NAME) && ~isempty(obj.bindingModel)
				% We don't have this parameter, so try binding model
				oldVal = obj.bindingModel.setParameterValue(param, obj.nBoundStates, newVal);
				return;
			end

			offset = 0;
			if (param.SENS_COMP ~= - 1)
				offset = offset + param.SENS_COMP;
				if (param.SENS_SECTION ~= -1)
					offset = offset + param.SENS_SECTION * obj.nComponents;
				end
			elseif (param.SENS_SECTION ~= -1)
				offset = offset + param.SENS_SECTION;
			end
			oldVal = obj.data.(param.SENS_NAME)(offset + 1);
			obj.data.(param.SENS_NAME)(offset + 1) = newVal;
		end

		function notifySync(obj)
			%NOTIFYSYNC Called after parameters have been sent to CADET
			%   NOTIFYSYNC() resets the HASCHANGED property.

			obj.notifySync@Model();
			if ~isempty(obj.bindingModel)
				obj.bindingModel.notifySync();
			end
		end
	end

	methods (Access = 'protected')

		function val = getHasChanged(obj)
			val = obj.getHasChanged@Model() || (~isempty(obj.bindingModel) && obj.bindingModel.hasChanged);
		end

		function loadobjInternal(obj, S)
			obj.loadobjInternal@Model(S);

			if ~isempty(S.bindingModel)
				ctor = str2func([S.bindingModelClass '.loadobj']);
				obj.bindingModel = ctor(S.bindingModel);
			end
		end

	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = LumpedRateModelWithoutPores();
				obj.loadobjInternal(S);
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
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
