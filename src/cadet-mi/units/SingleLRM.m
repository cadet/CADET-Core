
classdef SingleLRM < LumpedRateModelWithoutPores & SingleUnitOpSystem
	%SingleLRM Represents a system with one lumped rate model
	%   This class behaves like one single column (lumped rate model)
	%   fed by a piecewise cubic polynomial inlet profile.
	%
	%   The lumped rate model is assigned the id 0 and the implicit inlet model
	%   is assigned id 1.
	%
	%   This class is almost 'plug-in' compatible to the CADET 2.0 Matlab interface.
	%
	% See also LUMPEDRATEMODELWITHOUTPORES, SINGLEUNITOPSYSTEM, PIECEWISECUBICPOLYPROFILE
	
	% Copyright: (C) 2008-2020 The CADET Authors
	%            See the license note at the end of the file.

	properties
		inlet; % Piecewise cubic polynomial inlet profile
	end

	properties (Dependent, Transient)
		% Constant coefficient for each component and section in [mol / m^3_IV];
		% each row corresponds to a section, each column to a component such that the
		% format is nSections x nComponents.
		constant;
		% Linear coefficient for each component and section in [mol / (m^3_IV * s)];
		% each row corresponds to a section, each column to a component such that the
		% format is nSections x nComponents.
		linear;
		% Quadratic coefficient for each component and section in [mol / (m^3_IV * s^2)];
		% each row corresponds to a section, each column to a component such that the
		% format is nSections x nComponents.
		quadratic;
		% Cubic coefficient for each component and section in [mol / (m^3_IV * s^3)];
		% each row corresponds to a section, each column to a component such that the
		% format is nSections x nComponents.
		cubic;
	end

	methods
		
		function obj = SingleLRM()
			%SINGLELRM Constructs a SingleLRM object

			obj = obj@LumpedRateModelWithoutPores();
			obj@SingleUnitOpSystem();
			obj.inlet = PiecewiseCubicPolyProfile();

			obj.unitOpIdx = 0;
		end
		

		function S = saveobj(obj)
			S = obj.saveobj@LumpedRateModelWithoutPores();
			S.singleUnitOpSys = obj.saveobj@SingleUnitOpSystem();
			S.inletClass = class(obj.inlet);
			S.inlet = obj.inlet.saveobj();
		end


		function res = validate(obj, sectionTimes, subModels)
			%VALIDATE Validates the configuration of the main unit operation
			%   RES = VALIDATE(SECTIONTIMES) uses the section times of the simulator given in
			%   SECTIONTIMES to validate the LRM in this system, the inlet model, and the
			%   system itself. Returns true in RES if everything is fine and false otherwise.
			%
			%   RES = VALIDATE(..., SUBMODELS) determines whether the main unit operation and
			%   inlet model are also validated (SUBMODELS = true, default) or not (SUBMODELS = false).
			%
			% See also MODELSYSTEM.VALIDATE, MODEL.VALIDATE

			if (nargin <= 2) || isempty(subModels)
				subModels = true;
			end

			res = obj.validate@SingleUnitOpSystem(sectionTimes);
			if subModels
				res = res && obj.inlet.validate(sectionTimes, obj.nComponents) && ...
				   obj.validate@LumpedRateModelWithoutPores(sectionTimes);
			end
		end

		function res = assembleConfig(obj, unitOp)
			%ASSEMBLECONFIG Assembles the configuration according to the CADET file format spec
			%   RES = ASSEMBLECONFIG() returns a nested Matlab struct RES that represents the
			%   system (including LRM and inlet) as detailed in the CADET file format spec.
			%
			%   RES = ASSEMBLECONFIG(UNITOP) returns the configuration of the full system
			%   (UNITOP = 'all' or []), the system itself (UNITOP = 'system' or -1), the
			%   general rate model (UNITOP = 'lrm' or 0), or just the inlet (UNITOP = 'inlet'
			%   or 1) as detailed in the CADET file format spec.
			%
			% See also MODELSYSTEM.ASSEMBLECONFIG, MEXSIMULATOR.ASSEMBLECONFIG

			if (nargin <= 1)
				unitOp = [];
			end

			if ischar(unitOp)
				switch unitOp
					case 'all'
						unitOp = [];
					case 'lrm'
					case 'lrmp'
						unitOp = 0;
					case 'inlet'
						unitOp = 1;
					case 'system'
						unitOp = -1;
					otherwise
						error('CADET:funcParamError', 'Invalid option "%s" for unitOp.', unitOp);
				end
			end

			if isempty(unitOp)
				inletConfig = obj.inlet.assembleConfig();
				inletConfig.INLET_TYPE = obj.inlet.name;
				res = obj.assembleConfigBase(obj.assembleConfig@LumpedRateModelWithoutPores(), inletConfig);
			else
				switch unitOp
					case -1
						res = obj.assembleConfigBase([], []);
					case 0
						res = obj.assembleConfig@LumpedRateModelWithoutPores();
					case 1
						res = obj.inlet.assembleConfig();
						res.INLET_TYPE = obj.inlet.name;
						res.UNIT_TYPE = 'INLET';
						res.NCOMP = int32(obj.nComponents);
					otherwise
						error('CADET:funcParamError', 'Invalid option %d for unitOp.', unitOp);
				end
			end
		end

		function res = assembleReturnConfig(obj)
			%ASSEMBLERETURNCONFIG Assembles the return configuration according to the CADET file format spec
			%   RES = ASSEMBLERETURNCONFIG() returns a nested Matlab struct RES that contains the return
			%   configuration (return group in the file format spec) of the full system including the general
			%   rate model and the inlet.
			%
			% See also MODEL.ASSEMBLERETURNCONFIG, MODELSYSTEM.ASSEMBLERETURNCONFIG, MEXSIMULATOR.ASSEMBLERETURNCONFIG

			res = obj.assembleReturnConfigBase(obj.assembleReturnConfig@LumpedRateModelWithoutPores());
		end

		function res = assembleInitialConditions(obj)
			%ASSEMBLEINITIALCONDITIONS Assembles the initial conditions according to the CADET file format spec
			%   RES = ASSEMBLEINITIALCONDITIONS() returns a nested Matlab struct RES that represents only the
			%   initial conditions part of the system, the general rate model, and the inlet as detailed in
			%   the (full configuration) CADET file format spec.
			%
			% See also MODEL.ASSEMBLEINITIALCONDITIONS, MODELSYSTEM.ASSEMBLEINITIALCONDITIONS

			res = obj.assembleInitialConditionsBase(obj.assembleInitialConditions@LumpedRateModelWithoutPores());
		end

		function val = get.inlet(obj)
			val = obj.inlet;
		end

		function set.inlet(obj, val)
			validateattributes(val, {'PiecewiseCubicPolyProfile'}, {'nonempty', 'scalar'}, '', 'val');
			obj.inlet = val;
			obj.hasInletChanged = true;
		end

		function val = get.constant(obj)
			val = obj.inlet.constant;
		end

		function set.constant(obj, val)
			obj.inlet.constant = val;
			obj.hasInletChanged = true;
		end

		function val = get.linear(obj)
			val = obj.inlet.linear;
		end

		function set.linear(obj, val)
			obj.inlet.linear = val;
			obj.hasInletChanged = true;
		end

		function val = get.quadratic(obj)
			val = obj.inlet.quadratic;
		end

		function set.quadratic(obj, val)
			obj.inlet.quadratic = val;
			obj.hasInletChanged = true;
		end

		function val = get.cubic(obj)
			val = obj.inlet.cubic;
		end

		function set.cubic(obj, val)
			obj.inlet.cubic = val;
			obj.hasInletChanged = true;
		end

		function val = getParameterValue(obj, param)
			%GETPARAMETERVALUE Retrieves a parameter value from the model or inlet
			%   VAL = GETPARAMETERVALUE(PARAM) searches for the (single) parameter identified by
			%   the struct PARAM with the fields SENS_NAME, SENS_UNIT, SENS_COMP, SENS_PARTYPE,
			%   SENS_BOUNDPHASE, SENS_REACTION, and SENS_SECTION (as returned by MAKESENSITIVITY()).
			%   The returned value VAL contains the current value of the parameter (on the Matlab side,
			%   not in the current CADET configuration) or NaN if the parameter could not be found.
			%
			% See also MODEL.GETPARAMETERVALUE, MEXSIMULATOR.GETPARAMETERVALUE, SINGLELRM.SETPARAMETERVALUE,
			%   MAKESENSITIVITY

			val = obj.getParameterValue@SingleUnitOpSystem(param);
			if ~isnan(val)
				return;
			end

			if param.SENS_UNIT == 0
				val = obj.getParameterValue@LumpedRateModelWithoutPores(param);
			elseif param.SENS_UNIT == 1
				val = obj.inlet.getParameterValue(param);
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
			% See also MODEL.SETPARAMETERVALUE, MEXSIMULATOR.SETPARAMETERVALUE, SINGLELRM.GETPARAMETERVALUE,
			%   MAKESENSITIVITY

			oldVal = obj.setParameterValue@SingleUnitOpSystem(param, newVal);
			if ~isnan(oldVal)
				return;
			end

			if param.SENS_UNIT == 0
				oldVal = obj.setParameterValue@LumpedRateModelWithoutPores(param, newVal);
			elseif param.SENS_UNIT == 1
				oldVal = obj.inlet.setParameterValue(param, newVal);
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

			obj.notifySync@SingleUnitOpSystem(systemOnly);

			if ~systemOnly
				obj.notifySync@LumpedRateModelWithoutPores();
				obj.inlet.notifySync();
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

			if obj.hasInletChanged || obj.inlet.hasChanged
				changedUnits = [changedUnits; 1];
			end

			if obj.hasSystemChanged
				changedUnits = [changedUnits; -1];
			end
		end

	end

	methods (Access = 'protected')

		function loadobjInternal(obj, S)
			obj.loadobjInternal@LumpedRateModelWithoutPores(S);
			obj.loadobjInternal@SingleUnitOpSystem(S.singleUnitOpSys);

			ctor = str2func([S.inletClass '.loadobj']);
			obj.inlet = ctor(S.inlet);
		end

	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = SingleLRM();
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
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
