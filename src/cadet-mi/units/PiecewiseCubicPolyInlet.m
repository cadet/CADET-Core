
classdef PiecewiseCubicPolyInlet < InletModel
	%PiecewiseCubicPolyInlet Uses piecewise cubic polynomials as inlet profile
	%   Piecewise cubic polynomials can be used as a cubic spline (if certain continuity 
	%   conditions are met at the interfaces of two adjacent pieces), but can also represent
	%   discontinuous data (such as step profiles).
	%
	% See also PIECEWISECUBICPOLYPROFILE
	
	% Copyright: (C) 2008-2024 The CADET Authors
	%            See the license note at the end of the file.
	
	properties (Access = 'protected', Transient)
		profile; % Piecewise cubic polynomial profile
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
		
		function obj = PiecewiseCubicPolyInlet()
			%PIECEWISECUBICPOLYINLET Constructs an INLETMODEL object using a piecewise cubic polynomial profile

			obj = obj@InletModel();
			obj.profile = PiecewiseCubicPolyProfile();
			obj.data.INLET_TYPE = obj.profile.name;
		end
		
		function S = saveobj(obj)
			S = obj.saveobj@InletModel();
			S.profileClass = class(obj.profile);
			S.profile = obj.profile.saveobj();
		end


		function res = validate(obj, sectionTimes)
			%VALIDATE Validates the configuration of the Model
			%   RES = VALIDATE(SECTIONTIMES) uses the section times of the simulator given in
			%   SECTIONTIMES to validate the Model object. Returns true in RES if everything
			%   is fine and false otherwise.
			%
			% See also MODEL.VALIDATE

			res = obj.validate@InletModel(sectionTimes) & obj.profile.validate(sectionTimes, obj.nComponents);
		end

		function res = assembleConfig(obj)
			%ASSEMBLECONFIG Assembles the configuration according to the CADET file format spec
			%   RES = ASSEMBLECONFIG() returns a nested Matlab struct RES that represents the
			%   model as detailed in the CADET file format spec.
			%
			% See also MODEL.ASSEMBLECONFIG, MEXSIMULATOR.ASSEMBLECONFIG

			res = obj.assembleConfig@InletModel();
			res = obj.profile.assembleConfig(res);
		end

		function pp = makePiecewisePoly(obj, sectionTimes, comp, diffOrder)
			%MAKEPIECEWISEPOLY Assembles a Matlab piecewise polynomial that can be used by ppval
			%   PP = MAKEPIECEWISEPOLY(SECTIONTIMES, COMP) returns a piecewise polynomial in Matlab
			%   format that mimics the PiecewiseCubicPolyProfile of one component and can be passed
			%   to PPVAL, UNMKPP, etc. For this single component, it is fully equivalent to the profile.
			%   SECTION_TIMES is a vector with section times (start, end) and COMP is the index of 
			%   the requested component.
			%
			%   PP = MAKEPIECEWISEPOLY(..., DIFFORDER) returns the DIFFORDER'th derivative of the
			%   profile. DIFFORDER is the order of the derivative and defaults to 0.
			%
			% See also PIECEWISECUBICPOLYPROFILE.MAKEPIECEWISEPOLY

			validateattributes(comp, {'numeric'}, {'>=', 1, '<=', obj.nComponents, 'finite', 'real'}, '', 'comp');
			if nargin <= 3
				diffOrder = 0;
			end
			pp = obj.profile.makePiecewisePoly(sectionTimes, comp, diffOrder);
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
			%   PIECEWISECUBICPOLYINLET.SETPARAMETERVALUE, MAKESENSITIVITY

			if param.SENS_UNIT ~= obj.unitOpIdx
				val = nan;
				return;
			end
			val = obj.profile.getParameterValue(param);
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
			%   PIECEWISECUBICPOLYINLET.GETPARAMETERVALUE, MAKESENSITIVITY

			if param.SENS_UNIT ~= obj.unitOpIdx
				oldVal = nan;
				return;
			end
			oldVal = obj.profile.setParameterValue(param, newVal);
		end

		function notifySync(obj)
			%NOTIFYSYNC Called after parameters have been sent to CADET
			%   NOTIFYSYNC() resets the HASCHANGED property.

			obj.notifySync@InletModel();
			obj.profile.notifySync();
		end

		function val = get.constant(obj)
			val = obj.profile.constant;
		end

		function set.constant(obj, val)
			obj.profile.constant = val;
			obj.hasChanged = true;
		end

		function val = get.linear(obj)
			val = obj.profile.linear;
		end

		function set.linear(obj, val)
			obj.profile.linear = val;
			obj.hasChanged = true;
		end

		function val = get.quadratic(obj)
			val = obj.profile.quadratic;
		end

		function set.quadratic(obj, val)
			obj.profile.quadratic = val;
			obj.hasChanged = true;
		end

		function val = get.cubic(obj)
			val = obj.profile.cubic;
		end

		function set.cubic(obj, val)
			obj.profile.cubic = val;
			obj.hasChanged = true;
		end

	end

	methods (Access = 'protected')

		function loadobjInternal(obj, S)
			obj.loadobjInternal@InletModel(S);

			ctor = str2func([S.profileClass '.loadobj']);
			obj.profile = ctor(S.profile);
		end

	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = PiecewiseCubicPolyInlet();
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
