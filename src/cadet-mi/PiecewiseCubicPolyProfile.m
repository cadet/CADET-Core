
classdef PiecewiseCubicPolyProfile < handle
	%PiecewiseCubicPolyProfile Uses piecewise cubic polynomials as a concentration profile
	%   Piecewise cubic polynomials can be used as a cubic spline (if certain continuity 
	%   conditions are met at the interfaces of two adjacent pieces), but can also represent
	%   discontinuous data (such as step profiles).
	
	% Copyright: (C) 2008-2017 The CADET Authors
	%            See the license note at the end of the file.
	
	properties (Transient)
		hasChanged; % Determines whether this object has changed after the last synchronization with CADET (C++ layer)
	end

	properties
		% Vector with strictly increasing elements which represent the start and end
		% of each interval, optional (ignored if object is not used as standalone)
		breaks;
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

	properties (Dependent, Transient)
		continuity; % Vector of logicals that determine whether a section transition is continuous
	end

	properties (Constant, Transient)
		name = 'PIECEWISE_CUBIC_POLY';
	end

	methods
		
		function obj = PiecewiseCubicPolyProfile()
			%PIECEWISECUBICPOLYPROFILE Constructs a piecewise cubic polynomial profile

			obj.breaks = [];
			obj.hasChanged = true;
		end
		

		function res = validate(obj, sectionTimes, nComponents, nonNegative)
			%VALIDATE Validates the configuration of the model
			%   RES = VALIDATE(SECTIONTIMES, NCOMPONENTS) validates the configuration of the profile and returns true if
			%   the settings are valid, otherwise false. SECTIONTIMES is the sections time vector configured in the
			%   Simulator and NCOMPONENTS is the number of components of the unit operation model.
			%
			%   RES = VALIDATE([], NCOMPONENTS) same as above but uses the local BREAKS property for SECTIONTIMES.
			%
			%   RES = VALIDATE(..., NONNEGATIVE) determines whether the profile is required to be nonnegative (default).

			if (nargin <= 3) || isempty(nonNegative)
				nonNegative = true;
			end

			if isempty(sectionTimes)
				validateattributes(obj.breaks, {'double'}, {'>=', 0.0, 'vector', 'increasing', 'nonempty', 'finite', 'real'}, '', 'breaks');
				if numel(obj.breaks) < 2
					error('CADET:invalidConfig', 'Expected breaks to have at least 2 elements.');
				end				
				sectionTimes = obj.breaks;
			end

			sectionTimes = sectionTimes(:);
			nSections = numel(sectionTimes) - 1;

			validateattributes(obj.constant, {'double'}, {'2d', 'nonempty', 'finite', 'real', '>=', 0.0, 'numel', nComponents * nSections}, '', 'constant');
			validateattributes(obj.linear, {'double'}, {'2d', 'nonempty', 'finite', 'real', 'numel', nComponents * nSections}, '', 'linear');
			validateattributes(obj.quadratic, {'double'}, {'2d', 'nonempty', 'finite', 'real', 'numel', nComponents * nSections}, '', 'quadratic');
			validateattributes(obj.cubic, {'double'}, {'2d', 'nonempty', 'finite', 'real', 'numel', nComponents * nSections}, '', 'cubic');

			if nonNegative
				% Determine minimum value of each polynomial on every piece

				% Polynomial f(x) = a * (x-s)^3 + b * (x-s)^2 + c * (x-s) + d on [s, t]
				% Extremum occurs at boundaries f(s) = d, f(t), or at z given by
				% f'(z) = 3 * a * (z-s)^2 + 2 * b * (z-s) + c == 0
				% Solve for z-s to obtain
				% (z-s)^2 + 2 * b / (3*a) * (z-s) + c / (3*a) == 0
				% z-s = -b / (3*a) +- sqrt( b^2 / (9 * a^2) - c / (3 * a) ) for a ~= 0
				% z-s = -c / (2*b) for a == 0
				sqrtDet = sqrt(obj.quadratic.^2 - 3 .* obj.linear .* obj.cubic);
				extremaA = (obj.cubic ~= 0) .* (-obj.quadratic - sqrtDet) ./ (3 .* obj.cubic) + (obj.cubic == 0) .* (-obj.linear ./ (2 .* obj.quadratic));
				extremaB = (obj.cubic ~= 0) .* (-obj.quadratic + sqrtDet) ./ (3 .* obj.cubic) + (obj.cubic == 0) .* (-obj.linear ./ (2 .* obj.quadratic));

				% For each component ...
				for j = 1:size(obj.constant, 2)
					pp = mkpp(sectionTimes, [obj.cubic(:,j), obj.quadratic(:,j), obj.linear(:,j), obj.constant(:,j)]);

					% For each interval / piece ...
					for i = 1:nSections
						% Check boundaries
						if any(ppval(pp, [sectionTimes(i), sectionTimes(i) + eps, sectionTimes(i+1), sectionTimes(i+1) - eps]) < 0)
							error('CADET:invalidConfig', 'Expected a nonnegative inlet profile (section %d, component %d).', i, j);
						end
						% Check extrema if they are valid (real, finite, inside current interval)
						if isreal(extremaA(i,j)) && isfinite(extremaA(i,j)) && (extremaA(i,j) >= 0) && (extremaA(i,j) <= sectionTimes(i+1) - sectionTimes(i))
							if ppval(pp, extremaA(i,j) + sectionTimes(i)) < 0
								error('CADET:invalidConfig', 'Expected a nonnegative inlet profile (section %d, component %d).', i, j);
							end
						end
						if isreal(extremaB(i,j)) && isfinite(extremaB(i,j)) && (extremaB(i,j) >= 0) && (extremaB(i,j) <= sectionTimes(i+1) - sectionTimes(i))
							if ppval(pp, extremaB(i,j) + sectionTimes(i)) < 0
								error('CADET:invalidConfig', 'Expected a nonnegative inlet profile (section %d, component %d).', i, j);
							end
						end
					end
				end
			end

			res = true;
		end

		function res = assembleConfig(obj, res)
			%ASSEMBLECONFIG Assembles the configuration according to the CADET file format spec
			%   RES = assembleConfig() Returns the configuration in a struct.
			%   RES = assembleConfig(RES) Adds the configuration to the given struct RES and returns it.

			if nargin < 1
				res = [];
			end

			for i = 1:size(obj.constant, 1)
				curPoly = [];
				curPoly.CONST_COEFF = double(obj.constant(i, :));
				curPoly.LIN_COEFF = double(obj.linear(i, :));
				curPoly.QUAD_COEFF = double(obj.quadratic(i, :));
				curPoly.CUBE_COEFF = double(obj.cubic(i, :));
				res.(sprintf('sec_%03d', i - 1)) = curPoly;
			end
		end

		function pp = makePiecewisePoly(obj, sectionTimes, comp, diffOrder)
			%MAKEPIECEWISEPOLY Assembles a Matlab piecewise polynomial that can be used by ppval
			%   PP = MAKEPIECEWISEPOLY(SECTIONTIMES, COMP) returns a piecewise polynomial in Matlab
			%   format that mimics the PiecewiseCubicPolyProfile of one component and can be passed
			%   to PPVAL, UNMKPP, etc. For this single component, it is fully equivalent to the profile.
			%   SECTIONTIMES is a vector with section times (start, end) and COMP is the index of 
			%   the requested component (1-based).
			%
			%   PP = MAKEPIECEWISEPOLY([], COMP) same as above, but uses the local BREAKS property
			%   instead of SECTIONTIMES.
			%
			%   PP = MAKEPIECEWISEPOLY(..., DIFFORDER) returns the DIFFORDER'th derivative of the
			%   profile. DIFFORDER is the order of the derivative and defaults to 0.

			if isempty(sectionTimes)
				sectionTimes = obj.breaks;
			end

			validateattributes(sectionTimes, {'double'}, {'nonnegative', 'vector', 'increasing', 'nonempty', 'finite', 'real'}, '', 'sectionTimes');
			if numel(sectionTimes) < 2
				error('CADET:invalidConfig', 'Expected sectionTimes to have at least 2 elements.');
			end
			validateattributes(comp, {'numeric'}, {'>=', 1, '<=', size(obj.constant, 2), 'scalar', 'nonempty', 'finite', 'real'}, '', 'comp');

			if nargin <= 3
				diffOrder = 0;
			end

			nSections = min(numel(sectionTimes) - 1, size(obj.constant, 1));

			switch diffOrder
				case 0
					pp = mkpp(sectionTimes(1:nSections + 1), [obj.cubic(:, comp), obj.quadratic(:, comp), obj.linear(:, comp), obj.constant(:, comp)]);
				case 1
					pp = mkpp(sectionTimes(1:nSections + 1), [3 .* obj.cubic(:, comp), 2 .* obj.quadratic(:, comp), obj.linear(:, comp)]);
				case 2
					pp = mkpp(sectionTimes(1:nSections + 1), [6 .* obj.cubic(:, comp), 2 .* obj.quadratic(:, comp)]);
				case 3
					pp = mkpp(sectionTimes(1:nSections + 1), [6 .* obj.cubic(:, comp)]);
				otherwise
					pp = mkpp(sectionTimes(1:nSections + 1), zeros(nSections + 1, 1));
			end
		end

		function S = saveobj(obj)
			S = [];
			S.breaks = obj.breaks;
			S.constant = obj.constant;
			S.linear = obj.linear;
			S.quadratic = obj.quadratic;
			S.cubic = obj.cubic;			
		end
		
		function val = getParameterValue(obj, param)
			%GETPARAMETERVALUE Returns the value of a parameter
			%   VAL = GETPARAMETERVALUE(PARAM) returns the value of the parameter identified by PARAM
			%   If the parameter is not found, NaN is returned. PARAM is assumed to be struct with the
			%   fields SENS_SECTION, SENS_COMP, and SENS_NAME. The first two are 0-based indices and 
			%   the latter is a string (one of 'CONST_COEFF', 'LIN_COEFF', 'QUAD_COEFF', and 'CUBIC_COEFF').
			%
			% See also PIECEWISECUBICPOLYPROFILE.SETPARAMETERVALUE, MAKESENSITIVITY.

			section = param.SENS_SECTION + 1;
			comp = param.SENS_COMP + 1;
			if (section < 1) || (comp < 1)
				val = nan;
				return;
			end

			switch param.SENS_NAME
				case 'CONST_COEFF'
					val = obj.constant(section, comp);
				case 'LIN_COEFF'
					val = obj.linear(section, comp);
				case 'QUAD_COEFF'
					val = obj.quadratic(section, comp);
				case 'CUBIC_COEFF'
					val = obj.cubic(section, comp);
				otherwise
					val = nan;
			end
		end

		function oldVal = setParameterValue(obj, param, newVal)
			%SETPARAMETERVALUE Sets a parameter to a given value
			%   OLDVAL = SETPARAMETERVALUE(PARAM, NEWVAL) sets the parameter identified by PARAM
			%   to the value NEWVAL and returns the old value in OLDVAL. If the parameter is not
			%   found, NaN is returned. PARAM is assumed to be struct with the fields SENS_SECTION,
			%   SENS_COMP, and SENS_NAME. The first two are 0-based indices and the latter is a
			%   string (one of 'CONST_COEFF', 'LIN_COEFF', 'QUAD_COEFF', and 'CUBIC_COEFF').
			%
			% See also PIECEWISECUBICPOLYPROFILE.GETPARAMETERVALUE, MAKESENSITIVITY.

			section = param.SENS_SECTION + 1;
			comp = param.SENS_COMP + 1;
			if (section < 1) || (comp < 1)
				val = nan;
				return;
			end

			switch param.SENS_NAME
				case 'CONST_COEFF'
					oldVal = obj.constant(section, comp);
					obj.constant(section, comp) = newVal;
				case 'LIN_COEFF'
					oldVal = obj.linear(section, comp);
					obj.linear(section, comp) = newVal;
				case 'QUAD_COEFF'
					oldVal = obj.quadratic(section, comp);
					obj.quadratic(section, comp) = newVal;
				case 'CUBIC_COEFF'
					oldVal = obj.cubic(section, comp);
					obj.cubic(section, comp) = newVal;
				otherwise
					oldVal = nan;
			end
		end

		function notifySync(obj)
			%NOTIFYSYNC Called after parameters have been sent to CADET
			%   NOTIFYSYNC() resets the HASCHANGED property.

			obj.hasChanged = false;
		end

		function set.constant(obj, val)
			obj.constant = val;
			obj.hasChanged = true;
		end

		function set.linear(obj, val)
			obj.linear = val;
			obj.hasChanged = true;
		end

		function set.quadratic(obj, val)
			obj.quadratic = val;
			obj.hasChanged = true;
		end

		function set.cubic(obj, val)
			obj.cubic = val;
			obj.hasChanged = true;
		end

		function val = get.continuity(obj)
			if isempty(obj.breaks) || (numel(obj.breaks) - 1 ~= size(obj.constant, 1)) || ...
				(numel(obj.breaks) - 1 ~= size(obj.linear, 1)) || (numel(obj.breaks) - 1 ~= size(obj.quadratic, 1)) || ...
				(numel(obj.breaks) - 1 ~= size(obj.cubic, 1))

				val = [];
				return;
			end

			% Threshold for continuity detection (maximal jump height that is still considered continuous)
			threshold = 1e-12;

			nComp = min([size(obj.constant, 2), size(obj.linear, 2), size(obj.quadratic, 2), size(obj.cubic, 2)]);
			val = true(numel(obj.breaks) - 2, 1);

			% Check continuity for each transition
			for i = 1:numel(obj.breaks) - 2
				
				% Transition from section i to i + 1
				% Section i = [breaks(i), breaks(i+1)]
				% Section i+1 = [breaks(i+1), breaks(i+2)]
				
				% A transition is continuous if and only if the section
				% transition of each component is continuous
				
				for comp = 1:nComp
					% Piecewise polynomials are defined on intervals
					% T_i = [t_i, t_{i+1}] by the formula
					%   P_i(t) = c_n * (t - t_{i})^n + ... + c_1 * (t - t_{i}) + c_0
					% Thus, if the transition from section i to i+1 is
					% continuous, it holds that
					%   L_i( t_{i+1} - t_{i} ) = L_{i+1}( 0 ),
					% where L_i( t ) = c_n * t^n + ... + c_1 * t + c_0
					% is the unshifted polynomial.

					% Build left polynomial, i.e., on section i
					leftCoefs = [obj.cubic(i, comp), ...
						obj.quadratic(i, comp), ...
						obj.linear(i, comp), ...
						obj.constant(i, comp)];

					% Instead of taking the right polynomial, i.e., on
					% section i+1, we take only the absolute value of this
					% polynomial. This is valid since the right polynomial
					% would be evaluated at 0 which only retains the
					% absolute value.
					rightVal = obj.constant(i+1, comp);
				
					% Compare
					val(i) = (abs(polyval(leftCoefs, obj.breaks(i+1) - obj.breaks(i)) - rightVal) <= threshold) && val(i);
					if ~val(i)
						% Early out if a component is not continuous
						break;
					end
				end
			end
		end
	
	end

	methods (Access = 'protected')

		function loadobjInternal(obj, S)
			obj.hasChanged = false;
			obj.breaks = S.breaks;
			obj.constant = S.constant;
			obj.linear = S.linear;
			obj.quadratic = S.quadratic;
			obj.cubic = S.cubic;
		end

	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = PiecewiseCubicPolyProfile();
				obj.loadobjInternal(S);
			end
		end

		function obj = fromPiecewisePolynomials(pps)
			%FROMPIECEWISEPOLYNOMIALS Creates a PiecewiseCubicPolyProfile from a cell array of piecewise (cubic) polynomials
			%   OBJ = FROMPIECEWISEPOLYNOMIALS(PPS) turns the cell array of piecewise polynomials PPS into a
			%   PiecewiseCubicPolyProfile. Polynomial coefficients of higher than cubic monomials are ignored. It is
			%   assumed that the pieces of the different polynomials in PPS are the same.
			%
			% See also PIECEWISECUBICPOLYPROFILE.FROMUNIFORMDATA, PIECEWISECUBICPOLYPROFILE.FROMRESULT.

			obj = PiecewiseCubicPolyProfile();

			validateattributes(pps, {'cell'}, {'nonempty'}, '', 'pps');
			[time, coefs] = unmkpp(pps{1});
			obj.breaks = time;

			% TODO: Check breaks of other piecewise polynomials and merge if necessary

			nComp = numel(pps);
			nSec = length(time) - 1;

			obj.constant = zeros(nSec, nComp);
			obj.linear = zeros(nSec, nComp);
			obj.quadratic = zeros(nSec, nComp);
			obj.cubic = zeros(nSec, nComp);

			for i = 1:nComp
				[~, coefs] = unmkpp(pps{i});

				obj.constant(:,i) = coefs(:,end);
				if size(coefs, 2) >= 2
					obj.linear(:,i) = coefs(:,end-1);
				end
				if size(coefs, 2) >= 3
					obj.quadratic(:,i) = coefs(:,end-2);
				end
				if size(coefs, 2) >= 4
					obj.cubic(:,i) = coefs(:,end-3);
				end
			end
		end

		function obj = fromUniformData(time, concentrations, method)
			%FROMUNIFORMDATA Creates a PiecewiseCubicPolyProfile from a time series
			%   OBJ = FROMUNIFORMDATA(TIME, CONCENTRATIONS) creates a PiecewiseCubicPolyProfile from a
			%   time series with the time points given in the vector TIME and the data given in the matrix
			%   CONCENTRATIONS. It is assumed that each row in CONCENTRATIONS refers to a time point in
			%   TIME. Here, 'uniform' refers to the fact that all columns in CONCENTRATIONS refer to the
			%   same TIME vector. The data is interpolated by a piecewise cubic hermite interpolating 
			%   polynomial (Matlab's pchip function) that preserves the shape.
			%
			%   OBJ = FROMUNIFORMDATA(..., METHOD) uses the given METHOD (function handle) for creating
			%   the piecewise interpolating polynomial instead of pchip.
			%
			% See also PCHIP, PIECEWISECUBICPOLYPROFILE.FROMPIECEWISEPOLYNOMIAL, PIECEWISECUBICPOLYPROFILE.FROMRESULT.

			obj = PiecewiseCubicPolyProfile();

			if (nargin <= 2) || isempty(method)
				method = @pchip;
			end

			validateattributes(time, {'double'}, {'nonnegative', 'increasing', 'vector', 'nonempty', 'finite', 'real'}, '', 'time');
			if numel(time) <= 1
				error('CADET:funcParamError', 'Expected time to have at least two entries.');
			end
			validateattributes(concentrations, {'double'}, {'2d', 'size', [length(time), NaN], 'finite', 'real'}, '', 'concentrations');

			nComp = size(concentrations, 2);
			nSec = length(time) - 1;
			obj.breaks = time;

			obj.constant = zeros(nSec, nComp);
			obj.linear = zeros(nSec, nComp);
			obj.quadratic = zeros(nSec, nComp);
			obj.cubic = zeros(nSec, nComp);

			% Create interpolating spline for each component
			for i = 1:nComp
				pp = method(time, concentrations(:,i));
				[~, coefs] = unmkpp(pp);

				obj.constant(:,i) = coefs(:,4);
				obj.linear(:,i) = coefs(:,3);
				obj.quadratic(:,i) = coefs(:,2);
				obj.cubic(:,i) = coefs(:,1);
			end
		end

		function obj = fromResult(res, unitOpIdx)
			%FROMRESULT Creates a PiecewiseCubicPolyProfile from a simulation result
			%   OBJ = FROMRESULT(RES) uses the given simulation results RES to construct an interpolating
			%   piecewise cubic polynomial that and returns the resulting PiecewiseCubicPolyProfile.
			%   By default, the outlet of the unit operation with id 0 is taken as source of the profile.
			%
			%   OBJ = FROMRESULT(..., UNITOPIDX) uses the outlet of unit operation UNITOPIDX (0-based index)
			%   as source.
			%
			% See also PIECEWISECUBICPOLYPROFILE.FROMPIECEWISEPOLYNOMIAL, PIECEWISECUBICPOLYPROFILE.FROMUNIFORMDATA.
			
			if isempty(res.solution.time)
				error('CADET:funcParamError', 'Expected res to have a nonempty solution.time field.');
			end

			if (nargin <= 1) || isempty(unitOpIdx)
				unitOpIdx = 0;
			end
			validateattributes(unitOpIdx, {'numeric'}, {'nonnegative', 'scalar', 'nonempty', 'finite', 'real'}, '', 'unitOpIdx');

			if isempty(res.solution.time)
				error('CADET:funcParamError', 'Expected nonempty res.solution.time field.');
			end

			obj.breaks = res.solution.time;
			if (length(res.solution.outlet) >= unitOpIdx + 1)
				obj = PiecewiseCubicPolyProfile.fromUniformData(res.solution.time, max(0.0, res.solution.outlet{unitOpIdx+1}));
			else
				if isempty(res.solution.outlet)
					error('CADET:funcParamError', 'Expected nonempty res.solution.outlet field.');
				else
					error('CADET:funcParamError', 'Expected unitOpIdx to be less than %d.', length(res.solution.outlet));
				end
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
