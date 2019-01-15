
classdef PiecewiseCubicPolyExtFun < ExternalFunction
	%PiecewiseCubicPolyExtFun External function based on piecewise cubic polynomials
	%   This external function takes into account time and axial position in the
	%   column, but ignores radial position in the bead. A quantity of interest
	%   is modeled by piecewise cubic polynomials at the column outlet. It is
	%   assumed that this quantity is transported inside the column with a known
	%   velocity. Based on the current time and the axial position inside the
	%   column, the correct position in the profile is determined and the function
	%   evaluated.
	%
	%   The coordinate system of the external profile begins at the column outlet
	%   and points backward to the column inlet.
	%
	%         1 / velocity                                                 0
	%      t <---|---------------------------------------------------------|
	%             _________________________________________________________
	%            |                                                         |
	%     Inlet  |  =>               Column                            =>  | Outlet
	%            |_________________________________________________________|
	%
	%   The sections of this piecewise cubic polynomial are independent of the
	%   simulator time sections. Note that the simulator time sections have to take
	%   discontinuities of the external functions into account.
	%
	% See also EXTERNALFUNCTION, MODEL
	
	% Copyright: (C) 2008-2019 The CADET Authors
	%            See the license note at the end of the file.
	
	properties (Constant)
		name = 'PIECEWISE_CUBIC_POLY'; % Name of the external function according to CADET file format specs
	end

	properties (Dependent, Transient)
		% Constant coefficient for each section [-];
		constant;
		CONST_COEFF;
		% Linear coefficient for each section [1 / s];
		linear;
		LIN_COEFF;
		% Quadratic coefficient for each section [1 / s^2];
		quadratic;
		QUAD_COEFF;
		% Cubic coefficient for each section [1 / s^3];
		cubic;
		CUBE_COEFF;
		% Normalized transport velocity in [1 / s]
		velocity;
		VELOCITY;
		% Time sections for polynomial pieces [s]
		sectionTimes;
		SECTION_TIMES;
	end

	methods

		function obj = PiecewiseCubicPolyExtFun()
			%PIECEWISECUBICPOLYEXTFUN Constructs a piecewise cubic polynomial external function object

			obj = obj@ExternalFunction();
		end

		function res = validate(obj, sectionTimes)
			%VALIDATE Validates the configuration of the external function
			%   RES = VALIDATE(SECTIONTIMES) uses the section times of the simulator given in
			%   SECTIONTIMES to validate the external function object. Returns true in RES if
			%   everything is fine and false otherwise.

			validateattributes(obj.sectionTimes, {'double'}, {'nonnegative', 'vector', 'increasing', 'nonempty', 'finite', 'real'}, '', 'sectionTimes');
			if numel(obj.sectionTimes) < 2
				error('CADET:invalidConfig', 'Expected sectionTimes to have at least 2 elements.');
			end

			nSections = numel(obj.sectionTimes) - 1;

			validateattributes(obj.constant, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nSections}, '', 'constant');
			validateattributes(obj.linear, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nSections}, '', 'linear');
			validateattributes(obj.quadratic, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nSections}, '', 'quadratic');
			validateattributes(obj.cubic, {'double'}, {'vector', 'nonempty', 'finite', 'real', 'numel', nSections}, '', 'cubic');
			validateattributes(obj.velocity, {'double'}, {'scalar', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'velocity');

			res = obj.validate@ExternalFunction(sectionTimes);
		end

		function pp = makePiecewisePoly(obj, diffOrder)
			%MAKEPIECEWISEPOLY Assembles a Matlab piecewise polynomial that can be used by ppval
			%   PP = MAKEPIECEWISEPOLY() returns a piecewise polynomial in Matlab format that
			%   mimics the PiecewiseCubicPolyExtFun and can be passed to PPVAL, UNMKPP, etc.
			%
			%   PP = MAKEPIECEWISEPOLY(DIFFORDER) returns the DIFFORDER'th derivative of the
			%   profile. DIFFORDER is the order of the derivative and defaults to 0.

			if nargin <= 1
				diffOrder = 0;
			end

			nSections = min(numel(obj.sectionTimes) - 1, numel(obj.constant));

			switch diffOrder
				case 0
					pp = mkpp(obj.sectionTimes(1:nSections + 1), [obj.cubic(:), obj.quadratic(:), obj.linear(:), obj.constant(:)]);
				case 1
					pp = mkpp(obj.sectionTimes(1:nSections + 1), [3 .* obj.cubic(:), 2 .* obj.quadratic(:), obj.linear(:)]);
				case 2
					pp = mkpp(obj.sectionTimes(1:nSections + 1), [6 .* obj.cubic(:), 2 .* obj.quadratic(:)]);
				case 3
					pp = mkpp(obj.sectionTimes(1:nSections + 1), [6 .* obj.cubic(:)]);
				otherwise
					pp = mkpp(obj.sectionTimes(1:nSections + 1), zeros(nSections + 1, 1));
			end
		end

		function val = evaluate(obj, t, z, r)
			%EVALUATE Evaluates the external function
			%   VAL = EVALUATE(T, Z, R) evaluates the external function at time T,
			%   normalized axial position Z (in [0, 1]), and normalized radial position
			%   R (in [0, 1]).

			validateattributes(t, {'double'}, {'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'time');
			validateattributes(z, {'double'}, {'nonempty', 'size', size(t), '>=', 0.0, '<=', 1.0, 'finite', 'real'}, '', 'z');

			% Calculate transformed time point
			transT = (1.0 - z) ./ obj.velocity + t;

			% Calculate output
			pp = obj.makePiecewisePoly();
			val = ppval(pp, transT);
		end

		function val = get.constant(obj)
			val = obj.data.CONST_COEFF;
		end

		function set.constant(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'constant');
			obj.data.CONST_COEFF = val;
			obj.hasChanged = true;
		end

		function val = get.linear(obj)
			val = obj.data.LIN_COEFF;
		end

		function set.linear(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'linear');
			obj.data.LIN_COEFF = val;
			obj.hasChanged = true;
		end

		function val = get.quadratic(obj)
			val = obj.data.QUAD_COEFF;
		end

		function set.quadratic(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'quadratic');
			obj.data.QUAD_COEFF = val;
			obj.hasChanged = true;
		end

		function val = get.cubic(obj)
			val = obj.data.CUBE_COEFF;
		end

		function set.cubic(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'cubic');
			obj.data.CUBE_COEFF = val;
			obj.hasChanged = true;
		end

		function val = get.velocity(obj)
			val = obj.data.VELOCITY;
		end

		function set.velocity(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'velocity');
			obj.data.VELOCITY = val;
			obj.hasChanged = true;
		end

		function val = get.sectionTimes(obj)
			val = obj.data.SECTION_TIMES;
		end

		function set.sectionTimes(obj, val)
			validateattributes(val, {'double'}, {'nonnegative', 'vector', 'increasing', 'nonempty', 'finite', 'real'}, '', 'sectionTimes');
			if numel(val) < 2
				error('CADET:invalidConfig', 'Expected sectionTimes to have at least 2 elements.');
			end
			obj.data.SECTION_TIMES = val;
			obj.hasChanged = true;
		end


		function val = get.CONST_COEFF(obj)
			val = obj.constant;
		end
		function set.CONST_COEFF(obj, val)
			obj.constant = val;
		end
		function val = get.LIN_COEFF(obj)
			val = obj.linear;
		end
		function set.LIN_COEFF(obj, val)
			obj.linear = val;
		end
		function val = get.QUAD_COEFF(obj)
			val = obj.quadratic;
		end
		function set.QUAD_COEFF(obj, val)
			obj.quadratic = val;
		end
		function val = get.CUBE_COEFF(obj)
			val = obj.cubic;
		end
		function set.CUBE_COEFF(obj, val)
			obj.cubic = val;
		end
		function val = get.VELOCITY(obj)
			val = obj.velocity;
		end
		function set.VELOCITY(obj, val)
			obj.velocity = val;
		end
		function val = get.SECTION_TIMES(obj)
			val = obj.sectionTimes;
		end
		function set.SECTION_TIMES(obj, val)
			obj.sectionTimes = val;
		end

	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = PiecewiseCubicPolyExtFun();
				obj.loadobjInternal(S);
			end
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
