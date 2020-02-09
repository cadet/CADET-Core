
classdef LinearInterpolationExtFun < ExternalFunction
	%LinearInterpolationExtFun Linear interpolating external function
	%   This external function takes into account time and axial position in the
	%   column, but ignores radial position in the bead. A quantity of interest
	%   is measured at the column outlet and recorded in a (time, value) like list.
	%   It is assumed that this quantity is transported inside the column with a
	%   known velocity. Based on the current time and the axial position inside
	%   the column, the correct time interval (subject to the velocity) of the
	%   data points is chosen and the measurements used for linear interpolation.
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
	% See also EXTERNALFUNCTION, MODEL
	
	% Copyright: (C) 2008-2020 The CADET Authors
	%            See the license note at the end of the file.
	
	properties (Constant)
		name = 'LINEAR_INTERP_DATA'; % Name of the external function according to CADET file format specs
	end

	properties (Dependent, Transient)
		% Profile values at time points [-]
		profile;
		DATA;
		% Time points in [s]
		time;
		TIME;
		% Normalized transport velocity in [1 / s]
		velocity;
		VELOCITY;
	end

	methods

		function obj = LinearInterpolationExtFun()
			%LINEARINTERPOLATIONEXTFUN Constructs a linear interpolating external function object

			obj = obj@ExternalFunction();
		end

		function res = validate(obj, sectionTimes)
			%VALIDATE Validates the configuration of the external function
			%   RES = VALIDATE(SECTIONTIMES) uses the section times of the simulator given in
			%   SECTIONTIMES to validate the external function object. Returns true in RES if
			%   everything is fine and false otherwise.

			validateattributes(obj.time, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'increasing', 'finite', 'real'}, '', 'time');
			validateattributes(obj.profile, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'data');
			validateattributes(obj.velocity, {'double'}, {'scalar', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'velocity');

			res = obj.validate@ExternalFunction(sectionTimes);
		end

		function val = evaluate(obj, t, z, r)
			%EVALUATE Evaluates the external function
			%   VAL = EVALUATE(T, Z, R) evaluates the external function at time T,
			%   normalized axial position Z (in [0, 1]), and normalized radial position
			%   R (in [0, 1]).

			validateattributes(t, {'double'}, {'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'time');
			validateattributes(z, {'double'}, {'nonempty', 'size', size(t), '>=', 0.0, '<=', 1.0, 'finite', 'real'}, '', 'z');

			% Calculate transformed time point
			transT = (1.0 - z(:)) ./ obj.velocity + t(:);

			% Cache values
			val = zeros(numel(t), 1);
			dataTime = obj.time;
			dataValues = obj.profile;
			
			% Determine index of section by binning
			if (exist('discretize') == 2)
				idx = discretize(transT, dataTime);
			else
				idx = nan(numel(transT), 1);
				for i = 1:numel(transT)
					if (transT(i) <= dataTime(1)) || (transT(i) >= dataTime(end))
						continue;
					end
					idx(i) = find((transT(i) >= dataTime(1:end-1)) & (transT(i) <= dataTime(2:end)), 1, 'first');
				end
			end

			% Extrapolate with constant values outside of time domain
			val(transT <= dataTime(1)) = dataValues(1);
			val(transT >= dataTime(end)) = dataValues(end);

			% Perform linear interpolation in the middle
			mask = isnan(idx);
			idx(mask) = [];
			val(~mask) = dataValues(idx) + (dataValues(idx + 1) - dataValues(idx)) .* (transT(~mask).' - dataTime(idx)) ./ (dataTime(idx+1) - dataTime(idx));

			val = reshape(val, size(t));
		end


		function val = get.time(obj)
			val = obj.data.TIME;
		end

		function set.time(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'increasing', 'finite', 'real'}, '', 'time');
			obj.data.TIME = val;
			obj.hasChanged = true;
		end

		function val = get.profile(obj)
			val = obj.data.DATA;
		end

		function set.profile(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', 'finite', 'real'}, '', 'profile');
			obj.data.DATA = val;
			obj.hasChanged = true;
		end

		function val = get.velocity(obj)
			val = obj.data.VELOCITY;
		end

		function set.velocity(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', '>', 0.0, 'finite', 'real'}, '', 'velocity');
			obj.data.VELOCITY = val;
			obj.hasChanged = true;
		end


		function val = get.TIME(obj)
			val = obj.time;
		end
		function set.TIME(obj, val)
			obj.time = val;
		end
		function val = get.DATA(obj)
			val = obj.profile;
		end
		function set.DATA(obj, val)
			obj.profile = val;
		end
		function val = get.VELOCITY(obj)
			val = obj.velocity;
		end
		function set.VELOCITY(obj, val)
			obj.velocity = val;
		end

	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = LinearInterpolationExtFun();
				obj.loadobjInternal(S);
			end
		end
		
	end
end

% =============================================================================
%  CADET - The Chromatography Analysis and Design Toolkit
%  
%  Copyright (C) 2008-2020: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. obj program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies obj distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
