function [pp, dF, dD] = piecewiseCubicHermitePoly(x, f, d)
%PIECEWISECUBICHERMITEPOLY Uses piecewise cubic Hermite interpolation to interpolate a function and its derivative
%
%   Creates a piecewise cubic interpolating polynomial that simultaneously
%   interpolates function values and derivatives at some points X.
%
%   PP = PIECEWISECUBICHERMITEPOLY(X, F, D) constructs a piecewise cubic Hermite
%   interpolating polynomial using abscissae X, function values F at X and
%   function derivatives D at X. The piecewise polynomial PP is returned in
%   Matlab format.
%
%   [PP, DF] = PIECEWISECUBICHERMITEPOLY(...) also returns a cell array DF
%   with derivatives of PP with respect to each function value F. This
%   means that DF has as many cells as the vector F has entries and that
%   each entry of DF is a piecewise cubic polynomial containing the
%   derivative of the interpolant with respect to the corresponding
%   function value.
%
%   [PP, DF, DD] = PIECEWISECUBICHERMITEPOLY(...) also returns a cell array
%   DD with derivatives of PP with respect to each function derivative D.
%   This means that DD has as many cells as the vector D has entries and
%   that each entry of DD is a piecewise cubic polynomial containing the
%   derivative of the interpolant with respect to the corresponding
%   function derivative.
%
% See also PCHIP, SPLINE, PPVAL.

% Copyright: (C) 2008-2016 The CADET Authors
%            See the license note at the end of the file.

	% Length of each segment
	h = diff(x);

	% Construct interpolation on each segment
	coeffs = zeros(length(h), 4);
	for i = 1:length(h)
		curH = h(i);
		
		% Calculate Hermite interpolating polynomial on [0, curH] using
		% function values f(i), f(i+1) and function derivative d(i), d(i+1)
		coeffs(i, :) = [(2 * (f(i) - f(i+1)) / curH + d(i) + d(i+1)) / curH^2, ...
			(3 * (-f(i) + f(i+1)) / curH - 2  * d(i) - d(i+1)) / curH, ...
			d(i), f(i)];
	end
	
	pp = mkpp(x, coeffs);

	% Calculate derivative of PP wrt. F and D
	if nargout > 1
		dF = cell(length(f), 1);
		dD = cell(length(f), 1);
		
		% First point (only appears in first segment as left end)
		coeffs = zeros(length(h), 4);
		curH = h(1);
		
		coeffs(1, :) = [2 / curH^3, -3 / curH^2, 0, 1];
		dF{1} = mkpp(x, coeffs);
		
		coeffs(1, :) = [1 / curH^2, -2 / curH, 1, 0];
		dD{1} = mkpp(x, coeffs);
		
		% Inner points (appear in segment j [left end] and j-1 [right end])
		for j = 2:length(f)-1
			coeffs = zeros(length(h), 4);
			curH = h(j);
			
			coeffs(j-1, :) = [-2 / curH^3, 3 / curH^2, 0, 0];
			coeffs(j, :) = [2 / curH^3, -3 / curH^2, 0, 1];
			dF{j} = mkpp(x, coeffs);
			
			coeffs(j-1, :) = [1 / curH^2, -1 / curH, 0, 0];
			coeffs(j, :) = [1 / curH^2, -2 / curH, 1, 0];
			dD{j} = mkpp(x, coeffs);
		end
		
		% Last point (only appears in last segment as right end)
		coeffs = zeros(length(h), 4);
		curH = h(end);
		
		coeffs(end, :) = [-2 / curH^3, 3 / curH^2, 0, 0];
		dF{end} = mkpp(x, coeffs);
		
		coeffs(end, :) = [1 / curH^2, -1 / curH, 0, 0];
		dD{end} = mkpp(x, coeffs);
	end
	
end

% =============================================================================
%  CADET - The Chromatography Analysis and Design Toolkit
%  
%  Copyright (C) 2008-2016: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
