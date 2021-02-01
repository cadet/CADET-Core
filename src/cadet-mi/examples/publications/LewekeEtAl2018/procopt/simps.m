function v = simps(x, f)
%SIMPS Calculates a definite integral using equidistant points and Simpson's rule
%
%   A quadrature rule for numerical integration is constructed using the
%   given equidistant points X. Simpson's rule exploits a quadratic
%   approximation of the integrand on slices of the interval to accomplish
%   fourth order convergence and integrate polynomials of up to third
%   degree exactly.
%
%   However, due to the fact that Simpson's rule requires an uneven number
%   of points to form slices of the domain that are later summed up, the
%   order is only guaranteed for an uneven number of points. In case of an
%   even number n, Simspon's rule is applied to the first n-1 points and a
%   special rule is employed for the last segment. Here, the two points to
%   the left of the segment are added to the calculation for maintaining
%   the order of convergence.
%
%   V = SIMPS(X) constructs a quadrature rule V for the equidistant points
%   X. X is a vector with abscissae and V is a vector with coefficients of
%   the linear form that computes the integral. In order to get the final
%   value, the dot-product of V and F(X) has to be calculated.
%   
%   V = SIMPS(X, F) calculates a definite integral using abscissae X and
%   corresponding function evaluations F that are both expected to be
%   vectors. The result is returned in V.
%
%   See also TRAPZ.

% Copyright: (C) 2008-2016 The CADET Authors
%            See the license note at the end of the file.

	if mod(numel(x), 2) == 0
		% Even number of points
		
		% Simpson rule on all but the last point (to get uneven number)		
		v = [ordinarySimpsonRule(x(1:end-1)); 0];
		
		% Add rule for last segment that uses points to the left of the
		% segment for maintaining order
		h = x(end) - x(end-1);
		
		v(end-3) = v(end-3) + h / 24;
		v(end-2) = v(end-2) - 5/24 * h;
		v(end-1) = v(end-1) + 19 / 24 * h;
		v(end) = 3 / 8 * h;
	else
		% Uneven number of points
		v = ordinarySimpsonRule(x);
	end

	if nargin > 1
		v = v.' * f(:);
	end
end

function v = ordinarySimpsonRule(x)
	% Uneven number of points
	v = zeros(length(x), 1);
	N = (numel(x) - 1) / 2;
	h = (x(end) - x(1)) / N;

	v(1) = 1;
	v(end) = 1;
	v(3:2:end-2) = 2;
	v(2:2:end-1) = 4;

	v = h / 6 .* v;
end

% =============================================================================
%  CADET
%  
%  Copyright (C) 2008-2016: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
