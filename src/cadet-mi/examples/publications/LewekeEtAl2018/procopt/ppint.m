function pp = ppint(pp, c)
%PPINT Construct anti-derivative of piecewise polynomial
%
%   The constant of the anti-derivative is set as C. In order to calculate
%   the correct anti-derivative of the i-th piece, we need to add the sum
%   of the definite integrals of all previous pieces. This is easy, since
%   the constant of the previous polynomial contains the sum of definite
%   integrals of its previous ones, and so on.
%
%   PP = PPINT(PP) calculates the anti-derivative of PP and returns it. The
%   constant of the anti-derivative is set as 0.
%
%   PP = PPINT(PP, C) calculates the anti-derivative of PP and returns it.
%   The constant of the anti-derivative is taken as C.
%
% See also PCHIP, SPLINE, PPVAL, PPDER.

% Copyright: (C) 2008-2016 The CADET Authors
%            See the license note at the end of the file.

	if (nargin <= 1) || isempty(c)
		c = 0;
	end

	% Note that the polynomial pieces of splines are shifted to their left
	% boundary.

	[breaks, coefs] = unmkpp(pp);
	% Take anti-derivative of each polynomial independently
	coefs = [coefs ./ repmat(size(coefs, 2):-1:1, size(coefs, 1), 1), zeros(size(coefs, 1), 1)];
	% Calculate constant term of each piece
	coefs(1, end) = c;
	for i = 2:size(coefs, 1)
		% Add definite integral of previous polynomial (note the coordinate
		% shift) to the sum of all of its previous pieces (saved in its
		% constant term).
		coefs(i, end) = diff(polyval(coefs(i-1, :), breaks(i-1:i) - breaks(i-1))) + coefs(i-1, end);
	end
	% Rebuild piecewise polynomial
	pp = mkpp(breaks, coefs);
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
