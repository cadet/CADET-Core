function pp = ppder(pp)
%PPDER Construct derivative of piecewise polynomial
%
%   Calculates the derivative of a piecewise polynomial as a piecewise
%   polynomial.
%
%   PP = PPDER(PP) calculates the derivative of PP and returns it
%
% See also PCHIP, SPLINE, PPVAL, PPINT.

% Copyright: (C) 2008-2016 The CADET Authors
%            See the license note at the end of the file.

	[breaks, coefs] = unmkpp(pp);

	% Take derivative of each polynomial independently
	coefs = coefs(:, 1:end-1) .* repmat(size(coefs, 2)-1:-1:1, size(coefs, 1), 1);

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
