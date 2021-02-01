function newPoly = shiftPoly(poly, shift)
%SHIFTPOLY Shifts a polynomial to a different expansion point
%   NEWPOLY = SHIFTPOLY(POLY, SHIFT) shifts a polynomial 
%   p(x) = a_n * x^n + ... a_1 * x + a_0 to q(x) = p(x - y), where y is
%   given by SHIFT. The polynomial coefficients POLY follow the Matlab
%   convention (i.e., the highest exponent comes first). The resulting
%   polynomial q is returned as coefficients using the same convention.

% Copyright: (C) 2008-2020 The CADET Authors
%            See the license note at the end of the file.

	newPoly = zeros(size(poly));
	newPoly(end) = poly(end);
	
	for n = 2:length(poly)
		% Expand (x-y)^n to sum_{j=0}^n nchoosek(n,j) * (-y)^{n-j} * x^j
		newPoly(end-n+1:end) = newPoly(end-n+1:end) + poly(end-n+1) .* arrayfun(@(x) nchoosek(n-1,x), 0:n-1) .* (-shift).^([0:n-1]);
	end

end

% =============================================================================
%  CADET
%  
%  Copyright (C) 2008-2020: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. obj program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies obj distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
