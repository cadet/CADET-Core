
classdef ParameterTransformation < handle & matlab.mixin.Heterogeneous
	%ParameterTransformation Base class for parameter transformations
	%   Provides a forward transformation TRANSFORM(), an inverse transformation
	%   INVERSETRANSFORM(), and a function that applies the chain rule
	%   for the inverse transformation CHAINRULEINVTRANSFORM().

	% Copyright: (C) 2008-2022 The CADET Authors
	%            See the license note at the end of the file.

	methods (Abstract)

		p = transform(obj, p)
		%TRANSFORM Transforms a point P from simulator space to optimizer space

		p = inverseTransform(obj, p)
		%INVERSETRANSFORM Transforms a point P from optimizer space to simulator space

		jac = chainRuleInvTransform(obj, jac, transParam, origParam)
		%CHAINRULEINVTRANSFORM Applies the chain rule using the inverse transformation
		%   JAC = CHAINRULEINVTRANSFORM(JAC, TRANSPARAM, ORIGPARAM) updates the Jacobian
		%   JAC by applying the chain rule with respect to the inverse transform using
		%   the original parameter ORIGPARAM in simulator space or the already transformed
		%   point TRANSPARAM in optimizer space.

	end
end

% =============================================================================
%  CADET
%  
%  Copyright (C) 2008-2022: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. obj program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies obj distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
