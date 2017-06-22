
classdef LogParameterTransformation < ParameterTransformation
	%LogParameterTransformation Logarithmic parameter transformation
	%   Performs the transformation y = ln(x) or y = ln(-x) with automatic
	%   detection of negative parameter domains.

	% Copyright: (C) 2008-2017 The CADET Authors
	%            See the license note at the end of the file.

	properties
		signs; % Stores the signs of the parameters
		subset; % The parameter subset to transform
	end

	methods

		function obj = LogParameterTransformation(baseParams, paramSubset)
			%LOGPARAMETERTRANSFORMATION Creates an object of the logarithmic parameter transformation
			%   LOGPARAMETERTRANSFORMATION(BASEPARAMS) uses the given expected / typical parameter values
			%   BASEPARAMS to detect which parameters have positive or negative domains.
			%
			%   LOGPARAMETERTRANSFORMATION(BASEPARAMS, PARAMSUBSET) additionally specifies a subset of
			%   parameters to be transformed. The other parameters are left unchanged. PARAMSUBSET can
			%   either be a logical array of the same size as BASEPARAMS, or a vector with indices.

			narginchk(1, 2);

			if any(baseParams == 0)
				error('CADET:funcParamError', 'Expected non-zero parameters.');
			end
			obj.signs = sign(baseParams(:));

			if (nargin <= 1) || isempty(paramSubset)
				obj.subset = true(size(obj.signs));
			else
				if ~islogical(paramSubset)
					validateattributes(paramSubset, {'numeric'}, {'finite', 'real', '>=', 1, '<=', numel(obj.signs)}, '', 'paramSubset');

					obj.subset = false(size(obj.signs));
					obj.subset(paramSubset) = true;
				else
					validateattributes(paramSubset, {'logical'}, {'vector', 'numel', numel(obj.signs)}, '', 'paramSubset');
					obj.subset = paramSubset(:);
				end
			end
		end

		function p = transform(obj, p)
			%TRANSFORM Applies the forward transformation from simulator to optimizer space y = ln(x)

			s = size(p);
			p = p(:);
			p(obj.subset) = log(obj.signs(obj.subset) .* p(obj.subset));
			reshape(p, s);
		end

		function p = inverseTransform(obj, p)
			%INVERSETRANSFORM Applies the backward transformation from optimizer to simulator space x = exp(y) or x = -exp(y)

			s = size(p);
			p = p(:);
			p(obj.subset) = obj.signs(obj.subset) .* exp(p(obj.subset));
			reshape(p, s);
		end

		function jac = chainRuleInvTransform(obj, jac, transParam, origParam)
			%CHAINRULEINVTRANSFORM Applies the chain rule to the given Jacobian with respect to the inverse transformation
			%   JAC = CHAINRULEINVTRANSFORM(JAC, TRANSPARAM, ORIGPARAM) updates the Jacobian
			%   JAC by applying the chain rule with respect to the inverse transform using
			%   the original parameter ORIGPARAM in simulator space or the already transformed
			%   point TRANSPARAM in optimizer space.
			%
			%   Since for g( exp(y) ) we have J = J_g * exp(y) and for g( -exp(y) ) the Jacobian
			%   is updated as J = J_g * (-exp(y)).
			
			% While origParam contains exp(y) or -exp(y), transParam contains y
			jac(:, obj.subset) = jac(:, obj.subset) .* repmat(origParam(obj.subset).', size(jac, 1), 1);
		end

	end
end

% =============================================================================
%  CADET - The Chromatography Analysis and Design Toolkit
%  
%  Copyright (C) 2008-2017: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. obj program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies obj distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
