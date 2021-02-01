
classdef CustomParameterTransformation < ParameterTransformation
	%CustomParameterTransformation Custom parameter transformation using function handles
	%   Performs a custom parameter transformation that is driven by
	%   user-defined functions.
	%
	%   P = TRANSFORMFCN(P) applies the forward transformation from
	%   simulator to optimizer space to the parameter vector P.
	%
	%   P = INVTRANSFORMFCN(P) applies the inverse transformation from
	%   optimizer to simulator space to the parameter vector P.
	%
	%   JAC = CHAINRULEINVFCN(JAC, TRANSPARAM, ORIGPARAM) updates the
	%   Jacobian JAC by applying the chain rule with respect to the
	%   inverse transform using the original parameter ORIGPARAM in
	%   simulator space or the already transformed point TRANSPARAM
	%   in optimizer space.

	% Copyright: (C) 2008-2021 The CADET Authors
	%            See the license note at the end of the file.

	properties
		transformFcn; % Function handle to function mapping simulator -> optimizer space
		invTransformFcn; % Function handle to function mapping optimizer -> simulator space
		chainRuleInvFcn; % Function handle that applies chain rule to given Jacobian with respect to inverse transformation
	end

	methods

		function obj = CustomParameterTransformation(transform, invTransform, chainRuleInv)
			%CUSTOMPARAMETERTRANSFORMATION Creates an object of the custom parameter transformation
			%   CUSTOMPARAMETERTRANSFORMATION() creates an empty custom parameter
			%   transformation that is only operable if its properties are set later.
			%
			%   CUSTOMPARAMETERTRANSFORMATION(TRANSFORM, INVTRANSFORM, CHAINRULEINV) sets
			%   TRANSFORM as the forward transformation, INVTRANSFORM as its inverse, and
			%   CHAINRULEINV as the chain rule applying function.

			if (nargin >= 1) && ~isempty(transform)
				validateattributes(transform, {'function_handle'}, {}, '', 'transform');
				obj.transformFcn = transform;
			else
				obj.transformFcn = @(p) p;
			end

			if (nargin >= 2) && ~isempty(invTransform)
				validateattributes(invTransform, {'function_handle'}, {}, '', 'invTransform');
				obj.invTransformFcn = invTransform;
			else
				obj.invTransformFcn = @(p) p;
			end

			if (nargin >= 3) && ~isempty(chainRuleInv)
				validateattributes(chainRuleInv, {'function_handle'}, {}, '', 'chainRuleInv');
				obj.chainRuleInvFcn = chainRuleInv;
			else
				obj.chainRuleInvFcn = @(jac, transParam, origParam) jac;
			end
		end

		function p = transform(obj, p)
			%TRANSFORM Applies the forward transformation from simulator to optimizer space
			p = obj.transformFcn(p);
		end

		function p = inverseTransform(obj, p)
			%INVERSETRANSFORM Applies the backward transformation from optimizer to simulator space
			p = obj.invTransformFcn(p);
		end

		function jac = chainRuleInvTransform(obj, jac, transParam, origParam)
			%CHAINRULEINVTRANSFORM Applies the chain rule to the given Jacobian with respect to the inverse transformation
			%   JAC = CHAINRULEINVTRANSFORM(JAC, TRANSPARAM, ORIGPARAM) updates the Jacobian
			%   JAC by applying the chain rule with respect to the inverse transform using
			%   the original parameter ORIGPARAM in simulator space or the already transformed
			%   point TRANSPARAM in optimizer space.
			
			jac = obj.chainRuleInvFcn(jac, transParam, origParam);
		end

		function set.transformFcn(obj, val)
			if ~isempty(val)
				validateattributes(val, {'function_handle'}, {}, '', 'transformFcn');
				obj.transformFcn = val;
			else
				obj.transformFcn = [];
			end
		end

		function set.invTransformFcn(obj, val)
			if ~isempty(val)
				validateattributes(val, {'function_handle'}, {}, '', 'invTransformFcn');
				obj.invTransformFcn = val;
			else
				obj.invTransformFcn = [];
			end
		end

		function set.chainRuleInvFcn(obj, val)
			if ~isempty(val)
				validateattributes(val, {'function_handle'}, {}, '', 'chainRuleInvFcn');
				obj.chainRuleInvFcn = val;
			else
				obj.chainRuleInvFcn = [];
			end
		end
	end
end

% =============================================================================
%  CADET
%  
%  Copyright (C) 2008-2021: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. obj program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies obj distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
