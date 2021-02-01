function [p, factor] = extractParam(param, idx)
%EXTRACTPARAM Extracts a single parameter from a joined parameter
%   EXTRACTPARAM(PARAM) extracts the first parameter from a joined parameter PARAM consisting
%   of multiple single parameters.
%
%   EXTRACTPARAM(PARAM, IDX) extracts the single parameter at position / index IDX of
%   the joined parameter PARAM.
%
%   P = EXTRACTPARAM(...) returns the single parameter.
%
%   [P, FACTOR] = EXTRACTPARAM(...) also returns the linear factor that determines the contribution
%   of the single parameter to the joined parameter sensitivity.
%
% See also MAKESENSITIVITY.

% Copyright: (C) 2008-2020 The CADET Authors
%            See the license note at the end of the file.

	if (nargin <= 1) || isempty(idx)
		idx = 1;
	end

	validateattributes(param, {'struct'}, {'nonempty', 'numel', 1}, mfilename(), 'param');
	validateattributes(idx, {'double'}, {'nonempty', 'scalar'}, mfilename(), 'idx');

	p.SENS_UNIT = param.SENS_UNIT(idx);
	p.SENS_NAME = param.SENS_NAME(idx);
	p.SENS_COMP = param.SENS_COMP(idx);
	p.SENS_PARTYPE = param.SENS_PARTYPE(idx);
	p.SENS_REACTION = param.SENS_REACTION(idx);
	p.SENS_BOUNDPHASE = param.SENS_BOUNDPHASE(idx);
	p.SENS_SECTION = param.SENS_SECTION(idx);

	if isfield(param, 'SENS_ABSTOL')
		p.SENS_ABSTOL = param.SENS_ABSTOL;
	end
	if isfield(param, 'SENS_FACTOR')	
		p.SENS_FACTOR = param.SENS_FACTOR(idx);
		factor = p.SENS_FACTOR;
	else
		factor = 1.0;
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
