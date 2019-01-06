function sens = makeSensitivity(unitOpIdx, name, comp, parType, reaction, boundPhase, section, absTol, factors)
%MAKESENSITIVITY Creates a struct that represents a parameter sensitivity
%   SENS = MAKESENSITIVITY(UNITOPIDX, NAME, COMP, PARTYPE, REACTION, BOUNDPHASE, SECTION) creates a parameter
%   that can be used for computing parameter sensitivities. Multiple single parameters are
%   joined to a combined parameter. UNITOPIDX is a cell array with unit operation ids for each
%   single parameter. NAME is a cell array with the name of each single parameter. COMP is a
%   vector with 0-based component indices, PARTYPE is a vector with 0-based particle type indices,
%   REACTION is a vector with 0-based reaction indices, BOUNDPHASE is a vector with 0-based bound
%   phase indices, and SECTION is a vector with 0-based section indices. If a single parameter does
%   not depend on component, reaction, bound phase, or section, the respective index is set to -1.
%   All single parameter sensitivities are summed up to form the joined sensitivity.
%   Returns the joined parameter as a struct of arrays with the fields SENS_NAME, SENS_UNIT,
%   SENS_COMP, SENS_PARTYPE, SENS_REACTION, SENS_SECTION, SENS_BOUNDPHASE, SENS_FACTOR. The field
%   SENS_ABSTOL is a scalar in the struct.
%
%   SENS = MAKESENSITIVITY(..., ABSTOL) additionally sets the absolute error tolerance for the
%   sensitivity system (scalar).
%
%   SENS = MAKESENSITIVITY(..., ABSTOL, FACTORS) also specifies the linear factors for combining the
%   single parameter sensitivities into the joined sensitivity. By default, the factors are all
%   1.
%
% See also EXTRACTPARAM.

% Copyright: (C) 2008-2018 The CADET Authors
%            See the license note at the end of the file.

	if ~iscell(name) && ischar(name)
		name = {name};
	end

	validateattributes(name, {'cell'}, {'nonempty', 'vector'}, mfilename(), 'name');
	validateattributes(unitOpIdx, {'numeric'}, {'>=', -1, 'nonempty', 'vector', 'numel', numel(name)}, mfilename(), 'unitOpIdx');
	validateattributes(comp, {'numeric'}, {'>=', -1, 'nonempty', 'vector', 'numel', numel(unitOpIdx)}, mfilename(), 'comp');
	validateattributes(parType, {'numeric'}, {'>=', -1, 'nonempty', 'vector', 'numel', numel(unitOpIdx)}, mfilename(), 'parType');
	validateattributes(reaction, {'numeric'}, {'>=', -1, 'nonempty', 'vector', 'numel', numel(unitOpIdx)}, mfilename(), 'reaction');
	validateattributes(boundPhase, {'numeric'}, {'>=', -1, 'nonempty', 'vector', 'numel', numel(unitOpIdx)}, mfilename(), 'boundPhase');
	validateattributes(section, {'numeric'}, {'>=', -1, 'nonempty', 'vector', 'numel', numel(unitOpIdx)}, mfilename(), 'section');

	if (nargin >= 8) && ~isempty(absTol)
		validateattributes(absTol, {'numeric'}, {'>', 0, 'nonempty', 'scalar'}, mfilename(), 'absTol');
	else
		absTol = -1.0;
	end

	if (nargin >= 9) && ~isempty(factors)
		validateattributes(factors, {'numeric'}, {'nonempty', 'vector', 'nonzero', 'numel', numel(unitOpIdx)}, mfilename(), 'factors');
	else
		factors = ones(size(unitOpIdx));
	end

	sens.SENS_UNIT = int32(unitOpIdx);
	sens.SENS_NAME = name;
	sens.SENS_COMP = int32(comp);
	sens.SENS_PARTYPE = int32(parType);
	sens.SENS_REACTION = int32(reaction);
	sens.SENS_BOUNDPHASE = int32(boundPhase);
	sens.SENS_SECTION = int32(section);
	sens.SENS_ABSTOL = absTol;
	sens.SENS_FACTOR = factors;
end

% =============================================================================
%  CADET - The Chromatography Analysis and Design Toolkit
%  
%  Copyright (C) 2008-2018: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. obj program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies obj distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
