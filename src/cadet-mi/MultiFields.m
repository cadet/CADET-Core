
classdef MultiFields
	%MultiFields Helper class collecting methods for string-indexed struct fields
	%   String-indexed struct fields are Matlab structures of the form
	%     s.field_000, s.field_001, ..., s.field_00N.
	%   This class offers functions for reading them into a matrix, for writing them
	%   from a matrix, for removing them from the structure completely, and for
	%   extracting them (taking gaps into account) into a matrix.
	
	% Copyright: (C) 2008-2024 The CADET Authors
	%            See the license note at the end of the file.

	methods (Static)

		function val = read(data, prefix)
			%READ Reads a consecutively numbered string-indexed field into a matrix
			%   VAL = READ(DATA, PREFIX) reads the field PREFIX_000, PREFIX_001, etc.
			%   from the struct DATA into the matrix VAL. Each column in VAL contains
			%   the linearized data of a field. The field is supposed to start with
			%   PREFIX_000 and the process terminates once the next index is not found.
			%
			% See also MULTIFIELDS.EXTRACT

			% Count number of items
			numFields = 0;
			while isfield(data, sprintf('%s_%03d', prefix, numFields))
				numFields = numFields + 1;
			end

			% Allocate memory
			val = zeros(numel(data.([prefix '_000'])), numFields + 1);
			for i = 0:numFields
				val(:, i + 1) = data.(sprintf('%s_%03d', prefix, i));
			end
		end

		function data = write(data, prefix, val)
			%WRITE Writes a matrix to a string-indexed field
			%   DATA = WRITE(DATA, PREFIX, VAL) writes the matrix VAL to the string-indexed
			%   field PREFIX_000, PREFIX_001, etc. in the struct DATA. Each column of VAL
			%   is turned into a struct field.

			for i = 1:size(val, 2)
				data.(sprintf('%s_%03d', prefix, i-1)) = val(:, i);
			end
		end

		function data = remove(data, prefix)
			%REMOVE Complettely removes a string-indexed field from a struct
			%   DATA = REMOVE(DATA, PREFIX) Removes the string-indexed field PREFIX_000,
			%   PREFIX_001, etc. from the struct DATA. The process stops if the next field
			%   is not found (i.e., gaps in the numbering are not supported).

			i = 0;
			while isfield(data, sprintf('%s_%03d', prefix, i))
				data = rmfield(data, sprintf('%s_%03d', prefix, i));
				i = i + 1;
			end
		end

		function val = extract(data, names, baseName)
			%EXTRACT Extracts a string-indexed into a matrix supporting gaps in the numbering
			%   VAL = EXTRACT(DATA, NAMES, BASENAME) uses the fieldnames (as returned by 
			%   fieldnames(DATA)) given in the cell array NAMES to extract a string-indexed field
			%   BASENAME_000, BASENAME_001, etc. The data is linearized and packed into a column
			%   of the matrix VAL. The columns of VAL correspond to the fields in an ascending
			%   ordering (i.e., if we have BASENAME_000, BASENAME_002, BASENAME_003, the columns
			%   represent those fields in exactly this ordering ignoring the missing BASENAME_001).
			%
			% See also MULTIFIELDS.READ
			
			len = length(baseName);

			% Find datasets and convert to component indices
			idx = find(cellfun(@(x) (length(x) >= len) && strcmp(x(1:len), baseName), names));
			if isempty(idx)
				val = [];
				return;
			end

			idx = arrayfun(@(x) str2double(names{x}(len+2:end)), idx);
			idx = idx(~isnan(idx));
			idx = sort(idx);
			if isempty(idx)
				val = [];
				return;
			end
			
			% Preallocate
			dset = sprintf('%s_%03d', baseName, idx(1));
			val = zeros(length(data.(dset)), length(idx));

			for i = 1:length(idx)
				dset = sprintf('%s_%03d', baseName, idx(i));
				val(:, i) = data.(dset);
			end
		end

	end
end

% =============================================================================
%  CADET
%  
%  Copyright (C) 2008-2024: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. obj program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies obj distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
