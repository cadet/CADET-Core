function data = switchStorageOrdering(data, dims)
%SWITCHSTORAGEORDERING Changes the storage order of an array from column-major (Matlab) to row-major (CADET) and vice versa
%   DATA = SWITCHSTORAGEORDERING(DATA) infers the dimensions of the given array DATA and
%   shuffles them around in order to switch between column-major (Matlab) and row-major (CADET).
%
%   DATA = SWITCHSTORAGEORDERING(DATA, DIMS) takes the given array DATA and dimensions DIMS (vector
%   as returned by calling size(DATA)) and shuffles them around in order to switch between 
%   column-major (Matlab) and row-major (CADET).

% Copyright: (C) 2008-2018 The CADET Authors
%            See the license note at the end of the file.

	if (nargin <= 1) || isempty(dims)
		dims = size(data);
	end
	reverseInds = length(dims):-1:1;
	data = permute(reshape(data(:), dims(reverseInds)), reverseInds);

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
