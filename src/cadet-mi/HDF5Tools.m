
classdef HDF5Tools < handle
	%HDF5Tools Provides tools and utilities for working with HDF5 files
	%   Matlab structs are written to HDF5 files and converted to groups and
	%   datasets. The process also works in the other direction, where an
	%   HDF5 file is read in and converted to Matlab structs and arrays.
	
	% Copyright: (C) 2008-2018 The CADET Authors
	%            See the license note at the end of the file.

	methods (Static)

		function struct2hdf(fileName, strct, path, nodelete, legacy_hdf5)
			%STRUCT2HDF Writes a (nested) struct to an HDF5 file
			%   STRUCT2HDF(FILENAME, STRCT) writes the given
			%   Matlab struct STRCT to the HDF5 file identified by FILENAME.
			%
			%   STRUCT2HDF(..., PATH) uses the root node given by PATH in the
			%   HDF5 file (defaults to '/').
			%
			%   STRUCT2HDF(..., PATH, NODELETE) if set to true, the file is
			%   deleted before it is written to (default). Otherwise, data is
			%   appended.
			%
			%   STRUCT2HDF(..., PATH, NODELETE, LEGACY_HDF5) uses legacy HDF5
			%   functions if set to true (default).
			
			% If not provided, set default input values
			if (nargin <= 4) || isempty(legacy_hdf5)
				legacy_hdf5 = true;
			end
			if (nargin <= 3) || isempty(nodelete)
				nodelete = 0;
			end
			if (nargin <= 2) || isempty(path)
				path = '';
			end
			
			% Remove file if existent
			if ~isempty(dir(fileName)) && ~nodelete
				delete(fileName);
			end
			
			fields = fieldnames(strct);
			
			% Loop over all fields in the struct
			for f = fields'
				newpath = [path '/' f{1}];
				curfield = strct.(f{1});
				
				% If current field is a struct, call this function recursively
				if isstruct(curfield)
					HDF5Tools.struct2hdf(fileName, curfield, newpath, 1, legacy_hdf5);
					
					% If current field is a string, write it
				elseif ischar(curfield)
					mode = 'append';
					if isempty(dir(fileName)); mode = 'overwrite'; end
					hdf5write(fileName, newpath, curfield, 'WriteMode', mode);
					
					% If current field is an integer, write it as int32
				elseif isinteger(curfield)
					if legacy_hdf5
						mode = 'append';
						if isempty(dir(fileName)); mode = 'overwrite'; end
						hdf5write(fileName, newpath, int32(curfield), 'WriteMode', mode);
					else
						if isempty(curfield)
							curfield = [int32(0)]; % Dummy
						end
						h5create(fileName, newpath, size(curfield), 'Datatype', 'int32');
						h5write(fileName, newpath, curfield);
					end
						
					% If current field is a double, write it as extendible dataset
				else % isdouble
					if legacy_hdf5
						mode = 'append';
						if isempty(dir(fileName)); mode = 'overwrite'; end
						hdf5write(fileName, newpath, curfield, 'WriteMode', mode);
					else
						h5create(fileName, newpath, [Inf Inf], 'Datatype', 'double', 'ChunkSize', size(curfield));
						h5write(fileName, newpath, curfield, [1 1], size(curfield));
					end
				end
			end
			
		end
		
		function h = hdf2struct(filename, path)
			%HDF2STRUCT Converts an HDF5 file to a nested Matlab struct
			%   H = HDF2STRUCT(FILENAME) reads the HDF5 file identified by FILENAME
			%   and returns it in the nested Matlab struct H.
			%
			%   H = HDF2STRUCT(FILENAME, PATH) uses the given root node PATH
			%   in the HDF5 file to read the data.
			
			if nargin <= 1 || isempty(path)
				path = '/'; 
			end;
			
			h = [];
			h.hdf = [];

			% Open the HDF5 file
			fid = H5F.open(filename,'H5F_ACC_RDONLY','H5P_DEFAULT');
			
			% Open the group in path
			gid = H5G.open(fid, path);
			
			% Iterate over every entry in the HDF5 file and call 'operate' on it
			[status, h] = H5O.visit(gid, 'H5_INDEX_NAME', 'H5_ITER_NATIVE', @operate, h);
			
			% Close the group
			H5G.close(gid);
			
			% Close the HDF5 file
			H5F.close(fid);

			assert(status ~= 1, 'CADET:hdf2structConversion', 'Error in hdf2struct.');
			h = h.hdf;
		end
	end
end

function [status, h] = operate(obj, name, h)
%OPERATE Converts every dataset into a Matlab struct
	
	% Quick return on root group
	if strcmpi(name, '.')
		status = 0;
		return
	end

	% Check for dataset or group by case - DIRTY HACK!
	% But checking object info was not successful :-(
	last_name = regexp(name, '\w*$', 'match', 'once');
	if strcmp(last_name, upper(last_name))
		typestr = 'Dataset';
	else
		typestr = 'Group';
	end

	switch typestr
		case 'Group'

			% Store group name
			h.group = strrep(name, '/', '.');
			
		case 'Dataset'
			
			% Extract the name of the dataset
			dset_name = regexp(name, '\w*$', 'match', 'once');
			
			dataset = [];
			try
				% Open the dataset
				dataset = H5D.open(obj, name, 'H5P_DEFAULT');

				% Read data and store in handles.hdf structure
				path = splitstring(h.group, '.');
				path = struct('type', repmat({'.'}, length(path) + 1, 1), 'subs', [path(:); {dset_name}]);
				h.hdf = subsasgn(h.hdf, path, H5D.read(dataset, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT'));
		   
				H5D.close(dataset);
			catch
				% Close dataset if still open
				if ~isempty(dataset)
					H5D.close(dataset);
				end
				
				% Return error and stop visiting
				status = 1;
				return;
			end
			
		otherwise
			error('CADET:hdf2structConversion', 'Non-reachable code has been reached');
	end
	status = 0;
end

% =============================================================================
%  CADET - The Chromatography Analysis and Design Toolkit
%  
%  Copyright (C) 2008-2018: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
