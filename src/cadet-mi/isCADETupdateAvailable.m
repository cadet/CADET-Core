function updateAvailable = isCADETupdateAvailable(quiet)
%ISCADETUPDATEAVAILABLE Checks GitHub for an updated CADET version
%   UPDATEAVAILABLE = ISCADETUPDATEAVAILABLE() checks GitHub for an updated
%   CADET version and returns true if one is available and false otherwise.
%   A notification is sent to the command window that informs the user whether
%   an update is available or if the checking process failed for some reason.
%
%   UPDATEAVAILABLE = ISCADETUPDATEAVAILABLE(QUIET) determines whether any
%   textual output is sent to the command window (QUIET = false, default) 
%   or not (QUIET = true).

% Copyright: (C) 2008-2020 The CADET Authors
%            See the license note at the end of the file.

	if (nargin <= 0) || isempty(quiet)
		quiet = false;
	end
	
	[updateAvailable, oldVer, newVer, link] = checkStable();
	
	if updateAvailable && ~quiet
		fprintf('There is an updated version (%s > %s) of CADET available.', newVer, oldVer);
		disp(['You can download the release suitable for your machine from ' link '.']);
	end
	if ~updateAvailable && ~quiet && isempty(newVer)
		disp(['Update check failed because the GitHub servers could not be reached.']);
	end
end

function [updateAvailable, oldVer, newVer, link] = checkStable()
%CHECKSTABLE Checks the tags on GitHub for a more recent CADET version
	
	updateAvailable = false;
	link = '<a href="https://github.com/modsim/CADET/releases">GitHub</a>';

	version = getInstalledVersion();
	if isempty(version)
		oldVer = '';
		newVer = '';
		return;
	end
	
	oldVer = version;
	version = splitstring(version, '.');
	
	% Get most recent stable version from GitHub
	stableVersion = [];
	try
		stableVersion = urlread('https://raw.githubusercontent.com/modsim/CADET/master/version.txt');
	end
	
	if isempty(stableVersion)
		newVer = '';
		return;
	end
	
	% Remove line breaks
	stableVersion = regexprep(stableVersion, '\r\n|\n|\r', '');

	% Compare versions from left to right
	newVer = stableVersion;
	stableVersion = splitstring(stableVersion, '.');
	
	for i = 1:min(length(stableVersion), length(version))
		if str2double(stableVersion(i)) > str2double(version(i))
			updateAvailable = true;
			break;
		elseif str2double(stableVersion(i)) < str2double(version(i))
			% Local version is newer than remote version
			return;
		end
	end
	
	% All digits line up, so check if there is an additional digit in the
	% GitHub version
	if (~updateAvailable) && (length(stableVersion) > length(version))
		updateAvailable = str2double(stableVersion(length(version)+1)) > 0.0;
	end
end

function [version, commit] = getInstalledVersion()
%GETINSTALLEDVERSION Returns the installed CADET version and commit hash

	version = [];
	commit = [];
	
	[version, commit] = Simulator.getVersion();
	if isempty(version)
		% No working interface available
		disp('Error in CADET installation. Have you run installCADET.m?');
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
