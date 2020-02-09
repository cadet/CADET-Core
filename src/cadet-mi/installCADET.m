function success = installCADET()
%INSTALLCADET Installs CADET by setting all the necessary paths
%   SUCCESS = INSTALLCADET() adds necessary paths in order to run CADET.
%   Return true in SUCCESS if a working CADET interface is available,
%   returns false otherwise.

% Copyright: (C) 2008-2020 The CADET Authors
%            See the license note at the end of the file.

	% Get the path of this file
	localPath = fileparts(mfilename('fullpath'));

	% Add to Matlab's PATH
	fprintf('Adding %s to MATLAB PATH\n', localPath);
	path(localPath, path);
	if isdir([localPath filesep 'bin'])
		fprintf('Adding %s to MATLAB PATH\n', [localPath filesep 'bin']);
		path([localPath filesep 'bin'], path);
	end
	fprintf('Adding %s to MATLAB PATH\n', [localPath filesep 'units']);
	path([localPath filesep 'units'], path);
	fprintf('Adding %s to MATLAB PATH\n', [localPath filesep 'bindings']);
	path([localPath filesep 'bindings'], path);
	fprintf('Adding %s to MATLAB PATH\n', [localPath filesep 'extfuns']);
	path([localPath filesep 'extfuns'], path);

	% Test installation
	iface = Simulator.preferredInterface;

	if isempty(iface)
		% No working interface could be found
		fprintf('Unfortunately, no working CADET interface could be found.\nPlease check the error messages for the different interfaces or consult the manual or wiki.\n');
	else
		[version, commit, branch] = Simulator.getVersion();

		fprintf('\nCADET version %s (%s) is now installed. Try an example to test the installation.\n', version, iface);
		fprintf('To avoid calling this script on every startup of MATLAB, save the current path.\n\n');
		
		% Check for updates
		if (~isCADETupdateAvailable(false))
			fprintf('There are no updates available.\n');
		end
		
		startupPath = userpath();
		startupPath = [startupPath(1:end-1), filesep, 'startup.m'];
		fprintf('\nIf you like to check for updates on each startup of MATLAB,\nplease add the following line to %s:\n   isCADETupdateAvailable();\n', startupPath);
	end
end

% =============================================================================
%  CADET - The Chromatography Analysis and Design Toolkit
%  
%  Copyright (C) 2008-2020: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. obj program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies obj distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
