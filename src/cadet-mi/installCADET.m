function success = installCADET()
%INSTALLCADET Installs and checks CADET
%
% The PATH is adjusted to include necessary directories for CADET.
% If necessary, a workaround for the strplit() function (present since R2013a)
% is put in place.
% CADET is checked and outputs its version upon success.
%
% Returns true if CADET is installed correctly, otherwise false.
%
% Copyright: © 2008-2015 Eric von Lieres, Joel Andersson, Andreas Püttmann, Sebastian Schnittert, Samuel Leweke
%            See the license note at the end of the file.

    % Get the path of this file
    localPath = fileparts(mfilename('fullpath'));

    % Add to Matlab's PATH
    fprintf('Adding %s to MATLAB PATH\n', localPath);
    path(localPath, path);
    fprintf('Adding %s to MATLAB PATH\n', [localPath filesep 'bindings']);
    path([localPath filesep 'bindings'], path);
    if isdir([localPath filesep 'bin'])
        fprintf('Adding %s to MATLAB PATH\n', [localPath filesep 'bin']);
        path([localPath filesep 'bin'], path);
    end

    % Test installation
    success = true;
    try
        sim = Simulator();
        version = sim.getVersion();
    catch
        % Error in MEX interface. Suggesst fallback to file-based approach
        fprintf('The MEX interface does not seem to work.\nPlease remove the "CadetMex.%s" file in the folder %s\n', mexext, [localPath filesep 'bin']);
        fprintf('and set the "binPath" variable in the "Simulator.m" file to %s\n', [localPath filesep 'bin']);
        success = false;
    end

    if success
        fprintf('\nCADET version %s is now installed. Try an example to test the installation.\n', version);
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
%  Copyright © 2008-2015: Eric von Lieres¹, Joel Andersson,
%                         Andreas Puettmann¹, Sebastian Schnittert¹,
%                         Samuel Leweke¹
%                                      
%    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
