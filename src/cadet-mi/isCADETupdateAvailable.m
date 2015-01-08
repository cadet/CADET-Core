function updateAvailable = isCADETupdateAvailable(quiet)
%ISCADETUPDATEAVAILABLE Checks GitHub for an updated CADET version
%
% Parameters:
%   - quiet: Optional. Set to true to suppress output to the command
%       window. Set to false to issue a notification if an update is
%       available. Defaults to false.
%
% Returns true if an update is available, otherwise false.
%
% Copyright: © 2008-2015 Eric von Lieres, Joel Andersson, Andreas Püttmann, Sebastian Schnittert, Samuel Leweke
%            See the license note at the end of the file.

    if (nargin <= 0) || isempty(quiet)
        quiet = false;
    end
    
    [updateAvailable, oldVer, newVer, link] = checkStable();
    
    if updateAvailable && ~quiet
        disp(['There is an updated version (' newVer ' > ' oldVer ') of CADET available.']);
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

    % Compare versions from left to right
    newVer = stableVersion;
    stableVersion = splitstring(stableVersion, '.');
    
    for i = 1:min(length(stableVersion), length(version))
        if str2double(stableVersion(i)) > str2double(version(i))
            updateAvailable = true;
            break;
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
    
    try
        sim = Simulator();
        [version, commit] = sim.getVersion();
    catch
        % Error in MEX interface
        disp('Error in CADET installation. Have you run installCADET.m?');
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
