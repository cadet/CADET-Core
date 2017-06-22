function splitted = splitstring(str, delim)
%SPLITSTRING Splits a string using a given delimiter
%   SPLITTED = SPLITSTRING(STR, DELIM) splits the string STR into parts using
%   the delimiter DELIM. The result is a cell array with the splitted string
%   parts.

% Copyright: (C) 2008-2017 The CADET Authors
%            See the license note at the end of the file.

    splitted = textscan(str,'%s','delimiter',delim);
    splitted = splitted{1};
end

% =============================================================================
%  CADET - The Chromatography Analysis and Design Toolkit
%  
%  Copyright (C) 2008-2017: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
