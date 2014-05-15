function splitted = splitstring( str, delim )
%SPLITSTRING Splits a string using a given delimiter
%
% Parameters:
%	- str: String to split
%	- delim: Delimiter
%
% Returns cell array with splitted parts.
%
% Copyright: © 2008-2014 Eric von Lieres, Joel Andersson, Andreas Püttmann, Sebastian Schnittert, Samuel Leweke
%            See the license note at the end of the file.

    splitted = textscan(str,'%s','delimiter',delim);
    splitted = splitted{1};
end

% =============================================================================
%  CADET - The Chromatography Analysis and Design Toolkit
%  
%  Copyright © 2008-2014: Eric von Lieres¹, Joel Andersson,
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
