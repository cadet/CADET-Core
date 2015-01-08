
classdef Helpers
    %HELPERS Contains various helper functions
    %
    % Copyright: © 2008-2015 Eric von Lieres, Joel Andersson, Andreas Püttmann, Sebastian Schnittert, Samuel Leweke
    %            See the license note at the end of the file.
    
    methods (Static, Access = 'public')

        function [exists, val] = getValueFromNestedStruct(strct, fieldPath)
            % getValueFromNestedStruct Fetches the value of a field in a (deeply) nested struct
            %
            % Parameters:
            %   - strct: (Nested) Struct
            %   - fieldPath: Cell array with field names
            %
            % Returns:
            %   - exists: True if the requested field exists, otherwise false
            %   - val: Value of the field or [] if it does not exist
            %
            % Example:
            %   strct.first.level.part.A = 7;
            %   getValueFromNestedStruct(strct, [{'first'}, {'level'}, {'part'}, {'A'}])
            %   % Returns exists = true, val = 7;
            
            if ~isfield(strct, fieldPath{1})
                exists = false;
                val = [];
                return;
            end
            
            found = false;
            while ~found
                if ~isfield(strct, fieldPath{1})
                    % Not found
                    break;
                end
                strct = strct.(fieldPath{1});
                
                if length(fieldPath) == 1
                    found = true;
                    val = strct;
                    exists = true;
                    break;
                end
                
                fieldPath = fieldPath(2:end);
            end
            
            if ~found
                exists = false;
                val = [];
            end
        end
        
        function strct = setValueToNestedStruct(strct, fieldPath, val)
            % setValueToNestedStruct Sets a value to a path in a nested struct
            %
            % Parameters:
            %   - strct: Struct to change
            %   - fieldPath: Path to the field to change (f.e., 'test.A.b')
            %   - val: Value to assign
            %
            % Returns: The changed struct
            
            access = struct('type', {}, 'subs', {});
            
            for i = 1:length(fieldPath)
                if ~isfield(strct, fieldPath{i})
                    % Create field if it does not exist
%                    strct.(fieldPath{i}) = [];
                end
                
                access(i).type = '.';
                access(i).subs = fieldPath{i};
            end
            
            % Write
            strct = subsasgn(strct, access, val);
        end
        
        function ok = checkScalar(val, name)
            ok = true;
            if ~isscalar(val)
                ok = false;
                disp(['Error: ' name ' has to be scalar']);
            end
        end

        function ok = checkPositiveVector(val, len, name)
            ok = true;
            if ~any(length(val(:)) == len)
                ok = false;
                disp(['Error: ' name ' has to be a vector of length ' num2str(len)]);
            end
            if any(val <= 0)
                ok = false;
                disp(['Error: Each entry of ' name ' has to be positive, i.e., greater than 0']);
            end
        end

        function ok = checkNonnegativeVector(val, len, name)
            ok = true;
            if ~any(length(val(:)) == len)
                ok = false;
                disp(['Error: ' name ' has to be a vector of length ' num2str(len)]);
            end
            if any(val < 0)
                ok = false;
                disp(['Error: Each entry of ' name ' has to be nonnegative, i.e., greater than or equal 0']);
            end
        end
        
        function ok = checkPositiveScalar(val, name)
            ok = true;
            if val <= 0
                ok = false;
                disp(['Error: ' name ' has to be positive, i.e., greater than 0']);
            end
            if ~isscalar(val)
                ok = false;
                disp(['Error: ' name ' has to be scalar']);
            end
        end
        
        function ok = checkNonnegativeScalar(val, name)
            ok = true;
            if val < 0
                ok = false;
                disp(['Error: ' name ' has to be positive, i.e., greater than or equal 0']);
            end
            if ~isscalar(val)
                ok = false;
                disp(['Error: ' name ' has to be scalar']);
            end
        end
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
