classdef Helpers
    %HELPERS Contains various helper functions
    
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

