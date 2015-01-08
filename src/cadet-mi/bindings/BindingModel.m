
classdef BindingModel < handle
    %BINDINGMODEL Base class for all binding models
    %
    % Copyright: © 2008-2015 Eric von Lieres, Joel Andersson, Andreas Püttmann, Sebastian Schnittert, Samuel Leweke
    %            See the license note at the end of the file.
    
    properties
        name;       % Name of the isotherm as given in the file format specification
        parameters; % A cell array of structs which define parameters of this isotherm
    end
    
    methods
        
        function obj = BindingModel()
            %BINDINGMODEL Constructs an Isotherm object
            obj.name = 'BINDING';
            obj.parameters = [];
        end
                
        function ok = checkParameters(obj, params, nComponents)
            %CHECKPARAMETERS Checks the values of the given parameters for plausibility
            %
            % Parameters: 
            %   - params: Struct with parameters as fields
            %
            % Returns: True if the values are ok, otherwise false

            ok = true;
        end
        
        function params = convertParameters(obj, params)
            %CONVERTPARAMETERS Converts the values of the given parameters to the required types specified in the file format
            %
            % Parameters: 
            %   - params: Struct with parameters as fields
            %
            % Returns: The struct with converted parameter values

            % Nothing to do here
        end
        
        function [result, paramName] = hasParameter(obj, paramName)
            % hasParameter Checks whether the given paramter exists in the model
            %
            % Parameters:
            %   - paramName: Name of the parameter according to CADET's
            %       file format specs
            %
            % Returns: 
            %   - result: True if the parameter exists, otherwise false
            %   - paramName: Actual name of the parameter
            
            idx = find(cellfun(@(val) strcmpi(val.displayName, paramName) || strcmpi(val.name, paramName), obj.parameters));
            result = ~isempty(idx);
            if result
                paramName = obj.parameters{idx(1)}.name;
            else
                paramName = '';
            end
        end
        
    end
    
    methods (Access = 'protected', Static)
        
        function ok = checkNonnegativeVector(val, len, name)
            % checkNonnegativeVector Checks whether the given vector is nonnegative and has a certain length
            %
            % Parameters:
            %   - val: Vector to check
            %   - len: Valid length
            %   - name: Name of the setting for the error message
            %
            % Returns: True if the vector is valid, otherwise false
            
            ok = true;
            if length(val(:)) ~= len
                ok = false;
                disp(['Error: ' name ' has to be a vector of length ' num2str(len)]);
            end
            if any(val < 0.0)
                ok = false;
                disp(['Error: Each entry of ' name ' has to be nonnegative, i.e., greater than or equal 0']);
            end
        end
        
        function ok = checkPositiveVector(val, len, name)
            % checkPositiveVector Checks whether the given vector is positive and has a certain length
            %
            % Parameters:
            %   - val: Vector to check
            %   - len: Valid length
            %   - name: Name of the setting for the error message
            %
            % Returns: True if the vector is valid, otherwise false
            
            ok = true;
            if length(val(:)) ~= len
                ok = false;
                disp(['Error: ' name ' has to be a vector of length ' num2str(len)]);
            end
            if any(val <= 0.0)
                ok = false;
                disp(['Error: Each entry of ' name ' has to be positive, i.e., greater than 0']);
            end
        end
        
        function param = createParameter(name, displayName, perComponent)
            %CREATEPARAMETER Creates a default parameter structure with given base values
            %
            % Parameters:
            %   - name: Name of the parameter. Has to match the name in the file format specification.
            %   - displayName: Optional. Name of the parameter when displayed in Matlab (f.e., in a message). Default: Use name.
            %   - perComponent: Optional. Boolean which determines whether the parameter is a vector with nComponents entries. Default: true.
            %
            % Returns: The created parameter structure.

            if (nargin <= 1) || isempty(displayName)
                displayName = name;
            end

            if (nargin <= 2) || isempty(perComponent)
                perComponent = true;
            end

            param.name = name;
            param.displayName = displayName;
            param.perComponent = perComponent;
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
