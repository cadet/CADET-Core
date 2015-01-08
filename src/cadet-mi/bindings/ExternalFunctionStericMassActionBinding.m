
classdef ExternalFunctionStericMassActionBinding < BindingModel
    %EXTERNALFUNCTIONSTERICMASSACTIONBINDING Represents the External Function Steric Mass Action isotherm by providing its parameters and helper functions
    % Provides the parameters of the External Function Steric Mass Action isotherm along with the methods to validate and convert them to
    % the types specified in the CADET file format.
    %
    % Copyright: © 2008-2015 Eric von Lieres, Joel Andersson, Andreas Püttmann, Sebastian Schnittert, Samuel Leweke
    %            See the license note at the end of the file.
    
    methods

        function obj = ExternalFunctionStericMassActionBinding()
            %EXTERNALFUNCTIONSTERICMASSACTIONBINDING Constructs an ExternalFunctionStericMassActionBinding object
            obj.name = 'EXTERNAL_STERIC_MASS_ACTION';
            obj.parameters = cell(12,1);

            obj.addParameterGroup('KA', 0, true);
            obj.addParameterGroup('KD', 4, true);
            obj.addParameterGroup('NU', 8, true);
            obj.addParameterGroup('SIGMA', 12, true);
            obj.addParameterGroup('LAMBDA', 16, false);
        end
        
        function ok = checkParameters(obj, params, nComponents)
            %CHECKPARAMETERS Checks the values of the given parameters for plausibility
            %
            % Parameters: 
            %   - params: Struct with parameters as fields
            %
            % Returns: True if the values are ok, otherwise false

            ok = true;
            for i = 1:16
                par = obj.parameters{i};
                if length(params.(par.name)) ~= nComponents
                    ok = false;
                    disp(['Error: ' par.name ' has to be a vector of length ' num2str(nComponents)]);
                end
            end

            for i = 17:length(obj.parameters)
                par = obj.parameters{i};
                if length(params.(par.name)) ~= 1
                    ok = false;
                    disp(['Error: ' par.name ' has to be a scalar']);
                end
            end
        end
    end
    
    methods (Access = 'private')

        function addParameterGroup(obj, namePrefix, offset, perComponent)
            %ADDPARAMETERGROUP Adds a group of parameters with the given prefix
            %
            % Parameters:
            %   - namePrefix: Name of the group
            %   - offset: Offset in the parameters property
            %   - perComponent: Set to true if each component has a parameter
            %
            % Example: addParameterGroup('test', 2, true) will populate
            %   parameters(3:7) with parameters named test, test_T, test_TT, and test_TTT.

            for i = 1:4
                name = namePrefix;
                if i > 1
                    name = [name '_' repmat('T', 1, i-1)];
                end
                obj.parameters{offset + i} = BindingModel.createParameter(['EXTSMA_' name], name, perComponent);
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
