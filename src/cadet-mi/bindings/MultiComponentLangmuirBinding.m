
classdef MultiComponentLangmuirBinding < BindingModel
    %MULTICOMPONENTLANGMUIRBINDING Represents the Multi Component Langmuir isotherm by providing its parameters and helper functions
    % Provides the parameters of the Multi Component Langmuir isotherm along with the methods to validate and convert them to
    % the types specified in the CADET file format.
    %
    % Copyright: © 2008-2015 Eric von Lieres, Joel Andersson, Andreas Püttmann, Sebastian Schnittert, Samuel Leweke
    %            See the license note at the end of the file.
    
    methods
        
        function obj = MultiComponentLangmuirBinding()
            %MULTICOMPONENTLANGMUIRADSORPTION Constructs a MultiComponentLangmuir object
            obj.name = 'MULTI_COMPONENT_LANGMUIR';
            obj.parameters = cell(3,1);

            obj.parameters{1} = BindingModel.createParameter('MCL_KA', 'KA', true);
            obj.parameters{2} = BindingModel.createParameter('MCL_KD', 'KD', true);
            obj.parameters{3} = BindingModel.createParameter('MCL_QMAX', 'QMAX', true);
        end
                
        function ok = checkParameters(obj, params, nComponents)
            %CHECKPARAMETERS Checks the values of the given parameters for plausibility
            %
            % Parameters: 
            %   - params: Struct with parameters as fields
            %
            % Returns: True if the values are ok, otherwise false

            ok = BindingModel.checkPositiveVector(params.MCL_QMAX, nComponents, 'MCL_QMAX');
            ok = BindingModel.checkNonnegativeVector(params.MCL_KA, nComponents, 'MCL_KA') && ok;
            ok = BindingModel.checkNonnegativeVector(params.MCL_KD, nComponents, 'MCL_KD') && ok;
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
