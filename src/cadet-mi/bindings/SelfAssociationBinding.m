
classdef SelfAssociationBinding < BindingModel
    %SELFASSOCIATIONBINDING Represents the Self Associtation isotherm by providing its parameters and helper functions
    % Provides the parameters of the Self Associtation isotherm along with the methods to validate and convert them to
    % the types specified in the CADET file format.
    %
    % Copyright: © 2008-2015 Eric von Lieres, Joel Andersson, Andreas Püttmann, Sebastian Schnittert, Samuel Leweke
    %            See the license note at the end of the file.
    
    methods
        
        function obj = SelfAssociationBinding()
            %SELFASSOCIATIONADSORPTION Constructs a SelfAssociation object
            obj.name = 'SELF_ASSOCIATION';
            obj.parameters = cell(6,1);

            obj.parameters{1} = BindingModel.createParameter('SAI_KA1', 'KA1', true);
            obj.parameters{2} = BindingModel.createParameter('SAI_KA2', 'KA2', true);
            obj.parameters{3} = BindingModel.createParameter('SAI_KD', 'KD', true);
            obj.parameters{4} = BindingModel.createParameter('SAI_NU', 'NU', true);
            obj.parameters{5} = BindingModel.createParameter('SAI_SIGMA', 'SIGMA', true);
            obj.parameters{6} = BindingModel.createParameter('SAI_LAMBDA', 'LAMBDA', false);
        end
                
        function ok = checkParameters(obj, params, nComponents)
            %CHECKPARAMETERS Checks the values of the given parameters for plausibility
            %
            % Parameters: 
            %   - params: Struct with parameters as fields
            %
            % Returns: True if the values are ok, otherwise false

            ok = BindingModel.checkNonnegativeVector(params.SAI_KA1, nComponents, 'SAI_KA1');
            ok = BindingModel.checkNonnegativeVector(params.SAI_KA2, nComponents, 'SAI_KA2') && ok;
            ok = BindingModel.checkNonnegativeVector(params.SAI_KD, nComponents, 'SAI_KD') && ok;
            ok = BindingModel.checkNonnegativeVector(params.SAI_NU, nComponents, 'SAI_NU') && ok;
            ok = BindingModel.checkNonnegativeVector(params.SAI_SIGMA, nComponents, 'SAI_SIGMA') && ok;
            ok = BindingModel.checkNonnegativeVector(params.SAI_LAMBDA, 1, 'SAI_LAMBDA') && ok;
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
