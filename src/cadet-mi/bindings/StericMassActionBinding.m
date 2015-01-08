
classdef StericMassActionBinding < BindingModel
    %STERICMASSACTIONBINDING Represents the Steric Mass Action isotherm by providing its parameters and helper functions
    % Provides the parameters of the Steric Mass Action isotherm along with the methods to validate and convert them to
    % the types specified in the CADET file format.
    %
    % Copyright: © 2008-2015 Eric von Lieres, Joel Andersson, Andreas Püttmann, Sebastian Schnittert, Samuel Leweke
    %            See the license note at the end of the file.
    
    methods
        
        function obj = StericMassActionBinding()
            %STERICMASSACTIONADSORPTION Constructs a StericMassAction object
            obj.name = 'STERIC_MASS_ACTION';
            obj.parameters = cell(5,1);

            obj.parameters{1} = BindingModel.createParameter('SMA_KA', 'KA', true);
            obj.parameters{2} = BindingModel.createParameter('SMA_KD', 'KD', true);
            obj.parameters{3} = BindingModel.createParameter('SMA_NU', 'NU', true);
            obj.parameters{4} = BindingModel.createParameter('SMA_SIGMA', 'SIGMA', true);
            obj.parameters{5} = BindingModel.createParameter('SMA_LAMBDA', 'LAMBDA', false);
        end
                
        function ok = checkParameters(obj, params, nComponents)
            %CHECKPARAMETERS Checks the values of the given parameters for plausibility
            %
            % Parameters: 
            %   - params: Struct with parameters as fields
            %
            % Returns: True if the values are ok, otherwise false

            ok = BindingModel.checkNonnegativeVector(params.SMA_KA, nComponents, 'SMA_KA');
            ok = BindingModel.checkNonnegativeVector(params.SMA_KD, nComponents, 'SMA_KD') && ok;
            ok = BindingModel.checkNonnegativeVector(params.SMA_NU, nComponents, 'SMA_NU') && ok;
            ok = BindingModel.checkNonnegativeVector(params.SMA_SIGMA, nComponents, 'SMA_SIGMA') && ok;
            ok = BindingModel.checkNonnegativeVector(params.SMA_LAMBDA, 1, 'SMA_KA') && ok;
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
