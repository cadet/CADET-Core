
classdef LinearBinding < BindingModel
    %LINEARBINDING Represents the linear isotherm by providing its parameters and helper functions
    % Provides the parameters of the linear isotherm along with the methods to validate and convert them to
    % the types specified in the CADET file format.
    %
    % Copyright: © 2008-2015 Eric von Lieres, Joel Andersson, Andreas Püttmann, Sebastian Schnittert, Samuel Leweke
    %            See the license note at the end of the file.
    
    methods
        
        function obj = LinearBinding()
            %LINEARADSORPTION Constructs a LinearBinding object
            obj.name = 'LINEAR';
            obj.parameters = cell(2,1);

            obj.parameters{1} = BindingModel.createParameter('LIN_KA', 'KA', true);
            obj.parameters{2} = BindingModel.createParameter('LIN_KD', 'KD', true);
        end

        function ok = checkParameters(obj, params, nComponents)
            %CHECKPARAMETERS Checks the values of the given parameters for plausibility
            %
            % Parameters: 
            %   - params: Struct with parameters as fields
            %
            % Returns: True if the values are ok, otherwise false

            ok = BindingModel.checkNonnegativeVector(params.LIN_KA, nComponents, 'LIN_KA');
            ok = BindingModel.checkNonnegativeVector(params.LIN_KD, nComponents, 'LIN_KD') && ok;
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
