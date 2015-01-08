
classdef MultiComponentBiLangmuirBinding < BindingModel
    %MULTICOMPONENTBILANGMUIRBINDING Represents the Multi Component Bi-Langmuir isotherm by providing its parameters and helper functions
    % Provides the parameters of the Multi Component Bi-Langmuir isotherm along with the methods to validate and convert them to
    % the types specified in the CADET file format.
    %
    % Copyright: © 2008-2015 Eric von Lieres, Joel Andersson, Andreas Püttmann, Sebastian Schnittert, Samuel Leweke
    %            See the license note at the end of the file.
    
    methods
        
        function obj = MultiComponentBiLangmuirBinding()
            %MULTICOMPONENTBILANGMUIRBINDING Constructs a MultiComponentBiLangmuir object
            obj.name = 'MULTI_COMPONENT_BILANGMUIR';
            obj.parameters = cell(6,1);

            obj.parameters{1} = BindingModel.createParameter('MCBL_KA1', 'KA1', true);
            obj.parameters{2} = BindingModel.createParameter('MCBL_KD1', 'KD1', true);
            obj.parameters{3} = BindingModel.createParameter('MCBL_KA2', 'KA2', true);
            obj.parameters{4} = BindingModel.createParameter('MCBL_KD2', 'KD2', true);
            obj.parameters{5} = BindingModel.createParameter('MCBL_QMAX1', 'QMAX1', true);
            obj.parameters{6} = BindingModel.createParameter('MCBL_QMAX2', 'QMAX2', true);
        end
                
        function ok = checkParameters(obj, params, nComponents)
            %CHECKPARAMETERS Checks the values of the given parameters for plausibility
            %
            % Parameters: 
            %   - params: Struct with parameters as fields
            %
            % Returns: True if the values are ok, otherwise false

            if (mod(nComponents, 2) ~= 0)
                disp(['For emulated second bound state the number of components must be divisible by 2 (got ' num2str(nComponents) ')'])
                ok = false;
            else
                ok = true;
            end

            ok = BindingModel.checkPositiveVector(params.MCBL_QMAX1, nComponents / 2, 'MCBL_QMAX1') && ok;
            ok = BindingModel.checkPositiveVector(params.MCBL_QMAX2, nComponents / 2, 'MCBL_QMAX2') && ok;
            ok = BindingModel.checkNonnegativeVector(params.MCBL_KA1, nComponents / 2, 'MCBL_KA1') && ok;
            ok = BindingModel.checkNonnegativeVector(params.MCBL_KD1, nComponents / 2, 'MCBL_KD1') && ok;
            ok = BindingModel.checkNonnegativeVector(params.MCBL_KA2, nComponents / 2, 'MCBL_KA2') && ok;
            ok = BindingModel.checkNonnegativeVector(params.MCBL_KD2, nComponents / 2, 'MCBL_KD2') && ok;
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
