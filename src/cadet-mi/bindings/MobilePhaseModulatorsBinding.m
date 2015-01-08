
classdef MobilePhaseModulatorsBinding < BindingModel
    %MOBILEPHASEMODULATORSBINDING Represents the Mobile Phase Modulators isotherm by providing its parameters and helper functions
    % Provides the parameters of the Mobile Phase Modulators isotherm along with the methods to validate and convert them to
    % the types specified in the CADET file format.
    %
    % Copyright: © 2008-2015 Eric von Lieres, Joel Andersson, Andreas Püttmann, Sebastian Schnittert, Samuel Leweke
    %            See the license note at the end of the file.

    methods
        
        function obj = MobilePhaseModulatorsBinding()
            %MOBILEPHASEMODULATORSBINDING Constructs a MobilePhaseModulatorsBinding object
            obj.name = 'MOBILE_PHASE_MODULATORS';
            obj.parameters = cell(5,1);

            obj.parameters{1} = BindingModel.createParameter('MPM_KA', 'KA', true);
            obj.parameters{2} = BindingModel.createParameter('MPM_KD', 'KD', true);
            obj.parameters{3} = BindingModel.createParameter('MPM_QMAX', 'QMAX', true);
            obj.parameters{4} = BindingModel.createParameter('MPM_BETA', 'BETA', true);
            obj.parameters{5} = BindingModel.createParameter('MPM_GAMMA', 'GAMMA', true);
        end
                
        function ok = checkParameters(obj, params, nComponents)
            %CHECKPARAMETERS Checks the values of the given parameters for plausibility
            %
            % Parameters: 
            %   - params: Struct with parameters as fields
            %
            % Returns: True if the values are ok, otherwise false

            ok = BindingModel.checkNonnegativeVector(params.MPM_KA, nComponents, 'MPM_KA');
            ok = BindingModel.checkNonnegativeVector(params.MPM_KD, nComponents, 'MPM_KD') && ok;
            ok = BindingModel.checkNonnegativeVector(params.MPM_QMAX, nComponents, 'MPM_QMAX') && ok;
            ok = BindingModel.checkNonnegativeVector(params.MPM_BETA, nComponents, 'MPM_BETA') && ok;
            ok = BindingModel.checkNonnegativeVector(params.MPM_GAMMA, nComponents, 'MPM_GAMMA') && ok;
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
