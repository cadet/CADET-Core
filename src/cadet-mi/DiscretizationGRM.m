
classdef DiscretizationGRM < handle
    %DISCRETIZATIONGRM Represents discretization specific parameters of a GRM model
    % Holds all discretization specific parameters (f.e., number of cells)
    % for a CADET simulation of a GRM discretization.
    %
    % Copyright: © 2008-2015 Eric von Lieres, Joel Andersson, Andreas Püttmann, Sebastian Schnittert, Samuel Leweke
    %            See the license note at the end of the file.
    
    properties
        
        input;                  % The root element of the CADET file format specification's input section.
                                % Advanced users may work with this structure directly. The normal workflow
                                % is to use the properties of this class which are synchronized with the
                                % struct.
                                
        nCellsColumn;           % Number of cells in the column
        nCellsParticle;         % Number of cells in the particle

        particleDiscretizationType; % Type of particle discretization (f.e., equivolume)
        particleCellPosition;       % Positions of the cells in the particle
        
        reconstructionType;     % Type of reconstruction
        wenoBoundaryHandling;   % How WENO handles left and right boundary
        wenoEpsilon;            % Epsilon in the WENO method
        wenoOrder;              % Order of the WENO method
        
    end
        
    properties (Access = 'protected')
        trackingInfo;       % Array of structures with infos on how to match properties to CADET file format specs
    end

    methods
        
        function obj = DiscretizationGRM()
            % DiscretizationGRM Constructs a DiscretizationGRM object and inserts as much default values as possible

            obj.createTrackingInfo();

            obj.nCellsColumn = 32;
            obj.nCellsParticle = 4;
            obj.particleCellPosition = [0.0, 0.5, 1.0];
            obj.particleDiscretizationType = 'EQUIDISTANT_PAR';
            obj.reconstructionType = 'WENO';
            obj.wenoBoundaryHandling = 0;
            obj.wenoEpsilon = 1e-12;
            obj.wenoOrder = 3;
        end

        function ok = checkValues(obj)
            % checkValues Checks the values of the parameters for plausibility
            %
            % Returns: True if everything is fine, otherwise false
            
            ok = true;

            ok = Helpers.checkPositiveScalar(obj.nCellsColumn, 'nCellsColumn') && ok;
            
            if obj.nCellsColumn <= 1
                ok = false;
                disp('Error: Number of cells in the column (nCellsColumn) has to be at least 2');
            end
            
            ok = Helpers.checkPositiveScalar(obj.nCellsParticle, 'nCellsParticle') && ok;

            ok = Helpers.checkPositiveScalar(obj.wenoEpsilon, 'wenoEpsilon') && ok;
            ok = Helpers.checkPositiveScalar(obj.wenoOrder, 'wenoOrder') && ok;
            
            if (obj.wenoOrder < 1) || (obj.wenoOrder > 3)
                ok = false;
                disp('Error: wenoOrder has to be between 1 and 3 (inclusive)');
            end

            if ~isscalar(obj.wenoBoundaryHandling)
                ok = false;
                disp('Error: wenoBoundaryHandling has to be scalar');
            end
            
            if (obj.wenoBoundaryHandling < 0) || (obj.wenoBoundaryHandling > 3)
                ok = false;
                disp('Error: wenoBoundaryHandling has to be between 0 and 3 (inclusive)');
            end
            
            if ~strcmp(obj.reconstructionType, 'WENO')
                ok = false;
                disp('Error: reconstructionType has to be WENO');
            end
            
            if ~strcmp(obj.particleDiscretizationType, 'EQUIDISTANT_PAR') && ...
                    ~strcmp(obj.particleDiscretizationType, 'EQUIVOLUME_PAR') && ...
                    ~strcmp(obj.particleDiscretizationType, 'USER_DEFINED_PAR')
                ok = false;
                disp('Error: particleDiscretizationType has to be one of EQUIDISTANT_PAR, EQUIVOLUME_PAR, and USER_DEFINED_PAR');
            end

            if strcmp(obj.particleDiscretizationType, 'USER_DEFINED_PAR')
                if length(obj.particleCellPosition) ~= obj.nCellsParticle + 1
                    ok = false;
                    disp('Error: Length of particleCellPosition has to be nCellsParticle + 1');
                end

                if any(obj.particleCellPosition(2:end) < obj.particleCellPosition(1:end-1))
                    ok = false;
                    disp('Error: particleCellPosition has to be monotonically increasing');
                end

                if (obj.particleCellPosition(1) ~= 0.0) || (obj.particleCellPosition(end) ~= 1.0)
                    ok = false;
                    disp('Error: particleCellPosition has to begin with 0.0 and end with 1.0');
                end
            end

        end

        function synchToStruct(obj)
            % synchToStruct Synchronizes the parameters to the struct ("input" property)

            obj.input = [];

            % Handle tracked variables
            for i = 1:length(obj.trackingInfo)
                ti = obj.trackingInfo{i};
                
                if ~ti.optional || (ti.optional && ~isempty(obj.(ti.propName)))
                    obj.input = Helpers.setValueToNestedStruct(obj.input, splitstring(ti.structName, '.'), ti.convToStruct(obj.(ti.propName)));
                end
            end

        end

        function synchFromStruct(obj, input)
            % synchFromStruct Synchronizes the parameters by taking them from the given struct
            %
            % Parameters:
            %   - input: Optional. Structure with all necessary fields for
            %       the Model (see CADET file format specs). If left out,
            %       the property "input" is used to synchronize.
            
            if (nargin > 1) && ~isempty(input)
                obj.input = input;
            end
            
            % Handle tracked variables
            for i = 1:length(obj.trackingInfo)
                ti = obj.trackingInfo{i};
                
                [exists, value] = Helpers.getValueFromNestedStruct(obj.input, splitstring(ti.structName, '.'));
                if ~exists && ~ti.optional
                    warning('CADET:trackedFieldNotInStruct', 'The tracked field %s, which corresponds to the property %s, is not in the given struct', ti.structName, ti.propName);
                end
                obj.(ti.propName) = ti.convToProp(value);
            end
            
        end
        
    end
    
    methods (Access = 'protected')

        function createTrackingInfo(obj)
            % createTrackingInfo Populates the tracking list
            % Tracked properties are easily synchronized with their struct
            % pendants (see CADET's file format specs). All you need to
            % setup a tracking is to specify the name of the property and
            % the corresponding field (with path) of the struct according
            % to the file format specs.
            %
            % If you add more parameters to the model which need to be put
            % in the file format, best practice is to add tracking info for
            % them.
            %
            % Conversion from the datatype of the property to the datatype
            % of the struct and back is done by means of function handles.
            
            obj.addTracking('nCellsColumn',     'discretization.NCOL', @(x) int32(x), @(x) double(x));
            obj.addTracking('nCellsParticle',   'discretization.NPAR', @(x) int32(x), @(x) double(x));

            obj.addTracking('particleDiscretizationType',   'discretization.PAR_DISC_TYPE');
            obj.addTracking('particleCellPosition',         'discretization.PAR_DISC_VECTOR');
            obj.addTracking('reconstructionType',           'discretization.RECONSTRUCTION');
            obj.addTracking('wenoBoundaryHandling',         'discretization.weno.BOUNDARY_MODEL', @(x) int32(x), @(x) double(x));
            obj.addTracking('wenoEpsilon',                  'discretization.weno.WENO_EPS');
            obj.addTracking('wenoOrder',                    'discretization.weno.WENO_ORDER', @(x) int32(x), @(x) double(x));
        end

        function addTracking(obj, propName, structName, convToStruct, convToProp, optional)
            % addTracking Adds a tracking info to the list
            %
            % Parameters:
            %   - propName: Name of the property
            %   - structName: Path to the field in the nested structs
            %   - convToStruct: Optional. Function handle which converts
            %       the property value to a struct value
            %   - convToProp: Optional. Function handle which converts the
            %       struct value to a property value
            %   - optional: Optional. True if the item is optional and
            %       should only be created if it contains some data
            
            info.propName = propName;
            info.structName = structName;
            
            if nargin <= 3 || isempty(convToStruct)
                info.convToStruct = @(x) x;
            else
                info.convToStruct = convToStruct;
            end
            
            if nargin <= 4 || isempty(convToProp)
                info.convToProp = @(x) x;
            else
                info.convToProp = convToProp;
            end
            
            if nargin <= 5 || isempty(optional)
                info.optional = false;
            else
                info.optional = optional;
            end
            
            obj.trackingInfo{end+1} = info;
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
