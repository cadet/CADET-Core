
classdef ModelGRM < handle
    %MODELGRM Represents a model of a single column with all model specific parameters for a CADET simulation 
    % Holds all model specific parameters (f.e., interstitial velocity, dispersion, etc.) which are necessary
    % for a CADET simulation. This class is mainly concerned with the transport model. Binding is handled by an
    % binding model class.
    %
    % Adding new parameters to the Model class:
    %   1. Add a property to hold the value(s)
    %   2. Add a tracking info line to the createTrackingInfo() function
    %   3. Add validation code to checkValues() function
    %
    % Copyright: © 2008-2015 Eric von Lieres, Joel Andersson, Andreas Püttmann, Sebastian Schnittert, Samuel Leweke
    %            See the license note at the end of the file.
    
    properties
        
        input;                  % The root element of the CADET file format specification's input section.
                                % Advanced users may work with this structure directly. The normal workflow
                                % is to use the properties of this class which are synchronized with the
                                % struct.

        nComponents;            % Number of chemical components in the column

        % Initial conditions
        
        initialMobileConcentration;     % Initial concentration of the mobile phase for each component
        initialBeadConcentration;       % Initial concentration of the liquid phase in the beads for each component
        initialSolidConcentration;      % Initial concentration of the solid phase for each component
        initialState;                   % Initial state vector
        initialSensitivities;           % Initial state variables of all sensitive parameters

        % Binding model
        
        bindingParameters;     % Struct with parameters of the binding model

        kineticBindingModel;    % Boolean indicating if the binding happens "infinitely" fast (false) or
                                % on a comparable time scale than the transport (true)

        % Transport
        
        dispersionColumn;           % Dispersion coefficient of the mobile phase transport inside the column
        interstitialVelocity;       % Interstitial velocity of the mobile phase transport inside the column

        porosityColumn;             % Porosity of the column
        porosityParticle;           % Porosity of the particle

        filmDiffusion;              % Film diffusion coefficient
        diffusionParticle;          % Diffusion coefficient of the mobile phase components inside the particle
        diffusionParticleSurface;   % Diffusion coefficient of the mobile phase components on the surface of the particle

        % Geometry
        
        columnLength;               % Length of the column
        particleRadius;             % Radius of the particles

        % Inlet
        
        nInletSections;         % Number of sections in the inlet specification
        sectionTimes;           % Vector with start and end times of sections (length: nInletSections + 1)
        sectionContinuity;      % Vector with booleans indicating continuity of a section transition (length: nInletSections - 1)
        sectionConstant;        % Vector with constant coefficients of the piecewise polynomials defining the inlet profile
        sectionLinear;          % Vector with linear coefficients of the piecewise polynomials defining the inlet profile
        sectionQuadratic;       % Vector with quadratic coefficients of the piecewise polynomials defining the inlet profile
        sectionCubic;           % Vector with cubic coefficients of the piecewise polynomials defining the inlet profile
        
        % External data profile
        
        extFuncVelocity;        % External function velocity (compare EXT_VELOCITY in file format specs)
        extFuncProfile;         % External function profile data (compare EXT_PROFILE in file format specs)
        extFuncDelta;           % External function spacing between data points (compare EXT_PROF_DELTA in file format specs)
    end
    
    properties (Dependent = true)
        bindingModel;               % Object of the used binding model class
    end
    
    properties(Constant)
        modelType = 'GENERAL_RATE_MODEL'; % Type of the model according to CADET file format specs
    end
    
    properties (Access = 'protected')
        trackingInfo;       % Array of structures with infos on how to match properties to CADET file format specs
        
        bindingHandle;     % Storage for the actual object of the used BindingModel class
    end

    methods
        
        function obj = ModelGRM()
            % ModelGRM Constructs a ModelGRM object and inserts as much default values as possible

            obj.createTrackingInfo();

            obj.bindingHandle = [];
            obj.bindingParameters = [];
            obj.nComponents = 0;
            obj.extFuncDelta = [];
            obj.extFuncProfile = [];
            obj.extFuncVelocity = [];
        end

        function ok = checkValues(obj)
            % checkValues Checks the values of the parameters for plausibility
            %
            % Returns: True if everything is fine, otherwise false
            
            ok = true;

            % General
            ok = Helpers.checkPositiveScalar(obj.nComponents, 'nComponents') && ok;
            
            % Binding model
            if isempty(obj.bindingHandle)
                ok = false;
                disp('Error: Binding model is unset');
            end

            params = obj.bindingHandle.parameters;
            for i = 1:length(params)
                curPar = params{i};
                if ~isfield(obj.bindingParameters, curPar.name)
                    ok = false;
                    disp(['Error: bindingParameters has to contain the field "' curPar.name '" when binding model ' obj.bindingHandle.name ' is chosen']);
                end
            end

            if ~obj.bindingHandle.checkParameters(obj.bindingParameters, obj.nComponents)
                ok = false;
                disp(['Error: bindingParameters are not valid for binding model ' obj.bindingHandle.name]);
            end

            ok = Helpers.checkScalar(obj.kineticBindingModel, 'kineticBindingModel') && ok;

            % Initial conditions
            if (isempty(obj.initialState))
                ok = Helpers.checkNonnegativeVector(obj.initialMobileConcentration, obj.nComponents, 'initialMobileConcentration') && ok;
                ok = Helpers.checkNonnegativeVector(obj.initialSolidConcentration, obj.nComponents, 'initialSolidConcentration') && ok;
                if (~isempty(obj.initialBeadConcentration))
                    ok = Helpers.checkNonnegativeVector(obj.initialBeadConcentration, obj.nComponents, 'initialBeadConcentration') && ok;
                end
            end
            
            % Transport
            if any(obj.dispersionColumn) < 0
                ok = false;
                disp(['Error: Each entry of dispersionColumn has to be nonnegative, i.e., greater than or equal 0']);
            end
            if (length(obj.dispersionColumn) ~= length(obj.sectionTimes) - 1) && ~isscalar(obj.dispersionColumn)
                ok = false;
                disp(['Error: dispersionColumn has to be scalar or a vector of length ' num2str(length(obj.sectionTimes) - 1)]);
            end
            
            if any(obj.interstitialVelocity < 0)
                ok = false;
                disp(['Error: Each entry of interstitialVelocity has to be nonnegative, i.e., greater than or equal 0']);
            end
            if (length(obj.interstitialVelocity) ~= length(obj.sectionTimes) - 1) && ~isscalar(obj.interstitialVelocity)
                ok = false;
                disp(['Error: interstitialVelocity has to be scalar or a vector of length ' num2str(length(obj.sectionTimes) - 1)]);
            end

            ok = Helpers.checkNonnegativeScalar(obj.porosityColumn, 'porosityColumn') && ok;
            ok = Helpers.checkNonnegativeScalar(obj.porosityParticle, 'porosityParticle') && ok;
            
            ok = Helpers.checkNonnegativeVector(obj.filmDiffusion, obj.nComponents .* [1, length(obj.sectionTimes)-1], 'filmDiffusion') && ok;
            ok = Helpers.checkNonnegativeVector(obj.diffusionParticle, obj.nComponents .* [1, length(obj.sectionTimes)-1], 'diffusionParticle') && ok;
            ok = Helpers.checkNonnegativeVector(obj.diffusionParticleSurface, obj.nComponents .* [1, length(obj.sectionTimes)-1], 'diffusionParticleSurface') && ok;

            % Geometry
            ok = Helpers.checkPositiveScalar(obj.columnLength, 'columnLength') && ok;
            ok = Helpers.checkPositiveScalar(obj.particleRadius, 'particleRadius') && ok;

            % Check inlet
            ok = Helpers.checkPositiveScalar(obj.nInletSections, 'nInletSections') && ok;

            if obj.nInletSections + 1 ~= length(obj.sectionTimes)
                ok = false;
                disp('Error: Length of sectionTimes has to be nInletSections + 1');
            end

            if any(obj.sectionTimes(2:end) < obj.sectionTimes(1:end-1))
                ok = false;
                disp('Error: sectionTimes has to be monotonically increasing');
            end

            if any(obj.sectionTimes() < 0)
                ok = false;
                disp('Error: sectionTimes has to be nonnegative');
            end

            if obj.nInletSections - 1 ~= length(obj.sectionContinuity)
                ok = false;
                disp('Error: Length of sectionContinuity has to be nInletSections - 1');
            end

            if ~all([obj.nComponents, obj.nInletSections] == size(obj.sectionConstant))
                ok = false;
                disp('Error: Size of sectionConstant has to be nComponents x nInletSections');
            end

            if ~all([obj.nComponents, obj.nInletSections] == size(obj.sectionLinear))
                ok = false;
                disp('Error: Size of sectionLinear has to be nComponents x nInletSections');
            end

            if ~all([obj.nComponents, obj.nInletSections] == size(obj.sectionQuadratic))
                ok = false;
                disp('Error: Size of sectionQuadratic has to be nComponents x nInletSections');
            end

            if ~all([obj.nComponents, obj.nInletSections] == size(obj.sectionCubic))
                ok = false;
                disp('Error: Size of sectionCubic has to be nComponents x nInletSections');
            end

            % External function
            if ~isempty(obj.extFuncProfile)

                if ~isscalar(obj.extFuncVelocity)
                    ok = false;
                    disp('Error: extFuncVelocity has to be a scalar');
                end

                if length(obj.extFuncProfile) <= 1
                    ok = false;
                    disp('Error: extFuncDelta has to be a vector of length greater than 1');
                end

                if length(obj.extFuncProfile) ~= length(obj.extFuncDelta)
                    ok = false;
                    disp('Error: extFuncProfile and extFuncDelta have to have the same length');
                end

                if obj.extFuncDelta(1) ~= 0.0
                    ok = false;
                    disp('Error: First entry of extFuncDelta has to be 0.0');
                end

                if any(obj.extFuncDelta(2:end) <= 0.0)
                    ok = false;
                    disp('Error: All entries but the first of extFuncDelta have to be positive');
                end
            end
        end
        
        function val = get.bindingModel(obj)
            val = obj.bindingHandle;
        end
        
        function set.bindingModel(obj, newBinding)
            obj.bindingParameters = [];
            if ~isempty(newBinding) && isobject(newBinding)
                % Take parameters from binding model
                obj.bindingHandle = newBinding;

                % Add fields in struct for parameters
                params = obj.bindingHandle.parameters;
                for i = 1:length(params)
                    par = params{i};

                    if isfield(par, 'perComponent') && par.perComponent
                        value = zeros(obj.nComponents, 1);
                    else
                        value = 0;
                    end

                    obj.bindingParameters.(par.name) = value;
                end
            else
                obj.bindingHandle = [];
            end
        end

        function synchToStruct(obj)
            % synchToStruct Synchronizes the parameters to the struct ("input" property)

            obj.input = [];
            obj.input.CHROMATOGRAPHY_TYPE = obj.modelType;
            
            % Handle tracked variables
            for i = 1:length(obj.trackingInfo)
                ti = obj.trackingInfo{i};
                
                if ~ti.optional || (ti.optional && ~isempty(obj.(ti.propName)))
                    obj.input = Helpers.setValueToNestedStruct(obj.input, splitstring(ti.structName, '.'), ti.convToStruct(obj.(ti.propName)));
                end
            end

            obj.input.model.ADSORPTION_TYPE = obj.bindingHandle.name;
            
            % Copy binding parameters
            params = obj.bindingHandle.parameters;
            for i = 1:length(params)
                par = params{i};
                obj.input.model.adsorption.(par.name) = obj.bindingParameters.(par.name);
            end

            % Initial state
            if (~isempty(obj.initialState))
                obj.input.model.INIT_STATE = obj.initialState;
            end
            if (~isempty(obj.initialSensitivities))
                obj.input.model.INIT_SENS = obj.initialSensitivities;
            end
            if (~isempty(obj.initialBeadConcentration))
                obj.input.model.INIT_CP = obj.initialBeadConcentration;
            end

            % Inlet
            for i = 1:obj.nInletSections
                sec_str = sprintf('sec_%03d', i-1);
                if ~isfield(sec_str, obj.input.model.inlet)
                    obj.input.model.inlet.(sec_str) = [];
                end
                
                obj.input.model.inlet.(sec_str).CONST_COEFF = obj.sectionConstant(:, i);
                obj.input.model.inlet.(sec_str).LIN_COEFF = obj.sectionLinear(:, i);
                obj.input.model.inlet.(sec_str).QUAD_COEFF = obj.sectionQuadratic(:, i);
                obj.input.model.inlet.(sec_str).CUBE_COEFF = obj.sectionCubic(:, i);
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

            assert(strcmp(obj.input.CHROMATOGRAPHY_TYPE, obj.modelType), ...
                'CADET:invalidModelType', 'The ModelGRM class does not accept models of type %s', obj.input.CHROMATOGRAPHY_TYPE);
            
            % Handle tracked variables
            for i = 1:length(obj.trackingInfo)
                ti = obj.trackingInfo{i};
                
                [exists, value] = Helpers.getValueFromNestedStruct(obj.input, splitstring(ti.structName, '.'));
                if ~exists && ~ti.optional
                    warning('CADET:trackedFieldNotInStruct', ...
                        'The tracked field %s, which corresponds to the property %s, is not in the given struct', ti.structName, ti.propName);
                end
                obj.(ti.propName) = ti.convToProp(value);
            end

            % Try to find the binding model by looking through all classes in the "binding" folder
            obj.bindingHandle = [];
            pathToIsotherms = fullfile(fileparts(mfilename('fullpath')), 'bindings');
            classes = dir(fullfile(pathToIsotherms, '*.m'));
            for i = 1:length(classes)
                % Create instance of class
                [~, name, ~] = fileparts(classes(i).name);
                instance = eval([name '();']);

                if strcmpi(instance.name, obj.input.model.ADSORPTION_TYPE)
                    obj.bindingHandle = instance;
                    break;
                end
            end

            if isempty(obj.bindingHandle)
                warning('CADET:bindingNotFound', 'Could not find binding model class named %s', obj.input.model.ADSORPTION_TYPE);
            end

            % Copy binding model parameters
            params = obj.bindingHandle.parameters;
            for i = 1:length(params)
                par = params{i};

                % Determine default value
                if isfield(par, 'perComponent') && par.perComponent
                    value = zeros(obj.nComponents, 1);
                else
                    value = 0;
                end

                if isfield(obj.input.model.adsorption, par.name)
                    % Set given value
                    obj.bindingParameters.(par.name) = obj.input.model.adsorption.(par.name);
                else
                    % Set default
                    obj.bindingParameters.(par.name) = value;
                end
            end            
            
            % Inlet
            obj.sectionConstant = zeros(obj.nComponents, obj.nInletSections);
            obj.sectionLinear = zeros(obj.nComponents, obj.nInletSections);
            obj.sectionQuadratic = zeros(obj.nComponents, obj.nInletSections);
            obj.sectionCubic = zeros(obj.nComponents, obj.nInletSections);

            for i = 1:obj.nInletSections
                sec_str = sprintf('sec_%03d', i-1);
                obj.sectionConstant(:, i) = obj.input.model.inlet.(sec_str).CONST_COEFF;
                obj.sectionLinear(:, i) = obj.input.model.inlet.(sec_str).LIN_COEFF;
                obj.sectionQuadratic(:, i) = obj.input.model.inlet.(sec_str).QUAD_COEFF;
                obj.sectionCubic(:, i) = obj.input.model.inlet.(sec_str).CUBE_COEFF;
            end

            % Initial state
            obj.initialState = [];
            if (isfield(obj.input.model, 'INIT_STATE'))
                obj.initialState = obj.input.model.INIT_STATE;
            end
            obj.initialSensitivities = [];
            if (isfield(obj.input.model, 'INIT_SENS'))
                obj.initialSensitivities = obj.input.model.INIT_SENS;
            end
            obj.initialBeadConcentration = [];
            if (isfield(obj.input.model, 'INIT_CP'))
                obj.initialBeadConcentration = obj.input.model.INIT_CP;
            end
        end
        
        function setInlet(obj, pp, componentIndex, continuityThreshold)
            % setInlet Takes a piecewise polynomial and converts it to CADET's inlet spline format.
            %
            % Note that calling this function resets the sectionContinuity
            % since it tries to detect continuity automatically.
            %
            % Parameters:
            %   - pp: Piecewise polynomial in Matlab format (result of
            %       spline(...) or interp1(..., 'pp') )
            %   - componentIndex: Index of the component (1-based)
            %   - continuityThreshold: Optional. Threshold for detection of
            %       section continuity (default: 1e-10). A negative value
            %       disables continuity detection.

            if (nargin <= 3) || isempty(continuityThreshold)
                continuityThreshold = 1e-10;
            end

            [breaks, coefs] = unmkpp(pp);

            obj.nInletSections = length(breaks) - 1;
            obj.sectionTimes = breaks(:);
            obj.sectionContinuity = false(obj.nInletSections - 1, 1);

            % Each row in coefs corresponds to a polynomial
            % The first entry in a row is the coefficient of the highest monomial

            % Resize arrays if necessary
            if ~all(size(obj.sectionConstant) == [obj.nComponents, obj.nInletSections])
                obj.sectionConstant = zeros(obj.nComponents, obj.nInletSections);
                obj.sectionLinear = zeros(obj.nComponents, obj.nInletSections);
                obj.sectionQuadratic = zeros(obj.nComponents, obj.nInletSections);
                obj.sectionCubic = zeros(obj.nComponents, obj.nInletSections);
            end

            % Populate arrays
            for i = 1:obj.nInletSections
                obj.sectionConstant(componentIndex,i) = coefs(i,end);

                if size(coefs, 2) > 1
                    obj.sectionLinear(componentIndex,i) = coefs(i,end-1);
                else
                    obj.sectionLinear(componentIndex,i) = 0;
                end

                if size(coefs, 2) > 2
                    obj.sectionQuadratic(componentIndex,i) = coefs(i,end-2);
                else
                    obj.sectionQuadratic(componentIndex,i) = 0;
                end

                if size(coefs,2) > 3
                    obj.sectionCubic(componentIndex,i) = coefs(i,end-3);
                else
                    obj.sectionCubic(componentIndex,i) = 0;
                end
            end
            
            if continuityThreshold >= 0
                obj.detectContinuity(continuityThreshold);
            end
        end
        
        function setInletsFromData(obj, timePoints, dataPoints, componentIndex, extrapolation)
            % setInletsFromData Takes data in form of time-concentration pairs for each component and assigns an interpolating piecewise polynomial to the inlet
            %
            % Note that calling this function resets the sectionContinuity
            % since it tries to detect continuity automatically.
            %
            % Parameters:
            %   - timePoints: Cell array with vector of time points
            %   - dataPoints: Cell array with vectors of concentrations at the time points
            %   - componentIndex: Optional. Array with indices of the components (1-based)
            %   - extrapolation: Optional. Mode of extrapolation. One of
            %       'zero', 'last', 'interp1'. Mode 'zero' sets signal to
            %       zero outside of its time domain, 'last' uses the last
            %       value, and 'interp1' uses Matlab's interp1 function.
            %       Default: 'zero'

            if nargin <= 3 || isempty(componentIndex)
                componentIndex = 1:length(timePoints);
            end

            if nargin <= 4 || isempty(extrapolation)
                extrapolation = 'zero';
            end
            
            assert(length(timePoints) == length(dataPoints), 'CADET:invalidInletProfile', ...
                'The number of time vectors differs from the number of data vectors');
            assert(length(timePoints) == length(componentIndex), 'CADET:invalidInletProfile', ...
                'The number of time vectors differs from the length of the component indices');
            
            % First pass: Determine unique time points
            globalTime = [];
            for i = 1:length(timePoints)
                tp = timePoints{i};
                dp = dataPoints{i};
                
                assert(length(tp) == length(dp), 'CADET:invalidInletProfile', ...
                    'The number of data points differs from the number of time points in inlet %d', i);
                
                globalTime = unique([globalTime; tp(:)]);
            end
            
            % Second pass: Interpolate to get all signals to the same resolution
            globalTime = sort(globalTime);
            for i = 1:length(timePoints)
                tp = timePoints{i};
                dp = dataPoints{i};
                
                newDp = interp1(tp, dp, globalTime, 'cubic');
                
                % Handle extrapolation
                if ~strcmpi(extrapolation, 'interp1')
                    
                    if strcmpi(extrapolation, 'zero')
                        
                        idx = (globalTime > tp(end)) | (globalTime < tp(1));
                        newDp(idx) = 0;
                        
                    elseif strcmpi(extrapolation, 'last')
                        
                        idx = (globalTime > tp(end));
                        newDp(idx) = dp(end);

                        idx = (globalTime < tp(end));
                        newDp(idx) = dp(1);
                        
                    else
                        error('CADET:invalidExtrapolationMethod', 'Extrapolation method %s is not supported', extrapolation);
                    end
                end
                
                pp = interp1(globalTime, newDp, 'cubic', 'pp');
                obj.setInlet(pp, componentIndex(i), -1);
            end
            
            obj.detectContinuity(1e-10);
         end
         
         function setInletFromData(obj, timePoints, dataPoints, componentIndex)
            % setInletFromData Takes data in form of time-concentration pairs and assigns an interpolating piecewise polynomial to the inlet
            % By default a piecewise cubic polynomial is used.
            %
            % Note that calling this function resets the sectionContinuity
            % since it tries to detect continuity automatically.
            %
            % Parameters:
            %   - timePoints: Vector with time points
            %   - dataPoints: Vector with concentrations at the time points
            %   - componentIndex: Index of the component (1-based)
            
            assert(length(timePoints) == length(dataPoints), 'CADET:invalidInletProfile', ...
                'The number of time points differs from the number of data points');
            assert(all(timePoints(2:end) > timePoints(1:end-1)), 'CADET:invalidInletProfile', ...
                'The time points have to increase monotonously');
            assert(all(timePoints(2:end) > timePoints(1:end-1)), 'CADET:invalidInletProfile', ...
                'The time points have to increase monotonously');

            pp = interp1(timePoints, dataPoints, 'cubic', 'pp');
            obj.setInlet(pp, componentIndex);
        end
        
        function [result, path] = hasParameter(obj, paramName)
            % hasParameter Checks whether the given paramter exists in the model
            %
            % Parameters:
            %   - paramName: Name of the parameter according to CADET's
            %       file format specs
            %
            % Returns: 
            %   - result: True if the parameter exists, otherwise false
            %   - path: Path to the field in the nested struct
                        
            idxModel = find(cellfun(@checkParamExistence, obj.trackingInfo));
            if ~isempty(idxModel)
                result = true;
                path = obj.trackingInfo{idxModel(1)}.structName;
                return;
            end
            
            if strcmpi(paramName, 'CONST_COEFF') || strcmpi(paramName, 'sectionconstant')
                result = true;
                path = 'model.inlet.sec_%03d.CONST_COEFF';
                return;
            end
                
            if strcmpi(paramName, 'LIN_COEFF') || strcmpi(paramName, 'sectionlinear')
                result = true;
                path = 'model.inlet.sec_%03d.LIN_COEFF';
                return;
            end
            
            if strcmpi(paramName, 'QUAD_COEFF') || strcmpi(paramName, 'sectionquadratic')
                result = true;
                path = 'model.inlet.sec_%03d.QUAD_COEFF';
                return;
            end
            
            if strcmpi(paramName, 'CUBE_COEFF') || strcmpi(paramName, 'sectioncubic')
                result = true;
                path = 'model.inlet.sec_%03d.CUBE_COEFF';
                return;
            end
            
            if ~isempty(obj.bindingHandle)
                [result, paramName] = obj.bindingHandle.hasParameter(paramName);
                path = ['model.adsorption.' paramName];
                return;
            end
            
            path = '';
            
            function res = checkParamExistence(val)
                pos = find(val.structName == '.', 1, 'last');
                sName = val.structName(pos+1:end);
                res = strcmpi(val.propName, paramName) || strcmpi(sName, paramName);
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
            
            obj.addTracking('nComponents', 'model.NCOMP', @(x) int32(x), @(x) double(x));

            % Initial conditions
            obj.addTracking('initialMobileConcentration', 'model.INIT_C');
            obj.addTracking('initialSolidConcentration', 'model.INIT_Q');

            obj.addTracking('bindingParameters', 'model.adsorption');
            obj.addTracking('kineticBindingModel', 'model.adsorption.IS_KINETIC', @(x) int32(x), @(x) x == 1);
            
            % Transport
            obj.addTracking('dispersionColumn', 'model.COL_DISPERSION');
            obj.addTracking('interstitialVelocity', 'model.VELOCITY');

            obj.addTracking('porosityColumn', 'model.COL_POROSITY');
            obj.addTracking('porosityParticle', 'model.PAR_POROSITY');

            obj.addTracking('filmDiffusion', 'model.FILM_DIFFUSION');
            obj.addTracking('diffusionParticle', 'model.PAR_DIFFUSION');
            obj.addTracking('diffusionParticleSurface', 'model.PAR_SURFDIFFUSION');

            % Geometry
            obj.addTracking('columnLength', 'model.COL_LENGTH');
            obj.addTracking('particleRadius', 'model.PAR_RADIUS');

            % Inlet
            obj.addTracking('nInletSections', 'model.inlet.NSEC', @(x) int32(x), @(x) double(x));
            obj.addTracking('sectionTimes', 'model.inlet.SECTION_TIMES');
            obj.addTracking('sectionContinuity', 'model.inlet.SECTION_CONTINUITY', @(x) int32(logical(x)), @(x) x == 1);

            obj.addTracking('extFuncVelocity', 'model.external.EXT_VELOCITY', [], [], true);
            obj.addTracking('extFuncProfile', 'model.external.EXT_PROFILE', [], [], true);
            obj.addTracking('extFuncDelta', 'model.external.EXT_PROF_DELTA', [], [], true);
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
        
        function detectContinuity(obj, threshold)
            % detectContinuity Fills the sectionContinuity with appropriate values by detecting continuous section transitions
            %
            % Parameters:
            %   - threshold: Optional. Threshold for continuity detection (default: 1e-10)

            if (nargin <= 1) || isempty(threshold)
                threshold = 1e-10;
            end

            obj.sectionContinuity = true(obj.nInletSections - 1, 1);

            % Visit each section transition
            for i = 1:obj.nInletSections-1
                
                % Transition from section i to i + 1
                % Section i = [sectionTimes(i), sectionTimes(i+1)]
                % Section i+1 = [sectionTimes(i+1), sectionTimes(i+2)]
                
                % A transition is continuous if and only if the section
                % transition of each component is continuous
                
                for comp = 1:obj.nComponents
                    % Piecewise polynomials are defined on intervals
                    % T_i = [t_i, t_{i+1}] by the formula
                    %   P_i(t) = c_n * (t - t_{i})^n + ... + c_1 * (t - t_{i}) + c_0
                    % Thus, if the transition from section i to i+1 is
                    % continuous, it holds that
                    %   L_i( t_{i+1} - t_{i} ) = L_{i+1}( 0 ),
                    % where L_i( t ) = c_n * t^n + ... + c_1 * t + c_0
                    % is the unshifted polynomial.

                    % Build left polynomial, i.e., on section i
                    leftCoefs = [obj.sectionCubic(comp, i), ...
                        obj.sectionQuadratic(comp, i), ...
                        obj.sectionLinear(comp, i), ...
                        obj.sectionConstant(comp, i)];

                    % Instead of taking the right polynomial, i.e., on
                    % section i+1, we take only the absolute value of this
                    % polynomial. This is valid since the right polynomial
                    % would be evaluated at 0 which only retains the
                    % absolute value.
                    rightVal = obj.sectionConstant(comp, i+1);
                
                    % Compare
                    obj.sectionContinuity(i) = (abs(polyval(leftCoefs, obj.sectionTimes(i+1) - obj.sectionTimes(i)) - rightVal) <= threshold) && obj.sectionContinuity(i);
                    if ~obj.sectionContinuity(i)
                        % Early out if a component is not continuous
                        break;
                    end
                end
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
