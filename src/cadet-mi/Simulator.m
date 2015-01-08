
classdef Simulator < handle
    %SIMULATOR Represents the CADET Simulator
    % Simulates a given Model using a given Discretization and computes 
    % Sensitivities if necessary.
    %
    % There are two available options for calling CADET from MATLAB:
    %   1. MEX
    %   2. HDF5 file exchange
    % If the MEX file is present, it is used to run CADET from within
    % MATLAB. This is the fastest method and preferred over HDF5 file
    % exchange which is used as backup solution.
    % The HDF5 file exchange method takes all settings and writes them to
    % an HDF5 file. After CADET has processed the file (in a separate
    % process, i.e., MATLAB waits for the CADET process to finish), the
    % results (written back to the same HDF5 file) are read in.
    %
    % Copyright: © 2008-2015 Eric von Lieres, Joel Andersson, Andreas Püttmann, Sebastian Schnittert, Samuel Leweke
    %            See the license note at the end of the file.
    
    properties
        model;          % Model to simulate
        discretization; % Discretization of the model
        solverOptions;  % Options controlling some aspects of the solver
        sensitivities;  % Cell array of structs with sensitivities
        sensitivityMethod;  % Method for computation of sensitivies

        % HDF5 Backend settings
        % Path variables guessed by CMake, override if necessary.
        
        % Directory of cadet-cs binary
        binPath = 'bin';
        
        % Directory of generated input / output files (has to be writable). 
        %This will usually point to /tmp (Linux).
        ioPath  = fileparts(tempname());
        
        % Directories of the SUNDIALS, HDF5 and CADET libs (separated by a
        % hash tag '#'). These are only necessary if you are not using a
        % standalone version of CADET (i.e., if you did not build CADET
        % from source)
        libPath = 'lib#sundials/lib#HDF5/lib';
        
        % Suffix of the HDF5 files this class works on (including file
        % extension)
        fileSuffix = '_cadet.h5';
        
        % Name of the cadet-cs binary
        binaryName = 'cadet-cs';
    end

    properties (Dependent)
        solutionTimes; % Vector with time points at which the solution is returned
        nThreads; % Number of threads CADET uses (OpenMP)
    end
    
    properties (Access = 'protected')
        useFileBackend; % Determines whether MEX (false) or HDF5 (true) communication is used
        
        exchangeFile;   % Filename for HDF5 file exchange
        
        sensitivityAccess; % Cell array of structs which cache how to access a sensitive parameter
    end
    
    methods

        function obj = Simulator(model, discretization, solverOptions)
            %SIMULATOR Constructs a Simulator object
            %
            % Parameters:
            %   - model: Optional. Model object
            %   - discretization: Optional. Discretization object
            %   - solverOptions: Optional. Struct with solver options

            % Determine backend
            obj.useFileBackend = (exist(['CadetMex.' mexext], 'file') ~= 3);
            
            if obj.useFileBackend
                % Setup exchange file
                [~, fileName] = fileparts(tempname());
                obj.exchangeFile = fullfile(obj.ioPath, [fileName obj.fileSuffix]);

                % Check if binary exists
                if exist(fullfile(obj.binPath, obj.binaryName), 'file') ~= 2
                    error('CADET:cadetBinaryNotFound', 'The CADET simulator binary could not be found at %s', fullfile(obj.binPath, obj.binaryName));
                end
                
                % Set linker search path for libraries
                if ispc
                    % Prepend libpath only if not present already
                    curEnv = getenv('PATH');
                    obj.libPath = strrep(obj.libPath, '#', ';');
                    if (length(curEnv) < length(obj.libPath)) || any(curEnv(1:length(obj.libPath)) ~= obj.libPath)
                        % Windows uses semicolons (;) to separate the PATH items
                        setenv('PATH', [obj.libPath ';' curEnv]);
                    end
                elseif ismac
                    % Prepend libpath only if not present already
                    curEnv = getenv('DYLD_LIBRARY_PATH');
                    obj.libPath = strrep(obj.libPath, '#', ':');
                    if (length(curEnv) < length(obj.libPath)) || any(curEnv(1:length(obj.libPath)) ~= obj.libPath)
                        setenv('DYLD_LIBRARY_PATH', [obj.libPath ':' curEnv]);
                    end
                elseif isunix
                    % Prepend libpath only if not present already
                    curEnv = getenv('LD_LIBRARY_PATH');
                    obj.libPath = strrep(obj.libPath, '#', ':');
                    if (length(curEnv) < length(obj.libPath)) || any(curEnv(1:length(obj.libPath)) ~= obj.libPath)
                        setenv('LD_LIBRARY_PATH', [obj.libPath ':' curEnv]);
                    end
                end
            end
            
            if nargin >= 1 && ~isempty(model)
                obj.model = model;
            else
                obj.model = [];
            end
            
            if nargin >= 2 && ~isempty(discretization)
                obj.discretization = discretization;
            else
                obj.discretization = [];
            end
            
            if nargin >= 3 && ~isempty(solverOptions)
                obj.solverOptions = solverOptions;
            else
                % Use default values
                obj.solverOptions = Simulator.getDefaultSolverOptions();
            end
            
            obj.sensitivities = [];
            obj.sensitivityMethod = 'ad1';
            obj.sensitivityAccess = [];
        end
        
        function delete(obj)
            % delete Destructs a Simulator object
            
            % Delete the HDF5 exchange file
            if obj.useFileBackend && exist(obj.exchangeFile, 'file') == 2
                delete(obj.exchangeFile);
            end
        end

        function result = simulate(obj, skipValidation, skipSynch, skipExtraction)
            % simulate Simulates the currently set model and returns the results
            %
            % Parameters:
            %   - skipValidation: Optional. Set to true to skip parameter
            %       validation (default: false)
            %   - skipSynch: Optional. Set to true to skip synchronizing
            %       the structs of the Model and Discretization to their
            %       struct representation (default: false)
            %   - skipExtraction: Optional. Set to true to skip the parsing
            %       of the raw CADET output (default: false)
            %
            % Returns: Struct with the fields (see Simulator.extractResults)
            %   - solution
            %       o time: Time points of the solution
            %       o outlet: Matrix with outlet signals (column is
            %           component)
            %       o outletComponents: Vector indicating the component
            %           index of each column of the outlet matrix
            %       o inlet: Matrix with inlet signals (column is
            %           component)
            %       o inletComponents: Vector indicating the component
            %           index of each column of the inlet matrix
            %       o column: Tensor with the raw column solution
            %       o particle: Tensor with the raw particle solution
            %       o flux: Tensor with the raw fluxes
            %       o lastState: Vector with the last state of the solution
            %   - sensitivity
            %       o jacobian: Tensor with the derivatives of the
            %           components with respect to the different
            %           parameters. First dimension contains the
            %           derivative, second the parameter, and third the
            %           component, i.e. jacobian(:, i, j) is the derivative
            %           of component j with respect to parameter i.
            %       o jacParams: Index matrix indicating the parameter of a
            %           "column" in the jacobian tensor. jacParams(i, j) is
            %           the real parameter index of jacobian(:, i, j).
            %       o jacComponents: Index matrix indicating the component
            %           of a "column" in the jacobian tensor. 
            %           jacComponents(i, j) is the real component index of 
            %           jacobian(:, i, j).
            %       o column: Tensor with the raw column derivatives
            %       o particle: Tensor with the raw particle derivatives
            %       o flux: Tensor with the raw flux derivatives
            %       o lastState: Vector with the last state of all sensitivities

            if nargin <= 1 || isempty(skipValidation)
                skipValidation = false;
            end
            
            if nargin <= 2 || isempty(skipSynch)
                skipSynch = false;
            end
            
            if nargin <= 3 || isempty(skipExtraction)
                skipExtraction = false;
            end

            task = obj.prepareSimulation(skipValidation, skipSynch);
            
            % Run simulator
            result = obj.run(task);
            
            % Parse results
            if ~skipExtraction
                result = Simulator.extractResults(result);
            end
        end
        
        function task = prepareSimulation(obj, skipValidation, skipSynch)
            % prepareSimulation Prepares the current model for simulation
            %
            % The prepared simulation can be run by calling Simulator.run
            % with the returned struct.
            %
            % Parameters:
            %   - skipValidation: Optional. Set to true to skip parameter
            %       validation (default: false)
            %   - skipSynch: Optional. Set to true to skip synchronizing
            %       the structs of the Model and Discretization to their
            %       struct representation (default: false)
            %
            % Returns: Internal struct required for simulation
 
            if nargin <= 1 || isempty(skipValidation)
                skipValidation = false;
            end
            
            if nargin <= 2 || isempty(skipSynch)
                skipSynch = false;
            end
            
            % Validate
            if ~skipValidation
                assert(obj.model.checkValues(), 'CADET:invalidParameters', 'Model has invalid parameters or settings');
                assert(obj.discretization.checkValues(), 'CADET:invalidParameters', 'Discretization has invalid parameters or settings');
                assert(obj.checkSolverOptions(), 'CADET:invalidParameters', 'Invalid solver options');
                assert(obj.checkSensitivities(), 'CADET:invalidParameters', 'Invalid sensitivities');
            end
            
            % Combine all parts to one big struct mimicking the CADET file
            % format specs
            if ~skipSynch
                obj.model.synchToStruct();
                obj.discretization.synchToStruct();
            end

            % Assemble task
            task = obj.model.input;
            task.discretization = obj.discretization.input.discretization;
            task.solver = obj.convertSolverOptions();
            task.sensitivity = obj.convertSensitivities();
            
            % Set absolute tolerances for parameter sensitivities
            counter = 0;
            for i = 1:length(obj.sensitivities)
                acc = obj.sensitivityAccess{i};
                curSens = obj.sensitivities{i};

                if ~curSens.active || ~curSens.autoAbsTol
                    continue;
                end
                
                % Taylor path to value in struct
                if curSens.SENS_COMP ~= -1
                    % Attention: SENS_COMP is zero-based
                    s = struct('type', [repmat({'.'}, length(acc.path), 1); {'()'}], 'subs', [acc.path(:); {[{curSens.SENS_COMP + 1}]}]);
                else
                    s = struct('type', repmat({'.'}, length(acc.path), 1), 'subs', acc.path(:));
                end
                
                % Fetch parameter value
                param = subsref(task, s);
                
                % Update absolute tolerance
                pset = sprintf('param_%03d', counter);

                abstol = task.solver.time_integrator.ABSTOL;
                if param ~= 0
                    abstol = abstol ./ abs(param);
                end

                task.sensitivity.(pset).SENS_ABSTOL = abstol;
                counter = counter + 1;
            end
            
        end
        
        function result = runWithParameters(obj, task, params, skipExtraction)
            % runWithParameters Runs the simulator with the given struct and the given parameters
            %
            % Parameters:
            %   - task: Struct created by Simulator.prepareSimulation()
            %   - params: Optional. Vector with values for the parameters
            %       in the order of the sensitivities. If empty, the
            %       simulation is run with the parameters set in the model
            %   - skipExtraction: Optional. Set to true to skip parsing /
            %       extraction of CADET's results (default: false)
            %
            % Returns: Same struct as Simulator.simulate()

            if (nargin <= 2)
                params = [];
            end
            
            if (nargin <= 3) || isempty(skipExtraction)
                skipExtraction = false;
            end
            
            % Set parameters
            assert(isempty(params) || (length(params) == length(obj.sensitivities)), 'CADET:invalidNumberOfParameters', ...
                'The number of parameters does not match the number of configured sensitivities');

            if ~isempty(params)
                counter = 0;
                for i = 1:length(obj.sensitivities)
                    acc = obj.sensitivityAccess{i};
                    curSens = obj.sensitivities{i};

                    % Taylor path to value in struct
                    if curSens.SENS_COMP ~= -1
                        s = struct('type', repmat({'.'}, length(acc.path), 1), 'subs', acc.path(:));
                        
                        % Attention: SENS_COMP is zero-based
                        if (curSens.SENS_SECTION ~= -1) && (length(subsref(task, s)) >= task.model.NCOMP * task.model.inlet.NSEC)
                            % Attention: SENS_SECTION is zero-based
                            s = struct('type', [repmat({'.'}, length(acc.path), 1); {'()'}], 'subs', [acc.path(:); {[{curSens.SENS_SECTION * task.model.NCOMP + (curSens.SENS_COMP + 1)}]}]);
                        else
                            s = struct('type', [repmat({'.'}, length(acc.path), 1); {'()'}], 'subs', [acc.path(:); {[{curSens.SENS_COMP + 1}]}]);                        
                        end
                    else
                        s = struct('type', repmat({'.'}, length(acc.path), 1), 'subs', acc.path(:));
                        if (curSens.SENS_SECTION ~= -1) && (length(subsref(task, s)) >= task.model.inlet.NSEC)
                            % Attention: SENS_SECTION is zero-based
                            s = struct('type', [repmat({'.'}, length(acc.path), 1); {'()'}], 'subs', [acc.path(:); {[{curSens.SENS_SECTION + 1}]}]);
                        end
                    end

                    % Assign new parameter value
                    task = subsasgn(task, s, params(i));

                    if curSens.active

                        % Update absolute tolerance
                        if curSens.autoAbsTol
                            pset = sprintf('param_%03d', counter);

                            abstol = task.solver.time_integrator.ABSTOL;
                            if params(i) ~= 0
                                abstol = abstol ./ abs(params(i));
                            end

                            task.sensitivity.(pset).SENS_ABSTOL = abstol;
                        end

                        counter = counter + 1;
                    end
                end
            end
            
            % Run simulator
            result = obj.run(task);
            
            % Parse results
            if ~skipExtraction
                result = Simulator.extractResults(result);
            end
        end
        
        function result = run(obj, task)
            % run Runs the simulator with the given struct as input
            %
            % Parameters:
            %   - task: Struct mimicking the CADET file format specs
            %
            % Returns: Struct with the raw CADET output

            if obj.useFileBackend
                % HDF5 file exchange

                % Set number of OpenMP threads
                setenv('OMP_NUM_THREADS', num2str(obj.nThreads));

                % Create the HDF5 file
                Simulator.struct2hdf(obj.exchangeFile, task, '/input', [], []);

                % Call simulator
                [failed, msgText] = system([fullfile(obj.binPath, obj.binaryName) ' ' obj.exchangeFile]);
                if failed
                    fprintf('Simulation failed (file %s):\n%s', obj.exchangeFile, msgText);
                    error('CADET:simulationFailed', 'Simulation failed! Check your settings and try again.');
                end

                % Load output from the HDF5 file
                try
                    result = Simulator.hdf2struct(obj.exchangeFile, '/output');
                catch
                    % Failure reading the HDF5 file. Something in CADET has
                    % gone wrong, i.e., simulation failed
                    fprintf('Simulation failed (file %s):\n%s', obj.exchangeFile, msgText);
                    error('CADET:simulationFailed', 'Simulation failed! Check your settings and try again.');
                end
            else
                % MEX
                data.input = task;
                try
                    out = CadetMex(data);
                    result = out.output;
                catch e
                    % Something went wrong
                    error('CADET:simulationFailed', 'Simulation failed! Check your settings and try again.\n%s', e.message);
                end
            end
                        
        end

        function addParameter(obj, name, component, section, sensitive, abstol, fdDelta)
            % addParameter Adds a (sensitive) parameter to the simulation
            %
            % Note that the Model (and its Isotherm if the parameter
            % belongs to it) have to be set to the Simulator before adding
            % parameters.
            %
            % Parameters:
            %   - name: Name of the parameter
            %   - component: Index of the component of the parameter
            %       (one-based) or -1 if independent of component
            %   - section: Index of the section of the parameter
            %       (one-based) or -1 if independent of section
            %   - sensitive: Optional. Set to true to enable computation of
            %       sensitivities with respect to this parameter (default:
            %       true)
            %   - abstol: Absolute tolerance for error estimation. Set to
            %       empty vector to enable automatic tolerance calculation
            %       via the formula 
            %           sovlerOptions.time_integrator.ABSTOL / |paramValue|
            %   - fdDelta: Optional. Step size for finite difference
            %       calculations (default: 0.01)
            
            autoAbsTol = false;
            
            if (nargin <= 2) || isempty(component)
                component = -1;
            end
            
            if (nargin <= 3) || isempty(section)
                section = -1;
            end

            if (nargin <= 4) || isempty(sensitive)
                sensitive = true;
            end
            
            if (nargin <= 5) || isempty(abstol)
                abstol = obj.solverOptions.time_integrator.ABSTOL;
                autoAbsTol = true;
            end
            
            if (nargin <= 6) || isempty(fdDelta)
                fdDelta = 0.01;
            end
            
            % Convert to zero-based
            if component ~= -1
                component = component - 1;
            end
            
            if section ~= -1
                section = section - 1;
            end
            
            actualName = obj.addSensitivityAccess(name, component, section);

            s = [];
            s.SENS_NAME = actualName;
            s.SENS_COMP = component;
            s.SENS_SECTION = section;
            s.SENS_ABSTOL = abstol;
            s.SENS_FD_DELTA = fdDelta;
            s.active = sensitive;
            s.autoAbsTol = autoAbsTol;
            
            obj.sensitivities{end+1} = s;
        end
        
        function setParameters(obj, names, components, sections, sensitives, abstols, fdDeltas)
            % setParameters Sets (sensitive) parameters of the simulation
            %
            % Note that the Model (and its Isotherm if the parameter
            % belongs to it) have to be set to the Simulator before adding
            % parameters.
            %
            % Parameters:
            %   - names: Cell array with parameter names
            %   - components: Array with indices of the parameters'
            %       components (one-based) or -1 if independent of component
            %   - sections: Array with indices of the parameters' sections
            %       (one-based) or -1 if independent of section
            %   - sensitives: Optional. Array indicating whether to enable 
            %       computation of sensitivities with respect to the 
            %       parameter (default: true)
            %   - abstols: Array with absolute tolerances for error
            %       estimation. Set to empty vector or negative value to 
            %       enable automatic tolerance calculation via the formula
            %           sovlerOptions.time_integrator.ABSTOL / |paramValue|
            %   - fdDeltas: Optional. Array with step sizes for finite 
            %       difference calculations (default: 0.01)
            
            autoAbsTol = false(length(names), 1);
            
            if (nargin <= 2) || isempty(components)
                components = -1 * ones(length(names), 1);
            end
            
            if (nargin <= 3) || isempty(sections)
                sections = -1 * ones(length(names), 1);
            end
            
            if (nargin <= 4) || isempty(sensitives)
                sensitives = true(length(names), 1);
            end
            
            if (nargin <= 5) || isempty(abstols)
                abstols = obj.solverOptions.time_integrator.ABSTOL * ones(length(names), 1);
                autoAbsTol = true(length(names), 1);
            end
            
            if (nargin <= 6) || isempty(fdDeltas)
                fdDeltas = 0.01 * ones(length(names), 1);
            end

            assert(length(components) == length(names), 'CADET:invalidParameterSpecs', ...
                'When setting parameters the length of components has to match the length of names');
            assert(length(sections) == length(names), 'CADET:invalidParameterSpecs', ...
                'When setting parameters the length of sections has to match the length of names');
            assert(length(sensitives) == length(names), 'CADET:invalidParameterSpecs', ...
                'When setting parameters the length of sensitives has to match the length of names');
            assert(length(abstols) == length(names), 'CADET:invalidParameterSpecs', ...
                'When setting parameters the length of abstols has to match the length of names');
            assert(length(fdDeltas) == length(names), 'CADET:invalidParameterSpecs', ...
                'When setting parameters the length of fdDeltas has to match the length of names');
            assert(all(components ~= 0), 'CADET:invalidParameterSpecs', ...
                'Found a 0 although the component index is one-based');
            assert(all(sections ~= 0), 'CADET:invalidParameterSpecs', ...
                'Found a 0 although the section index is one-based');
            
            % Convert to zero-based
            components = -1 .* (components == -1) + (components - 1) .* (components >= 0);
            sections = -1 .* (sections == -1) + (sections - 1) .* (sections >= 0);
            
            obj.sensitivityAccess = [];
            obj.sensitivities = cell(length(names), 1);
            for i = 1:length(names)
               
                actualName = obj.addSensitivityAccess(names{i}, components(i), sections(i));

                s = [];
                s.SENS_NAME = actualName;
                s.SENS_COMP = components(i);
                s.SENS_SECTION = sections(i);
                s.SENS_ABSTOL = abstols(i);
                s.SENS_FD_DELTA = fdDeltas(i);
                s.active = sensitives(i);
                s.autoAbsTol = autoAbsTol(i) || (abstols(i) < 0);

                obj.sensitivities{i} = s;
            end
        end
        
        function clearParameters(obj)
            % clearParameters Clears all parameters
            
            obj.sensitivities = [];
            obj.sensitivityAccess = [];
        end
        
        function saveAsHDF5(obj, task, fileName)
            % saveAsHDF5 Saves the simulation task to a cadet-cs compatible HDF5 file
            Simulator.struct2hdf(fileName, task, '/input', [], []);
        end
        
        function val = get.solutionTimes(obj)
            val = obj.solverOptions.USER_SOLUTION_TIMES;
        end
        
        function set.solutionTimes(obj, val)
            assert(all(val >= 0) && all(val(2:end) > val(1:end-1)), 'CADET:invalidTimePoints', ...
                'The solution times have to be nonnegative monotonically increasing');
            
            if isempty(val)
                obj.solverOptions.USER_SOLUTION_TIMES = [];
                obj.solverOptions.WRITE_AT_USER_TIMES = false;
            else
                obj.solverOptions.WRITE_AT_USER_TIMES = true;
                obj.solverOptions.USER_SOLUTION_TIMES = val;
            end
        end
        
        function val = get.nThreads(obj)
            val = obj.solverOptions.NTHREADS;
        end

        function set.nThreads(obj, val)
            assert(val >= 1, 'CADET:invalidThreads', 'The number of threads should be at least 1');
            obj.solverOptions.NTHREADS = val;
        end
        
        function [version, commit] = getVersion(obj)
            % getVersion Returns the version of CADET and commit hash it was built from
            
            if obj.useFileBackend
                % Call simulator
                [~, msgText] = system(fullfile(obj.binPath, obj.binaryName));
                
                % Parse message
                parts = splitstring(msgText, '\n');
                parts = splitstring(parts{6}, ' ');
                
                version = parts{5};
                commit = parts{9};
            else
                % MEX
                [version, commit] = CadetMex();
            end
        end
        
    end
    
    methods (Access = 'protected')
        
        function ok = checkSolverOptions(obj)
            % checkSolverOptions Checks the solverOptions for errors
            ok = true;

            if obj.solverOptions.WRITE_AT_USER_TIMES
                if isempty(obj.solverOptions.USER_SOLUTION_TIMES)
                    ok = false;
                    disp('Error: WRITE_AT_USER_TIMES in solverOptions is set but USER_SOLUTION_TIMES in solverOptions / solutionTimes is empty');
                end
                
                if any(obj.solverOptions.USER_SOLUTION_TIMES(2:end) <= obj.solverOptions.USER_SOLUTION_TIMES(1:end-1))
                    ok = false;
                    disp('Error: USER_SOLUTION_TIMES in solverOptions / solutionTimes has to be monotonically increasing');
                end
                
                if max(obj.solverOptions.USER_SOLUTION_TIMES) > max(obj.model.sectionTimes)
                    ok = false;
                    disp('Error: USER_SOLUTION_TIMES in solverOptions / solutionTimes must not exceed inlet time specified in the model');
                end
            end
            
            if obj.solverOptions.schur_solver.GS_TYPE ~= 0 && obj.solverOptions.schur_solver.GS_TYPE ~= 1
                ok = false;
                disp('Error: GS_TYPE in solverOptions.schur_solver has to be either 0 or 1');
            end
            
            if obj.solverOptions.schur_solver.SCHUR_SAFETY < 0.0
                ok = false;
                disp('Error: SCHUR_SAFETY in solverOptions.schur_solver has to be nonnegative');
            end
            
            if (obj.solverOptions.schur_solver.MAX_KRYLOV < 0.0) || (obj.solverOptions.schur_solver.MAX_KRYLOV > obj.discretization.nCellsColumn)
                ok = false;
                disp(['Error: MAX_KRYLOV in solverOptions.schur_solver has to be between 0 and ' num2str(obj.discretization.nCellsColumn)]);
            end
            
            if obj.solverOptions.schur_solver.MAX_RESTARTS < 0.0
                ok = false;
                disp('Error: MAX_RESTARTS in solverOptions.schur_solver has to be nonnegative');
            end
            
            if obj.solverOptions.time_integrator.ABSTOL <= 0.0
                ok = false;
                disp('Error: ABSTOL in solverOptions.time_integrator has to be positive');
            end
            
            if obj.solverOptions.time_integrator.RELTOL < 0.0
                ok = false;
                disp('Error: RELTOL in solverOptions.time_integrator has to be nonnegative');
            end
            
            if obj.solverOptions.time_integrator.INIT_STEP_SIZE < 0.0
                ok = false;
                disp('Error: INIT_STEP_SIZE in solverOptions.time_integrator has to be nonnegative');
            end
            
            if obj.solverOptions.time_integrator.MAX_STEPS < 0.0
                ok = false;
                disp('Error: MAX_STEPS in solverOptions.time_integrator has to be nonnegative');
            end
        end
        
        function so = convertSolverOptions(obj)
            % convertSolverOptions Converts the solverOptions to a struct complying with CADET's file format specs
            
            so = obj.solverOptions;
            
            % Logicals
            so.PRINT_PROGRESS               = int32(logical(so.PRINT_PROGRESS));
            so.PRINT_STATISTICS             = int32(logical(so.PRINT_STATISTICS));
            so.PRINT_TIMING                 = int32(logical(so.PRINT_TIMING));
            so.PRINT_PARAMLIST              = int32(logical(so.PRINT_PARAMLIST));
            so.PRINT_CONFIG                 = int32(logical(so.PRINT_CONFIG));
            so.USE_ANALYTIC_JACOBIAN        = int32(logical(so.USE_ANALYTIC_JACOBIAN));
            so.WRITE_AT_USER_TIMES          = int32(logical(so.WRITE_AT_USER_TIMES));
            so.WRITE_SOLUTION_TIMES         = int32(logical(so.WRITE_SOLUTION_TIMES));
            so.WRITE_SOLUTION_COLUMN_OUTLET = int32(logical(so.WRITE_SOLUTION_COLUMN_OUTLET));
            so.WRITE_SOLUTION_COLUMN_INLET  = int32(logical(so.WRITE_SOLUTION_COLUMN_INLET));
            so.WRITE_SOLUTION_ALL           = int32(logical(so.WRITE_SOLUTION_ALL));
            so.WRITE_SOLUTION_LAST          = int32(logical(so.WRITE_SOLUTION_LAST));
            so.WRITE_SENS_COLUMN_OUTLET     = int32(logical(so.WRITE_SENS_COLUMN_OUTLET));
            so.WRITE_SENS_ALL               = int32(logical(so.WRITE_SENS_ALL));
            so.WRITE_SENS_LAST              = int32(logical(so.WRITE_SENS_LAST));
            so.NTHREADS                     = int32(so.NTHREADS);
            
            % Integers
            so.schur_solver.GS_TYPE         = int32(so.schur_solver.GS_TYPE);
            so.schur_solver.MAX_KRYLOV      = int32(so.schur_solver.MAX_KRYLOV);
            so.schur_solver.MAX_RESTARTS    = int32(so.schur_solver.MAX_RESTARTS);

            so.time_integrator.MAX_STEPS    = int32(so.time_integrator.MAX_STEPS);
        end
        
        function ok = checkSensitivities(obj, useModel)
            % checkSensitivities Checks the sensitivities for errors
            %
            % Parameters:
            %   - useModel: Optional. Set to true to use values from the
            %       Model for validation (default: true)
            %
            % Returns true if the sensitivities are valid
            
            if (nargin <= 1) || isempty(useModel)
                useModel = true;
            end
            
            ok = true;

            if ~strcmp(obj.sensitivityMethod, 'ad1') && ...
                    ~strcmp(obj.sensitivityMethod, 'fd1') && ...
                    ~strcmp(obj.sensitivityMethod, 'fd2') && ...
                    ~strcmp(obj.sensitivityMethod, 'fd4')
                ok = false;
                disp('Error: sensitivityMethod has to be one of "ad1", "fd1", "fd2", and "fd4"');
            end
            
            for i = 1:length(obj.sensitivities)
                curSens = obj.sensitivities{i};
                
                if ~isfield(curSens, 'SENS_NAME') || isempty(curSens.SENS_NAME) || ~ischar(curSens.SENS_NAME)
                    ok = false;
                    disp(['Error: In sensitivity ' num2str(i) ': Field SENS_NAME has to be a nonempty string']);
                end
                
                if ~isfield(curSens, 'SENS_SECTION') || isempty(curSens.SENS_SECTION) || ~isscalar(curSens.SENS_SECTION) || (curSens.SENS_SECTION < -1)
                    ok = false;
                    disp(['Error: In sensitivity ' num2str(i) ': Field SENS_SECTION has to be a scalar of value at least -1']);
                end
                
                if ~isfield(curSens, 'SENS_COMP') || isempty(curSens.SENS_COMP) || ~isscalar(curSens.SENS_COMP) || (curSens.SENS_COMP < -1)
                    ok = false;
                    disp(['Error: In sensitivity ' num2str(i) ': Field SENS_COMP has to be a scalar of value at least -1']);
                end
                
                if ~isfield(curSens, 'autoAbsTol') || isempty(curSens.autoAbsTol)
                    ok = false;
                    disp(['Error: In sensitivity ' num2str(i) ': Field autoAbsTol has to be a scalar logical']);
                end
                
                if ~isfield(curSens, 'active') || isempty(curSens.active)
                    ok = false;
                    disp(['Error: In sensitivity ' num2str(i) ': Field active has to be a scalar logical']);
                end
                
                if isfield(curSens, 'active') && curSens.active
                    % Only check fields if parameter is sensitive
                    
                    if ~isfield(curSens, 'SENS_ABSTOL') || isempty(curSens.SENS_ABSTOL) || ~isscalar(curSens.SENS_ABSTOL) || (curSens.SENS_ABSTOL <= 0.0)
                        ok = false;
                        disp(['Error: In sensitivity ' num2str(i) ': Field SENS_ABSTOL has to be a scalar of value greater than 0.0']);
                    end

                    if ~isfield(curSens, 'SENS_FD_DELTA') || isempty(curSens.SENS_FD_DELTA) || ~isscalar(curSens.SENS_FD_DELTA) || (curSens.SENS_FD_DELTA <= 0.0)
                        ok = false;
                        disp(['Error: In sensitivity ' num2str(i) ': Field SENS_FD_DELTA has to be a scalar of value greater than 0.0']);
                    end
                end
                
                if useModel && ~isempty(obj.model)
                    
                    if isfield(curSens, 'SENS_SECTION') && (curSens.SENS_SECTION >= obj.model.nInletSections)
                        ok = false;
                        disp(['Error: In sensitivity ' num2str(i) ': Field SENS_SECTION has to be smaller than nInletSections (' num2str(obj.model.nInletSections) ')']);
                    end
                    
                    if isfield(curSens, 'SENS_COMP') && (curSens.SENS_COMP >= obj.model.nComponents)
                        ok = false;
                        disp(['Error: In sensitivity ' num2str(i) ': Field SENS_COMP has to be smaller than nComponents (' num2str(obj.model.nComponents) ')']);
                    end
                    
                    if isfield(curSens, 'SENS_NAME') && ~obj.model.hasParameter(curSens.SENS_NAME)
                        ok = false;
                        disp(['Error: In sensitivity ' num2str(i) ': The parameter "' curSens.SENS_NAME '" does not exist']);
                    end
                    
                end
            end
            
            if ~isempty(obj.model)
                obj.rebuildSensitivityAccess();
            end
        end

        function sens = convertSensitivities(obj)
            % convertSensitivities Converts the sensitivities to a struct complying with CADET's file format specs
            
            sens = [];
            
            counter = 0;
            for i = 1:length(obj.sensitivities)
                
                if ~obj.sensitivities{i}.active
                    continue;
                end
                
                fn = sprintf('param_%03d', counter);
                sens.(fn) = obj.sensitivities{i};
                
                sens.(fn).SENS_COMP = int32(sens.(fn).SENS_COMP);
                sens.(fn).SENS_SECTION = int32(sens.(fn).SENS_SECTION);
                sens.(fn) = rmfield(sens.(fn), 'active');
                sens.(fn) = rmfield(sens.(fn), 'autoAbsTol');

                counter = counter + 1;
            end
            
            sens.SENS_METHOD = obj.sensitivityMethod;
            sens.NSENS = int32(counter);
        end
        
        function actualName = addSensitivityAccess(obj, name, compIdx, secIdx)
            % addSensitivityAccess Adds an entry to the sensitivityAccess list
            %
            % Parameters:
            %   - name: Name of the sensitive parameter
            %   - compIdx: Index of the component (zero-based)
            %   - secIdx: Index of the section (zero-based)
            %
            % Returns: The actual parameter name according to CADET's file
            %   format specs
            
            if isempty(obj.model)
                actualName = name;
                return;
            end
            
            [result, path] = obj.model.hasParameter(name);
            
            if result
                pos = find(path == '.', 1, 'last');
                if ~isempty(pos)
                    actualName = path(pos+1:end);
                else
                    actualName = path;
                end

                s = [];
                s.path = splitstring(sprintf(path, secIdx), '.');
                s.name = actualName;
            else
                s = [];
                s.path = path;
                s.name = name;
            end
            
            obj.sensitivityAccess{end+1} = s;
        end
        
        function rebuildSensitivityAccess(obj)
            % rebuildSensitivityAccess Rebuilds the sensitivityAccess list

            if isempty(obj.model)
                return;
            end
            
            obj.sensitivityAccess = [];
            for i = 1:length(obj.sensitivities)
                curSens = obj.sensitivities{i};
                actualName = obj.addSensitivityAccess(curSens.SENS_NAME, curSens.SENS_COMP, curSens.SENS_SECTION);
                obj.sensitivities{i}.SENS_NAME = actualName;
            end
        end
         
    end
    
    methods (Static = true, Access = 'public')

        function ext = extractResults(result)
            % extractResults Extracts the returned data from a CADET simulation and sorts it into a struct
            %
            % Parameters:
            %   - result: Struct with the output of CADET
            %
            % Returns: A struct with the parsed / converted output
            %   - solution
            %       o time: Time points of the solution
            %       o outlet: Matrix with outlet signals (column is
            %           component)
            %       o outletComponents: Vector indicating the component
            %           index of each column of the outlet matrix
            %       o inlet: Matrix with inlet signals (column is
            %           component)
            %       o inletComponents: Vector indicating the component
            %           index of each column of the inlet matrix
            %       o column: Tensor with the raw column solution
            %       o particle: Tensor with the raw particle solution
            %       o flux: Tensor with the raw fluxes
            %       o lastState: Vector with the last state of the solution
            %   - sensitivity
            %       o jacobian: Tensor with the derivatives of the
            %           components with respect to the different
            %           parameters. First dimension contains the
            %           derivative, second the parameter, and third the
            %           component, i.e. jacobian(:, i, j) is the derivative
            %           of component j with respect to parameter i.
            %       o jacParams: Index matrix indicating the parameter of a
            %           "column" in the jacobian tensor. jacParams(i, j) is
            %           the real parameter index of jacobian(:, i, j).
            %       o jacComponents: Index matrix indicating the component
            %           of a "column" in the jacobian tensor. 
            %           jacComponents(i, j) is the real component index of 
            %           jacobian(:, i, j).
            %       o column: Tensor with the raw column derivatives
            %       o particle: Tensor with the raw particle derivatives
            %       o flux: Tensor with the raw flux derivatives
            %       o lastState: Vector with the last state of all sensitivities
            
            ext = [];

            solNames = fieldnames(result.solution);
            idxOutlet = find(cellfun(@(x) (length(x) >= 22) && strcmp(x(1:22), 'SOLUTION_COLUMN_OUTLET'), solNames));
            idxInlet = find(cellfun(@(x) (length(x) >= 21) && strcmp(x(1:21), 'SOLUTION_COLUMN_INLET'), solNames));
            
            idxTime = find(cellfun(@(x) strcmp(x, 'SOLUTION_TIMES'), solNames), 1);
            idxColumn = find(cellfun(@(x) strcmp(x, 'SOLUTION_COLUMN'), solNames), 1);
            idxParticle = find(cellfun(@(x) strcmp(x, 'SOLUTION_PARTICLE'), solNames), 1);
            idxBoundary = find(cellfun(@(x) strcmp(x, 'SOLUTION_BOUNDARY'), solNames), 1);
            idxLastState = find(cellfun(@(x) strcmp(x, 'SOLUTION_LAST'), solNames), 1);
            
            ext.solution = [];
            ext.solution.time = [];
            
            ext.solution.outlet = [];
            ext.solution.outletComponents = [];
            ext.solution.inlet = [];
            ext.solution.inletComponents = [];
            
            ext.solution.column = [];
            ext.solution.particle = [];
            ext.solution.flux = [];
            ext.solution.lastState = [];

            if ~isempty(idxTime)
                ext.solution.time = result.solution.SOLUTION_TIMES;
            end

            if ~isempty(idxColumn)
                ext.solution.column = result.solution.SOLUTION_COLUMN;
            end

            if ~isempty(idxParticle)
                ext.solution.particle = result.solution.SOLUTION_PARTICLE;
            end
            
            if ~isempty(idxBoundary)
                ext.solution.flux = result.solution.SOLUTION_BOUNDARY;
            end
            
            if ~isempty(idxLastState)
                ext.solution.lastState = result.solution.SOLUTION_LAST;
            end

            % Extract outlet signals
            if ~isempty(idxOutlet)
                % Convert to component indices
                idxOutlet = arrayfun(@(x) str2double(solNames{x}(29:end)), idxOutlet);
                idxOutlet = sort(idxOutlet);
                
                % Preallocate
                dset = sprintf('SOLUTION_COLUMN_OUTLET_COMP_%03d', idxOutlet(1));
                ext.solution.outlet = zeros(length(result.solution.(dset)), length(idxOutlet));
                ext.solution.outletComponents = zeros(length(idxOutlet), 1);

                for comp = 1:length(idxOutlet)
                    dset = sprintf('SOLUTION_COLUMN_OUTLET_COMP_%03d', idxOutlet(comp));
                    ext.solution.outlet(:, comp) = result.solution.(dset);
                    ext.solution.outletComponents(comp) = idxOutlet(comp);
                end
            end

            % Extract inlet signals
            if ~isempty(idxInlet)
                % Convert to component indices
                idxInlet = arrayfun(@(x) str2double(solNames{x}(28:end)), idxInlet);
                idxInlet = sort(idxInlet);
                
                % Preallocate
                dset = sprintf('SOLUTION_COLUMN_INLET_COMP_%03d', idxInlet(1));
                ext.solution.inlet = zeros(length(result.solution.(dset)), length(idxInlet));
                ext.solution.inletComponents = zeros(length(idxInlet), 1);

                for comp = 1:length(idxInlet)
                    dset = sprintf('SOLUTION_COLUMN_INLET_COMP_%03d', idxInlet(comp));
                    ext.solution.inlet(:, comp) = result.solution.(dset);
                    ext.solution.inletComponents(comp) = idxInlet(comp);
                end
            end

            % Sensitivities
            ext.sensitivity = [];
            ext.sensitivity.jacobian = [];
            ext.sensitivity.jacParams = [];
            ext.sensitivity.jacComponents = [];

            ext.sensitivity.column = [];
            ext.sensitivity.particle = [];
            ext.sensitivity.flux = [];
            ext.sensitivity.lastState = [];
            
            if ~isfield(result, 'sensitivity')
                return;
            end
            
            solNames = fieldnames(result.sensitivity);

            idxColumn = find(cellfun(@(x) strcmp(x, 'SENS_COLUMN'), solNames), 1);
            idxParticle = find(cellfun(@(x) strcmp(x, 'SENS_PARTICLE'), solNames), 1);
            idxBoundary = find(cellfun(@(x) strcmp(x, 'SENS_BOUNDARY'), solNames), 1);
            idxLastState = find(cellfun(@(x) strcmp(x, 'SENS_LAST'), solNames), 1);
            
            if ~isempty(idxColumn)
                ext.sensitivity.column = result.sensitivity.SENS_COLUMN;
            end
            
            if ~isempty(idxParticle)
                ext.sensitivity.particle = result.sensitivity.SENS_PARTICLE;
            end
            
            if ~isempty(idxBoundary)
                ext.sensitivity.flux = result.sensitivity.SENS_BOUNDARY;
            end

            if ~isempty(idxLastState)
                ext.sensitivity.lastState = result.sensitivity.SENS_LAST;
            end
            
            idxParams = find(cellfun(@(x) (length(x) >= 5) && strcmp(x(1:5), 'param'), solNames));
            if ~isempty(idxParams)
                % Convert to param indices
                idxParams = arrayfun(@(x) str2double(solNames{x}(7:end)), idxParams);
                idxParams = sort(idxParams);
                
                pset = sprintf('param_%03d', idxParams(1));
                solNames = fieldnames(result.sensitivity.(pset));
                
                % Get component indices
                idxComps = cellfun(@(x) str2double(x(25:end)), solNames);
                idxComps = sort(idxComps);
                
                % Preallocate
                dset = sprintf('SENS_COLUMN_OUTLET_COMP_%03d', idxComps(1));
                ext.sensitivity.jacobian = zeros([length(result.sensitivity.(pset).(dset)), length(idxParams), length(idxComps)]);
                ext.sensitivity.jacParams = zeros(length(idxParams), length(idxComps));
                ext.sensitivity.jacComponents = zeros(length(idxParams), length(idxComps));
                
                for i = 1:length(idxParams)
                    pset = sprintf('param_%03d', idxParams(i));
                    for j = 1:length(idxComps)
                        dset = sprintf('SENS_COLUMN_OUTLET_COMP_%03d', idxComps(j));
                        ext.sensitivity.jacobian(:, i, j) = result.sensitivity.(pset).(dset);
                        ext.sensitivity.jacParams(i, j) = idxParams(i);
                        ext.sensitivity.jacComponents(i, j) = idxComps(j);
                    end
                end
            end
        end
        
    end
    
    methods (Static = true, Access = 'private')
        
        function so = getDefaultSolverOptions()
            % getDefaultSolverOptions Returns default solver options
            so = [];

            so.PRINT_PROGRESS               = false;
            so.PRINT_STATISTICS             = false;
            so.PRINT_TIMING                 = false;
            so.PRINT_PARAMLIST              = false;
            so.PRINT_CONFIG                 = false;
            so.USE_ANALYTIC_JACOBIAN        = true;
            so.WRITE_AT_USER_TIMES          = false;
            so.USER_SOLUTION_TIMES          = [];
            so.LOG_LEVEL                    = 'ERROR';
            so.WRITE_SOLUTION_TIMES         = true;
            so.WRITE_SOLUTION_COLUMN_OUTLET = true;
            so.WRITE_SOLUTION_COLUMN_INLET  = false;
            so.WRITE_SOLUTION_ALL           = false;
            so.WRITE_SOLUTION_LAST          = false;
            so.WRITE_SENS_COLUMN_OUTLET     = true;
            so.WRITE_SENS_ALL               = false;
            so.WRITE_SENS_LAST              = false;
            so.NTHREADS                     = 1;
             
            % Schur solver
            so.schur_solver = [];
            so.schur_solver.GS_TYPE                  = 1;
            so.schur_solver.MAX_KRYLOV               = 0;
            so.schur_solver.MAX_RESTARTS             = 0;
            so.schur_solver.SCHUR_SAFETY             = 1e-8;
            
            % Time integrator
            so.time_integrator.ABSTOL                = 1e-8;
            so.time_integrator.RELTOL                = 0.0;
            so.time_integrator.INIT_STEP_SIZE        = 1e-6;
            so.time_integrator.MAX_STEPS             = 10000;
        end
        
        function struct2hdf(fileName, strct, path, nodelete, legacy_hdf5)
            % struct2hdf Writes a (nested) struct to an HDF5 file
            %
            % Parameters:
            %   - fileName: Name of the HDF5 file
            %   - strct: Struct to write
            %   - path: Optional. Root node in the HDF5 file (default: '/')
            %   - nodelete: Optional. Controls whether the file is deleted 
            %       (false, default) before it is written to, or if data is
            %       appended (true)
            %   - legacy_hdf5: Optional. Determines whether to use legacy
            %       HDF5 functions (default: true)
            
            % If not provided, set default input values
            if (nargin <= 4) || isempty(legacy_hdf5)
                legacy_hdf5 = true;
            end
            if (nargin <= 3) || isempty(nodelete)
                nodelete = 0;
            end
            if (nargin <= 2) || isempty(path)
                path = '';
            end
            
            % Remove file if existent
            if ~isempty(dir(fileName)) && ~nodelete
                delete(fileName);
            end
            
            fields = fieldnames(strct);
            
            % Loop over all fields in the struct
            for f = fields'
                newpath = [path '/' f{1}];
                curfield = strct.(f{1});
                
                % If current field is a struct, call this function recursively
                if isstruct(curfield)
                    Simulator.struct2hdf(fileName, curfield, newpath, 1, legacy_hdf5);
                    
                    % If current field is a string, write it
                elseif ischar(curfield)
                    mode = 'append';
                    if isempty(dir(fileName)); mode = 'overwrite'; end
                    hdf5write(fileName, newpath, curfield, 'WriteMode', mode);
                    
                    % If current field is an integer, write it as int32
                elseif isinteger(curfield)
                    if legacy_hdf5
                        mode = 'append';
                        if isempty(dir(fileName)); mode = 'overwrite'; end
                        hdf5write(fileName, newpath, int32(curfield), 'WriteMode', mode);
                    else
                        if isempty(curfield)
                            curfield = [int32(0)]; % Dummy
                        end
                        h5create(fileName, newpath, size(curfield), 'Datatype', 'int32');
                        h5write(fileName, newpath, curfield);
                    end
                        
                    % If current field is a double, write it as extendible dataset
                else % isdouble
                    if legacy_hdf5
                        mode = 'append';
                        if isempty(dir(fileName)); mode = 'overwrite'; end
                        hdf5write(fileName, newpath, curfield, 'WriteMode', mode);
                    else
                        h5create(fileName, newpath, [Inf Inf], 'Datatype', 'double', 'ChunkSize', size(curfield));
                        h5write(fileName, newpath, curfield, [1 1], size(curfield));
                    end
                end
            end
            
        end
        
        function h = hdf2struct(filename, path)
            % hdf2struct Converts a HDF5 file to a nested Matlab struct
            %
            % Parameters:
            %   - filename: Name of the HDF5 file
            %   - path: Optional. Path to root node from which conversion
            %       starts (default: '/')
            %
            % Returns: (Nested) struct mimicking the HDF5 file
            
            if nargin <= 1 || isempty(path)
                path = '/'; 
            end;
            
            h = [];
            h.hdf = [];

            % Open the HDF5 file
            fid = H5F.open(filename,'H5F_ACC_RDONLY','H5P_DEFAULT');
            
            % Open the group in path
            gid = H5G.open(fid, path);
            
            % Iterate over every entry in the HDF5 file and call 'operate' on it
            [status, h] = H5O.visit(gid, 'H5_INDEX_NAME', 'H5_ITER_NATIVE', @operate, h);
            
            % Close the group
            H5G.close(gid);
            
            % Close the HDF5 file
            H5F.close(fid);

            assert(status ~= 1, 'CADET:hdf2structConversion', 'Error in hdf2struct.');
            h = h.hdf;
        end
        
    end
end


function [status, h] = operate(obj, name, h)
%OPERATE Converts every dataset into a matlab struct
    
    % Quick return on root group
    if strcmpi(name, '.')
        status = 0;
        return
    end

    % Check for dataset or group by case - DIRTY HACK!
    % But checking object info was not successful :-(
    last_name = regexp(name, '\w*$', 'match', 'once');
    if strcmp(last_name, upper(last_name))
        typestr = 'Dataset';
    else
        typestr = 'Group';
    end

    switch typestr
        case 'Group'

            % Store group name
            h.group = strrep(name, '/', '.');
            
        case 'Dataset'
            
            % Extract the name of the dataset
            dset_name = regexp(name, '\w*$', 'match', 'once');
            
            dataset = [];
            try
                % Open the dataset
                dataset = H5D.open(obj, name, 'H5P_DEFAULT');

                % Read data and store in handles.hdf structure
                path = splitstring(h.group, '.');
                path = struct('type', repmat({'.'}, length(path) + 1, 1), 'subs', [path(:); {dset_name}]);
                h.hdf = subsasgn(h.hdf, path, H5D.read(dataset, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT'));
           
                H5D.close(dataset);
            catch
                % Close dataset if still open
                if ~isempty(dataset)
                    H5D.close(dataset);
                end
                
                % Return error and stop visiting
                status = 1;
                return;
            end
            
        otherwise
            error('CADET:hdf2structConversion', 'Non-reachable code has been reached');
    end
    status = 0;
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
