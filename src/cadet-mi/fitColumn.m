function [outParams, res, success] = fitColumn(fitData, initParams, loBound, upBound, logScale, parameterTransform, quietMode, optOptions )
%FITCOLUMN Fits a model to the given datasets.
%
% [outParams, res, success] = fitColumn(fitData, initParams, loBound, upBound, logScale, quietMode )
% [outParams, res, success] = fitColumn(fitData, initParams, loBound, upBound, logScale, parameterTransform, quietMode )
% [outParams, res, success] = fitColumn(fitData, initParams, loBound, upBound, logScale, parameterTransform, quietMode, optOptions )
%
% Parameters:
%   - fitData: Array of structs, one for each experiment. Each struct has
%       the fields:
%           o tOut: Outlet time points
%           o outMeas: Outlet measurements, (n x 1)-vector
%           o idxSum: Indices (1-based) of summed / observed outlet signals
%           o sim: Simulator object ready for simulation
%           o task: Executable simulation task
%           o links: Cell array with vectors specifying the indices
%               (one-based!) of the fits in the fitData the parameters are
%               linked to. Linked parameters are estimated simultaneously, 
%               i.e., they have the same value for the linked fits.
%           o joins: Cell array of vectors indicating parameters which are 
%               joined to this parameter. Joined parameters are handled as 
%               one parameter within this fit. The vectors contain the
%               indices (one-based!) of the parameters joined to the 
%               current one.
%           o weight: Optional. Global weight of the experiments in a set
%               of experiments.
%           o weightComponent: Optional. Weight of each observed component
%               in this experiment.
%           o weightSignal: Optional. Weight for each data point of the
%               measurements. Defaults to 1.0 for each data point if
%               neglected. Has to have the same size as outMeas.
%   - initParams: Vector with initial parameters for the optimization. The
%       order of the parameters is in the order of the fitData's paramInfo.
%       Linked parameters are only given once at their first (possibly 
%       linked) occurrence.
%   - loBound: Vector with lower bounds of parameters for optimization. Can
%       be empty to indicate no lower bound. See initParams for ordering.
%   - upBound: Vector with upper bounds of parameters for optimization. Can
%       be empty to indicate no upper bound. See initParams for ordering.
%   - logScale: Optional. Vector of booleans determining whether a
%       parameter is transformed to log scale in optimization. Performance
%       and robustness of most optimizers are substantially improved if all
%       parameters have the same scale. Scaling can be enabled/disabled for
%       each parameter by setting logScale to true/false (scalar). See
%       initParams for ordering. Default is all true.
%   - parameterTransform: Optional. Struct which specifies a parameter
%       transformation. Fields:
%           o postLogTransform: Set to true to apply this transform after
%               the logarithmic transformation, i.e., the final transform
%               is
%                   p_opt = transform( log( p_sim ) ),
%               where p_sim is the parameter vector given to the simulator
%               and p_opt is the parameter vector used by the optimizer.
%               If set to false, the ordering is
%                   p_opt = log( transform( p_sim ) ).
%           o transform: Function handle to the function applying the
%               transformation from simulator parameter space to optimizer
%               parameter space
%           o invTransform: Function handle to the function applying the
%               transformation from optimizer parameter space to simulator
%               parameter space, i.e., the inverse of transform()
%           o chainRuleInvTransform: Function handle to the function
%               applying the chain rule to the Jacobian of the system. The
%               function is given an (N x nParams)-matrix and the current
%               point in simulator space (with respect to potential
%               logarithmization). It should then return a matrix of the
%               same size as the input matrix with the chain rule applied.
%               The ordering of the columns in the Jacobian and the current
%               parameters is the same as for initParams.
%           o transformBounds: Optional. Set to true to transform the upper
%               and lower bounds of the optimization problem. Defaults to
%               true.
%   - quietMode: Optional. Enables quiet mode, i.e. disables plots and
%       iteration monitoring. Defaults to false (quiet mode disabled).
%   - optOptions: Optional. Additional optimizer settings for the chosen
%       optimizer (MATLAB's lsqnonlin at the moment). The settings are
%       applied last and may overwrite options set by fitColumn().
%
% Links and joins: While links combine parameters across several fits,
%   joins combine parameters within one fit. Only parameters sharing the
%   same name, component, and section are linkable. On the contrary, any
%   two or more parameters in one fit can be joined. Joined parameters are
%   declared like normal parameters, but come after normal ("master")
%   parameters. The vectors in the joins cell array give the indices of the
%   parameters joined to the corresponding master parameter. Note that only
%   master parameters can be linked, i.e., a link to a joined ("slave") 
%   parameter is not recognised.
%
% Log Scaling: Log scaling works by taking the variables x and transforming
%   them to v = log(x). Therefore, the forward simulation function f is
%   transformed to F(v) = f(exp(v)) = f(x). This requires the incovation of
%   the chain rule to compute the Jacobian of F:
%       F'(v) = f'( exp(v) ) * exp(v).
%   Negative parameters are handled by taking v = log(-x) and, hence,
%       F(v) = f(-exp(v)) = f(x).
%   The negative parameters are indicated by a negative initial value.
%   Since the exp() function is non-negative, log scaling introduces a
%   natural boundary condition on the parameters, i.e., v >= 0 for positive
%   parameters and v <= 0 for negative parameters. Therefore, 0 is always a
%   lower or upper bound.
%
%   Important: Initial parameters and boundaries are always specified in
%       normal scale.
%
%   In future implementations there might be more transformations which are
%   more suited to the task at hand. For example, transforming with tan()
%   allows for sign changes during optimization and reduces the scale.
%   Other bound-aware transforms are given by 
%       log(x - lower), log(upper - x), logit( (x-lower) / (upper-lower)),
%       where logit(x) = log(x / (1-x)).
%
% Weighting: There are two different types of weights: Global weights
%   applied to an experiment in a set of experiments and local weights
%   applied to each observed component or data point in one single 
%   experiment. Both are optional and can be used separately and 
%   independent from each other. If no weights are given, weighting is 
%   disabled, i.e., every element receives the weight 1.
%
%   The weights are not normalized by default, this has to be done by the
%   user. Also note that the weights are multiplicative, i.e.,
%          weight * (simulated - measured)
%   is computed.
%
% Returns:
%   - outParams: Optimized parameters in the same order as initParams.
%   - res: Residual value returned by the optimizer or -1 in case of
%       failure.
%   - success: True if the optimizer returned successful, otherwise false.
%
% Copyright: © 2008-2015 Eric von Lieres, Joel Andersson, Andreas Püttmann, Sebastian Schnittert, Samuel Leweke
%            See the license note at the end of the file.

    if nargin <= 2
        loBound = [];
    end
    if nargin <= 3
        upBound = [];
    end

    if (nargin <= 4) || isempty(logScale)
        logScale = true(size(initParams));
    end
    if numel(logScale) == 1
        logScale = ones(size(initParams)) .* logScale;
    end
    
    customTransform = true;
    if (nargin <= 5) || isempty(parameterTransform)
        customTransform = false;
    end

    if (nargin <= 6) || isempty(quietMode)
        quietMode = false;
    end

    if (nargin == 6) && islogical(parameterTransform) && (numel(parameterTransform) == 1) && ~isstruct(parameterTransform)
        quietMode = parameterTransform;
        customTransform = false;
        parameterTransform = [];
    end

    if nargin <= 7
        optOptions = [];
    end

    enablePlot = ~quietMode;  % Enables or disables plots in each iteration
    enableWeightedPlot = enablePlot;
    
    if customTransform && (~isfield(parameterTransform, 'transformBounds') || isempty(parameterTransform.transformBounds))
        parameterTransform.transformBounds = true;
    end

    % Convert param vector from column (n x 1) to row (1 x n)
    if size(initParams, 1) > 1
        initParams = initParams';
    end

    % Check fitData for errors
    [hasError, fitData, numJoins] = checkInput(fitData);
    if hasError
        disp('Aborting due to errors.');
        outParams = initParams;
        res = -1;
        return;
    end

    % Prepare conversions between global and local parameters in two passes
    % First pass: Resolve linked parameters and build local to global map
    ctrId = 1;
    maxParams = max(arrayfun(@(idx) length(fitData{idx}.sim.sensitivities) - numJoins(idx), 1:length(fitData)));
    localToGlobal = zeros(length(fitData), maxParams);  % Maps (idxFit, idxLocalParam) -> idxGlobalParam
    for i = 1:length(fitData)
        pi = fitData{i}.sim.sensitivities;
        for j = 1:length(pi) - numJoins(i)
            % Get parameter signature
            comp = pi{j}.SENS_COMP;
            sec = pi{j}.SENS_SECTION;
            pn = pi{j}.SENS_NAME;

            % Check if an Id has been assigned to a linked parameter
            if localToGlobal(i, j) <= 0
                linkId = 0;
                if isfield(fitData{i}, 'links') && ~isempty(fitData{i}.links)
                    % Find indices of linked parameters
                    [linkFit, linkParam] = getIndicesOfLinkedParams(fitData, pn, comp, sec, numJoins, fitData{i}.links{j});

                    % Assign ID to linked parameters
                    for lk = 1:length(linkFit)
                        val = localToGlobal(linkFit(lk), linkParam(lk));
                        if val > 0
                            linkId = val;
                            break;
                        end
                    end
                end
                
                % If a linked parameter has an Id, use it
                if linkId > 0
                    localToGlobal(i, j) = linkId;
                end
            end
            
            if localToGlobal(i, j) <= 0
                % Assign an ID
                localToGlobal(i, j) = ctrId;
                ctrId = ctrId + 1;
            end
            
            if isfield(fitData{i}, 'links') && ~isempty(fitData{i}.links)
                % Find indices of linked parameters
                [linkFit, linkParam] = getIndicesOfLinkedParams(fitData, pn, comp, sec, numJoins, fitData{i}.links{j});

                if length(linkFit) ~= length(fitData{i}.links{j})
                    warning('Could not find all linked parameters of fit %d parameter %d (%s) in the other fits', i, j, pn);
                end
                
                % Assign ID to linked parameters
                for lk = 1:length(linkFit)
                    localToGlobal(linkFit(lk), linkParam(lk)) = localToGlobal(i, j);
                end
            end
        end
    end
    
    % Second pass: Build global to local map
    globalToLocalFit = zeros(ctrId - 1, 1);    % Maps idxGlobalParam -> idxFit
    globalToLocalParam = zeros(ctrId - 1, 1);  % Maps idxGlobalParam -> idxLocalParam
    globalIdOrder = [];
    for i = 1:length(fitData)
        
        % Map global to local indices
        for j = 1:length(fitData{i}.sim.sensitivities) - numJoins(i)
            if localToGlobal(i, j) <= 0
                error('Internal error: Detected unmapped parameter.');
            end
            
            if ~any(globalIdOrder == localToGlobal(i, j))
                % Add Id to mapped list
                globalIdOrder = [globalIdOrder, localToGlobal(i, j)];

                % Add mapping info
                globalToLocalFit(localToGlobal(i, j)) = i;
                globalToLocalParam(localToGlobal(i, j)) = j;
            end
        end
    end
    
    % Collect weights
    globalWeights = zeros(length(fitData), 1);
    localWeights = cell(length(fitData), 1);
    signalWeights = cell(length(fitData), 1);
    hasWeightError = false;
    for i = 1:length(fitData)
        curFit = fitData{i};

        if ~isfield(curFit, 'weight')
            globalWeights(i) = 1;
        else
            if (numel(curFit.weight) ~= 1) || any(curFit.weight < 0)
                disp(['Error in fit ' num2str(i) ' weights: Global weight is negative which is not allowed']);
                hasWeightError = true;
            else
                globalWeights(i) = curFit.weight;
            end
        end
        
        if ~isfield(curFit, 'weightComponent')
            localWeights{i} = ones(1, size(curFit.outMeas, 2));
        else
            localWeights{i} = curFit.weightComponent;
            if numel(curFit.weightComponent) ~= size(curFit.outMeas, 2)
                disp(['Error in fit ' num2str(i) ' weights: Number of weights in weightComponent vector (' num2str(numel(curFit.weightComponent)) ') does not match number of observed components (' num2str(size(curFit.outMeas, 2)) ')']);
                hasWeightError = true;
            end
            if any(curFit.weightComponent < 0)
                disp(['Error in fit ' num2str(i) ' weights: Some weights in weightComponent are negative which is not allowed']);
                hasWeightError = true;
            end
        end
        
        if ~isfield(curFit, 'weightSignal')
            signalWeights{i} = ones(size(curFit.outMeas));
        else
            signalWeights{i} = curFit.weightSignal;
            if any(size(curFit.weightSignal) ~= size(curFit.outMeas))
                disp(['Error in fit ' num2str(i) ' weights: Format of weights in weightSignal matrix (' ...
                    num2str(size(curFit.weightSignal,1)) 'x' num2str(size(curFit.weightSignal,2)) ') does not match format of data points (' ...
                    num2str(size(curFit.outMeas, 1)) 'x' num2str(size(curFit.outMeas, 2)) ')']);
                hasWeightError = true;
            end
            if any(curFit.weightSignal(:) < 0)
                disp(['Error in fit ' num2str(i) ' weights: Some weights in weightSignal are negative which is not allowed']);
                hasWeightError = true;
            end
        end
    end
    
    if any(~cellfun(@(x) isfield(x, 'weight'), fitData)) && any(cellfun(@(x) isfield(x, 'weight'), fitData))
        hasWeightError = true;
        disp('Error in fit setup: Experiments must either all have weights or none');
    end
    
    % Check parameters and bounds for errors
    [hasError] = checkParamSize(length(globalToLocalFit), length(initParams), length(loBound), length(upBound), length(logScale));
    if hasError || hasWeightError
        disp('Aborting due to errors.');
        outParams = initParams;
        res = -1;
        return;
    end

    % Transform to optimizer space
    if customTransform && ~parameterTransform.postLogTransform
        initParams = parameterTransform.transform(initParams);

        if parameterTransform.transformBounds
            if ~isempty(upBound)
                upBound = parameterTransform.transform(upBound);
            end
            if ~isempty(loBound)
                loBound = parameterTransform.transform(loBound);
            end
        end
    end

    % Convert parameters to log scale
    globalLogScale = zeros(1, length(globalToLocalFit));
    for i = 1:length(globalToLocalFit)
        idxFit = globalToLocalFit(i);
        if logScale(i)

            if initParams(i) == 0
                error('Error: Initial parameter %d must not be zero if using log scale for parameter %d in fit %d!', i, globalToLocalParam(i), idxFit);
            end

            globalLogScale(i) = 1 * (initParams(i) > 0) + (-1) * (initParams(i) < 0);
            initParams(i) = log(globalLogScale(i) .* initParams(i));

            % Convert bounds to log scale
            % Clip bounds to (0, inf) or (-inf, 0)
            if globalLogScale(i) == 1
                if ~isempty(loBound)
                    loBound(i) = max(loBound(i), 0);
                end
                if ~isempty(upBound)
                    upBound(i) = min(upBound(i), inf);
                end
            else
                if ~isempty(loBound)
                    tmp = loBound(i);
                    loBound(i) = max(-upBound(i), 0);
                end
                if ~isempty(upBound)
                    upBound(i) = min(-tmp, inf);
                end
            end

            if ~isempty(loBound)
                loBound(i) = log(loBound(i));
            end
            if ~isempty(upBound)
                upBound(i) = log(upBound(i));
            end
        end
    end
    
    % Transform to optimizer space
    if customTransform && parameterTransform.postLogTransform
        initParams = parameterTransform.transform(initParams);

        if parameterTransform.transformBounds
            if ~isempty(upBound)
                upBound = parameterTransform.transform(upBound);
            end
            if ~isempty(loBound)
                loBound = parameterTransform.transform(loBound);
            end
        end
    end

    % Fit
    opts = optimset('TolFun', 1e-8, 'TolX', 1e-8, 'MaxIter', 100, 'Diagnostics', 'off');
    options_monitor = optimset('Display', 'iter');
    options_lsq = optimset('Jacobian','on');
    
    opts =  optimset(opts, options_lsq);
    if ~quietMode
        opts = optimset(opts, options_monitor);
    end
    if ~isempty(optOptions)
        % Apply custom optimizer settings
        opts = optimset(opts, optOptions);
    end
    
    lastParamsTried = [];

    % Call the optimizer
    try
        [cadetparams, res, ~, exitflag] = lsqnonlin(@residual, initParams, loBound, upBound, opts);
        success = (exitflag > 0);
    catch exc
        disp('ERROR in lsqnonlin. Probably due to failed CADET simulation.');
        disp(exc.message);
        cadetparams = lastParamsTried;
        res = -1;
        success = false;
    end

    % Transform back to simulator space
    if customTransform && parameterTransform.postLogTransform
        cadetparams = parameterTransform.invTransform(cadetparams);
    end

    % Convert back to normal scale
    idxLog = globalLogScale ~= 0;
    outParams = cadetparams;
    outParams(idxLog) = globalLogScale(idxLog) .* exp(cadetparams(idxLog));

    % Transform back to simulator space
    if customTransform && ~parameterTransform.postLogTransform
        outParams = parameterTransform.invTransform(outParams);
    end

    function [globRes, globJac] = residual(cadetparams)
    %RESIDUAL Calculates the residual of the optimization problem
        
        % Save the current parameters 
        lastParamsTried = cadetparams;
        lastParamsTriedExp = cadetparams;

        nf = length(fitData);

        % Transform back to simulator space
        if customTransform && parameterTransform.postLogTransform
            cadetparams = parameterTransform.invTransform(cadetparams);
        end

        % Convert from log scale to normal scale
        idxLog = globalLogScale ~= 0;
        cadetparams(idxLog) = globalLogScale(idxLog) .* exp(cadetparams(idxLog));
        lastParamsTriedExp(idxLog) = globalLogScale(idxLog) .* exp(lastParamsTriedExp(idxLog));

        % Transform back to simulator space
        if customTransform && ~parameterTransform.postLogTransform
            cadetparams = parameterTransform.invTransform(cadetparams);
        end

        globRes = [];
        globJac = [];

        for k = 1:nf

            % Get local params from global params
            idxLocal = localToGlobal(k, :);
            idxLocal = idxLocal(idxLocal > 0);
            localParams = cadetparams(idxLocal);

            [~, idxGlobal] = sort(idxLocal);
            
            curFit = fitData{k};

            % Simulate with current parameters
            if nargout > 1
                [out, jac] = forwardSim(curFit.tOut, curFit.sim, curFit.task, curFit.joins, localParams);
                
                if size(curFit.outMeas, 2) > 1
                    % Concatenate Jacobians of selected components
                    jacSim = squeeze(mat2cell(jac(:,:,curFit.idxComp), size(jac, 1), size(jac, 2), ones(length(curFit.idxComp), 1)));
                    jacSim = vertcat(jacSim{:});
                else
                    % Sum simulated Jacobians since only sum is observed
                    jacSim = sum(jac(:,:,curFit.idxComp), 3);
                end
            else
                out = forwardSim(curFit.tOut, curFit.sim, curFit.task, curFit.joins, localParams);
            end

            weightSig = signalWeights{k};
            if size(curFit.outMeas, 2) > 1
                % Select observed components
                sim = out(:, curFit.idxComp);

                % Calculate local residual
                r = (sim - curFit.outMeas);
                locWeights = repmat(localWeights{k}, size(r,1), 1);
                r = r .* weightSig .* locWeights;
                r = globalWeights(k) * r(:);
                
                locWeights = locWeights(:);
            else
                % Sum simulated signals since only sum is observed
                sim = sum(out(:,curFit.idxComp), 2);

                % Calculate local residual
                r = globalWeights(k) .* weightSig .* (sim - curFit.outMeas);
                
                locWeights = ones(length(r), 1);
            end
            
            if nargout > 1
                % Copy local Jacobian to reordered globalized Jacobian
                jacSimGlobalized = zeros(size(jacSim,1), length(cadetparams));
                jacSimGlobalized(:, idxLocal) = jacSim;

                if customTransform && parameterTransform.postLogTransform
                    jacSimGlobalized = parameterTransform.chainRuleInvTransform(jacSimGlobalized, lastParamsTried);
                end

                % Adapt Jacobian to log scale by chain rule
                idxMask = false(size(idxLog));
                idxMask(idxLocal) = true;
                if any(idxLog & idxMask)
                    jacSimGlobalized(:, idxLog & idxMask) = jacSim(:, idxLog(idxLocal)) .* repmat(localParams(idxLog(idxLocal)), size(jacSim,1), 1);
                end

                if customTransform && ~parameterTransform.postLogTransform
                    jacSimGlobalized = parameterTransform.chainRuleInvTransform(jacSimGlobalized, lastParamsTriedExp);
                end

                % Apply weights
                jacSimGlobalized(:, idxLocal) = jacSimGlobalized(:, idxLocal) .* repmat(weightSig(:), 1, length(idxLocal)) .* repmat(locWeights, 1, length(idxLocal)) .* globalWeights(k);
                
                % Concatenate local jacobians to obtain global Jacobian
                globJac = [globJac; jacSimGlobalized];
            end

            % Concatenate local residuals to obtain global residual
            globRes = [globRes; r];

            % Plot
            if enablePlot
                
                numObserved = sum(curFit.idxComp);
                if ~all(islogical(curFit.idxComp))
                    numObserved = length(curFit.idxComp);
                end
                
                if (numObserved >= 2) && (size(curFit.outMeas, 2) <= 1)
                    % Plot simulated signals on the left
                    subplot(nf, 2, 2*k-1);

                    plot(curFit.tOut, out(:,curFit.idxComp));
                    grid on;

                    legNames = cell(numObserved, 1);
                    for num = 1:numObserved
                        legNames{num} = ['Comp ' num2str(num)];
                    end

                    handle = legend(legNames);
                    set(handle, 'Box', 'off')
                    set(handle, 'Location','NorthWest');

                    % Plot sum signal and fit on the right
                    subplot(nf, 2, 2*k);
                else
                    % Plot fits in a column from top to bottom
                    subplot(nf, 1, k);
                end
                
                if enableWeightedPlot && (size(curFit.outMeas, 2) > 1)
                    hdMeas = plot(curFit.tOut, repmat(localWeights{k}, size(curFit.outMeas, 1), 1) .* curFit.outMeas);
                    hold on;
                    hdSim = plot(curFit.tOut, repmat(localWeights{k}, size(sim, 1), 1) .* sim);
                else
                    hdMeas = plot(curFit.tOut, curFit.outMeas);
                    hold on;
                    hdSim = plot(curFit.tOut, sim);
                end
                hold off;
                grid on;
                
                if size(curFit.outMeas, 2) <= 1                
                    handle = legend('Meas', 'Sim Sum');
                    set(hdSim(1), 'Color', 'r');
                else
                    legNames = cell(size(curFit.outMeas, 2) * 2, 1);
                    for num = 1:size(curFit.outMeas, 2)
                        legNames{num} = ['Meas Comp ' num2str(num)];
                        legNames{num+size(curFit.outMeas, 2)} = ['Sim Comp ' num2str(num)];
                        
                        % Adapt line styles
                        set(hdSim(num), 'LineStyle', '--', 'LineWidth', 4.0);
                    end
                    
                    handle = legend(legNames);
                end
                
                set(handle, 'Box', 'off')
                set(handle, 'Location','NorthWest');
                
                str = sprintf('%g | ', localParams);
                title(str(1:end-3));
                drawnow;
            end
            
        end
    end
    
end

function [ sol, jac ] = forwardSim(tOut, sim, task, joins, localParams)
%FORWARDSIM Solve General Rate Model using CADET
% Parameters:
%   - tOut: Vector of outlet time points
%   - sim: Instance of the Simulator ready for simulation
%   - task: Executable task for the simulator
%   - joins: Cell array with vectors of joined parameters
%   - localParams: Current parameters given to CADET
%
% Returns:
%   - sol: Matrix with concentration profiles at the system outlet.
%       Each column gives the profile of a component.
%   - jac: Jacobian of the outlet concentrations with respect to
%       cadetParams. jac is a 3D matrix / rank 3 tensor in which the
%       first dimension corresponds to tOut, second to the parameter,
%       and the third to the component.

    if ~isempty(joins)
        % Assign value to joined parameters
        nTrueParams = length(localParams);
        localParams = [localParams, zeros(1, length(sim.sensitivities) - length(localParams))];
        for i = 1:length(joins)
            localParams(joins{i}) = localParams(i);
        end
    end
    
    % Simulate
    result = sim.runWithParameters(task, localParams);
    sol = result.solution.outlet;
    jac = result.sensitivity.jacobian;
    
    if ~isempty(joins)
        % Apply chain rule by adding joined parameters
        jacNew = jac(:, 1:nTrueParams, :);
        for i = 1:length(joins)
            jacNew(:, i, :) = jacNew(:, i, :) + sum(jac(:, joins{i}, :), 2);
        end
        jac = jacNew;
    end
end

function [linkFit, linkParam] = getIndicesOfLinkedParams(fitData, pn, comp, sec, numJoins, idxSearch)
%GETINDICESOFLINKEDPARAMS Finds the fit and parameter index of matching linked parameters
    linkFit = [];
    linkParam = [];
    for i = idxSearch
        pi = fitData{i}.sim.sensitivities;
        for j = 1:length(pi) - numJoins(i)
            % Check if parameter matches
            if strcmpi(pn, pi{j}.SENS_NAME) && (comp == pi{j}.SENS_COMP) && (sec == pi{j}.SENS_SECTION)
                % Add indices to list
                linkFit = [linkFit, i];
                linkParam = [linkParam, j];
            end
        end
    end
end

function [hasError, fitData, numJoins] = checkInput(fitData)
%CHECKINPUT Checks the fitData array for errors

    numJoins = zeros(length(fitData), 1);

    hasError = ((size(fitData,1) > 1) && (size(fitData,2) > 1));
    if hasError
        disp(['Error in fitData: No matrix supported, fitData has to be a vector of cells']);
    end
    for i = 1:length(fitData)

        % Check for necessary fields
        hasError = checkForField(fitData{i}, 'tOut', '', 'Field "tOut" is missing from struct', i) | hasError;
        hasError = checkForField(fitData{i}, 'outMeas', '', 'Field "outMeas" is missing from struct', i) | hasError;
        hasError = checkForField(fitData{i}, 'sim', '', 'Simulator is missing from struct', i) | hasError;
        hasError = checkForField(fitData{i}, 'idxComp', '', 'Indices of observed components is missing from struct', i) | hasError;

        % Check and create, if necessary, the task field
        if (~isfield(fitData{i}, 'task') || isempty(fitData{i}.task)) && isfield(fitData{i}, 'sim') && ~isempty(fitData{i}.sim)
            fitData{i}.task = fitData{i}.sim.prepareSimulation();
        else
            if ~isfield(fitData{i}, 'task') || isempty(fitData{i}.task)
                hasError = true;
                disp(['Error in fit ' num2str(i) ': Task field does not exist or is empty']);
            end
        end

        if ~isfield(fitData{i}, 'logScale')
            % Enable logScale by default
            fitData{i}.logScale = true;
        end
        
        if ~isfield(fitData{i}, 'joins')
            % No joins by default
            fitData{i}.joins = [];
        else
            
            % Total number of joined parameters in this fit
            numJoins(i) = sum(arrayfun(@(x) length(x{1}), fitData{i}.joins));
            
            % Check the joins
            for j = 1:length(fitData{i}.joins)
                curJoin = fitData{i}.joins{j};

                % Check if a join is in an invalid index range (joins are
                % the last parameters)
                if any((curJoin <= length(fitData{i}.sim.sensitivities) - numJoins(i)) | (curJoin > length(fitData{i}.sim.sensitivities)))
                    disp(['Error in fit ' num2str(i) ': Parameter indices of join ' num2str(j) ' are out of range']);
                    hasError = true;
                end
                
                % Check if joined parameters appear in other master
                % parameters
                for k = j+1:length(fitData{i}.joins)
                    if length(unique([curJoin(:); fitData{i}.joins{k}(:)])) ~= length(curJoin) + length(fitData{i}.joins{k})
                        disp(['Error in fit ' num2str(i) ': Joined parameters of parameter ' num2str(j) ' may not be joined to parameter' num2str(k)]);
                        hasError = true;
                    end
                end
            end
            
            % Check number of joins
            if length(fitData{i}.sim.sensitivities) - numJoins(i) ~= length(fitData{i}.joins)
                disp(['Error in fit ' num2str(i) ': Number of joins has to match number of parameters']);
                hasError = true;
            end
        end
        
        if ~isfield(fitData{i}, 'links')
            % No links by default
            fitData{i}.links = repmat({[]}, length(fitData{i}.sim.sensitivities) - numJoins(i), 1);
        else

            % Check number of links
            if length(fitData{i}.links) ~= length(fitData{i}.sim.sensitivities) - numJoins(i)
                disp(['Error in fit ' num2str(i) ': Number of links has to match number of parameters']);
                hasError = true;
            end

            for j = 1:length(fitData{i}.links)
                % Check if links have the correct range
                if any((fitData{i}.links{j} <= 0) | (fitData{i}.links{j} > length(fitData)))
                    disp(['Error in fit ' num2str(i) ': Link indices of parameter ' num2str(j) ' are out of range']);
                    hasError = true;
                end

                % Check if links include the current fit
                if any(fitData{i}.links{j} == i)
                    disp(['Error in fit ' num2str(i) ': Cannot link parameter ' num2str(j) ' to itself']);
                    hasError = true;
                end
            end
            
        end
        
        % Check for consistency
        % Time points should be strictly increasing
        hasError = checkTimePoints(fitData{i}.tOut, 'measurements', i) | hasError;

        if length(fitData{i}.tOut) ~= size(fitData{i}.outMeas, 1)
            disp(['Error in fit ' num2str(i) ' measurements: Number of time points and concentrations does not match']);
            hasError = true;
        end
        
        if size(fitData{i}.outMeas, 2) > 1
            if isfield('idxComp', fitData{i}) && (length(fitData{i}.idxComp) ~= size(fitData{i}.outMeas, 2))
                disp(['Error in fit ' num2str(i) ' measurements: Number of signals and number of indices of observed components do not match']);
                hasError = true;
            end
        end
    end
end

function hasError = checkParamSize(numParams, numInitParams, numLower, numUpper, numLogScale)
%CHECKPARAMSIZE Checks whether the initial params and bounds have the correct number of entries

    hasError = false;

    if numInitParams ~= numParams
        disp(['Error: Fit setup has ' num2str(numParams) ' parameters but ' num2str(numInitParams) ' initial parameters were given']);
        hasError = true;
    end
    
    if numLogScale ~= numParams
        disp(['Error: Fit setup has ' num2str(numParams) ' parameters but ' num2str(numLogScale) ' logScale settings were given']);
        hasError = true;
    end
    
    if (numLower > 0) && (numLower ~= numParams)
        disp(['Error: Fit setup has ' num2str(numParams) ' parameters but ' num2str(numLower) ' lower bounds were given']);
        hasError = true;
    end
    
    if (numUpper > 0) && (numUpper ~= numParams)
        disp(['Error: Fit setup has ' num2str(numParams) ' parameters but ' num2str(numUpper) ' upper bounds were given']);
        hasError = true;
    end
end

function isErroneous = checkForField(sct, fieldName, preText, text, noFit)
%CHECKFORFIELD Checks whether a given struct contains a field
    isErroneous = ~isfield(sct, fieldName);
    if isErroneous
        if ~isempty(preText) && (preText ~= '')
            preText = [' ' preText];
        end
        disp(['Error in fit ' num2str(noFit) preText ':' text]);
    end
end

function isErroneous = checkTimePoints(time, name, noFit)
%CHECKTIMEPOINTS Checks whether time points are monotonically increasing
    isErroneous = false;
    if any(time(2:end) - time(1:end-1) <= 0)
        idx = find(time(2:end) - time(1:end-1) <= 0);
        disp(['Error in fit ' num2str(noFit) ': Times of ' name ' are not strictly increasing. Look at the following indices (1-based):']);
        disp(idx+1);
        isErroneous = true;
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
