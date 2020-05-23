
classdef ParameterFit < handle
	%ParameterFit Describes a parameter fitting optimization problem holding models and data
	%
	%   A fit can contain multiple experiments (e.g., different loadings or gradients) which in 
	%   turn may consist of multiple wavelengths (i.e., linear combinbations of the components) 
	%   measured from different unit operations.
	%
	%   Measurements need not begin, nor end with the simulation. Each wavelength can provide
	%   measurements for only a slice of the global simulation time. This slice is indicated by
	%   a mask (i.e., index or logical vector) on the vector of simulation times. However, the
	%   solution has to be evaluated at the measurement time points, that is, the simulation
	%   time points have to include the measurement time points.
	%
	%   Weights can be applied on several levels. All residuals can be weighted inside one
	%   wavelength and the wavelengths can be weighted among each other. Finally, the experiments
	%   can be weighted.
	%
	%   Each experiment has certain parameters that are configured in the Simulator instance that
	%   is associated with the experiment. Parameters between different experiments can be linked
	%   together, such that the same value is assigned to those parameters. This is especially
	%   useful if parameters are shared during different experiments (e.g., when estimating binding
	%   parameters from multiple gradient experiments).
	%
	%   For evaluation of the residual, the order of the parameters is given by the order of the
	%   experiment and the ordering inside that experiment (i.e., experiment-parameter-major). For
	%   example, if experiment E1 has parameters E1P1, E1P2, and experiment E2 has parameters E2P1,
	%   E2P2, then the order of the parameters is E1P1, E1P2, E2P1, E2P2. Linked parameters appear
	%   only once as the first linked parameter in the first linked experiment. If, in the previous
	%   example, the parameters E1P1 and E2P2 were linked, then the ordering would be E1P1, E1P2,
	%   E2P1. Since E2P2 is linked to E1P1 which appears first, it is omitted.
	%
	%   The class offers the possibility of applying parameter transformations for aiding the
	%   optimizer. Multiple transformations are applied in the given order: Given the
	%   transformations [f, g] and parameter vector x, we work with
	%     y = g( f(x) ),
	%   which means that the simulator is called with
	%     x = invF( invG(y) ),
	%   where invF and invG denote the inverse transformations of f and g, respectively.
	%   This also explains why for the Jacobian the chain rule is invoked with the inverse
	%   parameter transformations.

	% Copyright: (C) 2008-2020 The CADET Authors
	%            See the license note at the end of the file.

	properties (Access = 'protected')
		models; % Models to be fitted
		data; % Data
		idxUnit; % Index of the unit operations to fit against
		idxComp; % Index of the components to fit against
		timeMask; % Index of the simulation time points that are compared to the data
		weightsExp; % Weight for each experiment
		weightsWave; % Weight for components
		weightsData; % Weight for data points
		links; % Parameter links between experiments
		localToGlobal; % Maps (idxExp, idxLocalParam) -> idxGlobalParam (uses 1 based indices)
		globalToLocal; % Maps idxGlobalParam -> (idxExp, idxLocalParam) (uses 1 based indices)
		paramOffset; % List with first index of parameter set of each experiment in the global parameter vector
		totalDataPoints; % Total number of data points
		names; % Names of the experiments
		waveNames; % Names of the wavelengths within the experiments
		fdFactors; % Finite differences step sizes
		lastResult; % Last simulation results
		lastParams; % Last parameters simulated successfully
	end

	properties (Dependent)
		nParameters; % Number of (global) parameters
		nExperiments; % Number of experiments
		lastSimulationResults; % Results of the last simulation performed
		lastParameters; % Parameter values applied for the last successful simulation
		simulators; % Simulators of the different experiments
	end

	properties
		parameterTransform; % Parameter transformation
	end

	methods

		function obj = ParameterFit()
			%PARAMETERFIT Creates an object of the ParameterFit class
			obj.models = [];
			obj.data = [];
			obj.parameterTransform = [];
			obj.idxUnit = [];
			obj.idxComp = [];
			obj.timeMask = [];
			obj.weightsExp = [];
			obj.weightsWave = [];
			obj.weightsData = [];
			obj.links = [];
			obj.localToGlobal = [];
			obj.globalToLocal = [];
			obj.paramOffset = [];
			obj.totalDataPoints = 0;
			obj.names = [];
			obj.waveNames = [];
			obj.fdFactors = [];
			obj.lastResult = [];
			obj.lastParams = [];
		end

		function addExperiment(obj, data, sim, idxUnit, idxComp, timeMask, weightExp, weightWave, weightData, name, waveNames, fdFactors)
			%ADDEXPERIMENT Adds an experiment to the fit
			%   ADDEXPERIMENT(DATA, SIM, IDXUNIT, IDXCOMP) adds a new experiment to the fit. An experiment can
			%   comprise of multiple wavelengths (i.e., linear combinations of components). 
			%   DATA is a cell array that holds the measurements for each wavelength as vector.
			%   SIM is a Simulator instance that is ready to run (i.e., it is configured) the model.
			%   IDXUNIT is a vector which contains the 0-based index of the unit operation a wavelength
			%   measurements corresponds to.
			%   IDXCOMP is a cell array of vectors that denote the coefficients of a linear combination of
			%   components (e.g., extinction coefficients).
			%
			%   ADDEXPERIMENT(..., TIMEMASK) employs a maks for each wavelength that can include or exclude
			%   simulation time points. TIMEMASK can be given as logical array (same size as simulation time
			%   points) or vector of indices (same size as data points). By default, all simulation time
			%   points between the beginning of the simulation and the last measurement are used.
			%
			%   ADDEXPERIMENT(..., TIMEMASK, WEIGHTEXP) additionally assigns a weight to this experiment,
			%   which is 1.0 by default (all experiments are equally important). WEIGHTEXP is expected to
			%   be a positive scalar.
			%
			%   ADDEXPERIMENT(..., TIMEMASK, WEIGHTEXP, WEIGHTWAVE) additionally assigns a weight to each wavelength
			%   of this experiment. By default, WEIGHTWAVE is set to a vector of ones (all wavelengths are
			%   equally important). WEIGHTWAVE is expected to be a vector with positive elements.
			%
			%   ADDEXPERIMENT(..., TIMEMASK, WEIGHTEXP, WEIGHTWAVE, WEIGHTDATA) additionally assigns a weight
			%   to each data point of this experiment. WEIGHTDATA is a cell array of vectors with non-negative
			%   elements. By default, all data points are assumed equally important, that is, they are all
			%   assigned a weight of 1.0.
			%
			%   ADDEXPERIMENT(..., TIMEMASK, WEIGHTEXP, WEIGHTWAVE, WEIGHTDATA, NAME) specifies a name for
			%   this experiment. This is only used for placing a caption on a plot. Defaults to a subsequent
			%   numbering in the style of 'Experiment 1', 'Experiment 2', etc. NAME is expected to be a string.
			%
			%   ADDEXPERIMENT(..., TIMEMASK, WEIGHTEXP, WEIGHTWAVE, WEIGHTDATA, NAME, WAVENAMES) specifies a name
			%   for each wavelength in this experiment. This is only used for legends in a plot. Defaults to a
			%   subsequent numbering in the style of 'Wavelength 1', 'Wavelength 2', etc. WAVENAMES is expected
			%   to be a cell array of strings.
			%
			%   ADDEXPERIMENT(..., TIMEMASK, WEIGHTEXP, WEIGHTWAVE, WEIGHTDATA, NAME, WAVENAMES, FDFACTORS)
			%   determines the relative step size for each parameter when using finite differences. The use of
			%   finite differences requires FDFACTORS to be all positive (if FDFACTORS is all negative, finite
			%   differences are disabled) and parameter sensitivities of SIM to be disabled. Instead, variable
			%   parameters are required. If parameter sensitivities are enabled, finite differences will not
			%   be used regardless of FDFACTORS. Defaults to 1e-6 for each parameter.

			% Check given data and model
			validateattributes(sim, Simulator.availableInterfaceClasses, {'nonempty', 'numel', 1}, '', 'sim');
			if ~sim.validate()
				error('CADET:funcParamError', 'Expected sim to be fully configured and valid.');
			end

			model = sim.model;
			numSimTime = numel(sim.solutionTimes);

			validateattributes(idxUnit, {'numeric'}, {'nonempty', 'finite', 'real', 'vector', '>=', 0, '<', model.numUnitOperations}, '', 'idxUnit');
			validateattributes(idxComp, {'cell'}, {'nonempty', 'vector', 'numel', numel(idxUnit)}, '', 'idxComp');
			validateattributes(data, {'cell'}, {'nonempty', 'vector', 'numel', length(idxUnit)}, '', 'data');

			if isa(model, 'SingleUnitOpSystem')
				arrayfun(@(x) validateattributes(idxComp{x}, {'numeric'}, {'nonempty', 'finite', 'real', 'vector', 'numel', model.nComponents}, '', sprintf('idxComp{%d}', x)), 1:length(idxUnit));
			else
				arrayfun(@(x) validateattributes(idxComp{x}, {'numeric'}, {'nonempty', 'finite', 'real', 'vector', 'numel', model.models(idxUnit(x) + 1).nComponents}, '', sprintf('idxComp{%d}', x)), 1:length(idxUnit));
			end
			arrayfun(@(x) validateattributes(data{x}, {'numeric'}, {'nonempty', 'finite', 'real', 'vector'}, '', sprintf('data{%d}', x)), 1:length(idxUnit));
			idxTooManyData = find(arrayfun(@(x) numel(data{x}) > numSimTime, 1:length(idxUnit)));
			if ~isempty(idxTooManyData)
				temp = sprintf('data{%d}, ', idxTooManyData);
				error('CADET:funcParamError', 'Expected %s to contain at most %d (length of sim.solutionTimes) points.', temp(1:end-2), numSimTime);
			end

			if (nargin <= 5) || isempty(timeMask)
				timeMask = cellfun(@(x) 1:size(x,1), data, 'UniformOutput', false);
			end
			validateattributes(timeMask, {'cell'}, {'nonempty', 'vector', 'numel', numel(idxUnit)}, '', 'timeMask');
			for i = 1:length(timeMask)
				if islogical(timeMask{i})
					if sum(timeMask{i}) ~= numel(data{i})
						error('CADET:funcParamError', 'Expected timeMask{%d} to select %d elements, but got %d.', i, numel(data{i}), sum(timeMask{i}));
					end
					if numel(timeMask{i}) ~= numSimTime
						error('CADET:funcParamError', 'Expected timeMask{%d} to have length %d, but got %d.', i, numSimTime, numel(timeMask{i}));
					end
				else
					if any((timeMask{i} < 1) | (timeMask{i} > numSimTime)) || (numel(timeMask{i}) ~= numel(data{i}))
						error('CADET:funcParamError', 'Expected timeMask{%d} to contain %d valid indices into the simulation results.', i, numel(data{i}));
					end
				end
			end

			arrayfun(@(x) validateattributes(timeMask{x}, {'numeric', 'logical'}, {'nonempty', 'finite', 'real', 'vector'}, '', sprintf('timeMask{%d}', x)), 1:length(idxUnit));

			if (nargin <= 6) || isempty(weightExp)
				weightExp = 1.0;
			end
			validateattributes(weightExp, {'numeric'}, {'nonempty', 'finite', 'real', 'scalar', '>', 0}, '', 'weightExp');

			if (nargin <= 7) || isempty(weightWave)
				weightWave = ones(length(idxUnit), 1);
			end
			validateattributes(weightWave, {'numeric'}, {'nonempty', 'finite', 'real', 'vector', '>', 0, 'numel', length(idxUnit)}, '', 'weightWave');

			if (nargin <= 8) || isempty(weightData)
				weightData = arrayfun(@(x) ones(size(data{x})), 1:length(idxUnit), 'UniformOutput', false);
			end
			validateattributes(weightData, {'cell'}, {'nonempty', 'vector', 'numel', length(idxUnit)}, '', 'weightData');
			arrayfun(@(x) validateattributes(weightData{x}, {'numeric'}, {'nonempty', '>=', 0, 'finite', 'real', 'vector', 'numel', numel(data{x})}, '', sprintf('weightData{%d}', x)), 1:length(idxUnit));

			if (nargin <= 9) || isempty(name)
				name = sprintf('Experiment %d', length(obj.models) + 1);
			end
			validateattributes(name, {'char'}, {'nonempty'}, '', 'name');

			if (nargin <= 10) || isempty(waveNames)
				waveNames = arrayfun(@(x) sprintf('Wavelength %d', x), 1:length(idxUnit), 'UniformOutput', false);
			end
			validateattributes(waveNames, {'cell'}, {'nonempty', 'numel', length(idxUnit)}, '', 'waveNames');
			arrayfun(@(x) validateattributes(waveNames{x}, {'char'}, {'nonempty'}, '', sprintf('waveNames{%d}', x)), 1:length(idxUnit));

			if (nargin <= 11) || isempty(fdFactors)
				fdFactors = 1e-6 .* ones(sim.nVariableParameters, 1);
			end
			validateattributes(fdFactors, {'numeric'}, {'nonempty', 'vector', 'numel', sim.nVariableParameters}, '', 'fdFactors');
			if ~(all(fdFactors > 0.0) || all(fdFactors <= 0.0))
				error('CADET:funcParamError', 'Expected fdFactors to contain only positive or only negative values.');
			end

			% Convert to standard format
			data = cellfun(@(x) x(:), data, 'UniformOutput', false);
			idxComp = cellfun(@(x) x(:), idxComp, 'UniformOutput', false);
			timeMask = cellfun(@(x) x(:), timeMask, 'UniformOutput', false);
			weightData = cellfun(@(x) x(:), weightData, 'UniformOutput', false);

			% Append to internal lists
			obj.data{end+1} = data;
			obj.models{end+1} = sim;
			obj.idxUnit{end+1} = idxUnit;
			obj.idxComp{end+1} = idxComp;
			obj.timeMask{end+1} = timeMask;
			obj.weightsExp(end+1) = weightExp;
			obj.weightsWave{end+1} = weightWave(:).';
			obj.weightsData{end+1} = weightData;
			obj.names{end+1} = name;
			obj.waveNames{end+1} = waveNames;
			obj.fdFactors{end+1} = fdFactors;

			obj.rebuildInternalDataStructures();
		end

		function removeExperiment(obj, idx)
			%REMOVEEXPERIMENT Removes an experiment from the fit
			%   REMOVEEXPERIEMNT(IDX) Removes the experiment at index IDX from the fit. Also removes all links
			%   that point to this experiment or originate from this experiment.

			validateattributes(idx, {'numeric'}, {'nonempty', 'vector', '>=', 1, '<=', length(obj.data)}, '', 'idx');

			% Remove links to experiments
			for i = 1:length(obj.links)
				curLink = obj.links{i};
				for j = 1:length(idx)
					idxLink = curLink(:, 1) == idx(j);
					curLink(idxLink, :) = [];
				end
				obj.links{i} = curLink;
			end

			obj.data{idx} = [];
			obj.models{idx} = [];
			obj.idxUnit{idx} = [];
			obj.idxComp{idx} = [];
			obj.timeMask{idx} = [];
			obj.weightsExp(idx) = [];
			obj.weightsWave{idx} = [];
			obj.weightsData{idx} = [];
			obj.names{idx} = [];
			obj.waveNames{idx} = [];
			obj.fdFactors{idx} = [];

			obj.rebuildInternalDataStructures();
		end

		function addLink(obj, idxExp, idxParam)
			%ADDLINK Adds a link between two (or more) parameters of two (or more) distinct experiments
			%   ADDLINK(IDXEXP) treats the parameters identified in IDXEXP as a single parameter.
			%   The parameters are identified by the matrix IDXEXP in which the first column denotes
			%   the 1-based index of the experiment and the second column contains the 1-based index
			%   of the parameter. The matrix has to have at least 2 rows (i.e., two parameters that
			%   are linked).
			%
			%   ADDLINK(IDXEXP, IDXPARAM) treats the parameters identified by the two vectors IDXEXP
			%   and IDXPARAM as a single parameter. Each pair IDXEXP(i), IDXPARAM(i) denotes a parameter
			%   identified by a 1-based index of its experiment and a 1-based index of the parameter
			%   inside that particular experiment. At least two parameters are required (i.e., IDXEXP
			%   and IDXPARAM must have at least two elements each).
			%   
			%   Note that links between parameters inside one single experiments are not possible and
			%   should be implemented by joint parameters.

			% Extract idxParam from idxExp if we are only given idxExp
			if (nargin <= 2) || isempty(idxParam)
				if (mod(numel(idxExp), 2) ~= 0) || (numel(idxExp) < 4) || isempty(idxExp)
					error('CADET:funcParamError', 'Expected idxExp to have a multiple of 2 and at least 4 elements.');
				else
					if (size(idxExp, 2) == 2) && (size(idxExp, 1) > 1)
						idxParam = idxExp(:, 2);
						idxExp = idxExp(:, 1);
					elseif (size(idxExp, 1) == 2) && (size(idxExp, 2) > 1)
						idxParam = idxExp(2, :);
						idxExp = idxExp(1, :);
					end
				end
			end

			% Check input
			validateattributes(idxExp, {'numeric'}, {'nonempty', 'vector', '>=', 1, '<=', length(obj.data)}, '', 'idxExp');
			if length(unique(idxExp)) ~= length(idxExp)
				error('CADET:funcParamError', 'Detected unsupported parameter link within one experiment. Use parameter joins instead.');
			end
			validateattributes(idxParam, {'numeric'}, {'nonempty', 'vector', '>=', 1, 'numel', length(idxExp)}, '', 'idxParam');

			for i = 1:length(idxExp)
				ie = idxExp(i);
				% Check if link table aready contains such a link
				for j = 1:length(idxExp)
					if i == j
						continue;
					end
					curLinks = obj.links{ie};
					if ~isempty(curLinks) && any(all(repmat([idxExp(j), idxParam(j)], size(curLinks, 1), 1) == curLinks(:, 2:end), 2))
						error('CADET:funcParamError', 'Link to parameter %d of experiment %d already exists in experiment %d.', idxParam(j), idxExp(j), ie);
					end
				end
			end

			% Add links to link tables
			for i = 1:length(idxExp)
				ie = idxExp(i);
				table = [ones(numel(idxExp), 1) .* idxParam(i), idxExp(:), idxParam(:)];
				table(i,:) = [];

				obj.links{ie} = [obj.links{ie}; table];
			end

			obj.rebuildInternalDataStructures();
		end

		function set.parameterTransform(obj, paramTransform)
			if ~isempty(paramTransform)
				validateattributes(paramTransform, {'ParameterTransformation'}, {'vector'}, '', 'paramTransform');
			end
			obj.parameterTransform = paramTransform;
		end

		function val = get.nExperiments(obj)
			val = length(obj.data);
		end

		function val = get.nParameters(obj)
			val = size(obj.globalToLocal, 1);
		end

		function val = get.lastSimulationResults(obj)
			val = obj.lastResult;
		end

		function set.lastSimulationResults(obj, val)
			if isempty(val)
				obj.lastResult = [];
			end
		end

		function val = get.lastParameters(obj)
			val = obj.lastParams;
		end

		function set.lastParameters(obj, val)
			if isempty(val)
				obj.lastParams = [];
			end
		end

		function val = get.simulators(obj)
			val = cellfun(@(x) x, obj.models);
		end

		function val = applyTransform(obj, p)
			%APPLYTRANSFORM Applies the forward parameter transformations
			%   VAL = APPLYTRANSFORM(P) applies all forward transformations to the point P.
			%   Given P in simulator space, VAL is returned in optimizer space.

			s = size(p);
			val = p(:);
			if ~isempty(obj.parameterTransform)
				for i = 1:length(obj.parameterTransform)
					val = obj.parameterTransform(i).transform(val(:));
				end
			end
			val = reshape(val, s);
		end

		function val = applyInverseTransform(obj, p)
			%APPLYINVERSETRANSFORM Applies the backward parameter transformations
			%   VAL = APPLYINVERSETRANSFORM(P) applies all backward transformations to the point P.
			%   Given P in optimizer space, VAL is returned in simulator space.

			s = size(p);
			val = p(:);
			if ~isempty(obj.parameterTransform)
				for i = length(obj.parameterTransform):-1:1
					val = obj.parameterTransform(i).inverseTransform(val(:));
				end
			end
			val = reshape(val, s);
		end

		function [res, jac] = residualVector(obj, p)
			%RESIDUALVECTOR Computes the vector of residuals for the given parameters
			%   RES = RESIDUALVECTOR(P) computes the residual vector RES using the given parameters P.
			%
			%   [RES, JAC] = RESIDUALVECTOR(P) computes the residual vector RES using the given parameters P.
			%   Also returns the Jacobian JAC of the residual with respect to the parameters.
			
			validateattributes(p, {'numeric'}, {'finite', 'real', 'vector', 'numel', size(obj.globalToLocal, 1)}, '', 'p');

			if ~isempty(obj.parameterTransform)
				trueP = obj.applyInverseTransform(p(:));
			else
				trueP = p(:);
			end

			res = zeros(obj.totalDataPoints, 1);
			if nargout > 1
				jac = zeros(obj.totalDataPoints, length(p));
			end

			resIdx = 1;
			newLastResult = cell(length(obj.models), 1);

			% Loop over all experiments
			for i = 1:length(obj.models)

				% Extract local parameters
				localParamIdx = obj.localToGlobal(obj.paramOffset(i):obj.paramOffset(i+1)-1, 3);
				localParams = trueP(localParamIdx);

				% Simulate				
				newLastResult{i} = obj.models{i}.runWithParameters(localParams, true);
				resSol = newLastResult{i}.solution.outlet(:, 1, :); % Format is [nTime, nComp]
				resJac = newLastResult{i}.sensitivity.jacobian; % Format is [nTime, nComp, nParam]

				curUnit = obj.idxUnit{i};
				curComp = obj.idxComp{i};
				curWeightWave = obj.weightsWave{i};
				curWeightData = obj.weightsData{i};
				curMask = obj.timeMask{i};
				curData = obj.data{i};
				fdSteps = obj.fdFactors{i};

				% Check if finite differences are required and enabled
				if (nargout > 1) && (isempty(resJac) || any(arrayfun(@(x) isempty(resJac{x}), curUnit + 1))) && all(fdSteps > 0.0)
					resJac = jacobianFiniteDifferences(obj.models{i}, localParams, fdSteps);
				end

				% Loop over all wavelengths
				for j = 1:length(curUnit)
					nPoints = numel(curData{j});
					unitTrace = curUnit(j) + 1;
					res(resIdx:resIdx+nPoints-1) = (obj.weightsExp(i) .* curWeightWave(j)) .* curWeightData{j} .* (resSol{unitTrace}(curMask{j}, :) * curComp{j} - curData{j});

					if (nargout > 1) && ~isempty(resJac{unitTrace})
						curJac = resJac{unitTrace}(:, 1, :, :); % Format is [nTime, nComp, nParam]

						% Loop over parameters
						for k = 1:size(curJac, 3)
							partialDiff = curJac(curMask{j}, :, k) * curComp{j};
							jac(resIdx:resIdx+nPoints-1, localParamIdx(k)) = (obj.weightsExp(i) .* curWeightWave(j)) .* curWeightData{j} .* partialDiff;
						end
					end

					resIdx = resIdx + nPoints;
				end
			end

			% Apply chain rule from parameter transformation
			if (nargout > 1) && ~isempty(obj.parameterTransform)
				curP = obj.parameterTransform(1).transform(trueP);
				for i = 1:length(obj.parameterTransform)
					jac = obj.parameterTransform(i).chainRuleInvTransform(jac, curP(:), trueP(:));
					curP = trueP(:);
					trueP = obj.parameterTransform(i).transform(trueP(:));
				end
			end

			obj.lastResult = newLastResult;
			obj.lastParams = p;
		end

		function [res, grad] = residualSumOfSquares(obj, p)
			%RESIDUALSUMOFSQUARES Computes the sum of squares of the residuals for the given parameters
			%   RES = RESIDUALSUMOFSQUARES(P) computes the residual sum of squares RES using the given parameters P.
			%
			%   [RES, GRAD] = RESIDUALSUMOFSQUARES(P) computes the residual vector RES using the given parameters P.
			%   Also returns the Jacobian JAC of the residual with respect to the parameters.

			if nargout > 1
				[res, jac] = obj.residualVector(p);
				res = sum(res.^2);
				grad = 2 .* jac.' * res;
			else
				res = obj.residualVector(p);
				res = sum(res.^2);
			end
		end

		function hd = plot(obj, hd, numIter, procData, useWeights, paramHist, paramLogPlot)
			%PLOT Either plots the current fit against the data and provides history of the iterands, or plots data only
			%   PLOT() if no simulation (e.g., call to PARAMETERFIT.RESIDUALSUMOFSQUARES or 
			%   PARAMETERFIT.RESIDUALVECTOR) has been performed yet, plots the data registered in the
			%   fit in a single figure using subplots. If a simluation has been performed, plots the
			%   current fit against the data.
			%
			%   PLOT(HD) updates the figures from a previous call to PARAMETERFIT.PLOT with
			%   the current fit.
			%
			%   PLOT(HD, NUMITER, PROCDATA) updates (or creates if HD is empty) figures with
			%   the current fit and the current number of iterations given in NUMITER and meta
			%   data given in PROCDATA. PROCDATA is a vector with two elements and contains the
			%   current residual (e.g., sum of squares) and the current step size (in this order).
			%
			%   PLOT(..., USEWEIGHTS) determines whether data and fit are scaled according to the
			%   weights registered in the ParameterFit object. This represents the way the optimizer
			%   sees the data. This scaling is enabled by default.
			%
			%   PLOT(..., USEWEIGHTS, PARAMHIST) additionally determines whether a history of the
			%   parameters over the course of the iteration is plotted. This is enabled by default.
			%
			%   PLOT(..., USEWEIGHTS, PARAMHIST, PARAMLOGPLOT) determines whether the history plot
			%   of the parameters (PARAMHIST = true) uses a logarithmic scale. Defaults to true.
			%
			%   HD = PLOT(...) The handles of the main plot components are returned in
			%   the struct HD and can be used for updating the figures in subsequent calls to 
			%   PARAMETERFIT.PLOT.
			%
			%   This function is explicitly build for calling it in an OutputFcn of the optimizer.
			%   See any of the parameter estimation examples for an instruction on how to call it.
			%
			%   See also SINGLEEXPSEPARATECOMPONENTS

			if isempty(obj.lastResult)
				% We have no simulations run yet, plot only data

				hd.figure = figure('Name', 'Parameter Fit');
				hd.sim = cell(length(obj.data), 1);
				hd.data = cell(length(obj.data), 1);
				hd.legend = cell(length(obj.data), 1);
				hd.title = cell(length(obj.data), 1);
				for i = 1:length(obj.data)
					subplot(length(obj.data), 1, i);
					
					curMask = obj.timeMask{i};
					curData = obj.data{i};
					curWeightData = obj.weightsData{i};
					hd.data{i} = zeros(length(curData));

					% Loop over wavelengths
					for j = 1:length(curData)
						plotData = curData{j};
						if useWeights
							plotData = obj.weightsExp(i) .* obj.weightsWave{i} .* curWeightData{j} .* curData{j};
						end

						hold on;
						hd.data{i}(j) = plot(obj.models{i}.solutionTimes(curMask{j}), plotData);
						hold off;
					end
					
					% Add legend
					hd.legend{i} = legend(obj.waveNames{i});
					set(hd.legend{i}, 'Box', 'off')
					set(hd.legend{i}, 'Location', 'NorthWest');

					hd.title{i} = title(sprintf('%s', obj.names{i}));
					grid on;
				end

				return;
			end

			if (nargin <= 4) || isempty(useWeights)
				useWeights = true;
			end
			validateattributes(useWeights, {'logical'}, {'scalar', 'nonempty'}, '', 'useWeights');

			if (nargin <= 5) || isempty(paramHist)
				paramHist = true;
			end
			validateattributes(paramHist, {'logical'}, {'scalar', 'nonempty'}, '', 'paramHist');

			if (nargin <= 6) || isempty(paramLogPlot)
				paramLogPlot = true;
			end
			validateattributes(paramLogPlot, {'logical'}, {'scalar', 'nonempty'}, '', 'paramLogPlot');

			if (nargin >= 3) && ~isempty(numIter)
				validateattributes(numIter, {'numeric'}, {'scalar', 'nonempty', '>=', 0.0}, '', 'numIter');
			else
				numIter = [];
			end

			if (nargin >= 4) && ~isempty(procData)
				validateattributes(procData, {'numeric'}, {'vector', 'numel', 2, 'nonempty', '>=', 0.0}, '', 'procData');
			else
				procData = [];
			end

			% Show process monitor if both procData and numIter are present
			procMonitor = ~isempty(procData) && ~isempty(numIter);

			if paramHist
				if ~isempty(obj.parameterTransform)
					trueP = obj.applyInverseTransform(obj.lastParams(:));
				else
					trueP = obj.lastParams(:);
				end
			end

			if isempty(hd)
				% Create new figure
				hd.figure = figure('Name', 'Parameter Fit');
				hd.sim = cell(length(obj.data), 1);
				hd.data = cell(length(obj.data), 1);
				hd.legend = cell(length(obj.data), 1);
				hd.title = cell(length(obj.data), 1);
				for i = 1:length(obj.data)
					subplot(length(obj.data), 1, i);

					localLegendNames = cell(2*length(obj.data{i}), 1);

					resSol = obj.lastResult{i}.solution.outlet;

					% Loop over all wavelengths to assemble simulation results
					curData = obj.data{i};
					curComp = obj.idxComp{i};
					curUnit = obj.idxUnit{i};
					curMask = obj.timeMask{i};
					curWeightData = obj.weightsData{i};
					curWeightWave = obj.weightsWave{i};
					for j = 1:length(curUnit)
						unitTrace = curUnit(j) + 1;
						curSol = resSol{unitTrace}(curMask{j}, :) * curComp{j};
						plotData = curData{j};

						if useWeights
							localWeights = obj.weightsExp(i) .* curWeightWave(j) .* curWeightData{j};
							plotData = curData{j} .* localWeights;
							curSol = curSol .* localWeights;
						end

						% Plot all wavelengths
						hold on;
						hdAll = plot(obj.models{i}.solutionTimes(curMask{j}), [plotData, curSol]);
						hold off;
						hd.data{i} = [hd.data{i}; hdAll(1)];
						hd.sim{i} = [hd.sim{i}; hdAll(2)];

						localLegendNames{2*j-1} = obj.waveNames{i}{j};
						localLegendNames{2*j} = sprintf('Simulated %s', obj.waveNames{i}{j});
					end

					% Add legend
					hd.legend{i} = legend(localLegendNames);
					set(hd.legend{i}, 'Box', 'off')
					set(hd.legend{i}, 'Location','NorthWest');

					hd.title{i} = title(sprintf('%s', obj.names{i}));
					grid on;
				end

				if procMonitor || paramHist
					hd.figureProcMon = figure('Name', 'Process Monitor');

					if procMonitor && paramHist
						subplot(1,2,1);
					end
					if procMonitor
						hd.res = semilogy(numIter, procData, 'x-');
						hdTemp = legend('Residual', 'Step size');
						set(hdTemp, 'Location','SouthWest');
						hd.resTitle = title(sprintf('Iter %d - Residual %g - Step %g', numIter, procData(1), procData(2)));
						grid on;
					end

					if procMonitor && paramHist
						subplot(1,2,2);
					end
					if paramHist
						if paramLogPlot
							hd.param = semilogy(numIter, abs(trueP), 'x-');
						else
							hd.param = plot(numIter, trueP, 'x-');
						end
						hd.paramTitle = title(sprintf('%g | ', trueP));
						grid on;
					end
				end
			else
				% Update data in existing figure
				for i = 1:length(obj.data)
					resSol = obj.lastResult{i}.solution.outlet;

					% Loop over all wavelengths to assemble simulation results
					curComp = obj.idxComp{i};
					curUnit = obj.idxUnit{i};
					curMask = obj.timeMask{i};
					curWeightData = obj.weightsData{i};
					curWeightWave = obj.weightsWave{i};
					for j = 1:length(curUnit)
						unitTrace = curUnit(j) + 1;
						curSol = resSol{unitTrace}(curMask{j}, :) * curComp{j};

						if useWeights
							curSol = curSol .* (obj.weightsExp(i) .* curWeightWave(j) .* curWeightData{j});
						end

						% Update data in plot
						set(hd.sim{i}(j), 'YData', curSol);
					end
				end

				if procMonitor
					set(hd.res(1), 'YData', [get(hd.res(1), 'YData'), procData(1)], 'XData', [get(hd.res(1), 'XData'), numIter]);
					set(hd.res(2), 'YData', [get(hd.res(2), 'YData'), procData(2)], 'XData', [get(hd.res(2), 'XData'), numIter]);
					set(hd.resTitle, 'String', sprintf('Iter %d - Residual %g - Step %g', numIter, procData(1), procData(2)));
				end

				if paramHist
					if paramLogPlot
						plotP = abs(trueP);
					else
						plotP = trueP;
					end
					
					for i = 1:length(trueP)
						set(hd.param(i), 'YData', [get(hd.param(i), 'YData'), plotP(i)], 'XData', [get(hd.param(i), 'XData'), numIter]);
					end
					set(hd.paramTitle, 'String', sprintf('%g | ', trueP));
				end
			end

			drawnow;
		end

	end

	methods (Access = 'protected')

		function rebuildInternalDataStructures(obj)
			%rebuildInternalDataStructures Rebuilds internal data structures due to removed or added experiments,
			% or a new parameter link.
			
			% Compute total number of data points
			obj.totalDataPoints = 0;
			for i = 1:length(obj.data)
				curData = obj.data{i};
				obj.totalDataPoints = obj.totalDataPoints + sum(cellfun(@(x) numel(x), curData));
			end

			numParams = cellfun(@(x) x.nVariableParameters, obj.models);
			obj.paramOffset = [0; cumsum(numParams(:))] + 1;
			maxTotalParams = sum(numParams);

			obj.localToGlobal = nan(maxTotalParams, 3); % Maps (idxExp, idxLocalParam) -> idxGlobalParam

			if isempty(obj.links)
				obj.links = cell(numel(obj.data), 1);
			elseif length(obj.links) < length(obj.data)
				obj.links = [obj.links(:); cell(length(obj.data) - length(obj.links), 1)];
			end

			% First pass: Assign IDs to global parameters
			curIdx = 1;
			lastId = 0;
			for i = 1:length(obj.models)
				curLinks = obj.links{i};
				for j = 1:numParams(i)

					% Check if there is a link to a parameter appearing in a previous experiment
					if ~isempty(curLinks)
						idxLine = find((curLinks(:, 1) == j) & (curLinks(:, 2) < i), 1);

						if ~isempty(idxLine)
							idxLine = find((obj.localToGlobal(:, 1) == curLinks(idxLine, 2)) & (obj.localToGlobal(:, 2) == curLinks(idxLine, 3)), 1);
							% Yes, so we look up its ID and assign it
							obj.localToGlobal(curIdx, :) = [i, j, obj.localToGlobal(idxLine, 3)];
						end
					end
					
					if isempty(curLinks) || isempty(idxLine)
						% No, so assign the next global id
						lastId = lastId + 1;
						obj.localToGlobal(curIdx, :) = [i, j, lastId];
					end

					curIdx = curIdx + 1;
				end
			end

			% Second pass: Build inverse map
			obj.globalToLocal = nan(lastId, 2); % Maps idxGlobalParam -> (idxExp, idxLocalParam)
			for i = 1:lastId
				idxLine = find(obj.localToGlobal(:, 3) == i, 1);
				obj.globalToLocal(i, :) = obj.localToGlobal(idxLine, 1:2);
			end
		end

	end
end

% =============================================================================
%  CADET - The Chromatography Analysis and Design Toolkit
%  
%  Copyright (C) 2008-2020: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. obj program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies obj distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
