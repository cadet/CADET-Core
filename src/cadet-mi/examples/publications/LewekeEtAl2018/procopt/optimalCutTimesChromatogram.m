function [cutPoints, optimalYield, optimalPurity] = optimalCutTimesChromatogram(chrom, idxTarget, injMass, mode, fast, quiet)
%OPTIMALCUTTIMESCHROMATOGRAM Calculates optimal cut times with respect to purity and yield constraints
%  
%   Determines optimal cut points in a chromatogram in order to collect a
%   target component subject to purity and yield constraints.
%  
%   There are two modes:
%     1. Maximize yield subject to a 95 % purity constraint
%     2. Maximize purity subject to a 70 % yield constraint
%  
%   Introducing the mass m_i of component i between the cut points t_1 and 
%   t_2 as
%             / t_2
%      m_i = |      c_i(t) dt,
%            / t_1
%  
%   yield y_i is defined by
%      y_i = m_i / injMass_i
%   and purity is given by
%      p_i = m_i / (Sum_j m_j).
%   For the two cases the optimization problems are formulated as follows:
%     1:  max  y_{idxTarget}
%         s.t. p_{idxTarget} >= 0.95
%              t_1 <= t_2
%     2:  max  p_{idxTarget}
%         s.t. y_{idxTarget} >= 0.7
%              t_1 <= t_2
%  
%   The given chromatogram is converted into a cubic spline. In order to
%   evaluate m_i, the spline, which is a piecewise polynomial, is
%   integrated analytically by calculating its anti-derivative.
%
%   OPTIMALCUTTIMESCHROMATOGRAM(CHROM, IDXTARGET, INJMASS) performs cut
%   time optimization for the given chromatogram CHROM, which is a matrix
%   whose first column denotes time and every other column is a
%   concentration at the given time point (i.e., a chromatogram of a 
%   component). The target component is given by the 1-based index
%   IDXTARGET. The array INJMASS is supposed to contain the injected masses
%   of the different components in the same order as in CHROM.
%
%   OPTIMALCUTTIMESCHROMATOGRAM(..., MODE) selects the objective of the 
%   optimization (1 for maximum yield [default], 2 for maximum purity).
%
%   OPTIMALCUTTIMESCHROMATOGRAM(..., MODE, FAST) additionally determines
%   whether loose tolerances and low maximum iteration numbers (TRUE) are
%   used during optimization or not (FALSE, default).
%
%   OPTIMALCUTTIMESCHROMATOGRAM(..., MODE, FAST, QUIET) additionally
%   determines whether the optimizer's iterations are visualized and
%   statistics printed (FALSE, default) or no output is generated during 
%   optimization (TRUE).
%
%   [CUTPOINTS, OPTIMALYIELD, OPTIMALPURITY] = OPTIMALCUTTIMESCHROMATOGRAM(...)
%   returns the cut times in CUTPOINTS, the achieved yield OPTIMALYIELD
%   and purity OPTIMALPURITY at the optimum.
%
% See also BILINEARGRADIENTPROCESSCUTTIMES

% Copyright: © 2015 Samuel Leweke, Eric von Lieres
%            See the license note at the end of the file.

	if (nargin <= 3) || isempty(mode)
		mode = 1;
	end
	
	if (nargin <= 4) || isempty(fast)
		fast = false;
	end
	
	if (nargin <= 5) || isempty(quiet)
		quiet = false;
	end
	
	% Minimum purity constraint (95 %)
	minPurity = 0.95;

	% Minimum yield constraint (60 %)
	minYield = 0.7;
		
	% Create splines and anti-derivatives from chromatograms
	chromSpl = cell(size(chrom, 2) - 1, 1);
	chromInts = cell(size(chrom, 2) - 1, 1);
	for i = 1:length(chromSpl)
		chromSpl{i} = pchip(chrom(:,1), chrom(:,1+i));
		chromInts{i} = ppint(chromSpl{i});
	end
	
	% Find time of target component's peak
	[~, initPos] = max(chrom(:, idxTarget + 1));
	initPos = chrom(initPos, 1);
	
	% Bounds and initial cut points
	lb = [0, 0];
	ub = chrom(end, 1) .* [1, 1];
	initPos = initPos .* [0.99, 1.01];
	
	% Storage for nested functions
	hdPlot = [];

	% Set options and perform optimization of the different scenarios
	options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Diagnostics', 'off', ...
		'GradObj', 'on', 'GradConstr', 'on', 'Display', 'off');

	if fast
		options = optimoptions(options, 'MaxIter', 50, 'TolCon', 1e-6, 'TolFun', 1e-6, 'TolX', 1e-8);
	else
		options = optimoptions(options, 'MaxIter', 500, 'TolCon', 1e-10, 'TolFun', 1e-10, 'TolX', 1e-12);
	end
	
	if ~quiet
		options = optimoptions(options, 'Display', 'iter', 'OutputFcn', @progressMonitor);
	end
	
	if mode == 1
		[cutPoints, fval] = fmincon(@yield, initPos, [1, -1], [0], [], [], lb, ub, @purity, options);
	else
		[cutPoints, fval] = fmincon(@purity, initPos, [1, -1], [0], [], [], lb, ub, @yield, options);
	end
	
	% Calculate optimal purity and yield
	optimalYield = -yield(cutPoints);
	optimalPurity = -purity(cutPoints);
	
	function varargout = yield(x)
		% Calculates the yield
		
		% Evaluate yield formula
		y = mass(x, chromInts, idxTarget) / injMass(idxTarget);
		% Compute gradient
		grady = [ppval(chromSpl{idxTarget}, x)] ./ injMass(idxTarget);
		% First element needs to be negated because it is t_1 (lower
		% integral boundary)
		grady(1) = -grady(1);
		
		if nargout == 1
			% Call as objective function => max yield
			varargout{1} = -y;
		elseif nargout == 2
			% Call as objective function => max yield, return [y, grady]
			varargout{1} = -y;
			varargout{2} = -grady;
		elseif nargout == 4
			% Call as constraint => return [y,eq,grady,gradeq]
			varargout{1} = minYield - y;
			varargout{2} = [];
			varargout{3} = -grady';
			varargout{4} = [];
		end
	end

	function varargout = purity(x)
		% Calculates the purity
		
		% Precompute masses
		masses = arrayfun(@(idx) mass(x, chromInts, idx), 1:length(chromSpl));
		mSum = sum(masses);
		mTarget = masses(idxTarget);

		% Evaluate purity formula
		y = mTarget / mSum;
		% First element needs to be negated because it is t_1 (lower
		% integral boundary). Note that the product / quotient rule and the
		% chain rule have to be applied.
		grady = [ppval(chromSpl{idxTarget}, x) ./ mSum - mTarget * sum(cell2mat(cellfun(@(s) ppval(s, x), chromSpl, 'UniformOutput', false)), 1) / mSum^2];
		grady(1) = -grady(1);
		
		if nargout == 1
			% Call as objective function => max purity
			varargout{1} = -y;
		elseif nargout == 2
			% Call as objective function => max purity, return [y, grady]
			varargout{1} = -y;
			varargout{2} = -grady;
		elseif nargout == 4
			% Call as constraint => return [y,eq,grady,gradeq]
			varargout{1} = minPurity - y;
			varargout{2} = [];
			varargout{3} = -grady';
			varargout{4} = [];
		end
	end

	function stop = progressMonitor(x, optimValues, state)
		%PROGRESSMONITOR Callback function for progress report invoked by fmincon

		stop = false;
		if ~strcmp(state, 'iter')
			return;
		end

		% Calculate yield and purity
		masses = arrayfun(@(idx) mass(x, chromInts, idx), 1:length(chromSpl));
		mSum = sum(masses);
		mTarget = masses(idxTarget);
		yield = mTarget / injMass(idxTarget);
		purity = mTarget / mSum;

		% Plot

		if isempty(hdPlot)			
			hdPlot.figureA = figure('Name', 'Optimization');

			hdPlot.ax = subplot(1, 3, 1);
			plot(chrom(:,1), chrom(:,2:end));
			grid on;

			% Plot lines for start and end point
			ylim = get(gca,'ylim');
			hold on;
			hdPlot.cutS = line([x(1), x(1)], [0, ylim(2)], 'LineStyle','-', 'Color','k');
			hdPlot.cutE = line([x(2), x(2)], [0, ylim(2)], 'LineStyle','-', 'Color','k');
			hold off;
			set(gca, 'ylim', [0, ylim(2)]);
			set(gca, 'xlim', [min(chrom(:,1)), max(chrom(:,1))]);

			hdPlot.titleLeft = title(sprintf('From %g to %g => Yield %g and Purity %g', [x(:); yield; purity]));

			subplot(1, 3, 2);
			hdPlot.opt = semilogy(optimValues.iteration, [-optimValues.fval, optimValues.constrviolation, optimValues.firstorderopt, optimValues.stepsize], 'x-');
			hdTemp = legend('Function value', 'Constraint violation', 'Optimality', 'Step size');
			set(hdTemp, 'Location','SouthWest');
			hdPlot.optTitle = title(sprintf('Iteration %d', optimValues.iteration));
			grid on;

			subplot(1, 3, 3);
			hdPlot.param = plot(optimValues.iteration, x, 'x-');
			hdPlot.paramTitle = title('Parameters');
			grid on;
		else
			set(hdPlot.cutS, 'XData', [x(1), x(1)]);
			set(hdPlot.cutE, 'XData', [x(2), x(2)]);
			set(hdPlot.titleLeft, 'String', sprintf('From %g to %g => Yield %g and Purity %g', [x(:); yield; purity]));

			set(hdPlot.opt(1), 'YData', [get(hdPlot.opt(1), 'YData'), -optimValues.fval], 'XData', [get(hdPlot.opt(1), 'XData'), optimValues.iteration]);
			set(hdPlot.opt(2), 'YData', [get(hdPlot.opt(2), 'YData'), optimValues.constrviolation], 'XData', [get(hdPlot.opt(2), 'XData'), optimValues.iteration]);
			set(hdPlot.opt(3), 'YData', [get(hdPlot.opt(3), 'YData'), optimValues.firstorderopt], 'XData', [get(hdPlot.opt(3), 'XData'), optimValues.iteration]);
			set(hdPlot.opt(4), 'YData', [get(hdPlot.opt(4), 'YData'), optimValues.stepsize], 'XData', [get(hdPlot.opt(4), 'XData'), optimValues.iteration]);
			set(hdPlot.optTitle, 'String', sprintf('Iteration %d', optimValues.iteration));

			for i = 1:numel(x)
				set(hdPlot.param(i), 'YData', [get(hdPlot.param(i), 'YData'), x(i)], 'XData', [get(hdPlot.param(i), 'XData'), optimValues.iteration]);
			end
		end
		drawnow;
	end
end

function m = mass(t, chromInts, idxTarget)
%MASS Evaluates the mass of the target component between the given cut points
% Evaluates the anti-derivative of the target component's chromatogram in
% order to calculate the definite integral.
	m = diff(ppval(chromInts{idxTarget}, t));
end

% =============================================================================
%  
%  Copyright © 2015: Samuel Leweke¹, Eric von Lieres¹
%                                      
%    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
