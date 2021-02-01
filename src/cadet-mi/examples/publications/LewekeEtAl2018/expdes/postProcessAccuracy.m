function [ pulseLength, dist, relDist, confRelDistLo, confRelDistUp ] = postProcessAccuracy( data, doPlot )
%POSTPROCESSACCURACY Postprocesses data created by the RUNWORKFLOW() script for accuracy
%
%   Calculates mean absolute and relative deviation from the true parameters,
%   and the pulse lengths from the results of the RUNWORKFLOW() script. 
%
%   Confidence bands are estimated for the mean relative deviation by
%   a resampling Monte Carlo bootstrap method.
%
%   See RUNWORKFLOW() for details. The model considered here describes
%   affinity chromatography of lysozyme on Cibacron Blue Sepharose CL-6B at 
%   pH 2.7. Model parameters are based on benchmark 1 of the publication:
%   A. Püttmann, S. Schnittert, U. Naumann & E. von Lieres (2013).
%   Fast and accurate parameter sensitivities for the general rate model of
%   column liquid chromatography.
%   Computers & Chemical Engineering, 56, 46–57.
%   doi:10.1016/j.compchemeng.2013.04.021
%
%   POSTPROCESSACCURACY(DATA) postprocesses the cell array DATA obtained from a call
%   of RUNWORKFLOW and plots the results.
%
%   POSTPROCESSACCURACY(..., DOPLOT) determines via DOPLOT if the results are plotted.
%
%   [PULSELENGTH, DIST, RELDIST, CONFRELDISTLO, CONFRELDISTUP] = POSTPROCESSACCURACY(...)
%   returns 
%      - pulseLength: A vector with ascending injection pulse lengths in 
%           seconds.
%      - dist: Matrix with mean absolute deviation of the estimated parameters
%           from the true parameters for each pulse length (row) and parameter
%           (columns). Note that jointly estimated parameters appear as separate
%           columns.
%      - relDist: Matrix with mean relative deviation of the estimated parameters
%           from the true parameters for each pulse length (row) and parameter
%           (columns). Note that jointly estimated parameters appear as separate
%           columns.
%      - confRelDistLo: Matrix with lower bounds of 95 % confidence band on
%           mean relative deviation for each pulse length (row) and parameter
%           (columns). Note that jointly estimated parameters appear as separate
%           columns.
%      - confRelDistUp: Matrix with upper bounds of 95 % confidence band on
%           mean relative deviation for each pulse length (row) and parameter
%           (columns). Note that jointly estimated parameters appear as separate
%           columns.
%
%   See also RUNWORKFLOW, POSTPROCESSCOVDET

% Copyright: (C) 2008-2016 The CADET Authors
%            See the license note at the end of the file.

	if (nargin <= 1) || isempty(doPlot)
		doPlot = true;
	end

	% Confidence level (95 %)
	conf = 0.95;

	pulseLength = data{1}.inTimes;

	% Preallocate space
	dist = zeros(length(pulseLength), 4);
	relDist = zeros(length(pulseLength), 4);
	confRelDistLo = zeros(length(pulseLength), 4);
	confRelDistUp = zeros(length(pulseLength), 4);
	for mod = 1:length(data)

		localParamVals = [1.14, 4.88];

		if mod == 2
			% Estimate only k_a and fix q_max
			localParamVals = localParamVals(1);
		elseif mod == 3
			% Estimate only q_max and fix k_a
			localParamVals = localParamVals(2);
		end

		for j = 1:length(pulseLength)

			% Filter invalid samples
			idxValid = data{mod}.samples(:, end, j) > 0;
			samples = data{mod}.samples(idxValid, 1:end-2, j);

			diff = abs(samples - repmat(localParamVals, size(samples, 1), 1));
			relDiff = abs(samples - repmat(localParamVals, size(samples, 1), 1)) ./ repmat(localParamVals, size(samples, 1), 1);

			if mod == 1
				dist(j, 1:2) = mean(diff, 1);
				relDist(j, 1:2) = mean(relDiff, 1);

				temp = confInterval(relDiff(:, 1), conf);
				confRelDistLo(j, 1) = temp(1);
				confRelDistUp(j, 1) = temp(2);

				temp = confInterval(relDiff(:, 2), conf);
				confRelDistLo(j, 2) = temp(1);
				confRelDistUp(j, 2) = temp(2);
			else
				dist(j, mod+1) = mean(diff);
				relDist(j, mod+1) = mean(relDiff);

				temp = confInterval(relDiff, conf);
				confRelDistLo(j, mod+1) = temp(1);
				confRelDistUp(j, mod+1) = temp(2);
			end
		end

	end

	% Do some plotting
	if doPlot
		subplot(1,2,1);
		errorbar(repmat(pulseLength, 1, 4), relDist, relDist - confRelDistLo, confRelDistUp - relDist);
		grid on;
		grid on;
		legend('Keq Comb', 'Qmax Comb', 'Keq', 'Qmax');
		xlabel('Pulse length [s]');
		ylabel('mean( abs(p - p_{true}) / p_{true} )');
		set(gca, 'YScale', 'log');

		subplot(1,2,2);
		semilogy(pulseLength, dist);
		grid on;
		legend('Keq Comb', 'Qmax Comb', 'Keq', 'Qmax');
		xlabel('Pulse length [s]');
		ylabel('mean( abs(p - p_{true}) )');
	end
end

function [confInt] = confInterval(data, confLevel)
	d = sort(data);
	quant = [0 (0.5:(length(data)-0.5))./length(data) 1]';
	vals = [d(1); d(:); d(end)];
	confInt = interp1(quant, vals, ((1-confLevel) / 2) + [0, confLevel]);
end

% =============================================================================
%  CADET
%  
%  Copyright (C) 2008-2016: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
