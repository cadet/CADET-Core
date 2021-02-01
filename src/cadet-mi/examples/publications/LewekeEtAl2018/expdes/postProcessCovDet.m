function [ pulseLength, covDets, lowerQuant, upperQuant, failRates ] = postProcessCovDet( data, rel, doPlot )
%POSTPROCESSCOVDET Postprocesses data created by the RUNWORKFLOW() script for covariance determinants
%
%   Extracts the determinants of the covariance matrices, the optimization
%   failure rate, and the pulse lengths from the results of the RUNWORKFLOW()
%   script. 
%
%   Confidence bands are estimated for the covariance determinants by
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
%   POSTPROCESSCOVDET(DATA) postprocesses the cell array DATA obtained from a call
%   of RUNWORKFLOW and plots the results.
%
%   POSTPROCESSCOVDET(..., REL) determines whether the results are scaled
%   to the nominal parameter values. Defaults to true.
%
%   POSTPROCESSCOVDET(..., REL, DOPLOT) determines via DOPLOT if the 
%   results are plotted. Defaults to true.
%
%   [PULSELENGTH, COVDETS, LOWERQUANT, UPPERQUANT, FAILRATES] = POSTPROCESSCOVDET(...)
%   returns 
%      - pulseLength: A vector with ascending injection pulse lengths in 
%           seconds.
%      - covDets: Matrix with absolute value of the determinant of the
%           covariance matrix of the estimated parameters for each pulse length
%           (row) and parameter set (columns).
%      - lowerQuant: Matrix with lower bounds of 95 % confidence band on
%           absolute value of the determinant of the covariance matrix of the 
%           estimated parameters for each pulse length (row) and parameter set
%           (columns).
%      - upperQuant: Matrix with upper bounds of 95 % confidence band on
%           absolute value of the determinant of the covariance matrix of the 
%           estimated parameters for each pulse length (row) and parameter set
%           (columns).
%      - failRates: Vector with optimization failure rate for each pulse
%           length.
%
%   See also RUNWORKFLOW, POSTPROCESSACCURACY

% Copyright: (C) 2008-2016 The CADET Authors
%            See the license note at the end of the file.

	if (nargin <= 1) || isempty(rel)
		rel = true;
	end

	if (nargin <= 2) || isempty(doPlot)
		doPlot = true;
	end

	% Number of bootstrap samples
	nBSsamples = 1000;
	% Confidence level (95 %)
	conf = 0.95;

	% Comment out the following line to avoid resetting the RNG
	rng(0, 'twister');

	% Preallocate space
	pulseLength = data{1}.inTimes;
	lowerQuant = zeros(length(pulseLength), length(data));
	upperQuant = zeros(length(pulseLength), length(data));
	failRates = zeros(length(pulseLength), length(data));
	covDets = zeros(length(pulseLength), length(data));

	for mod = 1:length(data)

		% Collect metrics
		failRates(:, mod) = data{mod}.failRate;
		covDets(:, mod) = data{mod}.covDets;

		% Remove dimensionality and normalize
		if (size(data{mod}.samples, 2) > 3)
			if rel
				covDets = covDets ./ (1.14 * 4.88)^2;
			end
			covDets(:, mod) = covDets(:, mod).^(1 / (size(data{mod}.samples, 2) - 2));
		elseif (mod == 2) && rel
			covDets(:, mod) = covDets(:, mod) ./ 1.14^2; 
		elseif (mod == 3) && rel
			covDets(:, mod) = covDets(:, mod) ./ 4.88^2; 
		end

		for j = 1:length(pulseLength)

			% Filter invalid samples
			idxValid = data{mod}.samples(:, end, j) > 0;
			samples = data{mod}.samples(idxValid, 1:end-2, j);

			% Bootstrapping
			covDetsSamples = zeros(1, nBSsamples);
			for i = 1:nBSsamples
				% Generate new samples by randomly drawing with replacement
				covDetsSamples(i) = abs(det(cov(samples(randi(size(samples, 1), size(samples, 1), 1), :))));
			end

			% Remove dimensionality and normalize
			if (size(samples, 2) > 1)
				if rel
					covDetsSamples = covDetsSamples ./ (1.14 * 4.88)^2;
				end
				covDetsSamples = covDetsSamples.^(1 / size(samples, 2));
			elseif (mod == 2) && rel
				covDetsSamples = covDetsSamples ./ 1.14^2; 
			elseif (mod == 3) && rel
				covDetsSamples = covDetsSamples ./ 4.88^2; 
			end

			% Compute quantiles of new samples to obtain confidence band
			d = sort(covDetsSamples);

			quant = [0 (0.5:(nBSsamples-0.5))./nBSsamples 1]';
			vals = [d(1); d(:); d(end)];
			temp = interp1(quant, vals, ((1-conf) / 2) + [0, conf]);
			lowerQuant(j, mod) = temp(1);
			upperQuant(j, mod) = temp(2);
		end

	end

	% Do some plotting
	if doPlot
		errorbar(repmat(pulseLength, 1, length(data)), covDets, covDets - lowerQuant, upperQuant - covDets);
		grid on;
		legend('Keq & Qmax', 'Keq', 'Qmax');
		title('Precision');
		xlabel('Pulse length [s]');
		if rel
			ylabel('( det( Cov( p ) ) / prod(pTrue^2) )^(1 / nParam)');
		else
			ylabel('( det( Cov( p ) ) )^(1 / nParam)');
		end
		set(gca, 'YScale', 'log');
	end
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
