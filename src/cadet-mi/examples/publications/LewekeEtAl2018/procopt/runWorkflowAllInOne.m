function [ data ] = runWorkflowAllInOne(nSamples, nMultiStart, saveToFile, colLength, elutionConstraint, algorithm)
%RUNWORKFLOWALLINONE Runs the whole workflow (sampling, optimization) of the combined (gradient, cut points) yield optimization
%  
%   Runs the whole process optimization workflow consisting of first sampling
%   initial parameters and then starting a gradient-based optimizer on the
%   best samples, thus, using a multi-start strategy.
%  
%   See the function BILINEARGRADIENTPROCESSYIELD for details on the parameters
%   and the underlying optimization scenario. See SAMPLEINITIALPOINTSALLINONE 
%   for details on the sampling.
%
%   RUNWORKFLOWALLINONE(NSAMPLES) runs the workflow and uses NSAMPLES samples
%   from which the best 10 are used for optimization (multi-start).
%
%   RUNWORKFLOWALLINONE(..., NMULTISTART) uses the best NMULTISTART samples for
%   optimization (multi-start).
%
%   RUNWORKFLOWALLINONE(..., NMULTISTART, SAVETOFILE) determines via SAVETOFILE 
%   whether final and intermediate results are saved to file. Defaults to true.
%
%   RUNWORKFLOWALLINONE(..., NMULTISTART, SAVETOFILE, COLLENGTH) additionally sets the column
%   length which defaults to 0.014m.
%
%   RUNWORKFLOWALLINONE(..., NMULTISTART, SAVETOFILE, COLLENGTH, ELUTIONCONSTRAINT)
%   additionally sets the constraint type that enforces complete elution. If
%   ELUTIONCONSTRAINT = 1, the concentration at the end of the profile is constrained to 1e-6 mM.
%   If ELUTIONCONSTRAINT = 2, the eluted mass has to be at least 99.5 % of the injected mass.
%
%   RUNWORKFLOWALLINONE(..., NMULTISTART, SAVETOFILE, COLLENGTH, ELUTIONCONSTRAINT, ALGORITHM) 
%   sets the algorithm that FMINCON uses for optimization (e.g., 'interior-point', 'sqp').
%
%   DATA = RUNWORKFLOWALLINONE(...) returns a cell array of structs with the fields
%      - optCut: Optimal cut points in s
%      - optGrad: Optimal gradient shape parameters
%      - optYield: Optimized yield
%      - optPurity: Purity at optimum
%      - optOverlaps: Overlaps of the components at the optimum
%      - optResult: Struct with additional meta info returned by the optimizer.
%      - initShape: Initial gradient parameters
%      - initCut: Initial cut points
%
%   See also SAMPLEINITIALPOINTSALLINONE, BILINEARGRADIENTPROCESSYIELD

% Copyright: © 2015 Samuel Leweke, Eric von Lieres
%            See the license note at the end of the file.

	% Set default values
	if (nargin <= 1) || isempty(nMultiStart)
		nMultiStart = 10;
	end

	if (nargin <= 2) || isempty(saveToFile)
		saveToFile = true;
	end

	if (nargin <= 3) || isempty(colLength)
		colLength = 0.014;
	end

	if (nargin <= 4)
		elutionConstraint = 1;
	end

	if (nargin <= 5)
		algorithm = 'interior-point';
	end

	% Minimum purity requirement
	minPurity = 0.95;

	% Generate samples
	optimalCutPointSampling = true;
	[gradientShapes, cutTimes, yields, purities, overlaps] = sampleInitialPointsAllInOne(nSamples, colLength, true, optimalCutPointSampling);

	if saveToFile
		save(['CombinedOpt-Samples.mat'], 'gradientShapes', 'cutTimes', 'yields', 'purities', 'overlaps', 'optimalCutPointSampling');
	end
	
	% Filter samples not meeting the purity restriction
	subset = purities >= minPurity;
	purities = purities(subset);
	yields = yields(subset);
	gradientShapes = gradientShapes(subset, :);
	overlaps = overlaps(subset, :);
	cutTimes = cutTimes(subset, :);

	% Sort by yield
	[~, idxSorted] = sort(-yields);
	purities = purities(idxSorted);
	yields = yields(idxSorted);
	gradientShapes = gradientShapes(idxSorted, :);
	overlaps = overlaps(idxSorted, :);
	cutTimes = cutTimes(idxSorted, :);
	
	% Take the best nMultiStart samples (or less)
	nSamples = min(length(purities), nMultiStart);
	data = cell(nSamples, 1);
	for i = 1:nSamples
		% Run the optimizer
		[optCut, optGrad, optYield, optPurity, optOverlaps, optResult] = bilinearGradientProcessYield(cutTimes(i, :), gradientShapes(i, :), colLength, elutionConstraint, algorithm);

		% Refine cut points
		[cutPoints, optimalYield, optimalPurity] = bilinearGradientProcessCutTimes(optGrad, [], colLength);

		% Save data
		res.optGrad = optGrad;
		res.optCut = optCut;
		res.optYield = optYield;
		res.optPurity = optPurity;
		res.optCutRef = cutPoints;
		res.optYieldRef = optimalYield;
		res.optPurityRef = optimalPurity;
		res.optResult = optResult;
		res.optOverlaps = optOverlaps;
		res.initShape = gradientShapes(i,:);
		res.initCut = cutTimes(i,:);
		data{i} = res;

		if saveToFile
			save(['comb-run' num2str(i) '.mat'], 'optCut', 'optGrad', 'optYield', 'optPurity', 'cutPoints', 'optimalYield', 'optimalPurity', 'optOverlaps', 'optResult', 'i');
		end

	end
	
	if saveToFile
		samples = [];
		samples.purities = purities;
		samples.yields = yields;
		samples.gradientShapes = gradientShapes;
		samples.overlaps = overlaps;
		samples.cutTimes = cutTimes;

		save(['CombinedOpt-Results.mat'], 'data', 'samples', 'nMultiStart', 'algorithm', 'elutionConstraint', 'optimalCutPointSampling');
	end
	
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
