function data = runWorkflowTwoStage(nSamples, nMultiStart, saveToFile, colLength, absMode, elutionConstraint, algorithm)
%RUNWORKFLOWTWOSTAGE Runs the whole process optimization workflow (sampling, optimization)
%
%   Runs the whole process optimization workflow consisting of first sampling
%   initial parameters and then starting a gradient-based optimizer on the
%   best samples, thus, using a multi-start strategy.
%  
%   See the functions BILINEARGRADIENTPROCESSOVERLAPSUM and
%   BILINEARGRADIENTPROCESSOVERLAPMAX for details on the parameters and
%   the underlying optimization scenario. See SAMPLEINITIALPOINTSTWOSTAGE
%   for details on the sampling and BILINEARGRADIENTPROCESSCUTTIMES for
%   cut point optimization.
%  
%   The optimization problem is separated into two sub problems:
%    1. Optimize gradient shape to obtain minimum peak overlaps
%    2. Use fixed chromatogram to find optimal cut points
%
%   RUNWORKFLOWTWOSTAGE(NSAMPLES) runs the whole workflow with NSAMPLES
%   samples using the best 10 samples for optimization.
%
%   RUNWORKFLOWTWOSTAGE(..., NMULTISTART) uses the best NMULTISTART samples
%   for optimization (multi-start).
%
%   RUNWORKFLOWTWOSTAGE(..., NMULTISTART, SAVETOFILE) determines via SAVETOFILE
%   if final and intermediate results are saved to file. Defaults to true.
%
%   RUNWORKFLOWTWOSTAGE(..., NMULTISTART, SAVETOFILE, COLLENGTH) additionally
%   sets the column length which defaults to 0.014m.
%
%   RUNWORKFLOWTWOSTAGE(..., NMULTISTART, SAVETOFILE, COLLENGTH, ABSMODE) additionally sets
%   the absolute value function (regularization) used.
%
%   RUNWORKFLOWTWOSTAGE(..., NMULTISTART, SAVETOFILE, COLLENGTH, ABSMODE, ELUTIONCONSTRAINT) 
%   additionally sets the constraint type that enforces complete elution. If
%   ELUTIONCONSTRAINT = 1, the concentration at the end of the profile is constrained
%   to 1e-6 mM. If ELUTIONCONSTRAINT = 2, the eluted mass has to be at least 99.5 %
%   of the injected mass.
%
%   RUNWORKFLOWTWOSTAGE(..., NMULTISTART, SAVETOFILE, COLLENGTH, ABSMODE, ELUTIONCONSTRAINT, ALGORITHM) 
%   additionally sets the algorithm used by FMINCON (e.g., 'interior-point', 'sqp').
%
%   DATA = RUNWORKFLOWTWOSTAGE(...) returns an NSAMPLES x 3 x 2 cell
%   array with structs containing the optimization results. The first
%   dimension specifies the initial sample used, second the absMode,
%   and third if maximum of overlaps (1) or sum of overlaps (2) is
%   minimized. The structs have the following fields:
%      - optCut: Optimal cut points in s
%      - optGrad: Optimal gradient shape parameters
%      - optYield: Optimized yield
%      - optPurity: Purity at optimum
%      - optOverlaps: Overlaps of the components at the optimum
%      - optResult: Struct with additional meta info returned by the optimizer.
%      - initShape: Initial gradient parameters
%
%   See also BILINEARGRADIENTPROCESSOVERLAPSUM, BILINEARGRADIENTPROCESSOVERLAPMAX,
%      RUNOPTSCENARIO, BILINEARGRADIENTPROCESSCUTTIMES, SAMPLEINITIALPOINTSTWOSTAGE

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
		absMode = 0;
	end

	if (nargin <= 5)
		elutionConstraint = 1;
	end

	if (nargin <= 6)
		algorithm = 'interior-point';
	end

	% Generate samples
	[gradientShapes, objVals] = sampleInitialPointsTwoStage(nSamples, colLength, true);
	
	if saveToFile
		save(['ThreeCompMultOpt-Samples.mat'], 'gradientShapes', 'objVals');
	end

	% Find best nMultistart samples for each objective function and store them in rows of samples
	nMultiStart = min(nSamples, nMultiStart);
	samples = zeros(2 * nMultiStart, size(gradientShapes, 2));

	% Handle sum objective
	[~, idxSorted] = sort(objVals(:, 1));
	samples(1:nMultiStart, :) = gradientShapes(idxSorted(1:10), :);

	% Handle max objective
	[~, idxSorted] = sort(objVals(:, 2));
	samples(nMultiStart+1:2*nMultiStart, :) = gradientShapes(idxSorted(1:10), :);
	
	% Remove all duplicates
	samples = unique(samples, 'rows');
	
	% Run all combinations of initial samples and sum / max objective
	data = cell(size(samples, 1), 2);
	for i = 1:size(samples, 1)
		% Max of all overlaps
		[gradientShape, overlaps, optResult] = runOptScenario(false, samples(i, :), absMode, elutionConstraint, colLength, algorithm, saveToFile, i);
		[cutPoints, optimalYield, optimalPurity] = bilinearGradientProcessCutTimes(gradientShape, 1, colLength, true);

		% Save data
		res.initShape = samples(i,:);
		res.optGrad = gradientShape;
		res.optOverlaps = overlaps;
		res.optResult = optResult;
		res.optCut = cutPoints;
		res.optYield = optimalYield;
		res.optPurity = optimalPurity;
		data{i, 1} = res;
		
		% Sum of all overlaps
		[gradientShape, overlaps, optResult] = runOptScenario(true, samples(i, :), absMode, elutionConstraint, colLength, algorithm, saveToFile, i);
		[cutPoints, optimalYield, optimalPurity] = bilinearGradientProcessCutTimes(gradientShape, 1, colLength, true);

		% Save data
		res.initShape = samples(i,:);
		res.optGrad = gradientShape;
		res.optOverlaps = overlaps;
		res.optResult = optResult;
		res.optCut = cutPoints;
		res.optYield = optimalYield;
		res.optPurity = optimalPurity;
		data{i, 2} = res;
	end
	
	if saveToFile
		save(['ThreeCompMultOpt-Results.mat'], 'data', 'samples', 'absMode', 'elutionConstraint', 'colLength', 'algorithm');
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
