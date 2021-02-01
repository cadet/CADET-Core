function runWorkflow(nSamples, nBest, transformed, saveToFile)
%RUNWORKFLOW Runs the whole parameter estimation workflow (sampling, optimization)
%
%   Runs the whole parameter estimationworkflow consisting of first sampling
%   initial parameters and then starting a gradient-based optimizer on the
%   best samples, thus, using a multi-start strategy.
%
%   See the function RUNFITTING for details on the parameters and the
%   underlying optimization scenario. See SAMPLEINITIALPOINTS for details
%   on the sampling.
%
%   RUNWORKFLOW(NSAMPLES, NBEST) runs the whole workflow with NSAMPLES
%   initial samples and takes the best NBEST samples for optimization.
%
%   RUNWORKFLOW(..., TRANSFORMED) determines via the flag TRANSFORMED
%   if a parameter transformation is applied for SMA_NU in the optimization
%   process (fixes ordering of components). Defaults to false.
%
%   RUNWORKFLOW(..., TRANSFORMED, SAVETOFILE) additionally controls wheter
%   the results are saved to file by the flag SAVETOFILE. Defaults to true.
%
%   See also SAMPLEINITIALPOINTS, RUNFITTING, RUNFITTINGTRANSFORMED

% Copyright: (C) 2008-2016 The CADET Authors
%            See the license note at the end of the file.

	if (nargin <= 2) || isempty(transformed)
		transformed = false;
	end

	if (nargin <= 3) || isempty(saveToFile)
		saveToFile = true;
	end
	
	% Generate samples
	[samples, objVals] = sampleInitialPoints(nSamples, transformed);
	if saveToFile
		if transformed
			fileName = 'ThreeCompParEst-Samples-trans.mat';
		else
			fileName = 'ThreeCompParEst-Samples.mat';
		end
		save(fileName, 'samples', 'objVals');
	end

	% Make sure we don't use more initial points than samples
	nBest = min(nSamples, nBest);
	
	% Find best samples
	[~, idx] = sort(objVals);
	samples = samples(idx(1:nBest), :);
	
	% Create storage
	estimatedParams = zeros(nBest, size(samples, 2));
	residuals = zeros(nBest, 1);
	
	% Run parameter fitting for each initial sample
	for i = 1:nBest
		initParams = samples(i,:);
		if transformed
			[optParams, residual] = runFittingTransformed(initParams);
		else
			[optParams, residual] = runFitting(initParams);
		end
		residuals(i) = residual;
		if residual == -1
			estimatedParams(i, :) = initParams;
		else
			estimatedParams(i, :) = optParams;
		end

		if saveToFile
			% Construct file name for results
			nameBase = ['ThreeCompParEst-init' num2str(i)];
			if transformed
				nameBase = [nameBase '-trans'];
			end
			% Save results in MAT file
			save([nameBase '.mat'], 'optParams', 'i', 'residual', 'initParams');
			% Save figure
			saveas(gcf, [nameBase '.fig']);
		end
	end

	if saveToFile
		if transformed
			fileName = 'ThreeCompParEst-Results-trans.mat';
		else
			fileName = 'ThreeCompParEst-Results.mat';
		end
		save(fileName, 'estimatedParams', 'residuals', 'samples');
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
