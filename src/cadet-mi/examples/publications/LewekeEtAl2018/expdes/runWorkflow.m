function [data] = runWorkflow(nSamples, saveToFile)
%RUNWORKFLOW Runs the experimental design workflow (Monte-Carlo, post-processing)
%
%   The motivation of this study is given by the question whether, given a
%   specific amount of protein, it is better to use long pulses with low
%   concentration or short pulses with high concentrations to estimate
%   isotherm parameters. The model considered here describes affinity
%   chromatography of lysozyme on Cibacron Blue Sepharose CL-6B at pH 2.7. 
%   Model parameters are based on benchmark 1 of the following publication:
%   A. Püttmann, S. Schnittert, U. Naumann & E. von Lieres (2013).
%   Fast and accurate parameter sensitivities for the general rate model of
%   column liquid chromatography.
%   Computers & Chemical Engineering, 56, 46–57.
%   doi:10.1016/j.compchemeng.2013.04.021
%  
%   A Monte-Carlo approach is employed to tackle this question. Different
%   base pulse shapes (each injecting the same amount of protein) are
%   considered and the effect of measurement uncertainty in injected mass and
%   pulse length are studied. Three different estimation scenarios are
%   considered: Estimation of adsorption constant (k_a) and capacity (q_max),
%   estimation of k_a only, and estimation of q_max only.
%   For each sampled pulse shape a parameter estimation is conducted to
%   recover the isotherm parameters using the original pulse shape as input.
%   The samples are used to calculate covariance matrices and precision
%   metrics for each base pulse shape.
%  
%   This function drives the Monte-Carlo sampling (sampleParameters function)
%   and uses the results for post-processing.
%
%   RUNWORKFLOW(NSAMPLES, SAVETOFILE) runs the full workflow with NSAMPLES samples
%   per pulse shape and estimation scenario. The flag SAVETOFILE controls whether
%   results are saved to file.
%
%   DATA = RUNWORKFLOW(...) returns all results in a cell array DATA with 3
%   elements which correspond to the three parameter estimation scenarios. Each
%   cell contains a struct with the following elements:
%       o samples: 3D-array with parameter samples obtained for each pulse
%           shape. Dimensions are nSamples x (nParameters + 2) x nShapes.
%           The additional two columns contain the residual and whether the
%           optimizer was successful.
%       o sampledShapes: 3D-array with length and height of the sampled
%           pulse shapes. Dimensions are nSamples x 2 x nShapes.
%       o covMats: 3D-Array containing a parameter covariance matrix for
%           each base pulse shape in the first two dimensions. Only valid
%           samples are used for computing the covariance matrix.
%       o covMatsAll: 3D-Array containing a parameter covariance matrix for
%           each base pulse shape in the first two dimensions. All samples
%           are used for computing the covariance matrix.
%       o stdParams: Matrix containing the standard deviations of the
%           parameters sampled for each base pulse shape in the rows. Only
%           valid samples are used for computing the standard deviation.
%       o stdParamsAll: Matrix containing the standard deviations of the
%           parameters sampled for each base pulse shape in the rows. All
%           samples are used for computing the standard deviation.
%       o stdRes: Vector containing the standard deviations of the
%           residuals for each base pulse shape. Only valid samples are
%           used for computing the standard deviation.
%       o stdResAll: Vector containing the standard deviations of the
%           residuals for each base pulse shape. All samples are used for
%           computing the standard deviation.
%       o meanRes: Vector containing the means of the residuals for each
%           base pulse shape. Only valid samples are used for computing the
%           means.
%       o meanResAll: Vector containing the means of the residuals for each
%           base pulse shape. All samples are used for computing the means.
%       o covDets: Vector containing the determinants of the covariance
%           matrices (see above). Only valid samples are used.
%       o covDetsAll: Vector containing the determinants of the covariance
%           matrices (see above). All samples are used.
%       o failRate: Vector containing the failure rate of the optimization
%           for each base pulse shape. Between 0 (all successful) and 1
%           (all failed).
%       o inTimes: Vector containing the length of the base pulse shape for
%           each run.
%       o optRes: Cell array with additional optimizer output (e.g., number
%           of iterations) per sample.
%
%   See also SAMPLEPARAMETERS, POSTPROCESS

% Copyright: (C) 2008-2016 The CADET Authors
%            See the license note at the end of the file.

	% Comment out the following line to avoid resetting the RNG
	rng(0, 'twister');

	% Configure Monte-Carlo sampling
	opts.nSamples = nSamples;
	opts.saveToFile = saveToFile;

	% Standard deviation of absolute and relative Gaussian noise on the
	% injected mass
	opts.stdInMassAbs = 0;
	opts.stdInMassRel = 0.05; % High noise: 0.1 = 10 %, low noise: 0.05 = 5 %

	% Standard deviation of absolute and relative Gaussian noise on the
	% pulse length
	opts.stdInTimeAbs = 1; % High noise: 5, low noise: 1
	opts.stdInTimeRel = 0;

	% Standard deviation of absolute and relative Gaussian noise on the
	% artificial measurements
	opts.measNoiseAbs = 5e-5;
	opts.measNoiseRel = 0;

	% Preallocate space
	data = cell(3, 1);
	nShapes = 18;
	
	% Loop through parameter estimation scenarios
	for opMode = 1:3

		% Number of estimated parameters
		nParams = 1;
		if opMode == 1
			nParams = 2;
		end

		stats = [];

		% Preallocate space for ...

		% Samples
		stats.samples = zeros(nSamples, nParams + 2, nShapes);
		stats.sampledShapes = zeros(nSamples, 2, nShapes);

		% Covariance matrices and their determinants
		stats.covMats = zeros(nParams, nParams, nShapes);
		stats.covMatsAll = zeros(nParams, nParams, nShapes);
		stats.covDets = zeros(nShapes,1);
		stats.covDetsAll = zeros(nShapes,1);
		% Standard deviations of parameters
		stats.stdParams = zeros(nShapes, nParams);
		stats.stdParamsAll = zeros(nShapes, nParams);
		% Standard deviations of residuals
		stats.stdRes = zeros(nShapes,1);
		stats.stdResAll = zeros(nShapes,1);
		% Means of residuals
		stats.meanRes = zeros(nShapes,1);
		stats.meanResAll = zeros(nShapes,1);
		% Determinants of covariance matrices
		stats.covDets = zeros(nShapes,1);
		stats.covDetsAll = zeros(nShapes,1);
		% Optimizer failure rate
		stats.failRate = zeros(nShapes, 1);
		% Length of base pulse in seconds
		stats.inTimes = zeros(nShapes,1);
		% Additional optimizer output
		stats.optResult = cell(nShapes, 1);

		% Loop through base pulse shapes
		for i = 1:nShapes

			[samples, sampledShape, inTime, optResult] = sampleParameters(i, opMode, opts);

			stats.samples(:, :, i) = samples;
			stats.sampledShapes(:, :, i) = sampledShape;

			stats.inTimes(i) = inTime;
			stats.optResult = optResult;

			% Determine indices of valid optimizations
			idx = samples(:, end) > 0;

			% Calculate failure rate
			stats.failRate(i) = sum(~idx) / size(samples,1);

			% Calculate covariance matrices using all or only valid samples
			stats.covMats(:, :, i) = cov(samples(idx, 1:end-2));
			stats.covMatsAll(:, :, i) = cov(samples(:, 1:end-2));

			% Calculate means of residuals using all or only valid samples
			stats.meanRes(i, :) = mean(samples(idx, end-1));
			stats.meanResAll(i, :) = mean(samples(:, end-1));

			% Calculate standard deviations of residuals and parameters 
			% using all or only valid samples
			stats.stdParams(i, :) = std(samples(idx, 1:end-2), 0, 1);
			stats.stdParamsAll(i, :) = std(samples(:, 1:end-2), 0, 1);
			stats.stdRes(i, :) = std(samples(idx, end-1));
			stats.stdResAll(i, :) = std(samples(:, end-1));

			% Compute the determinant of the estimated covariance matrices
			stats.covDets(i) = abs(det(stats.covMats(:, :, i)));
			stats.covDetsAll(i) = abs(det(stats.covMatsAll(:, :, i)));

			% Save intermediates
			if saveToFile
				save(['ExpDes-Results-Mode' num2str(opMode) '-Shape' num2str(i) '.mat'], 'stats', 'samples', 'sampledShape', 'inTime', 'opts', 'optResult');
			end
		end

		% Save intermediates
		if saveToFile
			save(['ExpDes-Results-Part' num2str(opMode) '.mat'], 'stats');
		end

		data{opMode} = stats;
	end

	if saveToFile
		save(['ExpDes-Results.mat'], 'data', 'opts');
	end
end

% =============================================================================
%  CADET - The Chromatography Analysis and Design Toolkit
%  
%  Copyright (C) 2008-2016: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
