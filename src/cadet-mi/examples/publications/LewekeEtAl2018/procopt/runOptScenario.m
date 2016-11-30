function [gradientShape, overlaps, optResult] = runOptScenario(isSum, initShape, absMode, elutionConstraint, colLength, algorithm, saveToFile, idxSample)
%RUNOPTSCENARIO Runs a process optimization with various options
%  
%   See the functions BILINEARGRADIENTPROCESSOVERLAPMAX and
%   BILINEARGRADIENTPROCESSOVERLAPSUM for details on the 
%   parameters and the underlying optimization scenario.
%
%   RUNOPTSCENARIO(ISSUM, INITSHAPE) runs an optimization scenario where
%   ISSUM determines whether the sum of the peak overlaps is optimized or
%   their maximum. The initial shape of the bilinear gradient is given in
%   the vector INITSHAPE that determines the start concentration, slope,
%   and length of the first gradient.
%
%   RUNOPTSCENARIO(..., ABSMODE) uses ABSMODE to choose the type of (soft)
%   abs() function that is used on the product of two concentrations.
%   Defaults to 0 (no abs() used at all).
%
%   RUNOPTSCENARIO(..., ABSMODE, ELUTIONCONSTRAINT) additionally sets the constraint
%   type that enforces complete elution. If ELUTIONCONSTRAINT = 1, the concentration
%   at the end of the profile is constrained to 1e-6 mM. If ELUTIONCONSTRAINT = 2,
%   the eluted mass has to be at least 99.5 % of the injected mass.
%
%   RUNOPTSCENARIO(..., ABSMODE, ELUTIONCONSTRAINT, COLLENGTH) additionally sets
%   the column length which defaults to 0.014m.
%
%   RUNOPTSCENARIO(..., ABSMODE, ELUTIONCONSTRAINT, COLLENGTH, ALGORITHM) additionally
%   sets the algorithm used by FMINCON (e.g., 'interior-point', 'sqp').
%
%   RUNOPTSCENARIO(..., ABSMODE, ELUTIONCONSTRAINT, COLLENGTH, ALGORITHM, SAVETOFILE)
%   controls via SAVETOFILE whether the results are saved to file (both figures
%   and data). Defaults to true. The base name of the files is chosen as
%   mult[sum | max]-abs<ABSMODE>.
%
%   RUNOPTSCENARIO(..., ABSMODE, ELUTIONCONSTRAINT, COLLENGTH, ALGORITHM, SAVETOFILE, IDXSAMPLE)
%   adds the sample index IDXSAMPLE to the file name base
%   mult[sum | max]-abs<ABSMODE>-init<IDXSAMPLE>.
%
%   [GRADIENTSHAPE, OVERLAPS, OPTRESULT] = RUNOPTSCENARIO(...) returns
%   the shape of the gradient at the optimum in GRADIENTSHAPE, the overlaps
%   of the different components in OVERLAPS and a struct with optimizer
%   information in OPTRESULT.
%
%   See also BILINEARGRADIENTPROCESSOVERLAPSUM, BILINEARGRADIENTPROCESSOVERLAPMAX

% Copyright: © 2015 Samuel Leweke, Eric von Lieres
%            See the license note at the end of the file.

	if (nargin <= 2)
		absMode = 0;
	end

	if (nargin <= 3)
		elutionConstraint = 1;
	end

	if (nargin <= 4) || isempty(colLength)
		colLength = 0.014;
	end

	if (nargin <= 5)
		algorithm = 'interior-point';
	end

	if (nargin <= 6) || isempty(saveToFile)
		saveToFile = true;
	end

	if isSum
		[gradientShape, overlaps, optResult] = bilinearGradientProcessOverlapSum(initShape, colLength, absMode, elutionConstraint, algorithm);
		namePostfix = 'sum';
	else
		[gradientShape, overlaps, optResult] = bilinearGradientProcessOverlapMax(initShape, colLength, absMode, elutionConstraint, algorithm);
		namePostfix = 'max';
	end

	if saveToFile
		% Construct file name for results
		if (nargin >= 6) && ~isempty(idxSample)
			nameBase = ['mult' namePostfix '-abs' num2str(absMode) '-init' num2str(idxSample)];
		else
			nameBase = ['mult' namePostfix '-abs' num2str(absMode)];
		end
		% Save results in MAT file
		save([nameBase '.mat'], 'initShape', 'gradientShape', 'optResult', 'overlaps', 'absMode', 'elutionConstraint', 'algorithm', 'idxSample', 'colLength', 'isSum');
		% Save figure
		saveas(gcf, [nameBase '.fig']);
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
