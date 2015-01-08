function exampleF()
%EXAMPLEF One fit with two components, artifical data, and external inlet 
% profile
%
% Copyright: © 2008-2015 Eric von Lieres, Joel Andersson, Andreas Püttmann, Sebastian Schnittert, Samuel Leweke
%            See the license note at the end of the file.

    tIn = cell(2,1);
    inMeas = cell(2,1);
    
    % Inlet of component 1, can be read in from external source (e.g.,
    % Excel, CSV, DAT, etc.)
    tIn{1} = [0, 2000, 2020, 2040, 2060, 2080, 2100, 2120, 2140, 2160, 2180, 2200, 50000];
    inMeas{1} = [1, 1, 1, 0.9, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0, 0].* 6.14e-3;
    
    % Inlet of component 2
    tIn{2} = [0, 2000, 2020, 2040, 2060, 2080,  3080, 3120, 3140, 4040,   4080, 4120, 4160, 4200, 4240, 4280, 4320, 4360, 4400, 50000];
    inMeas{2} = [1, 1, 1, 0.9, 0.8, 0.6, 0.6, 0.8, 1.0, 1.0,   0.9, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0, 0].* 5.14e-3;
    
    % Param Info
    params = [{'MCL_KA'}, {'MCL_KD'}, {'MCL_KA'}, {'MCL_KD'}];
    comps = [1 1 2 2];
    secs = [-1 -1 -1 -1];

    fitData = cell(1);

    % Fit data
    fit = [];
    fit.idxComp = [1, 2];
    fit.tOut = linspace(0, 50000, 2001);
    fit.sim = createModel(fit, tIn, inMeas);
    fit.outMeas = generateArtificialData(fit);
    fit.outMeas = sum(fit.outMeas(:, fit.idxComp), 2); % Sum components
    
    fitData{1} = fit;

    % Set parameters for all simulators
    for i = 1:length(fitData)
        fitData{i}.sim.setParameters(params, comps, secs, true(length(params), 1));
        fitData{i}.task = fitData{i}.sim.prepareSimulation();
    end
    
    % Fit the data
    quietMode = false;      % Disable quiet mode
    loBound = [];
    upBound = [10, 1e-2, 10, 1e-2];
    
    initParams = [2.5, 7.2e-3, 0.2, 2e-3];  % True values: [1.14, 0.002, 1.0, 0.005]
    logScale = true(length(initParams), 1); % Enable log scaling

    [params, residual] = fitColumn(fitData, initParams, loBound, upBound, logScale, quietMode);
end

function [sim] = createModel(fit, tIn, inMeas)
    model = ModelGRM();
    
    % General
    model.nComponents = 2;

    % Initial conditions
	model.initialMobileConcentration = [0.0 0.0];
    model.initialSolidConcentration = [0.0 0.0];
    
    % Adsorption
    model.kineticBindingModel = true;
    model.bindingModel = MultiComponentLangmuirBinding();
    model.bindingParameters.MCL_KA         = [1.14 1.0];
    model.bindingParameters.MCL_KD         = [0.002 0.005];
    model.bindingParameters.MCL_QMAX       = [4.88 3.46];
    
    % Transport
    model.dispersionColumn          = 5.75e-8;
    model.filmDiffusion             = [6.9e-6 6.9e-6];
    model.diffusionParticle         = [6.07e-11 6.07e-11];
    model.diffusionParticleSurface  = [0.0 0.0];
    model.interstitialVelocity      = 5.75e-4;

    % Geometry
    model.columnLength        = 0.014;
    model.particleRadius      = 4.5e-5;
    model.porosityColumn      = 0.37;
    model.porosityParticle    = 0.75;
    
    % Inlet
    model.setInletsFromData(tIn, inMeas);

    % Discretization
    disc = DiscretizationGRM();
    disc.nCellsColumn = 16;
    disc.nCellsParticle = 4;
    
    % Simulator
    sim = Simulator(model, disc);
    sim.solutionTimes = fit.tOut;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%             The following code generates artificial data                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = generateArtificialData(fit)
    % Run
    result = fit.sim.simulate();
    data = result.solution.outlet;
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
