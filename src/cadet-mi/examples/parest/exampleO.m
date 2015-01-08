function exampleO()
%EXAMPLEO One fit with artificial data. The peaks are weighted by their mass.
%
% Copyright: © 2008-2015 Eric von Lieres, Joel Andersson, Andreas Püttmann, Sebastian Schnittert, Samuel Leweke
%            See the license note at the end of the file.

    % Param Info
    params = [{'SMA_KA'}, {'SMA_KA'}, {'SMA_KA'}];
    comps = [2 3 4];
    secs = [-1 -1 -1];
    
    fitData = cell(1);
    
    % Gradient elution with observed sum signal
    fit = [];

    fit.idxComp = [2 3 4];  % Sum of all three components (ignore salt)
    fit.tOut = linspace(0, 1500, 1001);
    fit.sim = createModel(fit.tOut);
    fit.outMeas = generateArtificialData(fit);
    fit.outMeas = fit.outMeas(:, fit.idxComp); % Take components, ignore salt
    
    % Compute the mass of each peak and use it to weight them
    fit.weightComponent = 1./[trapz(fit.tOut, fit.outMeas(:,1)), trapz(fit.tOut, fit.outMeas(:,2)), trapz(fit.tOut, fit.outMeas(:,3))];

    fitData{1} = fit;
   
    % Set parameters for all simulators
    for i = 1:length(fitData)
        fitData{i}.sim.setParameters(params, comps, secs, true(length(params), 1));
        fitData{i}.task = fitData{i}.sim.prepareSimulation();
    end
    
    % Fit the data
    quietMode = false;      % Disable quiet mode
    loBound = [];
    upBound = [];
    
    initParams = [20.5 4.9 1.8];  % True values: [35.5 1.59 7.7]
    logScale = true(length(initParams), 1); % Enable log scaling

    [params, residual] = fitColumn(fitData, initParams, loBound, upBound, logScale, quietMode);
end

function [sim] = createModel(tOut)
    model = ModelGRM();
    
    % General
    model.nComponents = 4;

    % Initial conditions
	model.initialMobileConcentration = [50.0 0.0 0.0 0.0];
    model.initialSolidConcentration = [1.2e3 0.0 0.0 0.0];
    
    % Adsorption
    model.kineticBindingModel = false;
    model.bindingModel = StericMassActionBinding();
    model.bindingParameters.SMA_LAMBDA     = 1.2e3;
    model.bindingParameters.SMA_KA         = [0.0 35.5 1.59 7.7];
    model.bindingParameters.SMA_KD         = [0.0 1000 1000 1000];
    model.bindingParameters.SMA_NU         = [0.0 4.7 5.29 3.7];
    model.bindingParameters.SMA_SIGMA      = [0.0 11.83 10.6 10.0];
    
    % Transport
    model.dispersionColumn          = 5.75e-8;
    model.filmDiffusion             = [6.9e-6 6.9e-6 6.9e-6 6.9e-6];
    model.diffusionParticle         = [7e-10 6.07e-11 6.07e-11 6.07e-11];
    model.diffusionParticleSurface  = [0.0 0.0 0.0 0.0];
    model.interstitialVelocity      = 5.75e-4;

    % Geometry
    model.columnLength        = 0.014;
    model.particleRadius      = 4.5e-5;
    model.porosityColumn      = 0.37;
    model.porosityParticle    = 0.75;
    
    % Inlet
    model.nInletSections = 3;
    model.sectionTimes = [0.0 10.0 90.0 1500.0];
    model.sectionContinuity = false(2,1);
    
    model.sectionConstant       = zeros(model.nComponents, model.nInletSections);
    model.sectionLinear         = zeros(model.nComponents, model.nInletSections);
    model.sectionQuadratic      = zeros(model.nComponents, model.nInletSections);
    model.sectionCubic          = zeros(model.nComponents, model.nInletSections);

    % Sec 1
    model.sectionConstant(1,1)  = 50.0;  % component 1
    model.sectionConstant(2,1)  = 0.95;  % component 2
    model.sectionConstant(3,1)  = 1.1;   % component 3
    model.sectionConstant(4,1)  = 1.0;   % component 4

    % Sec 2
    model.sectionConstant(1,2)  = 50.0;  % component 1

    % Sec 3
    model.sectionConstant(1,3)  = 100;  % component 1
    model.sectionLinear  (1,3)  = 0.2;  % component 1

    % Discretization
    disc = DiscretizationGRM();
    disc.nCellsColumn = 16;
    disc.nCellsParticle = 4;
    
    % Simulator
    sim = Simulator(model, disc);
    sim.solutionTimes = tOut;
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
