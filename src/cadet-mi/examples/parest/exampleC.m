function exampleC()
%EXAMPLEC Two fits with artificial data but different flow rates.
%
% Copyright: © 2008-2015 Eric von Lieres, Joel Andersson, Andreas Püttmann, Sebastian Schnittert, Samuel Leweke
%            See the license note at the end of the file.

    % Param Info
    params = [{'MCL_KA'}, {'VELOCITY'}];
    comps = [1 -1];
    secs = [-1 -1];

    fitData = cell(2,1);
    
    % First fit with flow rate 5.75e-4
    fit = [];
    fit.links = [{[2]}, {[]}];    % Link first parameters to second fit
                                  % velocity is independent
    fit.idxComp = [1];
    fit.tOut = linspace(0, 10000, 1001);
    fit.sim = createModel(fit.tOut, 1);
    fit.outMeas = generateArtificialData(fit);
    
    fitData{1} = fit;
    
    % Second fit with same model but flow rate 4.25e-4
    fit = [];
    fit.links = [{[1]}, {[]}];  % Link first parameters to first fit,
                                % unnecessary in this case since
                                % one link (regardless of its position) 
                                % is sufficient

    fit.idxComp = [1];
    fit.tOut = linspace(0, 10000, 1001);
    fit.sim = createModel(fit.tOut, 2);
    fit.outMeas = generateArtificialData(fit);

    fitData{2} = fit;
    
    % Set parameters for all simulators
    for i = 1:length(fitData)
        fitData{i}.sim.setParameters(params, comps, secs, true(length(params), 1));
        fitData{i}.task = fitData{i}.sim.prepareSimulation();
    end

    % Fit the data
    quietMode = false;      % Disable quiet mode
    loBound = [];
    upBound = [];
    
    % Parameters in order of first appearance:
    % mcl_ka, velocity of model 1, velocity of model 2
    initParams = [2.5, 5e-4, 5e-4];  % True values: [1.14, X]
    logScale = true(length(initParams), 1); % Enable log scaling

    [params, residual] = fitColumn(fitData, initParams, loBound, upBound, logScale, quietMode);
end

function [sim] = createModel(tOut, k)
    model = ModelGRM();
    
    % General
    model.nComponents = 1;

    % Initial conditions
	model.initialMobileConcentration = [0.0];
    model.initialSolidConcentration = [0.0];
    
    % Adsorption
    model.kineticBindingModel = true;
    model.bindingModel = MultiComponentLangmuirBinding();
    model.bindingParameters.MCL_KA         = [1.14];
    model.bindingParameters.MCL_KD         = [0.002];
    model.bindingParameters.MCL_QMAX       = [4.88];
    
    % Transport
    model.dispersionColumn          = 5.75e-8;
    model.filmDiffusion             = [6.9e-6];
    model.diffusionParticle         = [6.07e-11];
    model.diffusionParticleSurface  = [0.0];
    model.interstitialVelocity      = 5.75e-4;

    if k == 2
        model.interstitialVelocity  = 4.25e-4;
    end

    % Geometry
    model.columnLength        = 0.014;
    model.particleRadius      = 4.5e-5;
    model.porosityColumn      = 0.37;
    model.porosityParticle    = 0.75;
    
    % Inlet
    model.nInletSections = 1;
    model.sectionTimes = [0.0 10000];
    model.sectionContinuity = [];
    
    model.sectionConstant       = zeros(model.nComponents, model.nInletSections);
    model.sectionLinear         = zeros(model.nComponents, model.nInletSections);
    model.sectionQuadratic      = zeros(model.nComponents, model.nInletSections);
    model.sectionCubic          = zeros(model.nComponents, model.nInletSections);

    % Sec 1
    model.sectionConstant(1,1)  = 7.14e-3;  % component 1

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
