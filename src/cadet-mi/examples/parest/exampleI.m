function exampleI()
%EXAMPLEI Three fits with artificial data and different linked parameters.
%
% Linked parameters:
%
%   Parameters    Model
% ----------------------
%     mcl_ka      1,2,3
%     mcl_kd       1,3
%    velocity      1,2
% col_dispersion   2,3
% ----------------------
% Total: 4 parameters
%
% Copyright: © 2008-2015 Eric von Lieres, Joel Andersson, Andreas Püttmann, Sebastian Schnittert, Samuel Leweke
%            See the license note at the end of the file.

    fitData = cell(3,1);
    
    % First fit with flow rate 7e-4 and dispersion 5e-8
    
    % Param Info
    params = [{'mcl_ka'}, {'mcl_kd'}, {'velocity'}];
    comps = [1 1 -1];
    secs = [-1 -1 -1];

    fit = [];
    fit.links = [{[2,3]}, {[3]}, {[2]}];    % Link mcl_ka to fit 2, 3
                                            % Link mcl_kd to fit 3
                                            % Link velocity to fit 2

    fit.idxComp = [1];
    fit.tOut = linspace(0, 10000, 1001);
    fit.sim = createModel(fit.tOut, 1);
    fit.outMeas = generateArtificialData(fit);

    fit.sim.setParameters(params, comps, secs, true(length(params), 1));
    fit.task = fit.sim.prepareSimulation();
    fitData{1} = fit;
    
    % Second fit with same model as first fit but kd = 0.02

    % Param Info
    params = [{'mcl_ka'}, {'velocity'}, {'col_dispersion'}];
    comps = [1 -1 -1];
    secs = [-1 -1 -1];

    fit = [];
    fit.links = [{[]}, {[]}, {[3]}];        % Link col_dispersion to fit 3
    fit.idxComp = [1];
    fit.tOut = linspace(0, 10000, 1001);
    fit.sim = createModel(fit.tOut, 2);
    fit.outMeas = generateArtificialData(fit);

    fit.sim.setParameters(params, comps, secs, true(length(params), 1));
    fit.task = fit.sim.prepareSimulation();
    fitData{2} = fit;

    % Third fit with same model as first fit but flow rate 5e-4

    % Param Info
    params = [{'mcl_ka'}, {'mcl_kd'}, {'col_dispersion'}];
    comps = [1 1 -1];
    secs = [-1 -1 -1];

    fit = [];
    fit.idxComp = [1];
    fit.tOut = linspace(0, 10000, 1001);
    fit.sim = createModel(fit.tOut, 3);
    fit.outMeas = generateArtificialData(fit);

    fit.sim.setParameters(params, comps, secs, true(length(params), 1));
    fit.task = fit.sim.prepareSimulation();
    fitData{3} = fit;
    
    % Fit the data
    quietMode = false;      % Disable quiet mode
    loBound = []; %[1.05, 1e-5, 1e-6, 1e-10];
    upBound = [];
    
    % Parameters in order of first appearance:
    % mcl_ka, mcl_kd, velocity, col_dispersion
    initParams = [2.5, 5e-3, 5e-3, 1e-7];  % True values: [1.14, 0.002, 7e-4, 5.75e-8]
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
    
    if k == 2
        model.bindingParameters.MCL_KD     = [0.02];
    end
    
    % Transport
    model.dispersionColumn          = 5.75e-8;
    model.filmDiffusion             = [6.9e-6];
    model.diffusionParticle         = [6.07e-11];
    model.diffusionParticleSurface  = [0.0];
    model.interstitialVelocity      = 7e-4;

    if k == 1
        model.dispersionColumn      = 5e-8;
    end
    
    if k == 3
        model.interstitialVelocity  = 5e-4;
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
