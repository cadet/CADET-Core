function parallelSimulations()
%PARALLELSIMULATIONS Does multiple simulations of the same model with different parameters in parallel
%
% The Parallel Computing Toolbox is employed to compute many forward
% simulations in parallel.
%
% Copyright: © 2008-2015 Eric von Lieres, Joel Andersson, Andreas Püttmann, Sebastian Schnittert, Samuel Leweke
%            See the license note at the end of the file.

    [model, disc] = createModelDiscretization();

    % Simulator
    sim = Simulator(model, disc);
    sim.solutionTimes = linspace(0, 10000, 1001);
    
    % Add parameters (without sensitivity computation)
    sim.addParameter('VELOCITY', -1, -1, false);
    sim.addParameter('CONST_COEFF', 1, 1, false);
    
    % Prepare for later execution
    task = sim.prepareSimulation();
    
    % Run the simulation 100 times in parallel
    result = cell(100,1);
    parfor i = 1:length(result)
        % Draw random values for velocity and inlet concentration
        % Since rand() produces uniform random values in [0,1], the
        % expression -1 + 2*rand() transforms it to [-1, 1].
        randVel = 5.75e-4 + (-1 + 2*rand())*1e-4;
        randInlet = 7.14e-3 + (-1 + 2*rand())*1e-3;
        result{i} = sim.runWithParameters(task, [randVel, randInlet]);
    end
    
    % Plot all curves in one figure
    hold on;
    for i = 1:length(result)
        plot(result{i}.solution.time, result{i}.solution.outlet);
    end
    hold off;
    grid on;

end

function [model, disc] = createModelDiscretization()

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
