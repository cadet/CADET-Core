function resumedSimulation()
%RESUMEDSIMULATION Simulates for some time and then resumes the simulation
%
% Copyright: © 2008-2015 Eric von Lieres, Joel Andersson, Andreas Püttmann, Sebastian Schnittert, Samuel Leweke
%            See the license note at the end of the file.

    % Param Info
    params = [{'MCL_KA'}, {'MCL_KD'}];
    comps = [1 1];
    secs = [-1 -1];

    % First run
    % =========
    % Create model
    sim = createModel(2);
    sim.setParameters(params, comps, secs, true(length(params), 1));
    
    % Enable output of last state
    sim.solverOptions.WRITE_SOLUTION_LAST = true;
    sim.solverOptions.WRITE_SENS_LAST = true;
    
    % Simulate
    resultPart1 = sim.simulate();
    
    % Second run
    % ==========
    sim = createModel(3);
    sim.setParameters(params, comps, secs, true(length(params), 1));
    
    % Use last state as initial state for next run
    sim.model.initialState = resultPart1.solution.lastState;
    sim.model.initialSensitivities = resultPart1.sensitivity.lastState;

    % Simulate
    resultPart2 = sim.simulate();
    
    % Collect and merge solutions
    totalTime = [resultPart1.solution.time(1:end-1); resultPart2.solution.time + resultPart1.solution.time(end)];
    totalConc = [resultPart1.solution.outlet(1:end-1); resultPart2.solution.outlet];
    % ... and sensitivities
    totalJac = [squeeze(resultPart1.sensitivity.jacobian(1:end-1,:)); squeeze(resultPart2.sensitivity.jacobian)];
    
    % Complete run
    % ============
    sim = createModel(1);
    sim.setParameters(params, comps, secs, true(length(params), 1));

    % Simulate
    resultComplete = sim.simulate();
    
    % Plots
    subplot(2,2,1);
    plot(totalTime, [totalConc, resultComplete.solution.outlet]);
    legend('Piecewise', 'Complete');
    grid on;
    title('Solution');

    subplot(2,2,2);
    plot(totalTime, [totalJac(:,1), resultComplete.sensitivity.jacobian(:,1)]);
    legend('Piecewise', 'Complete');
    grid on;
    title('Sens KA');

    subplot(2,2,3);
    plot(totalTime, [totalJac(:,2), resultComplete.sensitivity.jacobian(:,2)]);
    legend('Piecewise', 'Complete');
    grid on;
    title('Sens QMAX');
    
    subplot(2,2,4);
    semilogy(totalTime, [abs(totalConc - resultComplete.solution.outlet), abs(totalJac - resultComplete.sensitivity.jacobian)]);
    legend('Solution', 'Sens KA', 'Sens QMAX');
    grid on;
    title('Absolute error');
end

function sim = createModel(mode)
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
    
    if mode == 1
        % Full simulation
        % Inlet
        model.nInletSections = 2;
        model.sectionTimes = [0.0 1000 4000];
        model.sectionContinuity = [false];

        solTimes = linspace(0, 4000, 4001);

        model.sectionConstant       = zeros(model.nComponents, model.nInletSections);
        model.sectionLinear         = zeros(model.nComponents, model.nInletSections);
        model.sectionQuadratic      = zeros(model.nComponents, model.nInletSections);
        model.sectionCubic          = zeros(model.nComponents, model.nInletSections);

        % Sec 1
        model.sectionConstant(1,1)  = 1.0;  % component 1
    elseif mode == 2
        % First part
        % Inlet
        model.nInletSections = 2;
        model.sectionTimes = [0.0 1000 1053];
        model.sectionContinuity = [false];

        solTimes = linspace(0, 1053, 1054);

        model.sectionConstant       = zeros(model.nComponents, model.nInletSections);
        model.sectionLinear         = zeros(model.nComponents, model.nInletSections);
        model.sectionQuadratic      = zeros(model.nComponents, model.nInletSections);
        model.sectionCubic          = zeros(model.nComponents, model.nInletSections);

        % Sec 1
        model.sectionConstant(1,1)  = 1.0;  % component 1
    else
        % Second part
        % Inlet
        model.nInletSections = 1;
        model.sectionTimes = [0.0 2947];
        model.sectionContinuity = [];

        solTimes = linspace(0, 2947, 2948);

        model.sectionConstant       = zeros(model.nComponents, model.nInletSections);
        model.sectionLinear         = zeros(model.nComponents, model.nInletSections);
        model.sectionQuadratic      = zeros(model.nComponents, model.nInletSections);
        model.sectionCubic          = zeros(model.nComponents, model.nInletSections);
    end
    
    % Discretization
    disc = DiscretizationGRM();
    disc.nCellsColumn = 16;
    disc.nCellsParticle = 4;
    
    % Simulator
    sim = Simulator(model, disc);
    sim.solutionTimes = solTimes;
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
