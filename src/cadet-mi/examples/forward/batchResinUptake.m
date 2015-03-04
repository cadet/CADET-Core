function batchResinUptake()
%BATCHRESINUPTAKE Simple batch resin uptake simulation
%
% The column is degraded to a tank by setting interstitial velocity and 
% axial dispersion to 0.
%
% Since there is no inlet condition any more, one inlet section spanning
% the whole computational domain is enough.
%
% The underlying batch resin model is given by the following equations:
%
% dc_i / dt = -1 / beta_c * 3 / R_p * k_fi * [ c_i - c_pi( . , R_p)
% dc_pi / dt + 1 / beta_p * dq_i / dt = D_pi * [ d^2 c_pi / dr^2 + 2 / r * dc_pi / dr] + D_si * [ d^2 q_i / dr^2 + 2 / r * dq_i / dr]
% dq_i / dt = f_isotherm( c_p, q)
%
% Boundary conditions:
% k_fi * [ c_i - c_pi( . , R_p) = eps_p * D_pi * dc_pi / dr ( . , R_p)
%                         + (1 - eps_p) * D_si * dq_i / dr ( . , R_p)
% dc_pi / dr ( . , 0) = 0
%
% Copyright: © 2008-2015 Eric von Lieres, Joel Andersson, Andreas Püttmann, Sebastian Schnittert, Samuel Leweke
%            See the license note at the end of the file.

    simTime = 50; % Simulated time in sec

    model = ModelGRM();
    
    % General
    model.nComponents = 3;

    % Initial conditions
	model.initialMobileConcentration = [7.0 7.0 7.0]; % Initial concentration in the tank
	model.initialBeadConcentration = [0.0 0.0 0.0];   % Initial concentration in the liquid phase of the beads
    model.initialSolidConcentration = [0.0 0.0 0.0];  % Initial concentration in the solid phase of the beads
    
    % Adsorption
    model.kineticBindingModel = true;
    model.bindingModel = MultiComponentLangmuirBinding();
    model.bindingParameters.MCL_KA         = [1.14 1.6 1.3];
    model.bindingParameters.MCL_KD         = [0.002 0.002 0.002];
    model.bindingParameters.MCL_QMAX       = [6.88 8 7];
    
    % Transport
    model.filmDiffusion             = [6.9e-6 6.3e-6 6.3e-6];
    model.diffusionParticle         = [6.07e-11 5.8e-11 6.07e-11];
    model.diffusionParticleSurface  = [0.0 0.0 0.0];
    model.interstitialVelocity      = 0; % No transport in column
    model.dispersionColumn          = 0; % No transport in column

    % Geometry
    model.columnLength        = 1; % Irrelevant
    model.particleRadius      = 4.5e-5;
    model.porosityColumn      = 0.37;
    model.porosityParticle    = 0.75;
    
    % Inlet
    model.nInletSections = 1;
    model.sectionTimes = [0.0 simTime];
    model.sectionContinuity = [];

    % No feed required
    model.sectionConstant       = zeros(model.nComponents, model.nInletSections);
    model.sectionLinear         = zeros(model.nComponents, model.nInletSections);
    model.sectionQuadratic      = zeros(model.nComponents, model.nInletSections);
    model.sectionCubic          = zeros(model.nComponents, model.nInletSections);

    % Discretization
    disc = DiscretizationGRM();
    disc.nCellsParticle = 100; % This is the only relevant discretization parameter
    disc.nCellsColumn = 3;     % We require at least 3 cells because of the underlying implementation
    disc.wenoOrder = 1;        % We don't need WENO scheme (degrade to upwind)
    
    % Simulator
    sim = Simulator(model, disc);
    sim.solutionTimes = linspace(0, simTime, 2001);
    sim.solverOptions.time_integrator.INIT_STEP_SIZE = 1e-10;
    sim.solverOptions.time_integrator.ABSTOL = 1e-12;
    sim.solverOptions.time_integrator.RELTOL = 1e-8;
    sim.solverOptions.WRITE_SOLUTION_ALL = true;
    
    % Run
    result = sim.simulate();
    result = correctOrdering(result);

    figure;
    plot(result.solution.time, result.solution.outlet);
    legend('Comp 1', 'Comp 2', 'Comp 3');
    grid on;
    title('Tank concentration');
    xlabel('Time [s]');
    ylabel('Conc. [mM]');

    % Centers of radial cells in equidistant configuration (default)
    rad = linspace(model.particleRadius, 0, disc.nCellsParticle+1);
    rad = 0.5 * (rad(1:end-1) + rad(2:end));

    for idxComp = 1:model.nComponents
        figure('Name', ['Component ' num2str(idxComp)]);
        
        subplot(1,2,1);
        surf(result.solution.time, rad, squeeze(result.solution.particle(idxComp,1,:,1,:)), 'EdgeColor', 'none');
        title(['Liquid bead conc.']);
        xlabel('Time [s]');
        ylabel('Bead radial pos. [m]');
        zlabel('Conc. [mM]');

        subplot(1,2,2);
        surf(result.solution.time, rad, squeeze(result.solution.particle(idxComp,2,:,1,:)), 'EdgeColor', 'none');
        title(['Solid bead conc.']);
        xlabel('Time [s]');
        ylabel('Bead radial pos. [m]');
        zlabel('Conc. [mM]');
    end
end

function result = correctOrdering(result)
% Convert nD-arrays / tensors from row-major to Matlab's column-major
% storage format

    if ~isempty(result.solution.particle)
        temp = result.solution.particle(:);

        if ndims(result.solution.particle) == 4
            nTime = size(result.solution.particle,1);
            nCol = size(result.solution.particle,2);
            nPar = size(result.solution.particle,3);
            nPhase = size(result.solution.particle,4);
            nComp = 1;
        else
            nTime = size(result.solution.particle,1);
            nCol = size(result.solution.particle,2);
            nPar = size(result.solution.particle,3);
            nPhase = size(result.solution.particle,4);
            nComp = size(result.solution.particle,5);
        end
        
        result.solution.particle = zeros([nComp,nPhase,nPar,nCol,nTime]);
        step = nComp * nPhase * nPar * nCol;
        for i = 1:nComp
            for j = 1:nCol
                for k = 1:nPar
                    for l = 1:nPhase
                        idxOffset = (j-1) * nPhase * nPar * nComp + (k-1) * nPhase * nComp + (l-1)*nComp + (i-1);
                        result.solution.particle(i,l,k,j,:) = temp(idxOffset+[1:step:length(temp)]);
                    end
                end
            end
        end
    end
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
