function pulseInjectionLinearBinding()
%PULSEINJECTIONLINEARBINDING Pulse injection of two components using a linear binding model
%
% Copyright: © 2008-2015 Eric von Lieres, Joel Andersson, Andreas Püttmann, Sebastian Schnittert, Samuel Leweke
%            See the license note at the end of the file.

    model = ModelGRM();
    
    % General
    model.nComponents = 2;

    % Initial conditions
	model.initialMobileConcentration = [0.0 0.0];
    model.initialSolidConcentration = [0.0 0.0];
    
    % Adsorption
    model.kineticBindingModel = false;
    model.bindingModel = LinearBinding();
    model.bindingParameters.LIN_KA         = [1.14 0.98];
    model.bindingParameters.LIN_KD         = [0.02 0.01];
    
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
    model.nInletSections = 2;
    model.sectionTimes = [0.0 10.0 4000.0];
    model.sectionContinuity = false(1,1);
    
    model.sectionConstant       = zeros(model.nComponents, model.nInletSections);
    model.sectionLinear         = zeros(model.nComponents, model.nInletSections);
    model.sectionQuadratic      = zeros(model.nComponents, model.nInletSections);
    model.sectionCubic          = zeros(model.nComponents, model.nInletSections);

    % Sec 1
    model.sectionConstant(:,1)  = 7.14e-3;  % component 1

    % Discretization
    disc = DiscretizationGRM();
    disc.nCellsColumn = 16;
    disc.nCellsParticle = 4;
    
    % Simulator
    sim = Simulator(model, disc);
    sim.solutionTimes = linspace(0, 4000.0, 4001);
    
    % Run
    result = sim.simulate();
            
    plot(result.solution.time, result.solution.outlet);
    legend('Comp 1', 'Comp 2');
    grid on;
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
