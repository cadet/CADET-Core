
classdef GeneralRateModel2D < Model
	%GeneralRateModel2D Represents the 2D general rate model of liquid column chromatography
	%   Holds all model specific parameters (e.g., interstitial velocity, dispersion, etc.) which are necessary
	%   for a CADET simulation. This class is mainly concerned with the transport model. Binding is handled by a
	%   binding model class.
	%
	% See also MODEL, SINGLEGRM2D, MODELSYSTEM
	
	% Copyright: (C) 2008-2020 The CADET Authors
	%            See the license note at the end of the file.

	properties
		bindingModel; % Object of the used binding model class
		reactionModelBulk; % Object of the used reaction model class for the bulk volume
		reactionModelParticle; % Object of the used reaction model class for the particle volume
	end

	properties (Dependent)
		
		% Discretization

		nCellsColumn; % Number of axial cells in the column
		nCellsRadial; % Number of radial cells in the column
		nCellsParticle; % Number of radial cells in the particle
		nBoundStates; % Number of bound states for each component
		nParticleTypes; % Number of particle types

		radialDiscretizationType; % Type of radial discretization (e.g., equivolume)
		radialCellPosition; % Positions of the cells in the radial direction (between 0 and columnRadius)

		particleDiscretizationType; % Type of particle discretization (e.g., equivolume)
		particleCellPosition; % Positions of the cells in the particle (dimensionless, between 0 and 1)

		useAnalyticJacobian; % Determines whether Jacobian is calculated analytically or via AD

		% Initial values
		initialBulk; % Initial concentrations for each component in the bulk volume in [mol / m^3_IV]
		initialParticle; % Initial concentrations for each component in the bead liquid volume in [mol / m^3_MP]
		initialSolid; % Initial concentrations for each component in the bead solid volume in [mol / m^3_SP]
		initialState; % Initial concentrations for each degree of freedom

		% Transport
		
		dispersionColumn; % Axial dispersion coefficient of the mobile phase transport inside the column in [m^2_IV / s]
		dispersionColumnMultiplexMode; % Multiplex mode of dispersionColumn
		dispersionColumnRadial; % Radial dispersion coefficient of the mobile phase transport inside the column in [m^2_IV / s]
		dispersionColumnRadialMultiplexMode; % Multiplex mode of dispersionColumnRadial
		interstitialVelocity; % Indicates direction of flow [m_IV / s]
		interstitialVelocityMultiplexMode; % Multiplex mode of interstitialVelocity

		filmDiffusion; % Film diffusion coefficient in [m / s]
		filmDiffusionMultiplexMode; % Multiplex mode of film diffusion
		diffusionParticle; % Diffusion coefficient of the mobile phase components inside the particle in [m^2_MP / s]
		diffusionParticleMultiplexMode; % Multiplex mode of particle diffusion
		diffusionParticleSurface; % Diffusion coefficient of the mobile phase components on the surface of the particle in [m^2_SP / s]
		diffusionParticleSurfaceMultiplexMode; % Multiplex mode of particle surface diffusion

		% Geometry
		
		porosityColumn; % Porosity of the column
		porosityParticle; % Porosity of the particle
		poreAccessibility; % Pore accessibility
		poreAccessibilityMultiplexMode; % Multiplex mode of pore accessibility

		columnLength; % Length of the column in [m]
		columnRadius; % Radius of the column in [m]
		particleRadius; % Radius of the particles in [m]
		particleCoreRadius; % Core radius of the particle in [m]
		particleTypeVolumeFractions; % Volume fractions of particle types
		particleTypeVolumeFractionsMultiplexMode; % Multiplex mode of particleTypeVolumeFractions

		% Numerical method for advection
		
		reconstructionType; % Type of reconstruction
		wenoBoundaryHandling; % How WENO handles left and right boundary
		wenoEpsilon; % Epsilon in the WENO method (for continuity detection)
		wenoOrder; % Order of the WENO method

		% Numerical method for solving the schur complement system (options for GMRES)

		gramSchmidtType; % Type of Gram-Schmidt orthogonalization process
		maxKrylovSize; % Maximum size of the Krylov subspace
		maxRestarts; % Maximum number of restarts in GMRES
		schurSafetyTol; % Schur-complement safety error tolerance
		linearSolverBulk; % Linear solver used for the column bulk block
	end

	properties (Dependent, Transient)
		nInletPorts; % Number of inlet ports
		nOutletPorts; % Number of outlet ports
	end
	
	properties (Constant)
		name = 'GENERAL_RATE_MODEL_2D'; % Type of the model according to CADET file format specs
		hasInlet = true; % Determines whether the unit operation has an inlet
		hasOutlet = true; % Determines whether the unit operation has an outlet
	end

	properties (Constant, Access = 'protected')
		hasConsistencySolver = true; % Determines whether this unit operation model has a consistency solver
	end

	methods
		
		function obj = GeneralRateModel2D()
			%GENERALRATEMODEL2D Constructs a GeneralRateModel2D object and inserts as much default values as possible

			obj = obj@Model();
			obj.data.discretization.weno = [];

			% Set some default values
			obj.bindingModel = [];
			obj.reactionModelBulk = [];
			obj.reactionModelParticle = [];
			obj.particleCoreRadius = 0;
			obj.particleTypeVolumeFractions = 1;
			obj.nParticleTypes = 1;

			obj.useAnalyticJacobian = true;
			obj.particleCellPosition = [];
			obj.particleDiscretizationType = 'EQUIDISTANT_PAR';

			obj.radialCellPosition = [];
			obj.radialDiscretizationType = 'EQUIDISTANT';

			obj.reconstructionType = 'WENO';
			obj.wenoBoundaryHandling = 0; % Decrease order of WENO scheme at boundary
			obj.wenoEpsilon = 1e-12;
			obj.wenoOrder = 3; % Largest possible order

			obj.gramSchmidtType = 1; % Modified Gram-Schmidt (more stable)
			obj.maxKrylovSize = 0; % Use largest possible size
			obj.maxRestarts = 0;
			obj.schurSafetyTol = 1e-8;
		end
		

		% Discretization

		function val = get.nCellsColumn(obj)
			val = double(obj.data.discretization.NCOL);
		end

		function set.nCellsColumn(obj, val)
			validateattributes(val, {'numeric'}, {'scalar', 'nonempty', '>=', 1, 'finite', 'real'}, '', 'nCellsColumn');
			obj.data.discretization.NCOL = int32(val);
			obj.hasChanged = true;
		end

		function val = get.nCellsRadial(obj)
			val = double(obj.data.discretization.NRAD);
		end

		function set.nCellsRadial(obj, val)
			validateattributes(val, {'numeric'}, {'scalar', 'nonempty', '>=', 1, 'finite', 'real'}, '', 'nCellsRadial');
			obj.data.discretization.NRAD = int32(val);
			obj.hasChanged = true;
		end

		function val = get.nCellsParticle(obj)
			val = double(obj.data.discretization.NPAR);
		end

		function set.nCellsParticle(obj, val)
			validateattributes(val, {'numeric'}, {'nonempty', 'vector', '>=', 1, 'finite', 'real'}, '', 'nCellsParticle');
			obj.data.discretization.NPAR = int32(val);
			obj.hasChanged = true;
		end

		function val = get.nBoundStates(obj)
			val = double(obj.data.discretization.NBOUND);
		end

		function set.nBoundStates(obj, val)
			validateattributes(val, {'numeric'}, {'nonnegative', 'vector', 'nonempty', 'finite', 'real'}, '', 'nBoundStates');
			obj.data.discretization.NBOUND = int32(val);
			obj.hasChanged = true;
		end

		function val = get.nParticleTypes(obj)
			val = double(obj.data.discretization.NPARTYPE);
		end

		function set.nParticleTypes(obj, val)
			validateattributes(val, {'numeric'}, {'scalar', 'nonempty', '>=', 1, 'finite', 'real'}, '', 'nParticleTypes');
			obj.data.discretization.NPARTYPE = int32(val);
			obj.hasChanged = true;
		end

		function val = get.particleDiscretizationType(obj)
			val = obj.data.discretization.PAR_DISC_TYPE;
		end

		function set.particleDiscretizationType(obj, val)
			if iscell(val)
				obj.data.discretization.PAR_DISC_TYPE = cellfun(@(s) validatestring(s, {'EQUIDISTANT_PAR', 'EQUIVOLUME_PAR', 'USER_DEFINED_PAR'}, '', 'particleDiscretizationType'), val, 'UniformOutput', false);
			else
				obj.data.discretization.PAR_DISC_TYPE = {validatestring(val, {'EQUIDISTANT_PAR', 'EQUIVOLUME_PAR', 'USER_DEFINED_PAR'}, '', 'particleDiscretizationType')};
			end
			obj.hasChanged = true;
		end

		function val = get.particleCellPosition(obj)
			if ~isfield(obj.data.discretization, 'PAR_DISC_VECTOR')
				val = [];
			else
				val = obj.data.discretization.PAR_DISC_VECTOR;
			end
		end
		
		function set.particleCellPosition(obj, val)
			if ~isempty(val)
				validateattributes(val, {'double'}, {'vector', 'increasing', '>=', 0.0, '<=', 1.0, 'finite', 'real'}, '', 'particleCellPosition');
				obj.data.discretization.PAR_DISC_VECTOR = val;
			else
				if isfield(obj.data.discretization, 'PAR_DISC_VECTOR')
					obj.data.discretization = rmfield(obj.data.discretization, 'PAR_DISC_VECTOR');
				end
			end
			obj.hasChanged = true;
		end

		function val = get.radialDiscretizationType(obj)
			val = obj.data.discretization.RADIAL_DISC_TYPE;
		end

		function set.radialDiscretizationType(obj, val)
			if iscell(val)
				obj.data.discretization.RADIAL_DISC_TYPE = cellfun(@(s) validatestring(s, {'EQUIDISTANT', 'EQUIVOLUME', 'USER_DEFINED'}, '', 'radialDiscretizationType'), val, 'UniformOutput', false);
			else
				obj.data.discretization.RADIAL_DISC_TYPE = {validatestring(val, {'EQUIDISTANT', 'EQUIVOLUME', 'USER_DEFINED'}, '', 'radialDiscretizationType')};
			end
			obj.hasChanged = true;
		end

		function val = get.radialCellPosition(obj)
			if ~isfield(obj.data.discretization, 'RADIAL_COMPARTMENTS')
				val = [];
			else
				val = obj.data.discretization.RADIAL_COMPARTMENTS;
			end
		end
		
		function set.radialCellPosition(obj, val)
			if ~isempty(val)
				validateattributes(val, {'double'}, {'vector', 'increasing', '>=', 0.0, '<=', obj.columnRadius, 'finite', 'real'}, '', 'radialCellPosition');
				obj.data.discretization.RADIAL_COMPARTMENTS = val;
			else
				if isfield(obj.data.discretization, 'RADIAL_COMPARTMENTS')
					obj.data.discretization = rmfield(obj.data.discretization, 'RADIAL_COMPARTMENTS');
				end
			end
			obj.hasChanged = true;
		end

		function val = get.useAnalyticJacobian(obj)
			val = logical(obj.data.discretization.USE_ANALYTIC_JACOBIAN);
		end

		function set.useAnalyticJacobian(obj, val)
			validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, '', 'useAnalyticJacobian');
			obj.data.discretization.USE_ANALYTIC_JACOBIAN = int32(logical(val));
			obj.hasChanged = true;
		end

		% Initial values

		function val = get.initialBulk(obj)
			val = obj.data.INIT_C;
		end

		function set.initialBulk(obj, val)
			validateattributes(val, {'double'}, {'nonnegative', 'vector', 'nonempty', 'finite', 'real'}, '', 'initialBulk');
			obj.data.INIT_C = val;
			obj.hasChanged = true;
		end

		function val = get.initialParticle(obj)
			if isfield(obj.data, 'INIT_CP')
				val = obj.data.INIT_CP;
			else
				val = [];
			end
		end

		function set.initialParticle(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'INIT_CP');
			else
				validateattributes(val, {'double'}, {'nonnegative', 'vector', 'finite', 'real'}, '', 'initialParticle');
				obj.data.INIT_CP = val;
			end
			obj.hasChanged = true;
		end

		function val = get.initialSolid(obj)
			if ~isfield(obj.data, 'INIT_Q')
				val = [];
			else
				val = obj.data.INIT_Q;
			end
		end

		function set.initialSolid(obj, val)
			if ~isempty(val)
				validateattributes(val, {'double'}, {'nonnegative', 'vector', 'finite', 'real'}, '', 'initialSolid');
			end
			obj.data.INIT_Q = val;
			obj.hasChanged = true;
		end

		function val = get.initialState(obj)
			if isfield(obj.data, 'INIT_STATE')
				val = obj.data.INIT_STATE;
			else
				val = [];
			end
		end

		function set.initialState(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'INIT_STATE');
			else
				validateattributes(val, {'double'}, {'nonnegative', 'vector', 'nonempty', 'finite', 'real'}, '', 'initialState');
				obj.data.INIT_STATE = val;
			end
			obj.hasChanged = true;
		end

		% Transport
		
		function val = get.dispersionColumn(obj)
			val = obj.data.COL_DISPERSION;
		end

		function set.dispersionColumn(obj, val)
			validateattributes(val, {'double'}, {'nonnegative', 'vector', 'nonempty', 'finite', 'real'}, '', 'dispersionColumn');
			obj.data.COL_DISPERSION = val(:);
			obj.hasChanged = true;
		end

		function val = get.dispersionColumnMultiplexMode(obj)
			if isfield(obj.data, 'COL_DISPERSION_MULTIPLEX')
				val = double(obj.data.COL_DISPERSION_MULTIPLEX);
			else
				val = [];
			end
		end

		function set.dispersionColumnMultiplexMode(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'COL_DISPERSION_MULTIPLEX');
			else
				validateattributes(val, {'double'}, {'nonnegative', 'scalar', 'nonempty', 'finite', 'real', '>=', 0, '<=', 7}, '', 'dispersionColumnMultiplexMode');
				obj.data.COL_DISPERSION_MULTIPLEX = int32(val);
			end
			obj.hasChanged = true;
		end

		function val = get.dispersionColumnRadial(obj)
			val = obj.data.COL_DISPERSION_RADIAL;
		end

		function set.dispersionColumnRadial(obj, val)
			validateattributes(val, {'double'}, {'nonnegative', 'vector', 'nonempty', 'finite', 'real'}, '', 'dispersionColumnRadial');
			obj.data.COL_DISPERSION_RADIAL = val(:);
			obj.hasChanged = true;
		end

		function val = get.dispersionColumnRadialMultiplexMode(obj)
			if isfield(obj.data, 'COL_DISPERSION_RADIAL_MULTIPLEX')
				val = double(obj.data.COL_DISPERSION_RADIAL_MULTIPLEX);
			else
				val = [];
			end
		end

		function set.dispersionColumnRadialMultiplexMode(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'COL_DISPERSION_RADIAL_MULTIPLEX');
			else
				validateattributes(val, {'double'}, {'nonnegative', 'scalar', 'nonempty', 'finite', 'real', '>=', 0, '<=', 7}, '', 'dispersionColumnRadialMultiplexMode');
				obj.data.COL_DISPERSION_RADIAL_MULTIPLEX = int32(val);
			end
			obj.hasChanged = true;
		end

		function val = get.interstitialVelocity(obj)
			if isfield(obj.data, 'VELOCITY')
				val = obj.data.VELOCITY;
			else
				val = [];
			end
		end

		function set.interstitialVelocity(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'VELOCITY');
			else
				validateattributes(val, {'double'}, {'positive', 'vector', 'nonempty', 'finite', 'real'}, '', 'interstitialVelocity');
				obj.data.VELOCITY = val(:);
			end
			obj.hasChanged = true;
		end

		function val = get.interstitialVelocityMultiplexMode(obj)
			if isfield(obj.data, 'VELOCITY_MULTIPLEX')
				val = double(obj.data.VELOCITY_MULTIPLEX);
			else
				val = [];
			end
		end

		function set.interstitialVelocityMultiplexMode(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'VELOCITY_MULTIPLEX');
			else
				validateattributes(val, {'double'}, {'nonnegative', 'scalar', 'nonempty', 'finite', 'real', '>=', 0, '<=', 3}, '', 'interstitialVelocityMultiplexMode');
				obj.data.VELOCITY_MULTIPLEX = int32(val);
			end
			obj.hasChanged = true;
		end

		function val = get.filmDiffusion(obj)
			val = obj.data.FILM_DIFFUSION;
			if (numel(val) >= obj.nComponents * obj.nParticleTypes) && (mod(numel(val), obj.nComponents * obj.nParticleTypes) == 0)
				val = permute(reshape(val, obj.nComponents, obj.nParticleTypes, numel(val) / (obj.nComponents * obj.nParticleTypes)), [3,2,1]);
			elseif (numel(val) >= obj.nComponents) && (mod(numel(val), obj.nComponents) == 0)
				val = reshape(val, obj.nComponents, numel(val) / obj.nComponents).';
			end
		end

		function set.filmDiffusion(obj, val)
			validateattributes(val, {'double'}, {'nonnegative', '3d', 'nonempty', 'finite', 'real'}, '', 'filmDiffusion');
			if length(size(val)) == 3
				val = permute(val, [3,2,1]);
			else
				val = val.';
			end
			obj.data.FILM_DIFFUSION = val(:);
			obj.hasChanged = true;
		end

		function val = get.filmDiffusionMultiplexMode(obj)
			if isfield(obj.data, 'FILM_DIFFUSION_MULTIPLEX')
				val = double(obj.data.FILM_DIFFUSION_MULTIPLEX);
			else
				val = [];
			end
		end

		function set.filmDiffusionMultiplexMode(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'FILM_DIFFUSION_MULTIPLEX');
			else
				validateattributes(val, {'double'}, {'nonnegative', 'scalar', 'nonempty', 'finite', 'real', '>=', 0, '<=', 3}, '', 'filmDiffusionMultiplexMode');
				obj.data.FILM_DIFFUSION_MULTIPLEX = int32(val);
			end
			obj.hasChanged = true;
		end

		function val = get.diffusionParticle(obj)
			val = obj.data.PAR_DIFFUSION;
			if (numel(val) >= obj.nComponents * obj.nParticleTypes) && (mod(numel(val), obj.nComponents * obj.nParticleTypes) == 0)
				val = permute(reshape(val, obj.nComponents, obj.nParticleTypes, numel(val) / (obj.nComponents * obj.nParticleTypes)), [3,2,1]);
			elseif (numel(val) >= obj.nComponents) && (mod(numel(val), obj.nComponents) == 0)
				val = reshape(val, obj.nComponents, numel(val) / obj.nComponents).';
			end
		end

		function set.diffusionParticle(obj, val)
			validateattributes(val, {'double'}, {'nonnegative', '3d', 'nonempty', 'finite', 'real'}, '', 'diffusionParticle');
			if length(size(val)) == 3
				val = permute(val, [3,2,1]);
			else
				val = val.';
			end
			obj.data.PAR_DIFFUSION = val(:);
			obj.hasChanged = true;
		end

		function val = get.diffusionParticleMultiplexMode(obj)
			if isfield(obj.data, 'PAR_DIFFUSION_MULTIPLEX')
				val = double(obj.data.PAR_DIFFUSION_MULTIPLEX);
			else
				val = [];
			end
		end

		function set.diffusionParticleMultiplexMode(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'PAR_DIFFUSION_MULTIPLEX');
			else
				validateattributes(val, {'double'}, {'nonnegative', 'scalar', 'nonempty', 'finite', 'real', '>=', 0, '<=', 3}, '', 'diffusionParticleMultiplexMode');
				obj.data.PAR_DIFFUSION_MULTIPLEX = int32(val);
			end
			obj.hasChanged = true;
		end

		function val = get.diffusionParticleSurface(obj)
			if ~isfield(obj.data, 'PAR_SURFDIFFUSION')
				val = [];
				return;
			end
			
			val = obj.data.PAR_SURFDIFFUSION;
		end

		function set.diffusionParticleSurface(obj, val)
			validateattributes(val, {'double'}, {'nonnegative', 'vector', 'finite', 'real'}, '', 'diffusionParticleSurface');
			obj.data.PAR_SURFDIFFUSION = val(:);
			obj.hasChanged = true;
		end

		function val = get.diffusionParticleSurfaceMultiplexMode(obj)
			if isfield(obj.data, 'PAR_SURFDIFFUSION_MULTIPLEX')
				val = double(obj.data.PAR_SURFDIFFUSION_MULTIPLEX);
			else
				val = [];
			end
		end

		function set.diffusionParticleSurfaceMultiplexMode(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'PAR_SURFDIFFUSION_MULTIPLEX');
			else
				validateattributes(val, {'double'}, {'nonnegative', 'scalar', 'nonempty', 'finite', 'real', '>=', 0, '<=', 3}, '', 'diffusionParticleSurfaceMultiplexMode');
				obj.data.PAR_SURFDIFFUSION_MULTIPLEX = int32(val);
			end
			obj.hasChanged = true;
		end

		% Geometry
		
		function val = get.porosityColumn(obj)
			val = obj.data.COL_POROSITY;
		end

		function set.porosityColumn(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>', 0.0, '<=', 1.0, 'finite', 'real'}, '', 'porosityColumn');
			obj.data.COL_POROSITY = val(:);
			obj.hasChanged = true;
		end

		function val = get.porosityParticle(obj)
			val = obj.data.PAR_POROSITY;
		end

		function set.porosityParticle(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>', 0.0, '<=', 1.0, 'finite', 'real'}, '', 'porosityParticle');
			obj.data.PAR_POROSITY = val;
			obj.hasChanged = true;
		end
		
		function val = get.poreAccessibility(obj)
			if isfield(obj.data, 'PORE_ACCESSIBILITY')
				val = obj.data.PORE_ACCESSIBILITY;
				if (numel(val) >= obj.nComponents * obj.nParticleTypes) && (mod(numel(val), obj.nComponents * obj.nParticleTypes) == 0)
					val = permute(reshape(val, obj.nComponents, obj.nParticleTypes, numel(val) / (obj.nComponents * obj.nParticleTypes)), [3,2,1]);
				end
			else
				val = [];
			end
		end

		function set.poreAccessibility(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'PORE_ACCESSIBILITY');
			else
				validateattributes(val, {'double'}, {'3d', 'nonempty', '>=', 0.0, '<=', 1.0, 'finite', 'real'}, '', 'poreAccessibility');
				if length(size(val)) == 2
					val = val.';
				else
					val = permute(val, [3,2,1]);
				end
				obj.data.PORE_ACCESSIBILITY = val(:);
			end
			obj.hasChanged = true;
		end

		function val = get.poreAccessibilityMultiplexMode(obj)
			if isfield(obj.data, 'PORE_ACCESSIBILITY_MULTIPLEX')
				val = double(obj.data.PORE_ACCESSIBILITY_MULTIPLEX);
			else
				val = [];
			end
		end

		function set.poreAccessibilityMultiplexMode(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'PORE_ACCESSIBILITY_MULTIPLEX');
			else
				validateattributes(val, {'double'}, {'nonnegative', 'scalar', 'nonempty', 'finite', 'real', '>=', 0, '<=', 3}, '', 'poreAccessibilityMultiplexMode');
				obj.data.PORE_ACCESSIBILITY_MULTIPLEX = int32(val);
			end
			obj.hasChanged = true;
		end

		function val = get.columnLength(obj)
			val = obj.data.COL_LENGTH;
		end

		function set.columnLength(obj, val)
			validateattributes(val, {'double'}, {'positive', 'scalar', 'nonempty', 'finite', 'real'}, '', 'columnLength');
			obj.data.COL_LENGTH = val;
			obj.hasChanged = true;
		end

		function val = get.columnRadius(obj)
			val = obj.data.COL_RADIUS;
		end

		function set.columnRadius(obj, val)
			validateattributes(val, {'double'}, {'positive', 'scalar', 'nonempty', 'finite', 'real'}, '', 'columnRadius');
			obj.data.COL_RADIUS = val;
			obj.hasChanged = true;
		end

		function val = get.particleRadius(obj)
			val = obj.data.PAR_RADIUS;
		end

		function set.particleRadius(obj, val)
			validateattributes(val, {'double'}, {'positive', 'vector', 'nonempty', 'finite', 'real'}, '', 'particleRadius');
			obj.data.PAR_RADIUS = val;
			obj.hasChanged = true;
		end

		function val = get.particleCoreRadius(obj)
			val = obj.data.PAR_CORERADIUS;
		end

		function set.particleCoreRadius(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'particleCoreRadius');		
			obj.data.PAR_CORERADIUS = val;
			obj.hasChanged = true;
		end

		function val = get.particleTypeVolumeFractions(obj)
			val = obj.data.PAR_TYPE_VOLFRAC;
		end

		function set.particleTypeVolumeFractions(obj, val)
			validateattributes(val, {'double'}, {'vector', 'nonempty', '>=', 0.0, '<=', 1.0, 'finite', 'real'}, '', 'particleTypeVolumeFractions');
			obj.data.PAR_TYPE_VOLFRAC = val;
			obj.hasChanged = true;
		end

		function val = get.particleTypeVolumeFractionsMultiplexMode(obj)
			if isfield(obj.data, 'PAR_TYPE_VOLFRAC_MULTIPLEX')
				val = double(obj.data.PAR_TYPE_VOLFRAC_MULTIPLEX);
			else
				val = [];
			end
		end

		function set.particleTypeVolumeFractionsMultiplexMode(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'PAR_TYPE_VOLFRAC_MULTIPLEX');
			else
				validateattributes(val, {'double'}, {'nonnegative', 'scalar', 'nonempty', 'finite', 'real', '>=', 0, '<=', 3}, '', 'particleTypeVolumeFractionsMultiplexMode');
				obj.data.PAR_TYPE_VOLFRAC_MULTIPLEX = int32(val);
			end
			obj.hasChanged = true;
		end

		% Numerical method for advection
		
		function val = get.reconstructionType(obj)
			val = obj.data.discretization.RECONSTRUCTION;
		end

		function set.reconstructionType(obj, val)
			obj.data.discretization.RECONSTRUCTION = validatestring(val, {'WENO'}, '', 'reconstructionType');
			obj.hasChanged = true;
		end

		function val = get.wenoBoundaryHandling(obj)
			val = double(obj.data.discretization.weno.BOUNDARY_MODEL);
		end

		function set.wenoBoundaryHandling(obj, val)
			validateattributes(val, {'numeric'}, {'>=', 0, '<=', 3, 'scalar', 'nonempty', 'finite', 'real'}, '', 'wenoBoundaryHandling');
			obj.data.discretization.weno.BOUNDARY_MODEL = int32(val);
			obj.hasChanged = true;
		end

		function val = get.wenoEpsilon(obj)
			val = obj.data.discretization.weno.WENO_EPS;
		end

		function set.wenoEpsilon(obj, val)
			validateattributes(val, {'double'}, {'positive', 'scalar', 'nonempty', 'finite', 'real'}, '', 'wenoEpsilon');
			obj.data.discretization.weno.WENO_EPS = val;
			obj.hasChanged = true;
		end

		function val = get.wenoOrder(obj)
			val = double(obj.data.discretization.weno.WENO_ORDER);
		end

		function set.wenoOrder(obj, val)
			validateattributes(val, {'numeric'}, {'>=', 1, 'scalar', 'nonempty', '<=', 3, 'finite', 'real'}, '', 'wenoOrder');
			obj.data.discretization.weno.WENO_ORDER = int32(val);
			obj.hasChanged = true;
		end

		% Numerical method for solving the schur complement system (options for GMRES)

		function val = get.gramSchmidtType(obj)
			val = double(obj.data.discretization.GS_TYPE);
		end

		function set.gramSchmidtType(obj, val)
			validateattributes(val, {'numeric'}, {'>=', 0, '<=', 1, 'scalar', 'nonempty', 'finite', 'real'}, '', 'gramSchmidtType');
			obj.data.discretization.GS_TYPE = int32(val);
			obj.hasChanged = true;
		end

		function val = get.maxKrylovSize(obj)
			val = double(obj.data.discretization.MAX_KRYLOV);
		end

		function set.maxKrylovSize(obj, val)
			validateattributes(val, {'numeric'}, {'nonnegative', 'scalar', 'nonempty', 'finite', 'real'}, '', 'maxKrylovSize');
			obj.data.discretization.MAX_KRYLOV = int32(val);
			obj.hasChanged = true;
		end

		function val = get.maxRestarts(obj)
			val = double(obj.data.discretization.MAX_RESTARTS);
		end

		function set.maxRestarts(obj, val)
			validateattributes(val, {'numeric'}, {'nonnegative', 'scalar', 'nonempty', 'finite', 'real'}, '', 'maxRestarts');
			obj.data.discretization.MAX_RESTARTS = int32(val);
			obj.hasChanged = true;
		end

		function val = get.schurSafetyTol(obj)
			val = obj.data.discretization.SCHUR_SAFETY;
		end

		function set.schurSafetyTol(obj, val)
			validateattributes(val, {'double'}, {'positive', 'scalar', 'nonempty', 'finite', 'real'}, '', 'schurSafetyTol');
			obj.data.discretization.SCHUR_SAFETY = val;
			obj.hasChanged = true;
		end

		function val = get.linearSolverBulk(obj)
			val = obj.data.discretization.LINEAR_SOLVER_BULK;
		end

		function set.linearSolverBulk(obj, val)
			obj.data.discretization.LINEAR_SOLVER_BULK = validatestring(val, {'DENSE', 'UMFPACK'}, '', 'linearSolverBulk');
			obj.hasChanged = true;
		end

		function set.bindingModel(obj, val)
			for i = 1:numel(val)
				if ~isempty(val(i)) && ~isa(val(i), 'BindingModel')
					error('CADET:invalidConfig', 'Expected a valid binding model at index %d.', i);
				end
			end
			obj.bindingModel = val;
			obj.hasChanged = true;
		end

		function set.reactionModelBulk(obj, val)
			if numel(val) > 1
				error('CADET:invalidConfig', 'Expected a single reaction model instead of array.');
			end
			if ~isempty(val) && ~isa(val, 'ReactionModel')
				error('CADET:invalidConfig', 'Expected a valid reaction model.');
			end
			obj.reactionModelBulk = val;
			obj.hasChanged = true;
		end

		function set.reactionModelParticle(obj, val)
			for i = 1:numel(val)
				if ~isempty(val(i)) && ~isa(val(i), 'ReactionModel')
					error('CADET:invalidConfig', 'Expected a valid reaction model at index %d.', i);
				end
			end
			obj.reactionModelParticle = val;
			obj.hasChanged = true;
		end

		function val = get.nInletPorts(obj)
			val = double(obj.data.discretization.NRAD);
		end

		function val = get.nOutletPorts(obj)
			val = double(obj.data.discretization.NRAD);
		end

		function S = saveobj(obj)
			S = obj.saveobj@Model();

			S.bindingModel = [];
			S.bindingModelClass = [];

			if ~isempty(obj.bindingModel)
				S.bindingModelClass = cell(numel(obj.bindingModel), 1);
				S.bindingModel = cell(numel(obj.bindingModel), 1);
				for i = 1:numel(obj.bindingModel)
					S.bindingModelClass{i} = class(obj.bindingModel(i));
					S.bindingModel{i} = obj.bindingModel(i).saveobj();
				end
			end

			S.reactionModelBulk = [];
			S.reactionModelBulkClass = [];

			if ~isempty(obj.reactionModelBulk)
				S.reactionModelBulkClass = class(obj.reactionModelBulk);
				S.reactionModelBulk = obj.reactionModelBulk.saveobj();
			end

			S.reactionModelParticle = [];
			S.reactionModelParticleClass = [];

			if ~isempty(obj.reactionModelParticle)
				S.reactionModelParticleClass = cell(numel(obj.reactionModelParticle), 1);
				S.reactionModelParticle = cell(numel(obj.reactionModelParticle), 1);
				for i = 1:numel(obj.reactionModelParticle)
					S.reactionModelParticleClass{i} = class(obj.reactionModelParticle(i));
					S.reactionModelParticle{i} = obj.reactionModelParticle(i).saveobj();
				end
			end
		end


		function res = validate(obj, sectionTimes)
			%VALIDATE Validates the configuration of the Model
			%   RES = VALIDATE(SECTIONTIMES) uses the section times of the simulator given in
			%   SECTIONTIMES to validate the Model object. Returns true in RES if everything
			%   is fine and false otherwise.
			%
			% See also MODEL.VALIDATE

			res = obj.validate@Model(sectionTimes);
			nSections = numel(sectionTimes) - 1;
			if ~isfield(obj.data, 'COL_DISPERSION')
				error('CADET:invalidConfig', 'Property dispersionColumn must be set.');
			end
			if ~isfield(obj.data, 'COL_DISPERSION_RADIAL')
				error('CADET:invalidConfig', 'Property dispersionColumn must be set.');
			end
			if ~isfield(obj.data, 'FILM_DIFFUSION')
				error('CADET:invalidConfig', 'Property filmDiffusion must be set.');
			end
			if ~isfield(obj.data, 'PAR_DIFFUSION')
				error('CADET:invalidConfig', 'Property diffusionParticle must be set.');
			end
			if ~isfield(obj.data, 'PAR_SURFDIFFUSION')
				error('CADET:invalidConfig', 'Property diffusionParticleSurface must be set.');
			end
			if ~isfield(obj.data, 'COL_POROSITY')
				error('CADET:invalidConfig', 'Property porosityColumn must be set.');
			end
			if ~isfield(obj.data, 'PAR_POROSITY')
				error('CADET:invalidConfig', 'Property porosityParticle must be set.');
			end
			if ~isfield(obj.data, 'COL_LENGTH')
				error('CADET:invalidConfig', 'Property columnLength must be set.');
			end
			if ~isfield(obj.data, 'COL_RADIUS')
				error('CADET:invalidConfig', 'Property columnLength must be set.');
			end
			if ~isfield(obj.data, 'PAR_RADIUS')
				error('CADET:invalidConfig', 'Property particleRadius must be set.');
			end


			validateattributes(obj.nCellsColumn, {'numeric'}, {'scalar', 'nonempty', '>=', 1, 'finite', 'real'}, '', 'nCellsColumn');
			validateattributes(obj.nCellsRadial, {'numeric'}, {'scalar', 'nonempty', '>=', 1, 'finite', 'real'}, '', 'nCellsRadial');
			validateattributes(obj.nCellsParticle, {'numeric'}, {'vector', 'nonempty', '>=', 1, 'finite', 'real'}, '', 'nCellsParticle');
			validateattributes(obj.nParticleTypes, {'numeric'}, {'scalar', 'nonempty', '>=', 1, 'finite', 'real'}, '', 'nParticleTypes');

			nParType = obj.nParticleTypes;
			validateattributes(obj.nBoundStates, {'numeric'}, {'nonnegative', 'vector', 'numel', obj.nComponents * nParType, 'finite', 'real'}, '', 'nBoundStates');
			
			for i = 1:length(obj.particleDiscretizationType)
				if strcmp(obj.particleDiscretizationType{i}, 'USER_DEFINED_PAR') && ~isempty(obj.particleCellPosition)
					validateattributes(obj.particleCellPosition, {'double'}, {'vector', 'increasing', '>=', 0.0, '<=', 1.0, 'numel', obj.nCellsParticle + nParType, 'finite', 'real'}, '', 'particleCellPosition');
					if (obj.particleCellPosition(1) ~= 0) || (obj.particleCellPosition(end) ~= 1)
						error('CADET:invalidConfig', 'Expected particleCellPosition to start with 0 and end with 1.');
					end
				end
			end

			if strcmp(obj.radialDiscretizationType, 'USER_DEFINED')
				validateattributes(obj.radialCellPosition, {'double'}, {'vector', 'increasing', '>=', 0.0, '<=', obj.columnRadius, 'numel', obj.nCellsRadial + 1, 'finite', 'real'}, '', 'radialCellPosition');
				if (obj.radialCellPosition(1) ~= 0) || (obj.radialCellPosition(end) ~= obj.columnRadius)
					error('CADET:invalidConfig', 'Expected radialCellPosition to start with 0 and end with %g (columnRadius).', obj.columnRadius);
				end
			end

			validateattributes(obj.initialBulk, {'double'}, {'nonnegative', 'vector', 'finite', 'real'}, '', 'initialBulk');
			if ~any(numel(obj.initialBulk) == [obj.nCellsRadial * obj.nComponents, obj.nComponents])
				error('CADET:invalidConfig', 'Expected initialBulk to be of size %d or %d', obj.nComponents, obj.nCellsRadial * obj.nComponents);
			end
			if ~isempty(obj.initialParticle)
				validateattributes(obj.initialParticle, {'double'}, {'nonnegative', 'vector', 'finite', 'real'}, '', 'initialParticle');
				if ~any(numel(obj.initialParticle) == [obj.nCellsRadial * obj.nComponents, obj.nComponents, obj.nCellsRadial * obj.nComponents * nParType, obj.nComponents * nParType])
					error('CADET:invalidConfig', 'Expected initialParticle to be of size %d, %d, %d, or %d', obj.nComponents, obj.nCellsRadial * obj.nComponents, obj.nComponents * nParType, obj.nCellsRadial * obj.nComponents * nParType);
				end
			end
			if ~isempty(obj.initialSolid)
				validateattributes(obj.initialSolid, {'double'}, {'nonnegative', 'vector', 'finite', 'real'}, '', 'initialSolid');
				if ~any(numel(obj.initialSolid) == [obj.nCellsRadial * sum(obj.nBoundStates), sum(obj.nBoundStates)])
					error('CADET:invalidConfig', 'Expected initialParticle to be of size %d, or %d', sum(obj.nBoundStates), obj.nCellsRadial * sum(obj.nBoundStates));
				end
			end
			if ~isempty(obj.initialState)
				nDof = obj.nCellsColumn * obj.nCellsRadial * obj.nComponents * (1 + nParType) + obj.nCellsColumn * obj.nCellsRadial * (obj.nCellsParticle(:).' * (obj.nComponents .* ones(nParType, 1) + kron(eye(nParType), ones(1, obj.nComponents)) * obj.nBoundStates(:)));
				if (numel(obj.initialState) ~= nDof) && (numel(obj.initialState) ~= 2 * nDof)
					error('CADET:invalidConfig', 'Expected initialState to be of size %d or %d.', nDof, 2*nDof);
				end
			end

			nTotalBnd = sum(obj.nBoundStates);
			if ~any(numel(obj.dispersionColumn) == [1, nSections, obj.nComponents, obj.nCellsRadial, nSections * obj.nCellsRadial, nSections * obj.nComponents, nSections * obj.nCellsRadial * obj.nComponents])
				error('CADET:invalidConfig', 'Expected dispersionColumn to be of size %d, %d, %d, %d, %d, %d, or %d.', ...
					1, nSections, obj.nComponents, obj.nCellsRadial, nSections * obj.nCellsRadial, nSections * obj.nComponents, nSections * obj.nCellsRadial * obj.nComponents);
			end
			if ~any(numel(obj.dispersionColumnRadial) == [1, nSections, obj.nComponents, obj.nCellsRadial, nSections * obj.nCellsRadial, nSections * obj.nComponents, nSections * obj.nCellsRadial * obj.nComponents])
				error('CADET:invalidConfig', 'Expected dispersionColumn to be of size %d, %d, %d, %d, %d, %d, or %d.', ...
					1, nSections, obj.nComponents, obj.nCellsRadial, nSections * obj.nCellsRadial, nSections * obj.nComponents, nSections * obj.nCellsRadial * obj.nComponents);
			end
			if ~isempty(obj.interstitialVelocity)
				if ~any(numel(obj.interstitialVelocity) == [1, nSections, obj.nCellsRadial, nSections * obj.nCellsRadial])
					error('CADET:invalidConfig', 'Expected interstitialVelocity to be of size %d, %d, %d, or %d.', 1, nSections, obj.nCellsRadial, nSections * obj.nCellsRadial);
				end
			end
			if (numel(obj.filmDiffusion) ~= obj.nComponents) && (numel(obj.filmDiffusion) ~= obj.nComponents * nSections) && (numel(obj.filmDiffusion) ~= obj.nComponents * nParType) && (numel(obj.filmDiffusion) ~= nSections * obj.nComponents * nParType)
				error('CADET:invalidConfig', 'Expected filmDiffusion to be of size %d, %d, %d, or %d.', obj.nComponents, obj.nComponents * nSections, obj.nComponents * nParType, nSections * obj.nComponents * nParType);
			end
			if (numel(obj.diffusionParticle) ~= obj.nComponents) && (numel(obj.diffusionParticle) ~= obj.nComponents * nSections) && (numel(obj.diffusionParticle) ~= obj.nComponents * nParType) && (numel(obj.diffusionParticle) ~= nSections * obj.nComponents * nParType)
				error('CADET:invalidConfig', 'Expected diffusionParticle to be of size %d, %d, %d, or %d.', obj.nComponents, obj.nComponents * nSections, obj.nComponents * nParType, nSections * obj.nComponents * nParType);
			end
			if (numel(obj.diffusionParticleSurface) ~= nTotalBnd) && (numel(obj.diffusionParticleSurface) ~= nSections * nTotalBnd) && (numel(obj.diffusionParticleSurface) ~= nSections * sum(obj.nBoundStates(1:obj.nComponents))) && (numel(obj.diffusionParticleSurface) ~= sum(obj.nBoundStates(1:obj.nComponents)))
				error('CADET:invalidConfig', 'Expected diffusionParticleSurface to be of size %d, %d, %d, or %d.', nTotalBnd, nTotalBnd * nSections, sum(obj.nBoundStates(1:obj.nComponents)), nSections * sum(obj.nBoundStates(1:obj.nComponents)));
			end

			validateattributes(obj.porosityColumn, {'double'}, {'vector', 'nonempty', '>', 0.0, '<=', 1.0, 'finite', 'real'}, '', 'porosityColumn');
			if (numel(obj.porosityColumn) ~= 1) && (numel(obj.porosityColumn) ~= obj.nCellsRadial)
				error('CADET:invalidConfig', 'Expected porosityColumn to be of size %d or %d (number of radial zones).', 1, obj.nCellsRadial);
			end
			validateattributes(obj.porosityParticle, {'double'}, {'vector', 'nonempty', '>', 0.0, '<=', 1.0, 'finite', 'real'}, '', 'porosityParticle');
			if (numel(obj.porosityParticle) ~= 1) && (numel(obj.porosityParticle) ~= nParType)
				error('CADET:invalidConfig', 'Expected porosityParticle to be of size %d or %d (number of particle types).', 1, nParType);
			end
			validateattributes(obj.columnLength, {'double'}, {'positive', 'scalar', 'nonempty', 'finite', 'real'}, '', 'columnLength');
			validateattributes(obj.columnRadius, {'double'}, {'positive', 'scalar', 'nonempty', 'finite', 'real'}, '', 'columnRadius');
			validateattributes(obj.particleRadius, {'double'}, {'positive', 'vector', 'nonempty', 'finite', 'real'}, '', 'particleRadius');
			if (numel(obj.particleRadius) ~= 1) && (numel(obj.particleRadius) ~= nParType)
				error('CADET:invalidConfig', 'Expected particleRadius to be of size %d or %d (number of particle types).', 1, nParType);
			end
			validateattributes(obj.particleCoreRadius, {'double'}, {'vector', 'nonempty', '>=', 0.0, 'finite', 'real'}, '', 'particleCoreRadius');
			if (numel(obj.particleCoreRadius) ~= 1) && (numel(obj.particleCoreRadius) ~= nParType)
				error('CADET:invalidConfig', 'Expected particleCoreRadius to be of size %d or %d (number of particle types).', 1, nParType);
			end

			if any(obj.particleCoreRadius >= obj.particleRadius)
				error('CADET:invalidConfig', 'Expected particleCoreRadius to be smaller than particleRadius.');
			end
			if ~isempty(obj.poreAccessibility)
				validateattributes(obj.poreAccessibility, {'double'}, {'vector', '>', 0.0, '<=', 1.0, 'finite', 'real'}, '', 'poreAccessibility');
				if (numel(obj.poreAccessibility) ~= obj.nComponents) && (numel(obj.poreAccessibility) ~= obj.nComponents * nSections) && (numel(obj.poreAccessibility) ~= obj.nComponents * nParType) && (numel(obj.poreAccessibility) ~= nSections * obj.nComponents * nParType)
					error('CADET:invalidConfig', 'Expected poreAccessibility to be of size %d, %d, %d, or %d.', obj.nComponents, obj.nComponents * nSections, obj.nComponents * nParType, nSections * obj.nComponents * nParType);
				end
			end
			validateattributes(obj.particleTypeVolumeFractions, {'double'}, {'vector', 'nonempty', '>=', 0.0, '<=', 1.0, 'finite', 'real'}, '', 'particleTypeVolumeFractions');
			if ~any(numel(obj.particleTypeVolumeFractions) == [nParType, nParType * obj.nCellsRadial, nParType * obj.nCellsColumn, nParType * obj.nCellsRadial * obj.nCellsColumn])
				error('CADET:invalidConfig', 'Expected particleTypeVolumeFractions to be of size %d, %d, %d, or %d.', nParType, nParType * obj.nCellsRadial, nParType * obj.nCellsColumn, nParType * obj.nCellsRadial * obj.nCellsColumn);
			end
			sumFracs = arrayfun(@(idx) sum(obj.particleTypeVolumeFractions(((idx-1) * nParType + 1):(idx * nParType))), 1:(numel(obj.particleTypeVolumeFractions) / nParType));
			if any(abs(sumFracs - 1.0) >= 1e-10)
				error('CADET:invalidConfig', 'Expected particleTypeVolumeFractions to sum to 1.0 for each cell.');
			end

			for i = 1:numel(obj.bindingModel)
				if ~isempty(obj.bindingModel(i)) && ~isa(obj.bindingModel(i), 'BindingModel')
					error('CADET:invalidConfig', 'Expected a valid binding model.');
				end
				if isempty(obj.bindingModel(i)) && (sum(obj.nBoundStates(((i-1) * obj.nComponents + 1):(i * obj.nComponents))) ~= 0)
					error('CADET:invalidConfig', 'Expected no bound states when using no binding model.');
				end

				if ~isempty(obj.bindingModel(i))
					res = obj.bindingModel(i).validate(obj.nComponents, obj.nBoundStates(((i-1) * obj.nComponents + 1):(i * obj.nComponents))) && res;
				end
			end

			if ~isempty(obj.reactionModelBulk)
				if ~isa(obj.reactionModelBulk, 'ReactionModel')
					error('CADET:invalidConfig', 'Expected a valid reaction model in bulk volume.');
				end

				res = obj.reactionModelBulk.validate(obj.nComponents, obj.nBoundStates(((i-1) * obj.nComponents + 1):(i * obj.nComponents))) && res;

				if ~obj.reactionModelBulk.hasBulkPhaseReactions
					warning('CADET:unexptectedConfig', 'Bulk reaction model does not contain bulk reactions.');
				end
			end

			for i = 1:numel(obj.reactionModelParticle)
				if ~isempty(obj.reactionModelParticle(i)) && ~isa(obj.reactionModelParticle(i), 'ReactionModel')
					error('CADET:invalidConfig', 'Expected a valid reaction model in particle volume.');
				end

				if ~isempty(obj.reactionModelParticle(i))
					res = obj.reactionModelParticle(i).validate(obj.nComponents, obj.nBoundStates(((i-1) * obj.nComponents + 1):(i * obj.nComponents))) && res;

					if ~obj.reactionModelParticle(i).hasLiquidPhaseReactions && ~obj.reactionModelParticle(i).hasSolidPhaseReactions
						warning('CADET:unexptectedConfig', 'Particle reaction model %d does neither contain liquid nor solid phase reactions.', i);
					end
				end
			end
		end

		function res = assembleConfig(obj)
			%ASSEMBLECONFIG Assembles the configuration according to the CADET file format spec
			%   RES = ASSEMBLECONFIG() returns a nested Matlab struct RES that represents the
			%   model as detailed in the CADET file format spec.
			%
			% See also MODEL.ASSEMBLECONFIG, MEXSIMULATOR.ASSEMBLECONFIG

			res = obj.assembleConfig@Model();

			if isempty(obj.bindingModel)
				error('CADET:invalidConfig', 'Expected valid binding model.');
			end

			res.ADSORPTION_MODEL = cell(numel(obj.bindingModel), 1);
			for i = 1:length(obj.bindingModel)
				if isempty(obj.bindingModel(i))
					res.ADSORPTION_MODEL{i} = 'NONE';
				else
					res.ADSORPTION_MODEL{i} = obj.bindingModel(i).name;
					res.(sprintf('adsorption_%03d', i-1)) = obj.bindingModel(i).assembleConfig();
				end
			end

			if isempty(obj.reactionModelBulk)
				res.REACTION_MODEL = 'NONE';
			else
				res.REACTION_MODEL = obj.reactionModelBulk.name;
				res.reaction_bulk = obj.reactionModelBulk.assembleConfig();
			end

			if isempty(obj.reactionModelParticle)
				error('CADET:invalidConfig', 'Expected valid reaction model for particle volume.');
			end

			res.REACTION_MODEL_PARTICLES = cell(numel(obj.reactionModelParticle), 1);
			for i = 1:length(obj.reactionModelParticle)
				if isempty(obj.reactionModelParticle(i))
					res.REACTION_MODEL_PARTICLES{i} = 'NONE';
				else
					res.REACTION_MODEL_PARTICLES{i} = obj.reactionModelParticle(i).name;
					res.(sprintf('reaction_particle_%03d', i-1)) = obj.reactionModelParticle(i).assembleConfig();
				end
			end
		end

		function res = assembleInitialConditions(obj)
			%ASSEMBLEINITIALCONDITIONS Assembles the initial conditions according to the CADET file format spec
			%   RES = ASSEMBLEINITIALCONDITIONS() returns a nested Matlab struct RES that represents only the
			%   initial conditions part of the model as detailed in the (full configuration) CADET file format
			%   spec.
			%
			% See also MODEL.ASSEMBLEINITIALCONDITIONS, MODELSYSTEM.ASSEMBLEINITIALCONDITIONS

			res = obj.assembleInitialConditions@Model();

			res.INIT_C = obj.data.INIT_C;
			if isfield(obj.data, 'INIT_Q')
				res.INIT_Q = obj.data.INIT_Q;
			else
				res.INIT_Q = [];
			end

			if isfield(obj.data, 'INIT_CP')
				res.INIT_CP = obj.data.INIT_CP;
			end

			if isfield(obj.data, 'INIT_STATE')
				res.INIT_STATE = obj.data.INIT_STATE;
			end
		end

		function val = getParameterValue(obj, param)
			%GETPARAMETERVALUE Retrieves a parameter value from the model
			%   VAL = GETPARAMETERVALUE(PARAM) searches for the (single) parameter identified by
			%   the struct PARAM with the fields SENS_NAME, SENS_UNIT, SENS_COMP, SENS_PARTYPE,
			%   SENS_BOUNDPHASE, SENS_REACTION, and SENS_SECTION (as returned by MAKESENSITIVITY()).
			%   The returned value VAL contains the current value of the parameter (on the Matlab side,
			%   not in the current CADET configuration) or NaN if the parameter could not be found.
			%
			% See also MODEL.GETPARAMETERVALUE, MEXSIMULATOR.GETPARAMETERVALUE,
			%   GENERALRATEMODEL.SETPARAMETERVALUE, MAKESENSITIVITY

			if param.SENS_UNIT ~= obj.unitOpIdx
				% Wrong unit operation
				val = nan;
				return;
			end

			if ~isfield(obj.data, param.SENS_NAME)
				% We don't have this parameter

				if ~isempty(obj.bindingModel)
					% Try binding model
					if (param.SENS_PARTYPE >= 0)
						val = obj.bindingModel(param.SENS_PARTYPE+1).getParameterValue(param, obj.nBoundStates((param.SENS_PARTYPE * obj.nComponents + 1):((param.SENS_PARTYPE + 1) * obj.nComponents)));
					else
						val = obj.bindingModel(1).getParameterValue(param, obj.nBoundStates(1:obj.nComponents));
					end
					return;
				end

				if ~isempty(obj.reactionModelBulk)
					% Try reaction model
					val = obj.reactionModelBulk.getParameterValue(param, obj.nBoundStates(1:obj.nComponents));
					return;
				end

				if ~isempty(obj.reactionModelParticle)
					% Try reaction model
					if (param.SENS_PARTYPE >= 0)
						val = obj.reactionModelParticle(param.SENS_PARTYPE+1).getParameterValue(param, obj.nBoundStates((param.SENS_PARTYPE * obj.nComponents + 1):((param.SENS_PARTYPE + 1) * obj.nComponents)));
					else
						val = obj.reactionModelParticle(1).getParameterValue(param, obj.nBoundStates(1:obj.nComponents));
					end
					return;
				end

				val = nan;
				return;
			end
			
			val = obj.data.(param.SENS_NAME);
			offset = 0;
			if (strcmp(param.SENS_NAME, 'PAR_SURFDIFFUSION'))
				if (param.SENS_BOUNDPHASE ~= -1)
					offset = offset + param.SENS_BOUNDPHASE;
				end
				if (param.SENS_COMP > 0)
					if (param.SENS_PARTYPE ~= -1)
						nBS = obj.nBoundStates((param.SENS_PARTYPE * obj.nComponents + 1):((param.SENS_PARTYPE + 1) * obj.nComponents));
					else
						nBS = obj.nBoundStates(1:obj.nComponents);
					end
					offset = offset + sum(nBS(1:param.SENS_COMP));
				end
				if (param.SENS_PARTYPE > 0)
					offset = offset + sum(obj.nBoundStates(1:(param.SENS_PARTYPE * obj.nComponents)));
				end
				if (param.SENS_SECTION ~= -1)
					offset = offset + param.SENS_SECTION * sum(obj.nBoundStates);
				end
			elseif (strcmp(param.SENS_NAME, 'PAR_TYPE_VOLFRAC'))
				if (param.SENS_PARTYPE ~= -1)
					offset = offset + param.SENS_PARTYPE;
				end
				if ((param.SENS_SECTION >= 0) && (param.SENS_REACTION >= 0))
					offset = offset + param.SENS_SECTION * obj.nParticleTypes * obj.nCellsRadial + param.SENS_REACTION * obj.nParticleTypes;
				else
					if (param.SENS_SECTION >= 0)
						offset = offset + param.SENS_SECTION * obj.nParticleTypes;
					elseif (param.SENS_REACTION >= 0)
						offset = offset + param.SENS_REACTION * obj.nParticleTypes;
					end
				end
			elseif (strcmp(param.SENS_NAME, 'INIT_C'))
				if (param.SENS_REACTION ~= -1)
					offset = offset + param.SENS_REACTION * obj.nComponents;
				end
				if (param.SENS_COMP ~= -1)
					offset = offset + param.SENS_COMP;
				end
			elseif (strcmp(param.SENS_NAME, 'INIT_CP'))
				if ((param.SENS_REACTION ~= -1) && (param.SENS_PARTYPE ~= -1))
					offset = offset + param.SENS_REACTION * obj.nComponents * obj.nParticleTypes + param.SENS_PARTYPE * obj.nComponents;
				elseif (param.SENS_PARTYPE ~= -1)
					offset = offset + param.SENS_PARTYPE * obj.nComponents;
				elseif (param.SENS_REACTION ~= -1)
					offset = offset + param.SENS_REACTION * obj.nComponents;
				end
				if (param.SENS_COMP ~= -1)
					offset = offset + param.SENS_COMP;
				end
			elseif (strcmp(param.SENS_NAME, 'INIT_Q'))
				if (param.SENS_REACTION ~= -1)
					if (param.SENS_PARTYPE ~= -1)
						offset = offset + param.SENS_REACTION * sum(obj.nBoundStates) + sum(obj.nBoundStates(1:((param.SENS_PARTYPE + 1) * obj.nComponents) + param.SENS_COMP));
					else
						offset = offset + param.SENS_REACTION * sum(obj.nBoundStates) + sum(obj.nBoundStates(1:param.SENS_COMP));
					end
				else
					if (param.SENS_PARTYPE ~= -1)
						offset = offset + sum(obj.nBoundStates(1:((param.SENS_PARTYPE + 1) * obj.nComponents) + param.SENS_COMP));
					else
						offset = offset + sum(obj.nBoundStates(1:param.SENS_COMP));
					end
				end
				offset = offset + param.SENS_BOUNDPHASE;
			elseif (strcmp(param.SENS_NAME, 'COL_DISPERSION') || strcmp(param.SENS_NAME, 'COL_DISPERSION_RADIAL'))
				if (param.SENS_SECTION ~= -1)
					if (param.SENS_PARTYPE ~= -1)
						if (param.SENS_COMP ~= -1)
							offset = param.SENS_SECTION * obj.nComponents * obj.nCellsRadial + param.SENS_PARTYPE * obj.nComponents + param.SENS_COMP;
						else
							offset = param.SENS_SECTION * obj.nCellsRadial + param.SENS_PARTYPE;
						end
					else
						if (param.SENS_COMP ~= -1)
							offset = param.SENS_SECTION * obj.nComponents + param.SENS_COMP;
						else
							offset = param.SENS_SECTION;
						end
					end
				else
					if (param.SENS_PARTYPE ~= -1)
						if (param.SENS_COMP ~= -1)
							offset = param.SENS_PARTYPE * obj.nComponents + param.SENS_COMP;
						else
							offset = param.SENS_PARTYPE;
						end
					else
						if (param.SENS_COMP ~= -1)
							offset = param.SENS_COMP;
						else
							offset = 0;
						end
					end
				end
			elseif (strcmp(param.SENS_NAME, 'VELOCITY'))
				if (param.SENS_SECTION ~= -1)
					if (param.SENS_PARTYPE ~= -1)
						offset = param.SENS_SECTION * obj.nCellsRadial + param.SENS_PARTYPE;
					else
						offset = param.SENS_SECTION;
					end
				else
					if (param.SENS_PARTYPE ~= -1)
						offset = param.SENS_PARTYPE;
					else
						offset = 0;
					end
				end
			else
				if (param.SENS_SECTION ~= -1)
					% Depends on section
					if (param.SENS_PARTYPE ~= -1)
						% Depends on particle type
						if (param.SENS_COMP ~= 1)
							% Depends on component
							offset = param.SENS_SECTION * obj.nComponents * length(obj.nCellsParticle) + param.SENS_PARTYPE * obj.nComponents + param.SENS_COMP;
						else
							% Independent of component
							offset = param.SENS_SECTION * length(obj.nCellsParticle) + param.SENS_PARTYPE;
						end
					else
						% Independent of particle type
						if (param.SENS_COMP ~= 1)
							% Depends on component
							offset = param.SENS_SECTION * obj.nComponents + param.SENS_COMP;
						else
							% Independent of component
							offset = param.SENS_SECTION;
						end
					end
				else
					% Independent of section
					if (param.SENS_PARTYPE ~= -1)
						% Depends on particle type
						if (param.SENS_COMP ~= 1)
							% Depends on component
							offset = param.SENS_PARTYPE * obj.nComponents + param.SENS_COMP;
						else
							% Independent of component
							offset = param.SENS_PARTYPE;
						end
					else
						% Independent of particle type
						if (param.SENS_COMP >= 0)
							offset = param.SENS_COMP;
						end
					end
				end
			end
			val = val(offset + 1);
		end

		function oldVal = setParameterValue(obj, param, newVal)
			%SETPARAMETERVALUE Sets a parameter value in the model
			%   OLDVAL = SETPARAMETERVALUE(PARAM, NEWVAL) searches for the parameter
			%   identified by the struct PARAM with the fields SENS_NAME, SENS_UNIT,
			%   SENS_COMP, SENS_PARTYPE, SENS_BOUNDPHASE, SENS_REACTION, and SENS_SECTION
			%   (as returned by MAKESENSITIVITY()). The returned value OLDVAL contains the old
			%   value of the parameter (on the Matlab side, not in the current CADET configuration)
			%   or NaN if the parameter could not be found. The value of the parameter is
			%   set to NEWVAL. The changes are not propagated to the underlying CADET simulator.
			%
			% See also MODEL.SETPARAMETERVALUE, MEXSIMULATOR.SETPARAMETERVALUE,
			%   GENERALRATEMODEL.GETPARAMETERVALUE, MAKESENSITIVITY

			if param.SENS_UNIT ~= obj.unitOpIdx
				% Wrong unit operation
				oldVal = nan;
				return;
			end

			if ~isfield(obj.data, param.SENS_NAME)
				% We don't have this parameter

				if ~isempty(obj.bindingModel)
					% Try binding model
					if (param.SENS_PARTYPE >= 0)
						oldVal = obj.bindingModel(param.SENS_PARTYPE+1).setParameterValue(param, obj.nBoundStates((param.SENS_PARTYPE * obj.nComponents + 1):((param.SENS_PARTYPE + 1) * obj.nComponents)), newVal);
					else
						oldVal = obj.bindingModel(1).setParameterValue(param, obj.nBoundStates(1:obj.nComponents), newVal);
					end
					return;
				end

				if ~isempty(obj.reactionModelBulk)
					% Try reaction model
					oldVal = obj.reactionModelBulk.setParameterValue(param, obj.nBoundStates(1:obj.nComponents), newVal);
					return;
				end

				if ~isempty(obj.reactionModelParticle)
					% Try reaction model
					if (param.SENS_PARTYPE >= 0)
						oldVal = obj.reactionModelParticle(param.SENS_PARTYPE+1).setParameterValue(param, obj.nBoundStates((param.SENS_PARTYPE * obj.nComponents + 1):((param.SENS_PARTYPE + 1) * obj.nComponents)), newVal);
					else
						oldVal = obj.reactionModelParticle(1).setParameterValue(param, obj.nBoundStates(1:obj.nComponents), newVal);
					end
					return;
				end

				oldVal = nan;
				return;
			end

			offset = 0;
			if (strcmp(param.SENS_NAME, 'PAR_SURFDIFFUSION'))
				if (param.SENS_BOUNDPHASE ~= -1)
					offset = offset + param.SENS_BOUNDPHASE;
				end
				if (param.SENS_COMP > 0)
					if (param.SENS_PARTYPE ~= -1)
						nBS = obj.nBoundStates((param.SENS_PARTYPE * obj.nComponents + 1):((param.SENS_PARTYPE + 1) * obj.nComponents));
					else
						nBS = obj.nBoundStates(1:obj.nComponents);
					end
					offset = offset + sum(nBS(1:param.SENS_COMP));
				end
				if (param.SENS_PARTYPE > 0)
					offset = offset + sum(obj.nBoundStates(1:(param.SENS_PARTYPE * obj.nComponents)));
				end
				if (param.SENS_SECTION ~= -1)
					offset = offset + param.SENS_SECTION * sum(obj.nBoundStates);
				end
			elseif (strcmp(param.SENS_NAME, 'PAR_TYPE_VOLFRAC'))
				if (param.SENS_PARTYPE ~= -1)
					offset = offset + param.SENS_PARTYPE;
				end
				if ((param.SENS_SECTION >= 0) && (param.SENS_REACTION >= 0))
					offset = offset + param.SENS_SECTION * obj.nParticleTypes * obj.nCellsRadial + param.SENS_REACTION * obj.nParticleTypes;
				else
					if (param.SENS_SECTION >= 0)
						offset = offset + param.SENS_SECTION * obj.nParticleTypes;
					elseif (param.SENS_REACTION >= 0)
						offset = offset + param.SENS_REACTION * obj.nParticleTypes;
					end
				end
			elseif (strcmp(param.SENS_NAME, 'INIT_C'))
				if (param.SENS_REACTION ~= -1)
					offset = offset + param.SENS_REACTION * obj.nComponents;
				end
				if (param.SENS_COMP ~= -1)
					offset = offset + param.SENS_COMP;
				end
			elseif (strcmp(param.SENS_NAME, 'INIT_CP'))
				if ((param.SENS_REACTION ~= -1) && (param.SENS_PARTYPE ~= -1))
					offset = offset + param.SENS_REACTION * obj.nComponents * obj.nParticleTypes + param.SENS_PARTYPE * obj.nComponents;
				elseif (param.SENS_PARTYPE ~= -1)
					offset = offset + param.SENS_PARTYPE * obj.nComponents;
				elseif (param.SENS_REACTION ~= -1)
					offset = offset + param.SENS_REACTION * obj.nComponents;
				end
				if (param.SENS_COMP ~= -1)
					offset = offset + param.SENS_COMP;
				end
			elseif (strcmp(param.SENS_NAME, 'INIT_Q'))
				if (param.SENS_REACTION ~= -1)
					if (param.SENS_PARTYPE ~= -1)
						offset = offset + param.SENS_REACTION * sum(obj.nBoundStates) + sum(obj.nBoundStates(1:((param.SENS_PARTYPE + 1) * obj.nComponents) + param.SENS_COMP));
					else
						offset = offset + param.SENS_REACTION * sum(obj.nBoundStates) + sum(obj.nBoundStates(1:param.SENS_COMP));
					end
				else
					if (param.SENS_PARTYPE ~= -1)
						offset = offset + sum(obj.nBoundStates(1:((param.SENS_PARTYPE + 1) * obj.nComponents) + param.SENS_COMP));
					else
						offset = offset + sum(obj.nBoundStates(1:param.SENS_COMP));
					end
				end
				offset = offset + param.SENS_BOUNDPHASE;
			elseif (strcmp(param.SENS_NAME, 'COL_DISPERSION') || strcmp(param.SENS_NAME, 'COL_DISPERSION_RADIAL'))
				if (param.SENS_SECTION ~= -1)
					if (param.SENS_PARTYPE ~= -1)
						if (param.SENS_COMP ~= -1)
							offset = param.SENS_SECTION * obj.nComponents * obj.nCellsRadial + param.SENS_PARTYPE * obj.nComponents + param.SENS_COMP;
						else
							offset = param.SENS_SECTION * obj.nCellsRadial + param.SENS_PARTYPE;
						end
					else
						if (param.SENS_COMP ~= -1)
							offset = param.SENS_SECTION * obj.nComponents + param.SENS_COMP;
						else
							offset = param.SENS_SECTION;
						end
					end
				else
					if (param.SENS_PARTYPE ~= -1)
						if (param.SENS_COMP ~= -1)
							offset = param.SENS_PARTYPE * obj.nComponents + param.SENS_COMP;
						else
							offset = param.SENS_PARTYPE;
						end
					else
						if (param.SENS_COMP ~= -1)
							offset = param.SENS_COMP;
						else
							offset = 0;
						end
					end
				end
			elseif (strcmp(param.SENS_NAME, 'VELOCITY'))
				if (param.SENS_SECTION ~= -1)
					if (param.SENS_PARTYPE ~= -1)
						offset = param.SENS_SECTION * obj.nCellsRadial + param.SENS_PARTYPE;
					else
						offset = param.SENS_SECTION;
					end
				else
					if (param.SENS_PARTYPE ~= -1)
						offset = param.SENS_PARTYPE;
					else
						offset = 0;
					end
				end
			else
				if (param.SENS_SECTION ~= -1)
					% Depends on section
					if (param.SENS_PARTYPE ~= -1)
						% Depends on particle type
						if (param.SENS_COMP ~= 1)
							% Depends on component
							offset = param.SENS_SECTION * obj.nComponents * length(obj.nCellsParticle) + param.SENS_PARTYPE * obj.nComponents + param.SENS_COMP;
						else
							% Independent of component
							offset = param.SENS_SECTION * length(obj.nCellsParticle) + param.SENS_PARTYPE;
						end
					else
						% Independent of particle type
						if (param.SENS_COMP ~= 1)
							% Depends on component
							offset = param.SENS_SECTION * obj.nComponents + param.SENS_COMP;
						else
							% Independent of component
							offset = param.SENS_SECTION;
						end
					end
				else
					% Independent of section
					if (param.SENS_PARTYPE ~= -1)
						% Depends on particle type
						if (param.SENS_COMP ~= 1)
							% Depends on component
							offset = param.SENS_PARTYPE * obj.nComponents + param.SENS_COMP;
						else
							% Independent of component
							offset = param.SENS_PARTYPE;
						end
					else
						% Independent of particle type
						if (param.SENS_COMP >= 0)
							offset = param.SENS_COMP;
						end
					end
				end
			end
			oldVal = obj.data.(param.SENS_NAME)(offset + 1);
			obj.data.(param.SENS_NAME)(offset + 1) = newVal;
		end

		function notifySync(obj)
			%NOTIFYSYNC Called after parameters have been sent to CADET
			%   NOTIFYSYNC() resets the HASCHANGED property.

			obj.notifySync@Model();

			if ~isempty(obj.bindingModel)
				for i = 1:length(obj.bindingModel)
					obj.bindingModel(i).notifySync();
				end
			end

			if ~isempty(obj.reactionModelBulk)
				obj.reactionModelBulk.notifySync();
			end

			if ~isempty(obj.reactionModelParticle)
				for i = 1:length(obj.reactionModelParticle)
					obj.reactionModelParticle(i).notifySync();
				end
			end
		end
	end

	methods (Access = 'protected')

		function val = getHasChanged(obj)
			val = false;
			if obj.getHasChanged@Model()
				val = true;
				return;
			end

			% Check binding models
			if ~isempty(obj.bindingModel)
				for i = 1:length(obj.bindingModel)
					if obj.bindingModel(i).hasChanged
						val = true;
						return
					end
				end
			end

			% Check reaction models
			if ~isempty(obj.reactionModelBulk)
				if obj.reactionModelBulk.hasChanged
					val = true;
					return
				end
			end

			if ~isempty(obj.reactionModelParticle)
				for i = 1:length(obj.reactionModelParticle)
					if obj.reactionModelParticle(i).hasChanged
						val = true;
						return
					end
				end
			end
		end

		function loadobjInternal(obj, S)
			obj.loadobjInternal@Model(S);

			obj.bindingModel = BindingModel.empty();
			if ~isempty(S.bindingModel)
				for i = 1:length(S.bindingModel)
					ctor = str2func([S.bindingModelClass{i} '.loadobj']);
					obj.bindingModel(i) = ctor(S.bindingModel{i});
				end
			end

			if ~isempty(S.reactionModelBulk)
				ctor = str2func([S.reactionModelBulkClass '.loadobj']);
				obj.reactionModelBulk = ctor(S.reactionModelBulk);
			else
				obj.reactionModelBulk = ReactionModel.empty();
			end

			obj.reactionModelParticle = ReactionModel.empty();
			if ~isempty(S.reactionModelParticle)
				for i = 1:length(S.reactionModelParticle)
					ctor = str2func([S.reactionModelParticleClass{i} '.loadobj']);
					obj.reactionModelParticle(i) = ctor(S.reactionModelParticle{i});
				end
			end
		end

	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = GeneralRateModel2D();
				obj.loadobjInternal(S);
			end
		end
		
	end

end

% =============================================================================
%  CADET - The Chromatography Analysis and Design Toolkit
%  
%  Copyright (C) 2008-2020: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
