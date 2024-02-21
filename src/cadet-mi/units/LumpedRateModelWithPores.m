
classdef LumpedRateModelWithPores < Model
	%LumpedRateModelWithPores Represents the lumped rate model of liquid column chromatography including pores
	%   Holds all model specific parameters (e.g., interstitial velocity, dispersion, etc.) which are necessary
	%   for a CADET simulation. This class is mainly concerned with the transport model. Binding is handled by a
	%   binding model class.
	%
	% See also MODEL, SINGLELRMP, MODELSYSTEM
	
	% Copyright: (C) 2008-2024 The CADET Authors
	%            See the license note at the end of the file.

	properties
		bindingModel; % Object of the used binding model class
		reactionModelBulk; % Object of the used reaction model class for the bulk volume
		reactionModelParticle; % Object of the used reaction model class for the particle volume
	end

	properties (Dependent)
		
		% Discretization

		nCellsColumn; % Number of axial cells in the column
		nBoundStates; % Number of bound states for each component
		nParticleTypes; % Number of particle types

		useAnalyticJacobian; % Determines whether Jacobian is calculated analytically or via AD

		% Initial values
		initialBulk; % Initial concentrations for each component in the bulk volume in [mol / m^3_IV]
		initialParticle; % Initial concentrations for each component in the bead liquid volume in [mol / m^3_MP]
		initialSolid; % Initial concentrations for each component in the bead solid volume in [mol / m^3_SP]
		initialState; % Initial concentrations for each degree of freedom

		% Transport
		
		dispersionColumn; % Dispersion coefficient of the mobile phase transport inside the column in [m^2_IV / s]
		interstitialVelocity; % Interstitial velocity of the mobile phase transport inside the column in [m_IV / s]

		filmDiffusion; % Film diffusion coefficient in [m / s]
		filmDiffusionMultiplexMode; % Multiplex mode of film diffusion

		% Geometry
		
		porosityColumn; % Porosity of the column
		porosityParticle; % Porosity of the particle
		poreAccessibility; % Pore accessibility
		poreAccessibilityMultiplexMode; % Multiplex mode of pore accessibility

		columnLength; % Length of the column in [m]
		particleRadius; % Radius of the particles in [m]
		crossSectionArea; % Cross section area of the column in [m]
		particleTypeVolumeFractions; % Volume fractions of particle types
		particleGeometry; % Type of particle geometry (i.e., 'SPHERE', 'CYLINDER', 'SLAB')

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
	end

	properties(Dependent, Transient)
		nInletPorts; % Number of inlet ports
		nOutletPorts; % Number of outlet ports
	end

	properties (Constant)
		name = 'LUMPED_RATE_MODEL_WITH_PORES'; % Type of the model according to CADET file format specs
		hasInlet = true; % Determines whether the unit operation has an inlet
		hasOutlet = true; % Determines whether the unit operation has an outlet
	end

	properties (Constant, Access = 'protected')
		hasConsistencySolver = true; % Determines whether this unit operation model has a consistency solver
	end

	methods
		
		function obj = LumpedRateModelWithPores()
			%LUMPEDRATEMODELWITHPORES Constructs a LumpedRateModelWithPores object and inserts as much default values as possible

			obj = obj@Model();
			obj.data.discretization.weno = [];

			% Set some default values
			obj.bindingModel = [];
			obj.reactionModelBulk = [];
			obj.reactionModelParticle = [];

			obj.useAnalyticJacobian = true;
			obj.particleTypeVolumeFractions = 1;
			obj.nParticleTypes = 1;
			obj.particleGeometry = 'SPHERE';

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

		function val = get.particleGeometry(obj)
			val = obj.data.discretization.PAR_GEOM;
		end

		function set.particleGeometry(obj, val)
			if iscell(val)
				obj.data.discretization.PAR_GEOM = cellfun(@(s) validatestring(s, {'SPHERE', 'CYLINDER', 'SLAB'}, '', 'particleGeometry'), val, 'UniformOutput', false);
			else
				obj.data.discretization.PAR_GEOM = {validatestring(val, {'SPHERE', 'CYLINDER', 'SLAB'}, '', 'particleGeometry')};
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

		function val = get.filmDiffusion(obj)
			val = obj.data.FILM_DIFFUSION;
			if (numel(val) >= obj.nComponents * numel(obj.porosityParticle)) && (mod(numel(val), obj.nComponents * numel(obj.porosityParticle)) == 0)
				val = permute(reshape(val, obj.nComponents, numel(obj.porosityParticle), numel(val) / (obj.nComponents * numel(obj.porosityParticle))), [3,2,1]);
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

		% Geometry
		
		function val = get.porosityColumn(obj)
			val = obj.data.COL_POROSITY;
		end

		function set.porosityColumn(obj, val)
			validateattributes(val, {'double'}, {'scalar', 'nonempty', '>', 0.0, '<=', 1.0, 'finite', 'real'}, '', 'porosityColumn');
			obj.data.COL_POROSITY = val;
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
				if (numel(val) >= obj.nComponents * numel(obj.porosityParticle)) && (mod(numel(val), obj.nComponents * numel(obj.porosityParticle)) == 0)
					val = permute(reshape(val, obj.nComponents, numel(obj.porosityParticle), numel(val) / (obj.nComponents * numel(obj.porosityParticle))), [3,2,1]);
				elseif (numel(val) >= obj.nComponents) && (mod(numel(val), obj.nComponents) == 0)
					val = reshape(val, obj.nComponents, numel(val) / obj.nComponents).';
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
				if length(size(val)) == 3
					val = permute(val, [3,2,1]);
				else
					val = val.';
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

		function val = get.particleRadius(obj)
			val = obj.data.PAR_RADIUS;
		end

		function set.particleRadius(obj, val)
			validateattributes(val, {'double'}, {'positive', 'vector', 'nonempty', 'finite', 'real'}, '', 'particleRadius');
			obj.data.PAR_RADIUS = val;
			obj.hasChanged = true;
		end

		function val = get.crossSectionArea(obj)
			if isfield(obj.data, 'CROSS_SECTION_AREA')
				val = obj.data.CROSS_SECTION_AREA;
			else
				val = [];
			end
		end

		function set.crossSectionArea(obj, val)
			if isempty(val)
				obj.data = rmfield(obj.data, 'CROSS_SECTION_AREA');
			else
				validateattributes(val, {'double'}, {'positive', 'scalar', 'nonempty', 'finite', 'real'}, '', 'crossSectionArea');
				obj.data.CROSS_SECTION_AREA = val;
			end
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
			val = 1;
		end

		function val = get.nOutletPorts(obj)
			val = 1;
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
			if (~isfield(obj.data, 'VELOCITY')) && (~isfield(obj.data, 'CROSS_SECTION_AREA'))
				error('CADET:invalidConfig', 'Property interstitialVelocity or crossSectionArea must be set.');
			end
			if ~isfield(obj.data, 'FILM_DIFFUSION')
				error('CADET:invalidConfig', 'Property filmDiffusion must be set.');
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
			if ~isfield(obj.data, 'PAR_RADIUS')
				error('CADET:invalidConfig', 'Property particleRadius must be set.');
			end

			nParType = obj.nParticleTypes;
			validateattributes(obj.nCellsColumn, {'numeric'}, {'scalar', 'nonempty', '>=', 1, 'finite', 'real'}, '', 'nCellsColumn');
			validateattributes(obj.nBoundStates, {'numeric'}, {'nonnegative', 'vector', 'numel', obj.nComponents * nParType, 'finite', 'real'}, '', 'nBoundStates');
			validateattributes(obj.nParticleTypes, {'numeric'}, {'scalar', 'nonempty', '>=', 1, 'finite', 'real'}, '', 'nParticleTypes');

			nTotalBnd = sum(obj.nBoundStates);

			validateattributes(obj.initialBulk, {'double'}, {'nonnegative', 'vector', 'numel', obj.nComponents, 'finite', 'real'}, '', 'initialBulk');
			if ~isempty(obj.initialParticle)
				validateattributes(obj.initialParticle, {'double'}, {'nonnegative', 'vector', 'finite', 'real'}, '', 'initialParticle');
				if ~any(numel(obj.initialParticle) == [obj.nComponents, obj.nComponents * nParType])
					error('CADET:invalidConfig', 'Expected initialParticle to be of size %d, or %d', obj.nComponents, obj.nComponents * nParType);
				end
			end
			if ~isempty(obj.initialSolid)
				validateattributes(obj.initialSolid, {'double'}, {'nonnegative', 'vector', 'numel', nTotalBnd, 'finite', 'real'}, '', 'initialSolid');
			end
			if ~isempty(obj.initialState)
				nDof = obj.nCellsColumn * obj.nComponents + obj.nCellsColumn * (obj.nComponents * nParType + nTotalBnd) + obj.nCellsColumn * obj.nComponents * nParType;
				if (numel(obj.initialState) ~= nDof) && (numel(obj.initialState) ~= 2 * nDof)
					error('CADET:invalidConfig', 'Expected initialState to be of size %d or %d.', nDof, 2*nDof);
				end
			end

			if (numel(obj.dispersionColumn) ~= 1) && (numel(obj.dispersionColumn) ~= nSections) && (numel(obj.dispersionColumn) ~= nComponents) && (numel(obj.dispersionColumn) ~= nComponents * nSections)
				error('CADET:invalidConfig', 'Expected dispersionColumn to be of size %d, %d, %d, or %d.', 1, nComponents, nSections, nComponents * nSections);
			end
			if ~isempty(obj.interstitialVelocity)
				if (numel(obj.interstitialVelocity) ~= 1) && (numel(obj.interstitialVelocity) ~= nSections)
					error('CADET:invalidConfig', 'Expected interstitialVelocity to be of size %d or %d (number of time sections).', 1, nSections);
				end
			end
			if (numel(obj.filmDiffusion) ~= obj.nComponents) && (numel(obj.filmDiffusion) ~= obj.nComponents * nSections) && (numel(obj.filmDiffusion) ~= obj.nComponents * nParType) && (numel(obj.filmDiffusion) ~= nSections * obj.nComponents * nParType)
				error('CADET:invalidConfig', 'Expected filmDiffusion to be of size %d, %d, %d, or %d.', obj.nComponents, obj.nComponents * nSections, obj.nComponents * nParType, nSections * obj.nComponents * nParType);
			end

			validateattributes(obj.porosityColumn, {'double'}, {'scalar', 'nonempty', '>', 0.0, '<=', 1.0, 'finite', 'real'}, '', 'porosityColumn');
			validateattributes(obj.porosityParticle, {'double'}, {'vector', 'nonempty', '>', 0.0, '<=', 1.0, 'finite', 'real'}, '', 'porosityParticle');
			if (numel(obj.porosityParticle) ~= 1) && (numel(obj.porosityParticle) ~= nParType)
				error('CADET:invalidConfig', 'Expected porosityParticle to be of size %d or %d (number of particle types).', 1, nParType);
			end
			validateattributes(obj.columnLength, {'double'}, {'positive', 'scalar', 'nonempty', 'finite', 'real'}, '', 'columnLength');
			validateattributes(obj.particleRadius, {'double'}, {'positive', 'vector', 'nonempty', 'finite', 'real'}, '', 'particleRadius');
			if (numel(obj.particleRadius) ~= 1) && (numel(obj.particleRadius) ~= nParType)
				error('CADET:invalidConfig', 'Expected particleRadius to be of size %d or %d (number of particle types).', 1, nParType);
			end
			if ~isempty(obj.poreAccessibility)
				validateattributes(obj.poreAccessibility, {'double'}, {'vector', '>', 0.0, '<=', 1.0, 'finite', 'real'}, '', 'poreAccessibility');
				if (numel(obj.poreAccessibility) ~= obj.nComponents) && (numel(obj.poreAccessibility) ~= obj.nComponents * nSections) && (numel(obj.poreAccessibility) ~= obj.nComponents * nParType) && (numel(obj.poreAccessibility) ~= nSections * obj.nComponents * nParType)
					error('CADET:invalidConfig', 'Expected poreAccessibility to be of size %d, %d, %d, or %d.', obj.nComponents, obj.nComponents * nSections, obj.nComponents * nParType, nSections * obj.nComponents * nParType);
				end
			end
			if (numel(obj.particleGeometry) ~= 1) && (numel(obj.particleGeometry) ~= nParType)
				error('CADET:invalidConfig', 'Expected particleGeometry to be of size %d or %d (number of particle types).', 1, nParType);
			end
			validateattributes(obj.particleTypeVolumeFractions, {'double'}, {'vector', 'nonempty', '>=', 0.0, '<=', 1.0, 'finite', 'real'}, '', 'particleTypeVolumeFractions');
			if abs(sum(obj.particleTypeVolumeFractions) - 1.0) >= 1e-10
				error('CADET:invalidConfig', 'Expected particleTypeVolumeFractions to sum to 1.0.');
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
				res.REACTION_MODEL_PARTICLES = 'NONE';
			else
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
			%   LUMPEDRATEMODELWITHPORES.SETPARAMETERVALUE, MAKESENSITIVITY

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
			if (strcmp(param.SENS_NAME, 'INIT_Q'))
				if (param.SENS_PARTYPE ~= -1)
					offset = offset + sum(obj.nBoundStates(1:((param.SENS_PARTYPE + 1) * obj.nComponents) + param.SENS_COMP));
				else
					offset = offset + sum(obj.nBoundStates(1:param.SENS_COMP));
				end
				offset = offset + param.SENS_BOUNDPHASE;
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
			%   LUMPEDRATEMODELWITHPORES.GETPARAMETERVALUE, MAKESENSITIVITY

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
			if (strcmp(param.SENS_NAME, 'INIT_Q'))
				if (param.SENS_PARTYPE ~= -1)
					offset = offset + sum(obj.nBoundStates(1:((param.SENS_PARTYPE + 1) * obj.nComponents) + param.SENS_COMP));
				else
					offset = offset + sum(obj.nBoundStates(1:param.SENS_COMP));
				end
				offset = offset + param.SENS_BOUNDPHASE;
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
				obj = LumpedRateModelWithPores();
				obj.loadobjInternal(S);
			end
		end
		
	end

end

% =============================================================================
%  CADET
%  
%  Copyright (C) 2008-2024: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
