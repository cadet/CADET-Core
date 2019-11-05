
classdef MassActionLawReaction < ReactionModel
	%MassActionLawReaction Reaction model based on mass action law
	%
	% See also REACTIONMODEL

	% Copyright: (C) 2008-2019 The CADET Authors
	%            See the license note at the end of the file.

	properties(Constant)
		name = 'MASS_ACTION_LAW'; % Name of the reaction model according to CADET file format specs
	end

	properties (Dependent, Transient)
		% Forward reaction rates in bulk phase
		kFwdBulk;
		MAL_KFWD_BULK;
		% Backward reaction rates in bulk phase
		kBwdBulk;
		MAL_KBWD_BULK;
		% Forward reaction rates in particle liquid phase
		kFwdLiquid;
		MAL_KFWD_LIQUID;
		% Backward reaction rates in particle liquid phase
		kBwdLiquid;
		MAL_KBWD_LIQUID;
		% Forward reaction rates in particle solid phase
		kFwdSolid;
		MAL_KFWD_SOLID;
		% Backward reaction rates in particle solid phase
		kBwdSolid;
		MAL_KBWD_SOLID;
		% Stoichiometric matrix for bulk phase reactions
		Sbulk;
		MAL_STOICHIOMETRY_BULK;
		% Stoichiometric matrix for particle liquid phase reactions
		Sliquid;
		MAL_STOICHIOMETRY_LIQUID;
		% Stoichiometric matrix for particle solid phase reactions
		Ssolid;
		MAL_STOICHIOMETRY_SOLID;
		% Exponent matrix for bulk phase reactions in forward direction
		expBulkFwd;
		MAL_EXPONENTS_BULK_FWD;
		% Exponent matrix for bulk phase reactions in backward direction
		expBulkBwd;
		MAL_EXPONENTS_BULK_BWD;
		% Exponent matrix for particle liquid phase reactions in forward direction
		expLiquidFwd;
		MAL_EXPONENTS_LIQUID_FWD;
		% Exponent matrix for particle liquid phase reactions in backward direction
		expLiquidBwd;
		MAL_EXPONENTS_LIQUID_BWD;
		% Exponent matrix for modifier particle solid phase components in particle liquid phase forward reactions
		expLiquidFwdModSolid;
		MAL_EXPONENTS_LIQUID_FWD_MODSOLID;
		% Exponent matrix for modifier particle solid phase components in particle liquid phase backward reactions
		expLiquidBwdModSolid;
		MAL_EXPONENTS_LIQUID_BWD_MODSOLID;
		% Exponent matrix for particle solid phase reactions in forward direction
		expSolidFwd;
		MAL_EXPONENTS_SOLID_FWD;
		% Exponent matrix for particle solid phase reactions in backward direction
		expSolidBwd;
		MAL_EXPONENTS_SOLID_BWD;
		% Exponent matrix for modifier particle liquid phase components in particle solid phase forward reactions
		expSolidFwdModLiquid;
		MAL_EXPONENTS_SOLID_FWD_MODLIQUID;
		% Exponent matrix for modifier particle liquid phase components in particle solid phase backward reactions
		expSolidBwdModLiquid;
		MAL_EXPONENTS_SOLID_BWD_MODLIQUID;

		hasBulkPhaseReactions; % Determines whether bulk phase reactions are present in the model
		hasLiquidPhaseReactions; % Determines whether particle liquid phase reactions are present in the model
		hasSolidPhaseReactions; % Determines whether particle solid phase reactions are present in the model
	end

	methods

		function obj = MassActionLawReaction()
			%MASSACTIONLAWREACTION Constructs a MassActionLawReaction object
			%   OBJ = MASSACTIONLAWREACTION() creates an empty MassActionLawReaction model

			obj = obj@ReactionModel();
		end

		function res = validate(obj, nComponents, nBoundStates)
			%VALIDATE Validates the reaction model parameters
			%   RES = VALIDATE(NCOMPONENTS, NBOUNDSTATES) uses the number of
			%   components NCOMPONENTS of the model and the number of bound
			%   states of each of these components as given by the vector
			%   NBOUNDSTATES to validate the parameters of the reaction model.
			%   Returns true in RES if everything is fine and false otherwise.

			if (isfield(obj.data, 'MAL_STOICHIOMETRY_BULK') && ~isempty(obj.data.MAL_STOICHIOMETRY_BULK))
				nReactions = size(obj.data.MAL_STOICHIOMETRY_BULK, 2);
				validateattributes(obj.data.MAL_STOICHIOMETRY_BULK, {'double'}, {'2d', 'size', [nComponents, nReactions], 'finite', 'real'}, '', 'Sbulk');
				validateattributes(obj.data.MAL_KFWD_BULK, {'double'}, {'vector', 'numel', nReactions, '>=', 0.0, 'finite', 'real'}, '', 'kFwdBulk');
				validateattributes(obj.data.MAL_KBWD_BULK, {'double'}, {'vector', 'numel', nReactions, '>=', 0.0, 'finite', 'real'}, '', 'kBwdBulk');

				if (isfield(obj.data, 'MAL_EXPONENTS_BULK_FWD') && ~isempty(obj.data.MAL_EXPONENTS_BULK_FWD))
					validateattributes(obj.data.MAL_EXPONENTS_BULK_FWD, {'double'}, {'2d', 'size', [nComponents, nReactions], 'finite', 'real'}, '', 'expBulkFwd');
				end

				if (isfield(obj.data, 'MAL_EXPONENTS_BULK_BWD') && ~isempty(obj.data.MAL_EXPONENTS_BULK_BWD))
					validateattributes(obj.data.MAL_EXPONENTS_BULK_BWD, {'double'}, {'2d', 'size', [nComponents, nReactions], 'finite', 'real'}, '', 'expBulkBwd');
				end
			end

			if (isfield(obj.data, 'MAL_STOICHIOMETRY_LIQUID') && ~isempty(obj.data.MAL_STOICHIOMETRY_LIQUID))
				nReactions = size(obj.data.MAL_STOICHIOMETRY_LIQUID, 2);
				validateattributes(obj.data.MAL_STOICHIOMETRY_LIQUID, {'double'}, {'2d', 'size', [nComponents, nReactions], 'finite', 'real'}, '', 'Sliquid');
				validateattributes(obj.data.MAL_KFWD_LIQUID, {'double'}, {'vector', 'numel', nReactions, '>=', 0.0, 'finite', 'real'}, '', 'kFwdLiquid');
				validateattributes(obj.data.MAL_KBWD_LIQUID, {'double'}, {'vector', 'numel', nReactions, '>=', 0.0, 'finite', 'real'}, '', 'kBwdLiquid');

				if (isfield(obj.data, 'MAL_EXPONENTS_LIQUID_FWD') && ~isempty(obj.data.MAL_EXPONENTS_LIQUID_FWD))
					validateattributes(obj.data.MAL_EXPONENTS_LIQUID_FWD, {'double'}, {'2d', 'size', [nComponents, nReactions], 'finite', 'real'}, '', 'expLiquidFwd');
				end

				if (isfield(obj.data, 'MAL_EXPONENTS_LIQUID_BWD') && ~isempty(obj.data.MAL_EXPONENTS_LIQUID_BWD))
					validateattributes(obj.data.MAL_EXPONENTS_LIQUID_BWD, {'double'}, {'2d', 'size', [nComponents, nReactions], 'finite', 'real'}, '', 'expLiquidBwd');
				end

				if (isfield(obj.data, 'MAL_EXPONENTS_LIQUID_FWD_MODSOLID') && ~isempty(obj.data.MAL_EXPONENTS_LIQUID_FWD_MODSOLID))
					validateattributes(obj.data.MAL_EXPONENTS_LIQUID_FWD_MODSOLID, {'double'}, {'2d', 'size', [sum(nBoundStates), nReactions], 'finite', 'real'}, '', 'expLiquidFwdModSolid');
				end

				if (isfield(obj.data, 'MAL_EXPONENTS_LIQUID_BWD_MODSOLID') && ~isempty(obj.data.MAL_EXPONENTS_LIQUID_BWD_MODSOLID))
					validateattributes(obj.data.MAL_EXPONENTS_LIQUID_BWD_MODSOLID, {'double'}, {'2d', 'size', [sum(nBoundStates), nReactions], 'finite', 'real'}, '', 'expLiquidBwdModSolid');
				end
			end

			if (isfield(obj.data, 'MAL_STOICHIOMETRY_SOLID') && ~isempty(obj.data.MAL_STOICHIOMETRY_SOLID))
				nReactions = size(obj.data.MAL_STOICHIOMETRY_SOLID, 2);
				validateattributes(obj.data.MAL_STOICHIOMETRY_SOLID, {'double'}, {'2d', 'size', [sum(nBoundStates), nReactions], 'finite', 'real'}, '', 'Ssolid');
				validateattributes(obj.data.MAL_KFWD_SOLID, {'double'}, {'vector', 'numel', nReactions, '>=', 0.0, 'finite', 'real'}, '', 'kFwdSolid');
				validateattributes(obj.data.MAL_KBWD_SOLID, {'double'}, {'vector', 'numel', nReactions, '>=', 0.0, 'finite', 'real'}, '', 'kBwdSolid');

				if (isfield(obj.data, 'MAL_EXPONENTS_SOLID_FWD') && ~isempty(obj.data.MAL_EXPONENTS_SOLID_FWD))
					validateattributes(obj.data.MAL_EXPONENTS_SOLID_FWD, {'double'}, {'2d', 'size', [sum(nBoundStates), nReactions], 'finite', 'real'}, '', 'expSolidFwd');
				end

				if (isfield(obj.data, 'MAL_EXPONENTS_SOLID_BWD') && ~isempty(obj.data.MAL_EXPONENTS_SOLID_BWD))
					validateattributes(obj.data.MAL_EXPONENTS_SOLID_BWD, {'double'}, {'2d', 'size', [sum(nBoundStates), nReactions], 'finite', 'real'}, '', 'expSolidBwd');
				end

				if (isfield(obj.data, 'MAL_EXPONENTS_SOLID_FWD_MODLIQUID') && ~isempty(obj.data.MAL_EXPONENTS_SOLID_FWD_MODLIQUID))
					validateattributes(obj.data.MAL_EXPONENTS_SOLID_FWD_MODLIQUID, {'double'}, {'2d', 'size', [nComponents, nReactions], 'finite', 'real'}, '', 'expSolidFwdModLiquid');
				end

				if (isfield(obj.data, 'MAL_EXPONENTS_SOLID_BWD_MODLIQUID') && ~isempty(obj.data.MAL_EXPONENTS_SOLID_BWD_MODLIQUID))
					validateattributes(obj.data.MAL_EXPONENTS_SOLID_BWD_MODLIQUID, {'double'}, {'2d', 'size', [nComponents, nReactions], 'finite', 'real'}, '', 'expSolidBwdModLiquid');
				end
			end

			res = obj.validate@ReactionModel(nComponents, nBoundStates);
		end

		function res = assembleConfig(obj)
			%ASSEMBLECONFIG Assembles the configuration according to the CADET file format spec
			%   RES = ASSEMBLECONFIG() returns a Matlab struct RES that represents the
			%   reaction model as detailed in the CADET file format spec.
			%
			%   Reaction models are supposed to store their configuration (conforming to the
			%   CADET file format spec) in the data field of this base class. However, by
			%   overwriting this function, derived classes can customize the assembly process.
			%
			% See also MODEL.ASSEMBLECONFIG

			res = obj.assembleConfig@ReactionModel();

			% Convert matrices to row-major ordering
			if isfield(obj.data, 'MAL_STOICHIOMETRY_BULK') && ~isempty(obj.data.MAL_STOICHIOMETRY_BULK)
				t = obj.data.MAL_STOICHIOMETRY_BULK.';
				res.MAL_STOICHIOMETRY_BULK = t(:);
			end

			if isfield(obj.data, 'MAL_STOICHIOMETRY_LIQUID') && ~isempty(obj.data.MAL_STOICHIOMETRY_LIQUID)
				t = obj.data.MAL_STOICHIOMETRY_LIQUID.';
				res.MAL_STOICHIOMETRY_LIQUID = t(:);
			end

			if isfield(obj.data, 'MAL_STOICHIOMETRY_SOLID') && ~isempty(obj.data.MAL_STOICHIOMETRY_SOLID)
				t = obj.data.MAL_STOICHIOMETRY_SOLID.';
				res.MAL_STOICHIOMETRY_SOLID = t(:);
			end

			if isfield(obj.data, 'MAL_EXPONENTS_BULK_FWD') && ~isempty(obj.data.MAL_EXPONENTS_BULK_FWD)
				t = obj.data.MAL_EXPONENTS_BULK_FWD.';
				res.MAL_EXPONENTS_BULK_FWD = t(:);
			end

			if isfield(obj.data, 'MAL_EXPONENTS_BULK_BWD') && ~isempty(obj.data.MAL_EXPONENTS_BULK_BWD)
				t = obj.data.MAL_EXPONENTS_BULK_BWD.';
				res.MAL_EXPONENTS_BULK_BWD = t(:);
			end

			if isfield(obj.data, 'MAL_EXPONENTS_LIQUID_FWD') && ~isempty(obj.data.MAL_EXPONENTS_LIQUID_FWD)
				t = obj.data.MAL_EXPONENTS_LIQUID_FWD.';
				res.MAL_EXPONENTS_LIQUID_FWD = t(:);
			end

			if isfield(obj.data, 'MAL_EXPONENTS_LIQUID_BWD') && ~isempty(obj.data.MAL_EXPONENTS_LIQUID_BWD)
				t = obj.data.MAL_EXPONENTS_LIQUID_BWD.';
				res.MAL_EXPONENTS_LIQUID_BWD = t(:);
			end

			if isfield(obj.data, 'MAL_EXPONENTS_LIQUID_FWD_MODSOLID') && ~isempty(obj.data.MAL_EXPONENTS_LIQUID_FWD_MODSOLID)
				t = obj.data.MAL_EXPONENTS_LIQUID_FWD_MODSOLID.';
				res.MAL_EXPONENTS_LIQUID_FWD_MODSOLID = t(:);
			end

			if isfield(obj.data, 'MAL_EXPONENTS_LIQUID_BWD_MODSOLID') && ~isempty(obj.data.MAL_EXPONENTS_LIQUID_BWD_MODSOLID)
				t = obj.data.MAL_EXPONENTS_LIQUID_BWD_MODSOLID.';
				res.MAL_EXPONENTS_LIQUID_BWD_MODSOLID = t(:);
			end

			if isfield(obj.data, 'MAL_EXPONENTS_SOLID_FWD') && ~isempty(obj.data.MAL_EXPONENTS_SOLID_FWD)
				t = obj.data.MAL_EXPONENTS_SOLID_FWD.';
				res.MAL_EXPONENTS_SOLID_FWD = t(:);
			end

			if isfield(obj.data, 'MAL_EXPONENTS_SOLID_BWD') && ~isempty(obj.data.MAL_EXPONENTS_SOLID_BWD)
				t = obj.data.MAL_EXPONENTS_SOLID_BWD.';
				res.MAL_EXPONENTS_SOLID_BWD = t(:);
			end

			if isfield(obj.data, 'MAL_EXPONENTS_SOLID_FWD_MODLIQUID') && ~isempty(obj.data.MAL_EXPONENTS_SOLID_FWD_MODLIQUID)
				t = obj.data.MAL_EXPONENTS_SOLID_FWD_MODLIQUID.';
				res.MAL_EXPONENTS_SOLID_FWD_MODLIQUID = t(:);
			end

			if isfield(obj.data, 'MAL_EXPONENTS_SOLID_BWD_MODLIQUID') && ~isempty(obj.data.MAL_EXPONENTS_SOLID_BWD_MODLIQUID)
				t = obj.data.MAL_EXPONENTS_SOLID_BWD_MODLIQUID.';
				res.MAL_EXPONENTS_SOLID_BWD_MODLIQUID = t(:);
			end
		end


		function val = get.kFwdBulk(obj)
			val = obj.data.MAL_KFWD_BULK;
		end

		function set.kFwdBulk(obj, val)
			validateattributes(val, {'double'}, {'vector', '>=', 0.0, 'finite', 'real'}, '', 'kFwdBulk');
			obj.data.MAL_KFWD_BULK = val;
			obj.hasChanged = true;
		end

		function val = get.kBwdBulk(obj)
			val = obj.data.MAL_KBWD_BULK;
		end

		function set.kBwdBulk(obj, val)
			validateattributes(val, {'double'}, {'vector', '>=', 0.0, 'finite', 'real'}, '', 'kBwdBulk');
			obj.data.MAL_KBWD_BULK = val;
			obj.hasChanged = true;
		end

		function val = get.kFwdLiquid(obj)
			val = obj.data.MAL_KFWD_LIQUID;
		end

		function set.kFwdLiquid(obj, val)
			validateattributes(val, {'double'}, {'vector', '>=', 0.0, 'finite', 'real'}, '', 'kFwdLiquid');
			obj.data.MAL_KFWD_LIQUID = val;
			obj.hasChanged = true;
		end

		function val = get.kBwdLiquid(obj)
			val = obj.data.MAL_KBWD_LIQUID;
		end

		function set.kBwdLiquid(obj, val)
			validateattributes(val, {'double'}, {'vector', '>=', 0.0, 'finite', 'real'}, '', 'kBwdLiquid');
			obj.data.MAL_KBWD_LIQUID = val;
			obj.hasChanged = true;
		end

		function val = get.kFwdSolid(obj)
			val = obj.data.MAL_KFWD_SOLID;
		end

		function set.kFwdSolid(obj, val)
			validateattributes(val, {'double'}, {'vector', '>=', 0.0, 'finite', 'real'}, '', 'kFwdSolid');
			obj.data.MAL_KFWD_SOLID = val;
			obj.hasChanged = true;
		end

		function val = get.kBwdSolid(obj)
			val = obj.data.MAL_KBWD_SOLID;
		end

		function set.kBwdSolid(obj, val)
			validateattributes(val, {'double'}, {'vector', '>=', 0.0, 'finite', 'real'}, '', 'kBwdSolid');
			obj.data.MAL_KBWD_SOLID = val;
			obj.hasChanged = true;
		end

		function val = get.Sbulk(obj)
			val = obj.data.MAL_STOICHIOMETRY_BULK;
		end

		function set.Sbulk(obj, val)
			obj.data.MAL_STOICHIOMETRY_BULK = val;
			obj.hasChanged = true;
		end

		function val = get.Sliquid(obj)
			val = obj.data.MAL_STOICHIOMETRY_LIQUID;
		end

		function set.Sliquid(obj, val)
			obj.data.MAL_STOICHIOMETRY_LIQUID = val;
			obj.hasChanged = true;
		end

		function val = get.Ssolid(obj)
			val = obj.data.MAL_STOICHIOMETRY_SOLID;
		end

		function set.Ssolid(obj, val)
			obj.data.MAL_STOICHIOMETRY_SOLID = val;
			obj.hasChanged = true;
		end

		function val = get.expBulkFwd(obj)
			val = obj.data.MAL_EXPONENTS_BULK_FWD;
		end

		function set.expBulkFwd(obj, val)
			obj.data.MAL_EXPONENTS_BULK_FWD = val;
			obj.hasChanged = true;
		end

		function val = get.expBulkBwd(obj)
			val = obj.data.MAL_EXPONENTS_BULK_BWD;
		end

		function set.expBulkBwd(obj, val)
			obj.data.MAL_EXPONENTS_BULK_BWD = val;
			obj.hasChanged = true;
		end

		function val = get.expLiquidFwd(obj)
			val = obj.data.MAL_EXPONENTS_LIQUID_FWD;
		end

		function set.expLiquidFwd(obj, val)
			obj.data.MAL_EXPONENTS_LIQUID_FWD = val;
			obj.hasChanged = true;
		end

		function val = get.expLiquidBwd(obj)
			val = obj.data.MAL_EXPONENTS_LIQUID_BWD;
		end

		function set.expLiquidBwd(obj, val)
			obj.data.MAL_EXPONENTS_LIQUID_BWD = val;
			obj.hasChanged = true;
		end

		function val = get.expLiquidFwdModSolid(obj)
			val = obj.data.MAL_EXPONENTS_LIQUID_FWD_MODSOLID;
		end

		function set.expLiquidFwdModSolid(obj, val)
			obj.data.MAL_EXPONENTS_LIQUID_FWD_MODSOLID = val;
			obj.hasChanged = true;
		end

		function val = get.expLiquidBwdModSolid(obj)
			val = obj.data.MAL_EXPONENTS_LIQUID_BWD_MODSOLID;
		end

		function set.expLiquidBwdModSolid(obj, val)
			obj.data.MAL_EXPONENTS_LIQUID_BWD_MODSOLI = val;
			obj.hasChanged = true;
		end

		function val = get.expSolidFwd(obj)
			val = obj.data.MAL_EXPONENTS_SOLID_FWD;
		end

		function set.expSolidFwd(obj, val)
			obj.data.MAL_EXPONENTS_SOLID_FWD = val;
			obj.hasChanged = true;
		end

		function val = get.expSolidBwd(obj)
			val = obj.data.MAL_EXPONENTS_SOLID_BWD;
		end

		function set.expSolidBwd(obj, val)
			obj.data.MAL_EXPONENTS_SOLID_BWD = val;
			obj.hasChanged = true;
		end

		function val = get.expSolidFwdModLiquid(obj)
			val = obj.data.MAL_EXPONENTS_SOLID_FWD_MODLIQUID;
		end

		function set.expSolidFwdModLiquid(obj, val)
			obj.data.MAL_EXPONENTS_SOLID_FWD_MODLIQUID = val;
			obj.hasChanged = true;
		end

		function val = get.expSolidBwdModLiquid(obj)
			val = obj.data.MAL_EXPONENTS_SOLID_BWD_MODLIQUID;
		end

		function set.expSolidBwdModLiquid(obj, val)
			obj.data.MAL_EXPONENTS_SOLID_BWD_MODLIQUID = val;
			obj.hasChanged = true;
		end


		function val = get.MAL_KFWD_BULK(obj)
			val = obj.kFwdBulk;
		end

		function set.MAL_KFWD_BULK(obj, val)
			obj.kFwdBulk = val;
		end

		function val = get.MAL_KBWD_BULK(obj)
			val = obj.kBwdBulk;
		end

		function set.MAL_KBWD_BULK(obj, val)
			obj.kBwdBulk = val;
		end

		function val = get.MAL_KFWD_LIQUID(obj)
			val = obj.kFwdLiquid;
		end

		function set.MAL_KFWD_LIQUID(obj, val)
			obj.kFwdLiquid = val;
		end

		function val = get.MAL_KBWD_LIQUID(obj)
			val = obj.kBwdLiquid;
		end

		function set.MAL_KBWD_LIQUID(obj, val)
			obj.kBwdLiquid = val;
		end

		function val = get.MAL_KFWD_SOLID(obj)
			val = obj.kFwdSolid;
		end

		function set.MAL_KFWD_SOLID(obj, val)
			obj.kFwdSolid = val;
		end

		function val = get.MAL_KBWD_SOLID(obj)
			val = obj.kBwdSolid;
		end

		function set.MAL_KBWD_SOLID(obj, val)
			obj.kBwdSolid = val;
		end

		function val = get.MAL_STOICHIOMETRY_BULK(obj)
			val = obj.Sbulk;
		end

		function set.MAL_STOICHIOMETRY_BULK(obj, val)
			obj.Sbulk = val;
		end

		function val = get.MAL_STOICHIOMETRY_LIQUID(obj)
			val = obj.Sliquid;
		end

		function set.MAL_STOICHIOMETRY_LIQUID(obj, val)
			obj.Sliquid = val;
		end

		function val = get.MAL_STOICHIOMETRY_SOLID(obj)
			val = obj.Ssolid;
		end

		function set.MAL_STOICHIOMETRY_SOLID(obj, val)
			obj.Ssolid = val;
		end

		function val = get.MAL_EXPONENTS_BULK_FWD(obj)
			val = obj.expBulkFwd;
		end

		function set.MAL_EXPONENTS_BULK_FWD(obj, val)
			obj.expBulkFwd = val;
		end

		function val = get.MAL_EXPONENTS_BULK_BWD(obj)
			val = obj.expBulkBwd;
		end

		function set.MAL_EXPONENTS_BULK_BWD(obj, val)
			obj.expBulkBwd = val;
		end

		function val = get.MAL_EXPONENTS_LIQUID_FWD(obj)
			val = obj.expLiquidFwd;
		end

		function set.MAL_EXPONENTS_LIQUID_FWD(obj, val)
			obj.expLiquidFwd = val;
		end

		function val = get.MAL_EXPONENTS_LIQUID_BWD(obj)
			val = obj.expLiquidBwd;
		end

		function set.MAL_EXPONENTS_LIQUID_BWD(obj, val)
			obj.expLiquidBwd = val;
		end

		function val = get.MAL_EXPONENTS_LIQUID_FWD_MODSOLID(obj)
			val = obj.expLiquidFwdModSolid;
		end

		function set.MAL_EXPONENTS_LIQUID_FWD_MODSOLID(obj, val)
			obj.expLiquidFwdModSolid = val;
		end

		function val = get.MAL_EXPONENTS_LIQUID_BWD_MODSOLID(obj)
			val = obj.expLiquidBwdModSolid;
		end

		function set.MAL_EXPONENTS_LIQUID_BWD_MODSOLID(obj, val)
			obj.expLiquidBwdModSolid = val;
		end

		function val = get.MAL_EXPONENTS_SOLID_FWD(obj)
			val = obj.expSolidFwd;
		end

		function set.MAL_EXPONENTS_SOLID_FWD(obj, val)
			obj.expSolidFwd = val;
		end

		function val = get.MAL_EXPONENTS_SOLID_BWD(obj)
			val = obj.expSolidBwd;
		end

		function set.MAL_EXPONENTS_SOLID_BWD(obj, val)
			obj.expSolidBwd = val;
		end

		function val = get.MAL_EXPONENTS_SOLID_FWD_MODLIQUID(obj)
			val = obj.expSolidFwdModLiquid;
		end

		function set.MAL_EXPONENTS_SOLID_FWD_MODLIQUID(obj, val)
			obj.expSolidFwdModLiquid = val;
		end

		function val = get.MAL_EXPONENTS_SOLID_BWD_MODLIQUID(obj)
			val = obj.expSolidBwdModLiquid;
		end

		function set.MAL_EXPONENTS_SOLID_BWD_MODLIQUID(obj, val)
			obj.expSolidBwdModLiquid = val;
		end


		function val = get.hasBulkPhaseReactions(obj)
			val = (isfield(obj.data, 'MAL_STOICHIOMETRY_BULK') && ~isempty(obj.data.MAL_STOICHIOMETRY_BULK)) && ...
				(isfield(obj.data, 'MAL_KFWD_BULK') && ~isempty(obj.data.MAL_KFWD_BULK)) && ...
				(isfield(obj.data, 'MAL_KBWD_BULK') && ~isempty(obj.data.MAL_KBWD_BULK));
		end

		function val = get.hasLiquidPhaseReactions(obj)
			val = (isfield(obj.data, 'MAL_STOICHIOMETRY_LIQUID') && ~isempty(obj.data.MAL_STOICHIOMETRY_LIQUID)) && ...
				(isfield(obj.data, 'MAL_KFWD_LIQUID') && ~isempty(obj.data.MAL_KFWD_LIQUID)) && ...
				(isfield(obj.data, 'MAL_KBWD_LIQUID') && ~isempty(obj.data.MAL_KBWD_LIQUID));
		end

		function val = get.hasSolidPhaseReactions(obj)
			val = (isfield(obj.data, 'MAL_STOICHIOMETRY_SOLID') && ~isempty(obj.data.MAL_STOICHIOMETRY_SOLID)) && ...
				(isfield(obj.data, 'MAL_KFWD_SOLID') && ~isempty(obj.data.MAL_KFWD_SOLID)) && ...
				(isfield(obj.data, 'MAL_KBWD_SOLID') && ~isempty(obj.data.MAL_KBWD_SOLID));
		end


	end

	methods (Static)

		function obj = loadobj(S)
			if isstruct(S)
				obj = MassActionLawReaction();
				obj.loadobjInternal(S);
			end
		end
		
	end

end

% =============================================================================
%  CADET - The Chromatography Analysis and Design Toolkit
%  
%  Copyright (C) 2008-2019: The CADET Authors
%            Please see the AUTHORS and CONTRIBUTORS file.
%  
%  All rights reserved. obj program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies obj distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
