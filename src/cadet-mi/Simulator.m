
classdef Simulator
	%Simulator Helper class to create CADET simulators using different interfaces
	%   This class provides methods to instantiate CADET simulators that use
	%   a specific interface (e.g., Matlab MEX, HDF5).
	%
	% See also MEXSIMULATOR
	
	% Copyright: (C) 2008-2020 The CADET Authors
	%            See the license note at the end of the file.
	
	properties(Constant)
		preferredInterface = Simulator.getPreferredInterface(); % Preferred (working) CADET interface
		availableInterfaceClasses = {'MexSimulator'};
		availableInterfaces = {'mex'};
	end

	methods (Static, Access = 'public')

		function interface = getPreferredInterface()
			%GETPREFERREDINTERFACE Returns the preferred working CADET interface
			%   INTERFACE = GETPREFERREDINTERFACE() returns the most preferred working
			%   CADET interface as string in INTERFACE. Returns an empty array if
			%   no working interface could be found.

			interface = [];

			% Try MEX interface
			try
				% The next method fails if the MEX interface is not working
				MexSimulator.getVersion();
				interface = 'mex';
				return;
			catch
				interface = [];
			end

		end

		function sim = create(interface)
			%CREATE Creates a simulator (using some specific interface)
			%   SIM = CREATE() tries different interfaces in the order
			%   of preference (MEX > HDF5) and returns a simulator object
			%   in SIM for the first working interface. Returns an empty
			%   array if no interface works.
			%
			%   SIM = CREATE(INTERFACE) tries to create a simulator using
			%   the interface given in the string INTERFACE ('mex', 'hdf5').
			%   Returns an empty array in SIM if the interface does not work.

			sim = [];
			if nargin == 0
				% Create preferred interface if there is any
				if ~isempty(Simulator.preferredInterface)
					sim = Simulator.create(Simulator.preferredInterface);
				end
			else
				% Create specific interface
				switch (interface)
					case 'mex'
						try
							% The next method fails if the MEX interface is not working
							MexSimulator.getVersion();
							sim = MexSimulator();
						catch
							sim = [];
						end
					otherwise
						error('CADET', 'Interface %s is not available.', interface);
				end				
			end
		end

		function [version, commit, branch] = getVersion()
			%GETVERSION Returns the version of CADET, commit hash and branch it was built from
			%   [VERSION] = getVersion() returns the version string VERSION of the underlying CADET
			%   simulator. Returns an empty array if no working interface has been found.
			%
			%   [VERSION, COMMIT] = getVersion() returns version string VERSION and git commit hash
			%   COMMIT. Returns empty arrays if no working interface has been found.
			%
			%   [VERSION, COMMIT, BRANCH] = getVersion() returns version string VERSION, git commit
			%   hash COMMIT, and branch name BRANCH. Returns empty arrays if no working interface
			%   has been found.
			
			switch Simulator.preferredInterface
				case 'mex'
					[version, commit, branch] = MexSimulator.getVersion();
				otherwise
					version = [];
					commit = [];
					branch = [];
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
%  All rights reserved. obj program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies obj distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
