// =============================================================================
//  CADET
//  
//  Copyright © 2008-2021: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

/**
 * @file 
 * Provides several implementations of ISolutionRecorder.
 */

#ifndef LIBCADET_SOLUTIONRECORDER_IMPL_HPP_
#define LIBCADET_SOLUTIONRECORDER_IMPL_HPP_

#include <vector>
#include <sstream>
#include <algorithm>
#include <numeric>

#include "cadet/SolutionRecorder.hpp"

namespace cadet
{

namespace detail
{
	CADET_CONST_OR_CONSTEXPR unsigned int numDefaultRecorderTimesteps = 100u;
}

/**
 * @brief Stores pieces of the solution of one single unit operation in internal buffers
 * @details The pieces of stored solutions are selectable at runtime.
 * @todo Use better storage than std::vector (control growth, maybe chunked storage -> needs chunked writes)
 */
class InternalStorageUnitOpRecorder : public ISolutionRecorder
{
public:

	struct StorageConfig
	{
		bool storeBulk;
		bool storeParticle;
		bool storeSolid;
		bool storeFlux;
		bool storeOutlet;
		bool storeInlet;
		bool storeVolume;
		bool storeLast;
	};

	InternalStorageUnitOpRecorder() : InternalStorageUnitOpRecorder(UnitOpIndep) { }

	InternalStorageUnitOpRecorder(UnitOpIdx idx) : _cfgSolution({false, false, false, true, false, false, false}),
		_cfgSolutionDot({false, false, false, false, false, false, false}), _cfgSensitivity({false, false, false, true, false, false, false}),
		_cfgSensitivityDot({false, false, false, true, false, false, false}), _storeTime(false), _storeCoordinates(false), _splitComponents(true), _splitPorts(true),
		_singleAsMultiPortUnitOps(false), _curCfg(nullptr), _nComp(0), _nVolumeDof(0), _numTimesteps(0), _numSens(0), _unitOp(idx), _needsReAlloc(false),
		_axialCoords(0), _radialCoords(0), _particleCoords(0)
	{
	}

	virtual ~InternalStorageUnitOpRecorder() CADET_NOEXCEPT
	{
	}

	virtual void clear()
	{
		// Clear solution storage
		_time.clear();
		clear(_data);
		clear(_dataDot);

		// Clear all sensitivity storage
		for (std::size_t i = 0; i < _sens.size(); ++i)
		{
			clear(_sens[i]);
			clear(_sensDot[i]);
		}
	}

	virtual void prepare(unsigned int numDofs, unsigned int numSens, unsigned int numTimesteps)
	{
		_numTimesteps = numTimesteps;
		_numSens = numSens;

		// Allocate sensitivity storage
		_sens.resize(numSens);
		_sensDot.resize(numSens);

		_needsReAlloc = false;
	}

	virtual void notifyIntegrationStart(unsigned int numDofs, unsigned int numSens, unsigned int numTimesteps)
	{
		_needsReAlloc = (numSens != _numSens) || (numTimesteps > _numTimesteps) || (numTimesteps == 0);

		// Clear all data from memory
		clear();

		_numTimesteps = numTimesteps;
		
		if (numSens != _numSens)
		{
			// Allocate sensitivity storage
			_sens.resize(numSens);
			_sensDot.resize(numSens);

			_numSens = numSens;
		}
	}

	virtual void unitOperationStructure(UnitOpIdx idx, const IModel& model, const ISolutionExporter& exporter)
	{
		// Only record one unit operation
		if (idx != _unitOp)
			return;

		_nComp = exporter.numComponents();
		_nVolumeDof = exporter.numVolumeDofs();
		_nInletPorts = exporter.numInletPorts();
		_nOutletPorts = exporter.numOutletPorts();

		// Query particle type specific structure
		const unsigned int numParTypes = exporter.numParticleTypes();
		_nParShells.resize(numParTypes, 0u);
		_nBoundStates.resize(numParTypes, 0u);
		for (unsigned int i = 0; i < numParTypes; ++i)
		{
			_nParShells[i] = exporter.numParticleShells(i);
			_nBoundStates[i] = exporter.numBoundStates(i);
		}

		// Query structure
		unsigned int len = 0;
		StateOrdering const* order = exporter.concentrationOrdering(len);
		_bulkLayout.clear();
		_bulkLayout.reserve(len + 1); // First slot is time
		_bulkLayout.push_back(0);
		_bulkCount = 1;
		for (unsigned int i = 0; i < len; ++i)
		{
			switch (order[i])
			{
				case StateOrdering::Component:
					_bulkLayout.push_back(exporter.numComponents());
					break;
				case StateOrdering::AxialCell:
					_bulkLayout.push_back(exporter.numAxialCells());
					_bulkCount *= exporter.numAxialCells();
					break;
				case StateOrdering::RadialCell:
					_bulkLayout.push_back(exporter.numRadialCells());
					_bulkCount *= exporter.numRadialCells();
				case StateOrdering::ParticleType:
				case StateOrdering::ParticleShell:
				case StateOrdering::BoundState:
					break;
			}
		}

		order = exporter.mobilePhaseOrdering(len);
		_particleLayout.clear();
		_particleLayout.resize(numParTypes, std::vector<std::size_t>(len + 1, 0)); // First slot is time
		_particleCount = std::vector<unsigned int>(numParTypes, 1u);
		unsigned int idxLayout = 1;
		for (unsigned int i = 0; i < len; ++i)
		{
			switch (order[i])
			{
				case StateOrdering::Component:
				{
					for (unsigned int j = 0; j < numParTypes; ++j)
						_particleLayout[j][idxLayout] = exporter.numComponents();

					++idxLayout;
					break;
				}
				case StateOrdering::AxialCell:
				{
					for (unsigned int j = 0; j < numParTypes; ++j)
					{
						_particleLayout[j][idxLayout] = exporter.numAxialCells();
						_particleCount[j] *= exporter.numAxialCells();
					}

					++idxLayout;
					break;
				}
				case StateOrdering::RadialCell:
				{
					for (unsigned int j = 0; j < numParTypes; ++j)
					{
						_particleLayout[j][idxLayout] = exporter.numRadialCells();
						_particleCount[j] *= exporter.numRadialCells();
					}

					++idxLayout;
					break;
				}
				case StateOrdering::ParticleType:
					break;
				case StateOrdering::ParticleShell:
				{
					for (unsigned int j = 0; j < numParTypes; ++j)
					{
						_particleLayout[j][idxLayout] = _nParShells[j];
						_particleCount[j] *= _nParShells[j];
					}

					++idxLayout;
					break;
				}
				case StateOrdering::BoundState:
					break;
			}
		}

		for (unsigned int j = 0; j < numParTypes; ++j)
			_particleLayout[j].resize(idxLayout);

		order = exporter.solidPhaseOrdering(len);
		_solidLayout.clear();
		_solidLayout.resize(numParTypes, std::vector<std::size_t>(len + 1, 0)); // First slot is time
		_solidCount = std::vector<unsigned int>(numParTypes, 1u);
		idxLayout = 1;
		for (unsigned int i = 0; i < len; ++i)
		{
			switch (order[i])
			{
				case StateOrdering::Component:
					break;
				case StateOrdering::AxialCell:
				{
					for (unsigned int j = 0; j < numParTypes; ++j)
					{
						_solidLayout[j][idxLayout] = exporter.numAxialCells();
						_solidCount[j] *= exporter.numAxialCells();
					}

					++idxLayout;
					break;
				}
				case StateOrdering::RadialCell:
				{
					for (unsigned int j = 0; j < numParTypes; ++j)
					{
						_solidLayout[j][idxLayout] = exporter.numRadialCells();
						_solidCount[j] *= exporter.numRadialCells();
					}

					++idxLayout;
					break;
				}
				case StateOrdering::ParticleType:
					break;
				case StateOrdering::ParticleShell:
				{
					for (unsigned int j = 0; j < numParTypes; ++j)
					{
						_solidLayout[j][idxLayout] = _nParShells[j];
						_solidCount[j] *= _nParShells[j];
					}

					++idxLayout;
					break;
				}
				case StateOrdering::BoundState:
				{
					for (unsigned int j = 0; j < numParTypes; ++j)
						_solidLayout[j][idxLayout] = _nBoundStates[j];

					++idxLayout;
					break;
				}
			}
		}

		for (unsigned int j = 0; j < numParTypes; ++j)
			_solidLayout[j].resize(idxLayout);

		order = exporter.fluxOrdering(len);
		_fluxLayout.clear();
		_fluxLayout.reserve(len + 1); // First slot is time
		_fluxLayout.push_back(0);
		for (unsigned int i = 0; i < len; ++i)
		{
			switch (order[i])
			{
				case StateOrdering::Component:
					_fluxLayout.push_back(exporter.numComponents());
					break;
				case StateOrdering::AxialCell:
					_fluxLayout.push_back(exporter.numAxialCells());
					break;
				case StateOrdering::ParticleType:
					_fluxLayout.push_back(numParTypes);
					break;
				case StateOrdering::RadialCell:
					_fluxLayout.push_back(exporter.numRadialCells());
					break;
				case StateOrdering::ParticleShell:
				case StateOrdering::BoundState:
					break;
			}
		}

		// Obtain coordinates
		if (_storeCoordinates)
		{
			_axialCoords.resize(exporter.numAxialCells());
			_radialCoords.resize(exporter.numRadialCells());
			const unsigned int numShells = std::accumulate(_nParShells.begin(), _nParShells.end(), 0u);
			_particleCoords.resize(numShells);

			exporter.axialCoordinates(_axialCoords.data());
			exporter.radialCoordinates(_radialCoords.data());

			unsigned int offset = 0;
			for (unsigned int i = 0; i < numParTypes; ++i)
			{
				exporter.particleCoordinates(i, _particleCoords.data() + offset);
				offset += _nParShells[i];
			}
		}

		// Validate config
		validateConfig(exporter, _cfgSolution);
		validateConfig(exporter, _cfgSolutionDot);
		validateConfig(exporter, _cfgSensitivity);
		validateConfig(exporter, _cfgSensitivityDot);

		// Everything is ok, we have nothing to do
		if (!_needsReAlloc)
		{
			// Reset for counting the number of received time steps
			_numTimesteps = 0;
			return;
		}

		// Allocate space for solution
		beginSolution();
		allocateMemory(exporter);		
		endSolution();

		beginSolutionDerivative();
		allocateMemory(exporter);		
		endSolution();

		// Allocate space for sensitivities
		for (std::size_t i = 0; i < _sens.size(); ++i)
		{
			beginSensitivity(i);
			allocateMemory(exporter);
			endSolution();

			beginSensitivityDot(i);
			allocateMemory(exporter);
			endSolution();
		}

		if (_storeTime)
		{
			if (_numTimesteps == 0)
				_time.reserve(detail::numDefaultRecorderTimesteps);
			else
				_time.reserve(_numTimesteps);
		}

		// Reset for counting the number of received time steps
		_numTimesteps = 0;
	}

	virtual void beginTimestep(double t)
	{
		++_numTimesteps;
		if (_storeTime)
			_time.push_back(t);
	}

	virtual void beginUnitOperation(cadet::UnitOpIdx idx, const cadet::IModel& model, const cadet::ISolutionExporter& exporter)
	{
		// Only record one unit operation
		if ((idx != _unitOp) || !_curCfg)
			return;

		unsigned int stride = 0;

		if (_curCfg->storeOutlet)
		{
			std::vector<double>& v = _curStorage->outlet;
			for (unsigned int j = 0; j < _nOutletPorts; ++j)
			{
				double const* outlet = exporter.outlet(j, stride);
				for (unsigned int i = 0; i < _nComp; ++i)
					v.push_back(outlet[i * stride]);
			}
		}

		if (_curCfg->storeInlet)
		{
			std::vector<double>& v = _curStorage->inlet;
			for (unsigned int j = 0; j < _nInletPorts; ++j)
			{
				double const* inlet = exporter.inlet(j, stride);
				for (unsigned int i = 0; i < _nComp; ++i)
					v.push_back(inlet[i * stride]);
			}
		}

		if (_curCfg->storeBulk)
		{
			std::vector<double>& v = _curStorage->bulk;
			double const* data = exporter.concentration();
			stride = exporter.bulkMobilePhaseStride();
			const unsigned int blockSize = exporter.numBulkDofs() / _bulkCount;
			for (unsigned int i = 0; i < _bulkCount; ++i, data += stride)
			{
				v.insert(v.end(), data, data + blockSize);
			}
		}

		if (_curCfg->storeParticle)
		{
			for (std::size_t parType = 0; parType < _nParShells.size(); ++parType)
			{
				std::vector<double>& cp = _curStorage->particle[parType];
				double const* data = exporter.particleMobilePhase(parType);
				stride = exporter.particleMobilePhaseStride(parType);
				const unsigned int blockSize = exporter.numParticleMobilePhaseDofs(parType) / _particleCount[parType];
				for (unsigned int i = 0; i < _particleCount[parType]; ++i, data += stride)
				{
					cp.insert(cp.end(), data, data + blockSize);
				}
			}
		}

		if (_curCfg->storeSolid)
		{
			for (std::size_t parType = 0; parType < _nParShells.size(); ++parType)
			{
				std::vector<double>& cs = _curStorage->solid[parType];
				double const* data = exporter.solidPhase(parType);
				stride = exporter.solidPhaseStride(parType);
				const unsigned int blockSize = exporter.numSolidPhaseDofs(parType) / _solidCount[parType];
				for (unsigned int i = 0; i < _solidCount[parType]; ++i, data += stride)
				{
					cs.insert(cs.end(), data, data + blockSize);
				}
			}
		}

		if (_curCfg->storeFlux)
		{
			double const* const data = exporter.flux();
			_curStorage->flux.insert(_curStorage->flux.end(), data, data + exporter.numFluxDofs());
		}

		if (_curCfg->storeVolume)
		{
			double const* const data = exporter.volume();
			_curStorage->volume.insert(_curStorage->volume.end(), data, data + exporter.numVolumeDofs());
		}
	}

	virtual void endUnitOperation() { }

	virtual void endTimestep() { }

	virtual void beginSolution()
	{
		_curCfg = &_cfgSolution;
		_curStorage = &_data;
	}

	virtual void endSolution()
	{
		_curCfg = nullptr;
		_curStorage = nullptr;
	}

	virtual void beginSolutionDerivative()
	{
		_curCfg = &_cfgSolutionDot;
		_curStorage = &_dataDot;
	}

	virtual void endSolutionDerivative() { endSolution(); }

	virtual void beginSensitivity(const cadet::ParameterId& pId, unsigned int sensIdx)
	{
		beginSensitivity(sensIdx);
	}

	virtual void endSensitivity(const cadet::ParameterId& pId, unsigned int sensIdx)
	{
		endSolution();
	}
	
	virtual void beginSensitivityDerivative(const cadet::ParameterId& pId, unsigned int sensIdx)
	{
		beginSensitivityDot(sensIdx);
	}

	virtual void endSensitivityDerivative(const cadet::ParameterId& pId, unsigned int sensIdx)
	{
		endSolution();
	}

	template <typename Writer_t>
	void writeCoordinates(Writer_t& writer)
	{
		if (!_storeCoordinates)
			return;

		if (!_axialCoords.empty())
			writer.template vector<double>("AXIAL_COORDINATES", _axialCoords);
		if (!_radialCoords.empty())
			writer.template vector<double>("RADIAL_COORDINATES", _radialCoords);

		if (!_particleCoords.empty())
		{
			std::ostringstream oss;
			unsigned int offset = 0;
			for (std::size_t pt = 0; pt < _nParShells.size(); ++pt)
			{
				if (_nParShells[pt] == 0)
					continue;

				oss.str("");
				oss << "PARTICLE_COORDINATES_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << pt;

				writer.template vector<double>(oss.str(), _nParShells[pt], _particleCoords.data() + offset);

				offset += _nParShells[pt];
			}
		}
	}

	template <typename Writer_t>
	void writeSolution(Writer_t& writer)
	{
		std::ostringstream oss;

		if (_storeTime)
			writer.template vector<double>("SOLUTION_TIMES", _time.size(), _time.data());

		beginSolution();
		writeData(writer, "SOLUTION", oss);
		endSolution();

		beginSolutionDerivative();
		writeData(writer, "SOLDOT", oss);
		endSolution();
	}

	template <typename Writer_t>
	void writeSensitivity(Writer_t& writer)
	{
		std::ostringstream oss;

		for (unsigned int param = 0; param < _numSens; ++param)
		{
			oss.str("");
			oss << "param_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << param;
			writer.pushGroup(oss.str());

			beginSensitivity(param);
			writeData(writer, "SENS", oss);
			endSolution();

			beginSensitivityDot(param);
			writeData(writer, "SENSDOT", oss);
			endSolution();

			writer.popGroup();
		}
	}

	template <typename Writer_t>
	void writeSensitivity(Writer_t& writer, unsigned int param)
	{
		std::ostringstream oss;

		beginSensitivity(param);
		writeData(writer, "SENS", oss);
		endSolution();

		beginSensitivityDot(param);
		writeData(writer, "SENSDOT", oss);
		endSolution();
	}

	inline StorageConfig& solutionConfig() CADET_NOEXCEPT { return _cfgSolution; }
	inline const StorageConfig& solutionConfig() const CADET_NOEXCEPT { return _cfgSolution; }
	inline void solutionConfig(const StorageConfig& cfg) CADET_NOEXCEPT { _cfgSolution = cfg; }

	inline StorageConfig& solutionDotConfig() CADET_NOEXCEPT { return _cfgSolutionDot; }
	inline const StorageConfig& solutionDotConfig() const CADET_NOEXCEPT { return _cfgSolutionDot; }
	inline void solutionDotConfig(const StorageConfig& cfg) CADET_NOEXCEPT { _cfgSolutionDot = cfg; }

	inline StorageConfig& sensitivityConfig() CADET_NOEXCEPT { return _cfgSensitivity; }
	inline const StorageConfig& sensitivityConfig() const CADET_NOEXCEPT { return _cfgSensitivity; }
	inline void sensitivityConfig(const StorageConfig& cfg) CADET_NOEXCEPT { _cfgSensitivity = cfg; }

	inline StorageConfig& sensitivityDotConfig() CADET_NOEXCEPT { return _cfgSensitivityDot; }
	inline const StorageConfig& sensitivityDotConfig() const CADET_NOEXCEPT { return _cfgSensitivityDot; }
	inline void sensitivityDotConfig(const StorageConfig& cfg) CADET_NOEXCEPT { _cfgSensitivityDot = cfg; }

	inline bool storeTime() const CADET_NOEXCEPT { return _storeTime; }
	inline void storeTime(bool st) CADET_NOEXCEPT { _storeTime = st; }

	inline bool storeCoordinates() const CADET_NOEXCEPT { return _storeCoordinates; }
	inline void storeCoordinates(bool sc) CADET_NOEXCEPT { _storeCoordinates = sc; }

	inline bool splitComponents() const CADET_NOEXCEPT { return _splitComponents; }
	inline void splitComponents(bool st) CADET_NOEXCEPT { _splitComponents = st; }

	inline bool splitPorts() const CADET_NOEXCEPT { return _splitPorts; }
	inline void splitPorts(bool st) CADET_NOEXCEPT { _splitPorts = st; }

	inline bool treatSingleAsMultiPortUnitOps() const CADET_NOEXCEPT { return _singleAsMultiPortUnitOps; }
	inline void treatSingleAsMultiPortUnitOps(bool smp) CADET_NOEXCEPT { _singleAsMultiPortUnitOps = smp; }

	inline UnitOpIdx unitOperation() const CADET_NOEXCEPT { return _unitOp; }
	inline void unitOperation(UnitOpIdx idx) CADET_NOEXCEPT { _unitOp = idx; }

	inline unsigned int numDataPoints() const CADET_NOEXCEPT { return _numTimesteps; }
	inline unsigned int numComponents() const CADET_NOEXCEPT { return _nComp; }
	inline unsigned int numInletPorts() const CADET_NOEXCEPT { return _nInletPorts; }
	inline unsigned int numOutletPorts() const CADET_NOEXCEPT { return _nOutletPorts; }

	inline double const* time() const CADET_NOEXCEPT { return _time.data(); }
	inline double const* inlet() const CADET_NOEXCEPT { return _data.inlet.data(); }
	inline double const* outlet() const CADET_NOEXCEPT { return _data.outlet.data(); }
	inline double const* bulk() const CADET_NOEXCEPT { return _data.bulk.data(); }
	inline double const* particle(unsigned int parType = 0) const CADET_NOEXCEPT { return _data.particle[parType].data(); }
	inline double const* solid(unsigned int parType = 0) const CADET_NOEXCEPT { return _data.solid[parType].data(); }
	inline double const* flux() const CADET_NOEXCEPT { return _data.flux.data(); }
	inline double const* volume() const CADET_NOEXCEPT { return _data.volume.data(); }
	inline double const* inletDot() const CADET_NOEXCEPT { return _dataDot.inlet.data(); }
	inline double const* outletDot() const CADET_NOEXCEPT { return _dataDot.outlet.data(); }
	inline double const* bulkDot() const CADET_NOEXCEPT { return _dataDot.bulk.data(); }
	inline double const* particleDot(unsigned int parType = 0) const CADET_NOEXCEPT { return _dataDot.particle[parType].data(); }
	inline double const* solidDot(unsigned int parType = 0) const CADET_NOEXCEPT { return _dataDot.solid[parType].data(); }
	inline double const* fluxDot() const CADET_NOEXCEPT { return _dataDot.flux.data(); }
	inline double const* volumeDot() const CADET_NOEXCEPT { return _dataDot.volume.data(); }
	inline double const* sensInlet(unsigned int idx) const CADET_NOEXCEPT { return _sens[idx].inlet.data(); }
	inline double const* sensOutlet(unsigned int idx) const CADET_NOEXCEPT { return _sens[idx].outlet.data(); }
	inline double const* sensBulk(unsigned int idx) const CADET_NOEXCEPT { return _sens[idx].bulk.data(); }
	inline double const* sensParticle(unsigned int idx, unsigned int parType = 0) const CADET_NOEXCEPT { return _sens[idx].particle[parType].data(); }
	inline double const* sensSolid(unsigned int idx, unsigned int parType = 0) const CADET_NOEXCEPT { return _sens[idx].solid[parType].data(); }
	inline double const* sensFlux(unsigned int idx) const CADET_NOEXCEPT { return _sens[idx].flux.data(); }
	inline double const* sensVolume(unsigned int idx) const CADET_NOEXCEPT { return _sens[idx].volume.data(); }
	inline double const* sensInletDot(unsigned int idx) const CADET_NOEXCEPT { return _sensDot[idx].inlet.data(); }
	inline double const* sensOutletDot(unsigned int idx) const CADET_NOEXCEPT { return _sensDot[idx].outlet.data(); }
	inline double const* sensBulkDot(unsigned int idx) const CADET_NOEXCEPT { return _sensDot[idx].bulk.data(); }
	inline double const* sensParticleDot(unsigned int idx, unsigned int parType = 0) const CADET_NOEXCEPT { return _sensDot[idx].particle[parType].data(); }
	inline double const* sensSolidDot(unsigned int idx, unsigned int parType = 0) const CADET_NOEXCEPT { return _sensDot[idx].solid[parType].data(); }
	inline double const* sensFluxDot(unsigned int idx) const CADET_NOEXCEPT { return _sensDot[idx].flux.data(); }
	inline double const* sensVolumeDot(unsigned int idx) const CADET_NOEXCEPT { return _sensDot[idx].volume.data(); }
protected:

	struct Storage
	{
		std::vector<double> outlet;
		std::vector<double> inlet;
		std::vector<double> bulk;
		std::vector<std::vector<double>> particle;
		std::vector<std::vector<double>> solid;
		std::vector<double> flux;
		std::vector<double> volume;
	};

	inline void beginSensitivity(unsigned int sensIdx)
	{
		_curCfg = &_cfgSensitivity;
		_curStorage = &_sens[sensIdx];
	}

	inline void beginSensitivityDot(unsigned int sensIdx)
	{
		_curCfg = &_cfgSensitivityDot;
		_curStorage = &_sensDot[sensIdx];
	}

	inline void validateConfig(const ISolutionExporter& exporter, StorageConfig& cfg)
	{
		// Only store fields that really exist
		cfg.storeOutlet = (exporter.numOutletPorts() > 0) && cfg.storeOutlet;
		cfg.storeInlet = (exporter.numInletPorts() > 0) && cfg.storeInlet;
		cfg.storeParticle = exporter.hasParticleMobilePhase() && cfg.storeParticle;
		cfg.storeSolid = exporter.hasSolidPhase() && cfg.storeSolid;
		cfg.storeFlux = exporter.hasParticleFlux() && cfg.storeFlux;
		cfg.storeVolume = exporter.hasVolume() && cfg.storeVolume;
	}

	inline void allocateMemory(const ISolutionExporter& exporter)
	{
		const unsigned int nAllocTimesteps = std::max(_numTimesteps, detail::numDefaultRecorderTimesteps);

		if (_curCfg->storeOutlet)
			_curStorage->outlet.reserve(nAllocTimesteps * _nComp * _nOutletPorts);

		if (_curCfg->storeInlet)
			_curStorage->inlet.reserve(nAllocTimesteps * _nComp * _nInletPorts);

		if (_curCfg->storeBulk)
			_curStorage->bulk.reserve(nAllocTimesteps * exporter.numBulkDofs());
		
		if (_curCfg->storeParticle)
		{
			_curStorage->particle.resize(_nParShells.size());
			for (std::size_t i = 0; i < _nParShells.size(); ++i)
				_curStorage->particle[i].reserve(nAllocTimesteps * exporter.numParticleMobilePhaseDofs(i));
		}
		
		if (_curCfg->storeSolid)
		{
			_curStorage->solid.resize(_nParShells.size());
			for (std::size_t i = 0; i < _nParShells.size(); ++i)
				_curStorage->solid[i].reserve(nAllocTimesteps * exporter.numSolidPhaseDofs(i));
		}

		if (_curCfg->storeFlux)
			_curStorage->flux.reserve(nAllocTimesteps * exporter.numFluxDofs());
		
		if (_curCfg->storeVolume)
			_curStorage->volume.reserve(nAllocTimesteps * exporter.numVolumeDofs());
	}

	template <typename Writer_t>
	void writeData(Writer_t& writer, const char* prefix, std::ostringstream& oss)
	{
		if (_curCfg->storeOutlet)
		{
			if (_splitPorts)
			{
				if (_splitComponents)
				{
					for (unsigned int port = 0; port < _nOutletPorts; ++port)
					{
						for (unsigned int comp = 0; comp < _nComp; ++comp)
						{
							oss.str("");
							if ((_nOutletPorts == 1) && !_singleAsMultiPortUnitOps)
							{
								oss << prefix << "_OUTLET_COMP_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << comp;
							}
							else
							{
								oss << prefix << "_OUTLET_PORT_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << port 
									<<  "_COMP_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << comp;
							}

							writer.template vector<double>(oss.str(), _numTimesteps, _curStorage->outlet.data() + comp + port * _nComp, _nComp * _nOutletPorts);
						}
					}
				}
				else
				{
					for (unsigned int port = 0; port < _nOutletPorts; ++port)
					{
						oss.str("");
						if ((_nOutletPorts == 1) && !_singleAsMultiPortUnitOps)
							oss << prefix << "_OUTLET";
						else
							oss << prefix << "_OUTLET_PORT_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << port;

						writer.template matrix<double>(oss.str(), _numTimesteps, _nComp, _curStorage->outlet.data() + port * _nComp, _nOutletPorts * _nComp, _nComp);
					}
				}
			}
			else
			{
				if (_splitComponents)
				{
					for (unsigned int comp = 0; comp < _nComp; ++comp)
					{
						oss.str("");
						oss << prefix << "_OUTLET_COMP_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << comp;
						if ((_nOutletPorts == 1) && !_singleAsMultiPortUnitOps)
							writer.template vector<double>(oss.str(), _numTimesteps, _curStorage->outlet.data() + comp, _nComp);
						else
							writer.template matrix<double>(oss.str(), _numTimesteps, _nOutletPorts, _curStorage->outlet.data() + comp, _nComp);
					}
				}
				else
				{
					oss.str("");
					oss << prefix << "_OUTLET";
					if ((_nOutletPorts == 1) && !_singleAsMultiPortUnitOps)
					{
						const std::vector<std::size_t> layout = {_numTimesteps, _nComp};
						writer.template tensor<double>(oss.str(), layout.size(), layout.data(), _curStorage->outlet.data());
					}
					else
					{
						const std::vector<std::size_t> layout = {_numTimesteps, _nOutletPorts, _nComp};
						writer.template tensor<double>(oss.str(), layout.size(), layout.data(), _curStorage->outlet.data());
					}
				}
			}
		}

		if (_curCfg->storeInlet)
		{
			if (_splitPorts)
			{
				if (_splitComponents)
				{
					for (unsigned int port = 0; port < _nInletPorts; ++port)
					{
						for (unsigned int comp = 0; comp < _nComp; ++comp)
						{
							oss.str("");
							if ((_nInletPorts == 1) && !_singleAsMultiPortUnitOps)
							{
								oss << prefix << "_INLET_COMP_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << comp;
							}
							else
							{
								oss << prefix << "_INLET_PORT_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << port 
									<<  "_COMP_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << comp;
							}

							writer.template vector<double>(oss.str(), _numTimesteps, _curStorage->inlet.data() + comp + port * _nComp, _nComp * _nInletPorts);
						}
					}
				}
				else
				{
					for (unsigned int port = 0; port < _nInletPorts; ++port)
					{
						oss.str("");
						if ((_nInletPorts == 1) && !_singleAsMultiPortUnitOps)
							oss << prefix << "_INLET";
						else
							oss << prefix << "_INLET_PORT_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << port;

						writer.template matrix<double>(oss.str(), _numTimesteps, _nComp, _curStorage->inlet.data() + port * _nComp, _nInletPorts * _nComp, _nComp);
					}
				}
			}
			else
			{
				if (_splitComponents)
				{
					for (unsigned int comp = 0; comp < _nComp; ++comp)
					{
						oss.str("");
						oss << prefix << "_INLET_COMP_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << comp;
						if ((_nInletPorts == 1) && !_singleAsMultiPortUnitOps)
							writer.template vector<double>(oss.str(), _numTimesteps, _curStorage->inlet.data() + comp, _nComp);
						else
							writer.template matrix<double>(oss.str(), _numTimesteps, _nInletPorts, _curStorage->inlet.data() + comp, _nComp);
					}
				}
				else
				{
					oss.str("");
					oss << prefix << "_INLET";
					if ((_nInletPorts == 1) && !_singleAsMultiPortUnitOps)
					{
						const std::vector<std::size_t> layout = {_numTimesteps, _nComp};
						writer.template tensor<double>(oss.str(), layout.size(), layout.data(), _curStorage->inlet.data());
					}
					else
					{
						const std::vector<std::size_t> layout = {_numTimesteps, _nInletPorts, _nComp};
						writer.template tensor<double>(oss.str(), layout.size(), layout.data(), _curStorage->inlet.data());
					}
				}
			}
		}

		if (_curCfg->storeBulk)
		{
			oss.str("");
			oss << prefix << "_BULK";
			_bulkLayout[0] = _numTimesteps;
			writer.template tensor<double>(oss.str(), _bulkLayout.size(), _bulkLayout.data(), _curStorage->bulk.data());
		}

		if (_curCfg->storeParticle)
		{
			if (_nParShells.size() <= 1)
			{
				oss.str("");
				oss << prefix << "_PARTICLE";
				_particleLayout[0][0] = _numTimesteps;
				writer.template tensor<double>(oss.str(), _particleLayout[0].size(), _particleLayout[0].data(), _curStorage->particle[0].data());
			}
			else
			{
				for (std::size_t parType = 0; parType < _nParShells.size(); ++parType)
				{
					std::vector<std::size_t>& pl = _particleLayout[parType];
					oss.str("");
					oss << prefix << "_PARTICLE_PARTYPE_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << parType;
					pl[0] = _numTimesteps;
					writer.template tensor<double>(oss.str(), pl.size(), pl.data(), _curStorage->particle[parType].data());
				}
			}
		}

		if (_curCfg->storeSolid)
		{
			if (_nParShells.size() <= 1)
			{
				oss.str("");
				oss << prefix << "_SOLID";
				_solidLayout[0][0] = _numTimesteps;
				writer.template tensor<double>(oss.str(), _solidLayout[0].size(), _solidLayout[0].data(), _curStorage->solid[0].data());
			}
			else
			{
				for (std::size_t parType = 0; parType < _nParShells.size(); ++parType)
				{
					std::vector<std::size_t>& pl = _solidLayout[parType];
					oss.str("");
					oss << prefix << "_SOLID_PARTYPE_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << parType;
					pl[0] = _numTimesteps;
					writer.template tensor<double>(oss.str(), pl.size(), pl.data(), _curStorage->solid[parType].data());
				}
			}
		}

		if (_curCfg->storeFlux)
		{
			oss.str("");
			oss << prefix << "_FLUX";
			_fluxLayout[0] = _numTimesteps;
			writer.template tensor<double>(oss.str(), _fluxLayout.size(), _fluxLayout.data(), _curStorage->flux.data());
		}

		if (_curCfg->storeVolume)
		{
			oss.str("");
			oss << prefix << "_VOLUME";
			_fluxLayout[0] = _numTimesteps;
			writer.template matrix<double>(oss.str(), _numTimesteps, _nVolumeDof, _curStorage->volume.data(), 1);
		}
	}

	inline void clear(Storage& s)
	{
		s.outlet.clear();
		s.inlet.clear();
		s.bulk.clear();

		for (auto& v : s.particle)
			v.clear();

		for (auto& v : s.solid)
			v.clear();

		s.flux.clear();
		s.volume.clear();
	}

	StorageConfig _cfgSolution;
	StorageConfig _cfgSolutionDot;
	StorageConfig _cfgSensitivity;
	StorageConfig _cfgSensitivityDot;
	bool _storeTime;
	bool _storeCoordinates;
	bool _splitComponents;
	bool _splitPorts;
	bool _singleAsMultiPortUnitOps;

	StorageConfig const* _curCfg;
	Storage* _curStorage;

	std::vector<double> _time;
	Storage _data;
	Storage _dataDot;
	std::vector<Storage> _sens;
	std::vector<Storage> _sensDot;

	std::vector<std::size_t> _bulkLayout;
	std::vector<std::vector<std::size_t>> _particleLayout;
	std::vector<std::vector<std::size_t>> _solidLayout;
	std::vector<std::size_t> _fluxLayout;

	unsigned int _nComp;
	unsigned int _nVolumeDof;
	unsigned int _nInletPorts;
	unsigned int _nOutletPorts;
	std::vector<unsigned int> _nParShells;
	std::vector<unsigned int> _nBoundStates;
	unsigned int _numTimesteps;
	unsigned int _numSens;
	UnitOpIdx _unitOp;

	bool _needsReAlloc;

	unsigned int _bulkCount; //!< Number of bulk mobile phase DOF blocks
	std::vector<unsigned int> _particleCount; //!< Number of particle mobile phase DOF blocks per particle type
	std::vector<unsigned int> _solidCount; //!< Number of solid phase DOF blocks per particle type

	std::vector<double> _axialCoords;
	std::vector<double> _radialCoords;
	std::vector<double> _particleCoords;
};


/**
 * @brief Stores pieces of the solution of the whole model system in recorders of single unit operations
 * @details Maintains a collection of InternalStorageUnitOpRecorder objects that store individual unit operations.
 *          The individual unit operation recorders are owned by this object and destroyed upon its own
 *          destruction.
 */
class InternalStorageSystemRecorder : public ISolutionRecorder
{
public:

	InternalStorageSystemRecorder() : _numTimesteps(0), _numSens(0), _storeTime(true)
	{
	}

	virtual ~InternalStorageSystemRecorder() CADET_NOEXCEPT
	{
		deleteRecorders();
	}

	virtual void clear()
	{
		_time.clear();

		for (InternalStorageUnitOpRecorder* rec : _recorders)
			rec->clear();
	}

	virtual void prepare(unsigned int numDofs, unsigned int numSens, unsigned int numTimesteps)
	{
		_numSens = numSens;
		if (numTimesteps > 0)
			_time.reserve(numTimesteps);
		else
			_time.reserve(detail::numDefaultRecorderTimesteps);

		for (InternalStorageUnitOpRecorder* rec : _recorders)
			rec->prepare(numDofs, numSens, numTimesteps);
	}

	virtual void notifyIntegrationStart(unsigned int numDofs, unsigned int numSens, unsigned int numTimesteps)
	{
		_numSens = numSens;
		_time.clear();

		if (numTimesteps > 0)
			_time.reserve(numTimesteps);
		else
			_time.reserve(detail::numDefaultRecorderTimesteps);

		for (InternalStorageUnitOpRecorder* rec : _recorders)
			rec->notifyIntegrationStart(numDofs, numSens, numTimesteps);
	}

	virtual void unitOperationStructure(UnitOpIdx idx, const IModel& model, const ISolutionExporter& exporter)
	{
		for (InternalStorageUnitOpRecorder* rec : _recorders)
			rec->unitOperationStructure(idx, model, exporter);

		// Reset for counting actual number of time steps
		_numTimesteps = 0;
	}

	virtual void beginTimestep(double t)
	{
		++_numTimesteps;
		if (_storeTime)
			_time.push_back(t);

		for (InternalStorageUnitOpRecorder* rec : _recorders)
			rec->beginTimestep(t);
	}

	virtual void beginUnitOperation(cadet::UnitOpIdx idx, const cadet::IModel& model, const cadet::ISolutionExporter& exporter)
	{
		for (InternalStorageUnitOpRecorder* rec : _recorders)
			rec->beginUnitOperation(idx, model, exporter);
	}

	virtual void endUnitOperation()
	{
		for (InternalStorageUnitOpRecorder* rec : _recorders)
			rec->endUnitOperation();
	}

	virtual void endTimestep()
	{
		for (InternalStorageUnitOpRecorder* rec : _recorders)
			rec->endTimestep();
	}

	virtual void beginSolution()
	{
		for (InternalStorageUnitOpRecorder* rec : _recorders)
			rec->beginSolution();
	}

	virtual void endSolution()
	{
		for (InternalStorageUnitOpRecorder* rec : _recorders)
			rec->endSolution();
	}

	virtual void beginSolutionDerivative()
	{
		for (InternalStorageUnitOpRecorder* rec : _recorders)
			rec->beginSolutionDerivative();
	}

	virtual void endSolutionDerivative()
	{
		for (InternalStorageUnitOpRecorder* rec : _recorders)
			rec->endSolutionDerivative();
	}

	virtual void beginSensitivity(const cadet::ParameterId& pId, unsigned int sensIdx)
	{
		for (InternalStorageUnitOpRecorder* rec : _recorders)
			rec->beginSensitivity(pId, sensIdx);
	}

	virtual void endSensitivity(const cadet::ParameterId& pId, unsigned int sensIdx)
	{
		for (InternalStorageUnitOpRecorder* rec : _recorders)
			rec->endSensitivity(pId, sensIdx);
	}

	virtual void beginSensitivityDerivative(const cadet::ParameterId& pId, unsigned int sensIdx)
	{
		for (InternalStorageUnitOpRecorder* rec : _recorders)
			rec->beginSensitivityDerivative(pId, sensIdx);
	}

	virtual void endSensitivityDerivative(const cadet::ParameterId& pId, unsigned int sensIdx)
	{
		for (InternalStorageUnitOpRecorder* rec : _recorders)
			rec->endSensitivityDerivative(pId, sensIdx);
	}

	template <typename Writer_t>
	void writeCoordinates(Writer_t& writer)
	{
		std::ostringstream oss;

		for (InternalStorageUnitOpRecorder* rec : _recorders)
		{
			oss.str("");
			oss << "unit_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << static_cast<int>(rec->unitOperation());

			if (rec->storeCoordinates())
			{
				writer.pushGroup(oss.str());
				rec->writeCoordinates(writer);
				writer.popGroup();
			}
		}
	}
	
	template <typename Writer_t>
	void writeSolution(Writer_t& writer)
	{
		std::ostringstream oss;

		if (_storeTime)
			writer.template vector<double>("SOLUTION_TIMES", _time.size(), _time.data());

		for (InternalStorageUnitOpRecorder* rec : _recorders)
		{
			oss.str("");
			oss << "unit_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << static_cast<int>(rec->unitOperation());

			writer.pushGroup(oss.str());
			rec->writeSolution(writer);
			writer.popGroup();
		}
	}

	template <typename Writer_t>
	void writeSensitivity(Writer_t& writer)
	{
		std::ostringstream oss;

		for (unsigned int param = 0; param < _numSens; ++param)
		{
			oss.str("");
			oss << "param_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << param;
			writer.pushGroup(oss.str());

			for (InternalStorageUnitOpRecorder* rec : _recorders)
			{
				oss.str("");
				oss << "unit_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << static_cast<int>(rec->unitOperation());

				writer.pushGroup(oss.str());
				rec->writeSensitivity(writer, param);
				writer.popGroup();
			}

			writer.popGroup();
		}
	}

	inline bool storeTime() const CADET_NOEXCEPT { return _storeTime; }
	inline void storeTime(bool st) CADET_NOEXCEPT { _storeTime = st; }

	inline bool anyUnitStoresCoordinates() const CADET_NOEXCEPT
	{
		for (InternalStorageUnitOpRecorder const* rec : _recorders)
		{
			if (rec->storeCoordinates())
				return true;
		}

		return false;
	}

	inline unsigned int numDataPoints() const CADET_NOEXCEPT { return _numTimesteps; }

	inline void addRecorder(InternalStorageUnitOpRecorder* rec)
	{
		_recorders.push_back(rec);
	}

	inline unsigned int numRecorders() const CADET_NOEXCEPT { return _recorders.size(); }
	inline InternalStorageUnitOpRecorder* recorder(unsigned int idx) CADET_NOEXCEPT { return _recorders[idx]; }
	inline InternalStorageUnitOpRecorder* recorder(unsigned int idx) const CADET_NOEXCEPT { return _recorders[idx]; }

	inline InternalStorageUnitOpRecorder* unitOperation(UnitOpIdx idx) CADET_NOEXCEPT
	{
		for (InternalStorageUnitOpRecorder* rec : _recorders)
		{
			if (rec->unitOperation() == idx)
				return rec;
		}
		return nullptr;
	}
	inline InternalStorageUnitOpRecorder const* unitOperation(UnitOpIdx idx) const CADET_NOEXCEPT
	{
		for (InternalStorageUnitOpRecorder const* rec : _recorders)
		{
			if (rec->unitOperation() == idx)
				return rec;
		}
		return nullptr;
	}

	inline void deleteRecorders()
	{
		for (InternalStorageUnitOpRecorder* rec : _recorders)
			delete rec;
		_recorders.clear();
	}

	inline double const* time() const CADET_NOEXCEPT { return _time.data(); }

protected:

	std::vector<InternalStorageUnitOpRecorder*> _recorders;
	unsigned int _numTimesteps;
	unsigned int _numSens;
	std::vector<double> _time;
	bool _storeTime;
};


} // namespace cadet

#endif  // LIBCADET_SOLUTIONRECORDER_IMPL_HPP_
