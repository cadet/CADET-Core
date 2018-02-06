// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2017: The CADET Authors
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

#include "cadet/SolutionRecorder.hpp"

namespace cadet
{

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
	};

	InternalStorageUnitOpRecorder() : InternalStorageUnitOpRecorder(UnitOpIndep) { }

	InternalStorageUnitOpRecorder(UnitOpIdx idx) : _cfgSolution({false, false, false, true, false, false}),
		_cfgSolutionDot({false, false, false, false, false, false}), _cfgSensitivity({false, false, false, true, false, false}),
		_cfgSensitivityDot({false, false, false, true, false, false}), _storeTime(false), _splitComponents(true), _curCfg(nullptr),
		_nComp(0), _nVolumeDof(0), _numTimesteps(0), _numSens(0), _unitOp(idx), _needsReAlloc(false)
	{
	}

	virtual ~InternalStorageUnitOpRecorder() CADET_NOEXCEPT
	{
		for (unsigned int i = 0; i < _sensOutlet.size(); ++i)
		{
			delete _sensOutlet[i];
			delete _sensInlet[i];
			delete _sensBulk[i];
			delete _sensParticle[i];
			delete _sensSolid[i];
			delete _sensFlux[i];
			delete _sensVolume[i];

			delete _sensOutletDot[i];
			delete _sensInletDot[i];
			delete _sensBulkDot[i];
			delete _sensParticleDot[i];
			delete _sensSolidDot[i];
			delete _sensFluxDot[i];
			delete _sensVolumeDot[i];
		}
	}

	virtual void clear()
	{
		// Clear solution storage
		_time.clear();
		_outlet.clear();
		_inlet.clear();
		_bulk.clear();
		_particle.clear();
		_solid.clear();
		_flux.clear();
		_volume.clear();

		_outletDot.clear();
		_inletDot.clear();
		_bulkDot.clear();
		_particleDot.clear();
		_solidDot.clear();
		_fluxDot.clear();
		_volumeDot.clear();

		// Clear all sensitivity storage
		for (unsigned int i = 0; i < _sensOutlet.size(); ++i)
		{
			_sensOutlet[i]->clear();
			_sensInlet[i]->clear();
			_sensBulk[i]->clear();
			_sensParticle[i]->clear();
			_sensSolid[i]->clear();
			_sensFlux[i]->clear();
			_sensVolume[i]->clear();

			_sensOutletDot[i]->clear();
			_sensInletDot[i]->clear();
			_sensBulkDot[i]->clear();
			_sensParticleDot[i]->clear();
			_sensSolidDot[i]->clear();
			_sensFluxDot[i]->clear();
			_sensVolumeDot[i]->clear();
		}
	}

	virtual void prepare(unsigned int numDofs, unsigned int numSens, unsigned int numTimesteps)
	{
		_numTimesteps = numTimesteps;
		_numSens = numSens;

		// Allocate sensitivity storage
		_sensOutlet.resize(numSens, nullptr);
		_sensInlet.resize(numSens, nullptr);
		_sensBulk.resize(numSens, nullptr);
		_sensParticle.resize(numSens, nullptr);
		_sensSolid.resize(numSens, nullptr);
		_sensFlux.resize(numSens, nullptr);
		_sensVolume.resize(numSens, nullptr);

		_sensOutletDot.resize(numSens, nullptr);
		_sensInletDot.resize(numSens, nullptr);
		_sensBulkDot.resize(numSens, nullptr);
		_sensParticleDot.resize(numSens, nullptr);
		_sensSolidDot.resize(numSens, nullptr);
		_sensFluxDot.resize(numSens, nullptr);
		_sensVolumeDot.resize(numSens, nullptr);

		for (unsigned int i = 0; i < numSens; ++i)
		{
			_sensOutlet[i] = new std::vector<double>();
			_sensInlet[i] = new std::vector<double>();
			_sensBulk[i] = new std::vector<double>();
			_sensParticle[i] = new std::vector<double>();
			_sensSolid[i] = new std::vector<double>();
			_sensFlux[i] = new std::vector<double>();
			_sensVolume[i] = new std::vector<double>();

			_sensOutletDot[i] = new std::vector<double>();
			_sensInletDot[i] = new std::vector<double>();
			_sensBulkDot[i] = new std::vector<double>();
			_sensParticleDot[i] = new std::vector<double>();
			_sensSolidDot[i] = new std::vector<double>();
			_sensFluxDot[i] = new std::vector<double>();
			_sensVolumeDot[i] = new std::vector<double>();
		}

		_needsReAlloc = false;
	}

	virtual void notifyIntegrationStart(unsigned int numDofs, unsigned int numSens, unsigned int numTimesteps)
	{
		_needsReAlloc = (numSens != _numSens) || (numTimesteps > _numTimesteps);

		// Clear all data from memory
		clear();		

		_numTimesteps = numTimesteps;
		
		if (numSens != _numSens)
		{
			// Delete all sensitivity storage
			for (unsigned int i = 0; i < _sensOutlet.size(); ++i)
			{
				delete _sensOutlet[i];
				delete _sensInlet[i];
				delete _sensBulk[i];
				delete _sensParticle[i];
				delete _sensSolid[i];
				delete _sensFlux[i];
				delete _sensVolume[i];

				delete _sensOutletDot[i];
				delete _sensInletDot[i];
				delete _sensBulkDot[i];
				delete _sensParticleDot[i];
				delete _sensSolidDot[i];
				delete _sensFluxDot[i];
				delete _sensVolumeDot[i];
			}

			// Allocate sensitivity storage
			_sensOutlet.clear();
			_sensOutlet.resize(numSens, nullptr);
			_sensInlet.clear();
			_sensInlet.resize(numSens, nullptr);
			_sensBulk.clear();
			_sensBulk.resize(numSens, nullptr);
			_sensParticle.clear();
			_sensParticle.resize(numSens, nullptr);
			_sensSolid.clear();
			_sensSolid.resize(numSens, nullptr);
			_sensFlux.clear();
			_sensFlux.resize(numSens, nullptr);
			_sensVolume.clear();
			_sensVolume.resize(numSens, nullptr);

			_sensOutletDot.clear();
			_sensOutletDot.resize(numSens, nullptr);
			_sensInletDot.clear();
			_sensInletDot.resize(numSens, nullptr);
			_sensBulkDot.clear();
			_sensBulkDot.resize(numSens, nullptr);
			_sensParticleDot.clear();
			_sensParticleDot.resize(numSens, nullptr);
			_sensSolidDot.clear();
			_sensSolidDot.resize(numSens, nullptr);
			_sensFluxDot.clear();
			_sensFluxDot.resize(numSens, nullptr);
			_sensVolumeDot.clear();
			_sensVolumeDot.resize(numSens, nullptr);

			// Populate with empty vectors
			for (unsigned int i = 0; i < numSens; ++i)
			{
				_sensOutlet[i] = new std::vector<double>();
				_sensInlet[i] = new std::vector<double>();
				_sensBulk[i] = new std::vector<double>();
				_sensParticle[i] = new std::vector<double>();
				_sensSolid[i] = new std::vector<double>();
				_sensFlux[i] = new std::vector<double>();
				_sensVolume[i] = new std::vector<double>();

				_sensOutletDot[i] = new std::vector<double>();
				_sensInletDot[i] = new std::vector<double>();
				_sensBulkDot[i] = new std::vector<double>();
				_sensParticleDot[i] = new std::vector<double>();
				_sensSolidDot[i] = new std::vector<double>();
				_sensFluxDot[i] = new std::vector<double>();
				_sensVolumeDot[i] = new std::vector<double>();
			}

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
				case StateOrdering::BoundState:
					break;
			}
		}

		order = exporter.mobilePhaseOrdering(len);
		_particleLayout.clear();
		_particleLayout.reserve(len + 1); // First slot is time
		_particleLayout.push_back(0);
		_particleCount = 1;
		for (unsigned int i = 0; i < len; ++i)
		{
			switch (order[i])
			{
				case StateOrdering::Component:
					_particleLayout.push_back(exporter.numComponents());
					break;
				case StateOrdering::AxialCell:
					_particleLayout.push_back(exporter.numAxialCells());
					_particleCount *= exporter.numAxialCells();
					break;
				case StateOrdering::RadialCell:
					_particleLayout.push_back(exporter.numRadialCells());
					_particleCount *= exporter.numRadialCells();
					 break;
				case StateOrdering::BoundState:
					break;
			}
		}

		order = exporter.solidPhaseOrdering(len);
		_solidLayout.clear();
		_solidLayout.reserve(len + 1); // First slot is time
		_solidLayout.push_back(0);
		_solidCount = 1;
		for (unsigned int i = 0; i < len; ++i)
		{
			switch (order[i])
			{
				case StateOrdering::Component:
					break;
				case StateOrdering::AxialCell:
					_solidLayout.push_back(exporter.numAxialCells());
					_solidCount *= exporter.numAxialCells();
					break;
				case StateOrdering::RadialCell:
					_solidLayout.push_back(exporter.numRadialCells());
					_solidCount *= exporter.numRadialCells();
					 break;
				case StateOrdering::BoundState:
					_solidLayout.push_back(exporter.numBoundStates());
					break;
			}
		}

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
				case StateOrdering::RadialCell:
				case StateOrdering::BoundState:
					break;
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
		for (unsigned int i = 0; i < _sensOutlet.size(); ++i)
		{
			beginSensitivity(i);
			allocateMemory(exporter);
			endSolution();

			beginSensitivityDot(i);
			allocateMemory(exporter);
			endSolution();
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
			double const* outlet = exporter.outlet(stride);
			for (unsigned int i = 0; i < _nComp; ++i)
				_curOutlet->push_back(outlet[i * stride]);
		}

		if (_curCfg->storeInlet)
		{
			double const* inlet = exporter.inlet(stride);
			for (unsigned int i = 0; i < _nComp; ++i)
				_curInlet->push_back(inlet[i * stride]);
		}

		if (_curCfg->storeBulk)
		{
			double const* data = exporter.concentration();
			stride = exporter.bulkMobilePhaseStride();
			const unsigned int blockSize = exporter.numBulkDofs() / _bulkCount;
			for (unsigned int i = 0; i < _bulkCount; ++i, data += stride)
			{
				_curBulk->insert(_curBulk->end(), data, data + blockSize);
			}
		}

		if (_curCfg->storeParticle)
		{
			double const* data = exporter.mobilePhase();
			stride = exporter.particleMobilePhaseStride();
			const unsigned int blockSize = exporter.numParticleMobilePhaseDofs() / _particleCount;
			for (unsigned int i = 0; i < _particleCount; ++i, data += stride)
			{
				_curParticle->insert(_curParticle->end(), data, data + blockSize);
			}
		}

		if (_curCfg->storeSolid)
		{
			double const* data = exporter.solidPhase();
			stride = exporter.solidPhaseStride();
			const unsigned int blockSize = exporter.numSolidPhaseDofs() / _solidCount;
			for (unsigned int i = 0; i < _solidCount; ++i, data += stride)
			{
				_curSolid->insert(_curSolid->end(), data, data + blockSize);
			}
		}

		if (_curCfg->storeFlux)
		{
			double const* const data = exporter.flux();
			_curFlux->insert(_curFlux->end(), data, data + exporter.numFluxDofs());
		}

		if (_curCfg->storeVolume)
		{
			double const* const data = exporter.volume();
			_curVolume->insert(_curVolume->end(), data, data + exporter.numVolumeDofs());
		}
	}

	virtual void endUnitOperation() { }

	virtual void endTimestep() { }

	virtual void beginSolution()
	{
		_curCfg = &_cfgSolution;
		_curOutlet = &_outlet;
		_curInlet = &_inlet;
		_curBulk = &_bulk;
		_curParticle = &_particle;
		_curSolid = &_solid;
		_curFlux = &_flux;
		_curVolume = &_volume;
	}

	virtual void endSolution()
	{
		_curCfg = nullptr;
		_curOutlet = nullptr;
		_curInlet = nullptr;
		_curBulk = nullptr;
		_curParticle = nullptr;
		_curSolid = nullptr;
		_curFlux = nullptr;
		_curVolume = nullptr;
	}

	virtual void beginSolutionDerivative()
	{
		_curCfg = &_cfgSolutionDot;
		_curOutlet = &_outletDot;
		_curInlet = &_inletDot;
		_curBulk = &_bulkDot;
		_curParticle = &_particleDot;
		_curSolid = &_solidDot;
		_curFlux = &_fluxDot;
		_curVolume = &_volumeDot;
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

	inline bool splitComponents() const CADET_NOEXCEPT { return _splitComponents; }
	inline void splitComponents(bool st) CADET_NOEXCEPT { _splitComponents = st; }

	inline UnitOpIdx unitOperation() const CADET_NOEXCEPT { return _unitOp; }
	inline void unitOperation(UnitOpIdx idx) CADET_NOEXCEPT { _unitOp = idx; }

	inline unsigned int numDataPoints() const CADET_NOEXCEPT { return _numTimesteps; }
	inline unsigned int numComponents() const CADET_NOEXCEPT { return _nComp; }

	inline double const* time() const CADET_NOEXCEPT { return _time.data(); }
	inline double const* inlet() const CADET_NOEXCEPT { return _inlet.data(); }
	inline double const* outlet() const CADET_NOEXCEPT { return _outlet.data(); }
	inline double const* bulk() const CADET_NOEXCEPT { return _bulk.data(); }
	inline double const* particle() const CADET_NOEXCEPT { return _particle.data(); }
	inline double const* solid() const CADET_NOEXCEPT { return _solid.data(); }
	inline double const* flux() const CADET_NOEXCEPT { return _flux.data(); }
	inline double const* volume() const CADET_NOEXCEPT { return _volume.data(); }
	inline double const* inletDot() const CADET_NOEXCEPT { return _inletDot.data(); }
	inline double const* outletDot() const CADET_NOEXCEPT { return _outletDot.data(); }
	inline double const* bulkDot() const CADET_NOEXCEPT { return _bulkDot.data(); }
	inline double const* particleDot() const CADET_NOEXCEPT { return _particleDot.data(); }
	inline double const* solidDot() const CADET_NOEXCEPT { return _solidDot.data(); }
	inline double const* fluxDot() const CADET_NOEXCEPT { return _fluxDot.data(); }
	inline double const* volumeDot() const CADET_NOEXCEPT { return _volumeDot.data(); }
	inline double const* sensInlet(unsigned int idx) const CADET_NOEXCEPT { return _sensInlet[idx]->data(); }
	inline double const* sensOutlet(unsigned int idx) const CADET_NOEXCEPT { return _sensOutlet[idx]->data(); }
	inline double const* sensBulk(unsigned int idx) const CADET_NOEXCEPT { return _sensBulk[idx]->data(); }
	inline double const* sensParticle(unsigned int idx) const CADET_NOEXCEPT { return _sensParticle[idx]->data(); }
	inline double const* sensSolid(unsigned int idx) const CADET_NOEXCEPT { return _sensSolid[idx]->data(); }
	inline double const* sensFlux(unsigned int idx) const CADET_NOEXCEPT { return _sensFlux[idx]->data(); }
	inline double const* sensVolume(unsigned int idx) const CADET_NOEXCEPT { return _sensVolume[idx]->data(); }
	inline double const* sensInletDot(unsigned int idx) const CADET_NOEXCEPT { return _sensInletDot[idx]->data(); }
	inline double const* sensOutletDot(unsigned int idx) const CADET_NOEXCEPT { return _sensOutletDot[idx]->data(); }
	inline double const* sensBulkDot(unsigned int idx) const CADET_NOEXCEPT { return _sensBulkDot[idx]->data(); }
	inline double const* sensParticleDot(unsigned int idx) const CADET_NOEXCEPT { return _sensParticleDot[idx]->data(); }
	inline double const* sensSolidDot(unsigned int idx) const CADET_NOEXCEPT { return _sensSolidDot[idx]->data(); }
	inline double const* sensFluxDot(unsigned int idx) const CADET_NOEXCEPT { return _sensFluxDot[idx]->data(); }
	inline double const* sensVolumeDot(unsigned int idx) const CADET_NOEXCEPT { return _sensVolumeDot[idx]->data(); }
protected:

	inline void beginSensitivity(unsigned int sensIdx)
	{
		_curCfg = &_cfgSensitivity;
		_curOutlet = _sensOutlet[sensIdx];
		_curInlet = _sensInlet[sensIdx];
		_curBulk = _sensBulk[sensIdx];
		_curParticle = _sensParticle[sensIdx];
		_curSolid = _sensSolid[sensIdx];
		_curFlux = _sensFlux[sensIdx];
		_curVolume = _sensVolume[sensIdx];
	}

	inline void beginSensitivityDot(unsigned int sensIdx)
	{
		_curCfg = &_cfgSensitivityDot;
		_curOutlet = _sensOutletDot[sensIdx];
		_curInlet = _sensInletDot[sensIdx];
		_curBulk = _sensBulkDot[sensIdx];
		_curParticle = _sensParticleDot[sensIdx];
		_curSolid = _sensSolidDot[sensIdx];
		_curFlux = _sensFluxDot[sensIdx];
		_curVolume = _sensVolumeDot[sensIdx];
	}

	inline void validateConfig(const ISolutionExporter& exporter, StorageConfig& cfg)
	{
		// Only store fields that really exist
		cfg.storeParticle = exporter.hasParticleMobilePhase() && cfg.storeParticle;
		cfg.storeSolid = exporter.hasSolidPhase() && cfg.storeSolid;
		cfg.storeFlux = exporter.hasParticleFlux() && cfg.storeFlux;
		cfg.storeVolume = exporter.hasVolume() && cfg.storeVolume;
	}

	inline void allocateMemory(const ISolutionExporter& exporter)
	{
		if (_curCfg->storeOutlet)
			_curOutlet->reserve(std::max(_numTimesteps, 100u) * _nComp);

		if (_curCfg->storeInlet)
			_curInlet->reserve(std::max(_numTimesteps, 100u) * _nComp);

		if (_curCfg->storeBulk)
			_curBulk->reserve(std::max(_numTimesteps, 100u) * exporter.numBulkDofs());
		
		if (_curCfg->storeParticle)
			_curParticle->reserve(std::max(_numTimesteps, 100u) * exporter.numParticleMobilePhaseDofs());
		
		if (_curCfg->storeSolid)
			_curSolid->reserve(std::max(_numTimesteps, 100u) * exporter.numSolidPhaseDofs());

		if (_curCfg->storeFlux)
			_curFlux->reserve(std::max(_numTimesteps, 100u) * exporter.numFluxDofs());
		
		if (_curCfg->storeVolume)
			_curVolume->reserve(std::max(_numTimesteps, 100u) * exporter.numVolumeDofs());
	}

	template <typename Writer_t>
	void writeData(Writer_t& writer, const char* prefix, std::ostringstream& oss)
	{
		if (_curCfg->storeOutlet)
		{
			if (_splitComponents)
			{
				for (unsigned int comp = 0; comp < _nComp; ++comp)
				{
					oss.str("");
					oss << prefix << "_OUTLET_COMP_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << comp;
					writer.template vector<double>(oss.str(), _numTimesteps, _curOutlet->data() + comp, _nComp);
				}
			}
			else
			{
				oss.str("");
				oss << prefix << "_OUTLET";
				writer.template matrix<double>(oss.str(), _numTimesteps, _nComp, _curOutlet->data(), 1);
			}
		}

		if (_curCfg->storeInlet)
		{
			if (_splitComponents)
			{
				for (unsigned int comp = 0; comp < _nComp; ++comp)
				{
					oss.str("");
					oss << prefix << "_INLET_COMP_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << comp;
					writer.template vector<double>(oss.str(), _numTimesteps, _curInlet->data() + comp, _nComp);
				}
			}
			else
			{
				oss.str("");
				oss << prefix << "_INLET";
				writer.template matrix<double>(oss.str(), _numTimesteps, _nComp, _curInlet->data(), 1);
			}
		}

		if (_curCfg->storeBulk)
		{
			oss.str("");
			oss << prefix << "_BULK";
			_bulkLayout[0] = _numTimesteps;
			writer.template tensor<double>(oss.str(), _bulkLayout.size(), _bulkLayout.data(), _curBulk->data());
		}

		if (_curCfg->storeParticle)
		{
			oss.str("");
			oss << prefix << "_PARTICLE";
			_particleLayout[0] = _numTimesteps;
			writer.template tensor<double>(oss.str(), _particleLayout.size(), _particleLayout.data(), _curParticle->data());
		}

		if (_curCfg->storeSolid)
		{
			oss.str("");
			oss << prefix << "_SOLID";
			_solidLayout[0] = _numTimesteps;
			writer.template tensor<double>(oss.str(), _solidLayout.size(), _solidLayout.data(), _curSolid->data());
		}

		if (_curCfg->storeFlux)
		{
			oss.str("");
			oss << prefix << "_FLUX";
			_fluxLayout[0] = _numTimesteps;
			writer.template tensor<double>(oss.str(), _fluxLayout.size(), _fluxLayout.data(), _curFlux->data());
		}

		if (_curCfg->storeVolume)
		{
			oss.str("");
			oss << prefix << "_VOLUME";
			_fluxLayout[0] = _numTimesteps;
			writer.template matrix<double>(oss.str(), _numTimesteps, _nVolumeDof, _curVolume->data(), 1);
		}
	}

	StorageConfig _cfgSolution;
	StorageConfig _cfgSolutionDot;
	StorageConfig _cfgSensitivity;
	StorageConfig _cfgSensitivityDot;
	bool _storeTime;
	bool _splitComponents;

	StorageConfig const* _curCfg;
	std::vector<double>* _curOutlet;
	std::vector<double>* _curInlet;
	std::vector<double>* _curBulk;
	std::vector<double>* _curParticle;
	std::vector<double>* _curSolid;
	std::vector<double>* _curFlux;
	std::vector<double>* _curVolume;

	std::vector<double> _time;
	std::vector<double> _outlet;
	std::vector<double> _inlet;
	std::vector<double> _bulk;
	std::vector<double> _particle;
	std::vector<double> _solid;
	std::vector<double> _flux;
	std::vector<double> _volume;

	std::vector<double> _outletDot;
	std::vector<double> _inletDot;
	std::vector<double> _bulkDot;
	std::vector<double> _particleDot;
	std::vector<double> _solidDot;
	std::vector<double> _fluxDot;
	std::vector<double> _volumeDot;

	std::vector<std::vector<double>*> _sensOutlet;
	std::vector<std::vector<double>*> _sensInlet;
	std::vector<std::vector<double>*> _sensBulk;
	std::vector<std::vector<double>*> _sensParticle;
	std::vector<std::vector<double>*> _sensSolid;
	std::vector<std::vector<double>*> _sensFlux;
	std::vector<std::vector<double>*> _sensVolume;

	std::vector<std::vector<double>*> _sensOutletDot;
	std::vector<std::vector<double>*> _sensInletDot;
	std::vector<std::vector<double>*> _sensBulkDot;
	std::vector<std::vector<double>*> _sensParticleDot;
	std::vector<std::vector<double>*> _sensSolidDot;
	std::vector<std::vector<double>*> _sensFluxDot;
	std::vector<std::vector<double>*> _sensVolumeDot;
	
	std::vector<std::size_t> _bulkLayout;
	std::vector<std::size_t> _particleLayout;
	std::vector<std::size_t> _solidLayout;
	std::vector<std::size_t> _fluxLayout;

	unsigned int _nComp;
	unsigned int _nVolumeDof;
	unsigned int _numTimesteps;
	unsigned int _numSens;
	UnitOpIdx _unitOp;

	bool _needsReAlloc;

	unsigned int _bulkCount; //!< Number of bulk mobile phase DOF blocks
	unsigned int _particleCount; //!< Number of particle mobile phase DOF blocks
	unsigned int _solidCount; //!< Number of solid phase DOF blocks
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
		_time.reserve(numTimesteps);

		for (InternalStorageUnitOpRecorder* rec : _recorders)
			rec->prepare(numDofs, numSens, numTimesteps);
	}

	virtual void notifyIntegrationStart(unsigned int numDofs, unsigned int numSens, unsigned int numTimesteps)
	{
		_numSens = numSens;
		_time.clear();
		_time.reserve(numTimesteps);

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

	inline unsigned int numDataPoints() const CADET_NOEXCEPT { return _numTimesteps; }

	inline void addRecorder(InternalStorageUnitOpRecorder* rec)
	{
		_recorders.push_back(rec);
	}

	inline unsigned int numRecorders() const CADET_NOEXCEPT { return _recorders.size(); }
	inline InternalStorageUnitOpRecorder* recorder(unsigned int idx) CADET_NOEXCEPT { return _recorders[idx]; }
	inline InternalStorageUnitOpRecorder* const recorder(unsigned int idx) const CADET_NOEXCEPT { return _recorders[idx]; }

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
