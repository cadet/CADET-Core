// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2015: Eric von Lieres¹, Joel Andersson,
//                         Andreas Puettmann¹, Sebastian Schnittert¹,
//                         Samuel Leweke¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include <cmath>
#include <vector>
#include <algorithm>

#include "ParticleDiscretization.hpp"
#include "CadetLogger.hpp"

namespace cadet
{
ParticleDiscretization::ParticleDiscretization(int npar) :
    _npar                           (npar),
    _par_cell_size                  (std::vector<double>(_npar, 0.0)),
    _par_center_radius              (std::vector<double>(_npar, 0.0)),
    _par_outer_surf_area_per_volume (std::vector<double>(_npar, 0.0)),
    _par_inner_surf_area_per_volume (std::vector<double>(_npar, 0.0))
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    setEquidistantRadialDisc();  // Use as default

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}

void ParticleDiscretization::setEquidistantRadialDisc()
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    const double dr = 1.0 / _npar;
    _par_cell_size.assign(_npar, dr);

    for (int cell = 0; cell < _npar; cell++)
    {
        _par_center_radius.at(cell) = 1.0 - (0.5 + static_cast<double>(cell)) * dr;

        // compute denominator -> corresponding to cell volume
        double vol = pow(1.0 - static_cast<double>(cell) * dr, 3) - pow(1.0 - static_cast<double>(cell + 1) * dr, 3);

        _par_outer_surf_area_per_volume.at(cell) = 3.0 * pow(1.0 - static_cast<double>(cell) * dr, 2) / vol;
        _par_inner_surf_area_per_volume.at(cell) = 3.0 * pow(1.0 - static_cast<double>(cell + 1) * dr, 2) / vol;
    }

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}

void ParticleDiscretization::setEquivolumeRadialDisc()
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    double r_out = 1.0;
    double r_in = 0.0;

    for (int cell = 0; cell < _npar; ++cell)
    {
        if (cell != (_npar - 1)) r_in = pow(pow(r_out, 3) - 1.0 / _npar, (1.0 / 3.0));

        _par_cell_size.at(cell) = r_out - r_in;

        _par_center_radius.at(cell) = r_out - 0.5 * _par_cell_size.at(cell);

        double vol = 1.0 / _npar;

        _par_outer_surf_area_per_volume.at(cell) = 3.0 * pow(r_out, 2) / vol;
        _par_inner_surf_area_per_volume.at(cell) = 3.0 * pow(r_in, 2) / vol;

        // for the next cell: r_out==r_in of the current cell
        r_out = r_in;
    }

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}

void ParticleDiscretization::setUserdefinedRadialDisc(const std::vector<double>& cell_interfaces)
{
    log::emit<Trace1>() << CURRENT_FUNCTION << Color::cyan << ": Called!" << Color::reset << log::endl;

    // Care for the right ordering and include 0.0 / 1.0 if not already in the vector.
    std::vector<double> orderedInterfaces = cell_interfaces;

    if (std::find(orderedInterfaces.begin(), orderedInterfaces.end(), 0.0) == orderedInterfaces.end())
        orderedInterfaces.push_back(0.0);
    if (std::find(orderedInterfaces.begin(), orderedInterfaces.end(), 1.0) == orderedInterfaces.end())
        orderedInterfaces.push_back(1.0);

    std::sort(orderedInterfaces.begin(), orderedInterfaces.end());
    std::reverse(orderedInterfaces.begin(), orderedInterfaces.end());

    for (int cell = 0; cell < _npar; ++cell)
    {
        _par_cell_size.at(cell) = orderedInterfaces.at(cell) - orderedInterfaces.at(cell + 1);

        _par_center_radius.at(cell) = orderedInterfaces.at(cell) - 0.5 * _par_cell_size.at(cell);

        // compute denominator -> corresponding to cell volume
        double vol = pow(orderedInterfaces.at(cell), 3) - pow(orderedInterfaces.at(cell + 1), 3);

        _par_outer_surf_area_per_volume.at(cell) = 3.0 * pow(orderedInterfaces.at(cell), 2) / vol;
        _par_inner_surf_area_per_volume.at(cell) = 3.0 * pow(orderedInterfaces.at(cell + 1), 2) / vol;
    }

    log::emit<Trace1>() << CURRENT_FUNCTION << Color::green << ": Finished!" << Color::reset << log::endl;
}

}  // namespace cadet


