// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2016: Eric von Lieres¹, Joel Andersson¹,
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

#ifndef PARTICLEDISCRETIZATION_HPP_
#define PARTICLEDISCRETIZATION_HPP_

namespace cadet
{

class ParticleDiscretization
{
public:
    // Constructor
    ParticleDiscretization(int npar);

    inline const std::vector<double> & getParCellCoords(void)   const { return _par_center_radius; }
    inline double getCellSize(int cell)                         const { return _par_cell_size.at(cell); }
    inline double getCenterRadius(int cell)                     const { return _par_center_radius.at(cell); }
    inline double getOuterSurfAreaPerVolume(int cell)           const { return _par_outer_surf_area_per_volume.at(cell); }
    inline double getInnerSurfAreaPerVolume(int cell)           const { return _par_inner_surf_area_per_volume.at(cell); }

    void setEquidistantRadialDisc();
    void setEquivolumeRadialDisc();
    void setUserdefinedRadialDisc(const std::vector<double> & cell_interfaces);

private:
    int _npar; //!< Number of particle (radial) discretization cells

    std::vector<double> _par_cell_size; //!< Cell size
    std::vector<double> _par_center_radius; //!< Particle cell-centered position for each particle cell
    std::vector<double> _par_outer_surf_area_per_volume;
    std::vector<double> _par_inner_surf_area_per_volume;

}; // class ParticleDiscretization

} // namespace cadet


#endif // PARTICLEDISCRETIZATION_HPP_
