// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2015: Eric von Lieres¹, Joel Andersson¹,
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

#ifndef SPARSEMATRIXELEMENT_HPP_
#define SPARSEMATRIXELEMENT_HPP_

namespace cadet
{

class SparseMatrixElement
{
public:

    inline void row(int r) { _row = r; }
    inline int  row() const { return _row; }

    inline void col(int c) { _col = c; }
    inline int  col() const { return _col; }

    inline void   value(double v) { _value = v; }
    inline double value() const { return _value; }

    inline void setElement(int row, int col, double value)
    {
        _row    = row;
        _col    = col;
        _value  = value;
    }

private:
    int         _row;
    int         _col;
    double      _value;
};

} // namespace cadet

#endif /* SPARSEMATRIXELEMENT_HPP_ */
