/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2020.

This file is part of BDSIM.

BDSIM is free software: you can redistribute it and/or modify 
it under the terms of the GNU General Public License as published 
by the Free Software Foundation version 3 of the License.

BDSIM is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BDSIM.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef TRKQuadrupole_h
#define TRKQuadrupole_h

#include "TRKElement.hh"

/**
 * @brief Quadrupole. Strength in Tesla/metre.
 */

class TRKQuadrupole: public TRKElement
{
public:
  TRKQuadrupole(double        strength,
		std::string   name,
		double        length,
		TRKAperture  *aperture,
		TRKPlacement *placement);
  virtual ~TRKQuadrupole();

  virtual void Track(TRKParticle* bunch, double step, TRKStrategy* strategy);

protected:
  /// output stream
  virtual void Print(std::ostream& out) const;

private:
  TRKQuadrupole(); ///< not implemented
};

#endif
