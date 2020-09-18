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
#ifndef TRKRBend_h
#define TRKRBend_h

#include "TRKElement.hh"

/**
 * @brief bend
 * 
 * dipole tracking but reference s-position is changed
 */

class TRKRBend: public TRKElement
{
public:
  TRKRBend(double        angle,
	   std::string   name,
	   double        length,
	   TRKAperture  *aperture,
	   TRKPlacement *placement,
	   double        k1);
  virtual ~TRKRBend();
  
  inline double GetAngle() const {return angle;}
  inline double GetK1() const {return k1;}  
  
  virtual void Track(TRKParticle& particle, double step, TRKStrategy* strategy);
  
protected:
  /// output stream
  virtual void Print(std::ostream& out) const;

private:
  TRKRBend(); ///< not implemented

  /// bending angle in rad
  double angle;
  double k1;
};

#endif
