/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2022.

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

#ifndef TRKBACKTRACKER_HH
#define TRKBACKTRACKER_HH

#include "BDSHitApertureImpact.hh"

#include "TRKParticle.hh"

class TRKElement;
class TRKStrategy;

class TRKBacktracker {
public:
  TRKBacktracker(TRKStrategy *strategy,
		 bool backtracking,
                 double lossPrecision);
  TRKBacktracker() = delete;
  BDSHitApertureImpact GetApertureImpact(TRKParticle &lostParticle,
                                         TRKElement *element,
                                         double step) const;

private:
  int NBisectionSteps(double interval) const;
  void BacktrackToLossPoint(TRKParticle &lostParticle,
			    TRKElement *element,
			    double step) const;
  void ParticleToBDSHitAperImpact(TRKParticle const &p) const;

  TRKStrategy * strategy;
  bool backtracking;
  double lossPrecision;


};

#endif /* TRKBACKTRACKER_HH */
