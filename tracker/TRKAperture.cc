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
#include "TRKAperture.hh"
#include <iterator>
#include <vector>
#include <iostream>
#include <algorithm> //for swap
#include <utility>

// Loss output
#include "BDSDebug.hh"
#include "BDSOutputBase.hh"
#include "BDSGlobalConstants.hh"

extern BDSOutputBase* trkOutput;

std::ostream& operator<<(std::ostream& out, const TRKAperture& aperture)
{
  aperture.PrintDetails(out);
  return out;
}

TRKBunch *TRKAperture::CheckAperture(TRKBunch *bunch) {
  auto lostIt = std::partition(
      bunch->begin(), bunch->end(),
      [this](TRKParticle p) { return !this->OutsideAperture(p); });

  if (lostIt == bunch->end()) {
    return new TRKBunch();;
  }
  auto nlost = std::distance(lostIt, bunch->end());
  std::vector<TRKParticle> lostParticles;
  lostParticles.reserve(nlost);
  std::copy(lostIt, bunch->end(), lostParticles.begin());
  
  bunch->Erase(lostIt, bunch->end());
  return new TRKBunch(std::move(lostParticles));
}
