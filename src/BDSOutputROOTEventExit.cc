/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2018.

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
#include "BDSOutputROOTEventExit.hh"

#ifndef __ROOTBUILD__
#include "CLHEP/Units/SystemOfUnits.h"
#include "BDSVolumeExitHit.hh"
#endif

ClassImp(BDSOutputROOTEventExit)

BDSOutputROOTEventExit::BDSOutputROOTEventExit()
{
  Flush();
}

BDSOutputROOTEventExit::~BDSOutputROOTEventExit()
{;}

#ifndef __ROOTBUILD__
void BDSOutputROOTEventExit::Fill(const BDSVolumeExitHit* hit)
{
  n++;
  totalEnergy.push_back( (float &&) (hit->totalEnergy / CLHEP::GeV));
  postStepKineticEnergy.push_back( (float &&) (hit->postStepKineticEnergy / CLHEP::GeV));
  X.push_back( (float &&) (hit->X / CLHEP::m));
  Y.push_back( (float &&) (hit->Y / CLHEP::m));
  Z.push_back( (float &&) (hit->Z / CLHEP::m));
  T.push_back( (float &&) (hit->T / CLHEP::ns));
  partID.push_back(hit->partID);
  trackID.push_back(hit->trackID);
  parentID.push_back(hit->parentID);
  weight.push_back(hit->weight);
  turn.push_back(hit->turnsTaken);
}
#endif

void BDSOutputROOTEventExit::Flush()
{
  n = 0;
  totalEnergy.clear();
  postStepKineticEnergy.clear();
  X.clear();
  Y.clear();
  Z.clear();
  T.clear();
  partID.clear();  
  trackID.clear();
  parentID.clear();
  weight.clear();
  turn.clear();
}
