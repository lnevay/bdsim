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
#ifndef BDSVOLUMEEXITSD_H
#define BDSVOLUMEEXITSD_H

#include "BDSVolumeExitHit.hh"

#include "globals.hh"
#include "G4StepStatus.hh"
#include "G4VSensitiveDetector.hh"


class BDSAuxiliaryNavigator;

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

/**
 * @brief Generates BDSVolumeExitHits if a particle is leaving a volume.
 *
 * If the step results in exiting a volume, record the global coordinates.
 * Can do general volume exit or world volume exit (default). See G4StepStatus.
 */

class BDSVolumeExitSD: public G4VSensitiveDetector
{
public:
  BDSVolumeExitSD(G4String name,
		  G4bool   worldExit = true);

  virtual ~BDSVolumeExitSD(){;}

  virtual void Initialize(G4HCofThisEvent* HCE);
  virtual G4bool ProcessHits(G4Step* step,
			     G4TouchableHistory* th);
  
private:
  /// assignment and copy constructor not implemented nor used
  BDSVolumeExitSD& operator=(const BDSVolumeExitSD&);
  BDSVolumeExitSD(BDSVolumeExitSD&);
  BDSVolumeExitSD() = delete;

  G4String     colName;        ///< Name prepared for collection.
  G4int        HCIDve;         ///< Hits collection ID for volume exit.
  G4StepStatus statusToMatch;  ///< World or volume exit status cache.
  BDSVolumeExitHitsCollection* collection; ///< Hits collection.
};

#endif

