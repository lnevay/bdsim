/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2024.

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
#include "BDSParticleCoordsSpherical.hh"
#include "BDSHitSamplerSphere.hh"

#include "globals.hh"
#include "G4Allocator.hh"

G4Allocator<BDSHitSamplerSphere> BDSAllocatorSamplerSphere;

BDSHitSamplerSphere::BDSHitSamplerSphere(G4int samplerIDIn,
					 const BDSParticleCoordsSpherical& coordsIn,
					 G4double totalEnergyIn,
					 G4double weightIn,
					 G4double SIn,
					 G4double momentumIn,
					 G4double massIn,
					 G4double chargeIn,
					 G4double rigidityIn,
					 G4int    pdgIDIn,
					 G4int    parentIDIn,
					 G4int    trackIDIn,
					 G4int    turnsTakenIn,
					 G4int    beamlineIndexIn,
					 G4int    nElectronsIn):
  samplerID(samplerIDIn),
  coords(coordsIn),
  totalEnergy(totalEnergyIn),
  weight(weightIn),
  S(SIn),
  momentum(momentumIn),
  mass(massIn),
  charge(chargeIn),
  rigidity(rigidityIn),
  pdgID(pdgIDIn),
  parentID(parentIDIn),
  trackID(trackIDIn),
  turnsTaken(turnsTakenIn),
  beamlineIndex(beamlineIndexIn),
  nElectrons(nElectronsIn)
{;}

BDSHitSamplerSphere::~BDSHitSamplerSphere()
{;}
