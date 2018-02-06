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
#include "BDSDebug.hh"
#include "BDSIonDefinition.hh"
#include "BDSParticleDefinition.hh"
#include "BDSPhysicalConstants.hh"
#include "BDSUtilities.hh"

#include "G4ParticleDefinition.hh"

#include <cmath>
#include <limits>
#include <ostream>
#include <stdexcept>


BDSParticleDefinition::BDSParticleDefinition(G4ParticleDefinition* particleIn,
					     G4double              totalEnergyIn,
					     G4double              ffact,
					     BDSIonDefinition*     ionDefinitionIn):
  particle(particleIn),
  ionDefinition(ionDefinitionIn),
  name(particleIn->GetParticleName()),
  mass(particleIn->GetPDGMass()),
  totalEnergy(totalEnergyIn),
  brho(std::numeric_limits<double>::max())// if zero charge infinite magnetic rigidity
{
  charge = particle->GetPDGCharge();
  if (ionDefinition) // may be nullptr
    {
      if (ionDefinition->OverrideCharge()) // if override for ions
	{charge = ionDefinition->Charge();}
    }
  
  kineticEnergy = totalEnergy - mass;

  CalculateMomentum();
  CalculateRigidity(ffact);
}

BDSParticleDefinition::BDSParticleDefinition(G4String          nameIn,
					     G4double          massIn,
					     G4double          chargeIn,
					     G4double          totalEnergyIn,
					     G4double          ffact,
					     BDSIonDefinition* ionDefinitionIn):
  particle(nullptr),
  ionDefinition(ionDefinitionIn),
  name(nameIn),
  mass(massIn),
  charge(chargeIn),
  totalEnergy(totalEnergyIn),
  brho(std::numeric_limits<double>::max())// if zero charge infinite magnetic rigidity
{
  kineticEnergy = totalEnergy - mass;

  CalculateMomentum();
  CalculateRigidity(ffact);
}

BDSParticleDefinition::~BDSParticleDefinition()
{
  delete ionDefinition;
}

std::ostream& operator<<(std::ostream& out, const BDSParticleDefinition& def)
{
  out << "Particle:       \""<< def.name << "\"" << G4endl;
  out << "Mass:            " << def.mass/CLHEP::GeV              << " GeV" << G4endl;
  out << "Charge:          " << def.charge                       << " e"   << G4endl;
  out << "Total Energy:    " << def.totalEnergy/CLHEP::GeV       << " GeV" << G4endl;
  out << "Kinetic Energy:  " << def.kineticEnergy/CLHEP::GeV     << " GeV" << G4endl;
  out << "Momentum:        " << def.momentum/CLHEP::GeV          << " GeV" << G4endl;
  out << "Rigidity (Brho): " << def.brho/(CLHEP::tesla*CLHEP::m) << " T*m" << G4endl;
  return out;
}

void BDSParticleDefinition::CalculateMomentum()
{
  try
    {momentum = std::sqrt(std::pow(totalEnergy,2) - std::pow(mass,2));}
  catch (std::domain_error) // sqrt(-ve)
    {
      G4cerr << __METHOD_NAME__ << "Total energy insufficient to include mass or particle" << G4endl;
      exit(1);
    }
}

void BDSParticleDefinition::CalculateRigidity(const G4double& ffact)
{
  // magnetic rigidity (brho)
  // formula: B(Tesla)*rho(m) = p(GeV)/(0.299792458 * charge(e))
  // charge (in e units); rigidity (in T*m)
  if (BDS::IsFinite(charge))
    {
      brho = ffact * momentum / CLHEP::GeV / BDS::cOverGeV / charge;
      brho *= CLHEP::tesla*CLHEP::m; // rigidity (in Geant4 units)
    }
}
