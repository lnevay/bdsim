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
#include "G4Version.hh"
#if G4VERSION_NUMBER > 1039
#include "BDSPhysicsRadioactivation.hh"

#include "G4GenericIon.hh"
#include "G4EmParameters.hh"
#include "G4LossTableManager.hh"
#include "G4NuclideTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicsListHelper.hh"
#include "G4Radioactivation.hh"
#include "G4Types.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4VAtomDeexcitation.hh"

#include "G4ParticleTypes.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4IonConstructor.hh"

#include "BDSGlobalConstants.hh"

BDSPhysicsRadioactivation::BDSPhysicsRadioactivation(G4bool atomicRearrangementIn):
  G4VPhysicsConstructor("BDSPhysicsRadioactivation"),
  atomicRearrangement(atomicRearrangementIn)
{
    // mandatory for G4NuclideTable
    //
    G4NuclideTable::GetInstance()->SetThresholdOfHalfLife(0.1*CLHEP::picosecond);
    G4NuclideTable::GetInstance()->SetLevelTolerance(1.0*CLHEP::eV);

    //read new PhotonEvaporation data set
    //
    G4DeexPrecoParameters* deex =
            G4NuclearLevelData::GetInstance()->GetParameters();
    deex->SetCorrelatedGamma(false);
    deex->SetStoreAllLevels(true);
    deex->SetIsomerProduction(true);
    deex->SetMaxLifeTime(G4NuclideTable::GetInstance()->GetThresholdOfHalfLife()
                         /std::log(2.));}

BDSPhysicsRadioactivation::~BDSPhysicsRadioactivation()
{;}

void BDSPhysicsRadioactivation::ConstructParticle()
{
    // pseudo-particles
    G4Geantino::GeantinoDefinition();

    // gamma
    G4Gamma::GammaDefinition();

    // leptons
    G4Electron::ElectronDefinition();
    G4Positron::PositronDefinition();

    G4NeutrinoE::NeutrinoEDefinition();
    G4AntiNeutrinoE::AntiNeutrinoEDefinition();

    // baryons
    G4Proton::ProtonDefinition();
    G4Neutron::NeutronDefinition();

    // ions
    G4IonConstructor iConstructor;
    iConstructor.ConstructParticle();
}

void BDSPhysicsRadioactivation::ConstructProcess()
{
  if (Activated())
    {return;}

  G4Radioactivation* ra = new G4Radioactivation();

  // atomic rearrangement
  ra->SetARM(atomicRearrangement);
  ra->SetAnalogueMonteCarlo(BDSGlobalConstants::Instance()->AnalogueMC()); // if FALSE: means that BRBias is activated per default, NSplit = 1 and Time Biasing between [0,1] sec.

  if (!BDSGlobalConstants::Instance()->AnalogueMC())
  {
      ra->SetSplitNuclei(BDSGlobalConstants::Instance()->NSplit());
      ra->SetBRBias(BDSGlobalConstants::Instance()->BRBias());

      if (!BDSGlobalConstants::Instance()->DecayBiasFilename().empty()){
          ra->SetDecayBias(BDSGlobalConstants::Instance()->DecayBiasFilename());
      }
  }

  // initialise atomic deexcitation
  G4LossTableManager* man = G4LossTableManager::Instance();
  G4VAtomDeexcitation* ad = man->AtomDeexcitation();
  if(!ad)
    {
      G4EmParameters::Instance()->SetAugerCascade(true);
      G4EmParameters::Instance()->SetDeexcitationIgnoreCut(true);
      ad = new G4UAtomicDeexcitation();
      ad->InitialiseAtomicDeexcitation();
      man->SetAtomDeexcitation(ad);
    }
  
  G4PhysicsListHelper::GetPhysicsListHelper()->RegisterProcess(ra, G4GenericIon::GenericIon());

  SetActivated();
}

#endif
