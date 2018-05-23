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
#include "BDSGlobalConstants.hh"
#include "BDSIonDefinition.hh"
#include "BDSModularPhysicsList.hh"
#include "BDSParticleDefinition.hh"
#include "BDSPhysicalConstants.hh"
#include "BDSPhysicsCherenkov.hh"
#include "BDSPhysicsCutsAndLimits.hh"
#include "BDSPhysicsEMDissociation.hh"
#include "BDSPhysicsLaserWire.hh"
#include "BDSPhysicsMuon.hh"
#include "BDSPhysicsSynchRad.hh"
#include "BDSUtilities.hh"

#include "parser/fastlist.h"
#include "parser/physicsbiasing.h"

// general geant4
#include "globals.hh"
#include "G4GenericBiasingPhysics.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4Version.hh"

// physics processes / builders (assumed Geant4.10.0 and upwards)
#include "G4ChargeExchangePhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmLivermorePolarizedPhysics.hh"
#include "G4EmLowEPPhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysicsSS.hh"
#include "G4EmStandardPhysicsWVI.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4HadronDElasticPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronElasticPhysicsLEND.hh"
#include "G4HadronElasticPhysicsXS.hh"
#include "G4HadronHElasticPhysics.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4HadronPhysicsQGSP_BERT.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronPhysicsShielding.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4IonElasticPhysics.hh"
#include "G4IonINCLXXPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4IonQMDPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4OpticalPhysics.hh"
#include "G4OpticalProcessIndex.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4SynchrotronRadiation.hh"

// version specific physics lists
#if G4VERSION_NUMBER > 1009
#include "G4HadronElasticPhysicsPHP.hh"
#endif

#if G4VERSION_NUMBER > 1019
#include "G4EmStandardPhysicsGS.hh"
#endif

#if G4VERSION_NUMBER > 1020
#include "G4SpinDecayPhysics.hh"
#endif

#if G4VERSION_NUMBER > 1022
#include "G4IonPhysicsPHP.hh"
#endif

#if G4VERSION_NUMBER > 1029
#include "G4MuonicAtomDecayPhysics.hh"
#endif

#if G4VERSION_NUMBER > 1039
#include "G4HadronPhysicsShieldingLEND.hh"
#endif

// particles
#include "G4AntiNeutrinoE.hh"
#include "G4AntiNeutron.hh"
#include "G4AntiProton.hh"
#include "G4BaryonConstructor.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4IonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4NeutrinoE.hh"
#include "G4Neutron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4ShortLivedConstructor.hh"

#include <iomanip>
#include <iterator>
#include <limits>
#include <map>
#include <ostream>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

BDSModularPhysicsList::BDSModularPhysicsList(G4String physicsList):
  opticalPhysics(nullptr),
  emWillBeUsed(false),
  usingIons(false)
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif 
  globals = BDSGlobalConstants::Instance();
  
  SetVerboseLevel(1);
  
  physicsConstructors.insert(std::make_pair("charge_exchange",        &BDSModularPhysicsList::ChargeExchange));
  physicsConstructors.insert(std::make_pair("cherenkov",              &BDSModularPhysicsList::Cherenkov));
  physicsConstructors.insert(std::make_pair("cuts_and_limits",        &BDSModularPhysicsList::CutsAndLimits));
  physicsConstructors.insert(std::make_pair("decay",                  &BDSModularPhysicsList::Decay));
  physicsConstructors.insert(std::make_pair("decay_radioactive",      &BDSModularPhysicsList::DecayRadioactive));
  physicsConstructors.insert(std::make_pair("em",                     &BDSModularPhysicsList::Em));
  physicsConstructors.insert(std::make_pair("em_extra",               &BDSModularPhysicsList::EmExtra));
  physicsConstructors.insert(std::make_pair("em_livermore",           &BDSModularPhysicsList::EmLivermore));
  physicsConstructors.insert(std::make_pair("em_livermore_polarised", &BDSModularPhysicsList::EmLivermorePolarised));
  physicsConstructors.insert(std::make_pair("em_low_ep",              &BDSModularPhysicsList::EmLowEP));
  physicsConstructors.insert(std::make_pair("em_penelope",            &BDSModularPhysicsList::EmPenelope));
  physicsConstructors.insert(std::make_pair("em_ss",                  &BDSModularPhysicsList::EmSS));
  physicsConstructors.insert(std::make_pair("em_wvi",                 &BDSModularPhysicsList::EmWVI));
  physicsConstructors.insert(std::make_pair("em_1",                   &BDSModularPhysicsList::Em1));
  physicsConstructors.insert(std::make_pair("em_2",                   &BDSModularPhysicsList::Em2));
  physicsConstructors.insert(std::make_pair("em_3",                   &BDSModularPhysicsList::Em3));
  physicsConstructors.insert(std::make_pair("em_4",                   &BDSModularPhysicsList::Em4));
  physicsConstructors.insert(std::make_pair("ftfp_bert",              &BDSModularPhysicsList::FTFPBERT));
  physicsConstructors.insert(std::make_pair("ftfp_bert_hp",           &BDSModularPhysicsList::FTFPBERTHP));
  physicsConstructors.insert(std::make_pair("hadronic_elastic",       &BDSModularPhysicsList::HadronicElastic));
  physicsConstructors.insert(std::make_pair("hadronic_elastic_d",     &BDSModularPhysicsList::HadronicElasticD));
  physicsConstructors.insert(std::make_pair("hadronic_elastic_h",     &BDSModularPhysicsList::HadronicElasticH));
  physicsConstructors.insert(std::make_pair("hadronic_elastic_hp",    &BDSModularPhysicsList::HadronicElasticHP));
  physicsConstructors.insert(std::make_pair("hadronic_elastic_lend",  &BDSModularPhysicsList::HadronicElasticLEND));
  physicsConstructors.insert(std::make_pair("hadronic_elastic_xs",    &BDSModularPhysicsList::HadronicElasticXS));
  physicsConstructors.insert(std::make_pair("ion",                    &BDSModularPhysicsList::Ion));
  physicsConstructors.insert(std::make_pair("ion_binary",             &BDSModularPhysicsList::IonBinary));
  physicsConstructors.insert(std::make_pair("ion_elastic",            &BDSModularPhysicsList::IonElastic));
  physicsConstructors.insert(std::make_pair("ion_elastic_qmd",        &BDSModularPhysicsList::IonElasticQMD));
  physicsConstructors.insert(std::make_pair("ion_em_dissociation",    &BDSModularPhysicsList::IonEMDissociation));
  physicsConstructors.insert(std::make_pair("ion_inclxx",             &BDSModularPhysicsList::IonINCLXX));
  physicsConstructors.insert(std::make_pair("lw",                     &BDSModularPhysicsList::LaserWire));
  physicsConstructors.insert(std::make_pair("muon",                   &BDSModularPhysicsList::Muon));
  physicsConstructors.insert(std::make_pair("neutron_tracking_cut",   &BDSModularPhysicsList::NeutronTrackingCut));
  physicsConstructors.insert(std::make_pair("optical",                &BDSModularPhysicsList::Optical));
  physicsConstructors.insert(std::make_pair("qgsp_bert",              &BDSModularPhysicsList::QGSPBERT));
  physicsConstructors.insert(std::make_pair("qgsp_bert_hp",           &BDSModularPhysicsList::QGSPBERTHP));
  physicsConstructors.insert(std::make_pair("qgsp_bic",               &BDSModularPhysicsList::QGSPBIC));
  physicsConstructors.insert(std::make_pair("qgsp_bic_hp",            &BDSModularPhysicsList::QGSPBICHP));
  physicsConstructors.insert(std::make_pair("shielding",              &BDSModularPhysicsList::Shielding));
  physicsConstructors.insert(std::make_pair("stopping",               &BDSModularPhysicsList::Stopping));
  physicsConstructors.insert(std::make_pair("synch_rad",              &BDSModularPhysicsList::SynchRad));
#if G4VERSION_NUMBER > 1019
  physicsConstructors.insert(std::make_pair("em_gs",                  &BDSModularPhysicsList::EmGS));
#endif
#if G4VERSION_NUMBER > 1020
  physicsConstructors.insert(std::make_pair("decay_spin",             &BDSModularPhysicsList::DecaySpin));
#endif
#if G4VERSION_NUMBER > 1022
  physicsConstructors.insert(std::make_pair("ion_php",                &BDSModularPhysicsList::IonPHP));
#endif
#if G4VERSION_NUMBER > 1029
  physicsConstructors.insert(std::make_pair("decay_muonic_atom",      &BDSModularPhysicsList::DecayMuonicAtom));
#endif
#if G4VERSION_NUMBER > 1039
  physicsConstructors.insert(std::make_pair("shielding_lend",         &BDSModularPhysicsList::ShieldingLEND));
#endif

  // old names and aliases
  aliasToOriginal["cerenkov"]      = "cherenkov";
  aliasToOriginal["cutsandlimits"] = "cuts_and_limits";
  aliasToOriginal["em_low"]        = "em_penelope";
  aliasToOriginal["hadronic"]      = "ftfp_bert";
  aliasToOriginal["hadronic_hp"]   = "ftfp_bert_hp";
  aliasToOriginal["ionbinary"]     = "ion_binary";
  aliasToOriginal["ioninclxx"]     = "ion_inclxx";
  aliasToOriginal["ionphp"]        = "ion_php";
  aliasToOriginal["spindecay"]     = "decay_spin";
  aliasToOriginal["synchrad"]      = "synch_rad";
  
  // prepare vector of valid names for searching when parsing physics list string
  for (const auto& constructor : physicsConstructors)
    {
      physicsLists.push_back(constructor.first);
      physicsActivated[constructor.first] = false;
    }

  // setup a list of incompatible physics lists for each one - mostly based on experience
  // initialise all to empty vectors and specify only ones that have some incompatible physics lists
  for (const auto kv : physicsConstructors)
    {incompatible.insert(std::make_pair(kv.first, std::vector<G4String>()));}
  incompatible["em"]     = {"em_ss", "em_wvi", "em_1",   "em_2", "em_3", "em_4"};
  incompatible["em_ss"]  = {"em",    "em_wvi", "em_1",   "em_2", "em_3", "em_4"};
  incompatible["em_wvi"] = {"em",    "em_ss",  "em_1",   "em_2", "em_3", "em_4"};
  incompatible["em_1"]   = {"em",    "em_ss",  "em_wvi", "em_2", "em_3", "em_4"};
  incompatible["em_2"]   = {"em",    "em_ss",  "em_wvi", "em_1", "em_3", "em_4"};
  incompatible["em_3"]   = {"em",    "em_ss",  "em_wvi", "em_1", "em_2", "em_4"};
  incompatible["em_4"]   = {"em",    "em_ss",  "em_wvi", "em_1", "em_2", "em_3"};
  incompatible["em_livermore"] = {"em_livermore_polarised"};
  incompatible["ftfp_bert"]    = {"ftfp_bert_hp", "qgsp_bert", "qgsp_bert_hp", "qgsp_bic", "qgsp_bic_hp"};
  incompatible["ftfp_bert_hp"] = {"ftfp_bert",    "qgsp_bert", "qgsp_bert_hp", "qgsp_bic", "qgsp_bic_hp"};
  incompatible["hadronic_elastic"]      = {"hadronic_elastic_d", "hadronic_elastic_h", "hadronic_elastic_hp", "hadronic_elastic_lend", "hadronic_elastic_xs"};
  incompatible["hadronic_elastic_d"]    = {"hadronic_elastic",   "hadronic_elastic_h", "hadronic_elastic_hp", "hadronic_elastic_lend", "hadronic_elastic_xs"};
  incompatible["hadronic_elastic_h"]    = {"hadronic_elastic",   "hadronic_elastic_d", "hadronic_elastic_hp", "hadronic_elastic_lend", "hadronic_elastic_xs"};
  incompatible["hadronic_elastic_hp"]   = {"hadronic_elastic",   "hadronic_elastic_d", "hadronic_elastic_h",  "hadronic_elastic_lend", "hadronic_elastic_xs"};
  incompatible["hadronic_elastic_lend"] = {"hadronic_elastic",   "hadronic_elastic_d", "hadronic_elastic_h",  "hadronic_elastic_hp",   "hadronic_elastic_xs"};
  incompatible["hadronic_elastic_xs"]   = {"hadronic_elastic",   "hadronic_elastic_d", "hadronic_elastic_h",  "hadronic_elastic_hp",   "hadronic_elastic_lend"};
  incompatible["ion_elastic"] = {"ion_elastic_qmd"};
  incompatible["qgsp_bert"]    = {"ftfp_bert", "ftfp_bert_hp", "qgsp_bert_hp", "qgsp_bic",     "qgsp_bic_hp"};
  incompatible["qgsp_bert_hp"] = {"ftfp_bert", "ftfp_bert_hp", "qgsp_bert",    "qgsp_bic",     "qgsp_bic_hp"};
  incompatible["qgsp_bic"]     = {"ftfp_bert", "ftfp_bert_hp", "qgsp_bert",    "qgsp_bert_hp", "qgsp_bic_hp"};
  incompatible["qgsp_bic_hp"]  = {"ftfp_bert", "ftfp_bert_hp", "qgsp_bert",    "qgsp_bert_hp", "qgsp_bic"};

#if G4VERSION_NUMBER > 1019
  for (const auto& name : {"em", "em_ss", "em_wvi", "em_1", "em_2", "em_3", "em_4"})
    {incompatible[name].push_back("em_gs");}
  incompatible["em_gs"] = {"em",    "em_ss",  "em_wvi", "em_1", "em_2", "em_3", "em_4"};
#endif
#if G4VERSION_NUMBER > 1020
  incompatible["decay"].push_back("decay_spin"); // append for safety in future
  incompatible["decay_spin"] = {"decay"};
#endif
#if G4VERSION_NUMBER > 1039
  incompatible["shielding"].push_back("shielding_lend");
  incompatible["shielding_lend"] = {"shielding"};
#endif
  
  ParsePhysicsList(physicsList);
  ConfigurePhysics();

  // register the physics constructors with base class mechanics.
  for(auto physics : constructors)
    {RegisterPhysics(physics);}
  
#ifdef BDSDEBUG
  Print();
#endif
}

BDSModularPhysicsList::~BDSModularPhysicsList()
{;}

void BDSModularPhysicsList::ConstructParticle()
{
  ConstructMinimumParticleSet();
  G4VModularPhysicsList::ConstructParticle();
}

void BDSModularPhysicsList::ConstructProcess()
{
  G4VModularPhysicsList::ConstructProcess();
  DumpCutValuesTable(100);
}

void BDSModularPhysicsList::Print()
{
  for (const auto& physics : physicsActivated)
    {
      G4String result = (physics.second ? "active" : "inactive");
      G4cout << std::setw(27) << ("\"" + physics.first + "\": ") << result << G4endl;
    }
}

void BDSModularPhysicsList::PrintDefinedParticles() const
{
  G4cout << __METHOD_NAME__ << "Defined particles: " << G4endl;
  auto it = G4ParticleTable::GetParticleTable()->GetIterator();
  it->reset();
  while ((*it)())
    {G4cout <<  it->value()->GetParticleName() << " ";}
  G4cout << G4endl;
}

void BDSModularPhysicsList::PrintPrimaryParticleProcesses() const
{
  auto particleName = globals->ParticleName();
  G4cout << "Register physics processes by name for the primary particle \"" << particleName << "\":" << G4endl;

  auto particle = G4ParticleTable::GetParticleTable()->FindParticle(particleName);
  if (!particle)
    {// could be ion that isn't defined
#ifdef BDSDEBUG
      G4cout << __METHOD_NAME__ << "primary particle not defined yet - could be ion" << G4endl;
#endif
      return;
    }
  auto pl = particle->GetProcessManager()->GetProcessList();
  for (G4int i = 0; i < pl->length(); i++)
    { G4cout << "\"" << (*pl)[i]->GetProcessName() << "\"" << G4endl; }
}

void BDSModularPhysicsList::ParsePhysicsList(G4String physListName)
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << "Physics list string: \"" << physListName << "\"" << G4endl;
#endif
  // string stream to vector will take a single string that contains words
  // delimited by whitespace and split them on the whitespace
  std::stringstream ss(physListName);
  std::istream_iterator<std::string> begin(ss);
  std::istream_iterator<std::string> end;
  std::vector<std::string> physicsListNamesS(begin, end);

  // convert to G4String for lower case convenience
  std::vector<G4String> physicsListNames;
  for (auto physicsListName : physicsListNamesS)
    {
      G4String name = G4String(physicsListName); // convert string to G4String.
      name.toLower(); // change to lower case - physics lists are case insensitive

      // search aliases
      auto result = aliasToOriginal.find(name);
      if (result != aliasToOriginal.end())
	{
	  G4cout << __METHOD_NAME__ << "alias \"" << name << "\" forwarding to \""
		 << result->second << "\"" << G4endl;
	  name = result->second; // overwrite name with the correct one
	}
      physicsListNames.push_back(name);
    }

  // seach for em physics (could be any order) - needed for different construction of muon phyiscs
  if (std::find(physicsListNames.begin(), physicsListNames.end(), "em") != physicsListNames.end())
    {emWillBeUsed = true;}

  for (const auto name : physicsListNames)
    {
      auto result = physicsConstructors.find(name);
      if (result != physicsConstructors.end())
	{
	  G4cout << __METHOD_NAME__ << "Constructing \"" << result->first << "\" physics list" << G4endl;
	  CheckIncompatiblePhysics(name);
	  auto mem = result->second;
	  (this->*mem)(); // call the function pointer in this instance of the class
	}
      else
	{
	  G4cout << "\"" << name << "\" is not a valid physics list. Available ones are: " << G4endl;
	  for (auto listName : physicsLists)
	    {G4cout << "\"" << listName << "\"" << G4endl;}
	  exit(1);
	}
    }

  //Always load cuts and limits.
  CutsAndLimits();
}

void BDSModularPhysicsList::ConstructMinimumParticleSet()
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif
  // e-, e+, v_e, v_e(bar)
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();

  // p, pbar, neutron, anti-neutron
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();

  // photon
  G4Gamma::Gamma();
}

void BDSModularPhysicsList::ConstructAllLeptons()
{
  G4LeptonConstructor leptons;
  leptons.ConstructParticle();
}

void BDSModularPhysicsList::ConstructAllShortLived()
{
  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();
}

void BDSModularPhysicsList::ConstructAllMesons()
{
  G4MesonConstructor mConstructor;
  mConstructor.ConstructParticle();
}

void BDSModularPhysicsList::ConstructAllBaryons()
{
  G4BaryonConstructor bConstructor;
  bConstructor.ConstructParticle();
}

void BDSModularPhysicsList::ConstructAllIons()
{
  usingIons = true; // all physics lists that use ions call this function so put this here
  G4GenericIon::GenericIonDefinition();
  G4IonConstructor iConstructor;
  iConstructor.ConstructParticle();
}

void BDSModularPhysicsList::ConfigurePhysics()
{
  if(opticalPhysics)
    {ConfigureOptical();}
}

void BDSModularPhysicsList::ConfigureOptical()
{
  // cherenkov turned on with optical even if it's not on as separate list
  opticalPhysics->Configure(G4OpticalProcessIndex::kCerenkov, true);
  opticalPhysics->Configure(G4OpticalProcessIndex::kScintillation, true);                                ///< Scintillation process index
  opticalPhysics->Configure(G4OpticalProcessIndex::kAbsorption,    globals->TurnOnOpticalAbsorption());  ///< Absorption process index
  opticalPhysics->Configure(G4OpticalProcessIndex::kRayleigh,      globals->TurnOnRayleighScattering()); ///< Rayleigh scattering process index
  opticalPhysics->Configure(G4OpticalProcessIndex::kMieHG,         globals->TurnOnMieScattering());      ///< Mie scattering process index
  opticalPhysics->Configure(G4OpticalProcessIndex::kBoundary,      globals->TurnOnOpticalSurface());     ///< Boundary process index
  opticalPhysics->Configure(G4OpticalProcessIndex::kWLS,           true);                                ///< Wave Length Shifting process index
  opticalPhysics->SetScintillationYieldFactor(globals->ScintYieldFactor());
  G4long maxPhotonsPerStep = globals->MaximumPhotonsPerStep();
  if (maxPhotonsPerStep >= 0)
    {opticalPhysics->SetMaxNumPhotonsPerStep(maxPhotonsPerStep);}
}

void BDSModularPhysicsList::CheckIncompatiblePhysics(const G4String& singlePhysicsIn) const
{
  // no need to check if key is present as we control both maps in the constructor
  const std::vector<G4String>& forbidden = incompatible.at(singlePhysicsIn);
  
  for (const auto key : forbidden)
    {// for each forbidden physics list, check if it's activated
      if (physicsActivated.at(key))
	{
	  G4cerr << __METHOD_NAME__ << "Incompatible physics list \"" << singlePhysicsIn
		 << "\" being used with already used \"" << key << "\"" << G4endl;
	  G4cout << "\"" << singlePhysicsIn << "\" cannot be used with the following:" << G4endl;
	  for (const auto& v : forbidden)
	    {G4cout << "\"" << v << "\"" << G4endl;}
	  exit(1);
	}
    }
}

void BDSModularPhysicsList::SetCuts()
{
  // set default value
  SetDefaultCutValue(globals->DefaultRangeCut());

  // overwrite when explicitly set in options
  if (globals->ProdCutPhotonsSet())
    {SetCutValue(globals->ProdCutPhotons(),  "gamma");}
  if (globals->ProdCutElectronsSet())
    {SetCutValue(globals->ProdCutElectrons(),"e-");}
  if (globals->ProdCutPositronsSet())
    {SetCutValue(globals->ProdCutPositrons(),"e+");}
  if (globals->ProdCutProtonsSet())
    {SetCutValue(globals->ProdCutProtons(),  "proton");}

  G4cout << __METHOD_NAME__ << "Default production range cut  " << GetDefaultCutValue()  << " mm" << G4endl;
  G4cout << __METHOD_NAME__ << "Photon production range cut   " << GetCutValue("gamma")  << " mm" << G4endl;
  G4cout << __METHOD_NAME__ << "Electron production range cut " << GetCutValue("e-")     << " mm" << G4endl;
  G4cout << __METHOD_NAME__ << "Positron production range cut " << GetCutValue("e+")     << " mm" << G4endl;
  G4cout << __METHOD_NAME__ << "Proton production range cut   " << GetCutValue("proton") << " mm" << G4endl;

#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << "List of all constructed particles by physics lists" << G4endl;
  for (auto particle : *G4ParticleTable::fDictionary)
    {G4cout << particle.second->GetParticleName() << ", ";}
  G4cout << G4endl;
#endif
  
  DumpCutValuesTable(); 
}

BDSParticleDefinition* BDSModularPhysicsList::ConstructBeamParticle(G4String particleName,
								    G4double totalEnergy,
								    G4double ffact) const
{
  particleName.toLower();

  BDSParticleDefinition* particleDefB = nullptr; // result can be constructed in two ways
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  
  if (particleName.contains("ion"))
    {
      usingIons = true;
      G4GenericIon::GenericIonDefinition(); // construct general ion particle
      auto ionDef = new BDSIonDefinition(particleName); // parse the ion definition

      G4IonTable* ionTable = particleTable->GetIonTable();
      G4double mass   = ionTable->GetIonMass(ionDef->Z(), ionDef->A());
      G4double charge = ionDef->Charge(); // correct even if overridden
      particleDefB = new BDSParticleDefinition(particleName, mass, charge,
					       totalEnergy, ffact, ionDef);
      // this takes ownership of the ion definition
    }
  else
    {
      ConstructBeamParticleG4(particleName);
      auto particleDef = particleTable->FindParticle(particleName);
      if (!particleDef)
	{
	  G4cerr << "Particle \"" << particleName << "\" not found: quitting!" << G4endl;
	  exit(1);
	}
      particleDefB = new BDSParticleDefinition(particleDef, totalEnergy, ffact);
    }

  globals->SetBeamParticleDefinition(particleDefB); // export for bunch distribution
  return particleDefB;
}

void BDSModularPhysicsList::ConstructBeamParticleG4(G4String name) const
{
  if (name == "proton")
    {G4Proton::ProtonDefinition();}
  else if (name == "antiproton")
    {G4AntiProton::AntiProtonDefinition();}
  else if (name == "e-")
    {G4Electron::ElectronDefinition();}
  else if (name == "e+")
    {G4Positron::PositronDefinition();}
  else if (name == "pi-")
    {G4PionMinus::PionMinusDefinition();}
  else if (name == "pi+")
    {G4PionPlus::PionPlusDefinition();}
  else if (name == "neutron")
    {G4Neutron::NeutronDefinition();}
  else if (name == "photon" || name == "gamma")
    {G4Gamma::Gamma();}
  else if (name == "mu-")
    {G4MuonMinus::MuonMinusDefinition();}
  else if (name == "mu+")
    {G4MuonPlus::MuonPlusDefinition();}
  else
    {G4cerr << "Unknown particle type \"" << name << "\"" << G4endl;}
}

void BDSModularPhysicsList::ChargeExchange()
{
  ConstructAllIons();
  if (!physicsActivated["charge_exchange"])
    {
      constructors.push_back(new G4ChargeExchangePhysics());
      physicsActivated["charge_exchange"] = true;
    }
}

void BDSModularPhysicsList::Cherenkov()
{
  if (!physicsActivated["cherenkov"])
    {
      constructors.push_back(new BDSPhysicsCherenkov(BDSGlobalConstants::Instance()->MaximumPhotonsPerStep(),
						     BDSGlobalConstants::Instance()->MaximumBetaChangePerStep()));
      physicsActivated["cherenkov"] = true;
      if (!physicsActivated["em"])
	{Em();} // requires em physics to work (found empirically)
    }
}

void BDSModularPhysicsList::CutsAndLimits()
{
  if(!physicsActivated["cuts_and_limits"])
    {
      constructors.push_back(new BDSPhysicsCutsAndLimits());
      physicsActivated["cuts_and_limits"] = true;
    }
}

void BDSModularPhysicsList::Decay()
{
  if(!physicsActivated["decay"])
    {
      constructors.push_back(new G4DecayPhysics());
      physicsActivated["decay"] = true;
    }
}

void BDSModularPhysicsList::DecayRadioactive()
{
  if(!physicsActivated["decay_radioactive"])
    {
      constructors.push_back(new G4RadioactiveDecayPhysics());
      physicsActivated["decay_radioactive"] = true;
    }
}

void BDSModularPhysicsList::Em()
{
  ConstructAllLeptons();
  if (!physicsActivated["em"])
    {
      constructors.push_back(new G4EmStandardPhysics());
      physicsActivated["em"] = true;
    }
}

void BDSModularPhysicsList::EmExtra()
{
  ConstructAllLeptons();

  // These are required by GammaNuclear and MuonNuclear which
  // are activated by default in G4EmExtraPhysics.
  ConstructAllShortLived();
  ConstructAllBaryons();
  ConstructAllIons();
  ConstructAllMesons();

  if (!physicsActivated["em_extra"])
    {
      auto constructor = new G4EmExtraPhysics();
#if G4VERSION_NUMBER > 1012
      constructor->Synch(true); // introduced geant version 10.1
#endif
#if G4VERSION_NUMBER > 1039
      G4bool useLENDGammaNuclear = BDSGlobalConstants::Instance()->UseLENDGammaNuclear();
      if (useLENDGammaNuclear)
	{
	  BDS::CheckLowEnergyNeutronDataExists("em_extra");
	  constructor->LENDGammaNuclear(true);
	}
#endif
      constructors.push_back(constructor);
      physicsActivated["em_extra"] = true;
    }
}

void BDSModularPhysicsList::EmLivermore()
{
  ConstructAllLeptons();
  if (!physicsActivated["em_livermore"])
    {
      constructors.push_back(new G4EmLivermorePhysics());
      physicsActivated["em_livermore"] = true;
    }
}

void BDSModularPhysicsList::EmLivermorePolarised()
{
  ConstructAllLeptons();
  if (!physicsActivated["em_livermore_polarised"])
    {
      constructors.push_back(new G4EmLivermorePolarizedPhysics());
      physicsActivated["em_livermore_polarised"] = true;
    }
}

void BDSModularPhysicsList::EmLowEP()
{
  ConstructAllLeptons();
  if (!physicsActivated["em_low_ep"])
    {
      constructors.push_back(new G4EmLowEPPhysics());
      physicsActivated["em_low_ep"] = true;
    }
}
							  
void BDSModularPhysicsList::EmPenelope()
{
  ConstructAllLeptons();
  if (!physicsActivated["em_penelope"])
    {
      constructors.push_back(new G4EmPenelopePhysics());
      physicsActivated["em_penelope"] = true;
    }
}

void BDSModularPhysicsList::EmSS()
{
  ConstructAllLeptons();
  if (!physicsActivated["em_ss"])
    {
      constructors.push_back(new G4EmStandardPhysicsSS());
      physicsActivated["em_ss"] = true;
    }
}

void BDSModularPhysicsList::EmWVI()
{
  ConstructAllLeptons();
  if (!physicsActivated["em_wvi"])
    {
      constructors.push_back(new G4EmStandardPhysicsWVI());
      physicsActivated["em_wvi"] = true;
    }
}

void BDSModularPhysicsList::Em1()
{
  ConstructAllLeptons();
  if (!physicsActivated["em_1"])
    {
      constructors.push_back(new G4EmStandardPhysics_option1());
      physicsActivated["em_1"] = true;
    }
}

void BDSModularPhysicsList::Em2()
{
  ConstructAllLeptons();
  if (!physicsActivated["em_2"])
    {
      constructors.push_back(new G4EmStandardPhysics_option2());
      physicsActivated["em_2"] = true;
    }
}

void BDSModularPhysicsList::Em3()
{
  ConstructAllLeptons();
  if (!physicsActivated["em_3"])
    {
      constructors.push_back(new G4EmStandardPhysics_option3());
      physicsActivated["em_3"] = true;
    }
}

void BDSModularPhysicsList::Em4()
{
  ConstructAllLeptons();
  if (!physicsActivated["em_4"])
    {
      constructors.push_back(new G4EmStandardPhysics_option4());
      physicsActivated["em_4"] = true;
    }
}

void BDSModularPhysicsList::FTFPBERT()
{
  ConstructAllLeptons();
  HadronicElastic(); // has to be here to prevent G4 segfault
  if(!physicsActivated["ftfp_bert"])
    {
      constructors.push_back(new G4HadronPhysicsFTFP_BERT());
      physicsActivated["ftfp_bert"] = true;
    }
}

void BDSModularPhysicsList::FTFPBERTHP()
{
  ConstructAllLeptons();
  HadronicElastic(); // has to be here to prevent G4 segfault
  if(!physicsActivated["ftfp_bert_hp"])
    {
      constructors.push_back(new G4HadronPhysicsFTFP_BERT_HP());
      physicsActivated["ftfp_bert_hp"] = true;
    }
}

void BDSModularPhysicsList::HadronicElastic()
{
  ConstructAllLeptons();
  if (!physicsActivated["hadronic_elastic"])
    {
      constructors.push_back(new G4HadronElasticPhysics());
      physicsActivated["hadronic_elastic"] = true;
    }
}

void BDSModularPhysicsList::HadronicElasticD()
{
  ConstructAllLeptons();
  if (!physicsActivated["hadronic_elastic_d"])
    {
      constructors.push_back(new G4HadronDElasticPhysics());
      physicsActivated["hadronic_elastic_d"] = true;
    }
}

void BDSModularPhysicsList::HadronicElasticH()
{
  ConstructAllLeptons();
  if (!physicsActivated["hadronic_elastic_h"])
    {
      constructors.push_back(new G4HadronHElasticPhysics());
      physicsActivated["hadronic_elastic_h"] = true;
    }
}

void BDSModularPhysicsList::HadronicElasticHP()
{
  ConstructAllLeptons();
  if (!physicsActivated["hadronic_elastic_hp"])
    {
      constructors.push_back(new G4HadronElasticPhysicsHP());
      physicsActivated["hadronic_elastic_hp"] = true;
    }
}

void BDSModularPhysicsList::HadronicElasticLEND()
{
  BDS::CheckLowEnergyNeutronDataExists("hadronic_elastic_lend");
  ConstructAllLeptons();
  if (!physicsActivated["hadronic_elastic_lend"])
    {
      constructors.push_back(new G4HadronElasticPhysicsLEND());
      physicsActivated["hadronic_elastic_lend"] = true;
    }
}

void BDSModularPhysicsList::HadronicElasticXS()
{
  ConstructAllLeptons();
  if (!physicsActivated["hadronic_elastic_xs"])
    {
      constructors.push_back(new G4HadronElasticPhysicsXS());
      physicsActivated["hadronic_elastic_xs"] = true;
    }
}

void BDSModularPhysicsList::Ion()
{
  ConstructAllBaryons();
  ConstructAllIons();
  ConstructAllLeptons();
  ConstructAllMesons();
  ConstructAllShortLived();

  if (!physicsActivated["ion"])
    {
      constructors.push_back(new G4IonPhysics());
      physicsActivated["ion"] = true;
    }
}

void BDSModularPhysicsList::IonBinary()
{
  BDS::CheckHighPrecisionDataExists("ion_binary");
  ConstructAllBaryons();
  ConstructAllIons();
  ConstructAllLeptons();
  ConstructAllMesons();
  ConstructAllShortLived();
  
  if (!physicsActivated["ion_binary"])
    {
      constructors.push_back(new G4IonBinaryCascadePhysics());
      physicsActivated["ion_binary"] = true;
    }
}


void BDSModularPhysicsList::IonElastic()
{
  ConstructAllBaryons();
  ConstructAllIons();
  ConstructAllLeptons();
  ConstructAllMesons();
  ConstructAllShortLived();
  
  if (!physicsActivated["ion_elastic"])
    {
      constructors.push_back(new G4IonElasticPhysics());
      physicsActivated["ion_elastic"] = true;
    }
}


void BDSModularPhysicsList::IonElasticQMD()
{
  ConstructAllBaryons();
  ConstructAllIons();
  ConstructAllLeptons();
  ConstructAllMesons();
  ConstructAllShortLived();
  
  if (!physicsActivated["ion_elastic_qmd"])
    {
      constructors.push_back(new G4IonQMDPhysics());
      physicsActivated["ion_elastic_qmd"] = true;
    }
}

void BDSModularPhysicsList::IonEMDissociation()
{
  ConstructAllBaryons();
  ConstructAllIons();
  ConstructAllLeptons();
  ConstructAllMesons();
  ConstructAllShortLived();
  if (!physicsActivated["ion_em_dissociation"])
    {
      constructors.push_back(new BDSPhysicsEMDissociation());
      physicsActivated["ion_em_dissociation"] = true;
    }
}

void BDSModularPhysicsList::IonINCLXX()
{
  BDS::CheckHighPrecisionDataExists("ion_inclxx");
  ConstructAllBaryons();
  ConstructAllIons();
  ConstructAllLeptons();
  ConstructAllMesons();
  ConstructAllShortLived();
  
  if (!physicsActivated["ion_inclxx"])
    {
      constructors.push_back(new G4IonINCLXXPhysics());
      physicsActivated["ion_inclxx"] = true;
    }
}

void BDSModularPhysicsList::LaserWire()
{
  if(!physicsActivated["lw"])
    {
      constructors.push_back(new BDSPhysicsLaserWire());
      physicsActivated["lw"] = true;
    }
}							  
							  
void BDSModularPhysicsList::Muon()
{
  if(!physicsActivated["muon"])
    {
      constructors.push_back(new BDSPhysicsMuon(emWillBeUsed));
      physicsActivated["muon"] = true;
    }
}

void BDSModularPhysicsList::NeutronTrackingCut()
{
  if(!physicsActivated["neutron_tracking_cut"])
    {
      auto ntc = new G4NeutronTrackingCut();
      G4double timeLimit = BDSGlobalConstants::Instance()->NeutronTimeLimit();
      G4double eKinLimit = BDSGlobalConstants::Instance()->NeutronKineticEnergyLimit();
      G4cout << __METHOD_NAME__ << "Neutron time limit: " << timeLimit / CLHEP::s << " s" << G4endl;
      G4cout << __METHOD_NAME__ << "Neutron kinetic energy limit: " << eKinLimit / CLHEP::MeV << G4endl;
      ntc->SetTimeLimit(timeLimit);
      ntc->SetKineticEnergyLimit(eKinLimit);
      constructors.push_back(ntc);
      physicsActivated["neutron_tracking_cut"] = true;
    }
}
							  
void BDSModularPhysicsList::Optical()
{
  if(!physicsActivated["optical"])
    {
      opticalPhysics = new G4OpticalPhysics();		  
      constructors.push_back(opticalPhysics);
      physicsActivated["optical"] = true;
    }
}

void BDSModularPhysicsList::QGSPBERT()
{
  ConstructAllLeptons();
  if(!physicsActivated["qgsp_bert"])
    {
      constructors.push_back(new G4HadronPhysicsQGSP_BERT());
      physicsActivated["qgsp_bert"] = true;
    }
}

void BDSModularPhysicsList::QGSPBERTHP()
{
  ConstructAllLeptons();
  if(!physicsActivated["qgsp_bert_hp"])
    {
      constructors.push_back(new G4HadronPhysicsQGSP_BERT_HP());
      physicsActivated["qgsp_bert_hp"] = true;
    }
}

void BDSModularPhysicsList::QGSPBIC()
{
  ConstructAllLeptons();
  if(!physicsActivated["qgsp_bic"])
    {
      constructors.push_back(new G4HadronPhysicsQGSP_BIC());
      physicsActivated["qgsp_bic"] = true;
    }
}

void BDSModularPhysicsList::QGSPBICHP()
{
  ConstructAllLeptons();
  if(!physicsActivated["qgsp_bic_hp"])
    {
      constructors.push_back(new G4HadronPhysicsQGSP_BIC_HP());
      physicsActivated["qgsp_bic_hp"] = true;
    }
}

void BDSModularPhysicsList::Shielding()
{
  if(!physicsActivated["shielding"])
    {
      constructors.push_back(new G4HadronPhysicsShielding());
      physicsActivated["shielding"] = true;
    }
}

void BDSModularPhysicsList::Stopping()
{
  ConstructAllShortLived();
  ConstructAllIons();
  if(!physicsActivated["stopping"])
    {
      constructors.push_back(new G4StoppingPhysics());
      physicsActivated["stopping"] = true;
    }
}

void BDSModularPhysicsList::SynchRad()
{
  ConstructAllLeptons();
  if(!physicsActivated["synch_rad"])
    {
      constructors.push_back(new BDSPhysicsSynchRad());
      physicsActivated["synch_rad"] = true;
    }
}

#if G4VERSION_NUMBER > 1019
void BDSModularPhysicsList::EmGS()
{
  ConstructAllLeptons();
  if (!physicsActivated["em_gs"])
    {
      constructors.push_back(new G4EmStandardPhysicsGS());
      physicsActivated["em_gs"] = true;
    }
}
#endif

#if G4VERSION_NUMBER > 1020
void BDSModularPhysicsList::DecaySpin()
{
  if(!physicsActivated["decay_spin"])
    {// this will replace regular decay for various processes
      constructors.push_back(new G4SpinDecayPhysics());
      physicsActivated["decay_spin"] = true;
    }
}
#endif

#if G4VERSION_NUMBER > 1022
void BDSModularPhysicsList::IonPHP()
{
  BDS::CheckHighPrecisionDataExists("ion_php");
  ConstructAllBaryons();
  ConstructAllIons();
  ConstructAllLeptons();
  ConstructAllMesons();
  ConstructAllShortLived();
  
  if (!physicsActivated["ion_php"])
    {
      constructors.push_back(new G4IonPhysicsPHP());
      physicsActivated["ion_php"] = true;
    }
}
#endif

#if G4VERSION_NUMBER > 1029
void BDSModularPhysicsList::DecayMuonicAtom()
{
  ConstructAllLeptons();
  if (!physicsActivated["decay_muonic_atom"])
    {
      constructors.push_back(new G4MuonicAtomDecayPhysics());
      physicsActivated["decay_muonic_atom"] = true;
    }
}
#endif

#if G4VERSION_NUMBER > 1039
void BDSModularPhysicsList::ShieldingLEND()
{
  BDS::CheckLowEnergyNeutronDataExists("shielding_lend");
  if(!physicsActivated["shielding_lend"])
    {
      constructors.push_back(new G4HadronPhysicsShieldingLEND());
      physicsActivated["shielding_lend"] = true;
    }
}
#endif

void BDSModularPhysicsList::BuildAndAttachBiasWrapper(const GMAD::FastList<GMAD::PhysicsBiasing>& biases)
{
  // particles we know we can bias
  std::map<G4String, G4bool> particlesToBias =
    {
      {"e-"     , false},
      {"e+"     , false},
      {"gamma"  , false},
      {"proton" , false},
      {"mu-"    , false},
      {"mu+"    , false},
      {"pi-"    , false},
      {"pi+"    , false}
    };

  // iterate through bias structures and turn on biasing for that particle if it's in the
  // map of acceptable particle definitions.
  for (auto b : biases)
    {
      G4String name = G4String(b.particle);
      if (particlesToBias.find(name) != particlesToBias.end())
	{particlesToBias[name] = true;}
      else
	{
	  G4cerr << __METHOD_NAME__ << "Not possible to bias \"" << name << "\"" << G4endl;
	  exit(1);
	}
    }

  // check whether we need to construct or attach biasing at all
  typedef std::pair<const G4String, G4bool> mapvalue;
  G4bool anyBiases = std::any_of(particlesToBias.begin(),
				 particlesToBias.end(),
				 [](mapvalue i){return i.second;});

  if (!anyBiases)
    {return;}

  // there are biases
  G4GenericBiasingPhysics* physBias = new G4GenericBiasingPhysics();
  for (auto part : particlesToBias)
    {
      if (part.second)
        {
          G4cout << __METHOD_NAME__ << "wrapping \"" << part.first << "\" for biasing" << G4endl;
          physBias->Bias(part.first);
        }
    }
  RegisterPhysics(physBias);
}
