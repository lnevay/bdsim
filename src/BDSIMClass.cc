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
#include "BDSIMClass.hh"

#include "BDSExecOptions.hh"     // executable command line options 
#include "BDSGlobalConstants.hh" //  global parameters

#include <cstdlib>      // standard headers 
#include <cstdio>
#include <signal.h>

#include "G4EventManager.hh" // Geant4 includes
#include "G4GenericBiasingPhysics.hh"
#include "G4GeometryManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4ParallelWorldPhysics.hh"
#include "G4ParticleDefinition.hh"
#include "G4SteppingManager.hh"
#include "G4TrackingManager.hh"
#include "G4VModularPhysicsList.hh"

#include "BDSAcceleratorModel.hh"
#include "BDSBunch.hh"
#include "BDSBunchFactory.hh"
#include "BDSColours.hh"
#include "BDSDetectorConstruction.hh"   
#include "BDSEventAction.hh"
#include "BDSFieldFactory.hh"
#include "BDSFieldLoader.hh"
#include "BDSGeometryFactorySQL.hh"
#include "BDSGeometryWriter.hh"
#include "BDSMaterials.hh"
#include "BDSOutput.hh" 
#include "BDSOutputFactory.hh"
#include "BDSParallelWorldUtilities.hh"
#include "BDSParser.hh" // Parser
#include "BDSParticleDefinition.hh"
#include "BDSPhysicsUtilities.hh"
#include "BDSPrimaryGeneratorAction.hh"
#include "BDSRandom.hh" // for random number generator from CLHEP
#include "BDSRunAction.hh"
#include "BDSRunManager.hh"
#include "BDSSamplerRegistry.hh"
#include "BDSSDManager.hh"
#include "BDSSteppingAction.hh"
#include "BDSStackingAction.hh"
#include "BDSTemporaryFiles.hh"
#include "BDSTrackingAction.hh"
#include "BDSUtilities.hh"
#include "BDSVisManager.hh"

BDSIM::BDSIM():
  ignoreSIGINT(false),
  usualPrintOut(true),
  initialised(false),
  argcCache(0),
  argvCache(nullptr),
  parser(nullptr),
  bdsOutput(nullptr),
  bdsBunch(nullptr),
  runManager(nullptr)
{;}

BDSIM::BDSIM(int argc, char** argv, bool usualPrintOutIn):
  ignoreSIGINT(false),
  usualPrintOut(usualPrintOutIn),
  initialised(false),
  argcCache(argc),
  argvCache(argv),
  parser(nullptr),
  bdsOutput(nullptr),
  bdsBunch(nullptr),
  runManager(nullptr)
{
  Initialise();
}  

int BDSIM::Initialise()
{
  /// Print header & program information
  G4cout<<"bdsim : version @BDSIM_VERSION@"<<G4endl;
  G4cout<<"        (C) 2001-@CURRENT_YEAR@ Royal Holloway University London"<<G4endl;
  G4cout<<G4endl;
  G4cout<<"        Reference: https://arxiv.org/abs/1808.10745"<<G4endl;
  G4cout<<"        Website:   http://www.pp.rhul.ac.uk/bdsim"<<G4endl;
  G4cout<<G4endl;

  /// Initialize executable command line options reader object
  const BDSExecOptions* execOptions = new BDSExecOptions(argcCache,argvCache);
  if (usualPrintOut)
    {execOptions->Print();}
  ignoreSIGINT = execOptions->IgnoreSIGINT(); // different sig catching for cmake
  
#ifdef BDSDEBUG
  G4cout << __FUNCTION__ << "> DEBUG mode is on." << G4endl;
#endif

  /// Parse lattice file
  parser = BDSParser::Instance(execOptions->InputFileName());
  /// Update options generated by parser with those from executable options.
  parser->AmalgamateOptions(execOptions->Options());
  parser->AmalgamateBeam(execOptions->Beam(), execOptions->Options().recreate);
  /// Check options for consistency
  parser->CheckOptions();
  
  /// Explicitly initialise materials to construct required materials before global constants.
  BDSMaterials::Instance()->PrepareRequiredMaterials(execOptions->Options().verbose);

  /// No longer needed. Everything can safely use BDSGlobalConstants from now on.
  delete execOptions; 

  /// Force construction of global constants after parser has been initialised (requires
  /// materials first). This uses the options and beam from BDSParser.
  /// Non-const as we'll update the particle definition.
  BDSGlobalConstants* globalConstants = BDSGlobalConstants::Instance();

  /// Initialize random number generator
  BDSRandom::CreateRandomNumberGenerator();
  BDSRandom::SetSeed(); // set the seed from options

  /// Construct output
  bdsOutput = BDSOutputFactory::CreateOutput(globalConstants->OutputFormat(),
					     globalConstants->OutputFileName());

  /// Check geant4 exists in the current environment
  if (!BDS::Geant4EnvironmentIsSet())
    {G4cout << "No Geant4 environmental variables found - please source geant4.sh environment" << G4endl; return 1;}

  /// Construct mandatory run manager (the G4 kernel) and
  /// register mandatory initialization classes.
  runManager = new BDSRunManager;

  /// Register the geometry and parallel world construction methods with run manager.
  BDSDetectorConstruction* realWorld = new BDSDetectorConstruction();
  /// Here the geometry isn't actually constructed - this is called by the runManager->Initialize()
  auto samplerWorlds = BDS::ConstructAndRegisterParallelWorlds(realWorld);
  runManager->SetUserInitialization(realWorld);  

  /// For geometry sampling, phys list must be initialized before detector.
  /// BUT for samplers we use a parallel world and this HAS to be before the physcis
#ifdef BDSDEBUG 
  G4cout << __FUNCTION__ << "> Constructing physics processes" << G4endl;
#endif
  G4String physicsListName = parser->GetOptions().physicsList;

  // sampler physics process for parallel world tracking must be instantiated BEFORE
  // regular physics.
  // Note, we purposively don't create a parallel world process for the curvilinear
  // world as we don't need the track information from it - unreliable that way. We
  // query the geometry directly using our BDSAuxiliaryNavigator class.
  auto samplerPhysics = BDS::ConstructSamplerParallelPhysics(samplerWorlds);
  G4VModularPhysicsList* physList = BDS::BuildPhysics(physicsListName);
  
  // Construction of the physics lists defines the necessary particles and therefore
  // we can calculate the beam rigidity for the particle the beam is designed w.r.t. This
  // must happen before the geometry is constructed (which is called by
  // runManager->Initialize()).
  BDSParticleDefinition* beamParticle;
  beamParticle = BDS::ConstructBeamParticle(globalConstants->ParticleName(),
					    globalConstants->BeamTotalEnergy(),
					    globalConstants->FFact());
  globalConstants->SetBeamParticleDefinition(beamParticle);
  if (usualPrintOut)
    {G4cout << "main> Beam particle properties: " << G4endl << *beamParticle;}
  // update rigidity where needed
  realWorld->SetRigidityForConstruction(beamParticle->BRho());
  realWorld->SetBeta0ForConstruction(beamParticle->Beta());
  BDSFieldFactory::SetDefaultRigidity(beamParticle->BRho());       // used for field loading
  BDSGeometryFactorySQL::SetDefaultRigidity(beamParticle->BRho()); // used for sql field loading
  
  BDS::RegisterSamplerPhysics(samplerPhysics, physList);
  auto biasPhysics = BDS::BuildAndAttachBiasWrapper(parser->GetBiasing());
  if (biasPhysics)//could be nullptr and can't be passed to geant4 like this
    {physList->RegisterPhysics(biasPhysics);}
  runManager->SetUserInitialization(physList);

  /// Instantiate the specific type of bunch distribution.
  bdsBunch = BDSBunchFactory::CreateBunch(beamParticle,
					  parser->GetBeam(),
					  globalConstants->BeamlineTransform(),
					  globalConstants->BeamlineS(),
					  globalConstants->GeneratePrimariesOnly());
  
  /// Optionally generate primaries only and exit
  /// Unfortunately, this has to be here as we can't query the geant4 particle table
  /// until after the physics list has been constructed and attached a run manager.
  if (globalConstants->GeneratePrimariesOnly())
    {
      const G4int pdgID = beamParticle->ParticleDefinition()->GetPDGEncoding();
      const G4double charge = beamParticle->Charge(); // note this may be different for PDG charge for an ion
      
      // output creation is duplicated below but with this if loop, we exit so ok.
      bdsOutput->NewFile();
      const G4int nToGenerate = globalConstants->NGenerate();
      const G4int printModulo = globalConstants->PrintModuloEvents();
      bdsBunch->BeginOfRunAction(nToGenerate);
      for (G4int i = 0; i < nToGenerate; i++)
	{
	  if (i%printModulo == 0)
	    {G4cout << "\r Primary> " << std::fixed << i << " of " << nToGenerate << G4endl;}
	  auto coords = bdsBunch->GetNextParticle();
	  bdsOutput->FillEventPrimaryOnly(coords, charge, pdgID);
	}
      // Write options now file open.
      const GMAD::OptionsBase* ob = BDSParser::Instance()->GetOptionsBase();
      bdsOutput->FillOptions(ob);

      // Write beam
      const GMAD::BeamBase* bb = BDSParser::Instance()->GetBeamBase();
      bdsOutput->FillBeam(bb);

      bdsOutput->CloseFile();
      delete bdsBunch;
      delete bdsOutput;
      delete runManager;
      return 0;
    }
  
  /// Print the geometry tolerance
  G4GeometryTolerance* theGeometryTolerance = G4GeometryTolerance::GetInstance();
  if (usualPrintOut)
    {
      G4cout << __FUNCTION__ << "> Geometry Tolerances: "     << G4endl;
      G4cout << __FUNCTION__ << ">" << std::setw(22) << "Surface: " << std::setw(10) << theGeometryTolerance->GetSurfaceTolerance() << " mm"   << G4endl;
      G4cout << __FUNCTION__ << ">" << std::setw(22) << "Angular: " << std::setw(10) << theGeometryTolerance->GetAngularTolerance() << " rad"  << G4endl;
      G4cout << __FUNCTION__ << ">" << std::setw(22) << "Radial: "  << std::setw(10) << theGeometryTolerance->GetRadialTolerance()  << " mm"   << G4endl;
    }
  /// Set user action classes
#ifdef BDSDEBUG 
  G4cout << __FUNCTION__ << "> Registering user action - Run Action"<<G4endl;
#endif
  runManager->SetUserAction(new BDSRunAction(bdsOutput, bdsBunch, beamParticle->IsAnIon()));

#ifdef BDSDEBUG 
  G4cout << __FUNCTION__ << "> Registering user action - Event Action"<<G4endl;
#endif
  BDSEventAction* eventAction = new BDSEventAction(bdsOutput);
  runManager->SetUserAction(eventAction);

#ifdef BDSDEBUG 
  G4cout << __FUNCTION__ << "> Registering user action - Stepping Action"<<G4endl;
#endif
  // Only add steppingaction if it is actually used, so do check here (for performance reasons)
  if (globalConstants->VerboseStep())
    {
      G4int verboseEventNumber = globalConstants->VerboseEventNumber();
      runManager->SetUserAction(new BDSSteppingAction(true, verboseEventNumber));
    }
  
#ifdef BDSDEBUG 
  G4cout << __FUNCTION__ << "> Registering user action - Tracking Action"<<G4endl;
#endif
  runManager->SetUserAction(new BDSTrackingAction(globalConstants->Batch(),
						  globalConstants->StoreTrajectory(),
						  globalConstants->TrajNoTransportation(),
						  eventAction));

#ifdef BDSDEBUG 
  G4cout << __FUNCTION__ << "> Registering user action - Stacking Action"<<G4endl;
#endif
  runManager->SetUserAction(new BDSStackingAction(globalConstants));

#ifdef BDSDEBUG 
  G4cout << __FUNCTION__ << "> Registering user action - Primary Generator"<<G4endl;
#endif
  runManager->SetUserAction(new BDSPrimaryGeneratorAction(bdsBunch, beamParticle));

  /// Initialize G4 kernel
#ifdef BDSDEBUG 
  G4cout << __FUNCTION__ << "> Initialising Geant4 kernel"<<G4endl;
#endif
  runManager->Initialize();

  /// Implement bias operations on all volumes only after G4RunManager::Initialize()
  realWorld->BuildPhysicsBias();

  if (usualPrintOut && BDSGlobalConstants::Instance()->PhysicsVerbose())
    {
      BDS::PrintPrimaryParticleProcesses(beamParticle->Name());
      BDS::PrintDefinedParticles();
    }

  /// Set verbosity levels
  runManager->SetVerboseLevel(globalConstants->VerboseRunLevel());
  G4EventManager::GetEventManager()->SetVerboseLevel(globalConstants->VerboseEventLevel());
  G4EventManager::GetEventManager()->GetTrackingManager()->SetVerboseLevel(globalConstants->VerboseTrackingLevel());
  G4EventManager::GetEventManager()->GetTrackingManager()->GetSteppingManager()->SetVerboseLevel(globalConstants->VerboseSteppingLevel());
  
  /// Close the geometry in preparation for running - everything is now fixed.
  G4bool bCloseGeometry = G4GeometryManager::GetInstance()->CloseGeometry();
  if(!bCloseGeometry)
    { 
      G4cerr << __FUNCTION__ << "> error - geometry not closed." << G4endl;
      return 1;
    }

  if (globalConstants->ExportGeometry())
    {
      BDSGeometryWriter geometrywriter;
      geometrywriter.ExportGeometry(globalConstants->ExportType(),
				    globalConstants->ExportFileName());
    }

  initialised = true;
  return 0;
}

void BDSIM::BeamOn(int nGenerate)
{
  G4cout.precision(10);
  /// Catch aborts to close output stream/file. perhaps not all are needed.
  struct sigaction act;
  act.sa_handler = &BDS::HandleAborts;
  sigemptyset(&act.sa_mask);
  act.sa_flags = 0;
  if (!ignoreSIGINT)
    {sigaction(SIGINT,  &act, 0);}
  sigaction(SIGABRT, &act, 0);
  sigaction(SIGTERM, &act, 0);
  sigaction(SIGSEGV, &act, 0);
  
  /// Run in either interactive or batch mode
  if(!BDSGlobalConstants::Instance()->Batch())   // Interactive mode
    {
      BDSVisManager visManager;
      visManager.StartSession(argcCache, argvCache);
    }
  else
    {// batch mode
      if (nGenerate < 0)
	{runManager->BeamOn(BDSGlobalConstants::Instance()->NGenerate());}
      else
	{runManager->BeamOn(nGenerate);}
    }
}

BDSIM::~BDSIM()
{
  /// Termination & clean up.
  G4GeometryManager::GetInstance()->OpenGeometry();
  
#ifdef BDSDEBUG 
  G4cout << __FUNCTION__ << "> BDSOutput deleting..."<<G4endl;
#endif
  delete bdsOutput;
  
#ifdef BDSDEBUG 
  G4cout << __FUNCTION__ << "> instances deleting..."<<G4endl;
#endif
  // Order important here because of singletons relying on each other
  delete BDSAcceleratorModel::Instance();
  delete BDSTemporaryFiles::Instance();
  delete BDSFieldFactory::Instance(); // this uses BDSGlobalConstants which uses BDSMaterials
  delete BDSGlobalConstants::Instance();
  delete BDSMaterials::Instance();

  // instances not used in this file, but no other good location for deletion
  delete BDSColours::Instance();
  delete BDSFieldLoader::Instance();
  delete BDSSDManager::Instance();
  delete BDSSamplerRegistry::Instance();
  
#ifdef BDSDEBUG 
  G4cout<< __FUNCTION__ << "> BDSRunManager deleting..."<<G4endl;
#endif
  delete runManager;
  delete bdsBunch;
  delete parser;

  if (usualPrintOut)
    {G4cout << __FUNCTION__ << "> End of Run. Thank you for using BDSIM!" << G4endl;}
}
