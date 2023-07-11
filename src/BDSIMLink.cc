/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2023.

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
#include "BDSIMLink.hh"

#include <csignal>
#include <cstdlib>      // standard headers 
#include <cstdio>

#include "G4EventManager.hh" // Geant4 includes
#include "G4GeometryManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4Version.hh"
#include "G4VModularPhysicsList.hh"

#include "BDSAcceleratorModel.hh"
#include "BDSAperturePointsLoader.hh"
#include "BDSBeamPipeFactory.hh"
#include "BDSBunch.hh"
#include "BDSBunchFactory.hh"
#include "BDSBunchSixTrackLink.hh"
#include "BDSCavityFactory.hh"
#include "BDSColours.hh"
#include "BDSComponentFactoryUser.hh"
#include "BDSDebug.hh"
#include "BDSException.hh"
#include "BDSExecOptions.hh"
#include "BDSFieldFactory.hh"
#include "BDSFieldLoader.hh"
#include "BDSGeometryFactory.hh"
#include "BDSGeometryFactorySQL.hh"
#include "BDSGeometryWriter.hh"
#include "BDSGlobalConstants.hh"
#include "BDSLinkComponent.hh"
#include "BDSLinkDetectorConstruction.hh"
#include "BDSLinkEventAction.hh"
#include "BDSLinkPrimaryGeneratorAction.hh"
#include "BDSLinkRunAction.hh"
#include "BDSLinkRunManager.hh"
#include "BDSLinkStackingAction.hh"
#include "BDSLinkTrackingAction.hh"
#include "BDSMaterials.hh"
#include "BDSOutput.hh" 
#include "BDSOutputFactory.hh"
#include "BDSParallelWorldUtilities.hh"
#include "BDSParser.hh"
#include "BDSParticleExternal.hh"
#include "BDSParticleDefinition.hh"
#include "BDSPhysicsUtilities.hh"
#include "BDSRandom.hh"
#include "BDSSamplerRegistry.hh"
#include "BDSSDManager.hh"
#include "BDSTemporaryFiles.hh"
#include "BDSUtilities.hh"
#include "BDSVisManager.hh"

#include <map>
#include <set>

BDSIMLink::BDSIMLink(BDSBunch* bunchIn):
  ignoreSIGINT(false),
  usualPrintOut(true),
  initialised(false),
  initialisationResult(1),
  argcCache(0),
  argvCache(nullptr),
  parser(nullptr),
  bdsOutput(nullptr),
  bdsBunch(bunchIn),
  runManager(nullptr),
  construction(nullptr),
  runAction(nullptr),
  currentElementIndex(0),
  userPhysicsList(nullptr)
{;}

BDSIMLink::BDSIMLink(int argc, char** argv, bool usualPrintOutIn):
  ignoreSIGINT(false),
  usualPrintOut(usualPrintOutIn),
  initialised(false),
  initialisationResult(1),
  argcCache(argc),
  argvCache(argv),
  parser(nullptr),
  bdsOutput(nullptr),
  bdsBunch(nullptr),
  runManager(nullptr),
  construction(nullptr),
  runAction(nullptr),
  currentElementIndex(0),
  userPhysicsList(nullptr)
{
  initialisationResult = Initialise();
}

int BDSIMLink::Initialise(int argc,
    char** argv,
    bool   usualPrintOutIn,
    double minimumKineticEnergy,
    bool   protonsAndIonsOnly)
{
  argcCache = argc;
  argvCache = argv;
  usualPrintOut = usualPrintOutIn;
  initialisationResult = Initialise(minimumKineticEnergy, protonsAndIonsOnly);
  return initialisationResult;
}

int BDSIMLink::Initialise(double minimumKineticEnergy,
                          bool   protonsAndIonsOnly)
{
  minimumKineticEnergy *= CLHEP::GeV; // units
  G4bool trackerDebug = false;

  /// Initialize executable command line options reader object
  const BDSExecOptions* execOptions = new BDSExecOptions(argcCache,argvCache);
  if (usualPrintOut)
    {execOptions->Print();}
  ignoreSIGINT = execOptions->IgnoreSIGINT(); // different sig catching for cmake
  execOptions->PrintCopyright();
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << "DEBUG mode is on." << G4endl;
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
  BDSRandom::CreateRandomNumberGenerator(globalConstants->RandomEngine());
  BDSRandom::SetSeed(); // set the seed from options

  /// Construct output
  bdsOutput = BDSOutputFactory::CreateOutput(globalConstants->OutputFormat(),
                                             globalConstants->OutputFileName());

  /// Check geant4 exists in the current environment
  if (!BDS::Geant4EnvironmentIsSet())
    {
      G4cerr << "No Geant4 environmental variables found - please source geant4.sh environment" << G4endl;
      G4cout << "A common fault is the wrong Geant4 environment as compared to the one BDSIM was compiled with." << G4endl;
      return 1;
    }

  /// Construct mandatory run manager (the G4 kernel) and
  /// register mandatory initialization classes.
  runManager = new BDSLinkRunManager();

  /// Register the geometry and parallel world construction methods with run manager.
  construction = new BDSLinkDetectorConstruction();

  /// Here the geometry isn't actually constructed - this is called by the runManager->Initialize()
  auto parallelWorldsRequiringPhysics = BDS::ConstructAndRegisterParallelWorlds(construction, true, false);

  runManager->SetUserInitialization(construction);

  // Set filters used in sensitive detectors that transfer particles back
  BDSSDManager::Instance()->SetLinkMinimumEK(minimumKineticEnergy);
  BDSSDManager::Instance()->SetLinkProtonsAndIonsOnly(protonsAndIonsOnly);
  
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << "> Constructing physics processes" << G4endl;
#endif
  G4String physicsListName = parser->GetOptions().physicsList;

#if G4VERSION_NUMBER > 1049
  // from 10.5 onwards they have a looping particle killer that warnings and kills particles
  // deemed to be looping that are <100 MeV. This is unrelated to the primary energy so troublesome.
  // set to the 'low' limits here ~10keV. This must be done before any physics is created as the
  // parameters are copied into the transportation physics process for each particle and it's very
  // hard to sift through and fix afterwards
  G4PhysicsListHelper::GetPhysicsListHelper()->UseLowLooperThresholds();
#endif
  auto parallelWorldPhysics = BDS::ConstructParallelWorldPhysics(parallelWorldsRequiringPhysics);
  G4int physicsVerbosity = BDSGlobalConstants::Instance()->PhysicsVerbosity();
  G4VModularPhysicsList* physList;
  if (userPhysicsList)
    {
      G4cout << "Using externally registered user defined physics list" << G4endl;
      physList = userPhysicsList;
    }
  else
    {physList = BDS::BuildPhysics(physicsListName, physicsVerbosity);}
  
  BDS::RegisterSamplerPhysics(parallelWorldPhysics, physList);
  // Construction of the physics lists defines the necessary particles and therefore
  // we can calculate the beam rigidity for the particle the beam is designed w.r.t. This
  // must happen before the geometry is constructed (which is called by
  // runManager->Initialize()).
  BDSParticleDefinition* designParticle = nullptr;
  BDSParticleDefinition* beamParticle = nullptr;
  G4bool beamDifferentFromDesignParticle = false;
  BDS::ConstructDesignAndBeamParticle(BDSParser::Instance()->GetBeam(),
				      globalConstants->FFact(),
				      designParticle,
				      beamParticle,
				      beamDifferentFromDesignParticle);
  if (usualPrintOut)
    {
      G4cout << "Design particle properties: " << G4endl << *designParticle;
      if (beamDifferentFromDesignParticle)
	{G4cout << "Beam particle properties: " << G4endl << *beamParticle;}
    }
  // update rigidity where needed
  construction->SetDesignParticle(designParticle);
  //BDSFieldFactory::SetDesignParticle(designParticle);
  
  //auto biasPhysics = BDS::BuildAndAttachBiasWrapper(parser->GetBiasing());
  //if (biasPhysics)//could be nullptr and can't be passed to geant4 like this
  //  {physList->RegisterPhysics(biasPhysics);}
  runManager->SetUserInitialization(physList);

  /// Instantiate the specific type of bunch distribution.
  GMAD::Beam b;
  b.distrType = "sixtracklink";
  if (!bdsBunch)
    {
      bdsBunch = BDSBunchFactory::CreateBunch(beamParticle,
                                              parser->GetBeam(),
                                              globalConstants->BeamlineTransform(),
                                              globalConstants->BeamlineS(),
                                              globalConstants->GeneratePrimariesOnly());
    }
  else
    {
      bdsBunch->SetOptions(beamParticle, b, BDSBunchType::reference); // update particle definition only
    }
  /// We no longer need beamParticle so delete it to avoid confusion. The definition is
  /// held inside bdsBunch (can be updated dynamically).
  delete beamParticle;
  
  /// Print the geometry tolerance
  G4GeometryTolerance* theGeometryTolerance = G4GeometryTolerance::GetInstance();
  if (usualPrintOut)
    {
      G4cout << __METHOD_NAME__ << "Geometry Tolerances: "     << G4endl;
      G4cout << __METHOD_NAME__ << std::setw(12) << "Surface: " << std::setw(7) << theGeometryTolerance->GetSurfaceTolerance() << " mm"   << G4endl;
      G4cout << __METHOD_NAME__ << std::setw(12) << "Angular: " << std::setw(7) << theGeometryTolerance->GetAngularTolerance() << " rad"  << G4endl;
      G4cout << __METHOD_NAME__ << std::setw(12) << "Radial: "  << std::setw(7) << theGeometryTolerance->GetRadialTolerance()  << " mm"   << G4endl;
    }
  /// Set user action classes
  runAction = new BDSLinkRunAction();
  BDSLinkEventAction* eventAction = new BDSLinkEventAction(bdsOutput, runAction, trackerDebug);
  runManager->SetUserAction(eventAction);
  runManager->SetUserAction(runAction);
  G4int verboseSteppingEventStart = globalConstants->VerboseSteppingEventStart();
  G4int verboseSteppingEventStop  = BDS::VerboseEventStop(verboseSteppingEventStart,
                                                          globalConstants->VerboseSteppingEventContinueFor());
  runManager->SetUserAction(new BDSLinkTrackingAction(globalConstants->Batch(),
                                                      eventAction,
                                                      verboseSteppingEventStart,
                                                      verboseSteppingEventStop,
                                                      globalConstants->VerboseSteppingPrimaryOnly(),
                                                      globalConstants->VerboseSteppingLevel()));
  runManager->SetUserAction(new BDSLinkStackingAction(globalConstants, std::set<G4int>(), protonsAndIonsOnly, minimumKineticEnergy));
  
  /*
  // Only add steppingaction if it is actually used, so do check here (for performance reasons)
  G4int verboseSteppingEventStart = globalConstants->VerboseSteppingEventStart();
  G4int verboseSteppingEventStop  = BDS::VerboseEventStop(verboseSteppingEventStart,
							  globalConstants->VerboseSteppingEventContinueFor());
  if (globalConstants->VerboseSteppingBDSIM())
    {
      runManager->SetUserAction(new BDSSteppingAction(true,
						      verboseSteppingEventStart,
						      verboseSteppingEventStop));
    }
  */
  
  auto primaryGeneratorAction = new BDSLinkPrimaryGeneratorAction(bdsBunch, &currentElementIndex, construction, trackerDebug);
  construction->SetPrimaryGeneratorAction(primaryGeneratorAction);
  runManager->SetUserAction(primaryGeneratorAction);
  //BDSFieldFactory::SetPrimaryGeneratorAction(primaryGeneratorAction);

  /// Initialize G4 kernel
  runManager->Initialize();
  
  /// Implement bias operations on all volumes only after G4RunManager::Initialize()
  construction->BuildPhysicsBias();

  if (usualPrintOut && BDSGlobalConstants::Instance()->PhysicsVerbose())
    {
      BDS::PrintPrimaryParticleProcesses(bdsBunch->ParticleDefinition()->Name());
      BDS::PrintDefinedParticles();
    }

  /// Set verbosity levels at run and G4 event level. Per event and stepping are controlled
  /// in event, tracking and stepping action. These have to be done here due to the order
  /// of construction in Geant4.
  runManager->SetVerboseLevel(std::min(globalConstants->VerboseRunLevel(), globalConstants->PhysicsVerbosity()));
  G4EventManager::GetEventManager()->SetVerboseLevel(globalConstants->VerboseEventLevel());
  G4EventManager::GetEventManager()->GetTrackingManager()->SetVerboseLevel(globalConstants->VerboseTrackingLevel());
  
  /// Close the geometry in preparation for running - everything is now fixed.
  G4bool bCloseGeometry = G4GeometryManager::GetInstance()->CloseGeometry();
  if (!bCloseGeometry)
    { 
      G4cerr << __METHOD_NAME__ << "error - geometry not closed." << G4endl;
      return 1;
    }

  if (globalConstants->ExportGeometry())
    {
      BDSGeometryWriter geometrywriter;
      geometrywriter.ExportGeometry(globalConstants->ExportType(),
				    globalConstants->ExportFileName());
    }

  const auto& nameInds = construction->NameToElementIndex();
  nameToElementIndex.insert(nameInds.begin(), nameInds.end());
  
  if (bdsOutput)
    {
      bdsOutput->InitialiseGeometryDependent();
      bdsOutput->NewFile();
    }

  initialised = true;
  return 0;
}

void BDSIMLink::BeamOn(int nGenerate)
{
  if (initialisationResult > 1 || !initialised)
    {return;} // a mode where we don't do anything

  G4cout.precision(10);
  /// Catch aborts to close output stream/file. perhaps not all are needed.
  struct sigaction act;
  act.sa_handler = &BDS::HandleAborts;
  sigemptyset(&act.sa_mask);
  act.sa_flags = 0;
  if (!ignoreSIGINT)
    {sigaction(SIGINT,  &act, nullptr);}
  sigaction(SIGABRT, &act, nullptr);
  sigaction(SIGTERM, &act, nullptr);
  sigaction(SIGSEGV, &act, nullptr);
  
  /// Run in either interactive or batch mode
  try
    {
      if (!BDSGlobalConstants::Instance()->Batch())   // Interactive mode
	{
	  BDSVisManager visManager = BDSVisManager(BDSGlobalConstants::Instance()->VisMacroFileName(),
						   BDSGlobalConstants::Instance()->Geant4MacroFileName());
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
  catch (const BDSException& exception)
    {
      // don't do this for now in case it's dangerous and we try tracking with open geometry
      //G4GeometryManager::GetInstance()->OpenGeometry();
      G4cout << exception.what() << G4endl;
      exit(1);
    } 
}

BDSIMLink::~BDSIMLink()
{
  /// Termination & clean up.
  G4GeometryManager::GetInstance()->OpenGeometry();
    
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << "deleting..." << G4endl;
#endif
  if (bdsOutput)
    {bdsOutput->CloseFile();}
  delete bdsOutput;

  try
    {
      // order important here because of singletons relying on each other
      delete BDSBeamPipeFactory::Instance();
      delete BDSCavityFactory::Instance();
      delete BDSGeometryFactory::Instance();
      delete BDSAcceleratorModel::Instance();
      delete BDSTemporaryFiles::Instance();
      delete BDSFieldFactory::Instance(); // this uses BDSGlobalConstants which uses BDSMaterials
      delete BDSGlobalConstants::Instance();
      delete BDSMaterials::Instance();
      
      // instances not used in this file, but no other good location for deletion
      if (initialisationResult < 2)
	{
	  delete BDSColours::Instance();
	  delete BDSFieldLoader::Instance();
	  //delete BDSSDManager::Instance();
	  delete BDSSamplerRegistry::Instance();
    BDSAperturePointsCache::Instance()->ClearCachedFiles();
	}
    }
  catch (...)
    {;} // ignore any exception as this is a destructor
  
  delete runManager;
  delete parser;

  if (usualPrintOut)
    {G4cout << __METHOD_NAME__ << "End of Run. Thank you for using BDSIM!" << G4endl;}
}

int BDSIMLink::GetLinkIndex(const std::string& elementName) const
{
  int result = -1;
  auto search = nameToElementIndex.find(elementName);
  if (search != nameToElementIndex.end())
    {result = search->second;}
  return result;
}

const BDSLinkComponent* BDSIMLink::GetLinkComponent(int linkID) const
{
  const BDSBeamline* bl = construction->LinkBeamline();
  if (!bl)
    {return nullptr;}
  if (linkID > (int)bl->size())
    {return nullptr;}
  const auto rawAccComponent = bl->at(linkID)->GetAcceleratorComponent();
  const auto linkComponent = dynamic_cast<const BDSLinkComponent*>(rawAccComponent);
  return linkComponent;
}

double BDSIMLink::GetChordLengthOfLinkElement(int beamlineIndex) const
{
  const BDSLinkComponent* component = GetLinkComponent(beamlineIndex);
  if (!component)
    {return -1.0;} // play it safe
  return component->ComponentChordLength();
}

double BDSIMLink::GetChordLengthOfLinkElement(const std::string& elementName)
{
  int linkID = GetLinkIndex(elementName);
  int beamlineIndex = linkIDToBeamlineIndex[linkID];
  return GetChordLengthOfLinkElement(beamlineIndex);
}

double BDSIMLink::GetArcLengthOfLinkElement(int beamlineIndex) const
{
  const BDSLinkComponent* component = GetLinkComponent(beamlineIndex);
  if (!component)
    {return -1.0;} // play it safe
  return component->ComponentArcLength();
}

double BDSIMLink::GetArcLengthOfLinkElement(const std::string& elementName)
{
  int linkID = GetLinkIndex(elementName);
  int beamlineIndex = linkIDToBeamlineIndex[linkID];
  return GetArcLengthOfLinkElement(beamlineIndex);
}

void BDSIMLink::SelectLinkElement(const std::string& elementName, G4bool debug)
{
  if (debug)
    {G4cout << "Searching for " << elementName;}
  currentElementIndex = GetLinkIndex(elementName);
  if (debug)
    {G4cout << ", Index " << currentElementIndex << G4endl;}
}

void BDSIMLink::SelectLinkElement(int index, G4bool debug)
{
  if (debug)
    {G4cout << "Searching for " << index << G4endl;}
  currentElementIndex = index;
}

int BDSIMLink::AddLinkCollimatorJaw(const std::string& collimatorName,
				     const std::string& materialName,
				     double length,
				     double halfApertureLeft,
				     double halfApertureRight,
				     double rotation,
				     double xOffset,
				     double yOffset,
                     double jawTiltLeft,
                     double jawTiltRight,
				     bool   buildLeftJaw,
				     bool   buildRightJaw,
				     bool   isACrystal,
				     double crystalAngle,
				     bool   sampleIn)
{
  G4GeometryManager* gm = G4GeometryManager::GetInstance();
  if (gm->IsGeometryClosed())
    {gm->OpenGeometry();}

  G4int linkID = construction->AddLinkCollimatorJaw(collimatorName,
				     materialName,
				     length,
				     halfApertureLeft,
				     halfApertureRight,
				     rotation,
				     xOffset,
				     yOffset,
                     jawTiltLeft,
                     jawTiltRight,
				     buildLeftJaw,
				     buildRightJaw,
				     isACrystal,
				     crystalAngle,
				     sampleIn);
  // update this class's nameToElementIndex map
  nameToElementIndex = construction->NameToElementIndex();
  linkIDToBeamlineIndex = construction->LinkIDToBeamlineIndex();
  
  if (bdsOutput)
    {bdsOutput->UpdateSamplers();}

  /// Close the geometry in preparation for running - everything is now fixed.
  G4bool bCloseGeometry = gm->CloseGeometry();
  if (!bCloseGeometry)
    {throw BDSException(__METHOD_NAME__, "error - geometry not closed.");}
  return (int)linkID;
}

BDSHitsCollectionSamplerLink* BDSIMLink::SamplerHits() const
{
  return runAction ? runAction->SamplerHits() : nullptr;
}

int BDSIMLink::GetCurrentMaximumSixTrackParticleID() const
{
  return runAction ? runAction->MaximumExternalParticleID() : 0;
}

void BDSIMLink::SetCurrentMaximumExternalParticleID(int currentMaximumExternalParticleID)
{
  if (runAction)
    {runAction->SetMaximumExternalParticleID(currentMaximumExternalParticleID);}
}
