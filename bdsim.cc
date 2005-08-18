/** BDSIM code.   
*   Author: Grahame A. Blair, Royal Holloway, Univ. of London.
*   @Copyright (c) 2002 by G.A.Blair.  ALL RIGHTS RESERVED. 
*   Update 6.04.2005 Agapov
*/

//  bdsim code 


const int DEBUG = 1;

#include "BDSGlobalConstants.hh" // global parameters


#include "G4UImanager.hh"        // G4 session managers
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#include "Randomize.hh"

#ifdef G4VIS_USE
#include "BDSVisManager.hh"
#endif


#include "stdlib.h"      // standard headers 
#include <stdio.h>
#include <unistd.h>
#include <getopt.h>

// tmp
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4HadronInelasticProcess.hh"  // processes


#include "BDSDetectorConstruction.hh"   // Geant4 includes
#include "BDSPhysicsList.hh"
#include "BDSPrimaryGeneratorAction.hh"
#include "BDSRunAction.hh"
#include "BDSEventAction.hh"
#include "BDSSteppingAction.hh"
#include "BDSUserTrackingAction.hh"
#include "BDSSteppingVerbose.hh"
#include "BDSRunManager.hh"


#include "BDSOutput.hh"  // interface to result output
#include "BDSBunch.hh"

#include "parser/gmad.h"  // GMAD parser




//=======================================================


// Print program usage
static void usage()
{
  G4cout<<"Usage: bdsim [options]"<<G4endl;
  G4cout<<"Options:"<<G4endl;
  G4cout<<"--file=<filename>   : specify the lattice file "<<G4endl
	<<"--output=<fmt>      : output format (root|ascii), default ascii"<<G4endl
	<<"--outfile=<file>    : output file name. Will be appended with _N"<<G4endl
        <<"                      where N = 0, 1, 2, 3... etc."<<G4endl
	<<"--vis_mac=<file>    : file with the visualization macro script, default vis.mac"<<G4endl
	<<"--help              : display this message"<<G4endl
	<<"--verbose           : run in verbose mode"<<G4endl
	<<"--batch             : batch mode - no graphics"<<G4endl;

}


BDSGlobalConstants* BDSGlobals;  // global options instance
BDSOutput bdsOutput;                // output interface
BDSBunch theBunch;  // bunch information
G4int outputFormat=_ASCII;
G4String outputFilename="output";  //receives a .txt or .root in BDSOutput
G4String inputFilename= "optics.gmad"; // input file with gmad lattice description
G4String visMacroFile="vis.mac"; // visualization macro file

G4bool verbose = false;  // run options
G4bool isBatch = false;



//=======================================================


int main(int argc,char** argv) {

  // NOTE: Visualization crashes when the argument parsing is added 


  // Parse the command line options 
  
   static struct option LongOptions[] = {
     { "help" , 0, 0, 0 },
     { "verbose", 0, 0, 0 },
     { "file", 1, 0, 0 },
     { "vis_mac", 1, 0, 0 },
     { "output", 1, 0, 0},
     { "outfile", 1, 0, 0},
     { "batch", 0, 0, 0 },  
     { 0, 0, 0, 0 }
   };
  
  int OptionIndex = 0;
  int c;
  int ThisOptionId;
  
  for(;;)
    {
      
      ThisOptionId = optind ? optind : 1;
      OptionIndex = 0;
      
      c = getopt_long(argc, argv, "Vv",
		      LongOptions, &OptionIndex );
      
      if ( c == -1 ) // end of options list
  	break;
      
      switch (c) {
      case 0:
	
	if( !strcmp(LongOptions[OptionIndex].name , "help") )
	  {

	    usage();
	    return 1;
	  }
	
	if( !strcmp(LongOptions[OptionIndex].name , "batch") )
	  {
	    isBatch = true;
	  }
	
	if( !strcmp(LongOptions[OptionIndex].name , "verbose") )
	  {
	    verbose = true; 
	  }
	if( !strcmp(LongOptions[OptionIndex].name , "output") )
	  {
	    if(optarg) {
	      if(!strcmp(optarg,"ascii")) outputFormat=_ASCII;
	      else if (!strcmp(optarg,"root")) outputFormat=_ROOT;
	      else G4cerr<<"unknown output format "<<optarg<<G4endl;
	    }
	  }
	if( !strcmp(LongOptions[OptionIndex].name , "outfile") )
	  {
	    if(optarg) {
	      outputFilename=optarg;
	    }
	  }
	if( !strcmp(LongOptions[OptionIndex].name , "file") )
	  {
	    if(optarg) {
	      inputFilename=optarg;
	    }
	    else {
	      G4cout<<"please specify the lattice filename"<<G4endl;
	    }
	  }
	if( !strcmp(LongOptions[OptionIndex].name , "vis_mac") )
	  {
	    if(optarg) {
	      visMacroFile=optarg;
	    }
	    else {
	      G4cout<<"please specify the visualization macro file"<<G4endl;
	    }
	  }
	
	
	break;
      case -1:
	break;
      default:
	break;
      }
      
    }
  
  
  G4cout<<"Using input file: "<<inputFilename<<G4endl;
    
  if( gmad_parser(inputFilename) == -1)   // parse lattice file
    {
      G4cout<<"can't open input file "<<inputFilename<<G4endl;
      exit(1);
    }


  BDSGlobals = new BDSGlobalConstants(options);
  theBunch.SetOptions(options);
 

  bdsOutput.SetFormat(outputFormat);


  // set default output formats:
  G4cout.precision(10);

  // choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);
  HepRandom::setTheSeed(BDSGlobals->GetRandomSeed());

  G4cout<<" seed from bdsglobals="<<BDSGlobals->GetRandomSeed()<<G4endl;
  G4cout<<"Random Number SEED ="<<HepRandom::getTheSeed()<<G4endl;

 
  BDSRunManager * runManager = new BDSRunManager;

  // set mandatory initialization classes
 
  BDSDetectorConstruction* detector = new BDSDetectorConstruction;

  if(DEBUG) G4cout<<"detector construction done"<<G4endl;
 
  runManager->SetUserInitialization(detector);

  if(DEBUG) G4cout<<"constructing phys list"<<G4endl;
 
  BDSPhysicsList* PhysList=new BDSPhysicsList;
  
  if(DEBUG) G4cout<<"user init phys list"<<G4endl;
  runManager->SetUserInitialization(PhysList);
  if(DEBUG) G4cout<<"user init phys list done"<<G4endl;

  G4UIsession* session=0;

 //  char **tmp = new char*[2];
//   tmp[0] = new char[32];
//   tmp[1] = new char[32];

//   strcpy(tmp[0],"hui");
//   strcpy(tmp[1],"vam");
  
  if(!isBatch)
    {
#ifdef G4UI_USE_XM
      if(DEBUG) G4cout<<"New Xm session"<<G4endl;
      session = new G4UIXm(0,tmp);
      if(DEBUG) G4cout<<"done"<<G4endl;
#else
#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);
#else
      session = new G4UIterminal();
#endif    
#endif
      
    }  


  if(DEBUG) G4cout<<"initializing visual manager"<<G4endl;


#ifdef G4VIS_USE  // visualization manager
  G4VisManager* visManager;
  if(!isBatch);
  {
    visManager = new BDSVisManager;
    visManager->Initialize();
  }
#endif
  
  
  // set user action classes
  if(DEBUG) G4cout<<"user action - detector"<<G4endl;
  runManager->SetUserAction(new BDSPrimaryGeneratorAction(detector));
  if(DEBUG) G4cout<<"user action - runaction"<<G4endl;

  runManager->SetUserAction(new BDSRunAction);

  
  //  BDSEventAction* theBDSEventAction = new BDSEventAction();
  if(DEBUG) G4cout<<"user action - eventaction"<<G4endl;
  runManager->SetUserAction(new BDSEventAction());

  if(DEBUG) G4cout<<"user action - steppingaction"<<G4endl;
  runManager->SetUserAction(new BDSSteppingAction);

  if(DEBUG) G4cout<<"user action - trackingaction"<<G4endl;
  runManager->SetUserAction(new BDSUserTrackingAction);
  

  if(DEBUG) G4cout<<"init kernel"<<G4endl;
  //Initialize G4 kernel
  runManager->Initialize();
  if(DEBUG) G4cout<<"init kernel done"<<G4endl;
  
  //  PhysList->BDSAddTransportation();
  
  // bias inelastic events
  
  G4ProcessManager *pManager = 0;
  G4ProcessVector *procVec = 0;
  G4String pName;
  G4HadronInelasticProcess* HadInProc=0;
  
  pManager = G4Electron::Electron()->GetProcessManager();
  procVec=pManager->GetProcessList();
  G4int nProc,iProc;
  nProc=pManager->GetProcessListLength();
  
  if(DEBUG) G4cout<<"electron - Hadron inelastic"<<G4endl;

  for(iProc=0;iProc<nProc;iProc++)
    {
      pName=(*procVec)[iProc]->GetProcessName();
      if(pName=="ElectroNuclear") 
	{
	  HadInProc=
	    static_cast<G4HadronInelasticProcess*>((*procVec)[iProc]);
	  HadInProc->
	    BiasCrossSectionByFactor(BDSGlobals->
				     GetHadronInelasticScaleFactor());
	}
    }
  
  pManager = G4Positron::Positron()->GetProcessManager();
  procVec=pManager->GetProcessList();
  nProc=pManager->GetProcessListLength();

  if(DEBUG) G4cout<<"positron - Hadron inelastic"<<G4endl;

  for(iProc=0;iProc<nProc;iProc++)
    {
      pName=(*procVec)[iProc]->GetProcessName();
      if(pName=="PositronNuclear")
	{
	  HadInProc=
	    static_cast<G4HadronInelasticProcess*>((*procVec)[iProc]);
	  HadInProc->
	    BiasCrossSectionByFactor(BDSGlobals->
				     GetHadronInelasticScaleFactor());
	}
    }
  
  
  pManager =G4Gamma::Gamma()->GetProcessManager();
  procVec=pManager->GetProcessList();
  nProc=pManager->GetProcessListLength();
  
  for(iProc=0;iProc<nProc;iProc++)
    {
      pName=(*procVec)[iProc]->GetProcessName();
      if(pName=="PhotoNuclear")
	{
	  HadInProc=
	    static_cast<G4HadronInelasticProcess*>((*procVec)[iProc]);
	  HadInProc->
	    BiasCrossSectionByFactor(BDSGlobals->
				     GetHadronInelasticScaleFactor());
	}
    }


  // neutrons - redefine the transportation process
  // TODO : interface to MCNP and Geant models
  
  
  
  bdsOutput.Init(0); // activate the output - setting the first filename to 
                     // be appended with _0

  if (!isBatch)   // Define UI session for interactive mode.
    {
      // get the pointer to the User Interface manager 
      G4UImanager* UI = G4UImanager::GetUIpointer();  
      // G4UIterminal is a (dumb) terminal.
      UI->ApplyCommand("/control/execute " + visMacroFile);    
#ifdef G4UI_USE_XM
      // Customize the G4UIXm menubar with a macro file :
      UI->ApplyCommand("/control/execute gui.mac");
#endif
      session->SessionStart();
      delete session;
      
    }
  else           // Batch mode
    { 
      runManager->BeamOn(BDSGlobals->GetNumberToGenerate());
    }
  
  // job termination  
#ifdef G4VIS_USE
  if(!isBatch) delete visManager;
#endif
  
    
  
  delete BDSGlobals;
  delete runManager;
     
  return 0;
}
