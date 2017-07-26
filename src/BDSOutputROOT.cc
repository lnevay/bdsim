#include "BDSOutputROOT.hh"

#include "parser/options.h"
#include "BDSDebug.hh"
#include "BDSOutputROOTEventHistograms.hh"
#include "BDSOutputROOTEventInfo.hh"
#include "BDSOutputROOTEventLoss.hh"
#include "BDSOutputROOTEventModel.hh"
#include "BDSOutputROOTEventOptions.hh"
#include "BDSOutputROOTEventRunInfo.hh"
#include "BDSOutputROOTEventSampler.hh"
#include "BDSOutputROOTEventTrajectory.hh"

#include "TFile.h"
#include "TObject.h"
#include "TTree.h"

BDSOutputROOT::BDSOutputROOT(G4String fileName,
			     G4int    fileNumberOffset):
  BDSOutput(fileName, fileNumberOffset)
{;}

BDSOutputROOT::~BDSOutputROOT()
{
  if (theRootOutputFile && theRootOutputFile->IsOpen())
    {theRootOutputFile->Write(0,TObject::kOverwrite);}
}

void BDSOutputROOT::NewFile() 
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ <<G4endl;
#endif
  G4String filename = GetNextFileName();
  
  theRootOutputFile      = new TFile(filename,"RECREATE", "BDS output file");
  // root file - note this sets the current 'directory' to this file!
  theRootOutputFile->cd();
  // options data tree
  theOptionsOutputTree   = new TTree("Options","BDSIM options");
  // model data tree
  theModelOutputTree     = new TTree("Model","BDSIM model");
  // run histogram tree
  theRunOutputTree       = new TTree("Run","BDSIM run histograms/information");
  // event data tree
  theEventOutputTree     = new TTree("Event","BDSIM event");
  
  // Build options and write structure
  theOptionsOutputTree->Branch("Options.",      "BDSOutputROOTEventOptions",optionsOutput,32000,2);
  
  // Build model and write structure
  theModelOutputTree->Branch("Model.",          "BDSOutputROOTEventModel",modelOutput,32000);

  // Build run data tree
  theRunOutputTree->Branch("Histos.",           "BDSOutputROOTEventHistograms",runHistos,32000,1);
  theRunOutputTree->Branch("Info.",             "BDSOutputROOTEventRunInfo",runInfo,32000,1);

  // Event info output
  theEventOutputTree->Branch("Info.",           "BDSOutputROOTEventInfo",evtInfo,32000,1);

  // Build primary structures
  if (WritePrimaries())
    {theEventOutputTree->Branch("Primary.",     "BDSOutputROOTEventSampler",primary,32000,1);}

  // Build loss and hit structures
  theEventOutputTree->Branch("Eloss.",          "BDSOutputROOTEventLoss",eLoss,4000,1);
  theEventOutputTree->Branch("PrimaryFirstHit.","BDSOutputROOTEventLoss",pFirstHit,4000,2);
  theEventOutputTree->Branch("PrimaryLastHit.", "BDSOutputROOTEventLoss",pLastHit, 4000,2);
  theEventOutputTree->Branch("TunnelHit.",      "BDSOutputROOTEventLoss",tHit, 4000,2);

  // Build trajectory structures
  theEventOutputTree->Branch("Trajectory.",     "BDSOutputROOTEventTrajectory",traj,4000,2);

  // Build event histograms
  theEventOutputTree->Branch("Histos.",         "BDSOutputROOTEventHistograms",evtHistos,32000,1);

  // build sampler structures
  for (G4int i = 0; i < (G4int)samplerTrees.size(); ++i)
    {
      auto samplerTreeLocal = samplerTrees.at(i);
      auto samplerName      = samplerNames.at(i);
      // set tree branches
      theEventOutputTree->Branch((samplerName+".").c_str(),
                                 "BDSOutputROOTEventSampler",
                                 samplerTreeLocal,32000,0);
    }
}

void BDSOutputROOT::WriteOptions()
{
  theOptionsOutputTree->Fill();
}

void BDSOutputROOT::WriteModel()
{
  theModelOutputTree->Fill();
}

void BDSOutputROOT::WriteFileEventLevel()
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ <<G4endl;
#endif
  theRootOutputFile->cd();
  theEventOutputTree->Fill();
}

void BDSOutputROOT::WriteFileRunLevel()
{
  theRootOutputFile->cd();
  theRunOutputTree->Fill();

  if(theRootOutputFile && theRootOutputFile->IsOpen())
    {theRootOutputFile->Write(nullptr,TObject::kOverwrite);}
}

void BDSOutputROOT::CloseFile()
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ <<G4endl;
#endif
  if(theRootOutputFile && theRootOutputFile->IsOpen())
  {
    theRootOutputFile->cd();
    theRootOutputFile->Write(0,TObject::kOverwrite);
    theRootOutputFile->Close();
    delete theRootOutputFile;
    theRootOutputFile = nullptr;
  }
}
