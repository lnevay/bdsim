#ifndef BDSOutputNone_h
#define BDSOutputNone_h 

#include "BDSOutputBase.hh"

/**
 * @brief No output class
 * 
 * Write no output
 *
 * @author Jochem Snuverink <Jochem.Snuverink@rhul.ac.uk>
 */

class BDSOutputNone : public BDSOutputBase {

public: 

  BDSOutputNone():BDSOutputBase(){} ///< default constructor
  virtual ~BDSOutputNone(){}

  /// write sampler hit collection
  virtual void WriteHits(BDSSamplerHitsCollection*){}
  /// make energy loss histo
  virtual void WriteEnergyLoss(BDSEnergyCounterHitsCollection*){}
  /// write primary loss histo
  virtual void WritePrimaryLoss(BDSEnergyCounterHit*){}
  /// write primary hits histo
  virtual void WritePrimaryHit(BDSEnergyCounterHit*){}
  /// write tunnel hits
  virtual void WriteTunnelHits(BDSTunnelHitsCollection*){}
  /// write a trajectory
  virtual void WriteTrajectory(std::vector<BDSTrajectory*> &){}
  /// write primary hit
  virtual void WritePrimary(G4double,G4double,G4double,G4double,G4double,G4double,G4double,G4double,G4double,G4int,G4int,G4int){}
  /// write a histogram
  virtual void WriteHistogram(BDSHistogram1D*){}
  /// write event info
  virtual void WriteEventInfo(time_t,time_t,G4float){}

  virtual void WriteTrackerBunch(G4String, TRKBunch*, bool){}
  virtual void WriteTrackerPrimaryLoss(TRKBunch*, int){}
  
  virtual void FillEvent(){}///< fill event structure
  virtual void Initialise(){}; ///< open the file
  virtual void Write(){};      ///< write to file
  virtual void Close(){};      ///< close the file
};

extern BDSOutputBase* bdsOutput;
#endif
