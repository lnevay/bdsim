/* BDSIM code.    Version 1.0
   Author: Grahame A. Blair, Royal Holloway, Univ. of London.
   Last modified 24.7.2002
   Copyright (c) 2002 by G.A.Blair.  ALL RIGHTS RESERVED. 
*/

#ifndef BDSEnergyCounterSD_h
#define BDSEnergyCounterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "BDSEnergyCounterHit.hh"
#include "G4Navigator.hh"
#include "G4GFlashSpot.hh"
#include "G4VGFlashSensitiveDetector.hh"

class G4VProcess;
class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class BDSEnergyCounterSD : public G4VSensitiveDetector, public G4VGFlashSensitiveDetector
{

public:
  BDSEnergyCounterSD(G4String name);
  ~BDSEnergyCounterSD();

  virtual void Initialize(G4HCofThisEvent*HCE);
  virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
  virtual G4bool ProcessHits(G4GFlashSpot*aSpot ,G4TouchableHistory* ROhist);

  G4double GetSPositionOfStep(G4Step* aStep,       G4TouchableHistory* readOutTH);
  G4double GetSPositionOfSpot(G4GFlashSpot* aSpot, G4TouchableHistory* readOutTH);
  G4String GetName();

private:
  /// assignment and copy constructor not implemented nor used
  BDSEnergyCounterSD& operator=(const BDSEnergyCounterSD&);
  BDSEnergyCounterSD(BDSEnergyCounterSD&);

  G4bool   verbose;
  G4String itsName;
  BDSEnergyCounterHitsCollection* energyCounterCollection;
  BDSEnergyCounterHitsCollection* primaryCounterCollection;
  G4int    HCIDe;
  G4int    HCIDp;
  G4double enrg;
  G4double X,Y,Z,S; // global coordinates
  G4double x,y,z;   // local coordinates
};

inline G4String BDSEnergyCounterSD::GetName()
{return itsName;}


#endif

