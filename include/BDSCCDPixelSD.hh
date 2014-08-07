#ifndef BDSCCDPixelSD_h
#define BDSCCDPixelSD_h 1

#include "G4VSensitiveDetector.hh"
#include "BDSCCDPixelHit.hh"
#include "G4TransportationManager.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class BDSCCDPixelSD : public G4VSensitiveDetector
{
  
public:
  BDSCCDPixelSD(G4String name);
  ~BDSCCDPixelSD();
  
  void Initialize(G4HCofThisEvent*HCE);
  G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
  void EndOfEvent(G4HCofThisEvent*HCE);
  void clear();
  void DrawAll();
  void PrintAll();
  
private:
  BDSCCDPixelHitsCollection *CCDPixelCollection;
  G4bool                    StoreHit;
  G4String itsType;
  G4String _collectionName;
};
#endif

