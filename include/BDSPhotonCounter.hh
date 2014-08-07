/** BDSIM

L. Deacon
A class for counting the number of photons generated by scintillators.

**/

//==============================================================

#ifndef BDSPhotonCounter_h
#define BDSPhotonCounter_h 

#include "G4Track.hh"

class BDSPhotonCounter 
{
public:
  static BDSPhotonCounter* Instance();
  ~BDSPhotonCounter();
  
  void countPhoton(const G4Track* aTrack);
  G4double nPhotons() const;
  G4double energy() const;

private:
  BDSPhotonCounter();
  G4double _nPhotons;
  G4double _energy;
  static BDSPhotonCounter* _instance;

};

inline G4double BDSPhotonCounter::nPhotons() const{
  return _nPhotons;
}

inline G4double BDSPhotonCounter::energy() const{
  return _energy;
}

#endif

