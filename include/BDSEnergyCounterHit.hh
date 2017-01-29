#ifndef BDSENERGYCOUNTERHIT_H
#define BDSENERGYCOUNTERHIT_H

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

class BDSEnergyCounterHit;

typedef G4THitsCollection<BDSEnergyCounterHit> BDSEnergyCounterHitsCollection;
extern G4Allocator<BDSEnergyCounterHit> BDSEnergyCounterHitAllocator;

// namespace BDS
// {
//   BDSEnergyCounterHit* LowestSPosPrimaryHit (BDSEnergyCounterHitsCollection* HC);
//   BDSEnergyCounterHit* HighestSPosPrimaryHit(BDSEnergyCounterHitsCollection* HC);
// }

/**
 * @brief Information recorded for a single piece of energy deposition.
 *
 */

class BDSEnergyCounterHit: public G4VHit
{
public:
  /// Default (in effect) constructor for energy counter hit. The intention (by a lack of
  /// setter methods is that all information should be provided as instantiation time for
  /// an instance of this class.
  BDSEnergyCounterHit(G4int    nCopyIn           = 0,    // copy number of physical volume
		      G4double energyIn          = 0,    // energy in this 'hit'
		      G4double XIn               = 0,    // global x
		      G4double YIn               = 0,    // global y
		      G4double ZIn               = 0,    // global x
		      G4double SBeforeIn         = 0,    // s of pre step coordinate
		      G4double SAfterIn          = 0,    // s of post step coordinate
		      G4double SHitIn            = 0,    // assigned s of loss
		      G4double xIn               = 0,    // local x
		      G4double yIn               = 0,    // local y
		      G4double zIn               = 0,    // local z
		      G4String nameIn            = "",   // volume name
		      G4int    partIDIn          = 0,    // PDG id - particle type
		      G4int    trackID           = -1,   // Track ID
		      G4int    parentID          = -1,   // Parent ID
		      G4double weightIn          = 1,    // weight
		      G4int    turnsTakenIn      = 1,    // turns taken if circular
		      G4int    eventNoIn         = 0,    // event number
		      G4double stepLengthIn      = 0,
		      G4int    beamlineIndexIn   = -1,
		      G4int    geomFlag          = -1);
  
  virtual ~BDSEnergyCounterHit();
  
  inline void* operator new(size_t) ;
  inline void operator delete(void *aHit);

  inline G4int    GetCopyNumber()      const {return copyNumber;}
  inline G4double GetEnergy()          const {return energy;}
  inline void     SetEnergy(G4double energyIn) {energy = energyIn;}
  inline G4double GetX()               const {return X;} 
  inline G4double GetY()               const {return Y;}
  inline G4double GetZ()               const {return Z;}
  inline G4double GetSBefore()         const {return sBefore;}
  inline G4double GetSAfter()          const {return sAfter;}
  inline G4double GetSHit()            const {return sHit;}
  inline G4double Getx()               const {return x;} 
  inline G4double Gety()               const {return y;}
  inline G4double Getz()               const {return z;} 
  inline G4String GetName()            const {return name;}
  inline G4int    GetPartID()          const {return partID;}
  inline G4int    GetTrackID()         const {return trackID;}
  inline G4int    GetParentID()        const {return parentID;}
  inline G4double GetWeight()          const {return weight;} 
  inline G4int    GetTurnsTaken()      const {return turnsTaken;}
  inline G4int    GetEventNo()         const {return eventNo;}
  inline G4double GetStepLength()      const {return stepLength;}
  inline G4int    GetBeamlineIndex()   const {return beamlineIndex;}
  inline G4int    GetGeomFlag()        const {return geomFlag;}
  inline G4double GetEnergyWeighted()  const {return weight * energy;}
  
private:
  /// Private default constructor (not implemented) as the constructor
  /// provided as defaults for everything which in effect provides a default
  /// constructor
  BDSEnergyCounterHit();

  G4int    copyNumber;
  G4double energy;
  /// @{ Global coordinate
  G4double X;
  G4double Y;
  G4double Z;
  /// @}
  
  G4double sBefore; ///< s coordinate of pre step point
  G4double sAfter;  ///< s coordinate of post step point

  /// The position this energy deposition is assigned to. As we represent
  /// a Monte Carlo of discrete processes, this is randomly scattered along
  /// the length of the track this energy deposition corresponds to (uniform
  /// random distribution).
  G4double sHit;
  
  /// @{ Local coordinate
  G4double x;
  G4double y;
  G4double z;
  /// @}
  
  G4String name;
  G4int    partID;
  G4int    trackID;
  G4int    parentID;
  G4double weight;
  G4int    turnsTaken;
  G4int    eventNo;
  G4double stepLength;
  G4int    beamlineIndex;
  G4int    geomFlag;
};

inline void* BDSEnergyCounterHit::operator new(size_t)
{
  void* aHit;
  aHit=(void*) BDSEnergyCounterHitAllocator.MallocSingle();
  return aHit;
}

inline void BDSEnergyCounterHit::operator delete(void *aHit)
{
 BDSEnergyCounterHitAllocator.FreeSingle((BDSEnergyCounterHit*) aHit);
}

#endif

///
