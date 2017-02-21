#ifndef BDSFIELDMAGMUONSPOILER_H
#define BDSFIELDMAGMUONSPOILER_H

#include "BDSFieldMag.hh"

#include "globals.hh" // geant4 types / globals
#include "G4ThreeVector.hh"

class BDSMagnetStrength;

/** 
 * @brief Field of a Muon Spoiler.
 *
 */

class BDSFieldMagMuonSpoiler: public BDSFieldMag
{
public:
  BDSFieldMagMuonSpoiler(BDSMagnetStrength const* strength,
			 G4double          const  brho);

  virtual ~BDSFieldMagMuonSpoiler(){;}

  /// Access the field value.
  virtual G4ThreeVector GetField(const G4ThreeVector &position,
				 const G4double       t = 0) const;
  
private:
  G4double bField;
};

#endif
