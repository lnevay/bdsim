#ifndef BDSPHYSICALVOLUMEINFO_H
#define BDSPHYSICALVOLUMEINFO_H

#include "G4String.hh"
#include "globals.hh" // geant4 types / globals

/**
 * @brief A class holding any information pertaining to a particular physical
 * volume in a BDSIM geant4 model.  
 * 
 * Currently, this is used for the name and
 * s position coordinate, but can be extended to store any other pertinent 
 * information.
 * 
 * @author Laurie Nevay
 */

class BDSPhysicalVolumeInfo
{
public:
  explicit BDSPhysicalVolumeInfo(G4double SPosAtMiddleOfElement);
  BDSPhysicalVolumeInfo(G4String nameIn,
			G4String placementNameIn,
			G4double sPosIn,
			G4int    beamlineIndexIn = -1);
  ~BDSPhysicalVolumeInfo();

  /// output stream
  friend std::ostream& operator<< (std::ostream &out, BDSPhysicalVolumeInfo const &info);

  /// Get the name of the logical volume
  inline G4String GetName() const {return name;}

  /// Set the name of the logical volume
  inline void SetName(G4String NameIn) {name = NameIn;}

  /// Get the name of the logical volume
  inline G4String GetPlacementName() const {return placementName;}

  /// Set the name of the logical volume
  inline void SetPlacementName(G4String placementNameIn)
  {placementName = placementNameIn;}

  /// Get the s position coordinate of the logical volume
  G4double GetSPos() const {return spos;}

  /// Set the s position coordinate of the logical volume
  void SetSPos(G4double SPosIn) {spos = SPosIn;}

  /// Get beamline index
  G4int  GetBeamlineIndex() const {return beamlineIndex;}
  
private:
  BDSPhysicalVolumeInfo();

  G4String name;
  G4String placementName;
  G4double spos;
  G4int    beamlineIndex;
};


#endif
