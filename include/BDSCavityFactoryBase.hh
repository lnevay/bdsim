#ifndef BDSCAVITYFACTORYBASE_H
#define BDSCAVITYFACTORYBASE_H

#include "globals.hh"
#include "G4String.hh"
#include "G4RotationMatrix.hh"

#include <vector>

class BDSCavity;
class BDSCavityInfo;
class BDSExtent;

class G4LogicalVolume;
class G4Material;
class G4UserLimits;
class G4VisAttributes;
class G4VPhysicalVolume;
class G4VSolid;

/**
 * @brief Abstract base class for rf cavity geometry factories.
 *
 * @author Laurie Nevay
 */

class BDSCavityFactoryBase
{
public:
  virtual ~BDSCavityFactoryBase(){;}

  /// Method that calls protected virtual methods in order. Sets general
  /// recipe of construction - in order below.
  BDSCavity* CreateCavity(G4String             name,
			  G4double             totalChordLength,
			  const BDSCavityInfo* info,
			  G4Material*          vacuumMaterial);
  
protected:
  BDSCavityFactoryBase();

  /// Ensure all internal members are reset before reuse.
  virtual void CleanUp();

  /// Create vacuumSolid and cavitySolid. Must return the container radius.
  virtual G4double CreateSolids(G4String             name,
				G4double             totalChordLength,
				const BDSCavityInfo* info) = 0;

  /// Create vacuumLV, cavityLV and containerLV.
  virtual void CreateLogicalVolumes(G4String             name,
				    const BDSCavityInfo* info,
				    G4Material*          vacuumMaterial);

  /// If the length is shorter than the default user limits in global constants,
  /// create a new user limits object and apply to all lvs in allLogicalVolumes.
  virtual void SetUserLimits(G4double length);

  /// Set vis attributes for cavityLV, vacuumLV and containerLV.
  virtual void SetVisAttributes(G4String colourName = "rfcavity");

  /// Place the vacuum and cavity centrally in the container logical volume.
  virtual void PlaceComponents(G4String name);

  /// Construct output object and register all temporary objects from vectors.
  virtual BDSCavity* BuildCavityAndRegisterObjects(const BDSExtent& extent);
  
  /// @{ Holder for registration at end of construction.
  std::vector<G4VSolid*>          allSolids;
  std::vector<G4LogicalVolume*>   allLogicalVolumes;
  std::vector<G4LogicalVolume*>   allSensitiveVolumes;
  std::vector<G4VPhysicalVolume*> allPhysicalVolumes;
  std::vector<G4RotationMatrix*>  allRotationMatrices;
  std::vector<G4UserLimits*>      allUserLimits;
  std::vector<G4VisAttributes*>   allVisAttributes;
  /// @}

  /// @{ Cache of particular solid or lv for common functionality in this class.
  G4VSolid*        vacuumSolid;
  G4VSolid*        cavitySolid; 
  G4VSolid*        containerSolid;
  G4LogicalVolume* vacuumLV;
  G4LogicalVolume* cavityLV;
  G4LogicalVolume* containerLV;
  /// @}

  /// @{ Cache of variable from BDSGlobalConstants.
  G4double    lengthSafety;  
  G4double    nSegmentsPerCircle;
  G4Material* emptyMaterial;
  G4bool      checkOverlaps;
  /// @}
};

#endif
