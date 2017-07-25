#ifndef BDSACCELERATORMODEL_H
#define BDSACCELERATORMODEL_H

#include "BDSBeamlineSet.hh"

#include "globals.hh"         // geant4 globals / types

#include <map>
#include <vector>

class BDSBeamline;
class BDSFieldObjects;
class G4LogicalVolume;
class G4ProductionCuts;
class G4Region;
class G4VPhysicalVolume;
class G4VSolid;

/**
 * @brief A holder class for all representations of the
 * accelerator model created in BDSIM. 
 *
 * This can be extended to 
 * allow inspection of the model. Holds the readout geometry
 * physical world in a location independent of detector 
 * construction.
 * 
 * In future, there may be several more beamlines than just
 * a flat one, and perhaps grouped into a more hierarchical
 * version.  These can be contained here and this class can
 * be extended as required.
 * 
 * @author Laurie Nevay
 */

class BDSAcceleratorModel
{
public:
  static BDSAcceleratorModel* Instance();
  ~BDSAcceleratorModel();

  /// @{ Register consituent of world.
  inline void RegisterWorldPV(G4VPhysicalVolume* worldIn) {worldPV = worldIn;}
  inline void RegisterWorldLV(G4LogicalVolume*   worldIn) {worldLV = worldIn;}
  inline void RegisterWorldSolid(G4VSolid*       worldIn) {worldSolid = worldIn;}
  /// @}
  
  /// Access the physical volume of the world
  inline G4VPhysicalVolume* GetWorldPV() const {return worldPV;}

  /// Register the main beam line set.
  inline void RegisterBeamlineSetMain(BDSBeamlineSet setIn) {mainBeamlineSet = setIn;}

  /// Register a set of beam lines to be managed and cleared up at the end of the simulation.
  void RegisterBeamlineSetExtra(G4String       name,
				BDSBeamlineSet setIn);
  
  /// @{ Accessor.
  inline const BDSBeamlineSet& BeamlineSetMain() const {return mainBeamlineSet;}
  const BDSBeamlineSet& BeamlineSet(G4String name) const;
  inline const std::map<G4String, BDSBeamlineSet>& ExtraBeamlines() const {return extraBeamlines;}
  const BDSBeamline* BeamlineMain() const {return mainBeamlineSet.massWorld;}
  /// @}

  /// Register the  beam line of arbitrary placements.
  inline void RegisterPlacementBeamline(BDSBeamline* placementBLIn) {placementBeamline = placementBLIn;}

  /// Access the beam line of arbitrary placements.
  inline BDSBeamline* PlacementBeamline() const {return placementBeamline;}
  
  /// Register the beam line containing all the tunnel segments
  inline void RegisterTunnelBeamline(BDSBeamline* beamlineIn) {tunnelBeamline = beamlineIn;}

  /// Access the beam line containing all the tunnel segments
  inline BDSBeamline* GetTunnelBeamline() const {return tunnelBeamline;}
  
  /// Register all field objects
  inline void RegisterFields(std::vector<BDSFieldObjects*>& fieldsIn){fields = fieldsIn;}

  /// Register a region and associated production cut as G4Region doesn't seem to delete
  /// it - note, no checking for double registration.
  void RegisterRegion(G4Region* region, G4ProductionCuts* cut);

  /// Register a temporary file for possible deletion at the deletion of the accelerator model
  /// based on the option from BDSGlobalConstants.
  void RegisterTemporaryFile(G4String fileName);

  /// Access region information. Will exit if not found.
  G4Region*         Region(G4String name) const;

  /// Simpler accessor for production cuts vs regions.
  G4ProductionCuts* ProductionCuts(G4String name) {return cuts.at(name);}

private:
  BDSAcceleratorModel(); ///< Default constructor is private as singleton.

  static BDSAcceleratorModel* instance;

  G4bool removeTemporaryFiles;

  G4VPhysicalVolume* worldPV;              ///< Physical volume of the mass world.
  G4LogicalVolume*   worldLV;
  G4VSolid*          worldSolid;

  BDSBeamlineSet mainBeamlineSet;
  std::map<G4String, BDSBeamlineSet> extraBeamlines; ///< Extra beamlines.

  BDSBeamline* tunnelBeamline;            ///< Tunnel segments beam line.
  BDSBeamline* placementBeamline;         ///< Placement beam line
  
  std::vector<BDSFieldObjects*> fields;       ///< All field objects.
  std::map<G4String, G4Region*> regions;      ///< All regions.
  std::map<G4String, G4ProductionCuts*> cuts; ///< Cuts corresponding to the regions.
  std::vector<G4String> temporaryFiles;       ///< All temporary file paths.
};

#endif
