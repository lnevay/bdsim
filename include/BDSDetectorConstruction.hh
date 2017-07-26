#ifndef BDSDETECTORCONSTRUCTION_H
#define BDSDETECTORCONSTRUCTION_H 

#include "BDSExtent.hh"

#include "globals.hh" // geant4 types / globals
#include "G4Transform3D.hh"
#include "G4Version.hh"
#include "G4VUserDetectorConstruction.hh"

#include <list>
#include <string>
#include <vector>

class G4LogicalVolume;
class G4Region;
class G4VPhysicalVolume;

namespace GMAD {
  struct Element;
  template<typename T> class FastList;
  class Placement;
}

class BDSAcceleratorModel;
class BDSBeamline;
class BDSBeamlineSet;
class BDSFieldObjects;
class BDSShowerModel;

#if G4VERSION_NUMBER > 1009
class BDSBOptrMultiParticleChangeCrossSection;
#endif

/**
 * @brief Class that constructs a Geant4 model of an accelerator.
 *
 * Mandatory class that must be supplied for a valid Geant4 simulation.
 * As construction is based on user input (via the parser) for a potentially
 * large number of items, the construction process is factorised into key steps
 * with member functions and makes use of factories for components and tunnel segments.
 * 
 * Maintained by Laurie Nevay & Jochem Snuverink
 */

class BDSDetectorConstruction: public G4VUserDetectorConstruction
{
public:
  BDSDetectorConstruction();
  virtual ~BDSDetectorConstruction();

  /// Overridden Geant4 method that must be implemented. Constructs the Geant4 geometry
  /// and returns the finished world physical volume.
  virtual G4VPhysicalVolume* Construct();

  virtual void ConstructSDandField();

  /// Create biasing operations.
  void BuildPhysicsBias();

  /// Public access to the world extent.
  BDSExtent WorldExtent() const {return worldExtent;}

  /// Loop over a beam line and place elements in a container (world). If a sensitive
  /// detector is specified, this is applied to each component. If regions are desired,
  /// the element is looked up in the region definitions and that is set up. If
  /// registerInfo, physical volume info is created and placed in a pv info registry.
  /// Public and static so it can be used by parallel world constructors. Last argument
  /// is whether to use the placement transform for curvilinear coordinate geometry that's
  /// different in the case of tilted dipoles.
  static void PlaceBeamlineInWorld(BDSBeamline*          beamline,
				   G4VPhysicalVolume*    containerPV,
				   G4bool                checkOverlaps     = false,
				   G4VSensitiveDetector* sensitiveDetector = nullptr,
				   G4bool                setRegions        = false,
				   G4bool                registerInfo      = false,
				   G4bool                useCLPlacementTransform = false);

  /// Create a transform based on the information in the placement.
  static G4Transform3D CreatePlacementTransform(const GMAD::Placement& placement);
  
private:
  /// assignment and copy constructor not implemented nor used
  BDSDetectorConstruction& operator=(const BDSDetectorConstruction&);
  BDSDetectorConstruction(BDSDetectorConstruction&);

  /// Create and set parameters for various G4Regions
  void InitialiseRegions();
  
  /// Build the main beam line and then any other required beam lines.
  void BuildBeamlines();

  /// Convert a parser beamline_list to BDSAcceleratorComponents with help of
  /// BDSComponentFactory and put in a BDSBeamline container that calculates coordinates
  /// and extents of the beamline.
  BDSBeamlineSet BuildBeamline(const GMAD::FastList<GMAD::Element>& beamLine,
			       G4String             name,
			       const G4Transform3D& initialTransform   = G4Transform3D(),
			       G4bool               beamlineIsCircular = false);

  /// Build the tunnel around the already constructed flat beam line.
  void BuildTunnel();
  
  /// Build the world volume using the extent of the BDSBeamline instance created
  /// in BuildBeamline()
  G4VPhysicalVolume* BuildWorld();
  
  /// Place beam line, tunnel beam line, end pieces and placements in world.
  void ComponentPlacement(G4VPhysicalVolume* worldPV);

  /// Detect whether the first element has an angled face such that it might overlap
  /// with a previous element.  Only used in case of a circular machine.
  G4bool UnsuitableFirstElement(std::list<GMAD::Element>::const_iterator element);

#if G4VERSION_NUMBER > 1009
  /// Function that creates physics biasing cross section
  BDSBOptrMultiParticleChangeCrossSection* BuildCrossSectionBias(const std::list<std::string>& biasList,
								 G4String defaultBias,
								 G4String elementName);

  /// List of bias objects - for memory management
  std::vector<BDSBOptrMultiParticleChangeCrossSection*> biasObjects;
#endif

#ifdef BDSDEBUG
  bool debug = true;
#else
  bool debug = false;
#endif

  ///@{ Variable copied from global constants
  G4bool verbose;
  G4bool checkOverlaps;
  ///@}

  /// Accelerator model pointer
  BDSAcceleratorModel* acceleratorModel;

  /// All fields
  std::vector<BDSFieldObjects*> fields;

  G4bool       circular;    ///< Whether or not we're building a circular machine.
  BDSExtent    worldExtent; ///< Record of the world extent.
  BDSBeamline* placementBL; ///< Placement beam line.
};

#endif

