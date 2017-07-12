#ifdef USE_GDML

#ifndef BDSGEOMETRYFACTORYGDML_H
#define BDSGEOMETRYFACTORYGDML_H

#include "BDSGeometryFactoryBase.hh"

#include "globals.hh"

#include <map>
#include <vector>

class BDSGeometryExternal;
class G4Colour;
class G4LogicalVolume;
class G4VPhysicalVolume;

/**
 * @brief Interface to Geant4's GDML loader.
 *
 * @author Stewart Boogert, Lawrence Deacon, Laurie Nevay
 */

class BDSGeometryFactoryGDML: public BDSGeometryFactoryBase
{
public:
  /// Singleton accessor.
  static BDSGeometryFactoryGDML* Instance();

  virtual ~BDSGeometryFactoryGDML();

  /// Main interface overrides base class method to construct a file with
  /// optional colour mapping.
  virtual BDSGeometryExternal* Build(G4String componentName,
				     G4String fileName,
				     std::map<G4String, G4Colour*>* colourMapping = nullptr,
				     G4double suggestedLength        = 0,
				     G4double suggestedOuterDiameter = 0);

private:
  /// Private default constructor as singleton.
  BDSGeometryFactoryGDML();

  /// Create a temporary file in the current working directory (even though the geometry file
  /// may not exist there) and create a copy of the input geometry file line by line, but replacing
  /// a 'key' with 'replacement'. Returns the temporary file name created. Naming only allows one
  /// temporary file per component - so far, not a problem.
  G4String ReplaceStringInFile(G4String componentName,
			       G4String filename,
			       G4String key,
			       G4String replacement);
  
  /// Singleton instance
  static BDSGeometryFactoryGDML* instance;

  /// Recursively append pvs and lvs from a given mother volume.  Pass by reference
  /// the output vectors
  void GetAllLogicalAndPhysical(const G4VPhysicalVolume*         volume,
				std::vector<G4VPhysicalVolume*>& pvs,
				std::vector<G4LogicalVolume*>&   lvs);
};

#endif

#endif
