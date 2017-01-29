#ifndef BDSBEAMPIPEFACTORY_H
#define BDSBEAMPIPEFACTORY_H

#include "BDSBeamPipeType.hh"

#include "globals.hh"           // geant4 globals / types
#include "G4ThreeVector.hh"

class G4Material;

class BDSBeamPipe;
class BDSBeamPipeInfo;
class BDSBeamPipeFactoryBase;

/**
 * @brief The main interface for using the beam pipe factories. 
 * 
 * Each function in this class mirrors one in BDSBeamPipeFactoryBase 
 * but with the addition of a BDSBeamPipeType enum as an argument. This
 * factory dispatches the call to the correct factory and should be 
 * the only one used by the caller.
 * 
 * @author Laurie Nevay
 */

class BDSBeamPipeFactory
{
public:
  static BDSBeamPipeFactory* Instance(); ///< Singleton accessor.
  
  ~BDSBeamPipeFactory();

  BDSBeamPipe* CreateBeamPipe(G4String         name,
			      G4double         length,
			      BDSBeamPipeInfo* bpi);
  
  BDSBeamPipe* CreateBeamPipe(BDSBeamPipeType beamPipeTypeIn,            // aperture type
			      G4String        nameIn,                    // name
			      G4double        lengthIn,                  // length [mm]
			      G4double        aper1 = 0,                 // aperture parameter 1
			      G4double        aper2 = 0,                 // aperture parameter 2
			      G4double        aper3 = 0,                 // aperture parameter 3
			      G4double        aper4 = 0,                 // aperture parameter 4
			      G4Material*     vacuumMaterialIn = nullptr,// vacuum material
			      G4double        beamPipeThicknessIn = 0,   // beampipe thickness [mm]
			      G4Material*     beamPipeMaterialIn = nullptr); // beampipe material

  BDSBeamPipe* CreateBeamPipe(BDSBeamPipeType beamPipeType,
			      G4String        name,
			      G4double        length,
			      G4ThreeVector   inputFaceNormal,
			      G4ThreeVector   outputFaceNormal,
			      G4double        aper1,
			      G4double        aper2,
			      G4double        aper3,
			      G4double        aper4,
			      G4Material*     vacuumMaterial,
			      G4double        beamPipeThickness,
			      G4Material*     beamPipeMaterial);

private:
  BDSBeamPipeFactory(); ///< Private constructor as singleton pattern.
  static BDSBeamPipeFactory* instance; ///< Singleton instance pointer.

  /// Return the appropriate factory singleton pointer given a type.
  BDSBeamPipeFactoryBase* GetAppropriateFactory(BDSBeamPipeType beamPipeTypeIn);

};



#endif
