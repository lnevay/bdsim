#ifndef BDSBEAMPIPEFACTORYLHC_H
#define BDSBEAMPIPEFACTORYLHC_H

#include "BDSBeamPipeFactoryBase.hh"
#include "BDSBeamPipe.hh"

/**
 * @brief factory for simple lhc aperture model beam pipes
 * 
 * singleton pattern
 * 
 * @author Laurie Nevay <laurie.nevay@rhul.ac.uk>
 */

class BDSBeamPipeFactoryLHC: public BDSBeamPipeFactoryBase
{
public:
  static BDSBeamPipeFactoryLHC* Instance(); ///< singleton accessor
  
  virtual ~BDSBeamPipeFactoryLHC();

  virtual BDSBeamPipe* CreateBeamPipe(G4String    nameIn,                     // name
				      G4double    lengthIn,                   // length [mm]
				      G4double    aper1 = 0,                  // aperture parameter 1
				      G4double    aper2 = 0,                  // aperture parameter 2
				      G4double    aper3 = 0,                  // aperture parameter 3
				      G4double    aper4 = 0,                  // aperture parameter 4
				      G4Material* vacuumMaterialIn = NULL,    // vacuum material
				      G4double    beamPipeThicknessIn = 0,    // beampipe thickness [mm]
				      G4Material* beamPipeMaterialIn = NULL   // beampipe material
				      );
  
  virtual BDSBeamPipe* CreateBeamPipeAngledInOut(G4String    nameIn,
						 G4double    lengthIn,
						 G4double    angleInIn,  // the normal angle of the input face
						 G4double    angleOutIn, // the normal angle of the output face
						 G4double    aper1 = 0,
						 G4double    aper2 = 0,
						 G4double    aper3 = 0,
						 G4double    aper4 = 0,
						 G4Material* vacuumMaterialIn = NULL,
						 G4double    beamPipeThicknessIn = 0,
						 G4Material* beamPipeMaterialIn = NULL
						 );

private:
  BDSBeamPipeFactoryLHC(); ///< private default constructor - singleton pattern
  static BDSBeamPipeFactoryLHC* _instance;

  void TestInputParameters(G4Material*& vacuumMaterialIn,
			   G4double&    beamPipeThicknessIn,
			   G4Material*& beamPipeMaterialIn,
			   G4double&    aper1In,
			   G4double&    aper2In,
			   G4double&    aper3In);

  //abstract common build features to one function
  //use member variables unique to this factory to pass them around

  BDSBeamPipe* CommonFinalConstruction(G4String    nameIn,
				       G4Material* vacuumMaterialIn,
				       G4Material* beamPipeMaterialIn,
				       G4double    lengthIn,
				       G4double    widthIn,
				       G4double    heightIn);
  void CreateGeneralAngledSolids(G4String      nameIn,
				 G4double      lengthIn,
				 G4double      aper1In,
				 G4double      aper2In,
				 G4double      aper3In,
				 G4double      beamPipeThicknessIn,
				 G4ThreeVector inputfaceIn,
				 G4ThreeVector outputfaceIn);

  /// function to create the container subtraction solid and avoid passing aper1,2,3 around further
  void CreateContainerSubtractionSolid(G4String& nameIn,
				       G4double& lengthIn,
				       G4double& beamPipeThicknessIn,
				       G4double& aper1In,
				       G4double& aper2In,
				       G4double& aper3In);

};
  
#endif
