#include "BDSCollimatorRectangular.hh"
#include "BDSCollimatorBase.hh"

#include "globals.hh" // geant4 globals / types
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4VSolid.hh"

BDSCollimatorRectangular::BDSCollimatorRectangular(G4String nameIn,
						   G4double lengthIn,
						   G4double outerDiameterIn,
						   G4double xApertureIn,
						   G4double yApertureIn,
						   G4double xOutApertureIn,
						   G4double yOutApertureIn,
						   G4String collimatorMaterialIn,
						   G4String vacuumMaterialIn,
						   G4String colourIn):
  BDSCollimatorBase(nameIn, lengthIn, outerDiameterIn, "rcol",
		    xApertureIn, yApertureIn,xOutApertureIn, yOutApertureIn,
		    collimatorMaterialIn,
		    vacuumMaterialIn, colourIn)
{;}

void BDSCollimatorRectangular::BuildInnerCollimator()
{
  if(tapered)
    {
      // Make subtracted volume longer than the solid volume
      G4double xGradient = std::abs((xAperture - xOutAperture)) / chordLength;
      G4double yGradient = std::abs((yAperture - yOutAperture)) / chordLength;
      
      G4double deltam = 0.1 * chordLength;
      G4double deltax = xGradient * deltam;
      G4double deltay = yGradient * deltam;
      
      innerSolid  = new G4Trd(name + "_inner_solid",             // name
                              xAperture + deltax,                // X entrance half length
                              xOutAperture - deltax,             // X exit half length
                              yAperture + deltay,                // Y entrance half length
                              yOutAperture - deltay,             // Y exit half length
                              (chordLength + 2*deltam) * 0.5);   // Z half length
    
      vacuumSolid = new G4Trd(name + "_vacuum_solid",               // name
                              xAperture - lengthSafetyLarge,        // X entrance half length
                              xOutAperture - lengthSafetyLarge,     // X exit half length
                              yAperture - lengthSafetyLarge,        // Y entrance half length
                              yOutAperture - lengthSafetyLarge,     // Y exit half length
                              chordLength*0.5 - lengthSafetyLarge); // Z half length
    }
  else
    {
      innerSolid  = new G4Box(name + "_inner_solid",    // name
                              xAperture,                // x half width
                              yAperture,                // y half width
                              chordLength);             // z half length
    // z half length long for unambiguous subtraction

      vacuumSolid = new G4Box(name + "_vacuum_solid",   // name
                              xAperture - lengthSafety, // x half width
                              yAperture - lengthSafety, // y half width
                              chordLength*0.5);         // z half length
    }
    
  RegisterSolid(innerSolid);
  RegisterSolid(vacuumSolid);
}
