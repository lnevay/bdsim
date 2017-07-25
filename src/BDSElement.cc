#include "BDSDebug.hh"
#include "BDSExtent.hh"
#include "BDSElement.hh"
#include "BDSFieldBuilder.hh"
#include "BDSFieldFactory.hh"
#include "BDSGeometryExternal.hh"
#include "BDSGeometryFactory.hh"

#include "globals.hh" // geant4 globals / types

class BDSFieldInfo;

BDSElement::BDSElement(G4String      nameIn,
		       G4double      lengthIn,
		       G4double      outerDiameterIn,
		       G4String      geometryIn,
		       G4String      fieldNameIn):
  BDSAcceleratorComponent(nameIn, lengthIn, 0, "element"),
  outerDiameter(outerDiameterIn),
  geometryFileName(geometryIn),
  fieldName(fieldNameIn)
{;}

void BDSElement::BuildContainerLogicalVolume()
{
#ifdef BDSDEBUG 
  G4cout << __METHOD_NAME__ <<G4endl;
#endif

  // The outerDiameter here is a suggested outerDiameter for the factory. Each subfactory may treat this
  // differently.
  BDSGeometryExternal* geom = BDSGeometryFactory::Instance()->BuildGeometry(name, geometryFileName, nullptr,
									    chordLength, outerDiameter);
  
  if (!geom)
    {
      G4cerr << __METHOD_NAME__ << "Error loading geometry in component \""
	     << name << "\"" << G4endl;
      exit(1);
    }
  
  // We don't register the geometry as a daughter as the geometry factory retains
  // ownership of the geometry and will clean it up at the end.
  
  // make the beam pipe container, this object's container
  containerLogicalVolume = geom->GetContainerLogicalVolume();
  containerSolid         = geom->GetContainerSolid();

  // set placement offset from geom so it's placed correctly in the beam line
  SetPlacementOffset(geom->GetPlacementOffset());
  
  // update extents
  InheritExtents(geom);

  const BDSExtent geomExtent = geom->GetExtent();
  BDSExtent nominalExt = BDSExtent(outerDiameter*0.5, outerDiameter*0.5, chordLength*0.5);
  if (nominalExt.TransverselyGreaterThan(geomExtent))
    {SetExtent(nominalExt);}

  G4double extLength = GetExtent().DZ();
  if (extLength > chordLength)
    {
      G4cerr << "BDSElement> The loaded geometry is larger than the specified length"
	     << " of the element, which will cause overlaps!" << G4endl
	     << "Calculated extent along z: " << extLength << " mm, vs specified "
	     << chordLength << G4endl;
      exit(1);
    }

  // Get the field definition from the parser
  // Note, the field factory manages the deletion of this info instance.
  SetField(BDSFieldFactory::Instance()->GetDefinition(fieldName));
}
