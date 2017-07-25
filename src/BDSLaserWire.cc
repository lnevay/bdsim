#include "BDSBeamPipeInfo.hh"
#include "BDSGlobalConstants.hh" 
#include "BDSLaserWire.hh"
#include "BDSMaterials.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"

BDSLaserWire::BDSLaserWire(G4String      nameIn,
			   G4double      lengthIn,
			   G4double      wavelengthIn,
			   G4ThreeVector directionIn):
  BDSAcceleratorComponent(nameIn, lengthIn, 0, "laserwire"),
  itsLaserDirection(directionIn),
  itsLaserWavelength(wavelengthIn)
{;}

BDSLaserWire::~BDSLaserWire()
{;}

void BDSLaserWire::BuildContainerLogicalVolume()
{
  G4double beamPipeRadius = BDSGlobalConstants::Instance()->GetDefaultBeamPipeModel()->aper1;
  containerSolid = new G4Box(name +"_container_solid",
			     beamPipeRadius,
			     beamPipeRadius,
			     chordLength*0.5);
  containerLogicalVolume=new G4LogicalVolume(containerSolid,
					     BDSMaterials::Instance()->GetMaterial("LaserVac"),
					     name + "_container_lv");
  BDSGlobalConstants::Instance()->SetLaserwireWavelength(containerLogicalVolume->GetName(),itsLaserWavelength);
  BDSGlobalConstants::Instance()->SetLaserwireDir(containerLogicalVolume->GetName(),itsLaserDirection);

  SetAcceleratorVacuumLogicalVolume(containerLogicalVolume);
}
