#include "BDSScreen.hh"

#include "BDSMultilayerScreen.hh"
#include "BDSScreenLayer.hh"

#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <list>
#include <sstream>

BDSScreen::BDSScreen(G4String         nameIn,  
		     G4double         chordLengthIn,
		     BDSBeamPipeInfo* beamPipeInfoIn,
		     G4TwoVector      sizeIn, //X Y dimensions of screen
		     G4double         screenAngleIn):
  BDSDrift(nameIn,chordLengthIn,beamPipeInfoIn),
  size(sizeIn), 
  screenAngle(screenAngleIn),
  screenPos(G4ThreeVector()),
  nLayers(0)
{
  mlScreen = new BDSMultilayerScreen(size, nameIn+"_mlscreen");
  
  screenRot = new G4RotationMatrix();
  screenRot->rotateY(screenAngle);
}

BDSScreen::~BDSScreen()
{
  delete screenRot;
  delete mlScreen;
}

void BDSScreen::Build()
{
  BDSDrift::Build(); //Build the beam pipe geometry.
  
  // Make the beam pipe wireframe
  G4VisAttributes* VisAtt1 = new G4VisAttributes(G4Colour(0.4, 0.4, 0.4, 0.3));
  VisAtt1->SetForceWireframe(true);
  VisAtt1->SetVisibility(true);
  containerLogicalVolume->SetVisAttributes(VisAtt1);

  PlaceScreen(); //Place the screen in the beam pipe
}

void BDSScreen::AddScreenLayer(G4double thickness, G4String material, G4int isSampler)
{
  std::stringstream ss;
  ss << nLayers;
  G4String lNum = ss.str();
  G4String lName = name+"_"+lNum;
  mlScreen->AddScreenLayer(thickness,material,lName, isSampler);
  if(!isSampler)
    {RegisterSensitiveVolume(mlScreen->LastLayer()->GetLog());}
  nLayers++;
}

void BDSScreen::PlaceScreen()
{
  mlScreen->Build();//Build the screen.
  G4LogicalVolume* vacuumLV = GetAcceleratorVacuumLogicalVolume();
  mlScreen->Place(screenRot, screenPos, vacuumLV); //Place the screen in the beampipe centre.
}
