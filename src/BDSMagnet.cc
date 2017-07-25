#include "BDSBeamPipe.hh"
#include "BDSBeamPipeFactory.hh"
#include "BDSBeamPipeInfo.hh"
#include "BDSBeamPipeType.hh"
#include "BDSDebug.hh"
#include "BDSExecOptions.hh"
#include "BDSFieldBuilder.hh"
#include "BDSFieldInfo.hh"
#include "BDSGlobalConstants.hh"
#include "BDSMaterials.hh"
#include "BDSMagnetGeometryType.hh"
#include "BDSMagnetOuter.hh"
#include "BDSMagnetOuterInfo.hh"
#include "BDSMagnetOuterFactory.hh"
#include "BDSMagnetStrength.hh"
#include "BDSMagnetType.hh"
#include "BDSMagnet.hh"
#include "BDSUtilities.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4VPhysicalVolume.hh"

#include <cstdlib>
#include <cmath>
#include <string>

class G4Userlimits;


BDSMagnet::BDSMagnet(BDSMagnetType       typeIn,
		     G4String            nameIn,
		     G4double            lengthIn,
		     BDSBeamPipeInfo*    beamPipeInfoIn,
		     BDSMagnetOuterInfo* magnetOuterInfoIn,
		     BDSFieldInfo*       vacuumFieldInfoIn,
		     G4double            angleIn,
		     BDSFieldInfo*       outerFieldInfoIn):
  BDSAcceleratorComponent(nameIn, lengthIn, angleIn, typeIn.ToString(), beamPipeInfoIn),
  magnetType(typeIn),
  magnetOuterInfo(magnetOuterInfoIn),
  vacuumFieldInfo(vacuumFieldInfoIn),
  outerFieldInfo(outerFieldInfoIn),
  beampipe(nullptr),
  placeBeamPipe(false),
  magnetOuterOffset(G4ThreeVector(0,0,0)),
  outer(nullptr),
  beamPipePlacementTransform(G4Transform3D())
{
  outerDiameter   = magnetOuterInfoIn->outerDiameter;
  containerRadius = 0.5*outerDiameter;
  inputface       = G4ThreeVector(0,0,-1);
  outputface      = G4ThreeVector(0,0, 1);
  
  beampipe = nullptr;
  outer    = nullptr;

  placeBeamPipe = false;

  // It's not possible to build advanced outer geometry for a very thin magnet.
  if (lengthIn < 1e-4*CLHEP::m) // 100um minimum length for geometry
    {magnetOuterInfo->geometryType = BDSMagnetGeometryType::none;}
  // No beam pipe geometry for really short 'magnets'
  if (lengthIn < 1e-6*CLHEP::m)
    {GetBeamPipeInfo()->beamPipeType = BDSBeamPipeType::circularvacuum;}
}

G4String BDSMagnet::DetermineScalingKey(BDSMagnetType typeIn)
{
  G4String result = "none";
  switch (typeIn.underlying())
    {
    case BDSMagnetType::hkicker:
    case BDSMagnetType::vkicker:
    case BDSMagnetType::muonspoiler:
    case BDSMagnetType::sectorbend:
    case BDSMagnetType::rectangularbend:
      {result = "field"; break;}
    case BDSMagnetType::quadrupole:
      {result = "k1"; break;}
    case BDSMagnetType::sextupole:
      {result = "k2"; break;}
    case BDSMagnetType::octupole:
      {result = "k3"; break;}
    case BDSMagnetType::decapole:
      {result = "k5"; break;}
    default:
      {break;} // leave as none without complaint
    };

  return result;
}

void BDSMagnet::Build()
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif  
  BuildBeampipe();
  BuildVacuumField();
  BuildOuter();
  // Instead of BDSAcceleratorComponent::Build just call BuildContainerLogicalVolume
  // to control user limits in the case where there is no container and we just inherit
  // the beam pipe container
  BuildContainerLogicalVolume();
  BuildOuterField(); // must be done when the containerLV exists
  PlaceComponents(); // place things (if needed) in container
}

void BDSMagnet::BuildBeampipe()
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif
  beampipe = BDSBeamPipeFactory::Instance()->CreateBeamPipe(name,
							    chordLength - lengthSafety,
							    beamPipeInfo);

  beamPipePlacementTransform = beampipe->GetPlacementTransform().inverse();
  
  RegisterDaughter(beampipe);
  InheritExtents(beampipe);

  SetAcceleratorVacuumLogicalVolume(beampipe->GetVacuumLogicalVolume());

  /// Update record of normal vectors now beam pipe has been constructed.
  SetInputFaceNormal(BDS::RotateToReferenceFrame(beampipe->InputFaceNormal(), angle));
  SetOutputFaceNormal(BDS::RotateToReferenceFrame(beampipe->OutputFaceNormal(), -angle));
}

void BDSMagnet::BuildVacuumField()
{
  if (vacuumFieldInfo)
    {
      G4Transform3D newFieldTransform = vacuumFieldInfo->Transform() * beamPipePlacementTransform;
      vacuumFieldInfo->SetTransform(newFieldTransform);
      BDSFieldBuilder::Instance()->RegisterFieldForConstruction(vacuumFieldInfo,
								beampipe->GetContainerLogicalVolume(),
								true);
    }
}

void BDSMagnet::BuildOuter()
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ <<  *magnetOuterInfo << G4endl;
#endif
  G4double outerLength = chordLength - 2*lengthSafety;
  outer = BDSMagnetOuterFactory::Instance()->CreateMagnetOuter(magnetType,
							       magnetOuterInfo,
							       outerLength,
							       chordLength,
							       beampipe);

  if(outer)
    {
      // copy necessary bits out of BDSGeometryComponent that holds
      // container information for whole magnet object provided by
      // magnet outer factory.
      BDSGeometryComponent* container = outer->GetMagnetContainer();
      containerSolid    = container->GetContainerSolid()->Clone();
      G4ThreeVector contOffset = container->GetPlacementOffset();
      // set the main offset of the whole magnet which is placed w.r.t. the
      // zero coordinate of the container solid
      SetPlacementOffset(contOffset);

      RegisterDaughter(outer);
      InheritExtents(container); // update extents

      // Only clear after extents etc have been used
      outer->ClearMagnetContainer();
      
      endPieceBefore = outer->EndPieceBefore();
      endPieceAfter  = outer->EndPieceAfter();

      /// Update record of normal vectors now beam pipe has been constructed.
      SetInputFaceNormal(BDS::RotateToReferenceFrame(outer->InputFaceNormal(), angle));
      SetOutputFaceNormal(BDS::RotateToReferenceFrame(outer->OutputFaceNormal(), -angle));
    }
}

void BDSMagnet::BuildOuterField()
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif  
  if (outer && outerFieldInfo)
    {
      // determine key for this specific magnet instance
      G4String scalingKey = DetermineScalingKey(magnetType);
      
      G4LogicalVolume* vol = outer->GetContainerLogicalVolume();
      BDSFieldBuilder::Instance()->RegisterFieldForConstruction(outerFieldInfo,
								vol,
								true,
								vacuumFieldInfo->MagnetStrength(),
								scalingKey);
      // Attach to the container but don't propagate to daughter volumes. This ensures
      // any gap between the beam pipe and the outer also has a field.
      BDSFieldBuilder::Instance()->RegisterFieldForConstruction(outerFieldInfo,
								containerLogicalVolume,
								false,
								vacuumFieldInfo->MagnetStrength(),
								scalingKey);
    }
}

void BDSMagnet::BuildContainerLogicalVolume()
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif
  // note beam pipe is not optional!
  if (outer)
    {//build around that
      // container solid will have been updated in BuildOuter if the outer exists
      containerLogicalVolume = new G4LogicalVolume(containerSolid,
						   emptyMaterial,
						   name + "_container_lv");

      // user limits
      auto defaultUL = BDSGlobalConstants::Instance()->GetDefaultUserLimits();
      //copy the default and update with the length of the object rather than the default 1m
      G4UserLimits* ul = BDS::CreateUserLimits(defaultUL, std::max(chordLength, arcLength));
      if (ul != defaultUL) // if it's not the default register it
        {RegisterUserLimits(ul);}
      containerLogicalVolume->SetUserLimits(ul);
      containerLogicalVolume->SetVisAttributes(BDSGlobalConstants::Instance()->GetContainerVisAttr());
      
      placeBeamPipe = true;
    }
  else
    {
      // use beam pipe container volume as ours and no need to place beam pipe
      containerSolid         = beampipe->GetContainerSolid();
      containerLogicalVolume = beampipe->GetContainerLogicalVolume();
      InheritExtents(beampipe);
      placeBeamPipe = false;
    }
}

void BDSMagnet::PlaceComponents()
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif
  if (placeBeamPipe)
    {
      G4ThreeVector beamPipeOffset = -1*GetPlacementOffset();
      // place beampipe
      G4PVPlacement* beamPipePV = new G4PVPlacement(nullptr,                 // rotation
						    beamPipeOffset,          // position in container
						    beampipe->GetContainerLogicalVolume(),  // its logical volume
						    name + "_beampipe_pv",   // its name
						    containerLogicalVolume,  // its mother  volume
						    false,                   // no boolean operation
						    0,                       // copy number
                                                    checkOverlaps);
      
      RegisterPhysicalVolume(beamPipePV);
    }

  if (outer)
    {
      //G4ThreeVector placementOffset = magnetOuterOffset + outer->GetPlacementOffset();
      G4ThreeVector outerOffset = outer->GetPlacementOffset();
      
      // place outer volume
      G4PVPlacement* magnetOuterPV = new G4PVPlacement(nullptr,                // rotation
						       outerOffset,            // at normally (0,0,0)
						       outer->GetContainerLogicalVolume(), // its logical volume
						       name+"_outer_pv",       // its name
						       containerLogicalVolume, // its mother  volume
						       false,                  // no boolean operation
						       0,                      // copy number
                                                       checkOverlaps);

      RegisterPhysicalVolume(magnetOuterPV);
    }
}

void BDSMagnet::SetOuterField(BDSFieldInfo* outerFieldInfoIn)
{
  delete outerFieldInfo;
  outerFieldInfo = outerFieldInfoIn;
}

void BDSMagnet::SetVacuumField(BDSFieldInfo* vacuumFieldInfoIn)
{
  delete vacuumFieldInfo;
  vacuumFieldInfo = vacuumFieldInfoIn;
}

BDSMagnet::~BDSMagnet()
{
  delete magnetOuterInfo;
  delete vacuumFieldInfo;
  delete outerFieldInfo;
}
