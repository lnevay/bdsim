/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2020.

This file is part of BDSIM.

BDSIM is free software: you can redistribute it and/or modify 
it under the terms of the GNU General Public License as published 
by the Free Software Foundation version 3 of the License.

BDSIM is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BDSIM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "BDSBeamPipe.hh"
#include "BDSBeamPipeFactory.hh"
#include "BDSBeamPipeInfo.hh"
#include "BDSBeamPipeType.hh"
#include "BDSFieldBuilder.hh"
#include "BDSFieldInfo.hh"
#include "BDSDebug.hh"
#include "BDSException.hh"
#include "BDSMagnetGeometryType.hh"
#include "BDSMagnetOuter.hh"
#include "BDSMagnetOuterInfo.hh"
#include "BDSMagnetOuterFactory.hh"
#include "BDSMagnetStrength.hh"
#include "BDSAcceleratorComponent.hh"
#include "BDSMagnetType.hh"
#include "BDSUtilities.hh"
#include "BDSMagnet.hh"
#include "BDSLine.hh"
#include "BDSComponentFactory.hh"
#include "BDSBeamline.hh"
#include "BDSBeamlineElement.hh"
#include "BDSSimpleComponent.hh"
#include "BDSGlobalConstants.hh"
#include "BDSBendBuilder.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4Types.hh"
#include "G4VPhysicalVolume.hh"
#include "BDSIntegratorSet.hh"
#include "BDSPhysicalVolumeInfo.hh"
#include "BDSPhysicalVolumeInfoRegistry.hh"
#include "BDSAcceleratorComponentRegistry.hh"
#include "BDSGeometryFactoryGDML.hh"
#include "BDSGeometryExternal.hh"



#include "CLHEP/Units/SystemOfUnits.h"

#include "parser/element.h"
#include "parser/elementtype.h"

#include "BDSMagnetNoneSplitOuter.hh"

using namespace GMAD;

class G4Userlimits;

BDSMagnetNoneSplitOuter::BDSMagnetNoneSplitOuter(BDSMagnetType typeIn,
                                                 BDSBeamPipeInfo *beamPipeInfoIn,
                                                 BDSMagnetOuterInfo *magnetOuterInfoIn,
                                                 BDSFieldInfo *vacuumFieldInfoIn,
                                                 BDSFieldInfo *outerFieldInfoIn,
                                                 G4bool isThinIn,
                                                 const Element*          elementIn,
                                                 BDSMagnetStrength*      stIn,
                                                 G4double                brhoIn,
                                                 const BDSIntegratorSet* integratorSetIn,
                                                 G4double                incomingFaceAngleIn,
                                                 G4double                outgoingFaceAngleIn,
                                                 G4bool                  buildFringeFieldsIn,
                                                 const GMAD::Element*    prevElementIn,
                                                 const GMAD::Element*    nextElementIn)
                                                 :BDSMagnet(typeIn, elementIn->name, elementIn->l* CLHEP::m, beamPipeInfoIn, magnetOuterInfoIn, vacuumFieldInfoIn, -elementIn->angle, outerFieldInfoIn, isThinIn)
{
    element           = elementIn;
    st                = stIn;
    brho              = brhoIn;
    integratorSet     = integratorSetIn;
    incomingFaceAngle = incomingFaceAngleIn;
    outgoingFaceAngle = outgoingFaceAngleIn;
    buildFringeFields = buildFringeFieldsIn;
    prevElement       = prevElementIn;
    nextElement       = nextElementIn;
    beamPipeInfo      = beamPipeInfoIn;
    magnetOuterInfo   = magnetOuterInfoIn;
    vacuumFieldInfo   = vacuumFieldInfoIn;
    outerFieldInfo    = outerFieldInfoIn;
}


void BDSMagnetNoneSplitOuter::SBendWithSingleOuter(const G4String&         elementName,
                                                 const Element*          element,
                                                 BDSMagnetStrength*      st,
                                                 G4double                brho,
                                                 const BDSIntegratorSet* integratorSet,
                                                 G4double                incomingFaceAngle,
                                                 G4double                outgoingFaceAngle,
                                                 G4bool                  buildFringeFields,
                                                 const GMAD::Element*    prevElement,
                                                 const GMAD::Element*    nextElement)
{
    Element* el = new Element(*element);
    el->magnetGeometryType = "none";
    BDSAcceleratorComponent* pipeLine = BDS::BuildSBendLine(elementName,el,st,brho,integratorSet,incomingFaceAngle,outgoingFaceAngle,buildFringeFields,prevElement,nextElement);
    delete el;

    pipeLine->Initialise();

    /*BDSLine::iterator it = reinterpret_cast<BDSLine*>(pipeLine)->begin();

    if (buildFringeFields)
    {std::advance(it, 1);}*/

    BDSBeamPipeInfo* beamPipeInfo1 = new BDSBeamPipeInfo(*beamPipeInfo);
    //beamPipeInfo1->inputFaceNormal.set(0.1,0,-0.99498743710662);
    //beamPipeInfo1->outputFaceNormal.set(1,0,0.99);

    BDSBeamPipe* beamPipeTest = BDSBeamPipeFactory::Instance()->CreateBeamPipe(name+"_bp",
                                                   chordLength - 2*lengthSafety,
                                                                           beamPipeInfo1);

    //outer = BDSMagnetOuterFactory::Instance()->CreateMagnetOuter(BDSMagnetType::sectorbend, magnetOuterInfo, chordLength-2*lengthSafety, chordLength, reinterpret_cast<BDSMagnet*>(*it)->BeamPipe());
    outer = BDSMagnetOuterFactory::Instance()->CreateMagnetOuter(BDSMagnetType::sectorbend, magnetOuterInfo, chordLength-2*lengthSafety, chordLength, beamPipeTest);

    delete beamPipeTest;
    delete beamPipeInfo1;

    BDSGeometryComponent* container = outer->GetMagnetContainer();
    containerLogicalVolume = container->GetContainerLogicalVolume();
    containerSolid    = container->GetContainerSolid()->Clone();
    G4ThreeVector contOffset = container->GetPlacementOffset();
    // set the main offset of the whole magnet which is placed w.r.t. the
    // zero coordinate of the container solid
    SetPlacementOffset(contOffset);
    RegisterDaughter(outer);
    InheritExtents(container, contOffset); // update extents
    // Only clear after extents etc have been used
    outer->ClearMagnetContainer();
    endPieceBefore = outer->EndPieceBefore();
    endPieceAfter  = outer->EndPieceAfter();
    /// Update record of normal vectors now beam pipe has been constructed.
    SetInputFaceNormal(BDS::RotateToReferenceFrame(outer->InputFaceNormal(), angle));
    SetOutputFaceNormal(BDS::RotateToReferenceFrame(outer->OutputFaceNormal(), -angle));

    G4double offsetAngle = element->angle/2;
    G4double offsetLength = chordLength/2;

    G4ThreeVector     initialGlobalPosition = G4ThreeVector(0,0, -offsetLength) ;
    G4ThreeVector u = G4ThreeVector( std::cos(-offsetAngle), 0, std::sin(-offsetAngle));
    G4ThreeVector v = G4ThreeVector(0, 1,0);
    G4ThreeVector w = G4ThreeVector(-std::sin(-offsetAngle), 0, std::cos(-offsetAngle));
    G4RotationMatrix* initialGlobalRotation  =  new G4RotationMatrix(u, v, w);
    G4double          initialS = 0;

    BDSBeamline* beamline = new BDSBeamline(initialGlobalPosition, initialGlobalRotation,initialS);

    beamline->AddComponent(pipeLine);

    G4int i = 0;
    for (auto element : *beamline)
    {

        G4String placementName = element->GetPlacementName() + "_pv";
        G4Transform3D* placementTransform = element->GetPlacementTransform();
        G4int copyNumber = i;
        /*auto pv = new G4PVPlacement(*placementTransform,                  // placement transform
                                    element->GetContainerLogicalVolume(), // volume to be placed
                                    placementName,                        // placement name
                                    containerLogicalVolume,                          // volume to place it in
                                    false,                                // no boolean operation
                                    copyNumber,                           // copy number
                                    checkOverlaps);                       // overlap checking
*/
        i++; // for incremental copy numbers

        //RegisterPhysicalVolume(pv);

    }

    G4ThreeVector u1 = G4ThreeVector( 0, -1, 0);
    G4ThreeVector v1 = G4ThreeVector(1, 0,0);
    G4ThreeVector w1 = G4ThreeVector(0, 0, 1);
    G4RotationMatrix* initialGlobalRotation1  =  new G4RotationMatrix(u1, v1, w1);

    //G4ThreeVector outerOffset = outer->GetPlacementOffset();
    G4ThreeVector outerOffset = G4ThreeVector(0,0, 0);

    G4PVPlacement* magnetOuterPV = new G4PVPlacement(initialGlobalRotation1,                // rotation
                                                     outerOffset,            // at normally (0,0,0)
                                                     outer->GetContainerLogicalVolume(), // its logical volume
                                                     name+"_outer_pv",       // its name
                                                     containerLogicalVolume, // its mother  volume
                                                     false,                  // no boolean operation
                                                     0,                      // copy number
                                                     checkOverlaps);
    RegisterPhysicalVolume(magnetOuterPV);

}

void BDSMagnetNoneSplitOuter::Build()
{
   SBendWithSingleOuter(element->name, element, st, brho, integratorSet,
                         incomingFaceAngle, outgoingFaceAngle, buildFringeFields,
                         prevElement, nextElement);

    //BuildUserLimits();
    //BuildVacuumField();
    //BuildOuter();
    // Instead of BDSAcceleratorComponent::Build just call BuildContainerLogicalVolume
    // to control user limits in the case where there is no container and we just inherit
    // the beam pipe container
    //BuildContainerLogicalVolume();
    //BuildOuterField(); // must be done when the containerLV exists
    //PlaceComponents(); // place things (if needed) in container
}
