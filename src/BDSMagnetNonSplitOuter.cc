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
#include "BDSFieldBuilder.hh"
#include "BDSFieldInfo.hh"
#include "BDSDebug.hh"
#include "BDSException.hh"
#include "BDSMagnetOuterInfo.hh"
#include "BDSMagnetOuterFactory.hh"
#include "BDSMagnetStrength.hh"
#include "BDSAcceleratorComponent.hh"
#include "BDSMagnetType.hh"
#include "BDSUtilities.hh"
#include "BDSMagnet.hh"
#include "BDSBeamline.hh"
#include "BDSBeamlineElement.hh"
#include "BDSSimpleComponent.hh"
#include "BDSBendBuilder.hh"
#include "G4Box.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4Types.hh"
#include "BDSIntegratorSet.hh"
#include "BDSGeometryFactoryGDML.hh"
#include "BDSGeometryExternal.hh"

#include "parser/element.h"
#include "parser/elementtype.h"

#include "BDSMagnetNonSplitOuter.hh"

using namespace GMAD;

class G4Userlimits;

BDSMagnetNonSplitOuter::BDSMagnetNonSplitOuter(BDSMagnetType typeIn,
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
    namedVacuumVolumes = BDS::GetWordsFromString(G4String(element->namedVacuumVolumes));
}


void BDSMagnetNonSplitOuter::SBendWithSingleOuter(const G4String&         elementName,
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
    BDSAcceleratorComponent* pipeLine;

    // Creation of the beampipe by calling the BDS::BuildSBendLine with element magnetGeometryType to "none"
    if (element->apertureType != "none")
    {
        Element el = Element(*element);
        el.magnetGeometryType = "none";
        el.l -= lengthSafety;
        pipeLine = BDS::BuildSBendLine(elementName,&el,st,brho,integratorSet,incomingFaceAngle,outgoingFaceAngle,buildFringeFields,prevElement,nextElement);
        pipeLine->Initialise();
    }

    // Create of a beampipe for BDSMagnetOuterFactory::Instance()->CreateMagnetOuter() with the correct geometry information
    BDSBeamPipeInfo* beamPipeInfoTmp = new BDSBeamPipeInfo(*beamPipeInfo);
    std::pair<G4ThreeVector,G4ThreeVector> faces = BDS::CalculateFaces(element->angle/2 + incomingFaceAngle, element->angle/2 + outgoingFaceAngle);
    beamPipeInfoTmp->inputFaceNormal = faces.first;
    beamPipeInfoTmp->outputFaceNormal = faces.second;

    BDSBeamPipe* beamPipeTmp = BDSBeamPipeFactory::Instance()->CreateBeamPipe(name+"_bp",
                                                                               arcLength,
                                                                               beamPipeInfoTmp);

    // Create the Magnet outer
    outer = BDSMagnetOuterFactory::Instance()->CreateMagnetOuter(BDSMagnetType::sectorbend, magnetOuterInfo, chordLength, chordLength, beamPipeTmp);

    delete beamPipeTmp;
    delete beamPipeInfoTmp;

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

    if (!element->extractOuterContainer) // check if the beampipe and outer geometry elements must be placed in the container logical volume
    {

        // Create a beamline receiving the beampipe elements with the correct offset and angle
        G4double offsetAngle = element->angle/2;
        G4double offsetLength = chordLength/2;

        G4ThreeVector     initialGlobalPosition = G4ThreeVector(0,0, -offsetLength) ;
        G4ThreeVector u = G4ThreeVector( std::cos(-offsetAngle), 0, std::sin(-offsetAngle));
        G4ThreeVector v = G4ThreeVector(0, 1,0);
        G4ThreeVector w = G4ThreeVector(-std::sin(-offsetAngle), 0, std::cos(-offsetAngle));
        G4RotationMatrix* initialGlobalRotation  =  new G4RotationMatrix(u, v, w);
        G4double          initialS = 0;

        BDSBeamline* beamline = new BDSBeamline(initialGlobalPosition, initialGlobalRotation,initialS);

        if (element->apertureType != "none") // check if the beampipe elements must be placed in the container logical volume
        {
            beamline->AddComponent(pipeLine);

            G4int i = 0;
            for (auto element : *beamline)
            {

                G4String placementName = element->GetPlacementName() + "_pv";
                G4Transform3D* placementTransform = element->GetPlacementTransform();
                G4int copyNumber = i;
                auto vv = new G4PVPlacement(*placementTransform,                  // placement transform
                                            element->GetContainerLogicalVolume(), // volume to be placed
                                            placementName,                        // placement name
                                            containerLogicalVolume,               // volume to place it in
                                            false,                         // no boolean operation
                                            copyNumber,                           // copy number
                                            false);                       // overlap checking

                i++; // for incremental copy numbers

                if (vv->CheckOverlaps())
                {
                    throw BDSException(__METHOD_NAME__, "Overlapping detected for the beampipe elements");
                }

                RegisterPhysicalVolume(vv);

            }

        }

        G4ThreeVector outerOffset = outer->GetPlacementOffset();

        std::cout << "placing the content of the gdml file" << std::endl;
        auto gdml_world = outer->GetContainerLogicalVolume();

        if (element->includeGdmlWorldVolume) // Check if the outer logical volume (Gdml world volume) must be placed into the container logical volume
        {
            G4Transform3D* placementTransform = new G4Transform3D();

            auto vv = new G4PVPlacement(*placementTransform,                  // placement transform
                                        gdml_world,                           // volume to be placed
                                        gdml_world->GetName() + "_pv", // placement name
                                        containerLogicalVolume,               // volume to place it in
                                        false,                         // no boolean operation
                                        0,                           // copy number
                                        false);                      // overlap checking

            if (vv->CheckOverlaps())
            {
                throw BDSException(__METHOD_NAME__, "Overlapping detected for the outer elements");
            }

            RegisterPhysicalVolume(vv);

        }

        for (G4int j = 0; j < gdml_world->GetNoDaughters(); j++)
        {
            const auto& pv = gdml_world->GetDaughter(j);
            G4String placementName = pv->GetName() + "_pv";
            std::cout << "placing " << placementName << std::endl;
            G4int copyNumber = 1;

            if (!element->includeGdmlWorldVolume) // if the Gdml world volume has not been placed in the container logical volume, placed the outer geometry elements one by one instead
            {
                auto vv = new G4PVPlacement(pv->GetRotation(), pv->GetTranslation(), // placement transform
                                            pv->GetLogicalVolume(),                  // volume to be placed
                                            placementName,                           // placement name
                                            containerLogicalVolume,                  // volume to place it in
                                            false,                            // no boolean operation
                                            copyNumber,                             // copy number
                                            false);                         // overlap checking

                if (vv->CheckOverlaps())
                {
                    throw BDSException(__METHOD_NAME__, "Overlapping detected for the outer elements");
                }

                RegisterPhysicalVolume(vv);
            }

            // register the vacuum field to the external logical volumes registered in namedVacuumVolumes
            if (std::find(namedVacuumVolumes.begin(), namedVacuumVolumes.end(), pv->GetLogicalVolume()->GetName().substr(9).c_str()) != namedVacuumVolumes.end())
            {
                BDSFieldBuilder::Instance()->RegisterFieldForConstruction(vacuumFieldInfo,
                                                                          pv->GetLogicalVolume(),
                                                                          true);
            }

        }
    }

}


void BDSMagnetNonSplitOuter::Build()
{
   SBendWithSingleOuter(element->name, element, st, brho, integratorSet,
                         incomingFaceAngle, outgoingFaceAngle, buildFringeFields,
                         prevElement, nextElement);

   if (element->fieldOuter != "" and !(element->extractOuterContainer)) // check when to build the outer field
   {
       BuildOuterField(); // must be done when the containerLV exists
   }

}
