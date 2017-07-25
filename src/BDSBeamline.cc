#include "globals.hh" // geant4 globals / types
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"

#include "BDSDebug.hh"
#include "BDSAcceleratorComponent.hh"
#include "BDSBeamline.hh"
#include "BDSBeamlineElement.hh"
#include "BDSGlobalConstants.hh"
#include "BDSLine.hh"
#include "BDSOutputBase.hh"
#include "BDSSimpleComponent.hh"
#include "BDSTiltOffset.hh"
#include "BDSTransform3D.hh"
#include "BDSUtilities.hh"

#include <algorithm>
#include <iterator>
#include <ostream>
#include <utility>  // for std::pair
#include <vector>

G4double BDSBeamline::paddingLength = -1;

BDSBeamline::BDSBeamline(G4ThreeVector     initialGlobalPosition,
			 G4RotationMatrix* initialGlobalRotation)
{
  // initialise extents
  totalChordLength      = 0;
  totalArcLength        = 0;
  maximumExtentPositive = G4ThreeVector(0,0,0);
  maximumExtentNegative = G4ThreeVector(0,0,0);
  
  // initial rotation matrix
  if (initialGlobalRotation) // default is null
    {previousReferenceRotationEnd = initialGlobalRotation;}
  else
    {previousReferenceRotationEnd = new G4RotationMatrix();}

#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << "with initial position and rotation" << G4endl;
  G4cout << "Initial position: " << initialGlobalPosition << G4endl;
  G4cout << "Initial rotation: " << *previousReferenceRotationEnd << G4endl;
#endif

  // initial position
  previousReferencePositionEnd = initialGlobalPosition;

  // initial s coordinate
  previousSPositionEnd = 0;

  // gap between each element added to the beam line
  if (paddingLength <= 0)
    {paddingLength = 3 * BDSGlobalConstants::Instance()->LengthSafety();}
  //paddingLength = 3*CLHEP::mm;
}

BDSBeamline::BDSBeamline(G4Transform3D initialTransform):
  BDSBeamline(initialTransform.getTranslation(),
	      new G4RotationMatrix(initialTransform.getRotation()))
{;}

BDSBeamline::~BDSBeamline()
{
  for (iterator it = begin(); it != end(); ++it)
    {delete (*it);}
  // special case, if empty then previousReferenceRotationEnd is not used in the first element
  if (size()==0)
    {delete previousReferenceRotationEnd;}
  // components map goes out of scope - elements are already deleted - no need to
  // explicitly delete
}

void BDSBeamline::PrintAllComponents(std::ostream& out) const
{
  for (const auto& element : *this)
    {out << element;}
}

void BDSBeamline::PrintMemoryConsumption() const
{
  G4cout << __METHOD_NAME__ << "container size: " << sizeof(beamline) << G4endl;
  G4cout << __METHOD_NAME__ << "beamline element cumulative size: " << sizeof(BDSBeamlineElement) * beamline.size() << G4endl;
  G4cout << __METHOD_NAME__ << "full usage including components:  " << (sizeof(BDSBeamlineElement) + sizeof(BDSAcceleratorComponent)) * beamline.size() << G4endl;
}

std::ostream& operator<< (std::ostream& out, BDSBeamline const &bl)
{
  out << "BDSBeamline with "    << bl.size()           << " elements" << G4endl;
  out << "Total arc length:   " << bl.totalArcLength   << " mm"       << G4endl;
  out << "Total chord length: " << bl.totalChordLength << " mm"       << G4endl;

  return out;
}

void BDSBeamline::AddComponent(BDSAcceleratorComponent* component,
			       BDSTiltOffset*           tiltOffset,
			       BDSSamplerType           samplerType,
			       G4String                 samplerName)
{
  if (!component)
    {G4cerr << __METHOD_NAME__ << "invalid accelerator component " << samplerName << G4endl; exit(1);}

  // check the sampler name is allowed in the output
  if (BDSOutputBase::InvalidSamplerName(samplerName))
    {
      G4cerr << __METHOD_NAME__ << "invalid sampler name \"" << samplerName << "\"" << G4endl;
      BDSOutputBase::PrintProtectedNames(G4cerr);
      exit(1);
    }
  
  if (BDSLine* line = dynamic_cast<BDSLine*>(component))
    {
      G4int size = (G4int)line->size();
      for (G4int i = 0; i < size; ++i)
	{
	  if (i < size-1)
	    {AddSingleComponent((*line)[i], tiltOffset);}
	  else // only attach the desired sampler to the last one in the line
	    {AddSingleComponent((*line)[i], tiltOffset, samplerType, samplerName);}
	}
    }
  else
    {AddSingleComponent(component, tiltOffset, samplerType, samplerName);}
  // free memory - as once the rotations are calculated, this is no longer needed
  delete tiltOffset;
}

void BDSBeamline::AddSingleComponent(BDSAcceleratorComponent* component,
				     BDSTiltOffset*           tiltOffset,
				     BDSSamplerType           samplerType,
				     G4String                 samplerName)
{
#ifdef BDSDEBUG
  G4cout << G4endl << __METHOD_NAME__ << "adding component to beamline and calculating coordinates" << G4endl;
  G4cout << "component name:      " << component->GetName() << G4endl;
#endif
  
  // Test if it's a BDSTransform3D instance - this is a unique component that requires
  // rotation in all dimensions and can skip normal addition as isn't a real volume
  // that can be placed.  Apply the transform and skip the rest of this function by returning
  // This modifies the "end" coordinates, rotation and axes of the last element in the beamline
  if (BDSTransform3D* transform = dynamic_cast<BDSTransform3D*>(component))
    {
      ApplyTransform3D(transform);
      return;
    }

  // if it's not a transform3d instance, continue as normal
  // interrogate the item
  G4double      length   = component->GetChordLength();
  G4double      angle    = component->GetAngle();
  G4bool hasFiniteLength = BDS::IsFinite(length);
  G4bool hasFiniteAngle  = BDS::IsFinite(angle);
  G4bool hasFiniteTilt, hasFiniteOffset;
  G4ThreeVector offset;
  if (tiltOffset)
    {
      hasFiniteTilt   = tiltOffset->HasFiniteTilt();
      hasFiniteOffset = tiltOffset->HasFiniteOffset();
      offset          = tiltOffset->GetOffset(); // returns 3Vector
    }
  else
    {
      hasFiniteTilt   = false;
      hasFiniteOffset = false;
      offset          = G4ThreeVector();
    }
  G4ThreeVector eP       = component->GetExtentPositive() + offset;
  G4ThreeVector eN       = component->GetExtentNegative() + offset;
  G4ThreeVector placementOffset   = component->GetPlacementOffset();
  G4bool hasFinitePlacementOffset = BDS::IsFinite(placementOffset);
  G4ThreeVector oFNormal = component->InputFaceNormal();
  
#ifdef BDSDEBUG
  G4cout << "chord length                " << length      << " mm"         << G4endl;
  G4cout << "angle                       " << angle       << " rad"        << G4endl;
  if (tiltOffset)
    {G4cout << "tilt offsetX offsetY        " << *tiltOffset << " rad mm mm " << G4endl;}
  else
    {G4cout << "no tilt offset" << G4endl;}
  G4cout << "has finite length           " << hasFiniteLength              << G4endl;
  G4cout << "has finite angle            " << hasFiniteAngle               << G4endl;
  G4cout << "has finite tilt             " << hasFiniteTilt                << G4endl;
  G4cout << "has finite offset           " << hasFiniteOffset              << G4endl;
  G4cout << "extent positive             " << eP                           << G4endl;
  G4cout << "extent negative             " << eN                           << G4endl;
  G4cout << "object placement offset     " << placementOffset              << G4endl;
  G4cout << "has finite placement offset " << hasFinitePlacementOffset     << G4endl;
  G4cout << "output face normal          " << oFNormal                     << G4endl;
#endif

  // Check this won't overlap with any previous geometry. This is only done for elements
  // that aren't drifts as they should be built by the component factory to match any angles.
  if (!empty() && (component->GetType() != "drift") && (component->GetType() != "thinmultipole"))
    {// can only look back if there is an element - won't clash if no element; also add drifts always
      G4bool   keepGoing   = true;
      G4bool   checkFaces  = true;
      G4double zSeparation = 0;
      const BDSBeamlineElement* inspectedElement = back(); // remember we haven't added this new element yet
      // find previous non drift output face.
      G4ThreeVector iFNormal;
      G4String clasherName = "Unknown";
      while (keepGoing)
	{
	  if (inspectedElement) // valid element
	    {// decrement could return nullptr so have to check if valid element
	      if ((inspectedElement->GetType() == "drift")||(inspectedElement->GetType() == "thinmultipole")) // leave keepGoing true
		{
		  zSeparation += inspectedElement->GetChordLength();
		  inspectedElement = GetPrevious(inspectedElement); // decrement
		}
	      else
		{
		  keepGoing   = false; // found a non drift - stop here
		  iFNormal    = inspectedElement->GetAcceleratorComponent()->OutputFaceNormal();
		  clasherName = inspectedElement->GetAcceleratorComponent()->GetName();
		}
	    }
	  else
	    {
	      keepGoing  = false;
	      checkFaces = false; // got to the beginning with only drifts - don't check
	    }
	}
#ifdef BDSDEBUG
      G4cout << "input face normal           " << iFNormal << G4endl; // matches above debug formatting
#endif

      if (checkFaces)
	{
	  // now do checks
	  BDSExtent extOF = inspectedElement->GetAcceleratorComponent()->GetExtent(); // output face
	  BDSExtent extIF = component->GetExtent(); // input face
	  
	  G4bool willIntersect = BDS::WillIntersect(iFNormal, oFNormal, zSeparation, extIF, extOF);
	  if (willIntersect)
	    {
	      G4cout << __METHOD_NAME__ << "Error - angled faces of objects will cause overlap in beam line geometry" << G4endl;
	      G4cout << "\"" << component->GetName() << "\" will overlap with \""
		     << clasherName << "\"" << G4endl;
	      exit(1);
	    }
	}
    }
  
  // Calculate the reference placement rotation
  // rotations are done first as they're required to transform the spatial displacements.
  // if not the first element in the beamline, copy the rotation matrix (cumulative along line)
  // from end of last component, else use initial rotation matrix (no copy to prevent memory leak)
  G4RotationMatrix* referenceRotationStart;
  if (empty())
    {referenceRotationStart = previousReferenceRotationEnd;}
  else
    {
      previousReferenceRotationEnd = back()->GetReferenceRotationEnd();
      referenceRotationStart  = new G4RotationMatrix(*previousReferenceRotationEnd);
    }

  G4RotationMatrix* referenceRotationMiddle = new G4RotationMatrix(*referenceRotationStart);
  G4RotationMatrix* referenceRotationEnd    = new G4RotationMatrix(*referenceRotationStart);
  
  // if the component induces an angle in the reference trajectory, rotate the mid and end point
  // rotation matrices appropriately
  if (hasFiniteAngle)
    {
      // remember our definition of angle - +ve angle bends in -ve x direction in right
      // handed coordinate system
      // rotate about cumulative local y axis of beamline
      // middle rotated by half angle in local x,z plane
      G4ThreeVector rotationAxisOfBend = G4ThreeVector(0,1,0); // nominally about local unit Y
      G4ThreeVector rotationAxisOfBendEnd = rotationAxisOfBend; // a copy
      if (hasFiniteTilt)
	{
	  G4double tilt = tiltOffset->GetTilt();
	  G4RotationMatrix rotationAxisRM = G4RotationMatrix();
	  rotationAxisRM.rotateZ(tilt);
	  rotationAxisOfBend.transform(rotationAxisRM);
	  rotationAxisOfBendEnd.transform(rotationAxisRM);
	}
      referenceRotationMiddle->rotate(angle*0.5, rotationAxisOfBend.transform(*previousReferenceRotationEnd));
      // end rotated by full angle in local x,z plane
      referenceRotationEnd->rotate(angle, rotationAxisOfBendEnd.transform(*previousReferenceRotationEnd));
    }
  
  G4RotationMatrix* rotationStart  = new G4RotationMatrix(*referenceRotationStart);
  G4RotationMatrix* rotationMiddle = new G4RotationMatrix(*referenceRotationMiddle);
  G4RotationMatrix* rotationEnd    = new G4RotationMatrix(*referenceRotationEnd);
  // add the tilt to the rotation matrices (around z axis)
  if (hasFiniteTilt)
    {
      G4double tilt = tiltOffset->GetTilt();

      // transform a unit z vector with the rotation matrices to get the local axes
      // of rotation to apply the tilt.
      G4ThreeVector unitZ = G4ThreeVector(0,0,1);
      rotationStart ->rotate(tilt, unitZ.transform(*referenceRotationStart));
      unitZ = G4ThreeVector(0,0,1);
      rotationMiddle->rotate(tilt, unitZ.transform(*referenceRotationMiddle));
      unitZ = G4ThreeVector(0,0,1);
      rotationEnd   ->rotate(tilt, unitZ.transform(*referenceRotationEnd));
    }
  
  // calculate the reference placement position
  // if not the first item in the beamline, get the reference trajectory global position
  // at the end of the previous element
  if (!empty())
    {
      previousReferencePositionEnd = back()->GetReferencePositionEnd();
      // leave a small gap for unambiguous geometry navigation. Transform that length
      // to a unit z vector along the direction of the beam line before this component
      G4ThreeVector componentGap = G4ThreeVector(0,0,paddingLength).transform(*referenceRotationStart);
      previousReferencePositionEnd += componentGap;
    }
  
  G4ThreeVector referencePositionStart, referencePositionMiddle, referencePositionEnd;
  if (hasFiniteLength)
    {
      referencePositionStart = previousReferencePositionEnd;
      
      // calculate delta to mid point
      G4ThreeVector md = G4ThreeVector(0, 0, 0.5 * length);
      md.transform(*referenceRotationMiddle);
      referencePositionMiddle = referencePositionStart + md;
      // remember the end position is the chord length along the half angle, not the full angle
      // the particle achieves the full angle though by the end position.
      G4ThreeVector delta = G4ThreeVector(0, 0, length).transform(*referenceRotationMiddle);
      referencePositionEnd = referencePositionStart + delta;
    }
  else
    {
      // element has no finite size so all positions are previous end position
      // likely this is a transform3d or similar - but not hard coded just for transform3d
      referencePositionStart  = previousReferencePositionEnd;
      referencePositionMiddle = previousReferencePositionEnd;
      referencePositionEnd    = previousReferencePositionEnd;
    }
  
  // add the placement offset
  G4ThreeVector positionStart, positionMiddle, positionEnd;
  if (hasFiniteOffset || hasFinitePlacementOffset)
    {
      if (hasFiniteOffset && hasFiniteAngle) 
	{// do not allow x offsets for bends as this will cause overlaps
	  G4String name = component->GetName();
	  G4cout << __METHOD_NAME__ << "WARNING - element has x offset, but this will cause geometry"
		 << " overlaps: " << name << " - omitting x offset" << G4endl;
	  offset.setX(0.0);
	}
      // note the displacement is applied in the accelerator x and y frame so use
      // the reference rotation rather than the one with tilt already applied
      G4ThreeVector total = offset + placementOffset;
      G4ThreeVector displacement   = total.transform(*referenceRotationMiddle);
      positionStart  = referencePositionStart  + displacement;
      positionMiddle = referencePositionMiddle + displacement;
      positionEnd    = referencePositionEnd    + displacement;
    }
  else
    {
      positionStart  = referencePositionStart;
      positionMiddle = referencePositionMiddle;
      positionEnd    = referencePositionEnd;
    }
  
  // calculate the s position
  // if not the first element in the beamline, get the s position at the end of the previous element
  if (!empty())
    {previousSPositionEnd = back()->GetSPositionEnd();}
  
  G4double arcLength   = component->GetArcLength();
  G4double chordLength = component->GetChordLength();

  // integrate lengths
  totalChordLength += chordLength;
  totalArcLength   += arcLength;

  // advance s coordinate
  G4double sPositionStart, sPositionMiddle, sPositionEnd;
  sPositionStart  = previousSPositionEnd;
  sPositionMiddle = previousSPositionEnd + 0.5 * arcLength;
  sPositionEnd    = previousSPositionEnd + arcLength;

  // update the global constants
  BDSGlobalConstants::Instance()->SetSMax(sPositionEnd);

#ifdef BDSDEBUG
  // feedback about calculated coordinates
  G4cout << "calculated coordinates in mm and rad are " << G4endl;
  G4cout << "reference position start:  " << referencePositionStart   << G4endl;
  G4cout << "reference position middle: " << referencePositionMiddle  << G4endl;
  G4cout << "reference position end:    " << referencePositionEnd     << G4endl;
  G4cout << "reference rotation start:  " << *referenceRotationStart;
  G4cout << "reference rotation middle: " << *referenceRotationMiddle;
  G4cout << "reference rotation end:    " << *referenceRotationEnd;
  G4cout << "position start:            " << positionStart            << G4endl;
  G4cout << "position middle:           " << positionMiddle           << G4endl;
  G4cout << "position end:              " << positionEnd              << G4endl;
  G4cout << "rotation start:            " << *rotationStart;
  G4cout << "rotation middle:           " << *rotationMiddle;
  G4cout << "rotation end:              " << *rotationEnd;
#endif
  
  // construct beamline element
  BDSTiltOffset* tiltOffsetToStore = nullptr;
  if (tiltOffset)
    {tiltOffsetToStore = new BDSTiltOffset(*tiltOffset);} // copy as can be used multiple times
  
  BDSBeamlineElement* element;
  element = new BDSBeamlineElement(component,
				   positionStart,
				   positionMiddle,
				   positionEnd,
				   rotationStart,
				   rotationMiddle,
				   rotationEnd,
				   referencePositionStart,
				   referencePositionMiddle,
				   referencePositionEnd,
				   referenceRotationStart,
				   referenceRotationMiddle,
				   referenceRotationEnd,
				   sPositionStart,
				   sPositionMiddle,
				   sPositionEnd,
				   tiltOffsetToStore,
				   samplerType,
				   samplerName,
				   (G4int)beamline.size());

  // calculate extents for world size determination
  UpdateExtents(element);
  
  // append it to the beam line
  beamline.push_back(element);

  // register the s position at the end for curvilinear transform
  sEnd.push_back(sPositionEnd);

  // register it by name
  RegisterElement(element);

#ifdef BDSDEBUG
  G4cout << *element;
  G4cout << __METHOD_NAME__ << "component added" << G4endl;
#endif
}

void BDSBeamline::ApplyTransform3D(BDSTransform3D* component)
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << "- as it's a transform3d instance" << G4endl;
#endif
  // interrogate component
  G4double dx     = component->GetDX();
  G4double dy     = component->GetDY();
  G4double dz     = component->GetDZ();
  G4double dTheta = component->GetDTheta();
  G4double dPsi   = component->GetDPsi();
  G4double dPhi   = component->GetDPhi();
  
  // debug feedback
#ifdef BDSDEBUG
  G4cout << "dx     " << dx     << G4endl;
  G4cout << "dy     " << dy     << G4endl;
  G4cout << "dz     " << dz     << G4endl;
  G4cout << "dTheta " << dTheta << G4endl;
  G4cout << "dPsi   " << dPsi   << G4endl;
  G4cout << "dPhi   " << dPhi   << G4endl;
#endif

  // test validity for overlaps
  if (dz < 0)
    {
      G4cerr << __METHOD_NAME__ << "Problem with Transform3d: " << component->GetName() << G4endl;
      G4cerr << __METHOD_NAME__ << "dz = " << dz << " < 0 -> will overlap previous element" << G4endl;
    } 

  // if not the first element in the beamline, get information from the end of the
  // last element in the beamline
  if (!empty())
    {
      BDSBeamlineElement* last = back();
      previousReferenceRotationEnd = last->GetReferenceRotationEnd();
      previousReferencePositionEnd = last->GetReferencePositionEnd();
    }

  // apply position
  // transform the local dx,dy,dz displacement into the global frame then apply
  G4ThreeVector delta = G4ThreeVector(dx,dy,dz).transform(*previousReferenceRotationEnd);
  previousReferencePositionEnd = previousReferencePositionEnd + G4ThreeVector(dx,dy,dz);

#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << "existing end rotation matrix:" << *previousReferenceRotationEnd << G4endl;
#endif
  
  // apply rotation
  // euler angles must be applied in sequence about the cumulatively rotated axes
  // use unit vectors that are transformed to the current cumulative rotation of the beamline
  // as rotation axes for each angle.
  G4ThreeVector unitZ = G4ThreeVector(0,0,1);
  // transform to the current local z axis at the end of the beamline
  unitZ.transform(*previousReferenceRotationEnd);
  // apply the phi (alpha) euler angle about the z axis (step1)
  previousReferenceRotationEnd->rotate(dPhi, unitZ);
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << "after apply Phi:" << *previousReferenceRotationEnd << G4endl;
#endif
  // apply the theta (beta) euler angle about the N or x' axis (step 2)
  // transform a unit x to axes rotated by phi (called N or x')
  G4ThreeVector unitX = G4ThreeVector(1,0,0);
  unitX.transform(*previousReferenceRotationEnd);
  // rotate about N by dTheta
  previousReferenceRotationEnd->rotate(dTheta, unitX);
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << "after apply Phi & Theta:" << *previousReferenceRotationEnd << G4endl;
#endif
  // apply the psi (gamma) euler angle about the Z or z'' axis
  G4ThreeVector unitZPP = G4ThreeVector(0,0,1); //Z Prime Prime
  // get unit z to local axes already rotated by phi and theta hence zpp name
  unitZPP.transform(*previousReferenceRotationEnd);
  // rotate by psi about zpp
  previousReferenceRotationEnd->rotate(dPsi, unitZPP);
  
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << "new end rotation matrix:" << *previousReferenceRotationEnd << G4endl;
#endif
}

void BDSBeamline::AddBeamlineElement(BDSBeamlineElement* element)
{
  if (!element)
    {G4cerr << __METHOD_NAME__ << "invalid BDSBeamlineElement" << G4endl; exit(1);}
  if (!(element->GetAcceleratorComponent()))
    {G4cerr << __METHOD_NAME__ << "invalid BDSAcceleratorComponent" << G4endl; exit(1);}
  
  // update world extent for this beam line
  UpdateExtents(element);
  
  // append it to the beam line
  beamline.push_back(element);

  // register it by name
  RegisterElement(element);

  // no need to update any internal variables - that's done by AddSingleComponent()
}

G4ThreeVector BDSBeamline::GetMaximumExtentAbsolute() const
{
  G4ThreeVector mEA;
  for (int i=0; i<3; i++)
    {mEA[i] = std::max(std::abs(maximumExtentPositive[i]), std::abs(maximumExtentNegative[i]));}
  return mEA;
}

G4Transform3D BDSBeamline::GetGlobalEuclideanTransform(G4double s, G4double x, G4double y) const
{
  // check if s is in the range of the beamline
  if (s > totalArcLength)
    {
      G4cout << __METHOD_NAME__
	     << "s position \"" << s << "\" is beyond length of accelerator" << G4endl;
      G4cout << "Returning 0 transform" << G4endl;
      return G4Transform3D();
    }

  // find element that s position belongs to
  auto lower = std::lower_bound(sEnd.begin(), sEnd.end(), s);
  G4int index = lower - sEnd.begin(); // subtract iterators to get index
  const BDSBeamlineElement* element = beamline.at(index);
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
  G4cout << "S position requested: " << s     << G4endl;
  G4cout << "Index:                " << index << G4endl;
  G4cout << "Element: " << *element << G4endl;
#endif

  G4double dx = 0;
  // G4double dy = 0; // currently magnets can only bend in local x so avoid extra calculation

  // difference from centre of element to point in local coords)
  // difference in s from centre, normalised to arcLengh and scaled to chordLength
  // as s is really arc length but we must place effectively in chord length coordinates
  const BDSAcceleratorComponent* component = element->GetAcceleratorComponent();
  G4double arcLength   = component->GetArcLength();
  G4double chordLength = component->GetChordLength();
  G4double dS          = s - element->GetSPositionMiddle();
  G4double localZ      = dS * (chordLength / arcLength);
  G4double angle       = component->GetAngle();
  G4RotationMatrix rotation; // will be interpolated rotation
  G4RotationMatrix* rotMiddle = element->GetReferenceRotationMiddle();
  // find offset of point from centre of volume - 2 methods
  if (BDS::IsFinite(angle))
    {
      // finite bend angle - interpolate position and angle along arc due to change in angle
      // local unit z at start of element
      G4ThreeVector localUnitY = G4ThreeVector(0,1,0);
      localUnitY.transform(*(element->GetReferenceRotationStart()));
      // linearly interpolate angle -> angle * (s from beginning into component)/arcLength
      G4double partialAngle = angle * std::fabs(( (0.5*arcLength + dS) / arcLength));
      rotation = G4RotationMatrix(*(element->GetReferenceRotationStart())); // start rotation
      rotation.rotate(partialAngle, localUnitY); // rotate it by the partial angle about local Y
      dx = localZ*tan(partialAngle); // calculate difference of arc from chord at that point
    }
  else
    {rotation = G4RotationMatrix(*rotMiddle);}

  // note, magnets only bend in local x so no need to add dy as always 0
  G4ThreeVector dLocal    = G4ThreeVector(x + dx, y /*+ dy*/, localZ);
#ifdef BDSDEBUG
  G4cout << "Local offset from middle: " << dLocal << G4endl;
#endif
  // note, rotation middle is also the same as the coordinate frame of the g4 solid
  G4ThreeVector globalPos = element->GetReferencePositionMiddle() + dLocal.transform(*rotMiddle);
  // construct transform3d from global position and rotation matrix
  G4Transform3D result    = G4Transform3D(rotation, globalPos);
  
#ifdef BDSDEBUG
  G4cout << "Global offset from middle: " << dLocal    << G4endl;
  G4cout << "Resultant global position: " << globalPos << G4endl;
#endif
  return result;
}

const BDSBeamlineElement* BDSBeamline::GetPrevious(const BDSBeamlineElement* element) const
{
  // search for element
  auto result = find(beamline.begin(), beamline.end(), element);
  if (result != beamline.end())
    {// found
      return GetPrevious(result - beamline.begin());
    }
  else
    {return nullptr;}
}

const BDSBeamlineElement* BDSBeamline::GetPrevious(G4int index) const
{
  if (index < 1 || index > (G4int)(beamline.size()-1))
    {return nullptr;} // invalid index - inc beginning or end
  else
    {return beamline[index-1];}
}

const BDSBeamlineElement* BDSBeamline::GetNext(const BDSBeamlineElement* element) const
{
  // search for element
  auto result = find(beamline.begin(), beamline.end(), element);
  if (result != beamline.end())
    {// found
      return GetNext(result - beamline.begin());
    }
  else
    {return nullptr;}
}

const BDSBeamlineElement* BDSBeamline::GetNext(G4int index) const
{
  if (index < 0 || index > (G4int)(beamline.size()-2))
    {return nullptr;} // invalid index - inc beginning or end
  else
    {return beamline[index+1];}
}

void BDSBeamline::RegisterElement(BDSBeamlineElement* element)
{
  // check if base name already registered (can be single component placed multiple times)
  std::map<G4String, BDSBeamlineElement*>::iterator search = components.find(element->GetName());
  if (search == components.end())
    {// not registered
      components[element->GetPlacementName()] = element;
    }
}

BDSBeamlineElement* BDSBeamline::GetElement(G4String acceleratorComponentName,
					    G4int    i) const
{
  // build placement name based on acc component name and ith placement
  // matches construction in BDSBeamlineElement
  G4String placementName = acceleratorComponentName + "_" + std::to_string(i);
  const auto search = components.find(placementName);
  if (search == components.end())
    {return nullptr;} //wasn't found
  else
    {return search->second;}
}

G4Transform3D BDSBeamline::GetTransformForElement(G4String acceleratorComponentName,
						  G4int    i) const
{
  BDSBeamlineElement* result = GetElement(acceleratorComponentName, i);
  if (!result)
    {
      G4cerr << __METHOD_NAME__ << "No element named \""
	     << acceleratorComponentName << "\" found for placement number "
	     << i << G4endl;
      exit(1);
    }
  else
    {return G4Transform3D(*(result->GetRotationMiddle()), result->GetPositionMiddle());}
}

void BDSBeamline::UpdateExtents(BDSBeamlineElement* element)
{
  // calculate extents for world size determination
  // project size in global coordinates
  G4ThreeVector     referencePositionStart = element->GetReferencePositionStart();
  G4RotationMatrix* referenceRotationStart = element->GetReferenceRotationStart();
  G4ThreeVector     referencePositionEnd   = element->GetReferencePositionEnd();
  G4RotationMatrix* referenceRotationEnd   = element->GetReferenceRotationEnd();
  BDSAcceleratorComponent* component       = element->GetAcceleratorComponent();
  G4ThreeVector eP                         = component->GetExtentPositive();
  eP.setZ(0); // we get the z position from the start point, so only need the transverse bits
  G4ThreeVector eN                         = component->GetExtentNegative();
  eN.setZ(0); // we get the z position from the start point, so only need the transverse bits
  G4ThreeVector ePStart                    = G4ThreeVector(eP).transform(*referenceRotationStart);
  G4ThreeVector eNStart                    = G4ThreeVector(eN).transform(*referenceRotationStart);
  G4ThreeVector ePEnd                      = G4ThreeVector(eP).transform(*referenceRotationEnd);
  G4ThreeVector eNEnd                      = G4ThreeVector(eN).transform(*referenceRotationEnd);
  G4ThreeVector extentposStart             = referencePositionStart + ePStart;
  G4ThreeVector extentnegStart             = referencePositionStart + eNStart;
  G4ThreeVector extentposEnd               = referencePositionEnd   + ePEnd;
  G4ThreeVector extentnegEnd               = referencePositionEnd   + eNEnd;
 
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
  G4cout << "start position (global):       " << referencePositionStart << G4endl;
  G4cout << "end position (global):         " << referencePositionEnd   << G4endl;
  G4cout << "local extent +ve:              " << eP                     << G4endl;
  G4cout << "local extent -ve:              " << eN                     << G4endl;
  G4cout << "extent +ve at start in global: " << ePStart                << G4endl;
  G4cout << "extent -ve at start in global: " << eNStart                << G4endl;
  G4cout << "extent +ve at end in global:   " << ePEnd                  << G4endl;
  G4cout << "extent -ve at end in global:   " << eNEnd                  << G4endl;
  G4cout << "current global extent +ve:     " << maximumExtentPositive  << G4endl;
  G4cout << "current global extent -ve:     " << maximumExtentNegative  << G4endl;
#endif
  
  // loop over each size and compare to cumulative extent
  // do this at the start and end to be sure for long components
  // start
  for (int i=0; i<3; i++)
    {
      if (extentposStart[i] > maximumExtentPositive[i])
	{maximumExtentPositive[i] = extentposStart[i];}
      if (extentnegStart[i] < maximumExtentNegative[i])
	{maximumExtentNegative[i] = extentnegStart[i];}
    }
  // end
  for (int i=0; i<3; i++)
    {
      if (extentposEnd[i] > maximumExtentPositive[i])
	{maximumExtentPositive[i] = extentposEnd[i];}
      if (extentnegEnd[i] < maximumExtentNegative[i])
	{maximumExtentNegative[i] = extentnegEnd[i];}
    }
  // end comparing negative extents with positive world just in case
  for (int i=0; i<3; i++)
    {
      if (extentnegEnd[i] > maximumExtentPositive[i])
	{maximumExtentPositive[i] = extentnegEnd[i];}
    }
#ifdef BDSDEBUG
  G4cout << "new global extent +ve:         " << maximumExtentPositive << G4endl;
  G4cout << "new global extent -ve:         " << maximumExtentNegative << G4endl;
#endif
}

BDSBeamlineElement* BDSBeamline::ProvideEndPieceElementBefore(BDSSimpleComponent* endPiece,
							      G4int               index) const
{
  if (!IndexOK(index))
    {return nullptr;}

  G4double endPieceLength      = endPiece->GetChordLength();
  BDSBeamlineElement*  element = beamline[index];
  G4RotationMatrix* elRotStart = element->GetRotationMiddle();
  G4ThreeVector     elPosStart = element->GetPositionStart();
  G4ThreeVector positionMiddle = elPosStart - G4ThreeVector(0,0,endPieceLength*0.5).transform(*elRotStart);
  G4ThreeVector  positionStart = elPosStart - G4ThreeVector(0,0,endPieceLength).transform(*elRotStart);
  G4double         elSPosStart = element->GetSPositionStart();
  BDSTiltOffset*  elTiltOffset = element->GetTiltOffset();
  BDSTiltOffset*   forEndPiece = nullptr;
  if (elTiltOffset)
    {forEndPiece = new BDSTiltOffset(*elTiltOffset);}
  BDSBeamlineElement* result = new BDSBeamlineElement(endPiece,
						      positionStart,
						      positionMiddle,
						      elPosStart,
						      new G4RotationMatrix(*elRotStart),
						      new G4RotationMatrix(*elRotStart),
						      new G4RotationMatrix(*elRotStart),
						      positionStart,// for now the same
						      positionMiddle,
						      elPosStart,
						      new G4RotationMatrix(*elRotStart),
						      new G4RotationMatrix(*elRotStart),
						      new G4RotationMatrix(*elRotStart),
						      elSPosStart - endPieceLength,
						      elSPosStart - 0.5*endPieceLength,
						      elSPosStart,
                                                      forEndPiece);
  return result;
}

BDSBeamlineElement* BDSBeamline::ProvideEndPieceElementAfter(BDSSimpleComponent* endPiece,
							     G4int               index,
							     G4bool              flip) const
{
  if (!IndexOK(index))
    {return nullptr;}

  G4double endPieceLength      = endPiece->GetChordLength();
  BDSBeamlineElement*  element = beamline[index];
  G4RotationMatrix*   elRotEnd = new G4RotationMatrix(*(element->GetRotationMiddle()));
  G4ThreeVector       elPosEnd = element->GetPositionEnd();
  G4ThreeVector positionMiddle = elPosEnd + G4ThreeVector(0,0,endPieceLength*0.5).transform(*elRotEnd);
  G4ThreeVector    positionEnd = elPosEnd + G4ThreeVector(0,0,endPieceLength).transform(*elRotEnd);
  if (flip)
    {elRotEnd->rotateY(CLHEP::pi);}
  G4double           elSPosEnd = element->GetSPositionEnd();
  BDSTiltOffset*  elTiltOffset = element->GetTiltOffset();
  BDSTiltOffset*   forEndPiece = nullptr;
  if (elTiltOffset)
    {forEndPiece = new BDSTiltOffset(*elTiltOffset);}
  BDSBeamlineElement* result = new BDSBeamlineElement(endPiece,
						      elPosEnd,
						      positionMiddle,
						      positionEnd,
						      new G4RotationMatrix(*elRotEnd),
						      new G4RotationMatrix(*elRotEnd),
						      new G4RotationMatrix(*elRotEnd),
						      elPosEnd,
						      positionMiddle,
						      positionEnd,
						      new G4RotationMatrix(*elRotEnd),
						      new G4RotationMatrix(*elRotEnd),
						      new G4RotationMatrix(*elRotEnd),
						      elSPosEnd,
						      elSPosEnd + 0.5*endPieceLength,
						      elSPosEnd + endPieceLength,
						      forEndPiece);
  return result;
}

G4bool BDSBeamline::IndexOK(G4int index) const
{
  if (index < 0 || index > (G4int)(beamline.size()-1))
    {return false;}
  else
    {return true;}
}

std::vector<G4double> BDSBeamline::GetEdgeSPositions()const
{
  std::vector<G4double> sPos;
  sPos.reserve(beamline.size()+1);
  // add start position
  sPos.push_back(0.0);
  for (auto element : beamline)
    {sPos.push_back(element->GetSPositionEnd()/CLHEP::m);}
  return sPos;
}
