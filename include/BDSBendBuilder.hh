#ifndef BDSBENDBUILDER_H
#define BDSBENDBUILDER_H

#include "globals.hh" // geant4 globals / types

class BDSIntegratorSet;
class BDSLine;
class BDSMagnet;
class BDSMagnetStrength;

namespace GMAD
{
  struct Element;
}

namespace BDS
{  
  /// This calculates and constructs a BDSLine* of BDSMagnet*. The bend is split
  /// into multiple wedges. If a small or zero poleface angle is specified, the
  /// end wedges will have faces angled as appropriate, the remaining wedges will
  /// re-use a single identical 'central' wedge several times. For a stronger angled
  /// poleface, the faces of each wedge fade in/out from the poleface to the cental
  /// wedge in the middle. Thin fringefield elements are placed at the beginning and
  /// end of the beamline if required.

  BDSAcceleratorComponent* BuildSBendLine(const GMAD::Element*    element,
					  BDSMagnetStrength*      st,
					  const G4double          brho,
					  const BDSIntegratorSet* integratorSet);
  
  /// Construct beamline for an rbend.  A line is returned with a single
  /// magnet as the main dipole, but can have fringefield magnets placed
  /// either end if specified.
  BDSLine* BuildRBendLine(const GMAD::Element*    element,
			  const GMAD::Element*    prevElement,
			  const GMAD::Element*    nextElement,
			  const G4double          brho,
			  BDSMagnetStrength*      st,
			  const BDSIntegratorSet* integratorSet);

  /// Utility function to calculate the number of segments an sbend should be split into.
  /// Based on aperture error tolerance - default is 1mm.
  G4int CalculateNSBendSegments(const G4double length,
				const G4double angle,
				const G4double e1 = 0,
				const G4double e2 = 0,
				const G4double aperturePrecision = 1.0);

  /// Thin magnet for dipole fringe field.
  /// Is beampipe only, no outer magnet.
  BDSMagnet* BuildDipoleFringe(const GMAD::Element*     element,
			       G4double                 angleIn,
			       G4double                 angleOut,
			       G4String                 name,
			       const BDSMagnetStrength* st,
			       G4double                 brho,
			       const BDSIntegratorSet*  integratorSet);

  /// Function to return a single secotr bend section.
  BDSMagnet* BuildSingleSBend(const GMAD::Element*     element,
			      const G4String           name,
			      const G4double           arcLength,
			      const G4double           angle,
			      const G4double           angleIn,
			      const G4double           angleOut,
			      const BDSMagnetStrength* strength,
			      const G4double           brho,
			      const BDSIntegratorSet*  integratorSet,
			      const G4bool             yokeOnLeft);
  
  /// Function to calculate the value of the fringe field correction term.
  G4double CalculateFringeFieldCorrection(G4double rho,
					  G4double polefaceAngle,
					  G4double fint,
  					  G4double hgap);

  void UpdateSegmentAngles(const G4int    index,
			   const G4int    nSBends,
			   const G4double semiAngle,
			   const G4double e1,
			   const G4double e2,
			   G4double&      segmentAngleIn,
			   G4double&      segmentAngleOut);
}

#endif
