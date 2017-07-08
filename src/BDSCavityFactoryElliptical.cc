#include "BDSAcceleratorComponent.hh"
#include "BDSCavity.hh"
#include "BDSCavityFactoryElliptical.hh"
#include "BDSCavityInfo.hh"
#include "BDSColours.hh"
#include "BDSGlobalConstants.hh"

#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4Polycone.hh"
#include "G4SubtractionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4VSolid.hh"

#include "CLHEP/Units/SystemOfUnits.h"

#include <cmath>
#include <vector>

BDSCavityFactoryElliptical::BDSCavityFactoryElliptical()
{;}

BDSCavityFactoryElliptical::~BDSCavityFactoryElliptical()
{;}

BDSCavity* BDSCavityFactoryElliptical::CreateCavity(G4String             name,
						    G4double             totalChordLength,
						    const BDSCavityInfo* info,
						    G4Material*          vacuumMaterial)
{
  CleanUp();
  
  G4double lengthSafetyLarge = BDSAcceleratorComponent::lengthSafetyLarge;
  G4double chordLength = totalChordLength;
  // Elliptical Cavity Parameters
  // irisRSemiAxis    --> Semi-axis of the iris ellipse perpendicular to the length of the cavity.
  // irisZSemiAxis    --> Semi-axis of the iris ellipse along the length of the cavity.
  // equatorRSemiAxis --> Semi-axis of the equator ellipse perpendicular to the length of the cavity.
  // equatorZSemiAxis --> Semi-axis of the equator ellipse along the length of the cavity.
  // tangentAngle     --> Angle between the line vertical to the z-axis and the line tangent to both ellipses.
  // irisRadius       --> Radius of the iris.
  // noPoints         --> Number of points used to define the shape.  This is divided up between
  //                     the 3 constituent shapes, 1/4 in each iris ellipse and 1/2 in central dome.
  // length           --> Length of the cell.
  // equatorRadius    --> Radius at the equator.
  G4double equatorRadius    = info->equatorRadius;
  G4double irisRSemiAxis    = info->irisVerticalAxis;
  G4double irisZSemiAxis    = info->irisHorizontalAxis;     //iris ellipse horizontal semiaxis
  G4double equatorRSemiAxis = info->equatorEllipseSemiAxis ;//equator ellipse vertical semiaxis
  G4double equatorZSemiAxis = info->equatorEllipseSemiAxis; //equator ellipse horizontal semiaxis
  G4double tangentAngle     = info->tangentLineAngle;
  G4double irisRadius       = info->irisRadius;
  G4double thickness        = info->thickness;
  unsigned int noPoints     = info->numberOfPoints;

  // calculate cartesian coordinates (z, r) from parameters.
  // 2D cylindrical coordinates, z along the beamline:
  G4double zi = chordLength * 0.5;                // z coord of iris ellipse centre
  G4double ri = irisRadius + irisRSemiAxis;       // r coord of iris ellipse centre
  G4double ze = 0.0;                              // z coord of equator ellipse centre
  G4double re = equatorRadius - equatorRSemiAxis; // r coord of equator ellipse centre.
  // gradient of line connecting the ellipses.
  // Add a pi/2 because angle is defined from the vertical, clockwise.
  G4double m = tan(tangentAngle + 0.5*CLHEP::pi); 
   
  // Gradient from tangentAngle. Find the derivative of ellipses. Equate
  // and solve for the parameter.
  // atan finds solution in the first quadrant.
  G4double equatorParameterTangentPoint = atan(-equatorRSemiAxis/(m*equatorZSemiAxis));
  // Add pi to get desired solution (third quadrant)
  G4double irisParameterTangentPoint = atan(-irisRSemiAxis/(m*irisZSemiAxis)) + CLHEP::pi; 
 
  // Rounding down to a multiple of 4. This is so that number of points are
  // share equally and consistently between the constituent ellipses.
  noPoints = noPoints - (noPoints % 4);  

  // Vector Definitions:
  // Parametric equation parameters:
  // equatorParameter --> the values of a parameter  for generating the central ellipses
  // irisParameter    --> the values of a parameter  for generating the iris ellipse
  
  //The Cartesian coordinates to be calculated from above parameters:
  //rInnerCoord      --> the values of radius for the inner geometry
  //rOuterCoord      --> the values of radius for the outer geometry
  //zInnerCoord      --> the values of z for the inner geometry
  //zOuterCoord      --> the values of z for the outer goemetry
  
  //Vector declaration:
  std::vector<G4double> equatorParameter;
  std::vector<G4double> irisParameter;
  std::vector<G4double> rInnerCoord;
  std::vector<G4double> rOuterCoord;
  std::vector<G4double> zInnerCoord;
  std::vector<G4double> zOuterCoord; 
    
  // Generate values for the parameters used to define the ellipse shapes
    
  // Vector of parameter values for generating the iris ellipse.  This can then be
  // reflected to find the iris ellipse on the other side.
  for (unsigned int i = 0; i < noPoints/4; i++) //iris allocated one quarter of points.
    {irisParameter.push_back(1.5*CLHEP::pi + (1.5*CLHEP::pi - irisParameterTangentPoint )*i/((noPoints/4)-1));}

  //Vector of parameter values for generating the equator ellipse
  for (unsigned int i = 0; i < noPoints/2; i++) //equator ellipse allocated half of the points
    {
      G4double term1 = CLHEP::pi - equatorParameterTangentPoint;
      G4double term2 = CLHEP::pi - 2*equatorParameterTangentPoint;
      G4double term3 = (noPoints*0.5) - 1;
      G4double value = term1 - term2 * i  / term3;
      equatorParameter.push_back(value);
    }
  
  // Calculate the coordinates of the inner cavity shape, starting from -z through
  // to +z, in 3 sections: -z iris ellipse -> equator ellipse -> +z iris ellipse.
  // Geant requires the coordinates to be done in this (negative to positive)
  
  // -z iris ellipse:
  // noPoints/4 because of the total noPoints, the iris ellipse in -z is allocated a quarter of points.
  for (unsigned int i = 0; i < noPoints/4; i++) 
    {
      // An extra point at the very start for unambigious boolean subtraction later on.
      // Will not actually be part of the final geometry.
      if (i == 0) 
	{
	  //-10mm so boundaries don't lay upon each other when using boolean operation
	  zInnerCoord.push_back(- zi + irisZSemiAxis * cos(irisParameter[i]) - 10*CLHEP::mm); 
	  rInnerCoord.push_back(ri + irisRSemiAxis * sin(irisParameter[i]));
	}
      zInnerCoord.push_back(- zi + irisZSemiAxis * cos(irisParameter[i]));
      rInnerCoord.push_back(ri + irisRSemiAxis * sin(irisParameter[i]));
    }
  
  //equator ellipse (middle):
  for (unsigned int i = 0; i < noPoints/2; i++)  //central ellipse is allocated half the total points.  
    {
      zInnerCoord.push_back(ze + equatorZSemiAxis * cos(equatorParameter[i]));
      rInnerCoord.push_back(re + equatorRSemiAxis * sin(equatorParameter[i]));
    };
  
  //+z iris ellipse:
  for (unsigned int i = 0; i < noPoints/4; i++)  //the +z iris ellipse is allocated noPoints/4.
    {
      zInnerCoord.push_back(zi + irisZSemiAxis * cos(CLHEP::pi - irisParameter[noPoints/4 - 1 - i]));
      rInnerCoord.push_back(ri + irisRSemiAxis * sin(CLHEP::pi - irisParameter[noPoints/4 - 1 - i]));	 
      //If last point, add an extra point purely for unambigious boolean subtraction later.
      if (i == (noPoints/4 - 1))
	{
	  //+10mm so boundaries don't lay upon each other
	  zInnerCoord.push_back(zi + irisZSemiAxis * cos(CLHEP::pi - irisParameter[noPoints/4 - 1 - i]) + 10*CLHEP::mm); 
	  rInnerCoord.push_back(ri + irisRSemiAxis * sin(CLHEP::pi - irisParameter[noPoints/4 - 1 - i]));
	}
    }  
  
  // Generating points for the entire outer cavity shape by using the inner cavity
  // vector of coordinates--- This works by calculating the gradient at each point on the curve.
  // It then finds a unit vector pointing perpendicularly outwards from this point and multiplies
  // by the thickness.  This creates a second vector of points for the outer geometry.
  //starting from 1 and ending at n-1 so that the coordinates used for the boolean subtraction are excluded.
  for (unsigned int i = 1; i < zInnerCoord.size() - 1; i++)
    {
      //calculate the gradient at point i
      G4double gradientAt_i = (rInnerCoord[i-1]-rInnerCoord[i+1])/(zInnerCoord[i-1]-zInnerCoord[i+1]);

      // Since the calculation relies on using the two points either side, and the
      // ends have only one coordinate either side, take the ends as having gradient = 0.
      // So that the vector points vertically upwards.  Otherwise just do the calculation described above.
      if (i == 1 || i == (zInnerCoord.size() - 2))
	{
	  zOuterCoord.push_back(zInnerCoord[i]);
	  rOuterCoord.push_back(rInnerCoord[i] + thickness);
	}
      else
	{
	  zOuterCoord.push_back(- gradientAt_i*thickness * std::pow( (gradientAt_i * gradientAt_i + 1), -0.5) + zInnerCoord[i]);
	  rOuterCoord.push_back(+ thickness * std::pow( (gradientAt_i*gradientAt_i + 1), -0.5) + rInnerCoord[i]);
	}
    }

#ifdef BDSDEBUG
  G4cout << "Now printing the values of (zInnerCoord,rInnerCoord):" << G4endl;
  G4cout << "Length of zInnerCoord = " << zInnerCoord.size() << G4endl;
  G4cout << "Length of rInnerCoord = " << rInnerCoord.size() << G4endl;
  for (G4int i = 0; i < (G4int)zInnerCoord.size(); i++)
    {G4cout << "(" << zInnerCoord[i] << "," << rInnerCoord[i] << ")" << G4endl;};
  
  G4cout << "Now printing the values of (zOuterCoord,rOuterCoord):" << G4endl;
  G4cout << "Length of zOuterCoord = " << zOuterCoord.size() << G4endl;
  G4cout << "Length of rOuterCoord = " << rOuterCoord.size() << G4endl;
  for (G4int i = 0; i < (G4int)zOuterCoord.size(); i++)
    {G4cout << "(" << zOuterCoord[i] << "," << rOuterCoord[i] << ")" << G4endl;};
#endif

  // Array of inner r coordinates.  zeroes ensures the polycone will be solid.
  std::vector<G4double> solidArrayOuter(noPoints, 0.0); 
  // Array of inner r coordinates.  zeroes ensures the polycone will be solid. 1 extra
  // point either side for unambiguous  boolean subtraction.
  std::vector<G4double> solidArrayInner(noPoints+2, 0.0); 

  // Define the inner solid which is to be subtracted from the outer and also used to define the vacuum.
  G4VSolid* innerSolid = new G4Polycone(name + "_inner_solid", //name
					0.0,                   //start angle
					CLHEP::twopi,          //sweep angle
					zInnerCoord.size(),
					zInnerCoord.data(),
					solidArrayInner.data(),
					rInnerCoord.data());   //r coordinates

  // Define the outer solid from which the inner is subtracted.
  G4VSolid* outerSolid = new G4Polycone(name + "_outer_solid", // name
					0.0,                   // start angle
					CLHEP::twopi,          // sweep angle
					zOuterCoord.size(),    // number of corners
					zOuterCoord.data(),    // zOuterCoord array of points
					solidArrayOuter.data(),// array of zeroes for solid shape
					rOuterCoord.data());   // rOuterCoord array of points.
  allSolids.push_back(innerSolid);
  allSolids.push_back(outerSolid);
  
  // Do the subtraction
  G4VSolid* cavitySolid = new G4SubtractionSolid(name + "_cavity_solid",  // name
						 outerSolid,              // solid1
						 innerSolid);             // minus solid2
  
  cavityLV = new G4LogicalVolume(cavitySolid,
				 info->material,
				 name + "_cavity_lv");
  allLogicalVolumes.push_back(cavityLV);
 
  // Delete first and last elements of the InnerCoord Vectors, as these entries.
  // The reason we need to do this is the same vector is also used for making the
  // vacuum filling the cavity, but without the extra points for subtraction.
  rInnerCoord.erase (rInnerCoord.begin());
  rInnerCoord.pop_back();
  zInnerCoord.erase (zInnerCoord.begin());
  zInnerCoord.pop_back();
  
  //Take lengthsafety from all r points and the end z points, to avoid overlap with the cavity.
  for (unsigned int i = 0; i < noPoints; i++)
    {rInnerCoord[i] = rInnerCoord[i] - lengthSafetyLarge;}
  
  zInnerCoord[0]                    = zInnerCoord[0] + 2*lengthSafety;
  G4int lastInd = (G4int)zInnerCoord.size() - 1;
  zInnerCoord[lastInd] = zInnerCoord[lastInd] - 2*lengthSafety; // last element
  
  G4VSolid* vacuumSolid = new G4Polycone(name + "_inner_solid", //name
					 0.0,                   //start angle
					 CLHEP::twopi,          //sweep angle
					 zInnerCoord.size(),
					 zInnerCoord.data(),
					 solidArrayInner.data(),
					 rInnerCoord.data());   //r coordinates
  allSolids.push_back(vacuumSolid);
  
  vacuumLV = new G4LogicalVolume(vacuumSolid,
				 vacuumMaterial,
				 name + "_vacuum_lv");
  allLogicalVolumes.push_back(vacuumLV);

  G4double outerRadius = equatorRadius + thickness + lengthSafetyLarge;
  BuildContainerLogicalVolume(name, chordLength, outerRadius);
  
  SetUserLimits(chordLength, allLogicalVolumes);
  SetVisAttributes("srfcavity");
  PlaceComponents(name);

  BDSExtent ext     = BDSExtent(outerRadius, outerRadius,  chordLength*0.5);
  BDSCavity* result = BuildCavityAndRegisterObjects(ext);  
  return result;
}
