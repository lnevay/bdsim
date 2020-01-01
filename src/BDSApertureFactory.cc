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
#include "BDSAperture.hh"
#include "BDSApertureCircular.hh"
#include "BDSApertureFactory.hh"
#include "BDSApertureType.hh"

#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4Types.hh"

#include "G4Box.hh"
#include "G4CutTubs.hh"
#include "G4IntersectionSolid.hh"
#include "G4Tubs.hh"
#include "G4VSolid.hh"

#include "CLHEP/Units/SystemOfUnits.h"


BDSApertureFactory::BDSApertureFactory():
  intersectionRadiusRatio(1.3),
  productName(""),
  productLength(0),
  productApertureIn(nullptr),
  productApertureOut(nullptr),
  productLengthExtra(0),
  angledFaces(false)
{
  /*
  specialisations = {
  std::make_pair(BDSApertureType::*/
}

BDSApertureFactory::~BDSApertureFactory()
{;}

G4VSolid* BDSApertureFactory::CreateAperture(const G4String&    name,
					     G4double           length,
					     const BDSAperture* apertureIn,
					     const BDSAperture* apertureOut,
					     G4double           lengthExtraForBoolean)
{
  G4bool variedAperture = (G4bool)apertureOut; // ie valid pointer for shape out.

  productName        = name;
  productLength      = length;
  productApertureIn  = apertureIn;
  productApertureOut = apertureOut;
  productLengthExtra = lengthExtraForBoolean;
  productNormalIn    = G4ThreeVector();
  productNormalOut   = G4ThreeVector();
  angledFaces        = false;

  if (variedAperture)
    {return CreateDifferentEnds();}

  // else create non-tapering piece of geometry
  G4VSolid* product = nullptr;
  switch (apertureIn->apertureType.underlying())
    {
    case BDSApertureType::circular:
      {product = CreateCircular();    break;}
    case BDSApertureType::rectangular:
      {product = CreateRectangular(); break;}
    case BDSApertureType::elliptical:
      {product = CreateElliptical();  break;}
    case BDSApertureType::lhc:
      {product = CreateLHC();         break;}
    case BDSApertureType::rectellipse:
      {product = CreateRectEllipse(); break;}
    case BDSApertureType::racetrack:
      {product = CreateRaceTrack();   break;}
    case BDSApertureType::octagonal:
      {product = CreateOctagonal();   break;}
    case BDSApertureType::clicpcl:
      {product = CreateClicPCL();     break;}
    default:
      {break;}
    }
  return product;
}

G4VSolid* BDSApertureFactory::CreateCircular() const
{
  const BDSApertureCircular* ap = dynamic_cast<const BDSApertureCircular*>(productApertureIn);
  if (!ap)
    {return nullptr;}
  if (!angledFaces)
    {
      G4VSolid* product = new G4Tubs(productName,
				     0,
				     ap->radius,
				     0.5 * productLength + productLengthExtra,
				     0,
				     CLHEP::twopi);
      return product;
    }
  else
    {
      G4VSolid* product = new G4CutTubs(productName,
					0,
					ap->radius,
					0.5 * productLength + productLengthExtra,
					0,
					CLHEP::twopi,
					productNormalIn,
					productNormalOut);
      return product;
    }
}

G4VSolid* BDSApertureFactory::CreateRectangular() const
{return nullptr;}
/*
G4VSolid* BDSApertureFactory::CreateRectangular() const
{
  const BDSApertureRectangular* ap = dynamic_cast<const BDSApertureRectangular*>(productApertureIn);
  if (!ap)
    {return nullptr;}
  if (!angledFaces)
    {
      G4VSolid* product = new G4Box(productName,
				    ap->a,
				    ap->b,
				    0.5 * productLength + productLengthExtra);
      return product;
    }
  else
    {
      G4VSolid* box = new G4Box(productName + "_square",
				ap->a,
				ap->b,
				productLength + productLengthExtra); // factor 2 here
      G4double maxRadius = ap->RadiusToEncompass();
      G4double intersectionRadius = intersectionRadiusRatio * maxRadius;
      G4VSolid* cut = new G4CutTubs(productName + "_angled",
				    0,
				    intersectionRadius,
				    0.5 * productLength + productLengthExtra,
				    0,
				    CLHEP::twopi,
				    productNormalIn,
				    productNormalOut);
      G4VSolid* product = new G4IntersectionSolid(productName, box, cut);
      return product;
    }
}
*/

G4VSolid* BDSApertureFactory::CreateElliptical() const
{return CreateCircular();}
G4VSolid* BDSApertureFactory::CreateLHC() const
{return CreateCircular();}
G4VSolid* BDSApertureFactory::CreateRectEllipse() const
{return CreateCircular();}
G4VSolid* BDSApertureFactory::CreateRaceTrack() const
{return CreateCircular();}
G4VSolid* BDSApertureFactory::CreateOctagonal() const
{return CreateCircular();}
G4VSolid* BDSApertureFactory::CreateClicPCL() const
{return CreateCircular();}

G4VSolid* BDSApertureFactory::CreateDifferentEnds() const
{return CreateCircular();}
