/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2024.

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
#include "BDSBeamPipeFactory.hh"
#include "BDSBeamPipeFactoryBase.hh"
#include "BDSBeamPipeFactoryGeneral.hh"
#include "BDSBeamPipeFactoryLHCDetailed.hh"
#include "BDSBeamPipeInfo2.hh"
#include "BDSBeamPipeType.hh"
#include "BDSDebug.hh"
#include "BDSException.hh"

#include "G4String.hh"
#include "G4Types.hh"

BDSBeamPipeFactory::BDSBeamPipeFactory()
{
  general     = new BDSBeamPipeFactoryGeneral();
  lhcdetailed = new BDSBeamPipeFactoryLHCDetailed();
}

BDSBeamPipeFactory::~BDSBeamPipeFactory()
{
  delete general;
  delete lhcdetailed;
}

BDSBeamPipeFactoryBase* BDSBeamPipeFactory::GetAppropriateFactory(BDSBeamPipeType type)
{
  BDSBeamPipeFactoryBase* result;
  switch(type.underlying())
    {
    case BDSBeamPipeType::circular:
    case BDSBeamPipeType::elliptical:
    case BDSBeamPipeType::rectangular:
    case BDSBeamPipeType::lhc:
    case BDSBeamPipeType::rectellipse:
    case BDSBeamPipeType::racetrack:
    case BDSBeamPipeType::octagonal:
    case BDSBeamPipeType::circularvacuum:
    case BDSBeamPipeType::clicpcl:
    case BDSBeamPipeType::pointsfile:
    case BDSBeamPipeType::rhombus:
      {result = general; break;}
    case BDSBeamPipeType::lhcdetailed:
      {result = lhcdetailed; break;}
    default:
#ifdef BDSDEBUG
      G4cout << __METHOD_NAME__ << "unknown type \"" << type << "\" - circular beampipe factory by default" << G4endl;
#endif
      return circular;
      break;
    }
}

BDSBeamPipe* BDSBeamPipeFactory::CreateBeamPipeForVacuumIntersection(const G4String&  name,
								     G4double         length,
								     BDSBeamPipeInfo* bpi)
{
  BDSBeamPipeInfo copy = BDSBeamPipeInfo(*bpi);
  copy.ShrinkBy(BDSGlobalConstants::Instance()->LengthSafetyLarge());
  return CreateBeamPipe(BDSBeamPipeType::circularvacuum,
			name,
			length,
			bpi->aper1,
			bpi->aper2,
			bpi->aper3,
			bpi->aper4,
			bpi->vacuumMaterial,
			bpi->beamPipeThickness,
			bpi->beamPipeMaterial);
}

BDSBeamPipe* BDSBeamPipeFactory::CreateBeamPipe(const G4String&  name,
						G4double         length,
						BDSBeamPipeInfo* bpi)
{
  if ((bpi->inputFaceNormal.z() > -1) || (bpi->outputFaceNormal.z() < 1))
    {
      return CreateBeamPipe(bpi->beamPipeType,
			    name,
			    length,
			    bpi->inputFaceNormal,
			    bpi->outputFaceNormal,
			    bpi->aper1,
			    bpi->aper2,
			    bpi->aper3,
			    bpi->aper4,
			    bpi->vacuumMaterial,
			    bpi->beamPipeThickness,
			    bpi->beamPipeMaterial,
			    bpi->pointsFileName,
			    bpi->pointsUnit);
    }
  else
    {
      return CreateBeamPipe(bpi->beamPipeType,
			    name,
			    length,
			    bpi->aper1,
			    bpi->aper2,
			    bpi->aper3,
			    bpi->aper4,
			    bpi->vacuumMaterial,
			    bpi->beamPipeThickness,
			    bpi->beamPipeMaterial,
			    bpi->pointsFileName,
			    bpi->pointsUnit);
    }
  return result;
}

BDSBeamPipe* BDSBeamPipeFactory::CreateBeamPipe(const G4String&  name,
                                                G4double         length,
                                                BDSBeamPipeInfo2* bpi)
{
  BDSBeamPipeFactoryBase* factory = GetAppropriateFactory(bpi->beamPipeType);
  return factory->CreateBeamPipe(name, length, bpi);
}
