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
#include "BDSApertureType.hh"
#include "BDSDebug.hh"
#include "BDSException.hh"

#include "globals.hh"

#include <map>

// dictionary for BDSApertureType
template<>
std::map<BDSApertureType, std::string>* BDSApertureType::dictionary =
  new std::map<BDSApertureType, std::string> ({
   {BDSApertureType::circle,      "circular"},
   {BDSApertureType::ellipse,     "elliptical"},
   {BDSApertureType::rectangle,   "rectangular"},
   {BDSApertureType::rectcircle,  "rectcircle"},
   {BDSApertureType::rectellipse, "rectellipse"},
   {BDSApertureType::racetrack,   "racetrack"},
   {BDSApertureType::octagonal,   "octagonal"},
   {BDSApertureType::clicpcl,     "clicpcl"}
});	

BDSApertureType BDS::DetermineApertureType(G4String apertureType)
{
  std::map<G4String, BDSApertureType> types;
  types["circular"]       = BDSApertureType::circle;
  types["elliptical"]     = BDSApertureType::ellipse;
  types["rectangular"]    = BDSApertureType::rectangle;
  types["rectcircle"]     = BDSApertureType::rectcircle;
  types["rectellipse"]    = BDSApertureType::rectellipse;
  types["racetrack"]      = BDSApertureType::racetrack;
  types["octagonal"]      = BDSApertureType::octagonal;
  types["clicpcl"]        = BDSApertureType::clicpcl;

  // alternatives
  types["lhc"]            = BDSApertureType::rectcircle;
  
  apertureType.toLower();

  auto result = types.find(apertureType);
  if (result == types.end())
    {
      // it's not a valid key
      G4cerr << __METHOD_NAME__ << apertureType << " is not a valid apertureType" << G4endl;

      G4cout << "Available geometry types are:" << G4endl;
      for (auto it : types)
	{G4cout << "\"" << it.first << "\"" << G4endl;}
      throw BDSException(__METHOD_NAME__, "");
    }
  
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << "determined aperture type to be " << result->second  << G4endl;
#endif
  return result->second;
}
