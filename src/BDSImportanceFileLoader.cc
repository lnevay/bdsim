/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2019.

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

#include "BDSDebug.hh"
#include "BDSImportanceFileLoader.hh"

#include "globals.hh" // geant4 types / globals
#include "G4String.hh"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <exception>
#include <map>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#ifdef USE_GZSTREAM
#include "gzstream.h"
#endif

BDSImportanceFileLoader::BDSImportanceFileLoader()
{;}

BDSImportanceFileLoader::~BDSImportanceFileLoader()
{;}

std::map<G4String, G4double> BDSImportanceFileLoader::Load(G4String fileName)
{
  std::ifstream file;
  std::vector<G4String> volumes;
  std::vector<G4double> importanceValues;

  file.open(fileName);

  // test if file is valid
#ifdef USE_GZSTREAM
  bool validFile = file.rdbuf()->is_open();
#else
  bool validFile = file.is_open();
#endif

  if (!validFile)
    {
      G4cerr << "Invalid file name or no such file named \"" << fileName << "\"" << G4endl;
      exit(1);
    }
  else
    {G4cout << "BDSIM importance file - loading \"" << fileName << "\"" << G4endl;}

  std::string line;
  std::map<G4String, G4double> importance;

  while (std::getline(file, line))
    { // read a line only if it's not a blank one

      std::istringstream liness(line);
      std::string volume;
      G4double importanceValue;

      // Skip a line if it's only whitespace
      if (std::all_of(line.begin(), line.end(), isspace))
        {continue;}

      liness >> volume >> importanceValue;

      // importance world should be GDML import, modify PV name accordingly.
      G4String fullVolume = "importanceWorld_PREPEND" + volume + "_pv";
      volumes.push_back(fullVolume);
      importanceValues.push_back(importanceValue);
      importance[fullVolume] = importanceValue;
    }

  file.close();

  return importance;
}
