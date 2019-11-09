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
#include <G4SystemOfUnits.hh>
#include "scorer.h"

using namespace GMAD;

Scorer::Scorer()
{
  clear();
  PublishMembers();
}

void Scorer::clear()
{
  name                 = "";
  type                 = "";
  particleName         = "";
  particlePDGID        = 0;
  minimumEnergy        = 0;
  maximumEnergy        = 1e8; // 100 PeV - well beyond G4 physics for now
  minimumTime          = 0;
  maximumTime          = 0;
  conversionFactorFile = "";
  materialToInclude    = "";
  materialToExclude    = "";
  scoreWorldVolumeOnly = false;
}

void Scorer::PublishMembers()
{
  publish("name",                 &Scorer::name);
  publish("type",                 &Scorer::type);
  publish("particleName",         &Scorer::particleName);
  publish("particlePDGID",        &Scorer::particlePDGID);
  publish("minimumEnergy",        &Scorer::minimumEnergy);
  publish("maximumEnergy",        &Scorer::maximumEnergy);
  publish("minimumTime",          &Scorer::minimumTime);
  publish("maximumTime",          &Scorer::maximumTime);
  publish("conversionFactorFile", &Scorer::conversionFactorFile);
  publish("conversionFactorPath", &Scorer::conversionFactorPath);
  publish("materialToInclude",    &Scorer::materialToInclude);
  publish("materialToExclude",    &Scorer::materialToExclude);
  publish("scoreWorldVolumeOnly", &Scorer::scoreWorldVolumeOnly);
}

void Scorer::print()const
{
  std::cout << "field: "
	    << "name "                  << name                 << std::endl
	    << "type "                  << type                 << std::endl
	    << "particleName "          << particleName         << std::endl
	    << "particlePDGID "         << particlePDGID        << std::endl
	    << "minimumEnergy "         << minimumEnergy        << std::endl
	    << "maximumEnergy "         << maximumEnergy        << std::endl
	    << "minimumTime "           << minimumTime          << std::endl
	    << "maximumTime "           << maximumTime          << std::endl
	    << "conversionFactorFile "  << conversionFactorFile << std::endl
	    << "conversionFactorPath "  << conversionFactorPath << std::endl
	    << "materialToInclude "     << materialToInclude    << std::endl
	    << "materialToExclude "     << materialToExclude    << std::endl
	    << "scoreWorldVolumeOnly "  << scoreWorldVolumeOnly << std::endl;
}
