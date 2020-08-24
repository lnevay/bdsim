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
#include <iostream>

#include "TRKElement.hh"
#include "TRKElementLine.hh"

TRKElementLine::TRKElementLine(std::string nameIn) :
  TRKElement(nameIn, 0.0, nullptr, nullptr)
{;}

void TRKElementLine::AddElement(TRKElement* e) {
  line.push_back(e);
  length += e->GetLength();
}

/// output stream
std::ostream& operator<< (std::ostream &out, const TRKElementLine &line) {
  auto elIter = line.begin();
  auto elIterEnd = line.end();
  for (;elIter!=elIterEnd; ++elIter) {
    out << **elIter << std::endl;
  }
  return out;
}
