/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2023.

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
#ifndef BDSWARNING_H
#define BDSWARNING_H

#include "G4String.hh"

/**
 * @brief General warning that prints a message and hangs briefly.
 * 
 * @author Laurie Nevay
 */

namespace BDS
{
  void Warning(const G4String& message);

  /// This will only warn if activeThread is true in the multithreaded case.
  void Warning(const G4String& methodName,
	       const G4String& message,
	       bool activeThread = true);
}

#endif
