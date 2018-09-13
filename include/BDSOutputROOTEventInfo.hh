/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2018.

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
#ifndef BDSOUTPUTROOTEVENTINFO_H
#define BDSOUTPUTROOTEVENTINFO_H

#include "TROOT.h"
#include "TObject.h"

#include <ctime>

/**
 * @brief Information pertaining to an individual event.
 *
 * @author Stewart Boogert
 */

class BDSOutputROOTEventInfo: public TObject
{  
public:
  BDSOutputROOTEventInfo();

  virtual ~BDSOutputROOTEventInfo();
  void Flush();
  
  time_t startTime; ///< Time stamp at start of event.
  time_t stopTime;  ///< Time stamp at end of event.
  float  duration;  ///< Number of seconds event took to complete simulation (not writing out).
  std::string seedStateAtStart; ///< Seed state at the start of the event.
  int    index;             ///< Number of this event or run.
  bool   aborted;           ///< Whether the event was aborted or not.
  bool   primaryHitMachine; ///< Wether the primary particle hit the accelerator or not.
  double memoryUsageMb;     ///< memoryUsage (rusage.ru_maxrss)
 
  ClassDef(BDSOutputROOTEventInfo,4);
};

#endif
