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
#ifndef HISTOGRAMACCUMULATORV_H
#define HISTOGRAMACCUMULATORV_H

#include <vector>

class BDSOutputROOTEventHistograms;
class HistogramAccumulator;

template <class T>
class HistogramAccumulatorV
{
public:
  HistogramAccumulatorV(BDSOutputROOTEventHistograms* tempHistograms);
  ~HistogramAccumulatorV();
  void Accumulate(BDSOutputROOTEventHistograms* oneEventHistograms);
  BDSOutputROOTEventHistograms* Terminate();
  
private:
  std::vector<T*> accumulators1d;
  std::vector<T*> accumulators2d;
  std::vector<T*> accumulators3d;
};

#endif
