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
#ifndef BDSOUTPUTROOTEVENTHISTOGRAMS_H
#define BDSOUTPUTROOTEVENTHISTOGRAMS_H

#include <string>
#include <vector>

#include <boost/format.hpp>
#include <boost/histogram.hpp>

#ifndef __ROOTBUILD__
#include "globals.hh"
#endif

#include "Rtypes.h"
#include "TObject.h"

class TH1D;
class TH2D;
class TH3D;


/**
 * @brief Holder for a set of histograms to be stored.
 *
 * @author Stewart Boogert
 */

class BDSOutputROOTEventHistograms: public TObject
{
public:
  BDSOutputROOTEventHistograms();
  BDSOutputROOTEventHistograms(const BDSOutputROOTEventHistograms &h);
  BDSOutputROOTEventHistograms& operator=(const BDSOutputROOTEventHistograms&) = delete;
  BDSOutputROOTEventHistograms(std::vector<TH1D*>& histogram1DIn,
			       std::vector<TH2D*>& histogram2DIn,
			       std::vector<TH3D*>& histogram3DIn,
			       std::vector<boost::histogram::histogram<int>> histograms4DIn);
  virtual ~BDSOutputROOTEventHistograms();

  /// Interface function to create a 1D histogram using only standard types.
  int Create1DHistogramSTD(std::string name, std::string title,
			   int nbins, double xmin, double xmax);

#ifndef __ROOTBUILD__
  G4int Create1DHistogram(G4String name, G4String title,
                          G4int nbins, G4double xmin, G4double xmax);
  G4int Create1DHistogram(G4String name, G4String title,
                          std::vector<double>& edges);
  G4int Create2DHistogram(G4String name, G4String title,
                          G4int nxbins, G4double xmin, G4double xmax,
                          G4int nybins, G4double ymin, G4double ymax);
  G4int Create2DHistogram(G4String name, G4String title,
                          std::vector<double>& xedges,
                          std::vector<double>& yedges);
  G4int Create3DHistogram(G4String name, G4String title,
			  G4int nxbins, G4double xmin, G4double xmax,
			  G4int nybins, G4double ymin, G4double ymax,
			  G4int nzbins, G4double zmin, G4double zmax);
  G4int Create3DHistogram(G4String name, G4String title,
			  std::vector<double>& xedges,
			  std::vector<double>& yedges,
			  std::vector<double>& zedges);
  G4int Create4DHistogram(G4String name, G4String title,
                          G4int nxbins, G4double xmin, G4double xmax,
                          G4int nybins, G4double ymin, G4double ymax,
                          G4int nzbins, G4double zmin, G4double zmax,
                          G4int nebins, G4double emin, G4double emax);

  void Fill1DHistogram(G4int histoId, G4double value, G4double weight = 1.0);
  void Fill2DHistogram(G4int histoId, G4double xValue, G4double yValue, G4double weight = 1.0);
  void Fill3DHistogram(G4int histoId, G4double xValue, G4double yValue, G4double zValue, G4double weight = 1.0);
  void Fill4DHistogram(G4int histoId, G4double xValue, G4double yValue, G4double zvalue, G4double eValue);
  
  /// Set the value of a bin by (ROOT!!) global bin index. Note the TH3 function should
  /// be used to get ROOT's idea of a global bin index.
  void Set3DHistogramBinContent(G4int    histoId,
				G4int    globalBinID,
				G4double value);
  
  /// Add the values from one supplied 3D histogram to another. Uses TH3-Add().
  void AccumulateHistogram3D(G4int histoId,
			     TH3D* otherHistogram);
#endif
  /// Flush the contents.
  virtual void Flush();
  
  /// Copy (using the TH->Clone) method from another instance.
  void Fill(const BDSOutputROOTEventHistograms* rhs);

  /// Copy (without using the TH->Clone) method from another instance. (Quicker).
  void FillSimple(const BDSOutputROOTEventHistograms* rhs);

  /// @{ Accessors.
  std::vector<TH1D*>& Get1DHistograms() {return histograms1D;}
  std::vector<TH2D*>& Get2DHistograms() {return histograms2D;}
  std::vector<TH3D*>& Get3DHistograms() {return histograms3D;}
  TH1D* Get1DHistogram(int iHisto) const {return histograms1D[iHisto];}
  TH2D* Get2DHistogram(int iHisto) const {return histograms2D[iHisto];}
  TH3D* Get3DHistogram(int iHisto) const {return histograms3D[iHisto];}
  /// @}

private:
  std::vector<TH1D*> histograms1D;
  std::vector<TH2D*> histograms2D;
  std::vector<TH3D*> histograms3D;
  std::vector<boost::histogram::histogram<boost::histogram::axis::regular<double>>> histograms4D;

  ClassDef(BDSOutputROOTEventHistograms,3);
};

#endif
