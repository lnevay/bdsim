#ifndef EVENTANALYSIS_H
#define EVENTANALYSIS_H

#include <vector>
#include <map>

#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"

#include "Event.hh"
#include "SamplerAnalysis.hh"
#include "HistogramMerge.hh"

/**
 * @brief Event level analysis.
 * 
 * @author Stewart Boogert
 */

class EventAnalysis
{
public:
  EventAnalysis();
  EventAnalysis(Event* eventIn, TChain* chainIn);
  virtual ~EventAnalysis();

  /// Initialise each sampler analysis object in samplerAnalysis.
  void Initialise();
  void Process();
  void ProcessSamplers();

  /// Terminate each individual sampler analysis and append optical functions.
  void Terminate();

  /// Process histogram definitions from configuration instance.
  void SimpleHistograms();

  /// Create an individual histogram based on a string selection etc.
  void FillHistogram(std::string treeName, std::string histoName,
                     std::string nbins,    std::string binning,
                     std::string plot,     std::string selection);

  /// Write analysis including optical functions to an output file.
  void Write(TFile *outputFileName);

protected:
  Event  *event;
  TChain *chain;
  std::vector<SamplerAnalysis*> samplerAnalyses;
  std::vector<std::string>      histogramNames;  ///< Rebdsim generated histogram names
  std::map<std::string, TH1*>   histograms1D;    ///< Rebdsim 1d histogram
  std::map<std::string, TH2*>   histograms2D;    ///< Rebdsim 2d histograms
  HistogramMerge               *histoSum;        ///< Bdsim histograms
  std::vector<std::vector<std::vector<double>>> opticalFunctions; ///< optical functions from all samplers
  ClassDef(EventAnalysis,1);
};

#endif
