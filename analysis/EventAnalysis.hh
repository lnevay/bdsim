#ifndef EVENTANALYSIS_H
#define EVENTANALYSIS_H

#include <vector>

#include "Analysis.hh"
#include "TROOT.h"

class Event;
class SamplerAnalysis;
class TChain;
class TFile;

/**
 * @brief Event level analysis.
 *
 * @author Stewart Boogert
 */

class EventAnalysis: public Analysis
{
public:
  EventAnalysis();
  EventAnalysis(Event*  eventIn,
		TChain* chain,
		bool    processSamplersIn   = false,
		bool    debug               = false,
		double  printModuloFraction = 0.01,
		bool    emittanceOnTheFlyIn = false);
  virtual ~EventAnalysis();

  void SetPrintModuloFraction(double fraction);

  virtual void Process();

  /// Terminate each individual sampler analysis and append optical functions.
  virtual void Terminate();

  /// Write analysis including optical functions to an output file.
  virtual void Write(TFile *outputFileName);

protected:
  Event* event;
  std::vector<SamplerAnalysis*> samplerAnalyses;
  std::vector<std::vector<std::vector<double>>> opticalFunctions; ///< optical functions from all samplers
  ClassDef(EventAnalysis,1);

private:
  /// Initialise each sampler analysis object in samplerAnalysis.
  void Initialise();
  void ProcessSamplers(bool firstTime = false);

  int  printModulo;     ///< Cache of print modulo fraction
  bool processSamplers; ///< Whether to process samplers.
  bool emittanceOnTheFly; ///< Whether to calculate emittance fresh at each sampler.
};

#endif
