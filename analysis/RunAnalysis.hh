#ifndef RUNANALYSIS_H
#define RUNANALYSIS_H

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"

#include "BDSOutputROOTEventHistograms.hh"

#include "Run.hh"
#include "Analysis.hh"

/**
 * @brief Analysis of a run.
 *
 * @author Stuart Walker.
 */

class RunAnalysis : public Analysis
{
public:
  RunAnalysis();
  RunAnalysis(Run *r, TChain *c, bool debug = false);
  virtual ~RunAnalysis();

  virtual void Process();

protected:
  Run* run;

  ClassDef(RunAnalysis,1);
};

#endif
