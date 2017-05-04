#ifndef EVENTANALYSISORBIT_H
#define EVENTANALYSISORBIT_H

#include "EventAnalysis.hh"

class TFile;

class EventAnalysisOrbit: public EventAnalysis
{
public:
  using EventAnalysis::EventAnalysis;

  virtual ~EventAnalysisOrbit(){;}

  void Clear();
  void ExtractOrbit(int index);
  void WriteOrbit(TFile* f);

  std::vector<double> ss;
  std::vector<double> x;
  std::vector<double> xp;
  std::vector<double> y;
  std::vector<double> yp;
};

#endif
