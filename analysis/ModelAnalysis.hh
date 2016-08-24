#ifndef MODELANALYSIS_H
#define MODELANALYSIS_H

#include "TROOT.h"

#include "Analysis.hh"
#include "Model.hh"

class TChain;

class ModelAnalysis : public Analysis
{
public:
  ModelAnalysis();
  ModelAnalysis(Model *modelIn, TChain *chainIn);
  virtual ~ModelAnalysis(){};

  virtual void Process(){};
  virtual void Terminate(){};

protected:
  Model    *model;
  ClassDef(ModelAnalysis,1);
};


#endif
