#ifndef DATALOADER_H
#define DATALOADER_H

#include "TROOT.h"

#include "RebdsimTypes.hh"

#include <map>
#include <string>
#include <vector>

class Event;
class Options;
class Model;
class Run;
class TChain;

/**
 * @brief Loader for a ROOT file using classes used to generate the file.
 *
 * @author Stewart Boogert
 */

class DataLoader
{
public:
  DataLoader(std::string fileName,
	     bool        debugIn           = false,
	     bool        processSamplersIn = false,
	     bool        allBranchesOn     = true,
	     const RBDS::BranchMap* branchesToTurnOn = nullptr);
  virtual ~DataLoader();

  /// Create an instance of each class in the file to be overlaid by loading
  /// the ROOT file.
  void CommonCtor(std::string fileName);

  /// Build up the input file list.
  void BuildInputFileList(std::string inputPath);

  /// Open the first file in the file list and map the trees in it.
  void BuildTreeNameList();

  /// Inspect the first file (leaving it open...) and build a list of samplers
  /// by checking whether they have 'Sampler' in the name.
  void BuildEventBranchNameList();

  /// Create a tree for each sampler and add all the files to it.
  void ChainTrees();

  /// Map each chain to the member instance of each storage class in this class.
  void SetBranchAddress(bool allOn = true,
                        const RBDS::BranchMap* bToTurnOn = nullptr);

  /// @{ Accessor
  std::vector<std::string>   GetTreeNames()    {return treeNames;};
  std::vector<std::string>   GetBranchNames()  {return branchNames;}
  std::vector<std::string>   GetSamplerNames() {return samplerNames;}
  Options*                   GetOptions()      {return opt;}
  Model*                     GetModel()        {return mod;}
  Event*                     GetEvent()        {return evt;}
  Run*                       GetRun()          {return run;};
  TChain*                    GetOptionsTree()  {return optChain;}
  TChain*                    GetModelTree()    {return modChain;}
  TChain*                    GetEventTree()    {return evtChain;}
  TChain*                    GetRunTree()      {return runChain;}
  /// @}

private:
  DataLoader();
  
  bool debug;
  bool processSamplers;
  bool allBranchesOn;
  const RBDS::BranchMap* branchesToTurnOn;
  
  Options                      *opt;
  Model                        *mod;
  Event                        *evt;
  Run                          *run;

  std::vector<std::string>      fileNames;
  std::vector<std::string>      safeFileNames;

  std::vector<std::string>      treeNames;
  std::vector<std::string>      branchNames;  // non-sampler branch names
  std::vector<std::string>      samplerNames; // sampler branch names
  std::map<std::string, int>    samplerNameMap;

  TChain *optChain;
  TChain *modChain;
  TChain *evtChain;
  TChain *runChain;

  ClassDef(DataLoader,1);
};

#endif
