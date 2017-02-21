#include "Analysis.hh"
#include "Config.hh"
#include "rebdsim.hh"

#include "TChain.h"

Analysis::Analysis(std::string treeNameIn,
		   TChain*     chainIn,
		   std::string mergedHistogramNameIn,
		   bool        debugIn):
  treeName(treeNameIn),
  chain(chainIn),
  mergedHistogramName(mergedHistogramNameIn),
  histoSum(nullptr),
  debug(debugIn)
{;}

Analysis::~Analysis()
{
  delete chain;
  delete histoSum;
}

void Analysis::Execute()
{
  std::cout << "Analysis on \"" << treeName << "\" beginning" << std::endl;
  Process();
  SimpleHistograms();
  Terminate();
  std::cout << "Analysis on \"" << treeName << "\" complete" << std::endl;
}

void Analysis::SimpleHistograms()
{
  if(debug)
    {std::cout << __METHOD_NAME__ << std::endl;}

  // loop over histogram specifications and fill
  auto c = Config::Instance();
  if (c)
    {
      auto hd = Config::Instance()->GetHistoDefs();  // histogram definitions
      for(auto i : hd)
	{
	  if (i["treeName"] == treeName)
	    {
	      std::string histName = i["histName"];
	      std::cout << "Filling histogram \"" << histName << "\"" << std::endl;
	      FillHistogram(i["treeName"].data(), histName, i["nbins"], i["binning"], i["plot"], i["select"]);
	    }
	}
    }
}

void Analysis::Terminate()
{
  if (histoSum)
    {histoSum->Terminate();}
}

void Analysis::FillHistogram(std::string treeName, std::string histoName,
			     std::string nbins,    std::string binning,
			     std::string plot,     std::string selection)
{
  TH1::AddDirectory(kTRUE);
  TH2::AddDirectory(kTRUE);
  TH3::AddDirectory(kTRUE);
  
  if(debug)
    {std::cout << __METHOD_NAME__ << std::endl;}
  double xlow=0.0, xhigh=0.0;
  double ylow=0.0, yhigh=0.0;
  double zlow=0.0, zhigh=0.0;
  int ndim = Config::Dimension(nbins);
  int nxbin=0, nybin=0, nzbin=0;

  if(ndim == 1)
  {
    nxbin = Config::NBins(nbins,1);
    Config::Binning(binning,1,xlow,xhigh);
  }
  else if(ndim == 2)
  {
    nxbin = Config::NBins(nbins,1);
    nybin = Config::NBins(nbins,2);
    Config::Binning(binning,1,xlow,xhigh);
    Config::Binning(binning,2,ylow,yhigh);
  }
  else if(ndim == 3)
    {
      nxbin = Config::NBins(nbins,1);
      nybin = Config::NBins(nbins,2);
      nzbin = Config::NBins(nbins,3);
      Config::Binning(binning,1,xlow,xhigh);
      Config::Binning(binning,2,ylow,yhigh);
      Config::Binning(binning,3,zlow,zhigh);
    }
  else
    {std::cerr << "Invalid number of dimensions: \"" << ndim << "\"" << std::endl; exit(1);}

  if (!chain)
  {std::cout << __METHOD_NAME__ << "Error no tree found by name: " << treeName << std::endl; exit(1);}

  if(debug)
    {std::cout << __METHOD_NAME__ << treeName << " " << histoName << " " << plot << std::endl;}

  auto pltSav = histoName;
  auto pltCmd = plot+" >> "+pltSav;

  // create histograms
  if(ndim == 1)
    {new TH1D(pltSav.c_str(),pltSav.c_str(), nxbin, xlow, xhigh);}
  else if(ndim == 2)
    {new TH2D(pltSav.c_str(),pltSav.c_str(), nxbin, xlow, xhigh, nybin, ylow, yhigh);}
  else if(ndim == 3)
    {new TH3D(pltSav.c_str(),pltSav.c_str(), nxbin, xlow, xhigh, nybin, ylow, yhigh, nzbin, zlow, zhigh);}

  chain->Draw(pltCmd.c_str(),selection.c_str(),"goff");

  if(ndim == 1)
  {
    auto h = (TH1*)gDirectory->Get(pltSav.c_str());
    this->histogramNames.push_back(pltSav);
    this->histograms1D[pltSav] = h;
    // std::cout << h << std::endl;
  }
  else if(ndim == 2)
  {
    auto h = (TH2*)gDirectory->Get(pltSav.c_str());
    this->histogramNames.push_back(pltSav);
    this->histograms2D[pltSav] = h;
    // std::cout << h << std::endl;
  }
  else if (ndim == 3)
    {
      auto h = (TH3*)gDirectory->Get(pltSav.c_str());
      histogramNames.push_back(pltSav);
      histograms3D[pltSav] = h;
    }


  if(debug)
  {
    std::cout << __METHOD_NAME__ << "`" << pltSav << "'  `" << pltCmd << "' " << gDirectory->Get(pltSav.c_str()) << std::endl;
  }
}

void Analysis::Write(TFile* outputFile)
{
  //treeName typically has a "." at the end, deleting it here:
  std::string cleanedName = treeName.erase(treeName.size() - 1);
  std::string outputDirName = cleanedName + std::string("Histograms");
  TDirectory *rebdsimDir = outputFile->mkdir(outputDirName.c_str());
  rebdsimDir->cd();
  for(auto h : histograms1D)
    {h.second->Write();}
  for(auto h : histograms2D)
    {h.second->Write();}
  for (auto h : histograms3D)
    {h.second->Write();}
  outputFile->cd("/");
  
  // Merged Histograms for this analysis instance (could be run, event etc)
  if (histoSum)
    {
      std::cout << "Merging histograms from \"" << treeName << "\" analysis" << std::endl;
      TDirectory* bdsimDir = outputFile->mkdir(mergedHistogramName.c_str());
      bdsimDir->cd();
      histoSum->Write(outputFile);
    }
}
