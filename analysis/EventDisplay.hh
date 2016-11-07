#ifndef EVENTDISPLAY_H
#define EVENTDISPLAY_H

#include "TROOT.h"

#include "DataLoader.hh"
#include "Event.hh"
#include "TChain.h"

/**
 * @brief Event viewer using ROOT EVE framework.
 *
 * @author Stewart Boogert
 */

class EventDisplay
{
public:
  /// Singleton accessor.
  static EventDisplay* Instance();

  /// Singleton accessor / constructor w.r.t. a file.
  static EventDisplay* Instance(TString geoFileName,
				TString dataFileName);

  virtual ~EventDisplay();

  /// Load an entry from the model tree.
  void LoadModel(int iMod);

  /// Load an entry from the options tree.
  void LoadOptions(int iOpt);

  /// Load an entry in the event tree.
  void LoadData(int iEvt);

  /// Clear a currently displayed event.
  void ClearEvent();

  /// Draw all things loaded at this point.
  void Draw();

  /// Draw the geometry.
  void DrawModel();

  /// Draw energy loss deposits.
  void DrawElossHits();

  /// Draw tunnel energy deposits.
  void DrawTunnelHits();

  /// Draw hits on a sampler plane.
  void DrawSamplers();

  /// Draw loaded trajectories.
  void DrawTrajectories();

private:
  /// Private default constructor as singleton pattern.
  EventDisplay();

  /// Supplied constructor for loading a file.
  EventDisplay(TString geoFileNameIn,
	       TString dataFileNameIn);

  /// Load the geometry from gdml file.
  void LoadGeometry();

  TString geoFileName    = "";
  TString dataFileName   = "";

  DataLoader *dataLoader             = nullptr;
  Event      *event                  = nullptr;
  TChain     *eventTree              = nullptr;
  Options    *options                = nullptr;
  TChain     *optionsTree            = nullptr;
  Model      *model                  = nullptr;
  TChain     *modelTree              = nullptr;

  /// Singleton instance.
  static EventDisplay *_instance;

  ClassDef(EventDisplay,1);
};

#endif
