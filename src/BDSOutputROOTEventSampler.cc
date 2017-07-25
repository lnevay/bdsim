#include "BDSOutputROOTEventSampler.hh"

#include "TTree.h"

#ifndef __ROOTBUILD__ 
#include "CLHEP/Units/SystemOfUnits.h"
#include "BDSSamplerHit.hh"
#endif

templateClassImp(BDSOutputROOTEventSampler)

template <class T>
BDSOutputROOTEventSampler<T>::BDSOutputROOTEventSampler():samplerName("sampler")
{
  this->Flush();
}

template <class T>
BDSOutputROOTEventSampler<T>::BDSOutputROOTEventSampler(std::string samplerNameIn):
  samplerName(samplerNameIn)
{
  this->Flush();
}

template
<class T> BDSOutputROOTEventSampler<T>::~BDSOutputROOTEventSampler()
{;}

#ifndef __ROOTBUILD__ 
template <class T>
void BDSOutputROOTEventSampler<T>::Fill(G4double E,
					G4double xIn,
					G4double yIn,
					G4double zIn,
					G4double xpIn,
					G4double ypIn,
					G4double zpIn,
					G4double tIn,
					G4double weightIn,
					G4int    partIDIn,
					G4int    /*nEvent*/,
					G4int    TurnsTaken,
					G4int    beamlineIndex)
{
  this->n++;
  this->z = (T) (zIn / CLHEP::m);
  this->S = (T) (0 / CLHEP::m);

  this->energy.push_back((T &&) (E / CLHEP::GeV));
  this->x.push_back((T &&) (xIn / CLHEP::m));
  this->y.push_back((T &&) (yIn / CLHEP::m));

  this->xp.push_back((T &&) (xpIn / CLHEP::radian));
  this->yp.push_back((T &&) (ypIn / CLHEP::radian));
  this->zp.push_back((T &&) (zpIn / CLHEP::radian));
  this->t.push_back((T &&) (tIn / CLHEP::ns));
  this->weight.push_back((const T &) weightIn);
  this->partID.push_back(partIDIn);
  this->parentID.push_back(0);
  this->turnNumber.push_back(TurnsTaken);
  this->modelID = beamlineIndex;
}

template <class T>
void BDSOutputROOTEventSampler<T>::Fill(BDSSamplerHit *hit)
{
  // get single values
  this->n++;
  this->z = (T) (hit->GetZ() / CLHEP::m);
  this->S = (T) (hit->GetS() / CLHEP::m);

  this->energy.push_back((T &&) (hit->GetTotalEnergy() / CLHEP::GeV));
  this->x.push_back((T &&) (hit->GetX() / CLHEP::m));
  this->y.push_back((T &&) (hit->GetY() / CLHEP::m));

  this->xp.push_back((T &&) (hit->GetXPrime() / CLHEP::radian));
  this->yp.push_back((T &&) (hit->GetYPrime() / CLHEP::radian));
  this->zp.push_back((T &&) (hit->GetZPrime() / CLHEP::radian));

  this->t.push_back((T &&) (hit->GetT() / CLHEP::ns));
  this->modelID = hit->GetBeamlineIndex();

#if 0
  this->X.push_back(hit->GetGlobalX()           / CLHEP::m);
  this->Y.push_back(hit->GetGlobalY()           / CLHEP::m);
  this->Z = hit->GetGlobalZ()                   / CLHEP::m;  
  this->Xp.push_back(hit->GetGlobalXPrime()     / CLHEP::radian);
  this->Yp.push_back(hit->GetGlobalYPrime()     / CLHEP::radian);
  this->Zp.push_back(hit->GetGlobalZPrime()     / CLHEP::radian);
#endif

  this->weight.push_back((T &&) hit->GetWeight());
  this->partID.push_back(hit->GetPDGtype());
  this->parentID.push_back(hit->GetParentID());
  this->trackID.push_back(hit->GetTrackID());
  this->turnNumber.push_back(hit->GetTurnsTaken());

}
//#else
//void BDSOutputROOTEventSampler::SetBranchAddress(TTree *)
//{}
#endif

template <class T> void BDSOutputROOTEventSampler<T>::SetBranchAddress(TTree *)
{;}

template <class T> void BDSOutputROOTEventSampler<T>::Flush()
{
  this->n = 0;
  this->energy.clear();
  this->x.clear();
  this->y.clear();
  this->z = 0.0;
  this->xp.clear();
  this->yp.clear();
  this->zp.clear();
  this->t.clear();
#if 0
  this->X.clear();
  this->Y.clear();
  this->Z = 0.0;      
  this->Xp.clear();
  this->Yp.clear();
  this->Zp.clear();
#endif
  this->weight.clear();
  this->partID.clear();
  this->parentID.clear();
  this->trackID.clear();
  this->turnNumber.clear();
  this->S = 0.0;
  this->modelID = -1;
}

template class BDSOutputROOTEventSampler<float>;
template class BDSOutputROOTEventSampler<double>;
