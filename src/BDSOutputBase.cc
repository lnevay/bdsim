#include "BDSOutputBase.hh"

#include "BDSGlobalConstants.hh"

BDSOutputBase::BDSOutputBase():eventNumber(0),outputFileNumber(-1)
{
  numberEventPerFile = BDSGlobalConstants::Instance()->NumberOfEventsPerNtuple();
}

void BDSOutputBase::Commit(const time_t&  startTime,
			   const time_t&  stopTime,
			   const G4float& duration,
                           const std::string& seedStateAtStart)
{
  Write(startTime, stopTime, duration, seedStateAtStart);
  Close();
  Initialise();
}
