#ifndef TRKSAMPLERDATA_H
#define TRKSAMPLERDATA_H

#include <vector>

#include "TRKParticle.hh"


class TRKSamplerData
{
public:
  struct TRKSamplerDatum {
    TRKParticle particle;
    int turn;
    double s;
  };
  using SamplerParticles = std::vector<TRKSamplerDatum>;

  TRKSamplerData() = default;
  TRKSamplerData(int reserve);

  void RecordParticle(TRKParticle p, int turn, double s);

  using iterator = SamplerParticles::iterator;
  using const_iterator = SamplerParticles::const_iterator;
  iterator begin() { return particles.begin(); }
  iterator end() { return particles.end(); }
  const_iterator begin() const { return particles.begin(); }
  const_iterator end() const { return particles.end(); }
  const_iterator cbegin() const { return particles.cbegin(); }
  const_iterator cend() const { return particles.cend(); }

  std::pair<const_iterator, const_iterator> EventRange(int eventNumber);


private:
  SamplerParticles particles;
  bool sorted; // dirty bit
  void SortIfNotSorted();


};

#endif /* TRKSAMPLERDATA_H */
