#ifndef TRKTracker_h
#define TRKTracker_h

class TRKBunch;
class TRKLine;
class TRKStrategy;
class TRKElement;
class TRKParticle;
/**
 * @brief class that takes care of the tracker
 */
class TRKTracker {
private:
  TRKLine*  line;
  /// tracking type
  TRKStrategy* strategy;
  /// aperture check
  void CheckAperture(TRKBunch* bunch, TRKElement* element);
  /// maximum number of turns
  unsigned int maxTurns;
public :
  TRKTracker(TRKLine* line, TRKStrategy* strategy);
  ~TRKTracker();
  /// track bunch through beamline
  void Track(TRKBunch* bunch);
  /// check aperture method
  void Check();
  void Shift(TRKParticle& particle);
};

#endif
