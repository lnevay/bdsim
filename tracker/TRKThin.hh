#ifndef TRKThin_h
#define TRKThin_h

#include "TRKStrategy.hh"
#include "TRKElement.hh"

class TRKBunch;

/**
 * @brief Thin tracking
 */
class TRKThin : public TRKStrategy {
public:

  // todo, discuss, do we need tracking step for thin?
  TRKThin(int trackingSteps);
		     
  virtual ~TRKThin(); 

  /// returns tracking type
  //  TRKType trackingType()const {return type;}

  /**
   * Tracks a bunch through an element
   * @param[in]  el    The TRKElement.
   * @param[in]  bunch The bunch.
   */
  ///@{
  virtual void Track(TRKDrift*      el, TRKBunch* bunch);
  virtual void Track(TRKDipole*     el, TRKBunch* bunch);
  virtual void Track(TRKBend*       el, TRKBunch* bunch);
  virtual void Track(TRKQuadrupole* el, TRKBunch* bunch);
  virtual void Track(TRKSextupole*  el, TRKBunch* bunch);
  virtual void Track(TRKOctupole*   el, TRKBunch* bunch);
  virtual void Track(TRKSolenoid*   el, TRKBunch* bunch);
  ///@}

private:
  TRKThin(); ///< not implemented
};


#endif
