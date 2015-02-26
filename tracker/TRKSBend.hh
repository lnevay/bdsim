#ifndef TRKSBend_h
#define TRKSBend_h

#include "TRKElement.hh"

/**
 * @brief bend
 * 
 * dipole tracking but reference s-position is changed
 */
class TRKSBend : public TRKElement {
private: 
  /// bending angle in rad
  double angle;

public:
  TRKSBend(double angle, std::string name, double length, TRKAperture *aperture, TRKPlacement *placement);
  ~TRKSBend();
  
  double GetAngle()const{return angle;}
  
  virtual void Track(TRKBunch* bunch, TRKStrategy* strategy);
  
protected:
  /// output stream
  virtual void Print(std::ostream& out) const;

private:
  TRKSBend(); ///< not implemented
};

#endif
