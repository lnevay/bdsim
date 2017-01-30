#ifndef TRKElement_h
#define TRKElement_h

#include <string>
#include <iostream>

class TRKAperture;
class TRKBunch;
class TRKPlacement;
class TRKStrategy;
class TRKTiltOffset;

/**
 * @brief Base class element, also used for Drift
 */

class TRKElement
{ 
public:
  TRKElement(std::string   name,
	     double        length,
	     TRKAperture*  aperture,
	     TRKPlacement* placement,
	     double        strenght = 0);
  virtual ~TRKElement();

  /// track method, visitor pattern
  virtual void Track(TRKBunch* bunch, TRKStrategy* strategy);
  TRKBunch* CheckAperture(TRKBunch* bunch);
  // for when no aperture is used - we have to copy after state to before / confirm coordinates
  void ConfirmCoordinates(TRKBunch* bunch);

  inline std::string  GetName()     const {return name;}
  inline double       GetLength()   const {return length;}
  inline double       GetStrength() const {return strength;}
  inline TRKAperture* GetAperture() const {return aperture;}
  inline bool         Sample()      const {return sample;}
  
  /// output stream
  friend std::ostream& operator<< (std::ostream &out, const TRKElement &element);

  /// @{ Set and add offsets.
  void SetOffset(double x, double y);
  void SetTilt(double phi, double theta, double psi);
  void AddOffset(double x, double y);
  void AddTilt(double phi, double theta, double psi);
  ///@}

protected : 
  std::string    name;              ///< name of element -- do we need this? JS
  double         length;            ///< length of component [m]
  double         strength;          ///< strength parameter used differently in each derived class
  TRKAperture*   aperture;          ///< aperture of element
  TRKPlacement*  placement;         ///< location of element
  TRKTiltOffset* offsetIn;          ///< tilt and offset entrance of element
  TRKTiltOffset* offsetOut;         ///< tilt and offset exit of element
  bool           sample;            ///< whether to sample afterwards or not

  /// virtual print method for overloading operator<<. Virtual Friend Function Idiom
  virtual void Print(std::ostream& out) const;

private :
  TRKElement(); ///< not implemented
  
  /// global coordinates of local point - perhaps not needed
  double* LocalToGlobal(double /*vOut*/[]){return NULL;}
};

// declare drift
typedef TRKElement TRKDrift;

inline std::ostream& operator<< (std::ostream &out, const TRKElement &element)
{
  element.Print(out);
  return out;
}

#endif
