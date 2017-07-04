#ifndef BDSAUXILIARYNAVIGATOR_H
#define BDSAUXILIARYNAVIGATOR_H

#include "globals.hh" // geant4 types / globals
#include "G4AffineTransform.hh"
#include "G4Navigator.hh"
#include "G4ThreeVector.hh"

class BDSStep;
class G4Step;
class G4VPhysicalVolume;

/**
 * @brief Extra G4Navigator to get coordinate transforms.
 * 
 * All BDSIM steppers and magnetic fields require the ability 
 * to convert from global to local coordinates. The prescribed method 
 * to do this is by using a G4Navigator instance. There is of course the 
 * main navigator for tracking, but requesting information on a global
 * point changes the state of the navigator - ie the particle being 
 * tracked by the navigator is now that point. This is VERY dangerous.
 * 
 * This class provides a static auxiliary navigator that each derived 
 * class can use. Making the auxilliary navigator static is not done 
 * to reduce memory usage but because navigating from an unknown place 
 * to anywhere in the geometry is
 * much more costly than a relative move in the geometry. If we only
 * use one auxiliary navigator, it will always be relatively close
 * in the geometry even if a different stepper has been used and 
 * therefore more efficient. This is important as Geant4 may use the
 * steppers at least two or three times to estimate the best next step
 * and the stepper itself may make three steps (full, and two half) to
 * estimate the error in the integration.
 *
 * This class in fact provides access two static navigators. The intended
 * use is one for the real world and one for the read out geometry / world
 * for curvilinear coordinates.  All functions have an optional last argument
 * to select which navigator is required - the default is the curvilinear one.
 * 
 * @author Laurie Nevay
 */

class BDSAuxiliaryNavigator
{
public:
  BDSAuxiliaryNavigator();
  ~BDSAuxiliaryNavigator();

  /// Setup the navigator w.r.t. to a world volume - typically real world.
  static void AttachWorldVolumeToNavigator(G4VPhysicalVolume* worldPVIn)
  {auxNavigator->SetWorldVolume(worldPVIn); worldPV = worldPVIn;}

  /// Setup the navigator w.r.t. to the read out world / geometry to provide
  /// curvilinear coordinates.
  static void AttachWorldVolumeToNavigatorCL(G4VPhysicalVolume* curvilinearWorldPVIn)
  {auxNavigatorCL->SetWorldVolume(curvilinearWorldPVIn); curvilinearWorldPV = curvilinearWorldPVIn;}

  static void RegisterCurvilinearBridgeWorld(G4VPhysicalVolume* curvilinearBridgeWorldPVIn)
  {auxNavigatorCLB->SetWorldVolume(curvilinearBridgeWorldPVIn); curvilinearBridgeWorldPV = curvilinearBridgeWorldPVIn;}

  /// A wrapper for the underlying static navigator instance located within this class.
  G4VPhysicalVolume* LocateGlobalPointAndSetup(const G4ThreeVector& point,
					       const G4ThreeVector* direction = nullptr,
					       const G4bool pRelativeSearch   = true,
					       const G4bool ignoreDirection   = true,
					       G4bool useCurvilinear          = true) const;

  /// A safe way to locate and setup a point inside a volume. This uses the mid
  /// point of the step to locate the volume rather than the edges which may lie
  /// on a boundary and typically find the world or previous volume.
  G4VPhysicalVolume* LocateGlobalPointAndSetup(G4Step const* const step,
					       G4bool useCurvilinear = true) const;

  /// Calculate the local coordinates for both a pre and post step point. The mid point
  /// of the step is used for the volume (and therefore transform) lookup which should
  /// ensure the correct volume is found - avoiding potential boundary issues between
  /// parallel worlds. This uses LocateGlobalPointAndSetup.
  BDSStep ConvertToLocal(G4Step const* const step,
			 G4bool useCurvilinear = true) const;

  /// Calculate the local coordinates for a position and direction along a step
  /// length.  This is similar to the same function but for a G4Step but split
  /// apart. The direction vector can be used as the momentum vector without being
  /// a unit vector.  The 'post step' vector in the BDSStep instance will be the
  /// direction vector (of same magnitude) but rotated to the local frame.
  /// This uses LocateGlobalPointAndSetup. Can control default lookahead distance
  /// default (1 -> 1mm).
  BDSStep ConvertToLocal(const G4ThreeVector& globalPosition,
			 const G4ThreeVector& globalDirection,
			 const G4double       stepLength     = 0,
			 const G4bool&        useCurvilinear = true,
			 const G4double       marginLength   = 1) const;

  /// Convert back to global coordinates.  This DOES NOT update the transforms
  /// and uses the existing transforms inside this class - ie this relies on
  /// ConvertToLocal being used beforehand to initialise the transforms.  This is
  /// done as we'd need a look up point behind this point.
  BDSStep ConvertToGlobalStep(const G4ThreeVector& localPosition,
			      const G4ThreeVector& localDirection,
			      const G4bool&        useCurvilinear = true) const;

  /// Calculate the local coordinates of a global point.
  G4ThreeVector ConvertToLocal(const G4double globalPoint[3],
			       const G4bool&  useCurvilinear = true) const;

  /// Vector version
  G4ThreeVector ConvertToLocal(const G4ThreeVector& globalPosition,
			       const G4bool&        useCurvilinear = true) const;

  /// Calculate the local version of a global vector (axis). This is done
  /// w.r.t. a point - ie, the point is used to initialise the transform
  /// if not done already.
  G4ThreeVector ConvertAxisToLocal(const G4double globalPoint[3],
				   const G4double globalAxis[3],
				   const G4bool&  useCurvilinear = true) const;

  /// Vector version
  G4ThreeVector ConvertAxisToLocal(const G4ThreeVector& globalPoint,
				   const G4ThreeVector& globalAxis,
				   const G4bool&        useCurvilinear = true) const;

  /// Convert a vector (axis) from local to global coordinates. NOTE this
  /// function must only be used once the instance of this class has been initialised,
  /// setting up the transforms.  It is up to the developer to ensure this,
  /// otherwise you'll find a bad access.
  G4ThreeVector ConvertAxisToGlobal(const G4ThreeVector& localAxis,
				    const G4bool&        useCurvilinear = true) const;

  /// Convert a vector (axis) from local to global coordinates. Note this function
  /// must only be used once the instance of this class has been initialised, setting
  /// up the transforms. It is up to the developer to ensure this. This utility function
  /// operates on two threevectors in a pair.
  std::pair<G4ThreeVector, G4ThreeVector> ConvertAxisToGlobal(const std::pair<G4ThreeVector, G4ThreeVector>& localAxis,
							      const G4bool&        useCurvilinear = true) const;

  /// Convert a position in local coordinates to global coordinates. NOTE a
  /// similar caution to ConvertAxisToGlobal applies.
  G4ThreeVector ConvertToGlobal(const G4ThreeVector& localPosition,
				const G4bool&        useCurvilinear = true) const;

  /// Convert a vector (axis) from local to global coordinates. This uses a global
  /// position to ensure the transform is initialised.
  G4ThreeVector ConvertAxisToGlobal(const G4ThreeVector& globalPosition,
				    const G4ThreeVector& localAxis,
				    const G4bool&        useCurvilinear = true) const;

  /// Convert a position in local coordinates to global coordinates. This uses a
  /// global position to ensure the transform is initialised.
  G4ThreeVector ConvertToGlobal(const G4ThreeVector& globalPosition,
				const G4ThreeVector& localPosition,
				const G4bool&        useCurvilinear = true) const;
  
protected:
  mutable G4AffineTransform globalToLocal;
  mutable G4AffineTransform localToGlobal;
  mutable G4AffineTransform globalToLocalCL;
  mutable G4AffineTransform localToGlobalCL;
  
  /// Navgiator object for safe navigation in the real (mass) world without
  /// affecting tracking of the particle.
  static G4Navigator* auxNavigator;

  /// Navigator object for curvilinear world that contains simple cylinders
  /// for each element whose local coordinates represent the curvilinear coordinate
  /// system.
  static G4Navigator* auxNavigatorCL;

  /// Navigator object for bridge world. This contains briding volumes for the
  /// gaps in the curivlinear world. It therefore acts as a fall back if we find
  /// the world volume when we know we really shouldn't.
  static G4Navigator* auxNavigatorCLB;

private:
  /// Utility function to select appropriate navigator
  G4Navigator* Navigator(G4bool curvilinear) const;

  /// @{ Utility function to select appropriate transfrom.
  const G4AffineTransform& GlobalToLocal(G4bool curvilinear) const;
  const G4AffineTransform& LocalToGlobal(G4bool curvilinear) const;
  /// @}

  void InitialiseTransform(const G4bool& massworld        = true,
			   const G4bool& curvilinearWorld = true) const;
  
  /// Locate the supplied point the in the geometry and get and store
  /// the transform to that volume in the member variable. This function
  /// has to be const as it's called the first time in GetField which is
  /// a pure virtual const function from G4MagneticField that we have to
  /// implement and have to keep const. This function doesn't change the
  /// const pointer but does change the contents of what it points to.
  void InitialiseTransform(const G4ThreeVector& globalPosition) const;

  /// This is used to foricibly initialise the transforms using a position,
  /// momentum vector and step length. The free drift of the particle is
  /// calculated and the average of the two points is used to locate
  /// and initialise the transforms (in global coordinates).
  void InitialiseTransform(const G4ThreeVector& globalPosition,
			   const G4ThreeVector& globalMomentum,
			   const G4double       stepLength);
  
  /// Counter to keep track of when the last instance of the class is deleted
  /// and therefore when the navigators can be safely deleted without affecting
  static G4int numberOfInstances;
  
  /// @{ Cache of world PV to test if we're getting the wrong volume for the transform.
  static G4VPhysicalVolume* worldPV;
  static G4VPhysicalVolume* curvilinearWorldPV;
  static G4VPhysicalVolume* curvilinearBridgeWorldPV;
  /// @}
  
  /// Margin by which to advance the point along the step direction if the
  /// world volume is found for transforms. This is in an attempt to find a
  /// real curvilinear volume. Therefore, this should be greater than
  /// lengthSafety or the geometrical tolerance.
  G4double volumeMargin;
};


#endif
