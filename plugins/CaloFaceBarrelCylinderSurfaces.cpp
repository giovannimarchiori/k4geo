//==========================================================================
// k4geo
//--------------------------------------------------------------------------
//
// Installs a cylindrical helper surface at the inner face of a cylindrical
// calorimeter barrel, to be used by the tracking to extrapolate tracks to
// the calorimeter face.
//
// Adapted from DD4hep_CaloFaceBarrelSurfacePlugin (DD4hep/DDDetectors/src/
// CaloFaceBarrel_surfaces.cpp, author F.Gaede, DESY), which approximates
// the face with a polygon of planar surfaces.
//
// Arguments:
//   <DetElement name>  name of the calorimeter the surface is attached to
//   length=            full length of the cylinder along z
//   radius=            radius of the cylinder
//   systemID=          system ID written in the surface cellID
//                      (side = 0, layer = 0, module = 0, sensor = 0)
//
//==========================================================================

#include <cmath>
#include <string>

namespace {
struct UserData {
  double length = 0.;
  double radius = 0.;
  int systemID = 0;
};
} // namespace

#define SURFACEINSTALLER_DATA UserData
#define DD4HEP_USE_SURFACEINSTALL_HELPER k4geo_CaloFaceBarrelCylinderSurfacePlugin
#include "DD4hep/Printout.h"
#include "DD4hep/SurfaceInstaller.h"
#include "DDRec/Surface.h"
#include "DDSegmentation/BitField64.h"

namespace {

/// cylindrical surface bounded in z by its length (the default bounds would use the shape of the volume)
class CaloBarrelCylinderImpl : public dd4hep::rec::VolCylinderImpl {
  double _length = 0.;

public:
  CaloBarrelCylinderImpl(dd4hep::Volume vol, dd4hep::rec::SurfaceType typ, double thickness_inner,
                         double thickness_outer, dd4hep::rec::Vector3D origin_val)
      : dd4hep::rec::VolCylinderImpl(vol, typ, thickness_inner, thickness_outer, origin_val) {}

  void setData(double length) { _length = length; }
  void setID(dd4hep::CellID id_val) { _id = id_val; }

  bool insideBounds(const dd4hep::rec::Vector3D& point, double epsilon) const override {
    return (std::abs(distance(point)) < epsilon) && (std::abs(point.z() - origin().z()) < _length / 2.);
  }

  double length_along_u() const override { return 2. * M_PI * origin().rho(); }
  double length_along_v() const override { return _length; }

  /// circles at both ends of the cylinder plus a few lines along z, for display
  std::vector<std::pair<dd4hep::rec::Vector3D, dd4hep::rec::Vector3D>> getLines(unsigned nMax) override {
    std::vector<std::pair<dd4hep::rec::Vector3D, dd4hep::rec::Vector3D>> lines;
    const unsigned n = (nMax > 8 ? nMax / 2 : 4);
    const double r = origin().rho();
    const double z0 = origin().z();
    for (unsigned i = 0; i < n; ++i) {
      const double phi0 = 2. * M_PI * i / n;
      const double phi1 = 2. * M_PI * (i + 1) / n;
      for (double z : {z0 - _length / 2., z0 + _length / 2.}) {
        lines.emplace_back(dd4hep::rec::Vector3D(r * cos(phi0), r * sin(phi0), z),
                           dd4hep::rec::Vector3D(r * cos(phi1), r * sin(phi1), z));
      }
      if (i % (n / 4) == 0) {
        lines.emplace_back(dd4hep::rec::Vector3D(r * cos(phi0), r * sin(phi0), z0 - _length / 2.),
                           dd4hep::rec::Vector3D(r * cos(phi0), r * sin(phi0), z0 + _length / 2.));
      }
    }
    return lines;
  }
};

/// handle for the surface, analogous to dd4hep::rec::VolCylinder
class CaloBarrelCylinder : public dd4hep::rec::VolSurface {
public:
  CaloBarrelCylinder(dd4hep::Volume vol, dd4hep::rec::SurfaceType typ, double thickness_inner, double thickness_outer,
                     dd4hep::rec::Vector3D origin_val)
      : dd4hep::rec::VolSurface(new CaloBarrelCylinderImpl(vol, typ, thickness_inner, thickness_outer, origin_val)) {}

  CaloBarrelCylinderImpl* operator->() { return static_cast<CaloBarrelCylinderImpl*>(_surf); }
};

template <>
void Installer<UserData>::handle_arguments(int argc, char** argv) {
  for (int i = 0; i < argc; ++i) {
    char* ptr = ::strchr(argv[i], '=');
    if (ptr) {
      std::string name(argv[i], ptr);
      double value = dd4hep::_toDouble(++ptr);

      printout(dd4hep::DEBUG, "k4geo_CaloFaceBarrelCylinderSurfacePlugin", "argument[%d] = %s = %f", i, name.c_str(),
               value);

      if (name == "length")
        data.length = value;
      else if (name == "radius")
        data.radius = value;
      else if (name == "systemID")
        data.systemID = value;
      else {
        printout(dd4hep::WARNING, "k4geo_CaloFaceBarrelCylinderSurfacePlugin", "unknown parameter:  %s ", name.c_str());
      }
    }
  }
}

template <typename UserData>
void Installer<UserData>::install(dd4hep::DetElement component, dd4hep::PlacedVolume pv) {
  dd4hep::Volume comp_vol = pv.volume();

  if (data.radius <= 0. || data.length <= 0.) {
    printout(dd4hep::ERROR, "k4geo_CaloFaceBarrelCylinderSurfacePlugin",
             "invalid radius (%f) or length (%f) for %s: no surface installed", data.radius, data.length,
             component.name());
    stopScanning();
    return;
  }

  printout(dd4hep::INFO, "k4geo_CaloFaceBarrelCylinderSurfacePlugin",
           "install cylindrical tracking surface for :  %s  (radius = %f mm, length = %f mm)", component.name(),
           data.radius / dd4hep::mm, data.length / dd4hep::mm);

  // same encoding as used by DDKalTest (and by DD4hep_CaloFaceBarrelSurfacePlugin)
  dd4hep::DDSegmentation::BitField64 bf("system:5,side:-2,layer:9,module:8,sensor:8");
  bf["system"] = data.systemID;

  const double inner_thickness = 1e-6;
  const double outer_thickness = 1e-6;

  dd4hep::rec::Vector3D o(data.radius, 0., 0.);

  CaloBarrelCylinder surf(comp_vol, Type(Type::Cylinder, Type::Helper, Type::Sensitive), inner_thickness,
                          outer_thickness, o);
  surf->setData(data.length);
  surf->setID(bf.getValue());

  addSurface(component, surf);

  stopScanning();
}

} // namespace
