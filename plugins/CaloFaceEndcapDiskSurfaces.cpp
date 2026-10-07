//==========================================================================
// k4geo
//--------------------------------------------------------------------------
//
// Installs two planar helper surfaces (at +z and -z) at the inner face of a
// calorimeter endcap with a circular outer edge, to be used by the tracking
// to extrapolate tracks to the calorimeter face.
//
// Adapted from DD4hep_CaloFaceEndcapSurfacePlugin (DD4hep/DDDetectors/src/
// CaloFaceEndcap_surfaces.cpp, author F.Gaede, DESY), which bounds the
// surface by a polygon.
//
// Arguments:
//   <DetElement name>  name of the calorimeter the surfaces are attached to
//   zpos=              |z| of the surfaces
//   radius=            outer radius of the disks
//   systemID=          system ID written in the surface cellID
//                      (side = +1/-1, layer = 0, module = 0, sensor = 0)
//
//==========================================================================

#include <cmath>
#include <string>

namespace {
struct UserData {
  double zpos = 0.;
  double radius = 0.;
  int systemID = 0;
};
} // namespace

#define SURFACEINSTALLER_DATA UserData
#define DD4HEP_USE_SURFACEINSTALL_HELPER k4geo_CaloFaceEndcapDiskSurfacePlugin
#include "DD4hep/Printout.h"
#include "DD4hep/SurfaceInstaller.h"
#include "DDRec/Surface.h"
#include "DDSegmentation/BitField64.h"

namespace {

/// planar surface orthogonal to z, bounded by a circle of given radius around the z axis
class CaloEndcapDiskImpl : public dd4hep::rec::VolPlaneImpl {
  double _r = 0.;

public:
  CaloEndcapDiskImpl(dd4hep::rec::SurfaceType typ, double thickness_inner, double thickness_outer,
                     dd4hep::rec::Vector3D u_val, dd4hep::rec::Vector3D v_val, dd4hep::rec::Vector3D n_val,
                     dd4hep::rec::Vector3D o_val, dd4hep::Volume vol, int id_val)
      : dd4hep::rec::VolPlaneImpl(typ, thickness_inner, thickness_outer, u_val, v_val, n_val, o_val, vol, id_val) {}

  void setData(double radius) { _r = radius; }
  void setID(dd4hep::CellID id_val) { _id = id_val; }

  bool insideBounds(const dd4hep::rec::Vector3D& point, double epsilon) const override {
    return (std::abs(distance(point)) < epsilon) && (point.rho() < _r);
  }

  /// outer circle, for display
  std::vector<std::pair<dd4hep::rec::Vector3D, dd4hep::rec::Vector3D>> getLines(unsigned nMax) override {
    std::vector<std::pair<dd4hep::rec::Vector3D, dd4hep::rec::Vector3D>> lines;
    const unsigned n = (nMax > 4 ? nMax : 4);
    for (unsigned i = 0; i < n; ++i) {
      const double phi0 = 2. * M_PI * i / n;
      const double phi1 = 2. * M_PI * (i + 1) / n;
      lines.emplace_back(dd4hep::rec::Vector3D(_r * cos(phi0), _r * sin(phi0), origin().z()),
                         dd4hep::rec::Vector3D(_r * cos(phi1), _r * sin(phi1), origin().z()));
    }
    return lines;
  }
};

typedef dd4hep::rec::VolSurfaceHandle<CaloEndcapDiskImpl> CaloEndcapDisk;

template <>
void Installer<UserData>::handle_arguments(int argc, char** argv) {
  for (int i = 0; i < argc; ++i) {
    char* ptr = ::strchr(argv[i], '=');
    if (ptr) {
      std::string name(argv[i], ptr);
      double value = dd4hep::_toDouble(++ptr);

      printout(dd4hep::DEBUG, "k4geo_CaloFaceEndcapDiskSurfacePlugin", "argument[%d] = %s = %f", i, name.c_str(),
               value);

      if (name == "zpos")
        data.zpos = value;
      else if (name == "radius")
        data.radius = value;
      else if (name == "systemID")
        data.systemID = value;
      else {
        printout(dd4hep::WARNING, "k4geo_CaloFaceEndcapDiskSurfacePlugin", "unknown parameter:  %s ", name.c_str());
      }
    }
  }
}

template <typename UserData>
void Installer<UserData>::install(dd4hep::DetElement component, dd4hep::PlacedVolume pv) {
  dd4hep::Volume comp_vol = pv.volume();

  if (data.radius <= 0. || data.zpos <= 0.) {
    printout(dd4hep::ERROR, "k4geo_CaloFaceEndcapDiskSurfacePlugin",
             "invalid radius (%f) or zpos (%f) for %s: no surface installed", data.radius, data.zpos, component.name());
    stopScanning();
    return;
  }

  printout(dd4hep::INFO, "k4geo_CaloFaceEndcapDiskSurfacePlugin",
           "install disk tracking surfaces for :  %s  (|z| = %f mm, radius = %f mm)", component.name(),
           data.zpos / dd4hep::mm, data.radius / dd4hep::mm);

  // same encoding as used by DDKalTest (and by DD4hep_CaloFaceEndcapSurfacePlugin)
  dd4hep::DDSegmentation::BitField64 bf("system:5,side:-2,layer:9,module:8,sensor:8");
  bf["system"] = data.systemID;

  const double inner_thickness = 1e-6;
  const double outer_thickness = 1e-6;

  // as in the DD4hep plugin, the origin is shifted off the z axis to pick up air instead of vacuum
  dd4hep::rec::Vector3D u(1., 0., 0.), v(0., 1., 0.), n(0., 0., 1.), o(0., 0.5 * data.radius, data.zpos);

  CaloEndcapDisk surf_pz(comp_vol, Type(Type::Helper, Type::Sensitive), inner_thickness, outer_thickness, u, v, n, o);
  bf["side"] = 1;
  surf_pz->setData(data.radius);
  surf_pz->setID(bf.getValue());
  addSurface(component, surf_pz);

  CaloEndcapDisk surf_nz(comp_vol, Type(Type::Helper, Type::Sensitive), inner_thickness, outer_thickness, u, v, n,
                         -1 * o);
  bf["side"] = -1;
  surf_nz->setData(data.radius);
  surf_nz->setID(bf.getValue());
  addSurface(component, surf_nz);

  stopScanning();
}

} // namespace
