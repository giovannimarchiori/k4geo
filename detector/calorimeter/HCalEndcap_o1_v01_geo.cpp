
//====================================================================
//  Modified HCal Endcap Driver for the CLIC detector
//--------------------------------------------------------------------
//
//  Author     : N.Nikiforou
//  Adapted from PolyhedraBarrel Calorimeter by M. Frank
//
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"

using namespace std;
using namespace DD4hep;
using namespace DD4hep::Geometry;

static Ref_t create_detector(LCDD& lcdd, xml_h e, SensitiveDetector sens)  {
  xml_det_t   x_det     = e;
  int         det_id    = x_det.id();
  string      det_name  = x_det.nameStr();
  DetElement    sdet      (det_name,det_id);

  // --- create an envelope volume and position it into the world ---------------------
  
  Volume envelope = XML::createPlacedEnvelope( lcdd,  e , sdet ) ;
  
  if( lcdd.buildType() == BUILD_ENVELOPE ) return sdet ;
  
  //-----------------------------------------------------------------------------------
  

    
    
//   std::cout<<"Building HCal EndCap inside envelope."<<std::endl;
  xml_dim_t   dim       = x_det.dimensions();
  Material    air       = lcdd.air();
  int         nsides_inner = dim.nsides_inner();
  int         nsides_outer = dim.nsides_outer();
  double      rmin      = dim.rmin();
  double      rmax      = dim.rmax(); /// FIXME: IS THIS RIGHT?
  double      zmin      = dim.zmin();
  double      phi0      = dim.phi0();
  double       rcutout   = dim.hasAttr(_U(rmin2)) ? dim.rmin2() : 0;
  double       zcutout   = dim.hasAttr(_U(z2)) ? dim.z2() : 0;
  
  Layering    layering(x_det);
  double      totalThickness = layering.totalThickness();
  double hCal_cell_size      = lcdd.constant<double>("HCal_cell_size"); //Should try to obtain from segmentation
  
  PolyhedraRegular polyVolume(nsides_outer,rmin,rmax,totalThickness);
  PolyhedraRegular cutoutPolyVolume(nsides_inner,0,rmin+rcutout,zcutout);
  Position cutoutPos(0,0,(zcutout-totalThickness)/2.0);
  std::cout<<"Cutout z width will be  "<<zcutout<<std::endl; 
  
  Volume      endcapVol("endcap",SubtractionSolid(polyVolume,cutoutPolyVolume,cutoutPos),air);
  
  DetElement  endcapA(sdet,"endcap",det_id);
  Ref_t(endcapA)->SetName((det_name+"_A").c_str());
  
  int l_num = 1;
  int layerType   = 0;
  double layerZ   = -totalThickness/2;
  
  //Create caloData object to extend driver with data required for reconstruction
  DDRec::LayeredCalorimeterData* caloData = new DDRec::LayeredCalorimeterData ;
  caloData->layoutType = DDRec::LayeredCalorimeterData::EndcapLayout ;
  caloData->inner_symmetry = nsides_inner;
  caloData->outer_symmetry = nsides_outer; 
  caloData->phi0 = phi0;
  /// extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
  caloData->extent[0] = rmin ;
  caloData->extent[1] = rmax ; ///FIXME: CHECK WHAT IS NEEDED (EXSCRIBED?)
  caloData->extent[2] = zmin ;
  caloData->extent[3] = zmin + totalThickness;
  
  endcapVol.setAttributes(lcdd,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());
  
  for(xml_coll_t c(x_det,_U(layer)); c; ++c)  {
    xml_comp_t       x_layer  = c;
    double           l_thick  = layering.layer(l_num-1)->thickness();
    string           l_name   = _toString(layerType,"layer%d");
    int              l_repeat = x_layer.repeat();
    double            l_rcutout = x_layer.hasAttr(_U(gap)) ? x_layer.gap() : 0;
    
    std::cout<<"Number of layers in group "<<layerType<<" : "<<l_repeat<<std::endl; 
    
    Volume           l_vol(l_name,PolyhedraRegular(nsides_outer,rmin+l_rcutout,rmax,l_thick),air);
    
    vector<PlacedVolume> sensitives;
    
    int s_num = 1;
    double sliceZ = -l_thick/2;
    double totalAbsorberThickness=0.;
    for(xml_coll_t s(x_layer,_U(slice)); s; ++s)  {
      xml_comp_t x_slice = s;
      string     s_name  = _toString(s_num,"slice%d");
      double     s_thick = x_slice.thickness();
      Material   s_mat   = lcdd.material(x_slice.materialStr());
      Volume     s_vol(s_name,PolyhedraRegular(nsides_outer,rmin+l_rcutout,rmax,s_thick),s_mat);
      
      s_vol.setVisAttributes(lcdd.visAttributes(x_slice.visStr()));
      sliceZ += s_thick/2;
      PlacedVolume s_phv = l_vol.placeVolume(s_vol,Position(0,0,sliceZ));
      s_phv.addPhysVolID("slice",s_num);
      if ( x_slice.isSensitive() )  {
        sens.setType("calorimeter");
        s_vol.setSensitiveDetector(sens);
        sensitives.push_back(s_phv);
      }
      
      char val = x_slice.hasAttr(_U(radiator)) ? x_slice.attr < string > (_U(radiator))[0] : 'f';
      val = std::toupper(val);
      bool isAbsorber =  (val == 'T' || val == 'Y');
      
      if( isAbsorber ==true)
        totalAbsorberThickness+= s_thick;
      
      sliceZ += s_thick/2;
      s_num++;
    }
    l_vol.setVisAttributes(lcdd.visAttributes(x_layer.visStr()));
    if ( l_repeat <= 0 ) throw std::runtime_error(x_det.nameStr()+"> Invalid repeat value");
    for(int j=0; j<l_repeat; ++j) {
      string phys_lay = _toString(l_num,"layer%d");
      layerZ += l_thick/2;
      DetElement    layer_elt(endcapA, phys_lay, l_num);
      PlacedVolume  pv = endcapVol.placeVolume(l_vol,Position(0,0,layerZ));
      pv.addPhysVolID("layer", l_num);
      layer_elt.setPlacement(pv);
      for(size_t ic=0; ic<sensitives.size(); ++ic)  {
        PlacedVolume sens_pv = sensitives[ic];
        DetElement comp_elt(layer_elt,sens_pv.volume().name(),l_num);
        comp_elt.setPlacement(sens_pv);
      }
      
      ///FIXME: IS ORIENTATION RIGHT? WHICH SIDE DO WE NEED TO ADD TO STRUCTURE?
      DDRec::LayeredCalorimeterData::Layer caloLayer ;
      caloLayer.distance = zmin + layerZ;
      caloLayer.thickness = l_thick;
      caloLayer.absorberThickness = totalAbsorberThickness;
      caloLayer.cellSize0 = hCal_cell_size; //FIXME only temporary. Should get from Surfaces/Segmentation?
      caloLayer.cellSize1 = hCal_cell_size; //FIXME
      
      
      layerZ += l_thick/2;
      ++l_num;
    }
    ++layerType;
  }
  
  double z_pos = zmin+totalThickness/2;
  PlacedVolume pv;
  // Reflect it.

  DetElement  endcapB = endcapA.clone(det_name+"_B",x_det.id());
  
  //Removed rotations to align with envelope
  pv = envelope.placeVolume(endcapVol,Transform3D(RotationZYX(0,0,0),
                                                  Position(0,0,z_pos)));
  pv.addPhysVolID("barrel", 1);
  endcapA.setPlacement(pv);
  
  //Removed rotations
  pv = envelope.placeVolume(endcapVol,Transform3D(RotationZYX(0,M_PI,0),
                                                  Position(0,0,-z_pos)));
  pv.addPhysVolID("barrel", 2);
  endcapB.setPlacement(pv);
  
//   pv.addPhysVolID("system", det_id);
//   both_endcaps.setPlacement(pv);
//   sdet.add(sdetA);
  sdet.add(endcapB);
  
  sdet.addExtension< DDRec::LayeredCalorimeterData >( caloData ) ;
  
  return sdet;
  
}

DECLARE_DETELEMENT(HCalEndcap_o1_v01,create_detector)

