#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TMatrixT.h"
#include "XML/Utilities.h"
#include <DDRec/DetectorData.h>

namespace det {

namespace ECalEndcap_ParallelPlateTurbine_o1_v01 {
  unsigned ECalEndCapElementCounter = 0;

  double ngroup(double r, double dz, double alpha, double d) {
    
    // for an endcap with inner radius r and thickness dz compute the
    // maximum number of groups of thickness d that can be placed if the
    // tilt angle of the boards around the radial axis is alpha
    
    double phi1 = acos(dz/tan(alpha)/2/r);
    double phi2 = acos((dz/tan(alpha)/2+ d/sin(alpha))/r);

    double delta = std::abs(phi1-phi2);
    if ( delta > M_PI ) {
      delta = 2*M_PI-delta;
    }
    
    // number of groups times delta needs to be smaller than 2pi
    double ngroup_max = 2*M_PI/delta;

    return ngroup_max;
}

  void buildOneSide_Turbine(dd4hep::Detector& aLcdd, dd4hep::DetElement& caloDetElem,
                            dd4hep::SensitiveDetector& aSensDet, dd4hep::Volume& aEnvelope,
                            dd4hep::xml::Handle_t& aXmlElement, unsigned& iGroup) {

    dd4hep::xml::DetElement calo = aXmlElement.child(_Unicode(calorimeter));
    dd4hep::xml::Dimension caloDim(calo.dimensions());

    dd4hep::xml::DetElement group = calo.child(_Unicode(turbineGroup));
    std::string nobleLiquidMaterial = group.materialStr();

    dd4hep::xml::DetElement stack = group.child(_Unicode(stack));

    dd4hep::xml::DetElement xmlDetElem = aXmlElement;
    std::string nameDet = xmlDetElem.nameStr();

    dd4hep::xml::Dimension dim(aXmlElement.child(_Unicode(dimensions)));

    // build cryostat
    //  Retrieve cryostat data
    dd4hep::xml::DetElement cryostat = calo.child(_Unicode(cryostat));
    dd4hep::xml::Dimension cryoDim(cryostat.dimensions());
    double cryoThicknessFront = aLcdd.constant<float>("CryoEMECThicknessFront");
    double cryoThicknessBack = aLcdd.constant<float>("CryoEMECThicknessBack");
    double bathThicknessFront = aLcdd.constant<float>("BathThicknessFront");
    double bathThicknessBack = aLcdd.constant<float>("BathThicknessBack");

    dd4hep::xml::DetElement cryoFront = cryostat.child(_Unicode(front));
    dd4hep::xml::DetElement cryoBack = cryostat.child(_Unicode(back));
    dd4hep::xml::DetElement cryoInner = cryostat.child(_Unicode(inner));
    dd4hep::xml::DetElement cryoOuter = cryostat.child(_Unicode(outer));

    bool cryoFrontSensitive = cryoFront.isSensitive();
    bool cryoBackSensitive = cryoBack.isSensitive();
    bool cryoInnerSensitive = cryoInner.isSensitive();
    bool cryoOuterSensitive = cryoOuter.isSensitive();

    double bathRmin = cryoDim.rmin2(); // - margin for inclination
    double bathRmax = cryoDim.rmax1(); // + margin for inclination
    double bathDelZ = cryoDim.dz();
    dd4hep::Tube bathOuterShape(bathRmin, bathRmax, bathDelZ); // make it 4 volumes + 5th for detector envelope

    dd4hep::printout(dd4hep::INFO, "ECalEndcap_ParallelPlateTurbine_o1_v01", "Cryostat front thickness is %f", cryoDim.rmin2());

    float NLcenterZ = (cryoThicknessFront - cryoThicknessBack) / 2.;
    if (cryoThicknessFront > 0) {
      // 1. Create cryostat
      dd4hep::Tube cryoFrontShape(cryoDim.rmin1(), cryoDim.rmax2(), cryoThicknessFront / 2.);
      dd4hep::Tube cryoBackShape(cryoDim.rmin1(), cryoDim.rmax2(), cryoThicknessBack / 2.);
      dd4hep::Tube cryoInnerShape(cryoDim.rmin1(), cryoDim.rmin2(), cryoDim.dz());
      dd4hep::Tube cryoOuterShape(cryoDim.rmax1(), cryoDim.rmax2(), cryoDim.dz());
      dd4hep::printout(dd4hep::INFO, "ECalEndcap_ParallelPlateTurbine_o1_v01",
                       "ECAL endcap cryostat: front: rmin (cm) = %f rmax (cm) = %f dz (cm) = %f ", cryoDim.rmin1(),
                       cryoDim.rmin2(), cryoDim.dz());
      dd4hep::printout(dd4hep::INFO, "ECalEndcap_ParallelPlateTurbine_o1_v01",
                       "ECAL encdap cryostat: back: rmin (cm) =  %f rmax (cm) = %f dz (cm) = %f", cryoDim.rmax1(),
                       cryoDim.rmax2(), cryoDim.dz());
      dd4hep::printout(dd4hep::INFO, "ECalEndcap_ParallelPlateTurbine_o1_v01",
                       "ECAL endcap cryostat: side: rmin (cm) =  %f rmax (cm) = %f dz (cm) = %f", cryoDim.rmin2(),
                       cryoDim.rmax1(), cryoDim.dz() - caloDim.dz());
      dd4hep::printout(dd4hep::INFO, "ECalEndcap_ParallelPlateTurbine_o1_v01", "Cryostat is made out of %s",
                       cryostat.materialStr().c_str());

      dd4hep::Volume cryoFrontVol(cryostat.nameStr() + "_front", cryoFrontShape,
                                  aLcdd.material(cryostat.materialStr()));
      dd4hep::Volume cryoBackVol(cryostat.nameStr() + "_back", cryoBackShape, aLcdd.material(cryostat.materialStr()));
      dd4hep::Volume cryoInnerVol(cryostat.nameStr() + "_inner", cryoInnerShape,
                                  aLcdd.material(cryostat.materialStr()));
      dd4hep::Volume cryoOuterVol(cryostat.nameStr() + "_outer", cryoOuterShape,
                                  aLcdd.material(cryostat.materialStr()));

      dd4hep::Position cryoFrontPos(0, 0, NLcenterZ - cryoDim.dz() - cryoThicknessFront / 2.);
      dd4hep::PlacedVolume cryoFrontPhysVol = aEnvelope.placeVolume(cryoFrontVol, cryoFrontPos);
      dd4hep::Position cryoBackPos(0, 0, NLcenterZ + cryoDim.dz() + cryoThicknessBack / 2.);
      dd4hep::PlacedVolume cryoBackPhysVol = aEnvelope.placeVolume(cryoBackVol, cryoBackPos);
      dd4hep::Position cryoInnerOuterPos(0, 0, NLcenterZ);
      dd4hep::PlacedVolume cryoInnerPhysVol = aEnvelope.placeVolume(cryoInnerVol, cryoInnerOuterPos);
      dd4hep::PlacedVolume cryoOuterPhysVol = aEnvelope.placeVolume(cryoOuterVol, cryoInnerOuterPos);
      unsigned sidetype = 0x4; // probably not needed anymore...
      if (cryoFrontSensitive) {
        cryoFrontVol.setSensitiveDetector(aSensDet);
        cryoFrontPhysVol.addPhysVolID("cryo", 1);
        cryoFrontPhysVol.addPhysVolID("type", sidetype + 1);
        dd4hep::printout(dd4hep::INFO, "ECalEndcap_ParallelPlateTurbine_o1_v01", "Cryostat front volume set as sensitive");
      }
      if (cryoBackSensitive) {
        cryoBackVol.setSensitiveDetector(aSensDet);
        cryoBackPhysVol.addPhysVolID("cryo", 1);
        cryoBackPhysVol.addPhysVolID("type", sidetype + 2);
        dd4hep::printout(dd4hep::INFO, "ECalEndcap_ParallelPlateTurbine_o1_v01", "Cryostat back volume set as sensitive");
      }
      if (cryoInnerSensitive) {
        cryoInnerVol.setSensitiveDetector(aSensDet);
        cryoInnerPhysVol.addPhysVolID("cryo", 1);
        cryoInnerPhysVol.addPhysVolID("type", sidetype + 3);
        dd4hep::printout(dd4hep::INFO, "ECalEndcap_ParallelPlateTurbine_o1_v01", "Cryostat inner volume set as sensitive");
      }
      if (cryoOuterSensitive) {
        cryoOuterVol.setSensitiveDetector(aSensDet);
        cryoOuterPhysVol.addPhysVolID("cryo", 1);
        cryoOuterPhysVol.addPhysVolID("type", sidetype + 4);
        dd4hep::printout(dd4hep::INFO, "ECalEndcap_ParallelPlateTurbine_o1_v01", "Cryostat outer volume set as sensitive");
      }
      dd4hep::DetElement cryoFrontDetElem(caloDetElem, "cryo_front", 0);
      cryoFrontDetElem.setPlacement(cryoFrontPhysVol);
      dd4hep::DetElement cryoBackDetElem(caloDetElem, "cryo_back", 0);
      cryoBackDetElem.setPlacement(cryoBackPhysVol);
      dd4hep::DetElement cryoInnerDetElem(caloDetElem, "cryo_inner", 0);
      cryoInnerDetElem.setPlacement(cryoInnerPhysVol);
      dd4hep::DetElement cryoOuterDetElem(caloDetElem, "cryo_outer", 0);
      cryoOuterDetElem.setPlacement(cryoOuterPhysVol);
    }

    // 2. Create noble liquid bath
    dd4hep::Volume bathVol(nobleLiquidMaterial + "_bath", bathOuterShape, aLcdd.material(nobleLiquidMaterial));
    dd4hep::printout(dd4hep::INFO, "ECalEndcap_ParallelPlateTurbine_o1_v01",
                     "ECAL endcap bath: material = %s rmin (cm) = %f rmax (cm) = %f, dz (cm) = %f, thickness in front "
                     "of ECal (cm) = %f,  thickness behind ECal (cm) = %f",
                     nobleLiquidMaterial.c_str(), bathRmin, bathRmax, caloDim.dz(), caloDim.rmin() - cryoDim.rmin2(),
                     cryoDim.rmax1() - caloDim.rmax());

    //    dd4hep::Position bathPos(0, 0, (cryoThicknessFront - cryoThicknessBack)/2.);
    dd4hep::Position bathPos(0, 0, NLcenterZ);

    dd4hep::PlacedVolume bathPhysVol = aEnvelope.placeVolume(bathVol, bathPos);

    dd4hep::DetElement bathDetElem(caloDetElem, "bath", 1);

    bathDetElem.setPlacement(bathPhysVol);

    // 3. Create detector structure
    double length = dim.dz() * 2.;
    double zOffsetEnvelope = -length / 2.;

    dd4hep::xml::DetElement supportTubeElem = calo.child(_Unicode(supportTube));
    float supportTubeThickness = supportTubeElem.thickness();

    double rmin = bathRmin;
    double rmax = bathRmax;
    double ri = rmin + supportTubeThickness; // here the actual stacks start
    

    dd4hep::Tube supportTube(rmin, ri, bathDelZ);

    dd4hep::Volume supportTubeVol("supportTube", supportTube, aLcdd.material(supportTubeElem.materialStr()));
    if (supportTubeElem.isSensitive()) {
      supportTubeVol.setSensitiveDetector(aSensDet);
    }
    dd4hep::PlacedVolume supportTube_pv =
      bathVol.placeVolume(supportTubeVol, dd4hep::Position(0, 0, zOffsetEnvelope + dim.dz()));
    supportTube_pv.addPhysVolID("cryo", 1);
    dd4hep::DetElement supportTubeDetElem(bathDetElem, "supportTube", 0);
    supportTubeDetElem.setPlacement(supportTube_pv);


    /*
     * Build one group of nStackMin to nStackMax stacks
     *
     * Number of virtual wheels is nStackMax - nStackMin + 1.  The first,
     * innermost, virtual wheel has nStackMin stacks and each following
     * virtual wheel further out in radius has one stack more. The radii are
     * computed by scanning r in steps of 10 mm until the number of placeable
     * groups (including clearance) is less or equal the chosen number of
     * groups.
     * 
     */

    unsigned nStackMin = group.attr<int>(_Unicode(nStackMin));
    unsigned nStackMax = group.attr<int>(_Unicode(nStackMax));
    unsigned nGroups   = group.attr<int>(_Unicode(nGroups));

    /*
     * each stack is made of slices of different material 
     */
    
    std::vector<double> dSlice;        // full thickness of each slice
    //std::vector<bool> visSlice;        // visibility of each slice
    std::vector<bool> sensSlice;       // sensitivity of each slice
    std::vector<std::string> matSlice; // material of each slice
    
    double dStack = 0; // thickness of one stack

    for (dd4hep::xml::Collection_t l(stack, _Unicode(slice)); l; ++l) {
      double d = xml_comp_t(l).thickness();
      dSlice.push_back(d);
      dStack += d;
      sensSlice.push_back(xml_comp_t(l).isSensitive());
      matSlice.push_back(xml_comp_t(l).materialStr());
    }

    double dClearance = aLcdd.constant<float>("MinGroupClearance");

    unsigned nVirtualWheels = nStackMax-nStackMin+1;

    double alpha = group.attr<float>(_Unicode(angle));

    
    std::vector<double> rVirtWheel(nVirtualWheels);
    std::vector<dd4hep::Tube> virtualCyl(nVirtualWheels);

    // construct all radii
    double r=ri;
    double dz = bathDelZ * 2 - bathThicknessFront - bathThicknessBack;
    for(unsigned i=0;i<nVirtualWheels;i++) {
      bool foundR = false;
      double d = (nStackMin + i)*dStack + dClearance;
      while (!foundR) {
	double n = ngroup(r,dz,alpha,d);
	if ( n > nGroups ) {
	  foundR = true;
	  rVirtWheel[i] = r;
	  dd4hep::printout(dd4hep::DEBUG, "ECalEndcap_ParallelPlateTurbine_o1_v01", "Making virtual wheel with inner radius %f and %d stacks.", r, nStackMin+i);
	  virtualCyl[i] = dd4hep::Tube(r,rmax,dz/2,(90-45)*dd4hep::deg,(90+45)*dd4hep::deg); // one virtual cylinder quadrant from 45 degrees to 135 degrees 
	}
	else {
	  r += 10*dd4hep::mm;
	  if ( r > rmax ) {
	    dd4hep::printout(dd4hep::ERROR, "ECalEndcap_ParallelPlateTurbine_o1_v01",
			     "The requested number of groups (%d) can not be placed for %d stacks.",nGroups,nStackMin+i);
	    return;
	  }
	}
      }
    }

    /*
     * create basic slice boxes
     */
    
    std::vector<dd4hep::Solid> sliceBoxes;  // basic shapes prior to transform and intersection (nSlice)
    std::vector<dd4hep::Solid> sliceShapes; // shapes after transform and intersection (nStack x nSlice)
    std::vector<dd4hep::Volume> sliceVolumes; // volumes corresponding to shapes after transform and intersection (nStack x nSlice)
    double ds=0; // local deltaStack - grows with layers ...
    for(unsigned j=0;j<dSlice.size();j++) { // loop over number of slices per stack
      sliceBoxes.push_back(dd4hep::Box(dSlice[j]/2.,
				       1.9*(rmax-ri)/2,
				       1.5*(dz/2./sin(alpha)+(ds+dSlice[j])*cos(alpha))));
      ds += dSlice[j];
    }

    dd4hep::Solid groupShape; // shape of the entire group

    /*
     * create slice shapes
     */
    unsigned k=0;
    for(int i=0;i<nStackMax;i++) { // loop over number of stacks per group
      double dx=0;
      unsigned iCyl = (i<nStackMin?0:i-nStackMin+1);
      for(unsigned j=0;j<dSlice.size();j++) { // loop over number of slices per stack
	sliceShapes.push_back
	  (dd4hep::IntersectionSolid   // intersect virtual cylinder quadrant with transformed slice box
	   (virtualCyl[iCyl],          // virtual cylinder quadrant
	    sliceBoxes[j],             // basic slice shape 
	    dd4hep::Transform3D
	    (dd4hep::RotationZYX     
	     (0,-(M_PI/2. - alpha),0), // rotate slice around y-axis
	     dd4hep::Position
	     ((-i)*dStack-0.5*dSlice[j]-dx,  // offset normal to slice plane
	     (rmax+rVirtWheel[iCyl])/2,      // offset along slice plane in vertical direction
	      0))));
	
	sliceVolumes.push_back
	  (dd4hep::Volume
	   (std::string("slice_")+std::to_string(i)+"_"+std::to_string(j),
	    sliceShapes[k],aLcdd.material(matSlice[j])));
	if ( sensSlice[j] ) {
	  sliceVolumes[k].setSensitiveDetector(aSensDet);
	}
	dx += dSlice[j];
	if ( k == 0 ) {
	  // first solid is copied to groupShape
	  groupShape = sliceShapes[0];
	}
	else {
	  // build composite of group up to this shape + this shape
	  dd4hep::UnionSolid theUnion(groupShape,sliceShapes[k]);
	  groupShape = theUnion;
	}
	k++;
      }
    }

    dd4hep::Volume groupVolume("groupVol",groupShape,aLcdd.material("Air")); // volume of one entire group

    dd4hep::Position posZero(0, 0, 0);
    for(unsigned i=0;i<sliceVolumes.size();i++) { // loop over all slice volumes and add to group volume
      dd4hep::PlacedVolume sliceVol_pv = groupVolume.placeVolume(sliceVolumes[i], posZero);
      sliceVol_pv.addPhysVolID("stack",(int)(i/dSlice.size())); // stack index
      sliceVol_pv.addPhysVolID("slice",(i%dSlice.size())); // slice index
    }
	
    /*
     * now build all nGroup groups by rotating them in phi
     */
    
    for(unsigned i=0;i<nGroups;i++) {
      dd4hep::RotationZ rotZ(2*M_PI/nGroups*i-M_PI/2); // starting at phi=0 for i=0 means to rotate by -90 degrees
      dd4hep::PlacedVolume placedGroup = aEnvelope.placeVolume(groupVolume,dd4hep::Transform3D(rotZ, posZero));
      placedGroup.addPhysVolID("group", i);
    }

    iGroup = nGroups;
    dd4hep::printout(dd4hep::DEBUG, "ECalEndcap_ParallelPlateTurbine_o1_v01", "Total number of groups:  %d", nGroups);
    return;
  }

  static dd4hep::Ref_t createECalEndcapTurbine(dd4hep::Detector& aLcdd, dd4hep::xml::Handle_t aXmlElement,
                                               dd4hep::SensitiveDetector aSensDet) {

    dd4hep::xml::DetElement xmlDetElem = aXmlElement;
    std::string nameDet = xmlDetElem.nameStr();
    int idDet = xmlDetElem.id();
    dd4hep::xml::Dimension dim(xmlDetElem.dimensions());
    dd4hep::DetElement caloDetElem(nameDet, idDet);
    dd4hep::xml::Dimension sdType = xmlDetElem.child(_U(sensitive));
    aSensDet.setType(sdType.typeStr());

    unsigned numReadoutRhoLayers, numReadoutZLayers;
    // Create air envelope for one endcap (will be copied to make both endcaps)
    dd4hep::Tube endcapShape(dim.rmin1(), dim.rmax1(), dim.dz());

    dd4hep::Volume envelopeVol(nameDet + "_vol", endcapShape, aLcdd.material("Air"));

    dd4hep::printout(dd4hep::DEBUG, "ECalEndcap_ParallelPlateTurbine_o1_v01",
                     "Placing detector on the positive side: (cm) %f  with min, max radii %f %f and depth %f",
                     dim.z_offset(), dim.rmin1(), dim.rmax1(), dim.dz());

    unsigned iGroup = 0;
    buildOneSide_Turbine(aLcdd, caloDetElem, aSensDet, envelopeVol, aXmlElement, iGroup);

    dd4hep::Assembly endcapsAssembly("ECalEndcaps_turbine");

    // Place the envelope
    dd4hep::Transform3D envelopePositiveVolume_tr(dd4hep::RotationZYX(0, 0, 0),
                                                  dd4hep::Translation3D(0, 0, dim.z_offset()));
    dd4hep::PlacedVolume envelopePositivePhysVol = endcapsAssembly.placeVolume(envelopeVol, envelopePositiveVolume_tr);
    envelopePositivePhysVol.addPhysVolID("side", 1);

    // make another placement for the negative z endcap
    dd4hep::Transform3D envelopeNegativeVolume_tr(dd4hep::RotationZYX(0, 0, 180 * dd4hep::deg),
                                                  dd4hep::Translation3D(0, 0, -dim.z_offset()));
    dd4hep::PlacedVolume envelopeNegativePhysVol = endcapsAssembly.placeVolume(envelopeVol, envelopeNegativeVolume_tr);
    envelopeNegativePhysVol.addPhysVolID("side", -1);

    dd4hep::Volume motherVol = aLcdd.pickMotherVolume(caloDetElem);
    dd4hep::PlacedVolume envelopePhysVol = motherVol.placeVolume(endcapsAssembly);
    caloDetElem.setPlacement(envelopePhysVol);
    envelopePhysVol.addPhysVolID("system", idDet);

    // Create dummy caloData object for PandoraPFA
    // FIXME: fill calo and layer data information
    auto caloData = new dd4hep::rec::LayeredCalorimeterData;
    caloData->layoutType = dd4hep::rec::LayeredCalorimeterData::EndcapLayout;
    caloDetElem.addExtension<dd4hep::rec::LayeredCalorimeterData>(caloData);

    // save extent information
    caloData->extent[0] = dim.rmin1();
    caloData->extent[1] = dim.rmax1();
    caloData->extent[2] = dim.z_offset() - dim.dz();
    caloData->extent[3] = dim.z_offset() + dim.dz();

    // Set type flags
    dd4hep::xml::setDetectorTypeFlag(xmlDetElem, caloDetElem);

    return caloDetElem;
  }
} // namespace ECalEndcap_ParallelPlateTurbine_o1_v01
} // namespace det

DECLARE_DETELEMENT(ECalEndcap_ParallelPlateTurbine_o1_v01, det::ECalEndcap_ParallelPlateTurbine_o1_v01::createECalEndcapTurbine
)
