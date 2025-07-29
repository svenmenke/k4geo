#include "detectorSegmentations/FCCSWEndcapParallelPlateTurbine_k4geo.h"
#include "DD4hep/Detector.h"

namespace dd4hep {
namespace DDSegmentation {

  /// default constructor using an encoding string
  FCCSWEndcapParallelPlateTurbine_k4geo::FCCSWEndcapParallelPlateTurbine_k4geo(const std::string& cellEncoding) : Segmentation(cellEncoding) {

    commonSetup();
  }

  FCCSWEndcapParallelPlateTurbine_k4geo::FCCSWEndcapParallelPlateTurbine_k4geo(const BitFieldCoder* decoder) : Segmentation(decoder) {
    // define type and description

    commonSetup();
  }

  /// initialize variables, etc (needed for either version of the ctor)
  void FCCSWEndcapParallelPlateTurbine_k4geo::commonSetup() {
    _type = "FCCSWEndcapParallelPlateTurbine_k4geo";
    _description = "Parallel-Plate-Turbine-specific segmentation in the global coordinates";

    // register all necessary parameters
    registerParameter("offset_rho", "Offset in rho", m_offsetRho, 0.);

    registerIdentifier("identifier_rho", "Cell ID identifier for rho", m_rhoID, "rho");

    registerParameter("grid_size_rho", "Grid size in rho", m_gridSizeRho, 0.);
    registerParameter("grid_size_z", "Grid size in z", m_gridSizeZ, 0.);
    registerParameter("offset_z", "Offset in z1", m_offsetZ, 0.);
    registerParameter("offset_theta", "Angular offset in theta", m_offsetTheta, 0., SegmentationParameter::AngleUnit,
                      true);
    registerIdentifier("identifier_z", "Cell ID identifier for z", m_zID, "z");
    registerIdentifier("identifier_side", "Cell ID identifier for side", m_sideID, "side");
    registerIdentifier("identifier_group", "Cell ID identifier for group", m_groupID, "group");
    registerIdentifier("identifier_stack", "Cell ID identifier for stack", m_stackID, "stack");
    registerIdentifier("identifier_slice", "Cell ID identifier for slice", m_sliceID, "slice");
    dd4hep::Detector* dd4hepgeo = &(dd4hep::Detector::getInstance());

    try {
      m_groupAngle = dd4hepgeo->constant<double>("GroupAngle");
    } catch (...) {
      std::cout << "GroupAngle not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    
    try {
      m_nGroups = dd4hepgeo->constant<int>("nGroups");
    } catch (...) {
      std::cout << "nGroups not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
    
    try {
      m_numReadoutRhoLayers = dd4hepgeo->constant<int>("ECalEndcapNumReadoutRhoLayers");
    } catch (...) {
      std::cout << "ECalEndcapNumReadoutRhoLayers not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
  
    try {
      m_numReadoutZLayers = dd4hepgeo->constant<int>("ECalEndcapNumReadoutZLayers");
    } catch (...) {
      std::cout << "ECalEndcapNumReadoutZLayers not found in detector metadata, exiting..." << std::endl;
      exit(1);
    }
  }

  /// determine the local position based on the cell ID
  Vector3D FCCSWEndcapParallelPlateTurbine_k4geo::position(const CellID& cID) const {

    double rhoVal = rho(cID);
    double zVal = z(cID);
    double phiVal = phi(cID);
    Vector3D pos = PositionRhoZPhi(rhoVal, zVal, phiVal);
    // account for the fact that the -z endcap is mirrored wrt to the +z one
    if (pos.Z < 0.)
      pos.Y = -pos.Y;

    return pos;
  }

  /// determine the cell ID based on the position
  CellID FCCSWEndcapParallelPlateTurbine_k4geo::cellID(const Vector3D& /* localPosition */, const Vector3D& globalPosition,
                                          const VolumeID& vID) const {
    CellID cID = vID;
    CellID iGroup = _decoder->get(cID, m_groupID);
    CellID iStack = _decoder->get(cID, m_stackID);
    CellID iSlice = _decoder->get(cID, m_sliceID);

    double lRho = rhoFromXYZ(globalPosition);
    int iRho = positionToBin(lRho, m_gridSizeRho, m_offsetRho + m_gridSizeRho / 2.);
    if (iRho < 0) {
      iRho = 0;
    }
    if (iRho >= m_numReadoutRhoLayers) {
      iRho = m_numReadoutRhoLayers - 1;
    }
    _decoder->set(cID, m_rhoID, iRho);

    double lZ = TMath::Abs(globalPosition.Z);
    int iZ = positionToBin(lZ, m_gridSizeZ, m_offsetZ + m_gridSizeZ / 2.);
    if (iZ < 0) {
      iZ = 0;
    }
    if (iZ >= m_numReadoutZLayers) {
      iZ = m_numReadoutZLayers - 1;
    }
    _decoder->set(cID, m_zID, iZ);

    return cID;
  }

  /// determine rho based on the cell ID
  double FCCSWEndcapParallelPlateTurbine_k4geo::rho(const CellID& cID) const {
    CellID rhoValue = _decoder->get(cID, m_rhoID);

    return binToPosition(rhoValue, m_gridSizeRho, m_offsetRho) + m_gridSizeRho / 2.;
  }

  /// determine the azimuthal angle phi based on the cell ID
  double FCCSWEndcapParallelPlateTurbine_k4geo::phi(const CellID& cID) const {
    CellID iGroup = _decoder->get(cID, m_groupID);
    CellID iStack = _decoder->get(cID, m_stackID);

    double phiCent = twopi * (iGroup + 0.5) / (m_nGroups);
    double rhoLoc = rho(cID);

    double zdepth = m_numReadoutZLayers * m_gridSizeZ;

    double zLoc = TMath::Abs(z(cID)) - m_offsetZ - zdepth / 2;
    double x = zLoc / TMath::Tan(m_groupAngle);
    double y = TMath::Sqrt(rhoLoc * rhoLoc - x * x);
    // rotate about z axis by phiCent
    double xprime = x * TMath::Cos(phiCent) + y * TMath::Sin(phiCent);
    double yprime = y * TMath::Cos(phiCent) - x * TMath::Sin(phiCent);

    return TMath::ATan2(xprime, yprime);
  }

  /// determine local x in plane of blade based on the cell ID
  double FCCSWEndcapParallelPlateTurbine_k4geo::z(const CellID& cID) const {
    CellID zValue = _decoder->get(cID, m_zID);
    CellID sideValue = _decoder->get(cID, m_sideID);
    return ((long long int)sideValue) *
           (binToPosition(zValue, m_gridSizeZ, m_offsetZ) + m_gridSizeZ / 2.);
  }

} // namespace DDSegmentation
} // namespace dd4hep
