/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */

/*
Evolved version of TPCDigi that provides additional functionality to deal with background. Couple to the Mokka Sensitive Detector Driver TPCSD03.cc

SJA:FIXME: Still needs to be tidied up for production release.

Three cases can be consider in the treatment of SimTrackerHits
i)   A clean isolated hit; this will be smeared according to the parametric point resolution and converted to a TrackerHit
ii)  Two or Three hits which are considered to be closer than the double hit resolution and which therefore cannot be viewed as seperable hits. These will be merged and be assigned a large associated measurement error.
iii) A continuous set of hits within one pad row which cannot be resolved as single hits, these are condidered to be charaterisable as background hits created by extremely low pt charged particles (pt < 10MeV) and therefore are removed from the hit collection.

The Driver has been modified to take an additional collection of SimTrackerHits which are produced by the Mokka TPC Sensitive Driver TPCSD03.cc. These hits are produced for particles which have very low pt and often do not move outside of the dimensions of a single pad row. These hits need to be treated differently as they do not cross any geometric boundaries in a Padrow based TPC Geometry. This negates the need to voxalise the TPC in Geant4 which has proved in the past to be prohibitive in terms of processing time due to the vastly increased number of geometric volumes. 

Steve Aplin 26 June 2009 (DESY)

*/

#ifndef DDTPCDigiAlgorithm_h
#define DDTPCDigiAlgorithm_h 1

// Standard
#include <string>
#include <gsl/gsl_rng.h>
#include <vector>
#include <map>

// k4FWCore & EDM4HEP
#include <k4FWCore/Transformer.h>
#include <edm4hep/TrackerHitSimTrackerHitLinkCollection.h>
#include <edm4hep/MCParticle.h>
#include <edm4hep/SimTrackerHit.h>
#include <edm4hep/MutableTrackerHitPlane.h>
#include <UTIL/CellIDEncoder.h>


#include "DDRec/DetectorData.h"
#include "CLHEP/Vector/TwoVector.h"

class Voxel_tpc;

class TPCModularEndplate;

/** ====== DDTPCDigiProcessor ====== <br>
 *
 * This Processor depends on Circle.h from MarlinUtil
 * 
 * Caution: This digitiser presently does not process space-point like SimTrackerHits which have been flagged with CellIDs set to the negetive row number. This must be implemented in future. 
 *Produces TPC TrackerHit collection from SimTrackerHit collection, smeared in r-phi and z. 
 * Double hits are identified but are currently not added to the collection. This may be change 
 * at a later date when criteria for their seperation is defined. The resolutions are defined in 
 * the GEAR stearing file. 
 *
 * Resolution in r-phi is calculated according to the formular <br>
 * sigma_{point}^2 = sigma_0^2 + Cd^2/N_{eff} * L_{drift}
 * Cd^2/N_{eff}} = 25^2/(22/sin(theta)*h/6mm)
 * Cd = 25 ( microns / cm^(1/2) )
 * (this is for B=4T, h is the pad height = pad-row pitch in mm,
 * theta is the polar angle)       
 *
 * At the moment resolution in z assumed to be independent of drift length. <br>
 *
 * The type of TPC TrackerHit is set to 500 via method TrackerHitImpl::setType(int type) <br>
 * <h4>Input collections and prerequisites</h4> 
 * Processor requires collections of SimTrackerHits in TPC <br>
 * <h4>Output</h4>
 * Processor produces collection of digitized TrackerHits in TPC <br>
 * @param CollectionName The name of input SimTrackerHit collection <br>
 * (default name STpc01_TPC)
 * @param RejectCellID0 Whether or not to reject SimTrackerHits with Cell ID 0. Mokka drivers
 * TPC00-TPC03 encoded the pad row number in the cell ID, which should always be non-zero anyway.
 * Drivers TPC04 and TPC05 do not simulate pad rows and thus have the cell ID set to zero for all hits.
 * You will need to set RejectCellID0 to 0 in order to use this processor with these drivers, but note
 * that the implications for track reconstruction are not strictly defined. Mokka driver TPC06 uses
 * a mixed approach with one hit per pad row having non-zero cell ID, extra hits having 0 cell ID.
 * Typically, unless you use TPC04 or TPC05, you should not touch this parameter. <br>
 * (default value 1)
 * @param TPCTrackerHitsCol The name of output collection of TrackerHits <br>
 * (default name TPCTrackerHits) <br>
 * <br>
 * @authors F.Gaede, S. Aplin, DESY and A.Raspereza, MPI and S. Ferraro
 *
 * 06/2017 FG: replace Gear with DDRec
 *         FG: remove hits that are inside the module gaps on the endplate (see parameters TPCEndPlateModulexxxx)
 *             (set TPCEndPlateModuleGapPhi=0 to not remove any hits)
 * 
 * Changed 7/9/07 so that the const and diffusion resolution terms are taken as processor parameters rather than the gear file.
 * The parameters _pixZ and pixRP were also changed from gear parameters to processor parameters
 * clare.lynch@bristol.ac.uk
 *
 */
class DDTPCDigiAlgorithm : public k4FWCore::MultiTransformer<std::tuple<
	edm4hep::TrackerHitPlaneColection,
	edm4hep::TrackerHitSimTrackerHitLinkCollection>(
	const std::vector<const edm4hep::SimTrackerHitCollection*> &)>{
  
public:
  
  DDTPCDigiAlgorithm();

  ~DDTPCDigiAlgorithm();
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  StatusCode initialize();
  
  /** Called for every event - the working horse.
   */
  std::tuple<edm4hep::TrackerHitPlaneColection,
             edm4hep::TrackerHitSimTrackerHitLinkCollection> operator(
	     const std::vector<const edm4hep::SimTrackerHitCollection*>& inputCols) const; 
  
  /** Called after data processing for clean up.
   */
  StatusCode finalize();
  
  void writeVoxelToHit( Voxel_tpc* aVoxel, UTIL::LCRelationNavigator& hitSimHitNav) ;
  void writeMergedVoxelsToHit( std::vector <Voxel_tpc*>* hitList, UTIL::LCRelationNavigator& hitSimHitNav ) ;
  void plotHelixHitResidual(edm4hep::MCParticle *mcp, CLHEP::Hep3Vector *thisPointRPhi);
  double getPadPhi( CLHEP::Hep3Vector* thisPointRPhi, CLHEP::Hep3Vector* firstPointRPhi, CLHEP::Hep3Vector* middlePointRPhi, CLHEP::Hep3Vector* lastPointRPhi);
  double getPadTheta( CLHEP::Hep3Vector* firstPointRPhi, CLHEP::Hep3Vector* middlePointRPhi, CLHEP::Hep3Vector* lastPointRPhi );

protected:
  Gaudi::Property<bool> m_use_raw_hits_to_store_simhit_pointer{this, "UseRawHitsToStoreSimhitPointer", bool(false), "Store the pointer to the SimTrackerHits in RawHits (deprecated)."};
  
  Gaudi::Property<int> m_rejectCellID0{this, "RejectCellID0", (int)1, "Whether or not to use hits without proper cell ID (pad row)."};
  Gaudi::Propert<float> m_padWidth{};

  EVENT::MCParticle* m_mcp{};
  EVENT::MCParticle* m_previousMCP{};
  EVENT::MCParticle* m_nextMCP{};
  EVENT::MCParticle* m_nMinus2MCP{};
  EVENT::MCParticle* m_nPlus2MCP{};   

  SimTrackerHit* m_SimTHit{};
  SimTrackerHit* m_previousSimTHit{};
  SimTrackerHit* m_nextSimTHit{};
  SimTrackerHit* m_nPlus2SimHit{};
  SimTrackerHit* m_nMinus2SimHit{};

  // gsl random number generator
  gsl_rng * m_random {};

  Gaudi::Property<float> m_pointResoRPhi0{this, "PointResolutionRPhi", (float)0.050, "R-Phi Resolution constant in TPC."}; // Coefficient for RPhi point res independant of drift length 
  Gaudi::Property<float> m_pointResoPadPhi{this, "PointResolutionPadPhi", (float)0.900, "Pad Phi Resolution constant in TPC."}; // Coefficient for the point res dependance on relative phi angle to the pad verticle 
  Gaudi::Property<float> m_diffRPhi{this, "DiffusionCoeffRPhi", (float)0.025, "R-Phi Diffusion Coefficent in TPC."}; // Coefficient for the rphi point res dependance on diffusion 
  Gaudi::Property<int>   m_nEff{this, "N_eff", (int)22, "Number of Effective electrons per pad in TPC"}; // number of effective electrons 

  Gaudi::Property<float> m_pointResoZ0{this, "PointResolutionZ", (float)0.4, "TPC Z Resolution Coefficent independent of diffusion."}; // Coefficient Z point res independant of drift length 
  Gaudi::Property<float> m_diffZ{this, "DiffusionCoeffZ", (float)0.08,  "Z Diffusion Coefficent in TPC."}; // Coefficient for the Z point res dependance on diffusion 

  Gaudi::Property<float> m_binningZ{this, "HitSortingBinningZ", (float)5.0, "Defines spatial slice in Z."};
  Gaudi::Property<float> m_binningRPhi{this, "HitSortingBinningRPhi", (float)2.0, "Defines spatial slice in RP."};
  Gaudi::Property<float> m_doubleHitResZ{this, "DoubleHitResolutionZ", (float)5.0,  "Defines the minimum distance for two seperable hits in Z."};
  Gaudi::Property<float> m_doubleHitResRPhi{this, "DoubleHitResolutionRPhi", (float)2.0, "Defines the minimum distance for two seperable hits in RPhi."};
  Gaudi::Property<int>   m_maxMerge{this, "MaxClusterSizeForMerge", (int)3, "Defines the maximum number of adjacent hits which can be merged."};

  int m_nRechits{};

  std::vector< std::vector<Voxel_tpc *> > m_tpcRowHits{};
  std::map< Voxel_tpc *,SimTrackerHit *> m_tpcHitMap{};
  std::vector<float> m_length{};
  int lenpos{};

  edm4hep::TrackerHitPlaneCollection* m_trkhitVec{};
  CellIDEncoder<edm4hep::MutableTrackerHitPlane>* m_cellid_encoder {};

  int  m_NSimTPCHits{};
  int  m_NBackgroundSimTPCHits{};
  int  m_NPhysicsSimTPCHits{};
  int  m_NPhysicsAbove02GeVSimTPCHits{};
  int  m_NPhysicsAbove1GeVSimTPCHits{};
  int  m_NRecTPCHits{};
  
  int  m_NLostPhysicsTPCHits{};
  int  m_NLostPhysicsAbove02GeVPtTPCHits{};
  int  m_NLostPhysicsAbove1GeVPtTPCHits{};
  int  m_NRevomedHits{};

  const dd4hep::rec::FixedPadSizeTPCData*  m_tpc{};
  double m_bField{};
  
  std::vector<int>   tpcEPModNumExample = { 14, 18, 23, 28, 32, 37, 42, 46 };
  std::vector<float> tpcEPModPhi0Example = { 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07 }; 

  Gaudi::Property<std::vector<int>>    m_tpcEndPlateModuleNumbers{this, "TPCEndPlateModuleNumbers", tpcEPModNumExample, "Number of modules in the rings of the TPC endplate."};
  Gaudi::Property<std::vector<float>>  m_tpcEndPlateModulePhi0s{this, "TPCEndPlateModulePhi0s", tpcEPModPhi0Example, "Phi0s of modules in the rings of the TPC endplate."};
  Gaudi::Property<float>    m_tpcEndPlateModuleGapPhi{this, "TPCEndPlateModuleGapPhi", (float)1., "Gap size in mm of the gaps between the endplace modules in Phi."};
  Gaudi::Property<float>    m_tpcEndPlateModuleGapR{this, "TPCEndPlateModuleGapR", (float)1., "Gap size in mm of the gaps between the endplace modules in R."};

  TPCModularEndplate* m_tpcEP{} ;

  // Histogram
  TH1D * m_phiDiffHisto{};
  TH1D * m_thetaDiffHisto{};
  TH1D * m_phiRelHisto{};
  TH1D * m_thetaRelHisto{};

  TH1D * m_phiDistHisto{};
  TH1D * m_rPhiPullHisto{};
  TH1D * m_rPhiDiffHisto{};
  TH1D * m_zDiffHisto{};
  TH1D * m_zPullHisto{};
  TH2D * m_zSigmaVsZHisto{};
  TH1D * m_zSigmaHisto{};
  TH1D * m_rPhiSigmaHisto{};
  TH1D * m_radiusCheckHisto{};
  TH1D * m_ResidualsRPhiHisto{};

  TH1D * m_NSimTPCHitsHisto{};
  TH1D * m_NBackgroundSimTPCHitsHisto{};
  TH1D * m_NPhysicsSimTPCHitsHisto{};
  TH1D * m_NPhysicsAbove02GeVSimTPCHitsHisto{};
  TH1D * m_NPhysicsAbove1GeVSimTPCHitsHisto{};
  TH1D * m_NRecTPCHitsHisto{};

  TH1D * m_NLostPhysicsTPCHitsHisto{};
  TH1D * m_NLostPhysicsAbove02GeVPtTPCHitsHisto{};
  TH1D * m_NLostPhysicsAbove1GeVPtTPCHitsHisto{};
  TH1D * m_NRevomedHitsHisto{};

  TH1D * m_NKeptPhysicsTPCHitsHistoPercent{};
  TH1D * m_NKeptPhysicsAbove02GeVPtTPCHitsHistoPercent{};
  TH1D * m_NKeptPhysicsAbove1GeVPtTPCHitsHistoPercent{};

#endif
} ;
#endif



