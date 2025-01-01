#ifndef DDPlanarDigiAlgorithm_h
#define DDPlanarDigiAlgorithm_h 1

// k4FWCore & EDM4HEP
#include <k4FWCore/Transformer.h>
#include <k4Interface/IUniqueIDGenSvc.h>
#include <edm4hep/SimTrackerHitCollection.h>
#include <edm4hep/TrackerHitCollection.h>
#include <edm4hep/TrackerHitSimTrackerHitLinkCollection.h>
#include <edm4hep/EventHeaderCollection.h>

// Standard
#include <string>
#include <vector>
#include <map>

// DD4HEP
#include <gsl/gsl_rng.h>
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"

// ROOT
#include <TH1F.h>

/** ======= DDPlanarDigiProcessor ========== <br>
 * Creates TrackerHits from SimTrackerHits, smearing them according to the input parameters. 
 * The positions of "digitized" TrackerHits are obtained by gaussian smearing positions
 * of SimTrackerHits perpendicular and along the ladder according to the specified point resolutions. 
 * The geometry of the surface is retreived from DDRec::Surface associated to the hit via cellID.
 * 
 * 
 * <h4>Input collections and prerequisites</h4> 
 * Algorithm requires a collection of SimTrackerHits <br>
 * <h4>Output</h4>
 * Algorithm produces collection of smeared TrackerHits<br>
 * @param SimTrackHitCollectionName The name of input collection of SimTrackerHits <br>
 * (default name VXDCollection) <br>
 * @param TrackerHitCollectionName The name of output collection of smeared TrackerHits <br>
 * (default name VTXTrackerHits) <br>
 * @param ResolutionU resolution in direction of u (in mm) <br>
 * (default value 0.004) <br>
 * @param ResolutionV Resolution in direction of v (in mm) <br>
 * (default value 0.004) <br>
 * @param IsStrip whether the hits are 1 dimensional strip measurements <br>
 * (default value false)
 * @param Ladder_Number_encoded_in_cellID ladder number has been encoded in the cellID <br>
 * (default value false) <br>
 * @param Sub_Detector_ID ID of Sub-Detector using UTIL/ILDConf.h from lcio <br>
 * (default value lcio::ILDDetID::VXD) <br>
 * <br>
 * 
 * @author F.Gaede CERN/DESY, S. Aplin DESY, S. Ferraro
 * @date Dec 2014
 */
class DDPlanarDigiAlgorithm : public k4FWCore::MultiTransformer<std::tuple<edm4hep::TrackerHitCollection, 
                                                                           edm4hep::TrackerHitSimTrackerHitLinkCollection>(
                                                                     const edm4hep::SimTrackerHitCollection&,
                                                                     const edm4hep::EventHeaderCollection&)> {
public:
  
  DDPlanarDigiAlgorithm(const std::string& name, ISvcLocator* svcLoc);
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  StatusCode initialize();
  
  /** Called for every run.
   */
  std::tuple<edm4hep::TrackerHitCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection> operator(const edm4hep::SimTrackerHitCollection& inputSim
                                                                                                     const edm4hep::EventHeaderCollection& evHeader) const; 
  
  /** Called after data processing for clean up.
   */
  StatusCode finalize();
  
protected:
  std::vector<float> initRes{0.0040};
  std::vector<float> initTimeRes{-1.0};
  std::vector<float> initMinTime{-1e9};
  std::vector<float> initMaxTime{1e9};

  Gaudi::Property<std::vector<float>> m_resU{this, "ResolutionU", initRes, "Resolution in direction of u - either one per layer or one for all layers."};
  Gaudi::Property<std::vector<float>> m_resV{this, "ResolutionV", initRes, "Resolution in direction of u - either one per layer or one for all layers."};
  Gaudi::Property<std::vector<float>> m_resT{this, "ResolutionT", initTimeRes, "Resolution of time - either one per layer or one for all layers. If the single entry is negative, disable time smearing."};
    
  Gaudi::Property<bool> m_isStrip{this, "IsStrip", bool(false), "Whether hits are 1D strip hits."};
  Gaudi::Property<std::string> m_subDetName{this, "SubDetectorName", std::string("VXD"), "Name of sub detector."};

  Gaudi::Property<bool> m_forceHitsOntoSurface{this "ForceHitsOntoSurface", bool(false), "Project hits onto the surface in case they are not yet on the surface (default: false)."};
  Gaudi::Property<double> m_minEnergy{this, "MinimumEnergyPerHit", double(0.0), "Minimum Energy (in GeV!) to accept hits, other hits are ignored."};
  
  Gaudi::Property<bool> m_correctTimesForPropagation{this, "CorrectTimesForPropagation", bool(false),  "Correct hit time for the propagation: radial distance/c (default: false)."};
  Gaudi::Property<bool> m_useTimeWindow{this, "UseTimeWindow", bool(false), "Only accept hits with time (after smearing) within the specified time window (default: false)."};
  
  Gaudi::Property<std::vector<float>> m_timeWindow_min{this, "TimeWindowMin", initMinTime, "Minimum time a hit must have after smearing to be accepted [ns] - either one per layer or one for all layers."};
  Gaudi::Property<std::vector<float>> m_timeWindow_max{this, "TimeWindowMax", initMaxTime, "Maximum time a hit must have after smearing to be accepted [ns] - either one per layer or one for all layers."};

  
  SmartIF<IUniqueIDGenSvc> m_idGen;
  gsl_rng* m_rng ;

  const dd4hep::rec::SurfaceMap* m_map ;

  std::vector<TH1F*> m_h ;
  
} ;

#endif



