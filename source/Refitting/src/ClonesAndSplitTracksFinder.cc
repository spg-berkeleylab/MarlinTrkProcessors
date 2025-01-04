#include "ClonesAndSplitTracksFinder.h"
#include "HitsSorterAndDebugger.h"

#include <marlin/Exceptions.h>
#include <marlin/Global.h>
#include <marlin/VerbosityLevels.h>

/*
#include <MarlinTrk/Factory.h>
#include <MarlinTrk/IMarlinTrack.h>
#include <MarlinTrk/MarlinTrkUtils.h>
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "MarlinTrk/MarlinTrkDiagnostics.h"

#include <marlinutil/GeometryUtil.h>
*/

// edm4hep
#include <edm4hep/TrackerHit.h>
#include <edm4hep/MutableTrack.h>
#include <edm4hep/MutableTrackState.h>
#include "BitField64.h"
#include <DD4hep/Detector.h>

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DDRec/SurfaceManager.h"

#include <algorithm>


using namespace std;

DECLARE_COMPONENT(ClonesAndSplitTracksFinder)

ClonesAndSplitTracksFinder::ClonesAndSplitTracksFinder(const std::string& name, ISvcLocator* svcLoc) : Transformer(name, svcLoc,
              { KeyValues("InputTrackCollectionName", {"SiTracks"}) },
              { KeyValues("OutputTrackCollectionName", {"SiTracksMerged"}) }) {}

StatusCode ClonesAndSplitTracksFinder::initialize() {
  /*
  _trksystem = MarlinTrk::Factory::createMarlinTrkSystem("DDKalTest", nullptr, "");

  _magneticField = MarlinUtil::getBzAtOrigin();
  ///////////////////////////////

  _trksystem->setOption(MarlinTrk::IMarlinTrkSystem::CFG::useQMS, _MSOn);
  _trksystem->setOption(MarlinTrk::IMarlinTrkSystem::CFG::usedEdx, _ElossOn);
  _trksystem->setOption(MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing, _SmoothOn);
  _trksystem->init();

  */
  DD4hep::Detector lcdd = dd4hep::Detector::getInstance();
  const double position[3] = {0, 0, 0};
  double magField[3] = {0, 0, 0};
  lcdd.field().magneticField(position, magField);
  m_magneticField = magField[2] / dd4hep::tesla;
  m_encoder = std::make_shared<BitField64>( "system:5,side:-2,layer:6,module:11,sensor:8" );

  // Put default values for track fitting
  m_initialTrackError_d0    = 1.e6;
  m_initialTrackError_phi0  = 1.e2;
  m_initialTrackError_omega = 1.e-4;
  m_initialTrackError_z0    = 1.e6;
  m_initialTrackError_tanL  = 1.e2;
  m_maxChi2perHit           = 1.e2;

  return StatusCode::SUCCESS;
}

edm4hep::TrackCollection ClonesAndSplitTracksFinder::operator(const edm4hep::TrackCollection& input_track_col) const{
  // set the correct configuration for the tracking system for this event
  /*
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useQMS>       mson(_trksystem, _MSOn);
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::usedEdx>      elosson(_trksystem, _ElossOn);
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing> smoothon(_trksystem, _SmoothOn);
  */

  const int nTracks = input_track_col->getNumberOfElements();
  debug() << " >> ClonesAndSplitTracksFinder starts with " << nTracks << " tracks." << endmsg;

  // establish the track collection that will be created
  edm4hep::TrackCollection outputCollection;
  outputCollection.setSubsetCollection();
  m_encoder->reset();
  // if we want to point back to the hits we need to set the flag
  //LCFlagImpl trkFlag(0);
  //trkFlag.setBit(LCIO::TRBIT_HITS);
  //trackVec->setFlag(trkFlag.getFlag());

  //------------
  // FIRST STEP: REMOVE CLONES
  //------------

  std::vector<edm4hep::Track> tracksWithoutClones;
  removeClones(tracksWithoutClones, input_track_col);
  const int ntracksWithoutClones = tracksWithoutClones.size();
  debug() << " >> ClonesAndSplitTracksFinder found " << ntracksWithoutClones << " tracks without clones."
                        << endmsg;

  if (m_mergeSplitTracks && ntracksWithoutClones > 1) {
    debug() << " Try to merge tracks ..." << endmsg;

    //------------
    // SECOND STEP: MERGE TRACKS
    //------------

    mergeSplitTracks(outputCollection, input_track_col, tracksWithoutClones);
  } else {
    debug() << " Not even try to merge tracks ..." << endmsg;
    for (UInt_t iTrk = 0; iTrk < tracksWithoutClones.size(); iTrk++) {
      outputCollection.push_back(tracksWithoutClones.at(iTrk));
    }
  }

  return outputCollection;
}


StatusCode ClonesAndSplitTracksFinder::finalize() { return StatusCode::SUCCESS }


// Function to check if two KDtracks contain several hits in common
int ClonesAndSplitTracksFinder::overlappingHits(const edm4hep::Track& track1, const edm4hep::Track& track2) {
  int nHitsInCommon = 0;

  for (size_t itrackHit = 0; itrackHit < track1.trackerHits_size(); ++itrackHit) {
    auto it = std::find_if(track2.trackerHits_begin(), track2.trackerHits_end(), [&](const auto& hit2){
      return track1.getTrackerHits(itrackHit) == hits;
    });
    if (it != track2.trackerHits_end()) {
      nHitsInCommon++; // If it does overlap, increment counter
    }
  }

  return nHitsInCommon;
}


void ClonesAndSplitTracksFinder::removeClones(std::vector<edm4hep::Track>& tracksWithoutClones, const edm4hep::TrackCollection& input_track_col) {
  debug() << "ClonesAndSplitTracksFinder::removeClones " << endmsg;
  const int nTracks = input_track_col.size();

  // loop over the input tracks

  std::multimap<int, std::pair<int, edm4hep::Track*>> candidateClones;

  for (int iTrack = 0; iTrack < nTracks; ++iTrack) {  // first loop over tracks
    int countClones        = 0;
    edm4hep::Track track_i = input_track_col.at(iTrack);

    for (int jTrack = 0; jTrack < nTracks; ++jTrack) {  // second loop over tracks

      edm4hep::Track track_j = input_track_col.at(jTrack);
      if (track_i != track_j) {  // track1 != track2

        const unsigned int nOverlappingHits = overlappingHits(track_i, track_j);
        if (nOverlappingHits >= 2) {  // clones
          countClones++;
          edm4hep::Track* bestTrack;
          bestInClones(track_i, track_j, nOverlappingHits, bestTrack);
          candidateClones.insert(make_pair(iTrack, make_pair(jTrack, bestTrack)));
        } else {
          continue;
        }
      }

    }  // end second track loop

    if (countClones == 0) {
      tracksWithoutClones.push_back(track_i);
    }

  }  // end first track loop

  filterClonesAndMergedTracks(candidateClones, input_track_col, tracksWithoutClones, true);
}

void ClonesAndSplitTracksFinder::mergeSplitTracks(edm4hep::TrackCollection& outputCollection, const edm4hep::TrackCollection& input_track_col,
                                                  vector<edm4hep::Track>& tracksWithoutClones) {
  debug() << "ClonesAndSplitTracksFinder::mergeSplitTracks " << endmsg;

  std::multimap<int, std::pair<int, edm4hep::Track*>> mergingCandidates;
  std::set<int> iter_duplicates;

  for (UInt_t iTrack = 0; iTrack < tracksWithoutClones.size(); ++iTrack) {
    int    countMergingPartners  = 0;
    bool   toBeMerged            = false;
    edm4hep::Track track_i       = tracksWithoutClones.at(iTrack);

    double pt_i    = 0.3 * m_magneticField / (fabs(track_i.getOmega()) * 1000.);
    double theta_i = (M_PI / 2 - std::atan(track_i.getTanLambda())) * 180. / M_PI;
    double phi_i   = track_i.getPhi() * 180. / M_PI;

    //Merge only tracks with min pt
    //Try to avoid merging loopers for now
    if (pt_i < m_minPt) {
      debug() << " Track #" << iTrack << ": pt = " << pt_i << ", theta = " << theta_i << ", phi = " << phi_i
              << "\n Track #" << iTrack << " does not fulfil min pt requirement."
              << "\n TRACK STORED" << endmsg;

      outputCollection.push_back(track_i);
      continue;
    }

    for (UInt_t jTrack = iTrack + 1; jTrack < tracksWithoutClones.size(); ++jTrack) {
      edm4hep::Track track_j = tracksWithoutClones.at(jTrack);
      bool isCloseInTheta    = false, isCloseInPhi = false, isCloseInPt = false;

      if (track_j != track_i) {
        double pt_j    = 0.3 * m_magneticField / (fabs(track_j.getOmega() * 1000.));
        double theta_j = (M_PI / 2 - std::atan(track_j.getTanLambda())) * 180. / M_PI;
        double phi_j   = track_j.getPhi() * 180. / M_PI;
        debug() << " Track #" << iTrack << ": pt = " << pt_i << ", theta = " << theta_i << ", phi = " << phi_i
                << "\n Track #" << jTrack << ": pt = " << pt_j << ", theta = " << theta_j << ", phi = " << phi_j
                << endmsg;

        if (pt_j < _minPt) {
          debug() << " Track #" << jTrack << " does not fulfil min pt requirement. Skip. " << endmsg;
          continue;
        }

        double significanceTanLambda = calculateSignificanceTanLambda(track_i, track_j);
        double significancePhi       = calculateSignificancePhi(track_i, track_j);
        double significancePt        = calculateSignificancePt(track_i, track_j);

        debug() << " -> tanLambda significance = " << significanceTanLambda << " with cut at "
                              << m_maxSignificanceTheta << endmsg;
        if (significanceTanLambda < m_maxSignificanceTheta) {
          isCloseInTheta = true;
          debug() << " Tracks are close in theta " << endmsg;
        }

        debug() << " -> phi significance = " << significancePhi << " with cut at " << m_maxSignificancePhi << endmsg;
         
        if (significancePhi < m_maxSignificancePhi) {
          isCloseInPhi = true;
          debug() << " Tracks are close in phi " << endmsg;
        }
        debug() << " -> pt significance = " << significancePt << " with cut at " << m_maxSignificancePt << endmsg;
        
        if (significancePt < m_maxSignificancePt) {
          isCloseInPt = true;
          debug() << " Tracks are close in pt  " << endmsg;
        }

        debug() << " Track #" << iTrack << ": " << endsmg;
        debug() << " Track #" << jTrack << ": " << endmsg;

        toBeMerged = isCloseInTheta && isCloseInPhi && isCloseInPt;

        if (toBeMerged) {  // merging, refitting, storing in a container of mergingCandidates (multimap <*track1, pair<*track2,*trackMerged>>)
          edm4hep::Track* trkPtr = nullptr;
          mergeAndFit(track_i, track_j, trkPtr);
          if (not lcioTrkPtr) {
            continue;
          }
          mergingCandidates.insert(make_pair(iTrack, make_pair(jTrack, trkPtr)));
          countMergingPartners++;
          iter_duplicates.insert(iTrack);
          iter_duplicates.insert(jTrack);
        } else {  // no merging conditions met
          continue;
        }
      }

    }  // end loop on jTracks

    // Track was already found as duplicate
    const bool is_in = iter_duplicates.find(iTrack) != iter_duplicates.end();
    if (countMergingPartners == 0 && !is_in) {  // if track_i has no merging partner, store it in the output vec
      debug() << " Track #" << iTrack << " has no merging partners, so TRACK STORED." << endmsg;

      outputCollection.push_back(track_i);

    } else {
      debug() << " TRACK NOT STORED" << endmsg;
    }
    if (countMergingPartners != 0)
      debug() << " possible merging partners for track #" << iTrack << " are = " << countMergingPartners << endmsg;
  
  }  // end loop on iTracks

  std::vector<edm4hep::Track> finalTracks;
  filterClonesAndMergedTracks(mergingCandidates, input_track_col, finalTracks, false);

  for (UInt_t iTrk = 0; iTrk < finalTracks.size(); iTrk++) {
    debug() << " TRACK STORED" << endmsg;
    outputCollection.push_back(finalTracks.at(iTrk));
  }
}

double ClonesAndSplitTracksFinder::calculateSignificancePt(const edm4hep::Track first, const edm4hep::Track second) {
  float omegaFirst  =  first.getOmega();
  float omegaSecond = second.getOmega();

  double ptFirst  = 0.3 * m_magneticField / (fabs( first.getOmega() * 1000.));
  double ptSecond = 0.3 * m_magneticField / (fabs(second.getOmega() * 1000.));

  const float sigmaPOverPFirst  = sqrt( first.getCovMatrix()[5]) / fabs(omegaFirst);
  const float sigmaPOverPSecond = sqrt(second.getCovMatrix()[5]) / fabs(omegaSecond);
  const float sigmaPtFirst      = ptFirst  * sigmaPOverPFirst;
  const float sigmaPtSecond     = ptSecond * sigmaPOverPSecond;

  const double significance = calculateSignificance(ptFirst, ptSecond, sigmaPtFirst, sigmaPtSecond);

  return significance;
}

double ClonesAndSplitTracksFinder::calculateSignificancePhi(const edm4hep::Track first, const edm4hep::Track second) {
  float phiFirst  =  first.getPhi();
  float phiSecond = second.getPhi();
  float deltaPhi  = (M_PI - std::abs(std::abs(phiFirst - phiSecond) - M_PI));

  const float sigmaPhiFirst  = sqrt( first.getCovMatrix()[2]);
  const float sigmaPhiSecond = sqrt(second.getCovMatrix()[2]);

  const double significance = calculateSignificance(deltaPhi, 0.0, sigmaPhiFirst, sigmaPhiSecond);

  return significance;
}

double ClonesAndSplitTracksFinder::calculateSignificanceTanLambda(const edm4hep::Track first, const edm4hep::Track second) {
  float tanLambdaFirst  =  first.getTanLambda();
  float tanLambdaSecond = second.getTanLambda();

  const float sigmaTanLambdaFirst  = sqrt( first.getCovMatrix()[14]);
  const float sigmaTanLambdaSecond = sqrt(second.getCovMatrix()[14]);

  const double significance =
      calculateSignificance(tanLambdaFirst, tanLambdaSecond, sigmaTanLambdaFirst, sigmaTanLambdaSecond);

  return significance;
}

double ClonesAndSplitTracksFinder::calculateSignificance(const double firstPar, const double secondPar,
                                                         const double firstPar_sigma, const double secondPar_sigma) {
  const float delta      = fabs(firstPar - secondPar);
  const float sigmaDelta = sqrt(firstPar_sigma * firstPar_sigma + secondPar_sigma * secondPar_sigma);

  return delta / sigmaDelta;
}

void ClonesAndSplitTracksFinder::filterClonesAndMergedTracks(std::multimap<int, std::pair<int, edm4hep::Track*>>& candidates,
                                                             const edm4hep::TrackCollection& inputTracks, std::vector<edm4hep::Track> finalTracks,
                                                             bool clones) {
  std::vector<podio::RelationRange<edm4hep::TrackerHit>> savedHitVec;

  for (const auto& iter : candidates) {
    int    track_a_id                = iter.first;
    int    track_b_id                = iter.second.first;
    edm4hep::Track* track_final      = iter.second.second;
    int    countConnections          = candidates.count(track_a_id);
    bool   multiConnection           = (countConnections > 1);

    if (!multiConnection) {  // if only 1 connection

      if (clones) {  // clones: compare the track pointers
        auto it_trk = std::find(finalTracks.begin(), finalTracks.end(), track_final);
        if (it_trk != finalTracks.end()) {  // if the track is already there, do nothing
          continue;
        }
        finalTracks.push_back(track_final);
      } else {  // mergeable tracks: compare the sets of tracker hits

        podio::RelationRange<edm4hep::TrackerHit> track_final_hits = track_final.getTrackerHits();
        bool toBeSaved = true;

        for (const auto& hitsVec : savedHitVec) {
          if (std::equal(hitsVec.begin(), hitsVec.end(), track_final_hits.begin())) {
            toBeSaved = false;
            break;
          }
        }

        if (toBeSaved) {
          savedHitVec.push_back(track_final_hits);
          finalTracks.push_back(track_final);
        } else {
          delete track_final;
        }
      }

    } else {  // if more than 1 connection, clones and mergeable tracks have to be treated a little different

      if (clones) {  // clones

        //look at the elements with equal range. If their bestTrack is the same, store it (if not already in). If their bestTrack is different, don't store it
        auto ret = candidates.equal_range(
            track_a_id);  //a std::pair of iterators on the multimap [ std::pair<std::multimap<Track*,std::pair<Track*,Track*>>::iterator, std::multimap<Track*,std::pair<Track*,Track*>>::iterator> ]
        std::vector<edm4hep::Track*> bestTracksMultiConnections;
        for (std::multimap<int, std::pair<int, edm4hep::Track*>>::iterator it = ret.first; it != ret.second; ++it) {
          edm4hep::Track* track_best = it->second.second;
          bestTracksMultiConnections.push_back(track_best);
        }
        if (std::adjacent_find(bestTracksMultiConnections.begin(), bestTracksMultiConnections.end(),
                               std::not_equal_to<edm4hep::Track*>()) ==
            bestTracksMultiConnections.end()) {  //one best track with the same track key
          auto it_trk = find(finalTracks.begin(), finalTracks.end(), bestTracksMultiConnections.at(0));
          if (it_trk != finalTracks.end()) {  // if the track is already there, do nothing
            continue;
          }
          finalTracks.push_back(bestTracksMultiConnections.at(0));

        } else {  //multiple best tracks with the same track key
          continue;
        }

      }  // end of clones

      else {                 // mergeable tracks -- at the moment they are all stored (very rare anyways)
        delete track_final;  //not using the mergedTracks, so delete it

        edm4hep::Track track_a = inputTracks.at(track_a_id);
        edm4hep::Track track_b = inputTracks.at(track_b_id);

        auto trk1 = find(finalTracks.begin(), outputCollection.end(), track_a);

        if (trk1 != finalTracks.end()) {  // if the track1 is already there
          continue;
        }
        // otherwise store the two tracks
        finalTracks.push_back(track_a);
        finalTracks.push_back(track_b);

      }  // end of mergeable tracks
    }
  }
}

void ClonesAndSplitTracksFinder::mergeAndFit(edm4hep::Track track_i, edm4hep::Track track_j, edm4hep::Track*& trkPtr) {
  debug() << "ClonesAndSplitTracksFinder::mergeAndFit " << endmsg;
  podio::RelationRange<edm4hep::TrackerHit> trkHits_i = track_i.getTrackerHits();
  podio::RelationRange<edm4hep::TrackerHit> trkHits_j = track_j.getTrackerHits();

  podio::RelationRange<edm4hep::TrackerHit> trkHits;
  for (UInt_t iHits = 0; iHits < trkHits_i.size(); iHits++) {
    trkHits.push_back(trkHits_i.at(iHits));
  }
  //Remove common hits while filling for the second track
  for (UInt_t jHits = 0; jHits < trkHits_j.size(); jHits++) {
    if (std::find(trkHits.begin(), trkHits.end(), trkHits_j.at(jHits)) != trkHits.end()) {
      debug() << " This hit is already in the track" << endmsg;
      continue;
    } else {
      trkHits.push_back(trkHits_j.at(jHits));
    }
  }
  std::sort(trkHits.begin(), trkHits.end(), sort_by_radius);
  debug() << " Hits in track to be merged: " << endmsg;

  auto mergedTrack = std::unique_ptr<edm4hep::MutableTrack>(new edm4hep::MutableTrack);

  auto marlin_trk = std::unique_ptr<MarlinTrk::IMarlinTrack>(_trksystem->createTrack());

  // Make an initial covariance matrix with very broad default values
  EVENT::FloatVec covMatrix(15, 0);            // Size 15, filled with 0s
  covMatrix[0]  = (_initialTrackError_d0);     //sigma_d0^2
  covMatrix[2]  = (_initialTrackError_phi0);   //sigma_phi0^2
  covMatrix[5]  = (_initialTrackError_omega);  //sigma_omega^2
  covMatrix[9]  = (_initialTrackError_z0);     //sigma_z0^2
  covMatrix[14] = (_initialTrackError_tanL);   //sigma_tanl^2

  const bool direction = _extrapolateForward ? MarlinTrk::IMarlinTrack::forward : MarlinTrk::IMarlinTrack::backward;

  int fit_status = MarlinTrk::createFinalisedLCIOTrack(marlin_trk.get(), trkHits, mergedTrack.get(), direction, covMatrix,
                                                       _magneticField, _maxChi2perHit);

  if (fit_status != 0) {
    streamlog_out(DEBUG4) << "Fit failed with error status " << fit_status << std::endl;
    return;
  }
  streamlog_out(DEBUG8) << " >> Fit not failed ! " << std::endl;

  // fit finished - get hits in the fit
  std::vector<std::pair<EVENT::TrackerHit*, double>> hits_in_fit;
  std::vector<std::pair<EVENT::TrackerHit*, double>> outliers;

  // remember the hits are ordered in the order in which they were fitted

  marlin_trk->getHitsInFit(hits_in_fit);
  if (hits_in_fit.size() < 3) {
    streamlog_out(DEBUG4) << "Less than 3 hits in fit: Track discarded. Number of hits =  " << trkHits.size() << std::endl;
    return;
  }

  marlin_trk->getOutliers(outliers);

  std::vector<TrackerHit*> all_hits;
  all_hits.reserve(hits_in_fit.size() + outliers.size());

  for (unsigned ihit = 0; ihit < hits_in_fit.size(); ++ihit) {
    all_hits.push_back(hits_in_fit[ihit].first);
  }

  for (unsigned ihit = 0; ihit < outliers.size(); ++ihit) {
    all_hits.push_back(outliers[ihit].first);
  }

  UTIL::BitField64 encoder2(lcio::LCTrackerCellID::encoding_string());
  encoder2.reset();  // reset to 0
  MarlinTrk::addHitNumbersToTrack(mergedTrack.get(), all_hits, false, encoder2);
  MarlinTrk::addHitNumbersToTrack(mergedTrack.get(), hits_in_fit, true, encoder2);

  if (streamlog::out.write<streamlog::DEBUG5>()) {
    streamlog_out(DEBUG5) << " Merged track : " << std::endl;
  }

  lcioTrkPtr = mergedTrack.release();
}

void ClonesAndSplitTracksFinder::bestInClones(edm4hep::Track track_a, edm4hep::Track track_b, int nOverlappingHits, edm4hep::Track*& bestTrack) {
  // This function compares two tracks which have a certain number of overlapping hits and returns the best track
  // The best track is chosen based on length (in terms of number of hits) and chi2/ndf requirements
  // In general, the longest track is preferred. When clones have same length, the one with best chi2/ndf is chosen

  podio::RelationRange<edm4hep::TrackerHit> trackerHit_a = track_a.getTrackerHits();
  podio::RelationRange<edm4hep::TrackerHit> trackerHit_b = track_b.getTrackerHits();

  int trackerHit_a_size = trackerHit_a.size();
  int trackerHit_b_size = trackerHit_b.size();

  double b_chi2 = track_b.getChi2() / track_b.getNdf();
  double a_chi2 = track_a.getChi2() / track_a.getNdf();

  if (nOverlappingHits == trackerHit_a_size) {  // if the second track is the first track + segment
    bestTrack = &track_b;
  } else if (nOverlappingHits == trackerHit_b_size) {  // if the second track is a subtrack of the first track
    bestTrack = &track_a;
  } else if (trackerHit_b_size == trackerHit_a_size) {  // if the two tracks have the same length
    if (b_chi2 <= a_chi2) {
      bestTrack = &track_b;
    } else {
      bestTrack = &track_a;
    }
  } else if (trackerHit_b_size > trackerHit_a_size) {  // if the second track is longer
    bestTrack = &track_b;
  } else if (trackerHit_b_size < trackerHit_a_size) {  // if the second track is shorter
    bestTrack = &track_a;
  }
}
