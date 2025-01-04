#ifndef ClonesAndSplitTracksFinder_h
#define ClonesAndSplitTracksFinder_h 1

#include <k4FWCore/Transformer.h>

#include "BitField64.h"

#include <edm4hep/Track.h>
#include <edm4hep/MutableTrack.h>
#include <edm4hep/TrackCollection.h>

#include <cfloat>
#include <vector>

/*
namespace MarlinTrk {
  class IMarlinTrkSystem;
  class IMarlinTrack;
}
*/

class ClonesAndSplitTracksFinder : public k4FWCore::Transformer<edm4hep::TrackCollection(const edm4hep::TrackCollection&)> {
public:
  virtual k4FWCore::Transformer* newProcessor() { return new ClonesAndSplitTracksFinder; }

  ClonesAndSplitTracksFinder(const std::string& name, ISvcLocator* svcLoc);
  ClonesAndSplitTracksFinder(const ClonesAndSplitTracksFinder&) = delete;
  ClonesAndSplitTracksFinder& operator=(const ClonesAndSplitTracksFinder&) = delete;

  // Initialisation - run at the beginning to start histograms, etc.
  virtual StatusCode initialize();

  // Run over each event - the main algorithm
  virtual edm4hep::TrackCollection operator(const edm4hep::TrackCollection& input_track_col) const;

  // Called at the very end for cleanup, histogram saving, etc.
  virtual StatusCode finalize();

protected:
  // Checks for overlapping hits
  int overlappingHits(const edm4hep::Track*, const edm4hep::Track*);

  // Picks up the best track between two clones (based on chi2 and length requirements)
  void bestInClones(edm4hep::Track, edm4hep::Track, int, edm4hep::Track*&);

  // Service function to set the information from a Track* object to a MutableTrack* object
  void fromTrackToMutbleTrack(const edm4hep::Track*, edm4hep::MutableTrack*&);

  // Merges hits from two tracks in one and fits it
  void mergeAndFit(edm4hep::Track, edm4hep::Track, edm4hep::Track*&);

  // Removes doubles (from clone treatments and track merging) and filters multiple connections (clones and mergeable tracks treated differently)
  void filterClonesAndMergedTracks(std::multimap<int, std::pair<int, edm4hep::Track*>>&, edm4hep::TrackCollection&, vector<edm4hep::Track>&, bool);

  // Contains the whole merging procedure (calls filterClonesAndMergedTracks(bool false) and mergeAndFit)
  void mergeSplitTracks(edm4hep::TrackCollection&, const edm4hep::TrackCollection&, vector<edm4hep::Track>&);

  // Calculate significance in pt for two candidate clones
  double calculateSignificancePt(const edm4hep::Track, const edm4hep::Track);

  // Calculate significance in phi for two candidate clones
  double calculateSignificancePhi(const edm4hep::Track, const edm4hep::Track);

  // Calculate significance in tanLambda for two candidate clones
  double calculateSignificanceTanLambda(const edm4hep::Track, const edm4hep::Track);

  // Calculate significance for two candidate clones
  double calculateSignificance(const double firstPar, const double secondPar, const double firstPar_sigma,
                               const double secondPar_sigma);

  // Contains the whole clone skimming procedure (calls bestInClones and filterClonesAndMergedTracks(bool true))
  void removeClones(vector<edm4hep::Track>&, const edm4hep::TrackCollection&);

  //MarlinTrk::IMarlinTrkSystem* _trksystem = nullptr;

  Gaudi::Property<bool>   m_MSOn{this, "MultipleScatteringOn", true, "Use MultipleScattering in Fit"};
  Gaudi::Property<bool>   m_ElossOn{this, "EnergyLossOn", true, "Use Energy Loss in Fit"};
  Gaudi::Property<bool>   m_SmoothOn{this, "SmoothOn", false, "Smooth All Measurement Sites in Fit"};
  double m_magneticField = 0.0;
  Gaudi::Property<bool>   m_extrapolateForward{this, "extrapolateForward", true, "if true extrapolation in the forward direction in-out), otherwise backward (out-in)"};

  Gaudi::Property<double> m_minPt{this, "minTrackPt", 1.0, "minimum track pt for merging (in GeV/c)"};
  Gaudi::Property<double> m_maxSignificanceTheta{this, "maxSignificanceTheta", 0.0, "maximum significance separation in tanLambda"}; 
  Gaudi::Property<double> m_maxSignificancePhi{this, "maxSignificancePhi", 0.0, "maximum significance separation in phi"};
  Gaudi::Property<double> m_maxSignificancePt{this, "maxSignificancePt", 0.0, "maximum significance separation in pt"};

  Gaudi::Property<bool> m_mergeSplitTracks{this, "mergeSplitTracks", false, "if true, the merging of split tracks is performed"};

  // Track fit parameters
  double m_initialTrackError_d0;
  double m_initialTrackError_phi0;
  double m_initialTrackError_omega;
  double m_initialTrackError_z0;
  double m_initialTrackError_tanL;
  double m_maxChi2perHit;

  std::shared_ptr<BitField64> m_encoder{};
};

#endif
