/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */

#include "DDTPCDigiAlgorithm.h"

#include "FixedPadSizeDiskLayout.h"
#include "TPCModularEndplate.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

#include <gsl/gsl_randist.h>

#include "Circle.h"
#include "SimpleHelix.h"
#include "constants.h"
#include "LCCylinder.h"

//stl exception handler
#include <stdexcept>
#include "constants.h"
#include "voxel.h"

// EDM4HEP

// --- DD4hep ---
#include "DD4hep/Detector.h"
#include "DDRec/Vector3D.h"
#include "DD4hep/DD4hepUnits.h" 

using namespace constants;

DECLARE_COMPONENT(DDTPCDigiAlgorithm)

bool compare_phi( Voxel_tpc* a, Voxel_tpc* b)  { 
  return ( a->getPhiIndex() < b->getPhiIndex() ) ; 
} 

bool compare_z( Voxel_tpc* a, Voxel_tpc* b) { 
  return ( a->getZIndex() < b->getZIndex() ) ; 
} 


DDTPCDigiAlgorithm::~DDTPCDigiAlgorithm(){
  delete m_tpcEP ;
}

DDTPCDigiAlgorithm::DDTPCDigiAlgorithm(const std::string& name, ISvcLocator* svcLoc) : MultiTransformer(name, svcLoc,
      { KeyValues("TPCPadRowHitCollectionName", {"TPCCollection", "TPCSpacePointCollection", "TPCLowPtCollection"}))
        KeyValues("EventHeaderCollectionName", {"EventHeader"}) }, 
      { KeyValues("TPCTrackerHitsCol", {"TPCTrackerHits"}),
        KeyValues("SimTrkHitRelCollection", {"TPCTrackerHitRelations"}) } {}


StatusCode DDTPCDigiProcessor::initialize() 
{
  ITHistSvc* histSvc{nullptr};
  StatusCode sc1 = service("THistSvc", histSvc);
  if ( sc1.isFailure() ) {
    error() << "Could not locate HistSvc" << endmsg;
    return StatusCode::FAILURE;
  } 
  
  m_phiDiffHisto = new TH1D("Histograms/phi_diff",
                                         "Calculated Phi - Track Phi",
                                         201, -0.05, 0.05);
  
  m_thetaDiffHisto = new TH1D("Histograms/theta_diff",
                                           "Calculated Theta - Track Theta",
                                           201, -0.05, 0.05);
  
  m_phiRelHisto = new TH1D("Histograms/padPhi",
                                        "Phi Relative to the Pad",
                                        201, 0.0, 6.3);
  
  m_thetaRelHisto = new TH1D("Histograms/padtheta",
                                          "Theta Relative to the pad",
                                          201, 0.0, 6.3);
  
  m_rPhiDiffHisto = new TH1D("Histograms/rPhiDiff",
                                          "rPhi_rec - rPhi_sim",
                                          201, -1.0, 1.0);
  
  m_zDiffHisto = new TH1D("Histograms/zDiff",
                                       "Z_rec - Z_sim",
                                       201, -1.0, 1.0);
  
  m_zPullHisto = new TH1D("Histograms/zPull",
                                       "(z_rec - z_sim) / Sigma_z",
                                       201, -10.0, 10.0);
  
  m_phiDistHisto = new TH1D("Histograms/phiDist",
                                         "phi_rec - Phi_sim",
                                         201, -1.0, 1.0);
  
  m_rPhiPullHisto = new TH1D("Histograms/rPhiPull",
                                          "(rPhi_rec - rPhi_sim) / Sigma_rPhi",
                                          201, -10.0, 10.0);
  
  m_zSigmaVsZHisto = _HF->createHistogram2D("Histograms/zSigmaVsZ",
                                           "z Sigma vs Z ",
                                           3000, 0.0, 3000.0,
                                           201, -0.20, 5.20);
  
  m_zSigmaHisto = new TH1D("Histograms/zSigma",
                                        "z Sigma ",
                                        201, -0.20, 5.20);
  
  m_rPhiSigmaHisto = new TH1D("Histograms/rPhiSigma",
                                           "rPhi Sigma",
                                           201, -0.20, 0.20);
  
  m_radiusCheckHisto = new TH1D("Histograms/radiusCheck",
                                             "R_hit - TPC Rmin - ((RowIndex + 0.5 )* padheight)",
                                             201, -0.20, 0.20);
  
  m_ResidualsRPhiHisto = new TH1D("Histograms/ResidualsRPhi",
                                               "MC Track Phi - Hit Phi",
                                               50, -0.001, 0.001);
  
  m_NSimTPCHitsHisto = new TH1D("Histograms/SimTPCHits",
                                             "Number of SimTPC Hits",
                                             100, 0.0, 1000000.0);
  
  m_NBackgroundSimTPCHitsHisto = new TH1D("Histograms/NBackgroundSimTPCHits",
                                                       "Number of Background SimTPC Hits",
                                                       100, 0.0, 1000000.0);
  
  m_NPhysicsSimTPCHitsHisto = new TH1D("Histograms/NPhysicsSimTPCHits",
                                                    "Number of PhysicsSimTPC Hits",
                                                    100, 0.0, 100000.0);
  
  m_NPhysicsAbove02GeVSimTPCHitsHisto = new TH1D("Histograms/NPhysicsAbove02GeVTPCHits",
                                                              "Number of PhysicsSimTPC Hits above 0.2GeV pt",
                                                              100, 0.0, 100000.0);
  
  m_NPhysicsAbove1GeVSimTPCHitsHisto = new TH1D("Histograms/NPhysicsAbove1GeVPtTPCHits",
                                                             "Number of PhysicsSimTPC Hits above 1.0 GeV pt",
                                                             100, 0.0, 100000.0);
  
  m_NRecTPCHitsHisto = new TH1D("Histograms/NRecTPCHits",
                                             "Number of Rec TPC Hits",
                                             50, 0.0, 100000.0);
  
  m_NLostPhysicsTPCHitsHisto = new TH1D("Histograms/NLostPhysicsTPCHits",
                                                     "Number of PhysicsSimTPC Hits Lost",
                                                     100, 0.0, 5000.0);
  
  m_NLostPhysicsAbove02GeVPtTPCHitsHisto = new TH1D("Histograms/NLostPhysicsAbove02GeVPtTPCHits",
                                                                 "Number of PhysicsSimTPC Hits Lost above 0.2 GeV pt",
                                                                 100, 0.0, 5000.0);
  
  m_NLostPhysicsAbove1GeVPtTPCHitsHisto = new TH1D("Histograms/NLostPhysicsAbove1GeVPtTPCHits",
                                                                "Number of PhysicsSimTPC Hits Lost above 1.0 GeV pt",
                                                                100, 0.0, 1000.0);
  
  m_NRevomedHitsHisto = new TH1D("Histograms/NRevomedHits",
                                              "Number of Removed TPC hits",
                                              100, 0.0, 1000000.0);
  
  
  m_NKeptPhysicsTPCHitsHistoPercent = new TH1D("Histograms/NKeptPhysicsTPCHitsPercent",
                                                            "Number of PhysicsSimTPC Hits Kept",
                                                            303, 0.0, 1.01);
  
  m_NKeptPhysicsAbove02GeVPtTPCHitsHistoPercent = new TH1D("Histograms/NKeptPhysicsAbove02GeVPtTPCHitsPercent",
                                                                        "Number of PhysicsSimTPC Hits Kept above 0.2 GeV pt",
                                                                        303, 0.0, 1.01);
  
  m_NKeptPhysicsAbove1GeVPtTPCHitsHistoPercent = new TH1D("Histograms/NKeptPhysicsAbove1GeVPtTPCHitsPercent",
                                                                       "Number of PhysicsSimTPC Hits Kept above 1.0 GeV pt",
                                                                       303, 0.0, 1.01);
  
  
  //--- get the geometry data from dd4hep
  dd4hep::Detector& theDet = dd4hep::Detector::getInstance();
  dd4hep::DetElement tpcDE = theDet.detector("TPC") ;
  m_tpc = tpcDE.extension<dd4hep::rec::FixedPadSizeTPCData>() ;
  

  // fill the data for the TPC endplate
  m_tpcEP = new TPCModularEndplate( m_tpc ) ;

  if( m_tpcEndPlateModuleNumbers.size() != m_tpcEndPlateModulePhi0s.size() ){

    throw Exception(" DDTPCDigiProcessor: parameters tpcEndPlateModuleNumbers and tpcEndPlateModulePhi0s dont have the same number of elements ( module rings ) !! " ) ;
  }

  for(unsigned i=0, N=m_tpcEndPlateModuleNumbers.size() ; i<N ; ++i){

    m_tpcEP->addModuleRing(  m_tpcEndPlateModuleNumbers[i] , m_tpcEndPlateModulePhi0s[i] ) ;
  }
  m_tpcEP->initialize() ;

  // -----

  debug() << " initialized TPC geometry from TPCData: " <<  *m_tpc << endmsg; 

  double bfieldV[3] ;
  theDet.field().magneticField( { 0., 0., 0. }  , bfieldV  ) ;
  m_bField = bfieldV[2]/dd4hep::tesla ;
  //----
  
  //intialise random number generator 
  m_random = gsl_rng_alloc(gsl_rng_ranlxs2);
  m_idGen = serviceLocator()->service("UniqueIDGenSvc");
}

std::tuple<edm4hep::TrackerHitPlaneColection,
             edm4hep::TrackerHitSimTrackerHitLinkCollection> operator(
             const std::vector<const edm4hep::SimTrackerHitCollection*>& inputCols
             const edm4hep::EventHeaderCollection& evHeader) const{
  
  gsl_rng_set( m_random, m_idGen->getUniqueID(evHeader.eventNumber()[0], evHeader.runNumber()[0] );
  debug() << "seed set to " << m_idGen->getUniqueID(evHeader.eventNumber()[0], evHeader.runNumber()[0] << " for event number "<< evHeader.eventNumber()[0] << endmsg;

   int numberOfVoxelsCreated(0);
  
  m_NSimTPCHits = 0;
  m_NBackgroundSimTPCHits = 0;
  m_NPhysicsSimTPCHits = 0;
  m_NPhysicsAbove02GeVSimTPCHits = 0;
  m_NPhysicsAbove1GeVSimTPCHits = 0;
  m_NRecTPCHits = 0;
  
  m_NLostPhysicsTPCHits = 0;
  m_NLostPhysicsAbove02GeVPtTPCHits = 0;
  m_NLostPhysicsAbove1GeVPtTPCHits = 0;
  m_NRevomedHits = 0;
  
  static bool firstEvent = true;
  m_tpcHitMap.clear();
  m_tpcRowHits.clear();

  if(firstEvent==true) {
    if (! m_use_raw_hits_to_store_simhit_pointer ) {
      debug() << "SimTrackerHits are no longer stored in RawTrackerHits. Enable this deprecated feature by setting UseRawHitsToStoreSimhitPointer to true in steering file." << endmsg;
    }
    else{ 
        debug() << "SimTrackerHits will be stored in RawTrackerHits. This is a deprecated please use the relations." << endmsg;
    }
  }
  
  firstEvent = false ;
  

  m_padWidth = m_tpc->padWidth/dd4hep::mm; 
  // set size of row_hits to hold (n_rows) vectors
  m_tpcRowHits.resize( m_tpc->maxRow );
  
  // relations from created trackerhits to the SimTrackerHits that caused them
  edm4hep::TrackerHitSimTrackerHitLinkCollection hitSimHitCol;
  edm4hep::TrackerHitPlaneCollection hitCol;
  edm4hep::SimTrackerHitCollection STHcol = inputCols[0];

  float edep0=0.0;
    
    int n_sim_hits = STHcol.size();
    
    LCFlagImpl colFlag( STHcol->getFlag() ) ;    
    
    m_NSimTPCHits = n_sim_hits;
    
    debug() << "number of Pad-Row based SimHits = " << n_sim_hits << endmsg;
    
    
    // make sure that all the pointers are initialise to NULL
    m_mcp=NULL;         
    m_previousMCP=NULL; 
    m_nextMCP=NULL;     
    m_nMinus2MCP=NULL;  
    m_nPlus2MCP=NULL;   
    
    m_SimTHit=NULL;
    m_previousSimTHit=NULL;
    m_nextSimTHit=NULL;
    m_nMinus2SimHit=NULL;
    m_nPlus2SimHit=NULL;
    
    // loop over all the pad row based sim hits
    for(int i=0; i< n_sim_hits; i++){
      
      // this will used for nominaml smearing for very low pt rubish, so set it to zero initially
      double ptSqrdMC = 0;
      
      m_SimTHit = &(STHcol.at( i ));
      
      float edep;
      double padPhi(0.0);
      double padTheta (0.0);
      
      
      debug() << "processing hit " << i << "\n"
              << " address = " << m_SimTHit  
              << " x = "  << m_SimTHit->getPosition().x
              << " y = "  << m_SimTHit->getPosition().y
              << " z = "  << m_SimTHit->getPosition().z
              << endmsg; 
     

      CLHEP::Hep3Vector thisPoint(m_SimTHit->getPosition().x, m_SimTHit->getPosition().y, m_SimTHit->getPosition().z);
      double padheight = m_tpc->padHeight/dd4hep::mm ;
      
      // conversion constant. r = pt / (FCT*_bField)
      const double FCT = 2.99792458E-4;
      
      m_mcp = &(m_SimTHit.getParticle()); 
      
      // increase the counters for the different classification of simhits
      if(m_mcp){ 
        
        // get the pt of the MCParticle, this will used later to uses nominal smearing for low momentum rubish
        edm4hep::Vector3d momentumMC = m_mcp->getMomentum();
        ptSqrdMC = momentumMC.x*momentumMC.x+momentumMC.y*momentumMC.y ; 
        
        debug() << " mcp address = " << m_mcp
                << " px = "  << momentumMC.x
                << " py = "  << momentumMC.y
                << " pz = "  << momentumMC.z
                << endmsg;
        
        // SJA:FIXME: the fact that it is a physics hit relies on the fact that for overlay 
        // the pointer to the mcp is set to NULL. This distinction may not always be true ...
        ++m_NPhysicsSimTPCHits ;
        if( ptSqrdMC > (0.2*0.2) ) ++m_NPhysicsAbove02GeVSimTPCHits ;
        if( ptSqrdMC > 1.0 )  ++m_NPhysicsAbove1GeVSimTPCHits ;
        
        if(m_mcp) plotHelixHitResidual(m_mcp, &thisPoint);

      } else {
        ++m_NBackgroundSimTPCHits;
      }
      
      // if the hits contain the momentum of the particle use this to calculate the angles relative to the pad 
      if(colFlag.bitSet(LCIO::THBIT_MOMENTUM)) {
        
        edm4hep::Vector3d mcpMomentum = m_SimTHit->getMomentum();
        
        CLHEP::Hep3Vector mom(mcpMomentum.x, mcpMomentum.y, mcpMomentum.z);
        
        // const double pt = mom.perp();
        // const double radius = pt / (FCT*_bField);
        
        // const double tanLambda = mom.z()/pt;
        
        padPhi = fabs(thisPoint.deltaPhi(mom));
        padTheta = mom.theta();
        
      } 
      
      else { // LCIO::THBIT_MOMENTUM not set
        
        // as the momentum vector is not available from the hits use triplets of 
        // hits to fit a circle and calculate theta and phi relative to the pad        
        
        if (!(m_mcp) || (sqrt(ptSqrdMC) / (FCT*m_bField)) < ( padheight / (0.1 * twopi))) { 
          // if the hit has no record of it MCParticle then there is no way to know if this hit has consecutive hits from the same MCParticle
          // so just set nominal values theta=phi=90
          // here make a cut for particles which will suffer more than a 10 percent change in phi over the distance of the pad 
          // R > padheight/(0.1*2PI) 
          // in both cases set the angles to 90 degrees
          padTheta = twopi/4.0 ;
          padPhi = twopi/4.0 ;    
        }
        else{
          
          // if there is at least one more hit after this one, set the pointer to the MCParticle for the next hit
          if (i < (n_sim_hits-1) ) {
            m_nextSimTHit = &(STHcol.at( i+1 ));
            m_nextMCP     = &(m_nextSimTHit.getParticle());
          }
          else{ // set make sure that the pointers are set back to NULL so that the comparisons later hold
            m_nextSimTHit = NULL;
            m_nextMCP = NULL;
          }
          // if there is at least two more hits after this one, set the pointer to the MCParticle for the next but one hit
          if (i < (n_sim_hits-2) ) {
            m_nPlus2SimHit = &(STHcol.at( i+2 ));
            m_nPlus2MCP   = &(m_nPlus2SimHit.getParticle());            
          }
          else{ // set make sure that the pointers are set back to NULL so that the comparisons later hold
            m_nPlus2SimHit = NULL;
            m_nPlus2MCP = NULL;
          }
          
          if      ( m_mcp==m_previousMCP && m_mcp==m_nextMCP )    { // middle hit of 3 from the same MCParticle 
            
            CLHEP::Hep3Vector precedingPoint(m_previousSimTHit->getPosition().x, m_previousSimTHit->getPosition().y, m_previousSimTHit->getPosition().z);
            CLHEP::Hep3Vector followingPoint(m_nextSimTHit->getPosition().x, m_nextSimTHit->getPosition().y, m_nextSimTHit->getPosition().z);
            
            debug() << "address of m_previousSimTHit = " << m_previousSimTHit  
                    << " x = "  << m_previousSimTHit->getPosition().x 
                    << " y = "  << m_previousSimTHit->getPosition().y
                    << " z = "  << m_previousSimTHit->getPosition().z
                    << endmsg; 
            
            debug() << "address of m_nextSimTHit = " << m_nextSimTHit 
                    << " x = "  << m_nextSimTHit->getPosition().x
                    << " y = "  << m_nextSimTHit->getPosition().y
                    << " z = "  << m_nextSimTHit->getPosition().z
                    << endmsg;
            
            // get phi and theta using functions defined below
            padPhi = getPadPhi( &thisPoint, &precedingPoint, &thisPoint, &followingPoint);
            padTheta = getPadTheta(&precedingPoint, &thisPoint, &followingPoint);
            
          }
          else if ( m_mcp==m_nextMCP     && m_mcp==m_nPlus2MCP )  { // first  hit of 3 from the same MCParticle
            
            CLHEP::Hep3Vector followingPoint(m_nextSimTHit->getPosition().x, m_nextSimTHit->getPosition().y, m_nextSimTHit->getPosition().z);
            CLHEP::Hep3Vector nPlus2Point(m_nPlus2SimHit->getPosition().x, m_nPlus2SimHit->getPosition().y, m_nPlus2SimHit->getPosition().z);
            
            // get phi and theta using functions defined below
            padPhi = getPadPhi( &thisPoint, &thisPoint, &followingPoint, &nPlus2Point);
            padTheta = getPadTheta(&thisPoint, &followingPoint, &nPlus2Point);
            
          }
          else if ( m_mcp==m_previousMCP && m_mcp==m_nMinus2MCP ) { // last   hit of 3 from the same MCParticle
            
            CLHEP::Hep3Vector nMinus2Point(m_nMinus2SimHit->getPosition().x, m_nMinus2SimHit->getPosition().y, m_nMinus2SimHit->getPosition().z);
            CLHEP::Hep3Vector precedingPoint(m_previousSimTHit->getPosition().x, m_previousSimTHit->getPosition().y, m_previousSimTHit->getPosition().z);
            
            // get phi and theta using functions defined below
            padPhi = getPadPhi( &thisPoint, &nMinus2Point, &precedingPoint, &thisPoint);
            padTheta = getPadTheta(&nMinus2Point, &precedingPoint, &thisPoint);
            
          }
          else{ // the hit is isolated as either a single hit, or a pair of hits, from a single MCParticle  
            padTheta = twopi/4.0 ;
            padPhi = twopi/4.0 ;    
          }
        }
        
        if(colFlag.bitSet(LCIO::THBIT_MOMENTUM)) {
          
          edm4hep::Vector3d mcpMomentum = m_SimTHit->getMomentum() ;
          
          CLHEP::Hep3Vector mom(mcpMomentum.x, mcpMomentum.y, mcpMomentum.z);
          
          double trackPhi = mom.phi();
          
          if(trackPhi<0.0) trackPhi=trackPhi+twopi;
          if(trackPhi>twopi) trackPhi=trackPhi-twopi;
          if(trackPhi>twopi/2.0) trackPhi = trackPhi - twopi/2.0 ;
          
          double localPhi = thisPoint.phi() - padPhi;
          
          m_phiRelHisto->Fill(padPhi);
          m_phiDiffHisto->Fill((fabs(localPhi - trackPhi))/trackPhi);
          m_thetaRelHisto->Fill(padTheta);
          m_thetaDiffHisto->Fill( (sin(padTheta) - sin(mom.theta()))/sin(mom.theta()) );
          
          debug() << "track Phi = " << trackPhi * (360.0 / twopi) << "\n" 
                  << "localPhi = " << localPhi * (360.0 / twopi) << "\n" 
                  << "pad Phi = " << padPhi * (360.0 / twopi) << "\n" 
                  << "pad Phi from track mom = " << ( thisPoint.phi() - trackPhi ) * (360.0 / twopi) << "\n" 
                  << "padTheta = " << padTheta * (360.0 / twopi) << "\n" 
                  << "padTheta from track mom = " << mom.theta() * (360.0 / twopi) << endmsg; 
        }       
        
      }
      
      //      int pad = padLayout.getNearestPad(thisPoint.perp(),thisPoint.phi());
      int layerNumber = m_SimTHit->getCellID();
      
      if(m_rejectCellID && (layerNumber<1)) {
        continue;
      }
      
      edep = m_SimTHit->getEDep();
      
      // Calculate Point Resolutions according to Ron's Formula 
      
      // sigma_{RPhi}^2 = sigma_0^2 + Cd^2/N_{eff} * L_{drift}
      
      // sigma_0^2 = (50micron)^2 + (900micron*sin(phi))^2
      // Cd^2/N_{eff}} = 25^2/(22/sin(theta)*h/6mm)
      // Cd = 25 ( microns / cm^(1/2) )
      // (this is for B=4T, h is the pad height = pad-row pitch in mm,
      // theta is the polar angle)       
      
      // sigma_{z}^2 = (400microns)^2 + L_{drift}cm * (80micron/sqrt(cm))^2 
      
      double aReso =m_pointResoRPhi0*m_pointResoRPhi0 + (m_pointResoPadPhi*m_pointResoPadPhi * sin(padPhi)*sin(padPhi)) ;
      double driftLength = m_tpc->driftLength/dd4hep::mm - (fabs(thisPoint.z()));
      
      if (driftLength <0) { 
	debug() << "DDTPCDigiProcessor : Warning! driftLength < 0 " << driftLength << " --> wrong data in dd4hep::rec::FixedPadSizeTPCData ? " << "\n" 
                << "Setting driftLength to 0.1" << "\n"
                << "m_tpc->driftLength/dd4hep::mm = " << m_tpc->driftLength/dd4hep::mm << endmsg; 
        driftLength = 0.10;
      }
      
      padheight = m_tpc->padHeight/dd4hep::mm ; 
      
      //double bReso = ( (_diffRPhi * _diffRPhi) / _nEff ) * sin(padTheta) * ( 6.0 / (padheight) )  * ( 4.0 / _bField  ) ;
      // formula with new quadratic B-field correction term
      double bReso = ( (m_diffRPhi * m_diffRPhi) / m_nEff ) * sin(padTheta) * ( 6.0 / (padheight) ) * ( (4.0 * 4.0) / (m_bField * m_bField) ) ;

      double tpcRPhiRes = sqrt( aReso + bReso * (driftLength / 10.0) ); // driftLength in cm
      
      double tpcZRes  = sqrt(( m_pointResoZ0 * m_pointResoZ0 ) 
                             + 
                             ( m_diffZ * m_diffZ ) * (driftLength / 10.0) ); // driftLength in cm 
      
      FixedPadSizeDiskLayout padLayout( m_tpc );
      int padIndex = padLayout.getNearestPad(thisPoint.perp(),thisPoint.phi());
      
      double TPCPadPlaneRMin = m_tpc->rMinReadout/dd4hep::mm ;
      double TPCPadPlaneRMax = m_tpc->rMaxReadout/dd4hep::mm ;
      
      int iRowHit = padLayout.getRowNumber(padIndex);
      int iPhiHit = padLayout.getPadNumber(padIndex);

      int NBinsZ =  (int) ((2.0 * m_tpc->driftLength/dd4hep::mm) / m_binningZ);
      int iZHit = (int) ( (float) NBinsZ * ( m_tpc->driftLength/dd4hep::mm + thisPoint.z() ) / ( 2.0 * m_tpc->driftLength/dd4hep::mm ) ) ;
      
      if(iZHit<0) iZHit=0;
      if(iZHit>NBinsZ) iZHit=NBinsZ;
      
      // make sure that the hit lies at the middle of the pad ring
      thisPoint.setPerp(padLayout.getPadCenter(padIndex)[0]);
      
      if( (thisPoint.perp() < TPCPadPlaneRMin) || (thisPoint.perp() > TPCPadPlaneRMax) ) {
        debug() << "Hit R not in TPC " << "\n"
                << "R = " << thisPoint.perp() << "\n" 
                << "the tpc InnerRadius = " << TPCPadPlaneRMin << "\n"
                << "the tpc OuterRadius = " << TPCPadPlaneRMax << "\n"
                << "Hit Dropped " << endmsg;
        continue;
      }
      
      if( (fabs(thisPoint.z()) > _tpc->driftLength/dd4hep::mm) ) {
        debug() << "Hit Z not in TPC " << "\n"
                << "Z = " << thisPoint.z() << "\n" 
                << "the tpc Max Z = " << m_tpc->driftLength/dd4hep::mm << "\n"
                << "Hit Dropped " << endmsg;
        continue; 
      }
      
      //get energy deposit of this row
      edep=m_SimTHit->getEDep();

      // create a tpc voxel hit and store it for this row
      Voxel_tpc * atpcVoxel = new Voxel_tpc(iRowHit,iPhiHit,iZHit, thisPoint, edep, tpcRPhiRes, tpcZRes);
      
      m_tpcRowHits.at(iRowHit).push_back(atpcVoxel);
      ++numberOfVoxelsCreated;
      
      // store the simhit pointer for this tpcvoxel hit in the hit map
      m_tpcHitMap[atpcVoxel] = m_SimTHit; 
      
      // move the pointers on 
      m_nMinus2MCP = m_previousMCP;
      m_previousMCP = m_mcp ;
      m_nMinus2SimHit = m_previousSimTHit;
      m_previousSimTHit = m_SimTHit;
      
    }
  
  
  // now process the LowPt collection
  edm4hep::SimTrackerHitCollection STHcolLowPt = inputCols[1];
    
    int n_sim_hitsLowPt = STHcolLowPt.size()  ;
    
    m_NBackgroundSimTPCHits += n_sim_hitsLowPt;
    m_NSimTPCHits += n_sim_hitsLowPt;
    
    debug() << "number of LowPt hits:" << n_sim_hitsLowPt endmsg;
    
    // loop over the LowPt hit collection
    for(int i=0; i< n_sim_hitsLowPt; i++){  
      
      m_SimTHit = &( STHcolLowPt.a( i ) ) ;
      
      CLHEP::Hep3Vector thisPoint(m_SimTHit->getPosition().x, m_SimTHit->getPosition().y, m_SimTHit->getPosition().z);
      
      FixedPadSizeDiskLayout padLayout( m_tpc );
      const std::vector<double>& planeExt = padLayout.getPlaneExtent() ;
      double TPCPadPlaneRMin = planeExt[0] ;
      double TPCPadPlaneRMax = planeExt[1] ;
      
      int NBinsZ =  (int) ((2.0 * m_tpc->driftLength/dd4hep::mm) / m_binningZ);
      
      if( (thisPoint.perp() < TPCPadPlaneRMin) || (thisPoint.perp() > TPCPadPlaneRMax) ) {
        debug() << "Hit R not in TPC " << "\n"
                << "R = " << thisPoint.perp() << "\n" 
                << "the tpc InnerRadius = " << TPCPadPlaneRMin << "\n"
                << "the tpc OuterRadius = " << TPCPadPlaneRMax << "\n"
                << "Hit Dropped " << endmsg;
        continue;
      }
      
      if( (fabs(thisPoint.z()) > _tpc->driftLength/dd4hep::mm) ) {
        debug() << "Hit Z not in TPC " << "\n"
                << "Z = " << thisPoint.z() << "\n"
                << "the tpc Max Z = " << m_tpc->driftLength/dd4hep::mm << "\n"
                << "Hit Dropped " << endmsg;
        continue; 
      }
      
      int padIndex = padLayout.getNearestPad(thisPoint.perp(),thisPoint.phi());
      
      int iRowHit = padLayout.getRowNumber(padIndex);
      int iPhiHit = padLayout.getPadNumber(padIndex);
      int iZHit = (int) ( (float) NBinsZ * 
                         ( m_tpc->driftLength/dd4hep::mm + thisPoint.z() ) / ( 2.0 * m_tpc->driftLength/dd4hep::mm ) ) ;
      
      // shift the hit in r-phi to the nearest pad-row centre 
      thisPoint.setPerp(padLayout.getPadCenter(padIndex)[0]);
      
      // set the resolutions to the pads to digital like values
      double tpcRPhiRes = m_padWidth;
      double tpcZRes = m_binningZ;
      
      //get energy deposit of this hit
      edep0 =m_SimTHit->getEDep();

     // create a tpc voxel hit for this simhit and store it for this tpc pad row
      Voxel_tpc * atpcVoxel = new Voxel_tpc(iRowHit,iPhiHit,iZHit, thisPoint, edep0, tpcRPhiRes, tpcZRes);
      
      m_tpcRowHits.at(iRowHit).push_back(atpcVoxel);
      ++numberOfVoxelsCreated;      
      
      // store the simhit pointer for this voxel hit in a map
      m_tpcHitMap[atpcVoxel] = m_SimTHit; 
      
    }
  
  
  int number_of_adjacent_hits(0);
  
  debug() << "finished looping over simhits, number of voxels = " << numberOfVoxelsCreated << endmsg;
  
  int numberOfhitsTreated(0);
  
  vector <Voxel_tpc *> row_hits;
  
  // loop over the tpc rows containing hits and check for merged hits
  for (unsigned int i = 0; i<_tpcRowHits.size(); ++i){
    
    row_hits = m_tpcRowHits.at(i);
    std::sort(row_hits.begin(), row_hits.end(), compare_phi );
    
    // double loop over the hits in this row 
    for (unsigned int j = 0; j<row_hits.size(); ++j){
      
      ++numberOfhitsTreated;      
      
      for (unsigned int k = j+1; k<row_hits.size(); ++k){
        
        if(row_hits[k]->getPhiIndex() > (row_hits[j]->getPhiIndex())+2){ // SJA:FIXME: here we need an OR to catch the wrap around
          break; // only compare hits in adjacent phi bins
        }
        
        // look to see if the two hit occupy the same pad in phi or if not whether they are within the r-phi double hit resolution
        else if( row_hits[k]->getPhiIndex()==row_hits[j]->getPhiIndex() 
                || 
                ( (fabs(row_hits[k]->getHep3Vector().deltaPhi(row_hits[j]->getHep3Vector()))) * row_hits[j]->getR()) < m_doubleHitResRPhi ) {
          
          // if neighboring in phi then compare z
          map <Voxel_tpc*,SimTrackerHit*> ::iterator it;
          
          SimTrackerHit* Hit1 = NULL;
          SimTrackerHit* Hit2 = NULL;
          
          // search of the simhit pointers in the tpchit map
          it=m_tpcHitMap.find(row_hits[j]);
          if(it!= m_tpcHitMap.end()) {
            Hit1 = it->second ; // hit found 
          }
          
          it=_tpcHitMap.find(row_hits[k]);
          if(it!= m_tpcHitMap.end()) {
            Hit2 = it->second ; // hit found 
          }
          
          double pathlengthZ1(0.0);
          double pathlengthZ2(0.0);
          
          if( Hit1 && Hit2 ){ // if both sim hits were found
            
            // check if the track momentum has been stored for the hits
            bool momentum_set = true;
            
            if( STHcol != NULL ){
              LCFlagImpl colFlag( STHcol->getFlag() ) ;
              momentum_set = momentum_set && colFlag.bitSet(LCIO::THBIT_MOMENTUM) ;
            }            
            
            if( STHcolLowPt != NULL ){
              LCFlagImpl colFlag( STHcolLowPt->getFlag() ) ;
              momentum_set =  momentum_set && colFlag.bitSet(LCIO::THBIT_MOMENTUM) ;
            }            
            
            if( momentum_set ){
              
              edm4hep::Vector3d Momentum1 = Hit1->getMomentum() ;
              edm4hep::Vector3d Momentum2 = Hit2->getMomentum() ;
              
              CLHEP::Hep3Vector mom1(Momentum1.x, Momentum1.y, Momentum1.z);
              CLHEP::Hep3Vector mom2(Momentum2.x, Momentum2.y, Momentum2.z);
              
              pathlengthZ1 = fabs( Hit1->getPathLength() * mom1.cosTheta() );
              pathlengthZ2 = fabs( Hit2->getPathLength() * mom2.cosTheta() );
            } 
            else {
              pathlengthZ1 = m_doubleHitResZ ; // assume the worst i.e. that the track is moving in z 
              pathlengthZ2 = m_doubleHitResZ ; // assume the worst i.e. that the track is moving in z 
            }
            
            double dZ = fabs(row_hits[j]->getZ() - row_hits[k]->getZ());
            
            double spacial_coverage = 0.5*(pathlengthZ1 + pathlengthZ2) + m_binningZ; 
            
            if( (dZ - spacial_coverage) < m_doubleHitResZ ) {
              
              row_hits[j]->setAdjacent(row_hits[k]);
              row_hits[k]->setAdjacent(row_hits[j]);
              ++number_of_adjacent_hits;
              
            }
          } else {
            debug() << "Hit1=" << Hit1 << "Hit2=" << Hit2 << endmsg; 
          }
        }
      }
    }
    
    
    // now all hits have been checked for adjacent hits, go throught and write out the hits or merge
    
    for (unsigned int j = 0; j<row_hits.size(); ++j){
      
      Voxel_tpc* seed_hit = row_hits[j];
      
      if(seed_hit->IsMerged() || seed_hit->IsClusterHit()) { 
        continue;
      }
      
      if(seed_hit->getNumberOfAdjacent()==0){ // no adjacent hits so smear and write to hit collection
        edm4hep::MutableTrackerHitSimTrackerHitLink hitLink = hitSimHitCol.create();
        edm4hep::MutableTrackerHitPlane hit = hitCol.create();
        writeVoxelToHit(seed_hit, &hitLink, &hit);
      }
      
      else if(seed_hit->getNumberOfAdjacent() < (m_maxMerge)){ // potential 3-hit cluster, can use simple average merge. 
        
        vector <Voxel_tpc*>* hitsToMerge = new vector <Voxel_tpc*>;
        
        int clusterSize = seed_hit->clusterFind(hitsToMerge);
        
        if( clusterSize <= m_maxMerge ){ // merge cluster
          seed_hit->setIsMerged();
          std::vector<edm4hep::MutableTrackerHitSimTrackerHitLink*> hitLinks;
          edm4hep::MutableTrackerHitPlane* hit = hitCol.create();
          for (int i = 0; i < hitsToMerge->size(); ++i) {
            hitLinks.push_back(&(hitSimHitCol.create()));
          }
          writeMergedVoxelsToHit(hitsToMerge, hitLinks, &hit);
        }
        delete hitsToMerge;
      } 
    } 
  }
  
  int numberOfHits(0);
  // count up the number of hits merged or lost
  for (unsigned int i = 0; i<_tpcRowHits.size(); ++i){
    row_hits = m_tpcRowHits.at(i);
    for (unsigned int j = 0; j<row_hits.size(); ++j){
      numberOfHits++;
      Voxel_tpc* seed_hit = row_hits[j];
      if(seed_hit->IsMerged() || seed_hit->IsClusterHit() || seed_hit->getNumberOfAdjacent() > m_maxMerge ) { 
        ++m_NRevomedHits;
        m_mcp = &((m_tpcHitMap[ seed_hit ])->getParticle()); 
        if(m_mcp != NULL ) { 
          ++_NLostPhysicsTPCHits;        
          edm4hep::Vector3d mom= m_mcp->getMomentum() ;
          double ptSQRD = mom.x*mom.x+mom.y*mom.y ; 
          if( ptSQRD > (0.2*0.2) ) ++m_NLostPhysicsAbove02GeVPtTPCHits ;
          if( ptSQRD > 1.0 )  ++m_NLostPhysicsAbove1GeVPtTPCHits ;
        }
      }
    }
  }
  
  debug() << "the number of adjacent hits is " <<  number_of_adjacent_hits << "  m_doubleHitResZ " << m_doubleHitResZ << "\n"  
          << "number of rec_hits = "  << m_NRecTPCHits << "\n"
          << "finished row hits " << numberOfHits << " " << numberOfhitsTreated << endmsg;    
  
  // set the parameters to decode the type information in the collection
  // for the time being this has to be done manually
  // in the future we should provide a more convenient mechanism to 
  // decode this sort of meta information
  
  //StringVec typeNames ;
  //IntVec typeValues ;
  //typeNames.push_back( LCIO::TPCHIT ) ;
  //typeValues.push_back( 1 ) ;
  //_trkhitVec->parameters().setValues("TrackerHitTypeNames" , typeNames ) ;
  //_trkhitVec->parameters().setValues("TrackerHitTypeValues" , typeValues ) ;

  // delete voxels
  for (unsigned int i = 0; i<_tpcRowHits.size(); ++i){
    vector <Voxel_tpc *>* current_row = &m_tpcRowHits.at(i);  
    for (unsigned int j = 0; j<current_row->size(); ++j){
      delete current_row->at(j);
    }
  }
  
  m_NSimTPCHitsHisto->Fill(m_NSimTPCHits);
  m_NBackgroundSimTPCHitsHisto->Fill(m_NBackgroundSimTPCHits);
  m_NPhysicsSimTPCHitsHisto->Fill(m_NPhysicsSimTPCHits);
  m_NPhysicsAbove02GeVSimTPCHitsHisto->Fill(m_NPhysicsAbove02GeVSimTPCHits);
  m_NPhysicsAbove1GeVSimTPCHitsHisto->Fill(m_NPhysicsAbove1GeVSimTPCHits);
  m_NRecTPCHitsHisto->Fill(m_NRecTPCHits);
  
  m_NLostPhysicsTPCHitsHisto->Fill(m_NLostPhysicsTPCHits);
  m_NLostPhysicsAbove02GeVPtTPCHitsHisto->Fill(m_NLostPhysicsAbove02GeVPtTPCHits);
  m_NLostPhysicsAbove1GeVPtTPCHitsHisto->Fill(m_NLostPhysicsAbove1GeVPtTPCHits);
  m_NRevomedHitsHisto->Fill(m_NRevomedHits);
  
  m_NKeptPhysicsTPCHitsHistoPercent->Fill( (float)(m_NPhysicsSimTPCHits-m_NLostPhysicsTPCHits) / (float)m_NPhysicsSimTPCHits );
  m_NKeptPhysicsAbove02GeVPtTPCHitsHistoPercent->Fill( (float)(m_NPhysicsAbove02GeVSimTPCHits-m_NLostPhysicsAbove02GeVPtTPCHits) / (float)m_NPhysicsAbove02GeVSimTPCHits);
  m_NKeptPhysicsAbove1GeVPtTPCHitsHistoPercent->Fill( (float)(m_NPhysicsAbove1GeVSimTPCHits-m_NLostPhysicsAbove1GeVPtTPCHits) / (float)m_NPhysicsAbove1GeVSimTPCHits );

   debug() << "m_NSimTPCHits = " << m_NSimTPCHits << "\n"
           << "m_NBackgroundSimTPCHits = " << m_NBackgroundSimTPCHits << "\n"
           << "m_NPhysicsSimTPCHits = " << m_NPhysicsSimTPCHits << "\n"
           << "m_NPhysicsAbove02GeVSimTPCHits = " << m_NPhysicsAbove02GeVSimTPCHits << "\n"
           << "m_NPhysicsAbove1GeVSimTPCHits = " << m_NPhysicsAbove1GeVSimTPCHits << "\n"
           << "m_NRecTPCHits = " << m_NRecTPCHits<< "\n"
           << "m_NLostPhysicsTPCHits = " << m_NLostPhysicsTPCHits << "\n"
           << "m_NLostPhysicsAbove02GeVPtTPCHits = " << m_NLostPhysicsAbove02GeVPtTPCHits << "\n"
           << "m_NLostPhysicsAbove1GeVPtTPCHits = " << m_NLostPhysicsAbove1GeVPtTPCHits << "\n"
           << "m_NRevomedHits = " << m_NRevomedHits << endmsg;
 
  //Clear the maps and the end of the event.
  m_tpcHitMap.clear();
  m_tpcRowHits.clear();
  
  return std::make_tuple(std::move(hitCol), std::move(hitSimHitCol));
}


StatusCode DDTPCDigiAlgorithm::finalize()
{ 
  gsl_rng_free(m_random);
}

void DDTPCDigiAlgorithm::writeVoxelToHit( Voxel_tpc* aVoxel, edm4hep::MutableTrackerHitSimTrackerHitLink* hitLink, edm4hep::MutableTrackerHitPlane* trkHit ){
  Voxel_tpc* seed_hit  = aVoxel;
  
  //  if( seed_hit->getRowIndex() > 5 ) return ;
  
  //now the hit pos has to be smeared
  
  double tpcRPhiRes = seed_hit->getRPhiRes();
  double tpcZRes = seed_hit->getZRes();
 
  CLHEP::Hep3Vector point(seed_hit->getX(),seed_hit->getY(),seed_hit->getZ());
  

  //-------- fg: remove the hit if it lies within a module boundary in phi --------

  dd4hep::rec::Vector3D hitPos( seed_hit->getX(),seed_hit->getY(),seed_hit->getZ());

  double hitGapDist = m_tpcEP->computeDistanceRPhi( hitPos ) ;

  if( hitGapDist < m_tpcEndPlateModuleGapPhi / 2. ){

    debug() << " removing hit in endplate module gap : " << hitPos << " - distance : " <<  hitGapDist << endmsg;

    return ;
  }

  //fg: here we could add additional effects, e.g larger smearing of hits in the
  //    vicinity of the gaps ( need extra parameters )
  //    or larger smearing in the pad row next to a module ring boundary in r 
  //    ...
  //------------------------------------------------------------------------------

  double unsmearedPhi = point.phi();
  
  double randrp = gsl_ran_gaussian(m_random,tpcRPhiRes);
  double randz =  gsl_ran_gaussian(m_random,tpcZRes);
  
  point.setPhi( point.phi() + randrp/ point.perp() );
  point.setZ( point.z() + randz );

  // make sure the hit is not smeared beyond the TPC Max DriftLength
  if( fabs(point.z()) > m_tpc->driftLength/dd4hep::mm ) point.setZ( (fabs(point.z()) / point.z() ) * m_tpc->driftLength/dd4hep::mm );

#if 0   // fg: it turns out that is more correct to allow reconstructed hits to be outside of the drift volume due to the large uncertainties
  // make sure the hit is not smeared beyond the cathode:
  double dzCathode = _tpc->zMinReadout/dd4hep::mm ;
  if( point.z() > 0. && hitPos.z() < 0. ){
    point.setZ( -dzCathode ) ;
  }
  if( point.z() < 0. && hitPos.z() > 0. ) {
    point.setZ( dzCathode ) ;
  }
#endif

  edm4hep::Vector3d pos{point.x(),point.y(),point.z()}; 
  trkHit->setPosition(pos);
  trkHit->setEDep(seed_hit->getEDep());
  //  trkHit->setType( 500 );
  
//  int side = lcio::ILDDetID::barrel ;
//  
//  if( pos[2] < 0.0 ) side = 1 ;
 
  ACTSTracking::BitField64 bitField( "system:5,side:-2,layer:6,module:11,sensor:8" );
  bitField.setFieldValue("system", 4) ; //TODO: lcio::ILDDetID::TPC = 4...
  bitField.setFieldValue("layer", seed_hit->getRowIndex());
  bitField.setFieldValue("module", 0);
  bitField.setFieldValue("side", (pos[2] < 0 ? -1 :1 ));

  trkHit->setCellID(bitField.getValue());
  
  
  // check values for inf and nan
  if( std::isnan(unsmearedPhi) || std::isinf(unsmearedPhi) || std::isnan(tpcRPhiRes) || std::isinf(tpcRPhiRes) ) {
    std::stringstream errorMsg;
    errorMsg << "\nProcessor: DDTPCDigiProcessor \n" 
             << "unsmearedPhi = "
             <<  unsmearedPhi
             << " tpcRPhiRes = "
             <<  tpcRPhiRes 
             << "\n" ;
    throw Exception(errorMsg.str());
  }
  
  // For no error in R
  edm4hep::CovMatrix3f covMat({float(sin(unsmearedPhi)*sin(unsmearedPhi)*tpcRPhiRes*tpcRPhiRes),
                               float(-cos(unsmearedPhi)*sin(unsmearedPhi)*tpcRPhiRes*tpcRPhiRes),
                               float(cos(unsmearedPhi)*cos(unsmearedPhi)*tpcRPhiRes*tpcRPhiRes),
                               float(0.),
                               float(0.),
                               float(tpcZRes*tpcZRes) });
  
  trkHit->setCovMatrix(covMat);      
  
  if( m_tpcHitMap[seed_hit] == NULL ){
    std::stringstream errorMsg;
    errorMsg << "\nProcessor: DDTPCDigiProcessor \n" 
             << "SimTracker Pointer is NULL throwing exception\n"
             << "\n" ;
    throw Exception(errorMsg.str());
  }
  
  if(pos[0]*pos[0]+pos[1]*pos[1]>0.0){ 
    //    push back the SimTHit for this TrackerHit

    if (m_use_raw_hits_to_store_simhit_pointer) {
      trkHit->rawHits().push_back( m_tpcHitMap[seed_hit] );
    }                        
    hitLink->setTo(m_tpcHitMap[seed_hit]);
    hitLink->setFrom(trkHit);
    hitLink->setWeight(1.0);

    m_NRecTPCHits++;
  }
  
  SimTrackerHit* theSimHit = m_tpcHitMap[seed_hit];
  double rSimSqrd = theSimHit->getPosition().x*theSimHit->getPosition().x + theSimHit->getPosition().y*theSimHit->getPosition().y;
  
  double phiSim = atan2(theSimHit->getPosition().y,theSimHit->getPosition().x);
  
  double rPhiDiff = (point.phi() - phiSim)*sqrt(rSimSqrd);
  double rPhiPull = ((point.phi() - phiSim)*sqrt(rSimSqrd))/(sqrt((covMat[2])/(cos(point.phi())*cos(point.phi()))));
  
  double zDiff = point.getZ() - theSimHit->getPosition().z;
  double zPull = zDiff/sqrt(covMat[5]);
  
  
  m_rPhiDiffHisto->Fill(rPhiDiff);
  m_rPhiPullHisto->Fill(rPhiPull);
  m_phiDistHisto->Fill(point.phi() - phiSim);
  m_zDiffHisto->Fill(zDiff);
  m_zPullHisto->Fill(zPull);
  
  m_zSigmaVsZHisto->Fill(seed_hit->getZ(),sqrt(covMat[5]));
  m_rPhiSigmaHisto->Fill(sqrt((covMat[2])/(cos(point.phi())*cos(point.phi()))));
  m_zSigmaHisto->Fill(sqrt(covMat[5]));
}

void DDTPCDigiAlgorithm::writeMergedVoxelsToHit( vector <Voxel_tpc*>* hitsToMerge, std::vector<edm4hep::MutableTrackerHitSimTrackerHitLink*> hitLinks, edm4hep::MutableTrackerHitPlane* trkHit){
  
  double sumZ = 0;
  double sumPhi = 0;
  double sumEDep = 0;
  //  double R = 0;
  double lastR = 0;
  
  unsigned number_of_hits_to_merge = hitsToMerge->size();
  

  for(unsigned int ihitCluster = 0; ihitCluster < number_of_hits_to_merge; ++ihitCluster){
    
    sumZ += hitsToMerge->at(ihitCluster)->getZ();
    sumPhi += hitsToMerge->at(ihitCluster)->getPhi();
    sumEDep += hitsToMerge->at(ihitCluster)->getEDep();
    hitsToMerge->at(ihitCluster)->setIsMerged();
    lastR = hitsToMerge->at(ihitCluster)->getR();

    if (m_use_raw_hits_to_store_simhit_pointer) {
      trkHit->rawHits().push_back( m_tpcHitMap[hitsToMerge->at(ihitCluster)] );
    }                        
    
    hitLink[ihitCluster]->setTo(m_tpcHitMap[hitsToMerge->at(ihitCluster)]);
    hitLink[ihitCluster]->setFrom(trkHit);
    hitLink[ihitCluster]->setWeight(float(1.0/number_of_hits_to_merge));
  }
  
  double avgZ = sumZ/(hitsToMerge->size());
  double avgPhi = sumPhi/(hitsToMerge->size());

  //set deposit energy as mean of merged hits
  sumEDep=sumEDep/(double)number_of_hits_to_merge;
  
  CLHEP::Hep3Vector* mergedPoint = new CLHEP::Hep3Vector(1.0,1.0,1.0);
  mergedPoint->setPerp(lastR);
  mergedPoint->setPhi(avgPhi);
  mergedPoint->setZ(avgZ);
  
  //store hit variables

  // first the hit pos has to be smeared------------------------------------------------
  
  //FIXME: which errors should we use for smearing the merged hits ?
  //       this might be a bit large ....
  double tpcRPhiRes = m_padWidth;
  double tpcZRes = m_binningZ;
 
  CLHEP::Hep3Vector point( mergedPoint->getX(), mergedPoint->getY(), mergedPoint->getZ()  ) ;
  
//  double unsmearedPhi = point.phi();
  
  double randrp = gsl_ran_gaussian(m_random,tpcRPhiRes);
  double randz =  gsl_ran_gaussian(m_random,tpcZRes);
  
  point.setPhi( point.phi() + randrp/ point.perp() );
  point.setZ( point.z() + randz );
  
  // make sure the hit is not smeared beyond the TPC Max DriftLength
  if( fabs(point.z()) > m_tpc->driftLength/dd4hep::mm ) point.setZ( (fabs(point.z()) / point.z() ) * m_tpc->driftLength/dd4hep::mm );

#if 0   // fg: it turns out that is more correct to allow reconstructed hits to be outside of the drift volume due to the large uncertainties
  // make sure the hit is not smeared onto the other side of the cathode:
  double dzCathode = _tpc->zMinReadout/dd4hep::mm ;
  if( point.z() > 0. && mergedPoint->getZ() < 0. ){
    point.setZ( -dzCathode ) ;
  }
  if( point.z() < 0. && mergedPoint->getZ() > 0. ){
    point.setZ( dzCathode ) ;
  }
#endif

  edm4hep::Vector3d pos{point.x(),point.y(),point.z()}; 

  //---------------------------------------------------------------------------------
  trkHit->setPosition(pos);
  trkHit->setEDep(sumEDep);
  //  trkHit->setType( 500 );
  
  FixedPadSizeDiskLayout padLayout( m_tpc );
  int padIndex = padLayout.getNearestPad(mergedPoint->perp(),mergedPoint->phi());  
  int row = padLayout.getRowNumber(padIndex);  
 
  ACTSTracking::BitField64 bitField( "system:5,side:-2,layer:6,module:11,sensor:8" );
  bitField.setFieldValue("system", 4) ; //TODO: lcio::ILDDetID::TPC = 4...
  bitField.setFieldValue("layer", row);
  bitField.setFieldValue("module", 0);
  bitField.setFieldValue("side", (pos[2] < 0 ? -1 :1 ));

  trkHit->setCellID(bitField.getValue());  
  
  double phi = mergedPoint->getPhi();
  
  // check values for inf and nan
  if( std::isnan(phi) || std::isinf(phi) || std::isnan(tpcRPhiRes) || std::isinf(tpcRPhiRes) ) {
    std::stringstream errorMsg;
    errorMsg << "\nProcessor: DDTPCDigiProcessor \n" 
             << "phi = "
             <<  phi 
             << " tpcRPhiRes = "
             <<  tpcRPhiRes 
             << "\n" ;
    throw Exception(errorMsg.str());
  }
  
  // For no error in R
  edm4hep::CovMatrix3f covMat({float(sin(phi)*sin(phi)*tpcRPhiRes*tpcRPhiRes),
                               float(-cos(phi)*sin(phi)*tpcRPhiRes*tpcRPhiRes),
                               float(cos(phi)*cos(phi)*tpcRPhiRes*tpcRPhiRes),
                               float(0.),
                               float(0.),
                               float(tpcZRes*tpcZRes)});
  
  trkHit->setCovMatrix(covMat);      
  
  ++m_nRechits;
  
  delete mergedPoint; 
}


void DDTPCDigiAlgorithm::plotHelixHitResidual( edm4hep::MCParticle *mcp, CLHEP::Hep3Vector *thisPoint){
 /* 
  const double FCT = 2.99792458E-4;
  double charge = mcp->getCharge();
  edm4hep::Vector3d *mom = mcp->getMomentum();
  double pt = sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
  double radius = pt / (FCT*_bField);
  double tanLambda = mom[2]/pt;
  double omega = charge / radius;
  
  if(pt>1.0) {
    
    //FIXME SJA: this is only valid for tracks from the IP and should be done correctly for non prompt tracks
    double Z0 = 0.;      
    double D0  = 0.;
    
    edm4hep::Vector3d refPoint(0.,0.,0.);
    
    SimpleHelix* helix = new SimpleHelix(D0, 
                                         atan2(mom[1],mom[0]), 
                                         omega, 
                                         Z0, 
                                         tanLambda,
                                         refPoint);
    
    
    // an almost "infinite" cylinder in z        
    edm4hep::Vector3d startCylinder(0.,0.,-1000000.0);
    edm4hep::Vector3d endCylinder(0.,0.,1000000.0);
    bool endplane=true;
    
    LCCylinder cylinder(startCylinder,endCylinder,thisPoint->perp(),endplane);
    
    bool pointExists = false;
    
    double pathlength = helix->getIntersectionWithCylinder( cylinder, pointExists);
    
    LCErrorMatrix* errors = new LCErrorMatrix();
    
    if(pointExists){
      
      edm4hep::Vector3d intersection = helix->getPosition(pathlength, errors); 
      
      double intersectionPhi = atan2(intersection[1],intersection[0]);
      double residualRPhi = ((intersectionPhi-thisPoint->phi())) ;
      m_ResidualsRPhiHisto->Fill(residualRPhi);
      
    }
    
    delete errors;
    delete helix;
    
    
    int row = padLayout.getRowNumber(padLayout.getNearestPad(thisPoint->perp(),thisPoint->phi()));
    int pad = padLayout.getNearestPad(thisPoint->perp(),thisPoint->phi());
    
    double rHit_diff = thisPoint->perp()  
    - padLayout.getPlaneExtent()[0]  
    - (( row + 0.5 ) 
       * m_tpc->padHeight/dd4hep::mm ) ;
    
    m_radiusCheckHisto->Fill(rHit_diff);
    
    //      streamlog_out(MESSAGE) << "$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$" << endl;
    //      streamlog_out(MESSAGE) << "thisPoint->perp() = " << thisPoint->perp() << endl;
    //      streamlog_out(MESSAGE) << "TPC Sensitive rMin = " << padLayout.getPlaneExtent()[0] << endl;
    //      streamlog_out(MESSAGE) << "Row number + 0.5 = " <<  row + 0.5 << endl;
    //      streamlog_out(MESSAGE) << "Pad Height = " <<  padLayout.getPadHeight(pad) << endl;
    //      streamlog_out(MESSAGE) << "Row Height = " <<   padLayout.getRowHeight(row) << endl;
    //      streamlog_out(MESSAGE) << "R_hit - TPC Rmin - ((RowIndex + 0.5 )* padheight) = " << rHit_diff << endl;
    
  }
  return;*/
}



double DDTPCDigiAlgorithm::getPadPhi( CLHEP::Hep3Vector *thisPoint, CLHEP::Hep3Vector* firstPoint, CLHEP::Hep3Vector* middlePoint, CLHEP::Hep3Vector* lastPoint){
  
  CLHEP::Hep2Vector firstPointRPhi(firstPoint->x(),firstPoint->y());
  CLHEP::Hep2Vector middlePointRPhi(middlePoint->x(),middlePoint->y());
  CLHEP::Hep2Vector lastPointRPhi(lastPoint->x(),lastPoint->y());
  
  // check that the points are not the same, at least at the level of a tenth of a micron 
  if( (fabs( firstPointRPhi.x() - middlePointRPhi.x() ) < 1.e-05  && fabs( firstPointRPhi.y() - middlePointRPhi.y() ) < 1.e-05) 
     ||
     (fabs( middlePointRPhi.x() - lastPointRPhi.x() ) < 1.e-05  && fabs( middlePointRPhi.y() - lastPointRPhi.y() ) < 1.e-05) 
     ||
     (fabs( firstPointRPhi.x() - lastPointRPhi.x() ) < 1.e-05  && fabs( firstPointRPhi.y() - lastPointRPhi.y() ) < 1.e-05) 
     ) {
    
    warning() << " DDTPCDigiProcessor::getPadPhi "  
              << "2 of the 3 SimTracker hits passed to Circle Fit are the same hit taking pad phi as PI/2\n"
              << " firstPoint->x() "  << firstPoint->x() 
              << " firstPoint->y() "  << firstPoint->y() 
              << " firstPoint->z() "  << firstPoint->z() 
              << " middlePoint->x() "  << middlePoint->x() 
              << " middlePoint->y() "  << middlePoint->y() 
              << " middlePoint->z() "  << middlePoint->z() 
              << " lastPoint->x() "  << lastPoint->x() 
              << " lastPoint->y() "  << lastPoint->y() 
              << " lastPoint.z() "  << lastPoint->z() 
              << endmsg;
    
    return twopi/4.0 ;
    
  }
  
  
  
  Circle theCircle(&firstPointRPhi, &middlePointRPhi, &lastPointRPhi);
  
  double localPhi = atan2((thisPoint->y() - theCircle.GetCenter()->y()) ,(thisPoint->x() - theCircle.GetCenter()->x())) + (twopi/4.0) ;
  
  if(localPhi>twopi) localPhi=localPhi - twopi;
  if(localPhi<0.0) localPhi=localPhi + twopi;
  if(localPhi>twopi/2.0) localPhi = localPhi - twopi/2.0 ;
  
  double pointPhi = thisPoint->phi();
  
  if(pointPhi>twopi) pointPhi=pointPhi - twopi;
  if(pointPhi<0.0) pointPhi=pointPhi + twopi;
  if(pointPhi>twopi/2.0) pointPhi = pointPhi - twopi/2.0 ;
  
  double padPhi = fabs(pointPhi - localPhi);
  
  // check that the value returned is reasonable 
  if( std::isnan(padPhi) || std::isinf(padPhi) ) {
    std::stringstream errorMsg;
    errorMsg << "\nProcessor: DDTPCDigiProcessor \n" 
             << "padPhi = "
             <<  padPhi
             << "\n" ;
    throw Exception(errorMsg.str());
  }
  
  return padPhi;
  
}

double DDTPCDigiAlgorithm::getPadTheta(CLHEP::Hep3Vector* firstPoint, CLHEP::Hep3Vector* middlePoint, CLHEP::Hep3Vector* lastPoint){
  
  // Calculate thetaPad for current hit
  CLHEP::Hep2Vector firstPointRPhi(firstPoint->x(),firstPoint->y()) ;
  CLHEP::Hep2Vector middlePointRPhi(middlePoint->x(),middlePoint->y());
  CLHEP::Hep2Vector lastPointRPhi(lastPoint->x(),lastPoint->y());
  
  // check that the points are not the same, at least at the level of a tenth of a micron 
  if( (fabs( firstPointRPhi.x() - middlePointRPhi.x() ) < 1.e-05  && fabs( firstPointRPhi.y() - middlePointRPhi.y() ) < 1.e-05) 
     ||
     (fabs( middlePointRPhi.x() - lastPointRPhi.x() ) < 1.e-05  && fabs( middlePointRPhi.y() - lastPointRPhi.y() ) < 1.e-05) 
     ||
     (fabs( firstPointRPhi.x() - lastPointRPhi.x() ) < 1.e-05  && fabs( firstPointRPhi.y() - lastPointRPhi.y() ) < 1.e-05) 
     ) {
    
    warning() << " DDTPCDigiProcessor::getPadTheta "  
              << "2 of the 3 SimTracker hits passed to Circle Fit are the same hit taking pad phi as PI/2\n"
              << " firstPoint->x() "  << firstPoint->x() 
              << " firstPoint->y() "  << firstPoint->y() 
              << " firstPoint->z() "  << firstPoint->z() 
              << " middlePoint->x() "  << middlePoint->x() 
              << " middlePoint->y() "  << middlePoint->y() 
              << " middlePoint->z() "  << middlePoint->z() 
              << " lastPoint->x() "  << lastPoint->x() 
              << " lastPoint->y() "  << lastPoint->y() 
              << " lastPoint.z() "  << lastPoint->z() 
              << endmsg;
    
    return twopi/4.0 ;
    
  }
  
  
  Circle theCircle(&firstPointRPhi, &middlePointRPhi, &lastPointRPhi);
  
  double deltaPhi = firstPoint->deltaPhi(*lastPoint);
  
  double pathlength = fabs(deltaPhi) *  theCircle.GetRadius();
  
  double padTheta = atan ( pathlength / fabs(lastPoint->z() - firstPoint->z()) ) ;
  
  double pathlength1 = 2.0 * asin( ( sqrt (
                                           (middlePointRPhi.x() - firstPointRPhi.x()) * (middlePointRPhi.x()-firstPointRPhi.x())
                                           +
                                           (middlePointRPhi.y()-firstPointRPhi.y()) * (middlePointRPhi.y()-firstPointRPhi.y())
                                           ) / 2.0 ) / theCircle.GetRadius() ) * theCircle.GetRadius()  ;
  
  
  double pathlength2 = ( ( sqrt (  (lastPointRPhi.x()-middlePointRPhi.x()) * (lastPointRPhi.x()-middlePointRPhi.x())
                                   + (lastPointRPhi.y()-middlePointRPhi.y()) * (lastPointRPhi.y()-middlePointRPhi.y())
                                   ) / 2.0 ) / theCircle.GetRadius() >= 1.0 )  ?
    ( 2.0 * asin( 1.0 ) * theCircle.GetRadius() )
    :
    ( 2.0 * asin( ( sqrt (
                          (lastPointRPhi.x()-middlePointRPhi.x()) * (lastPointRPhi.x()-middlePointRPhi.x())
                          +
                          (lastPointRPhi.y()-middlePointRPhi.y()) * (lastPointRPhi.y()-middlePointRPhi.y())
                          ) / 2.0 ) / theCircle.GetRadius() ) * theCircle.GetRadius() ) ;
  
  padTheta = atan ((fabs(pathlength1 + pathlength2)) / (fabs(lastPoint->z() - firstPoint->z())) ) ;
  
  // check that the value returned is reasonable 
  if( std::isnan(padTheta) || std::isinf(padTheta) ) {
    std::stringstream errorMsg;
    errorMsg << "\nProcessor: DDTPCDigiProcessor \n" 
             << "padTheta = "
             <<  padTheta 
             << "\n" ;
    throw Exception(errorMsg.str());
  }
  
  return padTheta;
  
}
