/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "DDPlanarDigiAlgorithm.h"

// EDM4HEP
#include <edm4hep/Vector2f.h>
#include <edm4hep/MCParticle.h>
#include <edm4hep/SimTrackerHit.h>
#include <edm4hep/MutableTrackerHitPlane.h>
#include <edm4hep/MutableTrackerHitSimTrackerHitLink.h>

// Gaudi
#include <GaudiKernel/ITHistSvc.h>
#include "GaudiKernel/IRndmEngine.h"

// DD4HEP
#include "DD4hep/Detector.h"
#include "DDSegmentation/BitFieldCoder.h"
#include "DD4hep/DD4hepUnits.h"

// Standard and ROOT
#include <TMath.h>
#include <cmath>
#include <cfloat>
#include <sstream>
#include <climits>
#include <iostream>
#include <algorithm>

DECLARE_COMPONENT(DDPlanarDigiAlgorithm)

DDPlanarDigiAlgorithm::DDPlanarDigiAlgorithm(const std::string& name, ISvcLocator* svcLoc) : MultiTransformer(name, svcLoc,
       { KeyValues("SimTrackHitCollectionName", {"VXDCollection"}) },
       { KeyValues("TrackerHitCollectionName", {"VTXTrackerHits"}),
         KeyValues("SimTrkHitRelCollection", {"VTXTrackerHitRelations"}) }) 
{
  m_geoSvc = serviceLocator()->service("GeoSvc");  // important to initialize m_geoSvc
}

enum {
  hu = 0,
  hv,
  hT,
  hitE,
  hitsAccepted,
  diffu,
  diffv,
  diffT,
  hSize 
} ;

StatusCode DDPlanarDigiAlgorithm::initialize() {
  MsgStream log(msgSvc(), name());
  // Initialize random number gen
  if ( (m_rng.initialize(randSvc(), Rndm::Gauss(0., 1.))).isFailure() ) {
    log << MSG::ERROR << "Unable to initialize rng Gauss." << endmsg;
    return StatusCode::FAILURE;
  }

  m_h.resize( hSize );
  
  if( m_resU.size() !=  m_resV.size() ) {
    log << MSG::ERROR << "Inconsistent number of resolutions given for U and V coordinate: " 
        << "ResolutionU  :" <<   m_resU.size() << " != ResolutionV : " <<  m_resV.size() 
        << endmsg;
    return StatusCode::FAILURE;
  }

  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();


  //===========  get the surface map from the SurfaceManager ================

  dd4hep::rec::SurfaceManager& surfMan = *theDetector.extension<dd4hep::rec::SurfaceManager>() ;

  dd4hep::DetElement det = theDetector.detector( m_subDetName ) ;

  m_map = surfMan.map( det.name() ) ;

  if( ! m_map ) {   
    log << MSG::ERROR<< "Could not find surface map for detector: " 
        << m_subDetName << " in SurfaceManager " << endmsg;
    return StatusCode::FAILURE;
  }

  log << MSG::DEBUG << "DDPlanarDigiProcessor::init(): found " << m_map->size() 
      << " surfaces for detector:" <<  m_subDetName << endmsg;

  log << MSG::INFO << " *** DDPlanarDigiProcessor::init(): creating histograms" << endmsg;

  SmartIF<ITHistSvc> histSvc;
  histSvc = serviceLocator()->service("THistSvc");

  m_h[ hu ] = new TH1F( "hu" , "smearing u" , 50, -5. , +5. );
  m_h[ hv ] = new TH1F( "hv" , "smearing v" , 50, -5. , +5. );
  m_h[ hT ] = new TH1F( "hT" , "smearing time" , 50, -5. , +5. );

  (void)histSvc->regHist("/histos/digi_planar/hu", m_h[ hu ]);
  (void)histSvc->regHist("/histos/digi_planar/hv", m_h[ hv ]);
  (void)histSvc->regHist("/histos/digi_planar/hT", m_h[ hT ]);

  m_h[ diffu ] = new TH1F( "diffu" , "diff u" , 1000, -5. , +5. );
  m_h[ diffv ] = new TH1F( "diffv" , "diff v" , 1000, -5. , +5. );
  m_h[ diffT ] = new TH1F( "diffT" , "diff time" , 1000, -5. , +5. );

  (void)histSvc->regHist("/histos/digi_planar/diffu", m_h[ diffu ]);
  (void)histSvc->regHist("/histos/digi_planar/diffv", m_h[ diffv ]);
  (void)histSvc->regHist("/histos/digi_planar/diffT", m_h[ diffT ]);

  m_h[ hitE ] = new TH1F( "hitE" , "hitEnergy in keV" , 1000, 0 , 200 );
  m_h[ hitsAccepted ] = new TH1F( "hitsAccepted" , "Fraction of accepted hits [%]" , 201, 0 , 100.5 );
  
  (void)histSvc->regHist("/histos/digi_planar/hitE", m_h[ hitE ]);
  (void)histSvc->regHist("/histos/digi_planar/hitsAccepted", m_h[ hitsAccepted ]);
  
  return StatusCode::SUCCESS;
}

std::tuple<edm4hep::TrackerHitPlaneCollection, 
           edm4hep::TrackerHitSimTrackerHitLinkCollection> DDPlanarDigiAlgorithm::operator()(
     const edm4hep::SimTrackerHitCollection& inputSim) const{
    MsgStream log(msgSvc(), name());
   
    std::vector<long> seeds;
    (void)randSvc().get()->engine()->seeds(seeds);
    log << MSG::DEBUG << "seed set to ";
    for (long seed : seeds) {
      log << MSG::DEBUG << seed;
    }
    log << MSG::DEBUG << endmsg;

    unsigned nCreatedHits=0;
    unsigned nDismissedHits=0;
    
    edm4hep::TrackerHitPlaneCollection trkhitCol;

    // Relation collection TrackerHit, SimTrackerHit
    edm4hep::TrackerHitSimTrackerHitLinkCollection relCollection;
    
    std::string initString;  
    initString = m_geoSvc->constantAsString(m_encodingStringVariable.value());
    dd4hep::DDSegmentation::BitFieldCoder cellid_decoder(initString); 

    int nSimHits = inputSim.size();
    
    log << MSG::DEBUG << " processing collection with " <<  nSimHits  << " hits ... " << endmsg;
    
    for(int i=0; i< nSimHits; ++i){

      edm4hep::SimTrackerHit simTHit = inputSim.at( i );

      m_h[hitE]->Fill( simTHit.getEDep() * (dd4hep::GeV / dd4hep::keV) );

      if( simTHit.getEDep() < m_minEnergy ) {
        log << MSG::DEBUG << "Hit with insufficient energy " << simTHit.getEDep() * (dd4hep::GeV / dd4hep::keV) << " keV" << endmsg;
        continue;
      }
      
      const int cellID = simTHit.getCellID() ;
  
      //***********************************************************
      // get the measurement surface for this hit using the CellID
      //***********************************************************
      
      dd4hep::rec::SurfaceMap::const_iterator sI = m_map->find( cellID ) ;

      if( sI == m_map->end() ){
        log << MSG::ERROR << " DDPlanarDigiProcessor::processEvent(): no surface found for cellID : " 
            <<   cellid_decoder.fieldDescription() << "\n" << cellid_decoder.valueString(simTHit.getCellID()) << endmsg;
	      continue;
      }

      const dd4hep::rec::ISurface* surf = sI->second;
      int layer  = cellid_decoder.get(simTHit.getCellID(), "layer");

      dd4hep::rec::Vector3D oldPos( simTHit.getPosition().x, simTHit.getPosition().y, simTHit.getPosition().z );
      dd4hep::rec::Vector3D newPos;

      //************************************************************
      // Check if Hit is inside sensitive 
      //************************************************************
      
      if ( ! surf->insideBounds( dd4hep::mm * oldPos ) ) {
        log << MSG::DEBUG << "  hit at " << oldPos 
            << " " << cellid_decoder.fieldDescription() << "\n" << cellid_decoder.valueString(simTHit.getCellID())
            << " is not on surface "
            << *surf
            << " distance: " << surf->distance(  dd4hep::mm * oldPos )
            << endmsg;
        
        if( m_forceHitsOntoSurface ){
          
          dd4hep::rec::Vector2D lv = surf->globalToLocal( dd4hep::mm * oldPos  ) ;
          dd4hep::rec::Vector3D oldPosOnSurf = (1./dd4hep::mm) * surf->localToGlobal( lv ) ; 
          
          log << MSG::DEBUG << " moved to " << oldPosOnSurf << " distance "
              << (oldPosOnSurf-oldPos).r() << endmsg;        
          
          oldPos = oldPosOnSurf;

        } else {
          ++nDismissedHits;
          continue; 
        }
      }

      //***************************************************************
      // Smear time of the hit and apply the time window cut if needed
      //***************************************************************
      
      float hitT = simTHit.getTime();
      
      // Smearing time of the hit
      if (m_resT.value().size() and m_resT.value()[0] > 0.0) {
        float resT = m_resT.value().size() > 1 ? m_resT.value().at(layer) : m_resT.value().at(0);
        float tSmear = resT > 0.0 ? ( m_rng.shoot() *resT ) : 0.0;
        m_h[hT]->Fill( resT > 0.0 ? tSmear / resT : 0.0 );
        m_h[diffT]->Fill( tSmear );

        hitT += tSmear;
        log << MSG::DEBUG << "smeared hit at T: " << simTHit.getTime() << " ns to T: " << hitT
            << " ns according to resolution: " << resT << " ns" << endmsg;
      }
     
      // Correcting for the propagation time
      if (m_correctTimesForPropagation) {
        double dt = oldPos.r() / ( TMath::C() / 1e6 );
        hitT -= dt;
        log << MSG::DEBUG << "corrected hit at R: " << oldPos.r() << " mm by propagation time: "
            << dt << " ns to T: " << hitT << " ns" << endmsg;
      }
      
      // Skipping the hit if its time is outside the acceptance time window
      if (m_useTimeWindow) {
        float timeWindow_min = m_timeWindow_min.value().size() > 1 ? m_timeWindow_min.value().at(layer) : m_timeWindow_min.value().at(0);
        float timeWindow_max = m_timeWindow_max.value().size() > 1 ? m_timeWindow_max.value().at(layer) : m_timeWindow_max.value().at(0);
        if ( hitT < timeWindow_min || hitT > timeWindow_max ) {
          log << MSG::DEBUG << "hit at T: " << simTHit.getTime() << " smeared to: "
              << hitT << " is outside the time window: hit dropped"  << endmsg;
          ++nDismissedHits;
          continue;
        }
      }


      //*********************************************************************************
      // Try to smear the hit position but ensure the hit is inside the sensitive region
      //*********************************************************************************
      
      dd4hep::rec::Vector3D u = surf->u() ;
      dd4hep::rec::Vector3D v = surf->v() ;
      

      // get local coordinates on surface
      dd4hep::rec::Vector2D lv = surf->globalToLocal( dd4hep::mm * oldPos  ) ;
      double uL = lv[0] / dd4hep::mm ;
      double vL = lv[1] / dd4hep::mm ;

      bool accept_hit = false ;
      unsigned  tries   =  0 ;              
      static const unsigned MaxTries = 10 ; 
      
      float resU = ( m_resU.value().size() > 1 ?  m_resU.value().at(  layer )  : m_resU.value().at(0)  ) ;
      float resV = ( m_resV.value().size() > 1 ?  m_resV.value().at(  layer )  : m_resV.value().at(0)  ) ; 


      while( tries < MaxTries ) {
        
        if( tries > 0 ) {
          log << MSG::DEBUG << "retry smearing for " <<  cellid_decoder.fieldDescription() << "\n" << cellid_decoder.valueString(simTHit.getCellID()) << " : retries " << tries << endmsg;
        }

        double uSmear = m_rng.shoot() *resU ;
        double vSmear = m_rng.shoot() *resV ;

        
        // dd4hep::rec::Vector3D newPosTmp = oldPos +  uSmear * u ;  
        // if( ! _isStrip )  newPosTmp = newPosTmp +  vSmear * v ;  
        
        dd4hep::rec::Vector3D newPosTmp;
        if (m_isStrip){
            if (m_subDetName == "SET"){
                double xStripPos, yStripPos, zStripPos;
                //Find intersection of the strip with the z=centerOfSensor plane to set it as the center of the SET strip
                dd4hep::rec::Vector3D simHitPosSmeared = (1./dd4hep::mm) * ( surf->localToGlobal( dd4hep::rec::Vector2D( (uL+uSmear)*dd4hep::mm, 0.) ) );
                zStripPos = surf->origin()[2] / dd4hep::mm ;
                double lineParam = (zStripPos - simHitPosSmeared[2])/v[2];
                xStripPos = simHitPosSmeared[0] + lineParam*v[0];
                yStripPos = simHitPosSmeared[1] + lineParam*v[1];
                newPosTmp = dd4hep::rec::Vector3D(xStripPos, yStripPos, zStripPos);    
            } else {
                newPosTmp = (1./dd4hep::mm) * ( surf->localToGlobal( dd4hep::rec::Vector2D( (uL+uSmear)*dd4hep::mm, 0. ) ) );    
            }
        } else {
            newPosTmp = (1./dd4hep::mm) * ( surf->localToGlobal( dd4hep::rec::Vector2D( (uL+uSmear)*dd4hep::mm, (vL+vSmear)*dd4hep::mm ) ) );
        }

        log << MSG::DEBUG << " hit at    : " << oldPos
            << " smeared to: " << newPosTmp
            << " uL: "         << uL
            << " vL: "         << vL
            << " uSmear: "     << uSmear
            << " vSmear: "     << vSmear << endmsg;


        if ( surf->insideBounds( dd4hep::mm * newPosTmp ) ) {    
          
          accept_hit = true ;
          newPos     = newPosTmp ;

          m_h[hu]->Fill(  uSmear / resU ) ; 
          m_h[hv]->Fill(  vSmear / resV ) ; 

          m_h[diffu]->Fill( uSmear );
          m_h[diffv]->Fill( vSmear );

          break;  

        } else { 
          log << MSG::DEBUG << "  hit at " << newPosTmp 
              << " " << cellid_decoder.fieldDescription() << "\n" << cellid_decoder.valueString(simTHit.getCellID()) 
              << " is not on surface "                  
              << " distance: " << surf->distance( dd4hep::mm * newPosTmp ) 
              << endmsg;
        }
        
        ++tries;
      }
      
      if( accept_hit == false ) {
        log << MSG::DEBUG << "hit could not be smeared within ladder after " << MaxTries << "  tries: hit dropped"  << endmsg;
        ++nDismissedHits;
        continue; 
      } 
      
      //**************************************************************************
      // Store hit variables to TrackerHitPlaneImpl
      //**************************************************************************
      

      edm4hep::MutableTrackerHitPlane trkHit = trkhitCol.create();
                  
      //const int cellID1 = simTHit.getCellID1();
      trkHit.setCellID( cellID );
      //trkHit->setCellID1( cellID1 );
      
      trkHit.setPosition( newPos.const_array() );
      trkHit.setTime( hitT );
      trkHit.setEDep( simTHit.getEDep() );

      edm4hep::Vector2f u_direction;
      u_direction.a = u.theta();
      u_direction.b = u.phi();
      
      edm4hep::Vector2f v_direction;
      v_direction.a = v.theta();
      v_direction.b = v.phi();
      
      log << MSG::DEBUG << " U[0] = "<< u_direction.a << " U[1] = "<< u_direction.b 
                        << " V[0] = "<< v_direction.a << " V[1] = "<< v_direction.b
                        << endmsg;

      trkHit.setU( u_direction );
      trkHit.setV( v_direction );
      
      trkHit.setDu( resU );

      if( m_isStrip ) {

        // store the resolution from the length of the wafer - in case a fitter might want to treat this as 2d hit ....
        double stripRes = (surf->length_along_v() / dd4hep::mm ) / std::sqrt( 12. );
        trkHit.setDv( stripRes ); 

      } else {
        trkHit.setDv( resV );
      }

      if( m_isStrip ){
        ///TODO: UTIL::ILDTrkHitTrpyBit::ONEDIMENSIONAL = 29, couldn't find in edm4hep
        trkHit.setType( edm4hep::utils::setBit( trkHit.getType() ,  29, true ) );
      }

      //**************************************************************************
      // Set Relation to SimTrackerHit
      //**************************************************************************    

      // Set relation with LCRelationNavigator
      edm4hep::MutableTrackerHitSimTrackerHitLink rel = relCollection.create();
      rel.setFrom(trkHit);
      rel.setTo(simTHit);
      
      ++nCreatedHits;
      
      log << MSG::DEBUG << "-------------------------------------------------------" << endmsg;
      
    }
    
    // Filling the fraction of accepted hits in the event
    float accFraction = nSimHits > 0 ? float(nCreatedHits) / float(nSimHits) * 100.0 : 0.0;
    m_h[hitsAccepted]->Fill( accFraction );
    
    //**************************************************************************
    // Add collection to event
    //**************************************************************************    
        
    log << MSG::DEBUG << "Created " << nCreatedHits << " hits, " << nDismissedHits << " hits  dismissed\n" << endmsg;

    return std::make_tuple(std::move(trkhitCol), std::move(relCollection)); 
}


StatusCode DDPlanarDigiAlgorithm::finalize(){ 
  return StatusCode::SUCCESS;
}
