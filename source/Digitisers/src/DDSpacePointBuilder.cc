#include "DDSpacePointBuilder.h"

// EDM4HEP
#include <edm4hep/TrackerHit.h>
#include <edm4hep/TrackerHitPlane.h>
#include <edm4hep/SimTrackerHit.h>
#include <edm4hep/TrackerHitSimTrackerHitLink.h>
#include <edm4hep/MutableTrackerHit.h>
#include <edm4hep/MutableTrackerHitPlane.h>

// DD4HEP
#include "DDRec/DetectorData.h"

// CLHEP
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h"

// Standard
#include <cmath>
#include <sstream>



DECLARE_COMPONENT(DDSpacePointBuilder);


DDSpacePointBuilder::DDSpacePointBuilder(const std::string& name, ISvcLocator* svcLoc) : MultiTransformer(name, svcLoc, {
	KeyValues("TrackerHitCollection", {"FTDTrackerHits"}),
	KeyValues("TrackerHitSimHitRelCollection", {"FTDTrackerHitRelations"})}, {
	KeyValues("SpacePointsCollection", {"FTDSpacePoints"}),
	KeyValues("SimHitSpacePointRelCollection", {"FTDSimHitSpacepointRelations"})
}) {}


StatusCode DDSpacePointBuilder::initialize() { 
  /*
  MarlinTrk::IMarlinTrkSystem* trksystem =  MarlinTrk::Factory::createMarlinTrkSystem( "DDKalTest" , 0, "" ) ;
  
  
  if( trksystem == 0 ) {
    
    throw EVENT::Exception( std::string("  Cannot initialize MarlinTrkSystem of Type: ") + std::string("DDKalTest" )  ) ;
    
  }
  
  trksystem->init() ;  
  */
  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();
  //theDetector = dd4hep::Detector::getInstance();
  
  //===========  get the surface map from the SurfaceManager ================
  
  dd4hep::rec::SurfaceManager& surfMan = *theDetector.extension<dd4hep::rec::SurfaceManager>() ;
  dd4hep::DetElement det = theDetector.detector( m_subDetName ) ;
  //const dd4hep::rec::SurfaceMap *surfMap ;
  surfMap = surfMan.map( det.name() ) ;
  
 
  /*
  // alternative way to create a map of instances
  const dd4hep::rec::IGeometry& geom = dd4hep::rec::IGeometry::instance() ;
  const std::vector<const dd4hep::rec::ISurface*>& surfaces = geom.getSurfaces() ;
  // create map of surfaces
  std::map< long, const dd4hep::rec::ISurface* > surfMap ;
    for(std::vector<const dd4hep::rec::ISurface*>::const_iterator surf = surfaces.begin() ; surf != surfaces.end() ; ++surf){
    surfMap[(*surf)->id() ] = (*surf) ;
  }
  */

  return StatusCode::SUCCESS; 
}


std::tuple<edm4hep::TrackerHitCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection> operator(
        const edm4hep::TrackerHitCollection& inputHits,
        const edm4hep::TrackerHitSimTrackerHitLinkCollection& inputRels) const{
 
    unsigned createdSpacePoints = 0;
    unsigned rawStripHits = 0;
    unsigned possibleSpacePoints = 0;
    m_nOutOfBoundary = 0;
    m_nStripsTooParallel = 0;
    m_nPlanesNotParallel = 0;
    
    // Map TrackerHits to SimTrackerHits
    std::unordered_map<edm4hep::TrackerHit, std::vector<edm4hep::SimTrackerHit>> trackerHit2SimHit;
    for (const auto& hitRel : inputRels) {
      edm4hep::TrackerHit trackerHit = hitRel.getFrom();
      edm4hep::SimTrackerHit simTrackerHit = hitRel.getTo();
      trackerHit2SimHit[trackerHit].push_back(simTrackerHit);
    } 
   
    edm4hep::TrackerHitCollection spCol;    // output spacepoint collection

    // Relation navigator for creating SpacePoint - SimTrackerHit relations
    edm4hep::TrackerHitSimTrackerHitLinkCollection spRelCollection;

    unsigned nHits = inputHits.size();
    
    debug() << "Number of hits: " << nHits << endmsg;
    
    //store hits in map according to their CellID
    std::map< int , std::vector< edm4hep::TrackerHitPlane* > > map_cellID_hits;
    std::map< int , std::vector< edm4hep::TrackerHitPlane* > >::iterator it;
    
    for( unsigned i=0; i<nHits; i++){
      
      edm4hep::TrackerHitPlane trkHit = dynamic_cast<TrackerHitPlane>(inputHits.at( i ));

      if( trkHit != NULL) {
        debug() << "Add hit with CellID = " << trkHit.getCellID() << " " << getCellIDInfo( trkHit.getCellID() ) << endmsg;
        map_cellID_hits[ trkHit.getCellID() ].push_back( &trkHit );
      }
    }
    
    // now loop over all CellIDs
    for( it= map_cellID_hits.begin(); it!= map_cellID_hits.end(); it++ ){
      
      rawStripHits += it->second.size();
      
      std::vector< edm4hep::TrackerHitPlane* > hitsFront = it->second;
  
      int cellID = it->first;
     
      //get the CellIDs at the back of this sensor
      std::vector< int > cellIDsBack = getCellIDsAtBack( cellID );

      for( unsigned i=0; i< cellIDsBack.size(); i++ ){ 
        
        
        int cellIDBack = cellIDsBack[i];
        std::vector< TrackerHitPlane* > hitsBack = map_cellID_hits[ cellIDBack ];
        
	debug() << "strips: CellID " << cellID  << " " << getCellIDInfo( cellID )  << "(" << hitsFront.size()
		<< " hits) <---> CellID " << cellIDBack << getCellIDInfo( cellIDBack )
		<< "(" << hitsBack.size() << " hits)\n"
		<< "--> " << hitsFront.size() * hitsBack.size() << " possible combinations\n"
                << endmsg;
        
        possibleSpacePoints += hitsFront.size() * hitsBack.size();
        
        
        // Now iterate over all combinations and store those that make sense
        for( unsigned ifront=0; ifront<hitsFront.size(); ifront++ ){
          
          TrackerHitPlane* hitFront = hitsFront[ifront];
          
          for( unsigned j=0; j<hitsBack.size(); j++ ){
            
            
            TrackerHitPlane* hitBack = hitsBack[j];

            std::vector<edm4hep::SimTrackerHit> simHitsFront =  trackerHit2SimHit[ hitFront ];
            std::vector<edm4hep::SimTrackerHit> simHitsBack  = trackerHit2SimHit[ hitBack ];

            debug() << "attempt to create space point from:\n" << " front hit: " 
                    << hitFront << " no. of simhit = " << simHitsFront.size() << endmsg;
            if( simHitsFront.empty() == false ) { 
              edm4hep::SimTrackerHit simhit = simHitsFront.at(0);
              debug() << " first simhit = " << simhit << " mcp = "<< simhit.getParticle() 
                      << " ( " << simhit.getPosition().x << " " 
                               << simhit.getPosition().y << " " 
                               << simhitgetPosition().z << " ) " << endmsg; 
            }
            debug() << "  rear hit: " << hitBack << " no. of simhit = " << simHitsBack.size() ;
            if( simHitsBack.empty() == false ) { 
              edm4hep::SimTrackerHit simhit = simHitsBack.at(0);
              debug() << " first simhit = " << simhit << " mcp = "<< simhit.getParticle()
                      << " ( " << simhit.getPosition().x << " " 
                               << simhit.getPosition().y << " " 
                               << simhit.getPosition().z << " ) " << endmsg; 
            }
            
            bool ghost_hit = true;
            
            if (simHitsFront.size()==1 && simHitsBack.size() == 1) {
              debug() << "SpacePoint creation from two good hits:" << endmsg;
              ghost_hit = simHitsFront.at(0).getParticle() != simHitsBack.at(0).getMCParticle(); 
            }
            
            if ( ghost_hit == true ) {
              debug() << "SpacePoint Ghosthit!" << endmsg;
            }

            double strip_length_mm = 0;
	    strip_length_mm = m_striplength ;

            // add tolerence 
            strip_length_mm = strip_length_mm * (1.0 + m_striplength_tolerance);
            
	    edm4hep::MutableTrackerHitPlane spacePoint = createSpacePoint( hitFront, hitBack, strip_length_mm);

            if ( spacePoint != NULL ) {

              spacePoint.setCellID(cellID);
              
              // store the hits it's composed of:
              //spacePoint->rawHits().push_back( hitFront );
              //spacePoint->rawHits().push_back( hitBack );
              
              //spacePoint->setType( UTIL::set_bit( spacePoint->getType() ,  ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ) ) ;
              
              spCol.push_back( spacePoint ) ; 
              
              createdSpacePoints++;
              
              
              ///////////////////////////////
              // make the relations
              if( simHitsFront.size() == 1 ){
                
                edm4hep::SimTrackerHit simHit = simHitsFront[0];
                
                if( simHit != NULL ){
                  edm4hep::MutableTrackerHitSimTrackerHitLink link = spRelCollection.create();
                  link.setFrom(spacePoint);
                  link.setTo(simHit);
                  link.setWeight(0.5);
                }
              }
              if( simHitsBack.size() == 1 ){
                
                edm4hep::SimTrackerHit simHit = simHitsBack[0];
                
                if( simHit != NULL ){
                  edm4hep::MutableTrackerHitSimTrackerHitLink link = spRelCollection.create();
                  link.setFrom(spacePoint);
                  link.setTo(simHit);
                  link.setWeight(0.5);
                }
              }
            } else {
                 
              if ( ghost_hit == true ) {
                debug() << "Ghosthit correctly rejected" << endmsg;
              } else {
                debug() << "True hit rejected!" << endmsg;
              }
              
               //////////////////////////////////
            }
            
          }
          
        }
        
      }
      
    }
    
    debug() << "\nCreated " << createdSpacePoints
            << " space points ( raw strip hits: " << rawStripHits << ")\n\n"
            
            << "  There were " << rawStripHits << " strip hits available, giving " 
            << possibleSpacePoints << " possible space points\n\n"
    
            << "  " << _nStripsTooParallel << " space points couldn't be created, "
                                           << "because the strips were too parallel\n\n"
            << "  " << _nPlanesNotParallel << " space points couldn't be created, "
                    << "because the planes of the measurement surfaces where not parallel enough\n\n"
            << "  " << _nOutOfBoundary     << " space points couldn't be created, "
                                           << "because the result was outside the sensor boundary\n\n" 
            << endmsg;

   return std::make_tuple(std::move(spCol, std::move(spRelCollection)); 
}


StatusCode DDSpacePointBuilder::finalize(){
  return StatusCode::SUCCESS;
}

//TrackerHitImpl* DDSpacePointBuilder::createSpacePoint( TrackerHitPlane* a , TrackerHitPlane* b, double stripLength, const dd4hep::rec::SurfaceMap* surfMap ){
edm4hep::MutableTrackerHitPlane DDSpacePointBuilder::createSpacePoint( TrackerHitPlane* a , TrackerHitPlane* b, double stripLength ){  
  edm4hep::Vector3d pa = a->getPosition();
  CLHEP::Hep3Vector PA( pa.x, pa.y, pa.z );
  dd4hep::rec::Vector3D ddPA( pa.x * dd4hep::mm, pa.y * dd4hep::mm, pa.z * dd4hep::mm );
  float du_a = a->getdU();
  
  //const dd4hep::rec::ISurface* msA = surfMap[a->getCellID()];
  const dd4hep::rec::ISurface* msA = surfMap->find(a->getCellID())->second;
  debug() << " Do I find a surface " << *msA << endmsg;
  dd4hep::rec::Vector3D ddWA = msA->normal();
  dd4hep::rec::Vector3D ddUA = msA->u();
  dd4hep::rec::Vector3D ddVA = msA->v();
  CLHEP::Hep3Vector UA(ddUA.x() / dd4hep::mm, ddUA.y() / dd4hep::mm, ddUA.z() / dd4hep::mm);
  CLHEP::Hep3Vector VA(ddVA.x() / dd4hep::mm, ddVA.y() / dd4hep::mm, ddVA.z() / dd4hep::mm);
  CLHEP::Hep3Vector WA(ddWA.x() / dd4hep::mm, ddWA.y() / dd4hep::mm, ddWA.z() / dd4hep::mm); 
  
  edm4hep::Vector3d pb = b->getPosition();
  CLHEP::Hep3Vector PB( pb.x, pb.y, pb.z );
  dd4hep::rec::Vector3D ddPB( pb.x * dd4hep::mm, pb.y * dd4hep::mm, pb.z * dd4hep::mm );
  float du_b = b->getdU();
  
  const dd4hep::rec::ISurface* msB = surfMap->find(b->getCellID())->second;
  dd4hep::rec::Vector3D ddWB = msB->normal();
  dd4hep::rec::Vector3D ddUB = msB->u();
  dd4hep::rec::Vector3D ddVB = msB->v();

  CLHEP::Hep3Vector UB(ddUB.x() / dd4hep::mm, ddUB.y() / dd4hep::mm, ddUB.z() / dd4hep::mm);
  CLHEP::Hep3Vector VB(ddVB.x() / dd4hep::mm, ddVB.y() / dd4hep::mm, ddVB.z() / dd4hep::mm);
  CLHEP::Hep3Vector WB(ddWB.x() / dd4hep::mm, ddWB.y() / dd4hep::mm, ddWB.z() / dd4hep::mm);
  
  debug() << "\t ( " << pa.x << " " 
                     << pa.y << " " 
                     << pa.z << " ) <--> ( " 
                     << pb.x << " " 
                     << pb.y << " " 
                     << pb.z << " )" << endmsg;

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // First: check if the two measurement surfaces are parallel (i.e. the w are parallel or antiparallel)
  double angle = fabs(WB.angle(WA));

  double angleMax = 1.*M_PI/180.;
  if(( angle > angleMax )&&( angle < M_PI-angleMax )){
    
    m_nPlanesNotParallel++;
    debug() << "\tThe planes of the measurement surfaces are not parallel enough, "
            << "the angle between the W vectors is " << angle
            << " where the angle has to be smaller than " << angleMax 
            << " or bigger than " << M_PI-angleMax << endmsg;
    return NULL; //calculate the xing point and if that fails don't create a spacepoint
    
  }
  //////////////////////////////////////////////////////////////////////////////////////////////////////


  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Next: check if the angle between the strips is not 0
  angle = fabs(VB.angle(VA));
  double angleMin= 1.*M_PI/180.;
  if(( angle < angleMin )||( angle > M_PI-angleMin )){
    
    m_nStripsTooParallel++;
    debug() << "\tThe strips (V vectors) of the measurement surfaces are too parallel, "
            << "the angle between the V vectors is " << angle
            << " where the angle has to be between " << angleMax 
            << " or bigger than " << M_PI-angleMin << endmsg;
    return NULL; //calculate the xing point and if that fails don't create a spacepoint
    
  }
  //////////////////////////////////////////////////////////////////////////////////////////////////////
 

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Next we want to calculate the crossing point.
  
  CLHEP::Hep3Vector point;

  CLHEP::Hep3Vector vertex(0.,0.,0.);
  dd4hep::rec::Vector2D L1 = msA->globalToLocal(ddPA);
  dd4hep::rec::Vector2D L2 = msB->globalToLocal(ddPB);

  dd4hep::rec::Vector2D ddSL1, ddEL1, ddSL2, ddEL2;
  if (_subDetName == "SET"){
      ddSL1 = dd4hep::rec::Vector2D( L1.u(), L1.v() + (-stripLength * dd4hep::mm)/2.0 );
      ddEL1 = dd4hep::rec::Vector2D( L1.u(), L1.v() + (stripLength * dd4hep::mm)/2.0 );
      ddSL2 = dd4hep::rec::Vector2D( L2.u(), L2.v() + (-stripLength * dd4hep::mm)/2.0 );
      ddEL2 = dd4hep::rec::Vector2D( L2.u(), L2.v() + (stripLength * dd4hep::mm)/2.0 );
  }
  else{
      ddSL1 = dd4hep::rec::Vector2D( L1.u(), (-stripLength * dd4hep::mm)/2.0 );
      ddEL1 = dd4hep::rec::Vector2D( L1.u(), (stripLength * dd4hep::mm)/2.0 );
      ddSL2 = dd4hep::rec::Vector2D( L2.u(), (-stripLength * dd4hep::mm)/2.0 );
      ddEL2 = dd4hep::rec::Vector2D( L2.u(), (stripLength * dd4hep::mm)/2.0 );        
  }

  dd4hep::rec::Vector3D ddS1 = msA->localToGlobal(ddSL1);
  dd4hep::rec::Vector3D ddE1 = msA->localToGlobal(ddEL1);
  dd4hep::rec::Vector3D ddS2 = msB->localToGlobal(ddSL2);
  dd4hep::rec::Vector3D ddE2 = msB->localToGlobal(ddEL2);
  CLHEP::Hep3Vector S1 (ddS1.x() / dd4hep::mm, ddS1.y() / dd4hep::mm, ddS1.z() / dd4hep::mm);
  CLHEP::Hep3Vector E1 (ddE1.x() / dd4hep::mm, ddE1.y() / dd4hep::mm, ddE1.z() / dd4hep::mm);
  CLHEP::Hep3Vector S2 (ddS2.x() / dd4hep::mm, ddS2.y() / dd4hep::mm, ddS2.z() / dd4hep::mm);
  CLHEP::Hep3Vector E2 (ddE2.x() / dd4hep::mm, ddE2.y() / dd4hep::mm, ddE2.z() / dd4hep::mm);

  debug() << " stripLength = " << stripLength << "\n"
          << " S1 = " << S1 << "\n"
          << " E1 = " << E1 << "\n"
          << " S2 = " << S2 << "\n"
          << " E2 = " << E2 << endmsg;

  point.set(0.0, 0.0, 0.0);
  
  
  int valid_intersection = calculatePointBetweenTwoLines_UsingVertex( S1, E1, S2, E2, vertex, point );
  
  if (valid_intersection != 0) {
    debug() << "\tNo valid intersection for lines." << endmsg;
    return NULL;
  }
  
  debug() << "\tVertex: Position of space point (global) : ( " << point.x() << " " 
                                                               << point.y() << " " 
                                                               << point.z() << " )" << endmsg;
  
  // using dd4hep to check if hit within boundaries
  dd4hep::rec::Vector3D DDpoint( point.x() * dd4hep::mm, point.y() * dd4hep::mm, point.z() * dd4hep::mm );

  if ( !msA->insideBounds(DDpoint)){

    m_nOutOfBoundary++;
    debug() << " SpacePoint position lies outside the boundary of the layer." << endmsg;
    
    return NULL;
  }

  
  //Create the new TrackerHit
  edm4hep::MutableTrackerHitPlane spacePoint = new edm4hep::MutableTrackerHitPlane();
  
  edm4hep::Vector3d pos{point.x(), point.y(), point.z() };
  spacePoint.setPosition( pos );
  
  
  // set error treating the strips as stereo with equal and opposite rotation -- for reference see Karimaki NIM A 374 p367-370

  // first calculate the covariance matrix in the cartisian coordinate system defined by the sensor 
  // here we assume that du is the same for both sides
  
  if( fabs(du_a - du_b) > 1.0e-06 ){
    error() << "\tThe measurement errors of the two 1D hits must be equal" << endmsg;    
    assert( (fabs(du_a - du_b) > 1.0e-06) == false );
    return NULL; //measurement errors are not equal don't create a spacepoint
  }
 
  
  double du2 = du_a*du_a;
  
  // rotate the strip system back to double-layer wafer system
  CLHEP::Hep3Vector u_sensor = UA + UB;
  CLHEP::Hep3Vector v_sensor = VA + VB;
  CLHEP::Hep3Vector w_sensor = WA + WB;
  
  CLHEP::HepRotation rot_sensor( u_sensor, v_sensor, w_sensor );
  CLHEP::HepMatrix rot_sensor_matrix;
  rot_sensor_matrix = rot_sensor;
  
  double cos2_alpha = VA.cos2Theta(v_sensor) ; // alpha = strip angle   
  double sin2_alpha = 1 - cos2_alpha ; 
  
  CLHEP::HepSymMatrix cov_plane(3,0); // u,v,w
  
  cov_plane(1,1) = (0.5 * du2) / cos2_alpha;
  cov_plane(2,2) = (0.5 * du2) / sin2_alpha;
  
  CLHEP::HepSymMatrix cov_xyz= cov_plane.similarity(rot_sensor_matrix);
  
  debug() << "\t cov_plane  = " << cov_plane << "\n\n"  
          << "\tstrip_angle = " << VA.angle(VB)/(M_PI/180) / 2.0 << " degrees\n\n" 
          << "\t cov_xyz  = " << cov_xyz << endmsg;
  
  edm4hep::CovMatrix3f cov( 9 );
  
  for(int irow=0; irow<3; ++irow ){
    for(int jcol=0; jcol<irow+1; ++jcol){
      cov.setValue(cov_xyz[irow][jcol], irow, jcol);
    }
  }
  
  spacePoint.setCovMatrix(cov);

  const auto pointTime = std::min(a->getTime(), b->getTime());
  spacePoint.setTime(pointTime);

  debug() << "\tHit accepted" << endmsg;
  
  return spacePoint;
  
}



int DDSpacePointBuilder::calculatePointBetweenTwoLines_UsingVertex( 
                                                const CLHEP::Hep3Vector& PA, 
                                                const CLHEP::Hep3Vector& PB, 
                                                const CLHEP::Hep3Vector& PC, 
                                                const CLHEP::Hep3Vector& PD,
                                                const CLHEP::Hep3Vector& Vertex,
                                                CLHEP::Hep3Vector& point){

  
  // A general point on the line joining point PA to point PB is
  // x, where 2*x=(1+m)*PA + (1-m)*PB. Similarly for 2*y=(1+n)*PC + (1-n)*PD.
  // Suppose that v is the vertex. Requiring that the two 'general
  // points' lie on a straight through v means that the vector x-v is a 
  // multiple of y-v. This condition fixes the parameters m and n.
  // We then return the 'space-point' x, supposed to be the layer containing PA and PB. 
  // We require that -1<m<1, otherwise x lies 
  // outside the segment PA to PB; and similarly for n.
  
  bool ok = true;
  
//  debug() << " Vertex = " << Vertex << endmsg; 
//  
//  debug() << " PA = " << PA << endmsg;
//  debug() << " PB = " << PB << endmsg;
//  debug() << " PC = " << PC << endmsg;
//  debug() << " PD = " << PD << endmsg;
  
  CLHEP::Hep3Vector VAB(PA-PB);
  CLHEP::Hep3Vector VCD(PC-PD);

//  debug() << " VAB = " << VAB << endmsg;
//  debug() << " VCD = " << VCD << endmsg;
  
  CLHEP::Hep3Vector  s(PA+PB-2*Vertex);   // twice the vector from vertex to midpoint
  CLHEP::Hep3Vector  t(PC+PD-2*Vertex);   // twice the vector from vertex to midpoint

  CLHEP::Hep3Vector  qs(VAB.cross(s));  
  CLHEP::Hep3Vector  rt(VCD.cross(t));  

//  debug() << " s = " << s << endmsg;
//  debug() << " t = " << t << endmsg;
//  debug() << " qs = " << qs << endmsg;
//  debug() << " rt = " << rt << endmsg;
  
  
  double m = (-(s*rt)/(VAB*rt)); // ratio for first line
    
  double limit = 1.0;
  
  if (m>limit || m<-1.*limit) {
    
    debug() << "m' = " << m << endmsg;
    
    ok = false;
    
  } else {
    
    double n = (-(t*qs)/(VCD*qs)); // ratio for second line

	  if (n>limit || n<-1.*limit) {
  
      debug() << "n' = " << n << endmsg;
      
      ok = false;

    }
  }
  
  if (ok) {
    point = 0.5*(PA + PB + m*VAB);
  }
  
  return ok ? 0 : 1;
  
}



int DDSpacePointBuilder::calculatePointBetweenTwoLines( const CLHEP::Hep3Vector& P1, const CLHEP::Hep3Vector& V1, const CLHEP::Hep3Vector& P2, const CLHEP::Hep3Vector& V2, CLHEP::Hep3Vector& point ){
  
  // Richgungsvektor normal auf die anderen beiden:
  CLHEP::Hep3Vector n = V1.cross( V2 );
  
  // Now we want to rotate into a coordinate system, where n is parallel to the z axis
  // For this: first set phi to 0
  // then: set theta to 0 (we set phi to 0 first, so we can then rotate arount the y axis)
  CLHEP::HepRotation rot;
  rot.rotateZ( -n.phi() );
  CLHEP::Hep3Vector nPrime = rot * n; //now the phi of nPrime should be 0
  debug() << "phi of n' = " << nPrime.phi() << " (it should be 0!!!)" << endmsg;
  rot.rotateY( -n.theta() );
  nPrime = rot * n;
  debug() << "phi of n'' = " << nPrime.phi() << " (it should be 0!!!)" << endmsg;
  debug() << "theta of n'' = " << nPrime.theta() <<  " (it should be 0!!!)" << endmsg;
  
  // Now rotate all the vectors and points into this coordinatesystem.
  CLHEP::Hep3Vector P1prime = rot * P1;
  CLHEP::Hep3Vector V1prime = rot * V1;
  CLHEP::Hep3Vector P2prime = rot * P2;
  CLHEP::Hep3Vector V2prime = rot * V2;
  
  // What is the gain of rotating into this system?
  // A: 
  double x;
  double y;
  int res = calculateXingPoint( P1prime.x(), P1prime.y(), V1prime.x(), V1prime.y(), P2prime.x(), P2prime.y(), V2prime.x(), V2prime.y(), x, y );
  
  if ( res != 0 ) return 1;
  
  point.setX( x );
  point.setY( y );
  point.setZ( (P1prime.z() + P2prime.z())/2. );
  
  // Now transform back to the global coordinates
  point = rot.inverse() * point;
  
  
  return 0;
  
}


int DDSpacePointBuilder::calculateXingPoint( double x1, double y1, float ex1, float ey1, double x2, double y2, float ex2, float ey2, double& x, double& y ){


  float a = (x1*ey1 - y1*ex1) - (x2*ey1 - y2*ex1);
  float b = ex2*ey1 - ex1*ey2;

  const float epsilon = 0.00001;

  if( fabs(b) < epsilon ) return 1; // if b==0 the two directions e1 and e2 are parallel and there is no crossing!

  float t = a/b;

  x = x2 + t*ex2;
  y = y2 + t*ey2;

  return 0; 

}
 
std::vector< int > DDSpacePointBuilder::getCellIDsAtBack( int cellID ){
  
  std::vector< int > back;
  
  //find out detector, layer
  ACTSTracking::BitField64 cellIDer( "system:5,side:-2,layer:6,module:11,sensor:8" );
  cellIDer.setValue( cellID );
  

  int subdet = cellIDer[ "system" ];
  int layer  = cellIDer[ "layer" ];
  
  if (subdet != ILDDetID::FTD)  {
    
    //check if sensor is in front
    if( layer%2 == 0 ){ // even layers are front sensors
      
      cellIDer[ "layer" ] = layer + 1; 
      // it is assumed that the even layers are the front layers
      // and the following odd ones the back layers
      
      back.push_back( cellID );
      
    }
  }

  else{

    dd4hep::Detector & theDetector2 = dd4hep::Detector::getInstance();
    dd4hep::DetElement ftdDE = theDetector2.detector( m_subDetName);
    dd4hep::rec::ZDiskPetalsData* ft = ftdDE.extension<dd4hep::rec::ZDiskPetalsData>();

    int sensor = cellIDer[ "sensor" ];
    //int Nsensors = ft->layers.at(layer).petalNumber ; 
    int Nsensors = ft->layers.at(layer).sensorsPerPetal ;

    debug() << " layer " << layer << " sensors " << Nsensors << "\n" 
            << " so sensor " << sensor << " is connected with sensor " 
            << sensor + Nsensors/2 << endmsg;

    std::vector<dd4hep::rec::ZDiskPetalsStruct::SensorType> Sensors ;
    
    //check if sensor is in front
    //if(( Sensors.at(layer).DoubleSided ) && ( sensor <= Nsensors / 2 ) ){
    if (sensor <= Nsensors / 2 ) {
      
      cellIDer[ "sensor" ] = sensor + Nsensors / 2; 
      // it is assumed, that sensors 1 until n/2 will be on front
      // and sensor n/2 + 1 until n are at the back
      // so the sensor x, will have sensor x+n/2 at the back
      
      back.push_back( cellID);
      
    }  

  }
TrackerCellID::subdet
  return back; 
}


std::string DDSpacePointBuilder::getCellIDInfo( int cellID ){

  std::stringstream s;
  
  //find out layer, module, sensor
  ACTSTracking::BitField64  cellIDer( "system:5,side:-2,layer:6,module:11,sensor:8" );
  cellIDer.setValue( cellID );

  int subdet = cellID[ "system" ] ;
  int side   = cellID[ "side" ];
  int module = cellID[ "module" ];
  int sensor = cellID[ "sensor" ];
  int layer  = cellID[ "layer" ];
  
  s << "(su" << subdet 
    << ",si" << side 
    << ",la" << layer 
    << ",mo" << module 
    << ",se" << sensor 
    << ")";
  
  return s.str(); 
}
