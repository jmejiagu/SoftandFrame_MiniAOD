// -*- C++ -*-
//
// Package:    JPsiBc
// Class:      JPsiBc
// 

//=================================================
// original author:  Jhovanny Andres Mejia        |
//         created:  Fryday Sep 23                |
//         <jhovanny.andres.mejia.guisao@cern.ch> | 
//=================================================

// system include files
#include <memory>


// user include files
#include "myAnalyzers/JPsiKsPAT/src/Bc2SPAT.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
//#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Common/interface/TriggerResults.h"

//For kinematic fit:
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"            
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"


#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TH2F.h"

//
// constants, enums and typedefs
//

  typedef math::Error<3>::type CovarianceMatrix;

//
// static data member definitions
//

//
// constructors and destructor
//

Bc2SPAT::Bc2SPAT(const edm::ParameterSet& iConfig)
  :
  dimuon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
  //trakCollection_label(consumes<std::vector<pat::GenericParticle>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  trakCollection_label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  //trakCollection_label(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("Trak"))),
  //trakCollection_label_lowpt(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak_lowpt"))),
  primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  BSLabel_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bslabel"))),

  
  OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
  //upsilon_mass_(iConfig.getParameter<double>("upsilon_mass")),
  //triggerCuts_(iConfig.getParameter<uint32_t>("triggerCuts"))
 
  tree_(0), 

  tri_Dim25(0), tri_Dim20(0), tri_JpsiTk(0),

  // *******************************************************
 
  nB(0), nMu(0),

  Bc2s_mass(0), Bc2s_px(0), Bc2s_py(0), Bc2s_pz(0),
  Bc2s_pi2_px_track(0), Bc2s_pi2_py_track(0), Bc2s_pi2_pz_track(0), 
  Bc2s_pi3_px_track(0), Bc2s_pi3_py_track(0), Bc2s_pi3_pz_track(0),
  Bc2s_pi2_charge1(0), Bc2s_pi3_charge1(0),
  Bc2s_Prob(0),

  B_mass(0), B_px(0), B_py(0), B_pz(0), B_charge(0),
  B_k_charge1(0), B_k_px_track(0), B_k_py_track(0), B_k_pz_track(0),
  
  B_J_mass(0), B_J_px(0), B_J_py(0), B_J_pz(0),

  B_J_px1(0), B_J_py1(0), B_J_pz1(0),
  B_J_px2(0), B_J_py2(0), B_J_pz2(0), 
  B_J_charge1(0), B_J_charge2(0),
  
  nVtx(0),
  priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXE(0), priVtxYE(0), priVtxZE(0), priVtxCL(0),
  priVtxXYE(0), priVtxXZE(0), priVtxYZE(0),

  pVtxIPX(0),   pVtxIPY(0),   pVtxIPZ(0), pVtxIPXE(0),   pVtxIPYE(0),   pVtxIPZE(0), pVtxIPCL(0),
  pVtxIPXYE(0),   pVtxIPXZE(0),   pVtxIPYZE(0),
 
  // ************************ ****************************************************

  B_Prob(0), B_J_Prob(0), 
 
  B_DecayVtxX(0),     B_DecayVtxY(0),     B_DecayVtxZ(0),
  B_DecayVtxXE(0),    B_DecayVtxYE(0),    B_DecayVtxZE(0),
  B_DecayVtxXYE(0),   B_DecayVtxXZE(0),   B_DecayVtxYZE(0),

  B_J_DecayVtxX(0),   B_J_DecayVtxY(0),   B_J_DecayVtxZ(0),
  B_J_DecayVtxXE(0),  B_J_DecayVtxYE(0),  B_J_DecayVtxZE(0),
  B_J_DecayVtxXYE(0), B_J_DecayVtxXZE(0), B_J_DecayVtxYZE(0),

  trigger(0),
 
  run(0), event(0),
  lumiblock(0)


{
   //now do what ever initialization is needed
}


Bc2SPAT::~Bc2SPAT()
{

}


//
// member functions
//



// ------------ method called to for each event  ------------
void Bc2SPAT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;
  
  //*********************************
  // Get event content information
  //*********************************  

  // Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB; 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); 

  edm::Handle< View<pat::PackedCandidate> > thePATTrackHandle;
  iEvent.getByToken(trakCollection_label,thePATTrackHandle);

  edm::Handle< View<pat::Muon> > thePATMuonHandle;
  iEvent.getByToken(dimuon_Label,thePATMuonHandle);

// Get HLT results
 edm::Handle<edm::TriggerResults> triggerResults_handle;
 //iEvent.getByLabel(triggerResults_Label, triggerResults_handle);
 iEvent.getByToken(triggerResults_Label, triggerResults_handle);

    trigger = 0;
    if (triggerResults_handle.isValid()) {
      const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
      unsigned int NTRIGGERS = 5;
      std::string TriggersToTest[NTRIGGERS] = {
	 "HLT_Dimuon25_Jpsi","HLT_Dimuon20_Jpsi_Barrel_Seagulls",
	 "HLT_DoubleMu4_JpsiTrk_Displaced","HLT_DoubleMu4_JpsiTrkTrk_Displaced",
	 "HLT_DoubleMu4_3_Jpsi_Displaced"};
      
      for (unsigned int i = 0; i < NTRIGGERS; i++) {
	for (int version = 1; version < 19; version++) {
	  std::stringstream ss;
	  ss << TriggersToTest[i] << "_v" << version;
	  unsigned int bit = TheTriggerNames.triggerIndex(edm::InputTag(ss.str()).label());
	  if (bit < triggerResults_handle->size() && triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
	    trigger += (1<<i);
	    break;
	  }
	}
      }
    } else std::cout << "*** NO triggerResults found " << iEvent.id().run() << "," << iEvent.id().event() << std::endl;
    
    
  //*********************************
  //Now we get the primary vertex 
  //*********************************

  reco::Vertex bestVtx;
  //edm::Handle<std::vector<reco::Vertex> > primaryVertices_handle;
  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  // get primary vertex
  bestVtx = *(primaryVertices_handle->begin());

  priVtxX = bestVtx.x();
  priVtxY = bestVtx.y();
  priVtxZ = bestVtx.z();
  //priVtxXE = bestVtx.xError();
  //priVtxYE = bestVtx.yError();
  //priVtxZE = bestVtx.zError();
  priVtxXE = bestVtx.covariance(0, 0);
  priVtxYE = bestVtx.covariance(1, 1);
  priVtxZE = bestVtx.covariance(2, 2);
  priVtxXYE = bestVtx.covariance(0, 1);
  priVtxXZE = bestVtx.covariance(0, 2);
  priVtxYZE = bestVtx.covariance(1, 2);

  priVtxCL = ChiSquaredProbability((double)(bestVtx.chi2()),(double)(bestVtx.ndof())); 
  nVtx = primaryVertices_handle->size(); 
 
  lumiblock = iEvent.id().luminosityBlock();
  run = iEvent.id().run();
  event = iEvent.id().event();
  //int run1   =  iEvent.id().run();
  //int event1 =  iEvent.id().event();


  //*****************************************
  //Let's begin by looking for J/psi+pi^+

  unsigned int nMu_tmp = thePATMuonHandle->size();
  nMu = nMu_tmp;

 for(View<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) 
    {
      
      for(View<pat::Muon>::const_iterator iMuon2 = iMuon1+1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) 
	{
	  if(iMuon1==iMuon2) continue;
	  
	  //opposite charge 
	  if( (iMuon1->charge())*(iMuon2->charge()) == 1) continue;

	  TrackRef glbTrackP;	  
	  TrackRef glbTrackM;	  
	  
	  if(iMuon1->charge() == 1){ glbTrackP = iMuon1->track();}
	  if(iMuon1->charge() == -1){ glbTrackM = iMuon1->track();}
	  
	  if(iMuon2->charge() == 1) { glbTrackP = iMuon2->track();}
	  if(iMuon2->charge() == -1){ glbTrackM = iMuon2->track();}
	  
	  if( glbTrackP.isNull() || glbTrackM.isNull() ) 
	    {
	      //std::cout << "continue due to no track ref" << endl;
	      continue;
	    }

	  if(iMuon1->track()->pt()<4.0) continue;
	  if(iMuon2->track()->pt()<4.0) continue;

	  if(!(glbTrackM->quality(reco::TrackBase::highPurity))) continue;
	  if(!(glbTrackP->quality(reco::TrackBase::highPurity))) continue;

	  if( !(iMuon1->isSoftMuon(bestVtx)) ) continue;
	  if( !(iMuon2->isSoftMuon(bestVtx)) ) continue;	 
	  
	  reco::TransientTrack muon1TT((*theB).build(glbTrackP));
	  reco::TransientTrack muon2TT((*theB).build(glbTrackM));

	 // *****  Trajectory states to calculate DCA for the 2 muons *********************
	  FreeTrajectoryState mu1State = muon1TT.impactPointTSCP().theState();
	  FreeTrajectoryState mu2State = muon2TT.impactPointTSCP().theState();

	  if( !muon1TT.impactPointTSCP().isValid() || !muon2TT.impactPointTSCP().isValid() ) continue;

	  // Measure distance between tracks at their closest approach
	  ClosestApproachInRPhi cApp;
	  cApp.calculate(mu1State, mu2State);
	  if( !cApp.status() ) continue;
	  float dca = fabs( cApp.distance() );	  
	  if (dca < 0. || dca > 0.5) continue;
	  //cout<<" closest approach  "<<dca<<endl;

	  // *****  end DCA for the 2 muons *********************
	  
	  //The mass of a muon and the insignificant mass sigma 
	  //to avoid singularities in the covariance matrix.
	  ParticleMass muon_mass = 0.10565837; //pdg mass
	  ParticleMass psi_mass = 3.096916;
	  float muon_sigma = muon_mass*1.e-6;
	  //float psi_sigma = psi_mass*1.e-6;
	  
	  //Creating a KinematicParticleFactory
	  KinematicParticleFactoryFromTransientTrack pFactory;
	  VirtualKinematicParticleFactory vFactory;
		  
	  //initial chi2 and ndf before kinematic fits.
	  float chi = 0.;
	  float ndf = 0.;
	  vector<RefCountedKinematicParticle> muonParticles;
	  muonParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
	  muonParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
	 	  
	  KinematicParticleVertexFitter fitter;   	  
	  RefCountedKinematicTree psiVertexFitTree;
	  psiVertexFitTree = fitter.fit(muonParticles); 	  
	  if (!psiVertexFitTree->isValid()) 
	    {
	      //std::cout << "caught an exception in the psi vertex fit" << std::endl;
	      continue; 
	    }
	  
	  psiVertexFitTree->movePointerToTheTop();
	  
	  RefCountedKinematicParticle psi_vFit_noMC = psiVertexFitTree->currentParticle();
	  RefCountedKinematicVertex psi_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();
	  
	  if( psi_vFit_vertex_noMC->chiSquared() < 0 )
	    {
	      //std::cout << "negative chisq from psi fit" << endl;
	      continue;
	    }
	  	  
	  //if(psi_vFit_vertex_noMC->chiSquared()>50.) continue;
	  if(psi_vFit_noMC->currentState().mass()<2.9 || psi_vFit_noMC->currentState().mass()>3.3) continue;

	  double Omb_J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
	  if(Omb_J_Prob_tmp<0.01)continue;

	  //Now that we have a J/psi candidate, we look for pi^+ candidates
	  
	  for(View<pat::PackedCandidate>::const_iterator iTrack1 = thePATTrackHandle->begin(); 
		   iTrack1 != thePATTrackHandle->end(); ++iTrack1 ) 
		   {
	  
		   if(iTrack1->charge()==0) continue;
		   if(fabs(iTrack1->pdgId())!=211) continue;
		   if(iTrack1->pt()<2.0) continue;
		   if(!(iTrack1->trackHighPurity())) continue;
		   
		   if ( IsTheSame(*iTrack1,*iMuon1) || IsTheSame(*iTrack1,*iMuon2) ) continue;		    
	       
		   reco::TransientTrack pionTT((*theB).build(iTrack1->pseudoTrack()));

		   //ParticleMass kaon_mass = 0.493677;
		   //float kaon_sigma = kaon_mass*1.e-6;
		   ParticleMass pion_mass = 0.13957018;
		   float pion_sigma = pion_mass*1.e-6;

		   float chi = 0.;
		   float ndf = 0.;
		  
		   // ***************************
		   // Jpsipion invariant mass (before kinematic vertex fit)
		   // ***************************
		   TLorentzVector pion14V, Jpsi4V; 
		   pion14V.SetXYZM(iTrack1->px(),iTrack1->py(),iTrack1->pz(),pion_mass);
		   
		   Jpsi4V.SetXYZM(psi_vFit_noMC->currentState().globalMomentum().x(),psi_vFit_noMC->currentState().globalMomentum().y(),psi_vFit_noMC->currentState().globalMomentum().z(),psi_vFit_noMC->currentState().mass());

		   if ( (pion14V + Jpsi4V).M()<5.4 || (pion14V + Jpsi4V).M()>7.3 ) continue;
		   
		   //Now we are ready to combine!
		   // JPsi mass constraint is applied in the final BC fit,

		   vector<RefCountedKinematicParticle> vFitMCParticles;
		   vFitMCParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
		   vFitMCParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
		   vFitMCParticles.push_back(pFactory.particle(pionTT,pion_mass ,chi,ndf,pion_sigma));
		  		  
		   MultiTrackKinematicConstraint *  j_psi_c = new  TwoTrackMassKinematicConstraint(psi_mass);
		   KinematicConstrainedVertexFitter kcvFitter;
		   RefCountedKinematicTree vertexFitTree = kcvFitter.fit(vFitMCParticles, j_psi_c);
		   if (!vertexFitTree->isValid()) {
		     //std::cout << "caught an exception in the B vertex fit with MC" << std::endl;
		     continue;
		   }
		   vertexFitTree->movePointerToTheTop();
		 
		    RefCountedKinematicParticle bCandMC = vertexFitTree->currentParticle();
		    RefCountedKinematicVertex bDecayVertexMC = vertexFitTree->currentDecayVertex();
		    if (!bDecayVertexMC->vertexIsValid()){
		      // cout << "B MC fit vertex is not valid" << endl;
		      continue;
		    }
		    
		    if ( (bCandMC->currentState().mass() < 5.9) || (bCandMC->currentState().mass() > 6.8) ) {
		      // (debug) cout << "continue from bmass > 6.8 or < 5.9 = " << bCandMC->currentState().mass() << endl;
		      continue;
		    }
		  		    
		    double B_Prob_tmp       = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
		    if(B_Prob_tmp<0.01)continue;
		      		    
		    /*
		    // get children from final B fit

		    vertexFitTree->movePointerToTheFirstChild();
		    RefCountedKinematicParticle mu1CandMC = vertexFitTree->currentParticle();
		    
		    vertexFitTree->movePointerToTheNextChild();
		    RefCountedKinematicParticle mu2CandMC = vertexFitTree->currentParticle();		   

		    vertexFitTree->movePointerToTheNextChild();
		    RefCountedKinematicParticle kCandMC = vertexFitTree->currentParticle();
		  
		    KinematicParameters VCandKP = kCandMC->currentState().kinematicParameters();
		    */

		     // ********************* todos los vertices primarios y escogemos el de mejor pointing angle **************** 
		   reco::Vertex bestVtxIP;

		   //bestVtx = *(primaryVertices_handle->begin());


		   Double_t pVtxIPX_temp = -10000.0;
		   Double_t pVtxIPY_temp = -10000.0;
		   Double_t pVtxIPZ_temp = -10000.0;
		   Double_t pVtxIPXE_temp = -10000.0;
		   Double_t pVtxIPYE_temp = -10000.0;
		   Double_t pVtxIPZE_temp = -10000.0;
		   Double_t pVtxIPXYE_temp = -10000.0;
		   Double_t pVtxIPXZE_temp = -10000.0;
		   Double_t pVtxIPYZE_temp = -10000.0;
		   Double_t pVtxIPCL_temp = -10000.0;	
		   Double_t lip1 = -1000000.0;
		     for(size_t i = 0; i < primaryVertices_handle->size(); ++i) {
		       const Vertex &vtx = (*primaryVertices_handle)[i];
		       
		       Double_t dx1 = (*bDecayVertexMC).position().x() - vtx.x(); 
		       Double_t dy1 = (*bDecayVertexMC).position().y() - vtx.y();
		       Double_t dz1 = (*bDecayVertexMC).position().z() - vtx.z();
		       float cosAlphaXYb1 = ( bCandMC->currentState().globalMomentum().x() * dx1 + bCandMC->currentState().globalMomentum().y()*dy1 + bCandMC->currentState().globalMomentum().z()*dz1  )/( sqrt(dx1*dx1+dy1*dy1+dz1*dz1)* bCandMC->currentState().globalMomentum().mag() );

		       if(cosAlphaXYb1>lip1)
			 {
			   lip1 = cosAlphaXYb1 ;
			   pVtxIPX_temp = vtx.x();
			   pVtxIPY_temp = vtx.y();
			   pVtxIPZ_temp = vtx.z();
			   pVtxIPXE_temp = vtx.covariance(0, 0);
			   pVtxIPYE_temp = vtx.covariance(1, 1);
			   pVtxIPZE_temp = vtx.covariance(2, 2);
			   pVtxIPXYE_temp = vtx.covariance(0, 1);
			   pVtxIPXZE_temp = vtx.covariance(0, 2);
			   pVtxIPYZE_temp = vtx.covariance(1, 2);
			   pVtxIPCL_temp = (TMath::Prob(vtx.chi2(),(int)vtx.ndof()) );

			   bestVtxIP = vtx;
			 
			 }
                 
		     }

		     /*pVtxIPX->push_back( pVtxIPX_temp);
		     pVtxIPY->push_back(  pVtxIPY_temp);	    
		     pVtxIPZ->push_back(  pVtxIPZ_temp);
		     pVtxIPXE->push_back( pVtxIPXE_temp);
		     pVtxIPYE->push_back( pVtxIPYE_temp);	    
		     pVtxIPZE->push_back( pVtxIPZE_temp);
		     pVtxIPXYE->push_back( pVtxIPXYE_temp);
		     pVtxIPXZE->push_back( pVtxIPXZE_temp);	    
		     pVtxIPYZE->push_back( pVtxIPYZE_temp);
		     pVtxIPCL->push_back(  pVtxIPCL_temp);*/

		    //*****************************************
  //Let's begin by looking for two pions for Bc(2S)->Bcpi^+pi^-

		    for(View<pat::PackedCandidate>::const_iterator iTrack2 = iTrack1+1;
			iTrack2 != thePATTrackHandle->end(); ++iTrack2 ) 
		      {
			if(iTrack1==iTrack2) continue;
			if(iTrack2->charge()==0) continue;
			if(fabs(iTrack2->pdgId())!=211) continue;
			if(iTrack2->pt()<0.51) continue;
			if(!(iTrack2->trackHighPurity())) continue;

			if ( IsTheSame(*iTrack2,*iMuon1) || IsTheSame(*iTrack2,*iMuon2) ) continue;		    


			for(View<pat::PackedCandidate>::const_iterator iTrack3 = iTrack2+1;
			    iTrack3 != thePATTrackHandle->end(); ++iTrack3 ) 
			  {
			    if(iTrack3==iTrack2 ||  iTrack3==iTrack1) continue;
			    if(iTrack3->charge()==0) continue;
			    if(fabs(iTrack3->pdgId())!=211) continue;
			    if(iTrack3->pt()<0.51) continue;
			    if(!(iTrack3->trackHighPurity())) continue;

			    //opposite charge 
			    if(iTrack3->charge() == iTrack2->charge()) continue;			

			    if ( IsTheSame(*iTrack3,*iMuon1) || IsTheSame(*iTrack3,*iMuon2) ) continue;		    

			    //Now let's see if these two tracks make a vertex with Bc
			    reco::TransientTrack pion2TT((*theB).build(iTrack2->pseudoTrack()));
			    reco::TransientTrack pion3TT((*theB).build(iTrack3->pseudoTrack()));

			    // ***************************
			    // Bc(2S) invariant mass ()
			    // ***************************

			    //ParticleMass bc_mass = 6.2749;// 6.2749 +/- 0.0008   //http://pdglive.lbl.gov/Particle.action?init=0&node=S091
			    TLorentzVector Bc4V,pion24V,pion34V,Bc2S4V; 
			    
			    pion24V.SetXYZM(iTrack2->px(),iTrack2->py(),iTrack2->pz(),pion_mass);
			    pion34V.SetXYZM(iTrack3->px(),iTrack3->py(),iTrack3->pz(),pion_mass);
			    Bc4V.SetXYZM(bCandMC->currentState().globalMomentum().x(),bCandMC->currentState().globalMomentum().y(),bCandMC->currentState().globalMomentum().z(),bCandMC->currentState().mass());
			    
			    Bc2S4V=pion24V+pion34V+Bc4V;
			    //cout<<"mass Bc*: "<<Bc2S4V.M()<<endl;  // ATlas value Bc2S mass = 6842 +/- 4 +/- 5
			    //if ((Bc2S4V.M() - Bc4V.M() + bc_mass < 5.0) || (Bc2S4V.M() - Bc4V.M() + bc_mass > 6.0)) continue;
			    if ( (Bc2S4V.M() < 5.8) || (Bc2S4V.M()> 7.9) ) continue;
			    
			    // using VirtualKinematicParticleFactory vFactory for Bc
			    float Bc_dof  = bDecayVertexMC->degreesOfFreedom();
			    float Bc_chi2 = bDecayVertexMC->chiSquared();

			    vector<RefCountedKinematicParticle> vFitMCParticlesBcpipi;
			    vFitMCParticlesBcpipi.push_back(vFactory.particle(bCandMC->currentState(),Bc_chi2,Bc_dof,bCandMC));
			    vFitMCParticlesBcpipi.push_back(pFactory.particle(pion2TT,pion_mass,chi,ndf,pion_sigma));
			    vFitMCParticlesBcpipi.push_back(pFactory.particle(pion3TT,pion_mass,chi,ndf,pion_sigma));


			    KinematicParticleVertexFitter kcvFitterBcpipi;
			    RefCountedKinematicTree vertexFitTreeBcpipi = kcvFitterBcpipi.fit(vFitMCParticlesBcpipi);
			    if (!vertexFitTreeBcpipi->isValid()) {
			      //std::cout << "caught an exception in the Bc2S vertex fit" << std::endl;
			      continue;
			    }
			    
			    vertexFitTreeBcpipi->movePointerToTheTop();
			    RefCountedKinematicParticle bcpipiCandMC = vertexFitTreeBcpipi->currentParticle();
			    RefCountedKinematicVertex bcpipiDecayVertexMC = vertexFitTreeBcpipi->currentDecayVertex();
			    if (!bcpipiDecayVertexMC->vertexIsValid()){
			      //std::cout << "B MC fit vertex is not valid" << endl;
			      continue;
			    }

			    //if(bcpipiCandMC->currentState().mass()>7.5) continue;
			    if(bcpipiCandMC->currentState().mass()<6.3 || bcpipiCandMC->currentState().mass()>7.3) continue;			    
			    
			    // ************ fill candidate variables now

			    Bc2s_mass->push_back(bcpipiCandMC->currentState().mass());
			    Bc2s_px->push_back(bcpipiCandMC->currentState().globalMomentum().x());
			    Bc2s_py->push_back(bcpipiCandMC->currentState().globalMomentum().y());
			    Bc2s_pz->push_back(bcpipiCandMC->currentState().globalMomentum().z());
			    //Bc2s_charge->push_back(->currentState().particleCharge());
			    double Bc2s_Prob_tmp = TMath::Prob(bcpipiDecayVertexMC->chiSquared(),(int)bcpipiDecayVertexMC->degreesOfFreedom());
			    Bc2s_Prob->push_back(Bc2s_Prob_tmp);

			    Bc2s_pi2_px_track->push_back(iTrack2->px() );
			    Bc2s_pi2_py_track->push_back(iTrack2->py() );
			    Bc2s_pi2_pz_track->push_back(iTrack2->pz() );
			    Bc2s_pi2_charge1->push_back(iTrack2->charge() );

			    Bc2s_pi3_px_track->push_back(iTrack3->px() );
			    Bc2s_pi3_py_track->push_back(iTrack3->py() );
			    Bc2s_pi3_pz_track->push_back(iTrack3->pz() );
			    Bc2s_pi3_charge1->push_back(iTrack3->charge() );		    
			    
			    B_mass->push_back(bCandMC->currentState().mass());
			    B_px->push_back(bCandMC->currentState().globalMomentum().x());
			    B_py->push_back(bCandMC->currentState().globalMomentum().y());
			    B_pz->push_back(bCandMC->currentState().globalMomentum().z());
			    B_charge->push_back(bCandMC->currentState().particleCharge());
			    
			    B_k_px_track->push_back(iTrack1->px() );
			    B_k_py_track->push_back(iTrack1->py() );
			    B_k_pz_track->push_back(iTrack1->pz() );
			    B_k_charge1->push_back(iTrack1->charge() );
			    
			    //B_k_chi2->push_back(ngenT1);
			    
			    B_J_mass->push_back( psi_vFit_noMC->currentState().mass() );
			    B_J_px->push_back( psi_vFit_noMC->currentState().globalMomentum().x() );
			    B_J_py->push_back( psi_vFit_noMC->currentState().globalMomentum().y() );
			    B_J_pz->push_back( psi_vFit_noMC->currentState().globalMomentum().z() );
			    
			    B_J_px1->push_back(iMuon1->track()->px());
			    B_J_py1->push_back(iMuon1->track()->py());
			    B_J_pz1->push_back(iMuon1->track()->pz());		    
			    B_J_charge1->push_back(iMuon1->charge());
			    
			    B_J_px2->push_back(iMuon2->track()->px());
			    B_J_py2->push_back(iMuon2->track()->py());
			    B_J_pz2->push_back(iMuon2->track()->pz());
			    B_J_charge2->push_back(iMuon2->charge());	   
			    
			    //double B_Prob_tmp       = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
			    //double Omb_J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
			    B_Prob    ->push_back(B_Prob_tmp);
			    B_J_Prob  ->push_back(Omb_J_Prob_tmp);			    
			    
			    B_DecayVtxX ->push_back((*bDecayVertexMC).position().x());    
			    B_DecayVtxY ->push_back((*bDecayVertexMC).position().y());
			    B_DecayVtxZ ->push_back((*bDecayVertexMC).position().z());
			    
			    B_DecayVtxXE ->push_back(bDecayVertexMC->error().cxx());   
			    B_DecayVtxYE ->push_back(bDecayVertexMC->error().cyy());   
			    B_DecayVtxZE ->push_back(bDecayVertexMC->error().czz());
			    B_DecayVtxXYE ->push_back(bDecayVertexMC->error().cyx());
			    B_DecayVtxXZE ->push_back(bDecayVertexMC->error().czx());
			    B_DecayVtxYZE ->push_back(bDecayVertexMC->error().czy());		  
			    
			    B_J_DecayVtxX ->push_back( psi_vFit_vertex_noMC->position().x() );
			    B_J_DecayVtxY ->push_back( psi_vFit_vertex_noMC->position().y() );
			    B_J_DecayVtxZ ->push_back( psi_vFit_vertex_noMC->position().z() );
			    
			    B_J_DecayVtxXE ->push_back( psi_vFit_vertex_noMC->error().cxx() );
			    B_J_DecayVtxYE ->push_back( psi_vFit_vertex_noMC->error().cyy() );
			    B_J_DecayVtxZE ->push_back( psi_vFit_vertex_noMC->error().czz() );
			    B_J_DecayVtxXYE ->push_back( psi_vFit_vertex_noMC->error().cyx() );
			    B_J_DecayVtxXZE ->push_back( psi_vFit_vertex_noMC->error().czx() );
			    B_J_DecayVtxYZE ->push_back( psi_vFit_vertex_noMC->error().czy() );

			    pVtxIPX->push_back( pVtxIPX_temp);
			    pVtxIPY->push_back(  pVtxIPY_temp);	    
			    pVtxIPZ->push_back(  pVtxIPZ_temp);
			    pVtxIPXE->push_back( pVtxIPXE_temp);
			    pVtxIPYE->push_back( pVtxIPYE_temp);	    
			    pVtxIPZE->push_back( pVtxIPZE_temp);
			    pVtxIPXYE->push_back( pVtxIPXYE_temp);
			    pVtxIPXZE->push_back( pVtxIPXZE_temp);	    
			    pVtxIPYZE->push_back( pVtxIPYZE_temp);
			    pVtxIPCL->push_back(  pVtxIPCL_temp);
			    
			    // ********************* muon-trigger-machint**************** 
			    const pat::Muon* muon1 = &(*iMuon1);
			    const pat::Muon* muon2 = &(*iMuon2);
			    
			    int tri_Dim25_tmp = 0, tri_JpsiTk_tmp = 0,  tri_Dim20_tmp = 0;
			    
			    if(muon1->triggerObjectMatchByPath("HLT_Dimuon25_Jpsi_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_Dimuon25_Jpsi_v*")!=nullptr) tri_Dim25_tmp = 1;
			    if(muon1->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Displaced_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Displaced_v*")!=nullptr) tri_JpsiTk_tmp = 1;
			    if(muon1->triggerObjectMatchByPath("HLT_Dimuon20_Jpsi_Barrel_Seagulls_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_Dimuon20_Jpsi_Barrel_Seagulls_v*")!=nullptr) tri_Dim20_tmp = 1;
			    
			    tri_Dim25->push_back( tri_Dim25_tmp );	       
			    tri_JpsiTk->push_back( tri_JpsiTk_tmp );
			    tri_Dim20->push_back( tri_Dim20_tmp );			   			   
			    
			    /////////////////////////////////////////////////////
			    
			    nB++;	       
			    muonParticles.clear();
			    vFitMCParticles.clear();
			    vFitMCParticlesBcpipi.clear();
			  }//new pion
		      }//new pion
		   }
	}
    }
 
   
  if (nB > 0 ) 
    {

      //std::cout << "filling tree" << endl;
      tree_->Fill();
    }

   nB = 0; nMu = 0;
   trigger = 0;
    
   Bc2s_mass->clear();    Bc2s_px->clear();    Bc2s_py->clear();    Bc2s_pz->clear();

   Bc2s_pi2_px_track->clear(); Bc2s_pi2_py_track->clear(); Bc2s_pi2_pz_track->clear(); 
   Bc2s_pi3_px_track->clear(); Bc2s_pi3_py_track->clear(); Bc2s_pi3_pz_track->clear();
   Bc2s_pi2_charge1->clear(); Bc2s_pi3_charge1->clear(); Bc2s_Prob->clear();

   B_charge->clear();

   B_mass->clear();    B_px->clear();    B_py->clear();    B_pz->clear(); 
   B_k_charge1->clear(); B_k_px_track->clear(); B_k_py_track->clear(); B_k_pz_track->clear();

   B_J_mass->clear();  B_J_px->clear();  B_J_py->clear();  B_J_pz->clear();

   B_J_px1->clear();  B_J_py1->clear();  B_J_pz1->clear(), B_J_charge1->clear();
   B_J_px2->clear();  B_J_py2->clear();  B_J_pz2->clear(), B_J_charge2->clear();

   B_Prob->clear(); B_J_Prob->clear();

   B_DecayVtxX->clear();     B_DecayVtxY->clear();     B_DecayVtxZ->clear();
   B_DecayVtxXE->clear();    B_DecayVtxYE->clear();    B_DecayVtxZE->clear();
   B_DecayVtxXYE->clear();   B_DecayVtxXZE->clear();   B_DecayVtxYZE->clear();

   B_J_DecayVtxX->clear();   B_J_DecayVtxY->clear();   B_J_DecayVtxZ->clear();
   B_J_DecayVtxXE->clear();  B_J_DecayVtxYE->clear();  B_J_DecayVtxZE->clear();
   B_J_DecayVtxXYE->clear(); B_J_DecayVtxXZE->clear(); B_J_DecayVtxYZE->clear();

   nVtx = 0;
   priVtxX = 0;     priVtxY = 0;     priVtxZ = 0; 
   priVtxXE = 0;    priVtxYE = 0;    priVtxZE = 0; priVtxCL = 0;
   priVtxXYE = 0;   priVtxXZE = 0;   priVtxYZE = 0;

   pVtxIPX->clear();  pVtxIPY->clear();  pVtxIPZ->clear();
   pVtxIPXE->clear();  pVtxIPYE->clear();  pVtxIPZE->clear();  pVtxIPCL->clear();
   pVtxIPXYE->clear();  pVtxIPXZE->clear();  pVtxIPYZE->clear(); 
   
   tri_Dim25->clear(); tri_Dim20->clear(); tri_JpsiTk->clear();

}


bool Bc2SPAT::IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}


// ------------ method called once each job just before starting event loop  ------------


void
Bc2SPAT::beginJob()
{

  std::cout << "Beginning analyzer job with value of isMC= " << isMC_ << std::endl;

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple","Bs->J/psi f0 ntuple");

  tree_->Branch("nB",&nB,"nB/i");
  tree_->Branch("nMu",&nMu,"nMu/i");

  tree_->Branch("Bc2s_mass", &Bc2s_mass);
  tree_->Branch("Bc2s_px", &Bc2s_px);
  tree_->Branch("Bc2s_py", &Bc2s_py);
  tree_->Branch("Bc2s_pz", &Bc2s_pz);
  tree_->Branch("Bc2s_Prob", &Bc2s_Prob);

  tree_->Branch("Bc2s_pi2_px_track", &Bc2s_pi2_px_track);
  tree_->Branch("Bc2s_pi2_py_track", &Bc2s_pi2_py_track);
  tree_->Branch("Bc2s_pi2_pz_track", &Bc2s_pi2_pz_track);
  tree_->Branch("Bc2s_pi2_charge1", &Bc2s_pi2_charge1);

  tree_->Branch("Bc2s_pi3_px_track", &Bc2s_pi3_px_track);
  tree_->Branch("Bc2s_pi3_py_track", &Bc2s_pi3_py_track);
  tree_->Branch("Bc2s_pi3_pz_track", &Bc2s_pi3_pz_track);
  tree_->Branch("Bc2s_pi3_charge1", &Bc2s_pi3_charge1);

  tree_->Branch("B_charge", &B_charge);
  tree_->Branch("B_mass", &B_mass);
  tree_->Branch("B_px", &B_px);
  tree_->Branch("B_py", &B_py);
  tree_->Branch("B_pz", &B_pz);

  tree_->Branch("B_k_charge1", &B_k_charge1);
  tree_->Branch("B_k_px_track", &B_k_px_track);
  tree_->Branch("B_k_py_track", &B_k_py_track);
  tree_->Branch("B_k_pz_track", &B_k_pz_track);

  tree_->Branch("B_J_mass", &B_J_mass);
  tree_->Branch("B_J_px", &B_J_px);
  tree_->Branch("B_J_py", &B_J_py);
  tree_->Branch("B_J_pz", &B_J_pz);

  tree_->Branch("B_J_px1", &B_J_px1);
  tree_->Branch("B_J_py1", &B_J_py1);
  tree_->Branch("B_J_pz1", &B_J_pz1);
  tree_->Branch("B_J_charge1", &B_J_charge1);

  tree_->Branch("B_J_px2", &B_J_px2);
  tree_->Branch("B_J_py2", &B_J_py2);
  tree_->Branch("B_J_pz2", &B_J_pz2);
  tree_->Branch("B_J_charge2", &B_J_charge2);

  tree_->Branch("B_Prob",    &B_Prob);
  tree_->Branch("B_J_Prob",  &B_J_Prob);
       
  tree_->Branch("B_DecayVtxX",     &B_DecayVtxX);
  tree_->Branch("B_DecayVtxY",     &B_DecayVtxY);
  tree_->Branch("B_DecayVtxZ",     &B_DecayVtxZ);
  tree_->Branch("B_DecayVtxXE",    &B_DecayVtxXE);
  tree_->Branch("B_DecayVtxYE",    &B_DecayVtxYE);
  tree_->Branch("B_DecayVtxZE",    &B_DecayVtxZE);
  tree_->Branch("B_DecayVtxXYE",    &B_DecayVtxXYE);
  tree_->Branch("B_DecayVtxXZE",    &B_DecayVtxXZE);
  tree_->Branch("B_DecayVtxYZE",    &B_DecayVtxYZE);
 
  tree_->Branch("B_J_DecayVtxX",   &B_J_DecayVtxX);
  tree_->Branch("B_J_DecayVtxY",   &B_J_DecayVtxY);
  tree_->Branch("B_J_DecayVtxZ",   &B_J_DecayVtxZ);
  tree_->Branch("B_J_DecayVtxXE",  &B_J_DecayVtxXE);
  tree_->Branch("B_J_DecayVtxYE",  &B_J_DecayVtxYE);
  tree_->Branch("B_J_DecayVtxZE",  &B_J_DecayVtxZE);
  tree_->Branch("B_J_DecayVtxXYE",  &B_J_DecayVtxXYE);
  tree_->Branch("B_J_DecayVtxXZE",  &B_J_DecayVtxXZE);
  tree_->Branch("B_J_DecayVtxYZE",  &B_J_DecayVtxYZE);

  tree_->Branch("priVtxX",&priVtxX, "priVtxX/f");
  tree_->Branch("priVtxY",&priVtxY, "priVtxY/f");
  tree_->Branch("priVtxZ",&priVtxZ, "priVtxZ/f");
  tree_->Branch("priVtxXE",&priVtxXE, "priVtxXE/f");
  tree_->Branch("priVtxYE",&priVtxYE, "priVtxYE/f");
  tree_->Branch("priVtxZE",&priVtxZE, "priVtxZE/f");
  tree_->Branch("priVtxXYE",&priVtxXYE, "priVtxXYE/f");
  tree_->Branch("priVtxXZE",&priVtxXZE, "priVtxXZE/f");
  tree_->Branch("priVtxYZE",&priVtxYZE, "priVtxYZE/f");
  tree_->Branch("priVtxCL",&priVtxCL, "priVtxCL/f");

  tree_->Branch("pVtxIPX",     &pVtxIPX);
  tree_->Branch("pVtxIPY",     &pVtxIPY);
  tree_->Branch("pVtxIPZ",     &pVtxIPZ);
  tree_->Branch("pVtxIPXE",     &pVtxIPXE);
  tree_->Branch("pVtxIPYE",     &pVtxIPYE);
  tree_->Branch("pVtxIPZE",     &pVtxIPZE);
  tree_->Branch("pVtxIPXYE",     &pVtxIPXYE);
  tree_->Branch("pVtxIPXZE",     &pVtxIPXZE);
  tree_->Branch("pVtxIPYZE",     &pVtxIPYZE);
  tree_->Branch("pVtxIPCL",     &pVtxIPCL);

 
  tree_->Branch("nVtx",       &nVtx);
  tree_->Branch("run",        &run,       "run/I");
  tree_->Branch("event",        &event,     "event/I");
  tree_->Branch("lumiblock",&lumiblock,"lumiblock/I");

  tree_->Branch("trigger",            &trigger,            "trigger/i");
      
  // *************************

  tree_->Branch("tri_Dim25",&tri_Dim25);
  tree_->Branch("tri_Dim20",&tri_Dim20);
  tree_->Branch("tri_JpsiTk",&tri_JpsiTk);

}


// ------------ method called once each job just after ending the event loop  ------------
void Bc2SPAT::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(Bc2SPAT);
