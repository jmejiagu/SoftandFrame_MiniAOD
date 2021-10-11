// -*- C++ -*-
//
// Package:    JPsiphi
// Class:      JPsiphi
// 
//=================================================
// original author:  Jhovanny Andres Mejia        |
//         created:  November 2020                |
//         <jhovanny.andres.mejia.guisao@cern.ch> | 
//=================================================

// system include files
#include <memory>

// user include files
#include "myAnalyzers/JPsiKsPAT/src/JPsiphi.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/FWLite/interface/EventBase.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"
#include "RecoVertex/VertexTools/interface/VertexDistance.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/CLHEP/interface/Migration.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

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
JPsiphi::JPsiphi(const edm::ParameterSet& iConfig)
  :
  dimuon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
  trakCollection_label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  genCands_(consumes<reco::GenParticleCollection>(iConfig.getParameter < edm::InputTag > ("GenParticles"))), 
  packedGenToken_(consumes<pat::PackedGenParticleCollection>(iConfig.getParameter <edm::InputTag> ("packedGenParticles"))), 
  primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  BSLabel_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bslabel"))),

  OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),

  tree_(0), 

  mumC2(0), mumAngT(0), mumNHits(0), mumNPHits(0),
  mupC2(0), mupAngT(0), mupNHits(0), mupNPHits(0),
  mumdxy(0), mupdxy(0), mumdz(0), mupdz(0),
  muon_dca(0),

  tri_Dim25(0), tri_JpsiTk(0), tri_JpsiTkTk(0),

  mu1soft(0), mu2soft(0), mu1tight(0), mu2tight(0), 
  mu1PF(0), mu2PF(0), mu1loose(0), mu2loose(0),

  // *******************************************************
  
  nB(0), nMu(0),
  B_mass(0), B_px(0), B_py(0), B_pz(0),

  B_phi_mass(0), 
  B_phi_px1(0), B_phi_py1(0), B_phi_pz1(0), 
  B_phi_px2(0), B_phi_py2(0), B_phi_pz2(0),

  B_phi_px1_track(0), B_phi_py1_track(0), B_phi_pz1_track(0), 
  B_phi_px2_track(0), B_phi_py2_track(0), B_phi_pz2_track(0),

  B_phi_charge1(0), B_phi_charge2(0),
  k1dxy(0), k2dxy(0), k1dz(0), k2dz(0),
  k1dxy_e(0), k2dxy_e(0), k1dz_e(0), k2dz_e(0),
  k1InnerHits(0), k2InnerHits(0),

  B_J_mass(0), B_J_px(0), B_J_py(0), B_J_pz(0),
  //B_J_pt1(0),
  B_J_px1(0), B_J_py1(0), B_J_pz1(0), 
  //B_J_pt2(0), 
  B_J_px2(0), B_J_py2(0), B_J_pz2(0), 
  B_J_charge1(0), B_J_charge2(0),

  nVtx(0),
  priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXE(0), priVtxYE(0), priVtxZE(0), priVtxCL(0),
  priVtxXYE(0), priVtxXZE(0), priVtxYZE(0),
  
  // ************************ ****************************************************

  B_chi2(0), B_J_chi2(0), 
  B_Prob(0), B_J_Prob(0), 

  B_DecayVtxX(0),     B_DecayVtxY(0),     B_DecayVtxZ(0),
  B_DecayVtxXE(0),    B_DecayVtxYE(0),    B_DecayVtxZE(0),
  B_DecayVtxXYE(0),   B_DecayVtxXZE(0),   B_DecayVtxYZE(0),

  run(0), event(0),
  lumiblock(0)

{
   //now do what ever initialization is needed
}


JPsiphi::~JPsiphi()
{

}


// ------------ method called to for each event  ------------
void JPsiphi::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  edm::Handle<reco::GenParticleCollection> pruned;
  //edm::Handle<pat::PackedGenParticle> pruned; 
  iEvent.getByToken(genCands_, pruned);
  
  edm::Handle<pat::PackedGenParticleCollection> packed;
  iEvent.getByToken(packedGenToken_,packed);

  //*********************************
  // Get gen level information
  //*********************************  

  gen_b_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_phi_p4.SetPtEtaPhiM(0.,0.,0.,0.); 
  gen_kaon1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_kaon2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muon1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muon2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_b_vtx.SetXYZ(0.,0.,0.);
  gen_jpsi_vtx.SetXYZ(0.,0.,0.);
  gen_b_ct = -9999.;
  
  if ( (isMC_ || OnlyGen_) && pruned.isValid() ) {
    int foundit = 0;
    for (size_t i=0; i<pruned->size(); i++) {
      foundit = 0;
      const reco::Candidate *dau = &(*pruned)[i];
      if ( (abs(dau->pdgId()) == 531) ) { //&& (dau->status() == 2) ) {
	foundit++;
	gen_b_p4.SetPtEtaPhiM(dau->pt(),dau->eta(),dau->phi(),dau->mass());
	gen_b_vtx.SetXYZ(dau->vx(),dau->vy(),dau->vz());
	//int nkaon=0;
	for (size_t k=0; k<dau->numberOfDaughters(); k++) {
	  const reco::Candidate *gdau = dau->daughter(k);
	  if (gdau->pdgId()==443 ) { //&& gdau->status()==2) {
	    foundit++;
	    gen_jpsi_vtx.SetXYZ(gdau->vx(),gdau->vy(),gdau->vz());
	    gen_b_ct = GetLifetime(gen_b_p4,gen_b_vtx,gen_jpsi_vtx);
	    int nm=0;
	    for (size_t l=0; l<gdau->numberOfDaughters(); l++) {
	      const reco::Candidate *mm = gdau->daughter(l);
	      if (mm->pdgId()==13) { foundit++;
		if (mm->status()!=1) {
		  for (size_t m=0; m<mm->numberOfDaughters(); m++) {
		    const reco::Candidate *mu = mm->daughter(m);
		    if (mu->pdgId()==13 ) { //&& mu->status()==1) {
		      nm++;
		      gen_muon1_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
		      break;
		    }
		  }
		} else {
		  gen_muon1_p4.SetPtEtaPhiM(mm->pt(),mm->eta(),mm->phi(),mm->mass());
		  nm++;
		}
	      }
	      if (mm->pdgId()==-13) { foundit++;
		if (mm->status()!=1) {
		  for (size_t m=0; m<mm->numberOfDaughters(); m++) {
		    const reco::Candidate *mu = mm->daughter(m);
		    if (mu->pdgId()==-13 ) { //&& mu->status()==1) {
		      nm++;
		      gen_muon2_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
		      break;
		    }
		  }
		} else {
		  gen_muon2_p4.SetPtEtaPhiM(mm->pt(),mm->eta(),mm->phi(),mm->mass());
		  nm++;
		}
	      }
	    }
	    if (nm==2) gen_jpsi_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
	    else foundit-=nm;
	  }
	  
	  for (size_t lk=0; lk<packed->size(); lk++) {
	     const reco::Candidate * dauInPrunedColl = (*packed)[lk].mother(0);
	     int stable_id = (*packed)[lk].pdgId();
	     if (dauInPrunedColl != nullptr && isAncestor(gdau,dauInPrunedColl)) {
	       if(stable_id == 321) {foundit++;
		 gen_kaon1_p4.SetPtEtaPhiM((*packed)[lk].pt(),(*packed)[lk].eta(),(*packed)[lk].phi(),(*packed)[lk].mass());
	       }
	       if(stable_id == -321){ foundit++;
		 gen_kaon2_p4.SetPtEtaPhiM((*packed)[lk].pt(),(*packed)[lk].eta(),(*packed)[lk].phi(),(*packed)[lk].mass());
	       }
	     }
	   }
	 
	} // for (size_t k
      }   // if (abs(dau->pdgId())==531 )
      if (foundit>=6) break;
    } // for i
    if (foundit!=6) {
      gen_b_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      gen_phi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      gen_b_vtx.SetXYZ(0.,0.,0.);
      gen_jpsi_vtx.SetXYZ(0.,0.,0.);
      gen_b_ct = -9999.;
      std::cout << "Does not found the given decay " << run << "," << event << " foundit=" << foundit << std::endl; // sanity check
    }
  }
  
  //*********************************
  //Now we get the primary vertex 
  //*********************************

  reco::Vertex bestVtx;
  edm::Handle<std::vector<reco::Vertex> > recVtxs;
  iEvent.getByToken(primaryVertices_Label, recVtxs);

  bestVtx = *(recVtxs->begin());
  
  priVtxX = bestVtx.x();
  priVtxY = bestVtx.y();
  priVtxZ = bestVtx.z();
  priVtxXE = bestVtx.covariance(0, 0);
  priVtxYE = bestVtx.covariance(1, 1);
  priVtxZE = bestVtx.covariance(2, 2);
  priVtxXYE = bestVtx.covariance(0, 1);
  priVtxXZE = bestVtx.covariance(0, 2);
  priVtxYZE = bestVtx.covariance(1, 2);
  priVtxCL = ChiSquaredProbability((double)(bestVtx.chi2()),(double)(bestVtx.ndof())); 
  
  nVtx = recVtxs->size();  
 
  lumiblock = iEvent.id().luminosityBlock();
  run = iEvent.id().run();
  event = iEvent.id().event();

  unsigned int nMu_tmp = thePATMuonHandle->size();
  
  for(View<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) 
    {
      for(View<pat::Muon>::const_iterator iMuon2 = iMuon1+1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) 
	{
	  if(iMuon1==iMuon2) continue;
	  
	  //opposite charge 
	  if( (iMuon1->charge())*(iMuon2->charge()) == 1) continue;

	  const pat::Muon *patMuonP = 0;
	  const pat::Muon *patMuonM = 0;
	  TrackRef glbTrackP;	  
	  TrackRef glbTrackM;	  
	  
	  if(iMuon1->charge() == 1){ patMuonP = &(*iMuon1); glbTrackP = iMuon1->track();}
	  if(iMuon1->charge() == -1){patMuonM = &(*iMuon1); glbTrackM = iMuon1->track();}
	  
	  if(iMuon2->charge() == 1) {patMuonP = &(*iMuon2); glbTrackP = iMuon2->track();}
	  if(iMuon2->charge() == -1){patMuonM = &(*iMuon2); glbTrackM = iMuon2->track();}
	  
	  if( glbTrackP.isNull() || glbTrackM.isNull() ) 
	    {
	      //std::cout << "continue due to no track ref" << endl;
	      continue;
	    }

	  if(iMuon1->track()->pt()<4.0) continue;
	  if(iMuon2->track()->pt()<4.0) continue;
	  //if(fabs(iMuon1->eta())>2.2 || fabs(iMuon2->eta())>2.2) continue;

	  if(!(glbTrackM->quality(reco::TrackBase::highPurity))) continue;
	  if(!(glbTrackP->quality(reco::TrackBase::highPurity))) continue;
        
	  //Let's check the vertex and mass
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
	  //if (dca < 0. || dca > 0.5) continue;
	  //cout<<" closest approach  "<<dca<<endl;

	  
	  //The mass of a muon and the insignificant mass sigma 
	  //to avoid singularities in the covariance matrix.
	  ParticleMass muon_mass = 0.10565837; //pdg mass
	  ParticleMass psi_mass = 3.096916;
	  float muon_sigma = muon_mass*1.e-6;
	  //float psi_sigma = psi_mass*1.e-6;

	  //Creating a KinematicParticleFactory
	  KinematicParticleFactoryFromTransientTrack pFactory;
	  
	  //initial chi2 and ndf before kinematic fits.
	  float chi = 0.;
	  float ndf = 0.;
	  vector<RefCountedKinematicParticle> muonParticles;
	  try {
	    muonParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
	    muonParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
	  }
	  catch(...) { 
	    std::cout<<" Exception caught ... continuing 1 "<<std::endl; 
	    continue;
	  }

	  KinematicParticleVertexFitter fitter;   

	  RefCountedKinematicTree psiVertexFitTree;
	  try {
	    psiVertexFitTree = fitter.fit(muonParticles); 
	  }
	  catch (...) { 
	    std::cout<<" Exception caught ... continuing 2 "<<std::endl; 
	    continue;
	  }

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

	   //some loose cuts go here

	   if(psi_vFit_vertex_noMC->chiSquared()>26.) continue;
	   if(psi_vFit_noMC->currentState().mass()<2.9 || psi_vFit_noMC->currentState().mass()>3.3) continue;

	   double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
	   if(J_Prob_tmp<0.01)
	     {
	       continue;
	     }

	   //Now that we have a J/psi candidate, we look for phi candidates

	   for(View<pat::PackedCandidate>::const_iterator iTrack1 = thePATTrackHandle->begin();
	   iTrack1 != thePATTrackHandle->end(); ++iTrack1 )
	     {
	       //quality cuts track1
               if(iTrack1->charge()==0) continue;
	       if(fabs(iTrack1->pdgId())!=211) continue;
	       if(iTrack1->pt()<0.95) continue;
	       if(!(iTrack1->trackHighPurity())) continue;
	       if(iTrack1->numberOfPixelHits()<1)continue;
	       if(iTrack1->numberOfHits()<5)continue;

	       for(View<pat::PackedCandidate>::const_iterator iTrack2 = iTrack1+1;
	       iTrack2 != thePATTrackHandle->end(); ++iTrack2 ) 
		 {
		   //quality cuts track2
		   if(iTrack1==iTrack2) continue;
		   if(iTrack2->charge()==0) continue;
		   if(fabs(iTrack2->pdgId())!=211) continue;
		   if(iTrack2->pt()<0.95) continue;
		   if(!(iTrack2->trackHighPurity())) continue;
		   if(iTrack2->numberOfPixelHits()<1)continue;
		   if(iTrack2->numberOfHits()<5)continue;

		   if(iTrack1->charge() == iTrack2->charge()) continue;

		   //Now let's checks if our muons do not use the same tracks as we are using now
		   if ( IsTheSame(*iTrack1,*iMuon1) || IsTheSame(*iTrack1,*iMuon2) ) continue;
		   if ( IsTheSame(*iTrack2,*iMuon1) || IsTheSame(*iTrack2,*iMuon2) ) continue;

		   //Now let's see if these two tracks make a vertex
		   reco::TransientTrack pion1TT((*theB).build(iTrack1->pseudoTrack()));
		   reco::TransientTrack pion2TT((*theB).build(iTrack2->pseudoTrack()));

		   ParticleMass kaon_mass = 0.493677;
		   float kaon_sigma = kaon_mass*1.e-6;

		   // ***************************
		   // pipi invariant mass (before kinematic vertex fit)
		   // ***************************
		   TLorentzVector pion14V,pion24V,pipi4V, Jpsi4V; 
		   pion14V.SetXYZM(iTrack1->px(),iTrack1->py(),iTrack1->pz(),kaon_mass);
		   pion24V.SetXYZM(iTrack2->px(),iTrack2->py(),iTrack2->pz(),kaon_mass);

		   pipi4V=pion14V+pion24V;
		   if(pipi4V.M()<0.970 || pipi4V.M()>1.070) continue;
	       
		   //initial chi2 and ndf before kinematic fits.
		   float chi = 0.;
		   float ndf = 0.;
	
		   // ***************************
		   // Bs invariant mass (before kinematic vertex fit)
		   // ***************************

		   Jpsi4V.SetXYZM(psi_vFit_noMC->currentState().globalMomentum().x(),psi_vFit_noMC->currentState().globalMomentum().y(),psi_vFit_noMC->currentState().globalMomentum().z(),psi_vFit_noMC->currentState().mass());

		   if ( (pipi4V + Jpsi4V).M()<4.4 || (pipi4V + Jpsi4V).M()>6.3 ) continue;

		   vector<RefCountedKinematicParticle> vFitMCParticles;
		   vFitMCParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
		   vFitMCParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
		   vFitMCParticles.push_back(pFactory.particle(pion1TT,kaon_mass,chi,ndf,kaon_sigma));
		   vFitMCParticles.push_back(pFactory.particle(pion2TT,kaon_mass,chi,ndf,kaon_sigma));

                   // JPsi mass constraint is applied in the final Bs fit,                                                               
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
		     //std::cout << "B MC fit vertex is not valid" << endl;
		     continue;
		   }

		   if(bCandMC->currentState().mass()<5.0 || bCandMC->currentState().mass()>6.0) continue;
		   
		   if(bDecayVertexMC->chiSquared()<0 || bDecayVertexMC->chiSquared()>50 ) 
		     {
		       //std::cout << " continue from negative chi2 = " << bDecayVertexMC->chiSquared() << endl;
		       continue;
		     }		   

		   double B_Prob_tmp  = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
		   if(B_Prob_tmp<0.01)
		     {
		       continue;
		     }
		   
		   // get children from final B fit

		   vertexFitTree->movePointerToTheFirstChild();
                   RefCountedKinematicParticle mu1CandMC = vertexFitTree->currentParticle();

                   vertexFitTree->movePointerToTheNextChild();
                   RefCountedKinematicParticle mu2CandMC = vertexFitTree->currentParticle();

                   vertexFitTree->movePointerToTheNextChild();
                   RefCountedKinematicParticle T1CandMC = vertexFitTree->currentParticle();

                   vertexFitTree->movePointerToTheNextChild();
                   RefCountedKinematicParticle T2CandMC = vertexFitTree->currentParticle();

		   KinematicParameters psiMu1KP = mu1CandMC->currentState().kinematicParameters();
		   KinematicParameters psiMu2KP = mu2CandMC->currentState().kinematicParameters();		   

		   KinematicParameters phiPi1KP = T1CandMC->currentState().kinematicParameters();
		   KinematicParameters phiPi2KP = T2CandMC->currentState().kinematicParameters();

		   const reco::Muon *recoMuonM = patMuonM;
		   const reco::Muon *recoMuonP = patMuonP;
		  		   		   	       
		   // fill candidate variables now
		   
		   if(nB==0){
		     nMu  = nMu_tmp;
		     // cout<< "*Number of Muons : " << nMu_tmp << endl;
		   } // end nB==0

		   B_mass->push_back(bCandMC->currentState().mass());
		   B_px->push_back(bCandMC->currentState().globalMomentum().x());
		   B_py->push_back(bCandMC->currentState().globalMomentum().y());
		   B_pz->push_back(bCandMC->currentState().globalMomentum().z());

		   B_phi_mass->push_back( pipi4V.M() );

		   B_J_mass->push_back( psi_vFit_noMC->currentState().mass() );
		   B_J_px->push_back( psi_vFit_noMC->currentState().globalMomentum().x() );
		   B_J_py->push_back( psi_vFit_noMC->currentState().globalMomentum().y() );
		   B_J_pz->push_back( psi_vFit_noMC->currentState().globalMomentum().z() );

	           // You can get the momentum components (for muons and kaon) from the final B childrens or of the original Tracks. Here, a example for the kaons:
		   B_phi_px1->push_back(phiPi1KP.momentum().x());
		   B_phi_py1->push_back(phiPi1KP.momentum().y());
		   B_phi_pz1->push_back(phiPi1KP.momentum().z());
		   B_phi_px1_track->push_back(iTrack1->px());
		   B_phi_py1_track->push_back(iTrack1->py());
		   B_phi_pz1_track->push_back(iTrack1->pz());
		   B_phi_charge1->push_back(T1CandMC->currentState().particleCharge());

		   B_phi_px2->push_back(phiPi2KP.momentum().x());
		   B_phi_py2->push_back(phiPi2KP.momentum().y());
		   B_phi_pz2->push_back(phiPi2KP.momentum().z());
		   B_phi_px2_track->push_back(iTrack2->px());
		   B_phi_py2_track->push_back(iTrack2->py());
		   B_phi_pz2_track->push_back(iTrack2->pz());
		   B_phi_charge2->push_back(T2CandMC->currentState().particleCharge());

		   B_J_px1->push_back(psiMu1KP.momentum().x());
		   B_J_py1->push_back(psiMu1KP.momentum().y());
		   B_J_pz1->push_back(psiMu1KP.momentum().z());
		   B_J_charge1->push_back(mu1CandMC->currentState().particleCharge());

		   B_J_px2->push_back(psiMu2KP.momentum().x());
		   B_J_py2->push_back(psiMu2KP.momentum().y());
		   B_J_pz2->push_back(psiMu2KP.momentum().z());
		   B_J_charge2->push_back(mu2CandMC->currentState().particleCharge());

		   B_J_chi2->push_back(psi_vFit_vertex_noMC->chiSquared());
		   B_chi2->push_back(bDecayVertexMC->chiSquared());
             
		   //double B_Prob_tmp       = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
		   //double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
		   B_Prob    ->push_back(B_Prob_tmp);
		   B_J_Prob  ->push_back(J_Prob_tmp);

		   B_DecayVtxX ->push_back((*bDecayVertexMC).position().x());    
		   B_DecayVtxY ->push_back((*bDecayVertexMC).position().y());
		   B_DecayVtxZ ->push_back((*bDecayVertexMC).position().z());
		   B_DecayVtxXE ->push_back(bDecayVertexMC->error().cxx());   
		   B_DecayVtxYE ->push_back(bDecayVertexMC->error().cyy());   
		   B_DecayVtxZE ->push_back(bDecayVertexMC->error().czz());
		   B_DecayVtxXYE ->push_back(bDecayVertexMC->error().cyx());
		   B_DecayVtxXZE ->push_back(bDecayVertexMC->error().czx());
		   B_DecayVtxYZE ->push_back(bDecayVertexMC->error().czy());

		   // ********************* muon-trigger-machint ****************

		   const pat::Muon* muon1 = &(*iMuon1);
		   const pat::Muon* muon2 = &(*iMuon2);

		   int tri_Dim25_tmp = 0, tri_JpsiTk_tmp = 0,  tri_JpsiTkTk_tmp = 0;
		   
		   if (muon1->triggerObjectMatchByPath("HLT_Dimuon25_Jpsi_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_Dimuon25_Jpsi_v*")!=nullptr) tri_Dim25_tmp = 1;
		   if (muon1->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Displaced_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Displaced_v*")!=nullptr) tri_JpsiTk_tmp = 1;
		   if (muon1->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrkTrk_Displaced_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrkTrk_Displaced_v*")!=nullptr) tri_JpsiTkTk_tmp = 1;
		   
		   tri_Dim25->push_back( tri_Dim25_tmp );	       
		   tri_JpsiTk->push_back( tri_JpsiTk_tmp );
                   tri_JpsiTkTk->push_back( tri_JpsiTkTk_tmp );

	   // ************ Different muons Id, and other properties  ****************

		   mu1soft->push_back(iMuon1->isSoftMuon(bestVtx) );
		   mu2soft->push_back(iMuon2->isSoftMuon(bestVtx) );
		   mu1tight->push_back(iMuon1->isTightMuon(bestVtx) );
		   mu2tight->push_back(iMuon2->isTightMuon(bestVtx) );
		   mu1PF->push_back(iMuon1->isPFMuon());
		   mu2PF->push_back(iMuon2->isPFMuon());
		   mu1loose->push_back(muon::isLooseMuon(*iMuon1));
		   mu2loose->push_back(muon::isLooseMuon(*iMuon2));

		   mumC2->push_back( glbTrackM->normalizedChi2() );
		   mumAngT->push_back( muon::isGoodMuon(*recoMuonM,muon::TMOneStationTight) ); 
		   mumNHits->push_back( glbTrackM->numberOfValidHits() );
		   mumNPHits->push_back( glbTrackM->hitPattern().numberOfValidPixelHits() );	       
		   mupC2->push_back( glbTrackP->normalizedChi2() );
		   mupAngT->push_back( muon::isGoodMuon(*recoMuonP,muon::TMOneStationTight) ); 
		   mupNHits->push_back( glbTrackP->numberOfValidHits() );
		   mupNPHits->push_back( glbTrackP->hitPattern().numberOfValidPixelHits() );
                   mumdxy->push_back(glbTrackM->dxy(bestVtx.position()) );
		   mupdxy->push_back(glbTrackP->dxy(bestVtx.position()) );
		   mumdz->push_back(glbTrackM->dz(bestVtx.position()) );
		   mupdz->push_back(glbTrackP->dz(bestVtx.position()) );
		   muon_dca->push_back(dca);

		   k1dxy->push_back(iTrack1->dxy());
		   k2dxy->push_back(iTrack2->dxy());
		   k1dz->push_back(iTrack1->dz());
		   k2dz->push_back(iTrack2->dz());

		   k1dxy_e->push_back(iTrack1->dxyError());
		   k2dxy_e->push_back(iTrack2->dxyError());
		   k1dz_e->push_back(iTrack1->dzError());
		   k2dz_e->push_back(iTrack2->dzError());

		   k1InnerHits->push_back(iTrack1->lostInnerHits());
		   k2InnerHits->push_back(iTrack2->lostInnerHits());		   
		  		 		   		   
		   nB++;	       
		   
		   /////////////////////////////////////////////////
		   
		   //pionParticles.clear();
		   muonParticles.clear();
		   vFitMCParticles.clear();

		 }
	     }
	}
    }
 
   //fill the tree and clear the vectors
   if (nB > 0 ) 
     {
       //std::cout << "filling tree" << endl;
       tree_->Fill();
     }
   // *********

   nB = 0; nMu = 0;

   //triggersL = 0; 


   B_mass->clear();    B_px->clear();    B_py->clear();    B_pz->clear();
   B_phi_mass->clear(); 

   B_J_mass->clear();  B_J_px->clear();  B_J_py->clear();  B_J_pz->clear();

   B_phi_px1->clear(); B_phi_py1->clear(); B_phi_pz1->clear(); B_phi_charge1->clear(); 
   B_phi_px2->clear(); B_phi_py2->clear(); B_phi_pz2->clear(); B_phi_charge2->clear(); 

   B_phi_px1_track->clear(); B_phi_py1_track->clear(); B_phi_pz1_track->clear();
   B_phi_px2_track->clear(); B_phi_py2_track->clear(); B_phi_pz2_track->clear();

   B_J_px1->clear();  B_J_py1->clear();  B_J_pz1->clear(), B_J_charge1->clear();
   B_J_px2->clear();  B_J_py2->clear();  B_J_pz2->clear(), B_J_charge2->clear();

   B_chi2->clear(); B_J_chi2->clear();
   B_Prob->clear(); B_J_Prob->clear(); 

   B_DecayVtxX->clear();     B_DecayVtxY->clear();     B_DecayVtxZ->clear();
   B_DecayVtxXE->clear();    B_DecayVtxYE->clear();    B_DecayVtxZE->clear();
   B_DecayVtxXYE->clear();   B_DecayVtxXZE->clear();   B_DecayVtxYZE->clear();

   nVtx = 0; 
   priVtxX = 0;     priVtxY = 0;     priVtxZ = 0; 
   priVtxXE = 0;    priVtxYE = 0;    priVtxZE = 0; priVtxCL = 0;
   priVtxXYE = 0;   priVtxXZE = 0;   priVtxYZE = 0; 

   k1dxy->clear(); k2dxy->clear(); k1dz->clear(); k2dz->clear();
   k1dxy_e->clear(); k2dxy_e->clear(); k1dz_e->clear(); k2dz_e->clear();
   k1InnerHits->clear(); k2InnerHits->clear(); 

   mumC2->clear();
   mumAngT->clear(); mumNHits->clear(); mumNPHits->clear();
   mupC2->clear();
   mupAngT->clear(); mupNHits->clear(); mupNPHits->clear();
   mumdxy->clear(); mupdxy->clear(); mumdz->clear(); mupdz->clear(); muon_dca->clear();

   tri_Dim25->clear(); tri_JpsiTk->clear(); tri_JpsiTkTk->clear(); 

   mu1soft->clear(); mu2soft->clear(); mu1tight->clear(); mu2tight->clear();
   mu1PF->clear(); mu2PF->clear(); mu1loose->clear(); mu2loose->clear(); 
 
}

bool JPsiphi::IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

bool JPsiphi::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
    if (ancestor == particle ) return true;
    for (size_t i=0; i< particle->numberOfMothers(); i++) {
        if (isAncestor(ancestor,particle->mother(i))) return true;
    }
    return false;
}

double JPsiphi::GetLifetime(TLorentzVector b_p4, TVector3 production_vtx, TVector3 decay_vtx) {
   TVector3 pv_dv = decay_vtx - production_vtx;
   TVector3 b_p3  = b_p4.Vect();
   pv_dv.SetZ(0.);
   b_p3.SetZ(0.);
   Double_t lxy   = pv_dv.Dot(b_p3)/b_p3.Mag();
   return lxy*b_p4.M()/b_p3.Mag();
}

// ------------ method called once each job just before starting event loop  ------------

void 
JPsiphi::beginJob()
{
  std::cout << "Beginning analyzer job with value of isMC_ = " << isMC_ << std::endl;

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple","Bs->J/psi phi ntuple");

  tree_->Branch("nB",&nB,"nB/i");
  tree_->Branch("nMu",&nMu,"nMu/i");

  tree_->Branch("B_mass", &B_mass);
  tree_->Branch("B_px", &B_px);
  tree_->Branch("B_py", &B_py);
  tree_->Branch("B_pz", &B_pz);

  tree_->Branch("B_phi_mass", &B_phi_mass);

  tree_->Branch("B_J_mass", &B_J_mass);
  tree_->Branch("B_J_px", &B_J_px);
  tree_->Branch("B_J_py", &B_J_py);
  tree_->Branch("B_J_pz", &B_J_pz);

  tree_->Branch("B_phi_px1", &B_phi_px1);
  tree_->Branch("B_phi_py1", &B_phi_py1);
  tree_->Branch("B_phi_pz1", &B_phi_pz1);
  tree_->Branch("B_phi_px1_track", &B_phi_px1_track);
  tree_->Branch("B_phi_py1_track", &B_phi_py1_track);
  tree_->Branch("B_phi_pz1_track", &B_phi_pz1_track);
  tree_->Branch("B_phi_charge1", &B_phi_charge1); 
 
  tree_->Branch("B_phi_px2", &B_phi_px2);
  tree_->Branch("B_phi_py2", &B_phi_py2);
  tree_->Branch("B_phi_pz2", &B_phi_pz2);
  tree_->Branch("B_phi_px2_track", &B_phi_px2_track);
  tree_->Branch("B_phi_py2_track", &B_phi_py2_track);
  tree_->Branch("B_phi_pz2_track", &B_phi_pz2_track);
  tree_->Branch("B_phi_charge2", &B_phi_charge2);

  tree_->Branch("B_J_px1", &B_J_px1);
  tree_->Branch("B_J_py1", &B_J_py1);
  tree_->Branch("B_J_pz1", &B_J_pz1);
  tree_->Branch("B_J_charge1", &B_J_charge1);

  tree_->Branch("B_J_px2", &B_J_px2);
  tree_->Branch("B_J_py2", &B_J_py2);
  tree_->Branch("B_J_pz2", &B_J_pz2);
  tree_->Branch("B_J_charge2", &B_J_charge2);

  tree_->Branch("B_chi2",    &B_chi2);
  tree_->Branch("B_J_chi2",  &B_J_chi2);

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
  
  tree_->Branch("nVtx",       &nVtx);
  tree_->Branch("run",        &run,       "run/I");
  tree_->Branch("event",        &event,     "event/I");
  tree_->Branch("lumiblock",&lumiblock,"lumiblock/I");
    
  // *************************
  
  tree_->Branch("k1dxy",&k1dxy);
  tree_->Branch("k2dxy",&k2dxy);
  tree_->Branch("k1dz",&k1dz);
  tree_->Branch("k2dz",&k2dz);

  tree_->Branch("k1dxy_e",&k1dxy_e);
  tree_->Branch("k2dxy_e",&k2dxy_e);
  tree_->Branch("k1dz_e",&k1dz_e);
  tree_->Branch("k2dz_e",&k2dz_e);

  tree_->Branch("k1InnerHits",&k1InnerHits);
  tree_->Branch("k2InnerHits",&k2InnerHits);

  tree_->Branch("mumC2",&mumC2);  
  tree_->Branch("mumAngT",&mumAngT);
  tree_->Branch("mumNHits",&mumNHits);
  tree_->Branch("mumNPHits",&mumNPHits);
  tree_->Branch("mupC2",&mupC2);  
  tree_->Branch("mupAngT",&mupAngT);
  tree_->Branch("mupNHits",&mupNHits);
  tree_->Branch("mupNPHits",&mupNPHits);
  tree_->Branch("mumdxy",&mumdxy);
  tree_->Branch("mupdxy",&mupdxy);
  tree_->Branch("mumdz",&mumdz);
  tree_->Branch("mupdz",&mupdz);
  tree_->Branch("muon_dca",&muon_dca);

  tree_->Branch("tri_Dim25",&tri_Dim25);
  tree_->Branch("tri_JpsiTk",&tri_JpsiTk);
  tree_->Branch("tri_JpsiTkTk",&tri_JpsiTkTk); 

  tree_->Branch("mu1soft",&mu1soft);
  tree_->Branch("mu2soft",&mu2soft);
  tree_->Branch("mu1tight",&mu1tight);
  tree_->Branch("mu2tight",&mu2tight);
  tree_->Branch("mu1PF",&mu1PF);
  tree_->Branch("mu2PF",&mu2PF);
  tree_->Branch("mu1loose",&mu1loose);
  tree_->Branch("mu2loose",&mu2loose);

  // gen
  if (isMC_) {
     tree_->Branch("gen_b_p4",     "TLorentzVector",  &gen_b_p4);
     tree_->Branch("gen_phi_p4",   "TLorentzVector",  &gen_phi_p4);
     tree_->Branch("gen_kaon1_p4",  "TLorentzVector",  &gen_kaon1_p4);
     tree_->Branch("gen_kaon2_p4",  "TLorentzVector",  &gen_kaon2_p4);
     tree_->Branch("gen_jpsi_p4",   "TLorentzVector",  &gen_jpsi_p4);
     tree_->Branch("gen_muon1_p4",  "TLorentzVector",  &gen_muon1_p4);
     tree_->Branch("gen_muon2_p4",  "TLorentzVector",  &gen_muon2_p4);
     tree_->Branch("gen_b_vtx",    "TVector3",        &gen_b_vtx);
     tree_->Branch("gen_jpsi_vtx",  "TVector3",        &gen_jpsi_vtx);
     tree_->Branch("gen_b_ct",     &gen_b_ct,        "gen_b_ct/F");
  }

}


// ------------ method called once each job just after ending the event loop  ------------
void JPsiphi::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(JPsiphi);

