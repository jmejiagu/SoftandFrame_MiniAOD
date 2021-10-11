// -*- C++ -*-
//
// Package:    JPsiKs0
// Class:      JPsiKs0
// 
//=================================================
// original author:  Jhovanny Andres Mejia        |
//         created:  November 2020                |
//         <jhovanny.andres.mejia.guisao@cern.ch> | 
//=================================================

// system include files
#include <memory>

// user include files
#include "myAnalyzers/JPsiKsPAT/src/JPsiKs0.h"

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

#include "TFile.h"
#include "TTree.h"

#include <vector>
#include <utility>
#include <string>
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
JPsiKs0::JPsiKs0(const edm::ParameterSet& iConfig)
  :
  dimuon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
  trakCollection_label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  BSLabel_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bslabel"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  v0PtrCollection_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secundaryVerticesPtr"))),	       

  genParticles_ ( iConfig.getUntrackedParameter<std::string>("GenParticles",std::string("genParticles")) ),
  OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
  doMC_ ( iConfig.getUntrackedParameter<bool>("doMC",false) ),
  tree_(0), 

  mumC2(0), mumNHits(0), mumNPHits(0),
  mupC2(0), mupNHits(0), mupNPHits(0),
  mumdxy(0), mupdxy(0), mumdz(0), mupdz(0),
  muon_dca(0),

  tri_Dim25(0), tri_JpsiTk(0), tri_JpsiTkTk(0),

  mu1soft(0), mu2soft(0), mu1tight(0), mu2tight(0), 
  mu1PF(0), mu2PF(0), mu1loose(0), mu2loose(0),

  nVtx(0),
  priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXE(0), priVtxYE(0), priVtxZE(0), priVtxCL(0),
  priVtxXYE(0), priVtxXZE(0), priVtxYZE(0),
 
  // ************************ ****************************************************

  bDecayVtxX(0), bDecayVtxY(0), bDecayVtxZ(0), bDecayVtxXE(0), bDecayVtxYE(0), bDecayVtxZE(0),
  bDecayVtxXYE(0), bDecayVtxXZE(0), bDecayVtxYZE(0),
  
  VDecayVtxX(0), VDecayVtxY(0), VDecayVtxZ(0), VDecayVtxXE(0), VDecayVtxYE(0), VDecayVtxZE(0),
  VDecayVtxXYE(0), VDecayVtxXZE(0), VDecayVtxYZE(0), 

  // *******************************************************
  nB(0), nMu(0),
  B_mass(0), B_px(0), B_py(0), B_pz(0),
  
  B_Ks0_mass(0), B_Ks0_px(0), B_Ks0_py(0), B_Ks0_pz(0),
  B_Ks0_pt1(0), B_Ks0_px1(0), B_Ks0_py1(0), B_Ks0_pz1(0), 
  B_Ks0_pt2(0), B_Ks0_px2(0), B_Ks0_py2(0), B_Ks0_pz2(0), 

  B_Ks0_px1_track(0), B_Ks0_py1_track(0), B_Ks0_pz1_track(0), 
  B_Ks0_px2_track(0), B_Ks0_py2_track(0), B_Ks0_pz2_track(0), 

  pi1dxy(0), pi2dxy(0), pi1dz(0), pi2dz(0),
  pi1dxy_e(0), pi2dxy_e(0), pi1dz_e(0), pi2dz_e(0),
  B_Ks0_charge1(0), B_Ks0_charge2(0),

  B_J_mass(0), B_J_px(0), B_J_py(0), B_J_pz(0),
  B_J_pt1(0), B_J_px1(0), B_J_py1(0), B_J_pz1(0), 
  B_J_pt2(0), B_J_px2(0), B_J_py2(0), B_J_pz2(0), 
  B_J_charge1(0), B_J_charge2(0),

  B_Ks0_chi2(0), B_J_chi2(0), B_chi2(0),
  B_Prob(0), B_J_Prob(0), B_ks0_Prob(0),

  run(0), event(0),
  lumiblock(0)

{
   //now do what ever initialization is needed
}

JPsiKs0::~JPsiKs0()
{

}

// ------------ method called to for each event  ------------
void JPsiKs0::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;
  
  //*********************************
  // Get event content information
  //*********************************  

  edm::Handle<std::vector<reco::VertexCompositePtrCandidate>> theV0PtrHandle;
  iEvent.getByToken(v0PtrCollection_,  theV0PtrHandle);

// Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB; 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); 

  edm::Handle< View<pat::PackedCandidate> > thePATTrackHandle;
  iEvent.getByToken(trakCollection_label,thePATTrackHandle);

  edm::Handle< View<pat::Muon> > thePATMuonHandle;
  iEvent.getByToken(dimuon_Label,thePATMuonHandle);

  //*********************************
  //Now we get the primary vertex 
  //*********************************

  reco::Vertex bestVtx;
  reco::Vertex bestVtxBS;

  // get primary vertex
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

  //*****************************************
  //Let's begin by looking for J/psi->mu+mu-

  unsigned int nMu_tmp = thePATMuonHandle->size();
 
 for(View<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) 
    {
      
      for(View<pat::Muon>::const_iterator iMuon2 = iMuon1+1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) 
	{
  
	  if(iMuon1==iMuon2) continue;
	  
	  //opposite charge 
	  if( (iMuon1->charge())*(iMuon2->charge()) == 1) continue; // <-------------------------------------

	  TrackRef glbTrackP;	  
	  TrackRef glbTrackM;	  
	  
	  if(iMuon1->charge() == 1){glbTrackP = iMuon1->track();}
	  if(iMuon1->charge() == -1){glbTrackM = iMuon1->track();}
	  
	  if(iMuon2->charge() == 1) {glbTrackP = iMuon2->track();}
	  if(iMuon2->charge() == -1){glbTrackM = iMuon2->track();}
	  
	  if( glbTrackP.isNull() || glbTrackM.isNull() ) 
	    {
	      //std::cout << "continue due to no track ref" << endl;
	      continue;
	    }

	  if(iMuon1->track()->pt()<4.0) continue;
	  if(iMuon2->track()->pt()<4.0) continue;

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

	  // *****  end DCA for the 2 muons *********************

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
	  
	  RefCountedKinematicParticle psi_vFit_noMC = psiVertexFitTree->currentParticle();//masa del J/psi
	  RefCountedKinematicVertex psi_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();//vertice del J/psi
	  
	  if( psi_vFit_vertex_noMC->chiSquared() < 0 )
	    {
	      //std::cout << "negative chisq from psi fit" << endl;
	      continue;
	    }

	  double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
	  if(J_Prob_tmp<0.01)
	    {
	      continue;
	    }
	  
	   //some loose cuts go here

	   if(psi_vFit_vertex_noMC->chiSquared()>50.) continue;
	   if(psi_vFit_noMC->currentState().mass()<2.9 || psi_vFit_noMC->currentState().mass()>3.3) continue;

	   //  ***************  

	   if ( theV0PtrHandle->size()>0 && thePATMuonHandle->size()>=2 ) 
	     {
	       
	       for ( vector<VertexCompositePtrCandidate>::const_iterator iVee = theV0PtrHandle->begin();   iVee != theV0PtrHandle->end(); ++iVee )
		 {
		   //get Lam tracks from V0 candidate
		   vector<pat::PackedCandidate> v0daughters;
		     vector<Track> theDaughterTracks;
		     v0daughters.push_back( *(dynamic_cast<const pat::PackedCandidate *>(iVee->daughter(0))) );
		     v0daughters.push_back( *(dynamic_cast<const pat::PackedCandidate *>(iVee->daughter(1))) );
		     		     
		     for(unsigned int j = 0; j < v0daughters.size(); ++j)
		       {
			 theDaughterTracks.push_back(v0daughters[j].pseudoTrack());
		       }
		     			
		     // it does not have sences here. 
		     //if ( IsTheSame(*theDaughterTracks[0],*iMuon1) || IsTheSame(*theDaughterTracks[0],*iMuon2) ) continue;
		     //if ( IsTheSame(*theDaughterTracks[1],*iMuon1) || IsTheSame(*theDaughterTracks[1],*iMuon2) ) continue;
		     
		     //Now let's see if these two tracks make a vertex
		     reco::TransientTrack pion1TT((*theB).build(theDaughterTracks[0]));
		     reco::TransientTrack pion2TT((*theB).build(theDaughterTracks[1]));		     
		     
		     ParticleMass pion_mass = 0.13957018;
		     ParticleMass Ks0_mass = 0.497614;
		     float pion_sigma = pion_mass*1.e-6;
		     float Ks0_sigma = Ks0_mass*1.e-6;
		     
		     //initial chi2 and ndf before kinematic fits.
		     float chi = 0.;
		     float ndf = 0.;
		     vector<RefCountedKinematicParticle> pionParticles;
		     // vector<RefCountedKinematicParticle> muonParticles;
		     try {
		       pionParticles.push_back(pFactory.particle(pion1TT,pion_mass,chi,ndf,pion_sigma));
		       pionParticles.push_back(pFactory.particle(pion2TT,pion_mass,chi,ndf,pion_sigma));
		     }
		     catch(...) {
		       std::cout<<" Exception caught ... continuing 3 "<<std::endl;
		       continue;
		     }
		     
		     RefCountedKinematicTree Ks0VertexFitTree;
		     try{
		       Ks0VertexFitTree = fitter.fit(pionParticles); 
		     }
		     catch(...) {
		       std::cout<<" Exception caught ... continuing 4 "<<std::endl;                   
		       continue;
		     }
		     if (!Ks0VertexFitTree->isValid()) 
		       {
			 //std::cout << "invalid vertex from the Ks0 vertex fit" << std::endl;
			 continue; 
		       }
		     Ks0VertexFitTree->movePointerToTheTop();
		     
		     RefCountedKinematicParticle Ks0_vFit_noMC = Ks0VertexFitTree->currentParticle();
		     RefCountedKinematicVertex Ks0_vFit_vertex_noMC = Ks0VertexFitTree->currentDecayVertex();
		     
		     if( Ks0_vFit_vertex_noMC->chiSquared() < 0 )
		       { 
			 //std::cout << "negative chisq from ks fit" << endl;
			 continue;
		       }
		     
		     //some loose cuts go here
		     
		     if(Ks0_vFit_vertex_noMC->chiSquared()>50) continue;
		     if(Ks0_vFit_noMC->currentState().mass()<0.45 || Ks0_vFit_noMC->currentState().mass()>0.55) continue;
		     
		     Ks0VertexFitTree->movePointerToTheFirstChild();
		     RefCountedKinematicParticle T1CandMC = Ks0VertexFitTree->currentParticle();
		     
		     Ks0VertexFitTree->movePointerToTheNextChild();
		     RefCountedKinematicParticle T2CandMC = Ks0VertexFitTree->currentParticle();
		     
		     //  Ks0  mass constrain
		     // do mass constrained vertex fit
		     // creating the constraint with a small sigma to put in the resulting covariance 
		     // matrix in order to avoid singularities
		     // JPsi mass constraint is applied in the final B fit
		     
		     KinematicParticleFitter csFitterKs;
		     KinematicConstraint * ks_c = new MassKinematicConstraint(Ks0_mass,Ks0_sigma);
		     // add mass constraint to the ks0 fit to do a constrained fit:  
		     
		     Ks0VertexFitTree = csFitterKs.fit(ks_c,Ks0VertexFitTree);
		     if (!Ks0VertexFitTree->isValid()){
		       //std::cout << "caught an exception in the ks mass constraint fit" << std::endl;
		       continue; 
		     }
		     
		     Ks0VertexFitTree->movePointerToTheTop();
		     RefCountedKinematicParticle ks0_vFit_withMC = Ks0VertexFitTree->currentParticle();
		     
		     //Now we are ready to combine!
		     // JPsi mass constraint is applied in the final Bd fit,
		     
		     vector<RefCountedKinematicParticle> vFitMCParticles;
		     vFitMCParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
		     vFitMCParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
		     vFitMCParticles.push_back(ks0_vFit_withMC);
		     
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
		     
		     double B_Prob_tmp       = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
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
		   RefCountedKinematicParticle Ks0CandMC = vertexFitTree->currentParticle();
		   
		   KinematicParameters psiMu1KP = mu1CandMC->currentState().kinematicParameters();
		   KinematicParameters psiMu2KP = mu2CandMC->currentState().kinematicParameters();
		   KinematicParameters psiMupKP;
		   KinematicParameters psiMumKP;
	       
		   if ( mu1CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu1KP;
		   if ( mu1CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu1KP;
		   if ( mu2CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu2KP;
		   if ( mu2CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu2KP;	 

 		   GlobalVector Jp1vec(mu1CandMC->currentState().globalMomentum().x(),
				       mu1CandMC->currentState().globalMomentum().y(),
 				       mu1CandMC->currentState().globalMomentum().z());

 		   GlobalVector Jp2vec(mu2CandMC->currentState().globalMomentum().x(),
				       mu2CandMC->currentState().globalMomentum().y(),
 				       mu2CandMC->currentState().globalMomentum().z());

 		   GlobalVector Ks0p1vec(T1CandMC->currentState().globalMomentum().x(),
				        T1CandMC->currentState().globalMomentum().y(),
 				        T1CandMC->currentState().globalMomentum().z());

 		   GlobalVector Ks0p2vec(T2CandMC->currentState().globalMomentum().x(),
					T2CandMC->currentState().globalMomentum().y(),
					T2CandMC->currentState().globalMomentum().z());

		   KinematicParameters Ks0Pi1KP = T1CandMC->currentState().kinematicParameters();
		   KinematicParameters Ks0Pi2KP = T2CandMC->currentState().kinematicParameters();
		   KinematicParameters Ks0PipKP;
		   KinematicParameters Ks0PimKP;
	       
		   if ( T1CandMC->currentState().particleCharge() > 0 ) Ks0PipKP = Ks0Pi1KP;
		   if ( T1CandMC->currentState().particleCharge() < 0 ) Ks0PimKP = Ks0Pi1KP;
		   if ( T2CandMC->currentState().particleCharge() > 0 ) Ks0PipKP = Ks0Pi2KP;
		   if ( T2CandMC->currentState().particleCharge() < 0 ) Ks0PimKP = Ks0Pi2KP;	 

		   // fill candidate variables now
		   
		   if(nB==0){		    
		     nMu  = nMu_tmp;
		     // cout<< "*Number of Muons : " << nMu_tmp << endl;
		   } // end nB==0		     
		
		   B_mass->push_back(bCandMC->currentState().mass());
		   B_px->push_back(bCandMC->currentState().globalMomentum().x());
		   B_py->push_back(bCandMC->currentState().globalMomentum().y());
		   B_pz->push_back(bCandMC->currentState().globalMomentum().z());
		
		   B_Ks0_mass->push_back( Ks0_vFit_noMC->currentState().mass() );
		   B_Ks0_px->push_back( Ks0_vFit_noMC->currentState().globalMomentum().x() );
		   B_Ks0_py->push_back( Ks0_vFit_noMC->currentState().globalMomentum().y() );
		   B_Ks0_pz->push_back( Ks0_vFit_noMC->currentState().globalMomentum().z() );

		   B_J_mass->push_back( psi_vFit_noMC->currentState().mass() );
		   B_J_px->push_back( psi_vFit_noMC->currentState().globalMomentum().x() );
		   B_J_py->push_back( psi_vFit_noMC->currentState().globalMomentum().y() );
		   B_J_pz->push_back( psi_vFit_noMC->currentState().globalMomentum().z() );

		   B_Ks0_pt1->push_back(Ks0p1vec.perp());
		   B_Ks0_px1->push_back(Ks0Pi1KP.momentum().x());
		   B_Ks0_py1->push_back(Ks0Pi1KP.momentum().y());
		   B_Ks0_pz1->push_back(Ks0Pi1KP.momentum().z());
		   B_Ks0_px1_track->push_back(v0daughters[0].px());
		   B_Ks0_py1_track->push_back(v0daughters[0].py());
		   B_Ks0_pz1_track->push_back(v0daughters[0].pz());
		   B_Ks0_charge1->push_back(T1CandMC->currentState().particleCharge());

		   B_Ks0_pt2->push_back(Ks0p2vec.perp());
		   B_Ks0_px2->push_back(Ks0Pi2KP.momentum().x());
		   B_Ks0_py2->push_back(Ks0Pi2KP.momentum().y());
		   B_Ks0_pz2->push_back(Ks0Pi2KP.momentum().z());
		   B_Ks0_px2_track->push_back(v0daughters[1].px());
		   B_Ks0_py2_track->push_back(v0daughters[1].py());
		   B_Ks0_pz2_track->push_back(v0daughters[1].pz());
		   B_Ks0_charge2->push_back(T2CandMC->currentState().particleCharge());

		   B_J_pt1->push_back(Jp1vec.perp());
		   B_J_px1->push_back(psiMu1KP.momentum().x());
		   B_J_py1->push_back(psiMu1KP.momentum().y());
		   B_J_pz1->push_back(psiMu1KP.momentum().z());
		   B_J_charge1->push_back(mu1CandMC->currentState().particleCharge());

		   B_J_pt2->push_back(Jp2vec.perp());
		   B_J_px2->push_back(psiMu2KP.momentum().x());
		   B_J_py2->push_back(psiMu2KP.momentum().y());
		   B_J_pz2->push_back(psiMu2KP.momentum().z());
		   B_J_charge2->push_back(mu2CandMC->currentState().particleCharge());

		   B_Ks0_chi2->push_back(Ks0_vFit_vertex_noMC->chiSquared());
		   B_J_chi2->push_back(psi_vFit_vertex_noMC->chiSquared());
		   B_chi2->push_back(bDecayVertexMC->chiSquared());

		   //double B_Prob_tmp       = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
		   //double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
		   double ks0_Prob_tmp  = TMath::Prob(Ks0_vFit_vertex_noMC->chiSquared(),(int)Ks0_vFit_vertex_noMC->degreesOfFreedom());
		   B_Prob    ->push_back(B_Prob_tmp);
		   B_J_Prob  ->push_back(J_Prob_tmp);
		   B_ks0_Prob ->push_back(ks0_Prob_tmp);

	   // ************
		   bDecayVtxX->push_back((*bDecayVertexMC).position().x());
		   bDecayVtxY->push_back((*bDecayVertexMC).position().y());
		   bDecayVtxZ->push_back((*bDecayVertexMC).position().z());
		   bDecayVtxXE->push_back(bDecayVertexMC->error().cxx());
		   bDecayVtxYE->push_back(bDecayVertexMC->error().cyy());
		   bDecayVtxZE->push_back(bDecayVertexMC->error().czz());
		   bDecayVtxXYE->push_back(bDecayVertexMC->error().cyx());
		   bDecayVtxXZE->push_back(bDecayVertexMC->error().czx());
		   bDecayVtxYZE->push_back(bDecayVertexMC->error().czy());

		   VDecayVtxX->push_back( Ks0_vFit_vertex_noMC->position().x() );
		   VDecayVtxY->push_back( Ks0_vFit_vertex_noMC->position().y() );
		   VDecayVtxZ->push_back( Ks0_vFit_vertex_noMC->position().z() );
		   VDecayVtxXE->push_back( Ks0_vFit_vertex_noMC->error().cxx() );
		   VDecayVtxYE->push_back( Ks0_vFit_vertex_noMC->error().cyy() );
		   VDecayVtxZE->push_back( Ks0_vFit_vertex_noMC->error().czz() );
		   VDecayVtxXYE->push_back( Ks0_vFit_vertex_noMC->error().cyx() );
		   VDecayVtxXZE->push_back( Ks0_vFit_vertex_noMC->error().czx() );
		   VDecayVtxYZE->push_back( Ks0_vFit_vertex_noMC->error().czy() );

 // ********************* muon-trigger-machint**************** 
		   
		   const pat::Muon* muon1 = &(*iMuon1);
		   const pat::Muon* muon2 = &(*iMuon2);

		   int tri_Dim25_tmp = 0, tri_JpsiTk_tmp = 0,  tri_JpsiTkTk_tmp = 0;
		   
		   if (muon1->triggerObjectMatchByPath("HLT_Dimuon25_Jpsi_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_Dimuon25_Jpsi_v*")!=nullptr) tri_Dim25_tmp = 1;
		   if (muon1->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Displaced_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Displaced_v*")!=nullptr) tri_JpsiTk_tmp = 1;
		   if (muon1->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrkTrk_Displaced_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrkTrk_Displaced_v*")!=nullptr) tri_JpsiTkTk_tmp = 1;
		   
		   tri_Dim25->push_back( tri_Dim25_tmp );	       
		   tri_JpsiTk->push_back( tri_JpsiTk_tmp );
                   tri_JpsiTkTk->push_back( tri_JpsiTkTk_tmp );

 	   // ************
		  
		   mu1soft->push_back(iMuon1->isSoftMuon(bestVtx) );
		   mu2soft->push_back(iMuon2->isSoftMuon(bestVtx) );
		   mu1tight->push_back(iMuon1->isTightMuon(bestVtx) );
		   mu2tight->push_back(iMuon2->isTightMuon(bestVtx) );
		   mu1PF->push_back(iMuon1->isPFMuon());
		   mu2PF->push_back(iMuon2->isPFMuon());
		   mu1loose->push_back(muon::isLooseMuon(*iMuon1));
		   mu2loose->push_back(muon::isLooseMuon(*iMuon2));

		   mumC2->push_back( glbTrackM->normalizedChi2() );
		   mumNHits->push_back( glbTrackM->numberOfValidHits() );
		   mumNPHits->push_back( glbTrackM->hitPattern().numberOfValidPixelHits() );	       
		   mupC2->push_back( glbTrackP->normalizedChi2() );
		   mupNHits->push_back( glbTrackP->numberOfValidHits() );
		   mupNPHits->push_back( glbTrackP->hitPattern().numberOfValidPixelHits() );
                   mumdxy->push_back(glbTrackM->dxy(bestVtx.position()) );
		   mupdxy->push_back(glbTrackP->dxy(bestVtx.position()) );
		   mumdz->push_back(glbTrackM->dz(bestVtx.position()) );
		   mupdz->push_back(glbTrackP->dz(bestVtx.position()) );
		   muon_dca->push_back(dca);

		   pi1dxy->push_back(v0daughters[0].dxy());
		   pi2dxy->push_back(v0daughters[1].dxy());
		   pi1dz->push_back(v0daughters[0].dz());
		   pi2dz->push_back(v0daughters[1].dz());

		   pi1dxy_e->push_back(v0daughters[0].dxyError());
		   pi2dxy_e->push_back(v0daughters[1].dxyError());
		   pi1dz_e->push_back(v0daughters[0].dzError());
		   pi2dz_e->push_back(v0daughters[1].dzError());

		   // try refitting the primary without the tracks in the B reco candidate		   
		  
		  nB++;	       
		   
		   /////////////////////////////////////////////////
		   pionParticles.clear();
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

   B_mass->clear();    B_px->clear();    B_py->clear();    B_pz->clear();
   B_Ks0_mass->clear(); B_Ks0_px->clear(); B_Ks0_py->clear(); B_Ks0_pz->clear();

   B_J_mass->clear();  B_J_px->clear();  B_J_py->clear();  B_J_pz->clear();

   B_Ks0_pt1->clear(); B_Ks0_px1->clear(); B_Ks0_py1->clear(); B_Ks0_pz1->clear(); B_Ks0_charge1->clear(); 
   B_Ks0_pt2->clear(); B_Ks0_px2->clear(); B_Ks0_py2->clear(); B_Ks0_pz2->clear(); B_Ks0_charge2->clear(); 

   B_Ks0_px1_track->clear(); B_Ks0_py1_track->clear(); B_Ks0_pz1_track->clear(); 
   B_Ks0_px2_track->clear(); B_Ks0_py2_track->clear(); B_Ks0_pz2_track->clear(); 

   B_J_pt1->clear();  B_J_px1->clear();  B_J_py1->clear();  B_J_pz1->clear(), B_J_charge1->clear();
   B_J_pt2->clear();  B_J_px2->clear();  B_J_py2->clear();  B_J_pz2->clear(), B_J_charge2->clear();

   B_Ks0_chi2->clear(); B_J_chi2->clear(); B_chi2->clear();
   B_Prob->clear(); B_J_Prob->clear(); B_ks0_Prob->clear();

   // *********

   nVtx = 0;
   priVtxX = 0; priVtxY = 0; priVtxZ = 0; 
   priVtxXE = 0; priVtxYE = 0; priVtxZE = 0; priVtxCL = 0;
   priVtxXYE = 0;   priVtxXZE = 0;   priVtxYZE = 0;

   bDecayVtxX->clear(); bDecayVtxY->clear(); bDecayVtxZ->clear(); 
   bDecayVtxXE->clear(); bDecayVtxYE->clear(); bDecayVtxZE->clear(); 
   bDecayVtxXYE->clear(); bDecayVtxXZE->clear(); bDecayVtxYZE->clear();  

   VDecayVtxX->clear(); VDecayVtxY->clear(); VDecayVtxZ->clear();
   VDecayVtxXE->clear(); VDecayVtxYE->clear(); VDecayVtxZE->clear();
   VDecayVtxXYE->clear(); VDecayVtxXZE->clear(); VDecayVtxYZE->clear();
  
   pi1dxy->clear(); pi2dxy->clear(); pi1dz->clear(); pi2dz->clear();
   pi1dxy_e->clear(); pi2dxy_e->clear(); pi1dz_e->clear(); pi2dz_e->clear();

   mumC2->clear();
   mumNHits->clear(); mumNPHits->clear();
   mupC2->clear();
   mupNHits->clear(); mupNPHits->clear();
   mumdxy->clear(); mupdxy->clear(); mumdz->clear(); mupdz->clear(); muon_dca->clear();

   tri_Dim25->clear(); tri_JpsiTk->clear(); tri_JpsiTkTk->clear();
      
   mu1soft->clear(); mu2soft->clear(); mu1tight->clear(); mu2tight->clear();
   mu1PF->clear(); mu2PF->clear(); mu1loose->clear(); mu2loose->clear(); 

}

bool JPsiKs0::IsTheSame(const reco::Track& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

// ------------ method called once each job just before starting event loop  ------------

void 
JPsiKs0::beginJob()
{

  std::cout << "Beginning analyzer job with value of doMC = " << doMC_ << std::endl;

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple","Bs->J/psi Ks0 ntuple");

  tree_->Branch("nB",&nB,"nB/i");
  tree_->Branch("nMu",&nMu,"nMu/i");

  tree_->Branch("B_mass", &B_mass);
  tree_->Branch("B_px", &B_px);
  tree_->Branch("B_py", &B_py);
  tree_->Branch("B_pz", &B_pz);

  tree_->Branch("B_Ks0_mass", &B_Ks0_mass);
  tree_->Branch("B_Ks0_px", &B_Ks0_px);
  tree_->Branch("B_Ks0_py", &B_Ks0_py);
  tree_->Branch("B_Ks0_pz", &B_Ks0_pz);
 
  tree_->Branch("B_J_mass", &B_J_mass);
  tree_->Branch("B_J_px", &B_J_px);
  tree_->Branch("B_J_py", &B_J_py);
  tree_->Branch("B_J_pz", &B_J_pz);

  tree_->Branch("B_Ks0_pt1", &B_Ks0_pt1);
  tree_->Branch("B_Ks0_px1", &B_Ks0_px1);
  tree_->Branch("B_Ks0_py1", &B_Ks0_py1);
  tree_->Branch("B_Ks0_pz1", &B_Ks0_pz1);
  tree_->Branch("B_Ks0_px1_track", &B_Ks0_px1_track);
  tree_->Branch("B_Ks0_py1_track", &B_Ks0_py1_track);
  tree_->Branch("B_Ks0_pz1_track", &B_Ks0_pz1_track);
  tree_->Branch("B_Ks0_charge1", &B_Ks0_charge1); 
 
  tree_->Branch("B_Ks0_pt2", &B_Ks0_pt2);
  tree_->Branch("B_Ks0_px2", &B_Ks0_px2);
  tree_->Branch("B_Ks0_py2", &B_Ks0_py2);
  tree_->Branch("B_Ks0_pz2", &B_Ks0_pz2);
  tree_->Branch("B_Ks0_px2_track", &B_Ks0_px2_track);
  tree_->Branch("B_Ks0_py2_track", &B_Ks0_py2_track);
  tree_->Branch("B_Ks0_pz2_track", &B_Ks0_pz2_track);
  tree_->Branch("B_Ks0_charge2", &B_Ks0_charge2);

  tree_->Branch("B_J_pt1", &B_J_pt1);
  tree_->Branch("B_J_px1", &B_J_px1);
  tree_->Branch("B_J_py1", &B_J_py1);
  tree_->Branch("B_J_pz1", &B_J_pz1);
  tree_->Branch("B_J_charge1", &B_J_charge1);

  tree_->Branch("B_J_pt2", &B_J_pt2);
  tree_->Branch("B_J_px2", &B_J_px2);
  tree_->Branch("B_J_py2", &B_J_py2);
  tree_->Branch("B_J_pz2", &B_J_pz2);
  tree_->Branch("B_J_charge2", &B_J_charge2);

  tree_->Branch("B_chi2", &B_chi2);
  tree_->Branch("B_Ks0_chi2", &B_Ks0_chi2);
  tree_->Branch("B_J_chi2", &B_J_chi2);

  tree_->Branch("B_Prob",    &B_Prob);
  tree_->Branch("B_ks0_Prob", &B_ks0_Prob);
  tree_->Branch("B_J_Prob",  &B_J_Prob);
       
  // *************************

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

  tree_->Branch("bDecayVtxX",&bDecayVtxX);
  tree_->Branch("bDecayVtxY",&bDecayVtxY);
  tree_->Branch("bDecayVtxZ",&bDecayVtxZ);
  tree_->Branch("bDecayVtxXE",&bDecayVtxXE);
  tree_->Branch("bDecayVtxYE",&bDecayVtxYE);
  tree_->Branch("bDecayVtxZE",&bDecayVtxZE);
  tree_->Branch("bDecayVtxXYE",&bDecayVtxXYE);
  tree_->Branch("bDecayVtxXZE",&bDecayVtxXZE);
  tree_->Branch("bDecayVtxYZE",&bDecayVtxYZE);

  tree_->Branch("VDecayVtxX",&VDecayVtxX);
  tree_->Branch("VDecayVtxY",&VDecayVtxY);
  tree_->Branch("VDecayVtxZ",&VDecayVtxZ);
  tree_->Branch("VDecayVtxXE",&VDecayVtxXE);
  tree_->Branch("VDecayVtxYE",&VDecayVtxYE);
  tree_->Branch("VDecayVtxZE",&VDecayVtxZE);
  tree_->Branch("VDecayVtxXYE",&VDecayVtxXYE);
  tree_->Branch("VDecayVtxXZE",&VDecayVtxXZE);
  tree_->Branch("VDecayVtxYZE",&VDecayVtxYZE);

  tree_->Branch("pi1dxy",&pi1dxy);
  tree_->Branch("pi2dxy",&pi2dxy);
  tree_->Branch("pi1dz",&pi1dz);
  tree_->Branch("pi2dz",&pi2dz);

  tree_->Branch("pi1dxy_e",&pi1dxy_e);
  tree_->Branch("pi2dxy_e",&pi2dxy_e);
  tree_->Branch("pi1dz_e",&pi1dz_e);
  tree_->Branch("pi2dz_e",&pi2dz_e);

  tree_->Branch("mumC2",&mumC2);  
  tree_->Branch("mumNHits",&mumNHits);
  tree_->Branch("mumNPHits",&mumNPHits);
  tree_->Branch("mupC2",&mupC2);  
  tree_->Branch("mupNHits",&mupNHits);
  tree_->Branch("mupNPHits",&mupNPHits);
  tree_->Branch("mumdxy",&mumdxy);
  tree_->Branch("mupdxy",&mupdxy);
  tree_->Branch("muon_dca",&muon_dca);

  tree_->Branch("tri_Dim25",&tri_Dim25);
  tree_->Branch("tri_JpsiTk",&tri_JpsiTk);
  tree_->Branch("tri_JpsiTkTk",&tri_JpsiTkTk); 
 
  tree_->Branch("mumdz",&mumdz);
  tree_->Branch("mupdz",&mupdz);
  tree_->Branch("mu1soft",&mu1soft);
  tree_->Branch("mu2soft",&mu2soft);
  tree_->Branch("mu1tight",&mu1tight);
  tree_->Branch("mu2tight",&mu2tight);
  tree_->Branch("mu1PF",&mu1PF);
  tree_->Branch("mu2PF",&mu2PF);
  tree_->Branch("mu1loose",&mu1loose);
  tree_->Branch("mu2loose",&mu2loose);

}


// ------------ method called once each job just after ending the event loop  ------------
void JPsiKs0::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(JPsiKs0);

