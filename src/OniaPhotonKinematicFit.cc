// -*- C++ -*-
//
// Package:    OniaPhotonKinematicFit
// Class:      OniaPhotonKinematicFit
// 
/**

 Description: Kinematic Fit for Onia + Photon

 Implementation:
     Original work from Stefano Argiro, and the Torino group
     Adapter for MINIAOD by Alberto Sanchez-Hernandez
*/


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

///For kinematic fit:
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h> 
#include <DataFormats/PatCandidates/interface/PackedCandidate.h>

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
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"

#include <boost/foreach.hpp>
#include <string>

//
// class declaration
//

class OniaPhotonKinematicFit : public edm::EDProducer {
  public:
    explicit OniaPhotonKinematicFit(const edm::ParameterSet&);
    ~OniaPhotonKinematicFit() override {};
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    void produce(edm::Event&, const edm::EventSetup&) override;

  double upsilon_mass_;
  std::string product_name_;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> chi_Label;

  template<typename T>
    struct GreaterByVProb {
      typedef T first_argument_type;
      typedef T second_argument_type;
      bool operator()( const T & t1, const T & t2 ) const {
            return t1.userFloat("vProb") > t2.userFloat("vProb");
      }
    };
};

OniaPhotonKinematicFit::OniaPhotonKinematicFit(const edm::ParameterSet& iConfig) {
  chi_Label     = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("chi_cand"));
  upsilon_mass_ = iConfig.getParameter<double>("upsilon_mass");
  product_name_ = iConfig.getParameter<std::string>("product_name");
  produces<pat::CompositeCandidateCollection>(product_name_);
}

// ------------ method called to produce the data  ------------
void OniaPhotonKinematicFit::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  // Grab paramenters
  edm::Handle<pat::CompositeCandidateCollection> chiCandHandle;
  iEvent.getByToken(chi_Label, chiCandHandle);
  
  //Kinemati refit collection
  std::unique_ptr<pat::CompositeCandidateCollection> chicCompCandRefitColl(new pat::CompositeCandidateCollection);
  
  // Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB; 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); 
  
  int indexConversion=-1;

  for (pat::CompositeCandidateCollection::const_iterator chiCand=chiCandHandle->begin();chiCand!=chiCandHandle->end();++chiCand) { 
    indexConversion++;
    reco::TrackRef JpsiTk[2]={
      ( dynamic_cast<const pat::Muon*>(chiCand->daughter("dimuon")->daughter("muon1") ) )->innerTrack(),
      ( dynamic_cast<const pat::Muon*>(chiCand->daughter("dimuon")->daughter("muon2") ) )->innerTrack()
    };
      
    reco::TrackCollection convTracks;
    const reco::Track tk0=*(dynamic_cast<const pat::CompositeCandidate*>(chiCand->daughter("photon"))->userData<reco::Track>("track0"));
    const reco::Track tk1=*(dynamic_cast<const pat::CompositeCandidate*>(chiCand->daughter("photon"))->userData<reco::Track>("track1"));
    convTracks.push_back(tk0);
    convTracks.push_back(tk1);
    const reco::Vertex thePrimaryV=*(dynamic_cast<const pat::CompositeCandidate*>(chiCand->daughter("dimuon"))->userData<reco::Vertex>("PVwithmuons"));

    std::vector<reco::TransientTrack> MuMuTT;
    MuMuTT.push_back((*theB).build(&JpsiTk[0]));
    MuMuTT.push_back((*theB).build(&JpsiTk[1]));
      
    std::vector<reco::TransientTrack> EETT;
    EETT.push_back((*theB).build(convTracks[0]));
    EETT.push_back((*theB).build(convTracks[1]));
      
    const ParticleMass zero_mass(0);
    float zero_sigma = 1E-6;
      
    const ParticleMass eleMass(0.000511);
    float eleSigma = 1E-6;
      
    KinematicParticleFactoryFromTransientTrack pFactory;
    std::vector<RefCountedKinematicParticle> PhotonParticles;
    PhotonParticles.push_back(pFactory.particle(EETT[0],eleMass,float(0),float(0),eleSigma));
    PhotonParticles.push_back(pFactory.particle(EETT[1],eleMass,float(0),float(0),eleSigma));
      
    KinematicParticleVertexFitter fitter;
    RefCountedKinematicTree photonVertexFitTree;
    photonVertexFitTree = fitter.fit(PhotonParticles);
      
    if (!photonVertexFitTree->isValid()) { 
      edm::ParameterSet pSet;
      pSet.addParameter<double>("maxDistance", 3);
      pSet.addParameter<int>("maxNbrOfIterations", 10000);
      KinematicParticleVertexFitter fitter2(pSet);
      photonVertexFitTree = fitter2.fit(PhotonParticles);
    }
      
    if (photonVertexFitTree->isValid()) {
	  
      // now apply Photon mass constraint
      KinematicParticleFitter csFitterPhoton;
      KinematicConstraint * pho_c = new MassKinematicConstraint(zero_mass,zero_sigma);

      // add mass constraint to the photon fit to do a constrained fit:
      photonVertexFitTree->movePointerToTheTop();
      photonVertexFitTree = csFitterPhoton.fit(pho_c,photonVertexFitTree);
	  
      if (photonVertexFitTree->isValid()) {

	const ParticleMass muonMass(0.1056584);
	float muonSigma = muonMass*1E-6;
	      
	photonVertexFitTree->movePointerToTheTop();
	RefCountedKinematicParticle fitPhoton = photonVertexFitTree->currentParticle();
	      
	std::vector<RefCountedKinematicParticle> allChiBDaughters;
	allChiBDaughters.push_back(pFactory.particle (MuMuTT[0], muonMass, float(0), float(0), muonSigma));
	allChiBDaughters.push_back(pFactory.particle (MuMuTT[1], muonMass, float(0), float(0), muonSigma));
	allChiBDaughters.push_back(fitPhoton);
	      
	KinematicConstrainedVertexFitter constVertexFitter;
	      
	MultiTrackKinematicConstraint *upsilon_mtc = new  TwoTrackMassKinematicConstraint(upsilon_mass_);
	RefCountedKinematicTree ChiBTree = constVertexFitter.fit(allChiBDaughters,upsilon_mtc);

	if (!ChiBTree->isEmpty()) {
		  
	  ChiBTree->movePointerToTheTop();
          RefCountedKinematicParticle fitChiB = ChiBTree->currentParticle();
          RefCountedKinematicVertex ChiBDecayVertex = ChiBTree->currentDecayVertex();
		  
          if (fitChiB->currentState().isValid()) { //Get chib         

            float ChiBM_fit  = fitChiB->currentState().mass();
	    float ChiBPx_fit = fitChiB->currentState().kinematicParameters().momentum().x();
            float ChiBPy_fit = fitChiB->currentState().kinematicParameters().momentum().y();
	    float ChiBPz_fit = fitChiB->currentState().kinematicParameters().momentum().z();
            float ChiBVtxX_fit = ChiBDecayVertex->position().x();
            float ChiBVtxY_fit = ChiBDecayVertex->position().y();
            float ChiBVtxZ_fit = ChiBDecayVertex->position().z();
            float ChiBVtxP_fit = ChiSquaredProbability((double)(ChiBDecayVertex->chiSquared()),
                                                       (double)(ChiBDecayVertex->degreesOfFreedom()));
		      
            reco::CompositeCandidate recoChib(0, math::XYZTLorentzVector(ChiBPx_fit, ChiBPy_fit, ChiBPz_fit,
                                              sqrt(ChiBM_fit*ChiBM_fit + ChiBPx_fit*ChiBPx_fit + ChiBPy_fit*ChiBPy_fit +
                                              ChiBPz_fit*ChiBPz_fit)), math::XYZPoint(ChiBVtxX_fit,
                                              ChiBVtxY_fit, ChiBVtxZ_fit), 50551);
		      
	    pat::CompositeCandidate patChib(recoChib);
            patChib.addUserFloat("vProb",ChiBVtxP_fit);
            patChib.addUserInt("Index",indexConversion);  // this also holds the index of the current chiCand

            // lifetime using PV
            TVector3 vtx;
            TVector3 pvtx;
            VertexDistanceXY vdistXY;
            reco::Vertex myVertex = *ChiBDecayVertex;

            vtx.SetXYZ(ChiBVtxX_fit, ChiBVtxY_fit, 0);
            pvtx.SetXYZ(thePrimaryV.position().x(), thePrimaryV.position().y(), 0);
            TVector3 pperp(ChiBPx_fit, ChiBPy_fit, 0);
            AlgebraicVector3 vpperp(pperp.x(), pperp.y(), 0);

            TVector3 vdiff = vtx - pvtx;
            double cosAlpha = vdiff.Dot(pperp) / (vdiff.Perp() * pperp.Perp());
            Measurement1D distXY = vdistXY.distance(myVertex, thePrimaryV);
            double ctauPV = distXY.value() * cosAlpha * ChiBM_fit / pperp.Perp();
            GlobalError v1e = myVertex.error();
            GlobalError v2e = thePrimaryV.error();
            AlgebraicSymMatrix33 vXYe = v1e.matrix() + v2e.matrix();
            double ctauErrPV = sqrt(ROOT::Math::Similarity(vpperp, vXYe)) * ChiBM_fit / (pperp.Perp2());
            patChib.addUserFloat("ctauPV",ctauPV);
            patChib.addUserFloat("ctauErrPV",ctauErrPV);
            patChib.addUserFloat("cosAlpha",cosAlpha);

            //get first muon
            bool child = ChiBTree->movePointerToTheFirstChild();
            RefCountedKinematicParticle fitMu1 = ChiBTree->currentParticle();
            if (!child) break;

            float mu1M_fit  = fitMu1->currentState().mass();
            float mu1Q_fit  = fitMu1->currentState().particleCharge();
            float mu1Px_fit = fitMu1->currentState().kinematicParameters().momentum().x();
	    float mu1Py_fit = fitMu1->currentState().kinematicParameters().momentum().y();
            float mu1Pz_fit = fitMu1->currentState().kinematicParameters().momentum().z();
	    reco::CompositeCandidate recoMu1(mu1Q_fit, math::XYZTLorentzVector(mu1Px_fit, mu1Py_fit, mu1Pz_fit, 
                                             sqrt(mu1M_fit*mu1M_fit + mu1Px_fit*mu1Px_fit + mu1Py_fit*mu1Py_fit + 
                                             mu1Pz_fit*mu1Pz_fit)), math::XYZPoint(ChiBVtxX_fit, ChiBVtxY_fit, ChiBVtxZ_fit), 13);
	    pat::CompositeCandidate patMu1(recoMu1);
		      
            //get second muon
            child = ChiBTree->movePointerToTheNextChild();
	    RefCountedKinematicParticle fitMu2 = ChiBTree->currentParticle();
            if (!child) break;

	    float mu2M_fit  = fitMu2->currentState().mass();
            float mu2Q_fit  = fitMu2->currentState().particleCharge();
            float mu2Px_fit = fitMu2->currentState().kinematicParameters().momentum().x();
	    float mu2Py_fit = fitMu2->currentState().kinematicParameters().momentum().y();
	    float mu2Pz_fit = fitMu2->currentState().kinematicParameters().momentum().z();
            reco::CompositeCandidate recoMu2(mu2Q_fit, math::XYZTLorentzVector(mu2Px_fit, mu2Py_fit, mu2Pz_fit,
                                             sqrt(mu2M_fit*mu2M_fit + mu2Px_fit*mu2Px_fit + mu2Py_fit*mu2Py_fit + 
                                             mu2Pz_fit*mu2Pz_fit)), math::XYZPoint(ChiBVtxX_fit, ChiBVtxY_fit, ChiBVtxZ_fit), 13);
            pat::CompositeCandidate patMu2(recoMu2);
		      		  
            //Define Onia from two muons
            pat::CompositeCandidate ups;
            ups.addDaughter(patMu1,"muon1");
            ups.addDaughter(patMu2,"muon2");	
            ups.setP4(patMu1.p4()+patMu2.p4());
		  
            //get photon
            child = ChiBTree->movePointerToTheNextChild();
            RefCountedKinematicParticle fitGamma = ChiBTree->currentParticle();
	    if (!child) break;
		      
	    float gammaM_fit  = fitGamma->currentState().mass();
            float gammaPx_fit = fitGamma->currentState().kinematicParameters().momentum().x();
	    float gammaPy_fit = fitGamma->currentState().kinematicParameters().momentum().y();
            float gammaPz_fit = fitGamma->currentState().kinematicParameters().momentum().z();
	    reco::CompositeCandidate recoGamma(0, math::XYZTLorentzVector(gammaPx_fit, gammaPy_fit, gammaPz_fit, 
                                               sqrt(gammaM_fit*gammaM_fit + gammaPx_fit*gammaPx_fit + gammaPy_fit*gammaPy_fit +
                                               gammaPz_fit*gammaPz_fit)), math::XYZPoint(ChiBVtxX_fit, ChiBVtxY_fit, ChiBVtxZ_fit), 22);
            pat::CompositeCandidate patGamma(recoGamma);

            patChib.addDaughter(ups,"dimuon");
            patChib.addDaughter(patGamma,"photon");

	    chicCompCandRefitColl->push_back(patChib);
	  }
	}
      }
    }
  }  
  // End kinematic fit

  // ...ash not sorted, since we will use the best un-refitted candidate
  // now sort by vProb
  //OniaPhotonKinematicFit::GreaterByVProb<pat::CompositeCandidate> vPComparator;
  //std::sort(chicCompCandRefitColl->begin(),chicCompCandRefitColl->end(), vPComparator);

  iEvent.put(std::move(chicCompCandRefitColl),product_name_); 
  
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void OniaPhotonKinematicFit::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(OniaPhotonKinematicFit);
