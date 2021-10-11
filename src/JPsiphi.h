#ifndef _JPsiphi_h
#define _JPsiphi_h

// system include files
#include <memory>

// user include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"


#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // muy importante para MiniAOD
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "TFile.h"
#include "TTree.h"
#include <vector>
#include "TLorentzVector.h"
#include "TVector3.h"
#include <utility>
#include <string>



//
// class decleration
//

class JPsiphi : public edm::EDAnalyzer {
public:
  explicit JPsiphi(const edm::ParameterSet&);
  ~JPsiphi();
  void fillPsi(const reco::Candidate& genpsi);
  void fillV0(const reco::Candidate& genv0);
  int const getMuCat(reco::Muon const& muon) const;
  bool IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu);
  bool   isAncestor(const reco::Candidate*, const reco::Candidate*);
  double GetLifetime(TLorentzVector, TVector3, TVector3);
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void printout(const RefCountedKinematicVertex& myVertex) const;
  void printout(const RefCountedKinematicParticle& myParticle) const;
  void printout(const RefCountedKinematicTree& myTree) const;
  
  
    // ----------member data ---------------------------
  edm::EDGetTokenT<edm::View<pat::Muon>> dimuon_Label;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trakCollection_label;
  edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packedGenToken_;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
  edm::EDGetTokenT<reco::BeamSpot> BSLabel_;


  bool OnlyBest_;
  bool isMC_;
  bool OnlyGen_;
  // bool doMC_;
  //std::string genParticles_;

  TTree*      tree_;
  int mupCategory;
  int mumCategory;
  int mupME1Clean;
  int mumME1Clean;
  
  std::vector<float>       *mumC2;
  std::vector<int>         *mumAngT, *mumNHits, *mumNPHits; 
  std::vector<float>       *mupC2;
  std::vector<int>         *mupAngT, *mupNHits, *mupNPHits;
  std::vector<float>       *mumdxy, *mupdxy, *mumdz, *mupdz;
  std::vector<float>       *muon_dca;

  std::vector<int>         *tri_Dim25, *tri_JpsiTk, *tri_JpsiTkTk; 

  std::vector<bool>        *mu1soft, *mu2soft, *mu1tight, *mu2tight;  
  std::vector<bool>        *mu1PF, *mu2PF, *mu1loose, *mu2loose;  

  int                     muAcc, muTrig, weight;

  // *************************************
  unsigned int             nB;
  unsigned int             nMu;
  
  std::vector<float>       *B_mass, *B_px, *B_py, *B_pz;
  
  std::vector<float>       *B_phi_mass;
  std::vector<float>       *B_phi_px1, *B_phi_py1, *B_phi_pz1;
  std::vector<float>       *B_phi_px2, *B_phi_py2, *B_phi_pz2;

  std::vector<float>       *B_phi_px1_track, *B_phi_py1_track, *B_phi_pz1_track;
  std::vector<float>       *B_phi_px2_track, *B_phi_py2_track, *B_phi_pz2_track;

  std::vector<int>         *B_phi_charge1, *B_phi_charge2;
  std::vector<float>       *k1dxy, *k2dxy, *k1dz, *k2dz;
  std::vector<float>       *k1dxy_e, *k2dxy_e, *k1dz_e, *k2dz_e;
  std::vector<float>       *k1InnerHits, *k2InnerHits;

  std::vector<float>       *B_J_mass, *B_J_px, *B_J_py, *B_J_pz;
  //std::vector<float>     *B_J_pt1;
  std::vector<float>       *B_J_px1, *B_J_py1, *B_J_pz1;
  //std::vector<float>     *B_J_pt2;
  std::vector<float>       *B_J_px2, *B_J_py2, *B_J_pz2;
  std::vector<int>         *B_J_charge1, *B_J_charge2;

  unsigned int             nVtx;
  float                    priVtxX, priVtxY, priVtxZ, priVtxXE, priVtxYE, priVtxZE, priVtxCL;
  float                    priVtxXYE, priVtxXZE, priVtxYZE;
  
  // ********************************** ************************************************************************

  std::vector<float>       *B_chi2, *B_J_chi2;
  std::vector<float>       *B_Prob, *B_J_Prob;

  std::vector<float>       *B_DecayVtxX,  *B_DecayVtxY,  *B_DecayVtxZ;
  std::vector<double>      *B_DecayVtxXE, *B_DecayVtxYE, *B_DecayVtxZE;
  std::vector<double>      *B_DecayVtxXYE, *B_DecayVtxXZE, *B_DecayVtxYZE;

  int  run, event;
  int lumiblock;

  TLorentzVector gen_b_p4,gen_phi_p4,gen_kaon1_p4,gen_kaon2_p4,gen_jpsi_p4,gen_muon1_p4,gen_muon2_p4;
  TVector3       gen_b_vtx,gen_jpsi_vtx;
  float          gen_b_ct;

};
#endif
