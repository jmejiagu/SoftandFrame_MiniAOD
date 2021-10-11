#ifndef _Bc2SPAT_h
#define _Bc2SPAT_h

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

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // muy importante para MiniAOD

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

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


//
// class decleration
//

class Bc2SPAT : public edm::EDAnalyzer {
public:
  explicit Bc2SPAT(const edm::ParameterSet&);
  ~Bc2SPAT();
  void fillPsi(const reco::Candidate& genpsi);
  void fillV0(const reco::Candidate& genv0);
  //int const getMuCat(reco::Muon const& muon) const;
  bool IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu);
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void printout(const RefCountedKinematicVertex& myVertex) const;
  void printout(const RefCountedKinematicParticle& myParticle) const;
  void printout(const RefCountedKinematicTree& myTree) const;

  //void MatchMuonWithTriggers(const pat::Muon &iMuon, const std::vector<std::string>& TrigList, std::string &TrigListNameTmp);
  void CheckHLTTriggers(const std::vector<std::string>& TrigList);

  // ----------member data ---------------------------
  
  edm::EDGetTokenT<edm::View<pat::Muon>> dimuon_Label;
  //edm::EDGetTokenT<std::vector<pat::GenericParticle>> trakCollection_label;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trakCollection_label;
  //edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trakCollection_label_lowpt;
  //edm::EDGetTokenT<pat::PackedCandidateCollection> trakCollection_label;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
  edm::EDGetTokenT<reco::BeamSpot> BSLabel_;

 
  bool OnlyBest_;
  bool isMC_;
  bool OnlyGen_;

  TTree*      tree_;
  int mupCategory;
  int mumCategory;
  int mupME1Clean;
  int mumME1Clean;

  std::vector<int>         *tri_Dim25, *tri_Dim20, *tri_JpsiTk; 

  int                      muAcc, muTrig, weight;
  // *************************************
  unsigned int             nB;
  unsigned int             nMu;

  
  std::vector<float>       *Bc2s_mass, *Bc2s_px, *Bc2s_py, *Bc2s_pz;
  std::vector<float>       *Bc2s_pi2_px_track, *Bc2s_pi2_py_track, *Bc2s_pi2_pz_track;
  std::vector<float>       *Bc2s_pi3_px_track, *Bc2s_pi3_py_track, *Bc2s_pi3_pz_track;
  std::vector<int>         *Bc2s_pi2_charge1, *Bc2s_pi3_charge1;
  std::vector<float>       *Bc2s_Prob;

  std::vector<float>       *B_mass, *B_px, *B_py, *B_pz, *B_charge;
  std::vector<int>         *B_k_charge1; 
  std::vector<float>       *B_k_px_track, *B_k_py_track, *B_k_pz_track;

  std::vector<float>       *B_J_mass, *B_J_px, *B_J_py, *B_J_pz;

  std::vector<float>       *B_J_px1, *B_J_py1, *B_J_pz1;
  std::vector<float>       *B_J_px2, *B_J_py2, *B_J_pz2;
  std::vector<int>         *B_J_charge1, *B_J_charge2;

  // vertice primario CON mayor Pt
  unsigned int             nVtx;
  float                    priVtxX, priVtxY, priVtxZ, priVtxXE, priVtxYE, priVtxZE, priVtxCL;
  float                    priVtxXYE, priVtxXZE, priVtxYZE;

   // vertice primario CON mejor pointin-angle
   std::vector<float>          *pVtxIPX,  *pVtxIPY, *pVtxIPZ, *pVtxIPXE, *pVtxIPYE, *pVtxIPZE, *pVtxIPCL;
   std::vector<float>          *pVtxIPXYE,  *pVtxIPXZE, *pVtxIPYZE;
  
  // ********************************** ************************************************************************

  std::vector<float>       *B_Prob, *B_J_Prob;

  std::vector<float>       *B_DecayVtxX,  *B_DecayVtxY,  *B_DecayVtxZ;
  std::vector<double>      *B_DecayVtxXE, *B_DecayVtxYE, *B_DecayVtxZE;
  std::vector<double>      *B_DecayVtxXYE, *B_DecayVtxXZE, *B_DecayVtxYZE;
 
  std::vector<float>       *B_J_DecayVtxX,   *B_J_DecayVtxY,   *B_J_DecayVtxZ;
  std::vector<float>       *B_J_DecayVtxXE,  *B_J_DecayVtxYE,  *B_J_DecayVtxZE;
  std::vector<float>       *B_J_DecayVtxXYE, *B_J_DecayVtxXZE, *B_J_DecayVtxYZE;

  UInt_t trigger;

  int  run, event;
  int   lumiblock;

};
#endif
