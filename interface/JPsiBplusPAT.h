#ifndef _JPsiBplusPAT_h
#define _JPsiBplusPAT_h

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

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
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

class JPsiBplusPAT : public edm::EDAnalyzer {
public:
  explicit JPsiBplusPAT(const edm::ParameterSet&);
  ~JPsiBplusPAT();
  void fillPsi(const reco::Candidate& genpsi);
  void fillV0(const reco::Candidate& genv0);
  int const getMuCat(reco::Muon const& muon) const;
  bool const HasGoodME11(reco::Muon const& muon, double const dxdzCut) const;
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void printout(const RefCountedKinematicVertex& myVertex) const;
  void printout(const RefCountedKinematicParticle& myParticle) const;
  void printout(const RefCountedKinematicTree& myTree) const;

  int PdgIDatTruthLevel(reco::TrackRef Track, edm::Handle<reco::GenParticleCollection> genParticles, int &ParentID);
    
  void CheckL1Triggers(const edm::Event& iEvent, const edm::EventSetup& iSetup, std::string &TrigListNameL1Tmp);
  void MatchMuonWithTriggers(const pat::Muon &iMuon, const std::vector<std::string>& TrigList, std::string &TrigListNameTmp);
  void CheckHLTTriggers(const std::vector<std::string>& TrigList);
    
  void MatchMuonWithL1L2(const pat::Muon &iMuon, const std::vector<std::string>& TrigListL1L2, std::string &TrigListNameL1L2Tmp);
  float Myctau(const RefCountedKinematicParticle &CandMC, const RefCountedKinematicVertex &DecayVertexMC, 
	       const GlobalPoint &PVPtmp, const GlobalError &PVEtmp,float mass_tmp, 
	       float &ctau2Dtmp, float &ctauEtemp, float &ctauE2Dtemp );
    
    // ----------member data ---------------------------
  std::string hltTriggerResults_;
  std::string vtxSample;
  std::string genParticles_;
  std::string muonType;
  std::string muonTypeForPAT;
  bool doMC_;
  TTree*      tree_;
  int mupCategory;
  int mumCategory;
  int mupME1Clean;
  int mumME1Clean;
  
  
  
  std::vector<float>       *priRfVtxX, *priRfVtxY, *priRfVtxZ, *priRfVtxXE, *priRfVtxYE, *priRfVtxZE, *priRfVtxCL;
  std::vector<float>       *priRfVtxXYE, *priRfVtxXZE, *priRfVtxYZE;
  std::vector<int>         *priRfNTrkDif;
  //std::vector<float>       *bctau, *bctau2D; //*bctauBS, *bctauBS2D, *bctauRf, *bctau2DRf;
  //std::vector<float>       *bctauE, *bctau2DE;// *bctauBSE, *bctauBS2DE, *bctauRfE, *bctau2DRfE;
 
  std::vector<float>       *mumC2;
  std::vector<int>         *mumCat, *mumAngT, *mumNHits, *mumNPHits; 
  std::vector<float>       *mupC2;
  std::vector<int>         *mupCat, *mupAngT, *mupNHits, *mupNPHits;
  std::vector<float>       *mumdxy, *mupdxy, *mumdz, *mupdz;

  std::vector<bool>        *mu1soft, *mu2soft, *mu1tight, *mu2tight;  
  std::vector<bool>        *mu1PF, *mu2PF, *mu1loose, *mu2loose;  
 
  int                      muAcc, muTrig, weight;
  // *************************************
  unsigned int             nB;
  unsigned int             nMu;
    
  std::vector<float>       *B_mass, *B_px, *B_py, *B_pz, *B_charge;
  std::vector<float>       *B_k_px, *B_k_py, *B_k_pz,  *B_k_charge1; 
  std::vector<int>         *B_k_parentId1, *B_k_pId1;

  std::vector<float>       *B_J_mass, *B_J_px, *B_J_py, *B_J_pz;

  std::vector<float>       *B_J_pt1, *B_J_px1, *B_J_py1, *B_J_pz1;
  std::vector<float>       *B_J_pt2, *B_J_px2, *B_J_py2, *B_J_pz2;
  std::vector<int>         *B_J_charge1, *B_J_charge2;

 // vertice primario CON mayor Pt
  unsigned int             nVtx;
  float                    priVtxX, priVtxY, priVtxZ, priVtxXE, priVtxYE, priVtxZE, priVtxCL;
  float                    priVtxXYE, priVtxXZE, priVtxYZE;
 
  // vertice primario CON constrain de Beamspot y mayor Pt 
  float                    priVtxXBS, priVtxYBS, priVtxZBS, priVtxXBSE, priVtxYBSE, priVtxZBSE, priVtxCLBS;
  float                    priVtxXYBSE, priVtxXZBSE, priVtxYZBSE;

// todos los vertices primarios 
//  std::vector<float>       *pVtxX,  *pVtxY,  *pVtxZ,  *pVtxXE,  *pVtxYE,  *pVtxZE,  *pVtxCL;
//  std::vector<float>       *pVtxXYE, *pVtxXZE, *pVtxYZE;

// vertice primario CON mejor pointin-angle
 std::vector<float>          *pVtxIPX,  *pVtxIPY, *pVtxIPZ, *pVtxIPXE, *pVtxIPYE, *pVtxIPZE, *pVtxIPCL;
 std::vector<float>          *pVtxIPXYE,  *pVtxIPXZE, *pVtxIPYZE;


// todos los vertices primarios CON constrain de Beamspot
//  std::vector<float>       *pVtxBSX,  *pVtxBSY,  *pVtxBSZ,  *pVtxBSXE,  *pVtxBSYE,  *pVtxBSZE,  *pVtxBSCL;
//  std::vector<float>       *pVtxBSXYE, *pVtxBSXZE, *pVtxBSYZE;

 // vertice primario CON constrain de Beamspot y mejor pointin-angle
  std::vector<float>        *pVtxBSIPX,  *pVtxBSIPY,  *pVtxBSIPZ, *pVtxBSIPXE, *pVtxBSIPYE, *pVtxBSIPZE, *pVtxBSIPCL;
  std::vector<float>        *pVtxBSIPXYE,  *pVtxBSIPXZE,  *pVtxBSIPYZE;

// ************ esta es la informacion concerniente al beamspot *************           
  double                    PVXBS, PVYBS, PVZBS, PVXBSE, PVYBSE, PVZBSE;
  double                    PVXYBSE, PVXZBSE, PVYZBSE;
  
  // ********************************** ************************************************************************
 

  std::vector<float>       *B_chi2, *B_J_chi2;
  std::vector<float>       *B_Prob, *B_J_Prob;

  std::vector<float>       *B_DecayVtxX,  *B_DecayVtxY,  *B_DecayVtxZ;
  std::vector<double>      *B_DecayVtxXE, *B_DecayVtxYE, *B_DecayVtxZE;
  std::vector<double>      *B_DecayVtxXYE, *B_DecayVtxXZE, *B_DecayVtxYZE;
 
  std::vector<float>       *B_J_DecayVtxX,   *B_J_DecayVtxY,   *B_J_DecayVtxZ;
  std::vector<float>       *B_J_DecayVtxXE,  *B_J_DecayVtxYE,  *B_J_DecayVtxZE;
  std::vector<float>       *B_J_DecayVtxXYE, *B_J_DecayVtxXZE, *B_J_DecayVtxYZE;

  char triggersL[10000], triggersL1L[10000];

  char triggersMuP[10000],     triggersMuM[10000] ;
  char triggersL1L2_MuP[10000], triggersL1L2_MuM[10000];
  
  char triggersL1[10000];
  int  nTrgL, nTrgL1L,  nMuonTrgL,  nMuonPTrgL,        nMuonMTrgL;
  int  ntriggersL1L2_MuP, ntriggersL1L2_MuM;

  //std::vector<std::string>         *triggersL1L;
  std::vector<std::string>         *triggersMuPL, *triggersMuML;
  std::vector<std::string>         *triggersL1L2_MuPL, *triggersL1L2_MuML;
 
  int  run, event;


};
#endif
