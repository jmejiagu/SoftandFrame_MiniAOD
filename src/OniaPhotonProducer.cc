#include <myAnalyzers/JPsiKsPAT/src/OniaPhotonProducer.h>

OniaPhotonProducer::OniaPhotonProducer(const edm::ParameterSet& ps):
  dimuon_Label(consumes<pat::CompositeCandidateCollection>(ps.getParameter< edm::InputTag>("dimuons"))),
  photon_Label(consumes<pat::CompositeCandidateCollection>(ps.getParameter< edm::InputTag>("conversions"))),
  pi0OnlineSwitch_(ps.getParameter<bool>("pi0OnlineSwitch")),
  deltaMass_(ps.getParameter<std::vector<double> >("deltaMass")),
  dzMax_(ps.getParameter<double>("dzmax")),
  triggerMatch_(ps.getParameter<bool>("triggerMatch"))
{
  produces<pat::CompositeCandidateCollection>();
  candidates = 0;
  delta_mass_fail = 0;
  dz_cut_fail = 0;
  pizero_fail = 0;
}
 

void OniaPhotonProducer::produce(edm::Event& event, const edm::EventSetup& esetup){
  std::unique_ptr<pat::CompositeCandidateCollection> chiCandColl(new pat::CompositeCandidateCollection);

  edm::Handle<pat::CompositeCandidateCollection> dimuons;
  event.getByToken(dimuon_Label,dimuons);

  edm::Handle<pat::CompositeCandidateCollection> conversions;
  event.getByToken(photon_Label,conversions);

  // Note: since Dimuon cand are sorted by decreasing vertex probability then the first chi cand is the one associated with the "best" dimuon 
  for (pat::CompositeCandidateCollection::const_iterator  dimuonCand = dimuons->begin(); dimuonCand!= dimuons->end(); ++dimuonCand){

     // use only trigger-matched Jpsi or Upsilon if so requested 
     if (triggerMatch_){
         if (!dimuonCand->userInt("isTriggerMatched")) continue; 
     }

     // loop on conversion candidates, make chi cand
     for (pat::CompositeCandidateCollection::const_iterator conv = conversions->begin(); conv!= conversions->end(); ++conv){

	pat::CompositeCandidate chiCand = makeChiCandidate(*dimuonCand, *conv);
    
	if (!cutDeltaMass(chiCand,*dimuonCand)){
	   delta_mass_fail++;
	   continue;
	}
	const reco::Vertex *ipv = dimuonCand->userData<reco::Vertex>("commonVertex");
    	float dz = fabs(Getdz(*conv,ipv->position()));              // onia2mumu stores vertex as userData
	chiCand.addUserFloat("dz",dz);

	if (!cutdz(dz)){
	   dz_cut_fail++;	
	   continue;
	}

        int flags = (conv->userInt("flags")%32);
        bool pi0_fail = flags&8;
        if (pi0OnlineSwitch_ && pi0_fail) {
           pizero_fail++;
           continue;
        }

	chiCandColl->push_back(chiCand);
	candidates++;    
     }
  }
  event.put(std::move(chiCandColl));
}

float OniaPhotonProducer::Getdz(const pat::CompositeCandidate& c, const reco::Candidate::Point &p) {

  const reco::Candidate::LorentzVector& mom = c.p4();
  const reco::Candidate::Point& vtx = c.vertex();
  
  double dz = (vtx.Z()-p.Z()) - ((vtx.X()-p.X())*mom.X()+(vtx.Y()-p.Y())*mom.Y())/mom.Rho() * mom.Z()/mom.Rho();
  return (float) dz;  
  
}

void OniaPhotonProducer::endJob(){
  std::cout << "###########################" << std::endl;
  std::cout << "Chi Candidate producer report:" << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Delta mass fail: " << delta_mass_fail << std::endl;
  std::cout << "Dz fail:         " << dz_cut_fail << std::endl;
  std::cout << "Pi0 fail:        " << pizero_fail << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Found " << candidates << " Chi candidates." << std::endl;
  std::cout << "###########################" << std::endl;
}
  
const pat::CompositeCandidate OniaPhotonProducer::makeChiCandidate(const pat::CompositeCandidate& dimuon, 
				  const pat::CompositeCandidate& photon){
  pat::CompositeCandidate chiCand;
  chiCand.addDaughter(dimuon,"dimuon");
  chiCand.addDaughter(photon,"photon");
  const reco::Vertex *ipv = dimuon.userData<reco::Vertex>("commonVertex");
  chiCand.setVertex(ipv->position());
  reco::Candidate::LorentzVector vChic = dimuon.p4() + photon.p4();
  chiCand.setP4(vChic);
  return chiCand;
}

// check if the mass difference is in desired range
bool OniaPhotonProducer::cutDeltaMass(const pat::CompositeCandidate& chiCand,
				   const pat::CompositeCandidate& dimuonCand){
  float deltam = chiCand.p4().M() - dimuonCand.p4().M();
  float m1     = deltaMass_[0];
  float m2     = deltaMass_[1];
  return (deltam > m1 && deltam < m2);
}

//define this as a plug-in
DEFINE_FWK_MODULE(OniaPhotonProducer);
