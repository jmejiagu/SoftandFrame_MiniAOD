import FWCore.ParameterSet.Config as cms

### ==== This is our version of the patMuonsWithTrigger using MINIAOD

### unpack them
unpackedPatTrigger = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
  patTriggerObjectsStandAlone = cms.InputTag( 'slimmedPatTrigger' ),
  triggerResults              = cms.InputTag( 'TriggerResults::HLT' ),
  unpackFilterLabels          = cms.bool( True )
)

### then perform a match for all HLT triggers of interest
PATmuonTriggerMatchHLT = cms.EDProducer( "PATTriggerMatcherDRDPtLessByR",
    src     = cms.InputTag( "slimmedMuons" ),
    matched = cms.InputTag( "unpackedPatTrigger" ),
    matchedCuts = cms.string(""),
    maxDPtRel = cms.double( 0.5 ),
    maxDeltaR = cms.double( 0.5 ),
    resolveAmbiguities    = cms.bool( True ),
    resolveByMatchQuality = cms.bool( True )
)

PATmuonMatchHLTL2   = PATmuonTriggerMatchHLT.clone(matchedCuts = cms.string('coll("hltL2MuonCandidates")'), maxDeltaR = 0.3, maxDPtRel = 10.0)           #maxDeltaR Changed accordingly to Zoltan tuning. It was: 1.2
PATmuonMatchHLTL3   = PATmuonTriggerMatchHLT.clone(matchedCuts = cms.string('coll("hltIterL3FromL2MuonCandidates")'), maxDeltaR = 0.1, maxDPtRel = 10.0) #maxDeltaR Changed accordingly to Zoltan tuning. It was: 0.5
PATmuonMatchHLTL3v2 = PATmuonTriggerMatchHLT.clone(matchedCuts = cms.string('coll("hltIterL3MuonCandidates")'), maxDeltaR = 0.1, maxDPtRel = 10.0)       #needed because name changed
PATmuonMatchHLTL3T  = PATmuonTriggerMatchHLT.clone(matchedCuts = cms.string('coll("hltGlbTrkMuonCands")'), maxDeltaR = 0.1, maxDPtRel = 10.0)            #maxDeltaR Changed accordingly to Zoltan tuning. It was: 0.5
PATmuonMatchHLTTkMu = PATmuonTriggerMatchHLT.clone(matchedCuts = cms.string('coll("hltHighPtTkMuonCands")'), maxDeltaR = 0.1, maxDPtRel = 10.0)          #maxDeltaR Changed accordingly to Zoltan tuning. It was: 0.5

slimmedMuonsTriggerMatchers1Mu = cms.Sequence(
      PATmuonMatchHLTL2 +
      PATmuonMatchHLTL3 +
      PATmuonMatchHLTL3v2 +
      PATmuonMatchHLTL3T +
      PATmuonMatchHLTTkMu
)

slimmedMuonsTriggerMatchers1MuInputTags = [
    cms.InputTag('PATmuonMatchHLTL2'),
    cms.InputTag('PATmuonMatchHLTL3'),
    cms.InputTag('PATmuonMatchHLTL3v2'),
    cms.InputTag('PATmuonMatchHLTL3T'),
    cms.InputTag('PATmuonMatchHLTTkMu'),
]

### Embed them
slimmedMuonsWithTrigger = cms.EDProducer( "PATTriggerMatchMuonEmbedder",
    src     = cms.InputTag(  "slimmedMuons" ),
    matches = cms.VInputTag()
)
slimmedMuonsWithTrigger.matches += slimmedMuonsTriggerMatchers1MuInputTags

### Run the whole trigger Sequence 
slimmedMuonsWithTriggerSequence = cms.Sequence(
    unpackedPatTrigger *
    slimmedMuonsTriggerMatchers1Mu *
    slimmedMuonsWithTrigger
)
