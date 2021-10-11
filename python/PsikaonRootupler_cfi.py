import FWCore.ParameterSet.Config as cms

rootuple = cms.EDAnalyzer('JPsiKaon',
                          dimuons = cms.InputTag("slimmedMuons"),
                          Trak = cms.InputTag("packedPFCandidates"),
                          GenParticles = cms.InputTag("genParticles"),
                          packedGenParticles = cms.InputTag("packedGenParticles"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          bslabel = cms.InputTag("offlineBeamSpot"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          OnlyBest = cms.bool(False),
                          isMC = cms.bool(False),
                          OnlyGen = cms.bool(False),
                          )
