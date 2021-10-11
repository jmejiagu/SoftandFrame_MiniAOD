import FWCore.ParameterSet.Config as cms

rootuple = cms.EDAnalyzer('JPsiKs0',
                          dimuons = cms.InputTag("slimmedMuons"),
                          Trak = cms.InputTag("packedPFCandidates"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          secundaryVerticesPtr = cms.InputTag("slimmedKshortVertices"),
                          bslabel = cms.InputTag("offlineBeamSpot"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          OnlyBest = cms.bool(False),
                          isMC = cms.bool(False),
                          OnlyGen = cms.bool(False),
                          )
