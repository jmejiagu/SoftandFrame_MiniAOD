import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v24', '')# for 2018
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v20', '')# for 2017
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v27', '')# for 2016
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')

## Message Logger and Event range
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

        #MiniAOD UltraLegacy2018
        #'/store/data/Run2018A/Charmonium/MINIAOD/12Nov2019_UL2018_rsb-v1/10000/08F41CB9-8F1F-D44F-A5FC-D17E38328C4C.root',
        
        #MiniAOD UltraLegacy2017
        #'/store/data/Run2017F/Charmonium/MINIAOD/09Aug2019_UL2017-v1/20000/00BACB48-9B0F-8F48-A68B-2F08A3E9E681.root',
        
        #MiniAOD UltraLegacy2016
        #'/store/data/Run2016G/Charmonium/MINIAOD/21Feb2020_UL2016-v1/30000/00013A18-278D-5B48-9BEF-1083A8F5C9D7.root',

        #MiniAOD UltraLegacy BsToJpsiPhi MC sample 2017
        #'/store/mc/RunIISummer19UL17MiniAOD/BsToJpsiPhi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/106X_mc2017_realistic_v6-v1/10000/0C845AFD-DC22-CC4D-8BBC-32F7404B014C.root',

        #MiniAOD MC RelValPsi2SToJPsiPiPi_14
        'file:/eos/cms/store/relval/CMSSW_12_1_0_pre2/RelValPsi2SToJPsiPiPi_14/MINIAODSIM/121X_mcRun3_2021_realistic_v1-v3/10000/934b45ff-3f99-4f69-9648-bfd129ab7d90.root',
             
 )
)

process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
                                        triggerConditions = cms.vstring('HLT_Dimuon20_Jpsi_Barrel_Seagulls_v*',
                                                                        'HLT_Dimuon25_Jpsi_v*',
                                                                        'HLT_DoubleMu4_JpsiTrkTrk_Displaced_v*',
                                                                        'HLT_DoubleMu4_JpsiTrk_Displaced_v*',
                                                                        'HLT_DoubleMu4_Jpsi_Displaced_v*'                                   
                                                                       ),
                                        hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                        l1tResults = cms.InputTag( "" ),
                                        throw = cms.bool(False)
                                        )

process.rootuple = cms.EDAnalyzer('JPsiTrkTrk',
                          dimuons = cms.InputTag("slimmedMuons"),
                          Trak = cms.InputTag("packedPFCandidates"),
                          #Trak_lowpt = cms.InputTag("lostTracks"),
                          GenParticles = cms.InputTag("genParticles"),
                          packedGenParticles = cms.InputTag("packedGenParticles"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          bslabel = cms.InputTag("offlineBeamSpot"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          OnlyBest = cms.bool(False),
                          isMC = cms.bool(False),
                          OnlyGen = cms.bool(False),
                          #Trkmass           = cms.double(0.493677),
                          Trkmass           = cms.double(0.13957018),
                          TrkTrkMasscut     = cms.vdouble(0.350,4.0),        
                          BarebMasscut      = cms.vdouble(3.2,4.4),
                          bMasscut          = cms.vdouble(3.3,4.0),        
                          )


#process.load("myAnalyzers.JPsiKsPAT.PsiphiRootupler_cfi")
#process.rootuple.isMC = cms.bool(True) 
#process.rootuple.GenParticles = cms.InputTag("prunedGenParticles") 

process.TFileService = cms.Service("TFileService",
        #fileName = cms.string('Rootuple_BstoJpsiphi_MC2017UL_MiniAOD.root'),
        fileName = cms.string('Rootuple_Psi2StoJpsipipi_MCRelVal_MiniAOD.root'),                                  
)


process.mySequence = cms.Sequence(
                                   process.triggerSelection *
                                   process.rootuple
				   )
#process.p = cms.Path(process.mySequence)

#process.p = cms.Path(process.triggerSelection*process.rootuple)
process.p = cms.Path(process.rootuple)




