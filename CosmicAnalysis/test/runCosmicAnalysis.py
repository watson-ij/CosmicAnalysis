import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process('CosmicAnalysis',eras.Run3)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoMuon.TrackingTools.MuonServiceProxy_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '113X_dataRun3_Prompt_v2', '')
process.MessageLogger.cerr.FwkReport.reportEvery = 5000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
#process.maxEvents.input = cms.untracked.int32(10)
# Input source
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.skipEvents = cms.untracked.uint32(0)

from glob import glob
process.source.fileNames.extend(
    # ['/store/express/Commissioning2021/ExpressCosmics/FEVT/Express-v1/000/342/218/00000/007c63b2-8625-44ed-b1e6-a13ee77f2a42.root'],
    [f.replace('/eos/cms','') for f in glob('/eos/cms/store/express/Commissioning2021/ExpressCosmics/FEVT/Express-v1/000/342/218/00000/*')]
    # [('file:'+f) for f in glob('/eos/cms/store/express/Commissioning2021/ExpressCosmics/FEVT/Express-v1/000/342/218/00000/*')][:5]
)

process.options = cms.untracked.PSet()

process.TFileService = cms.Service("TFileService",fileName = cms.string("cosmics.root"))

process.CosmicAnalysis = cms.EDAnalyzer('CosmicAnalysis',
    process.MuonServiceProxy,
    gemRecHits = cms.InputTag("gemRecHits"),
    cscSegments = cms.InputTag("cscSegments"),
    muons = cms.InputTag("muons"),
)
process.p = cms.Path(process.CosmicAnalysis)
