import FWCore.ParameterSet.Config as cms


from CSCTFUtils.OfflineDQMCSCTF.offline_dqm_csctf_cfi import *


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
                                        SkipEvent = cms.untracked.vstring('ProductNotFound')
                                        )


##testing

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
# Other statements
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'GR_P_V56', '')

#globalTagValue = 'GR_H_V58'
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.connect = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'
#process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
#process.GlobalTag.globaltag = globalTagValue+'::All'
#es_prefer_GlobalTag = cms.ESPrefer('GlobalTag')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	'/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/248/003/00000/DCB92542-0714-E511-989A-02163E01438F.root',
    ),
    secondaryFileNames = cms.untracked.vstring(
	'/store/data/Run2015A/SingleMu/RAW/v1/000/248/003/00000/B424BBCB-6512-E511-8E78-02163E012852.root'
    ),
)

process.TFileService = cms.Service("TFileService",
				   fileName = cms.string("/tmp/rewang/csctf.root")
                                  )


#process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

#-------------------------------------
# sequences needed for L1 trigger DQM
#

# standard unpacking sequence
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")

process.RawToDigi = cms.Sequence(process.csctfDigis
                    )

process.RawToDigi.remove("siPixelDigis")
process.RawToDigi.remove("siStripDigis")
process.RawToDigi.remove("scalersRawToDigi")
process.RawToDigi.remove("castorDigis")



process.csctfDigis.producer = cms.InputTag("rawDataCollector")



process.p = cms.Path(process.RawToDigi * process.OfflineDQMCSCTF)


### local run
process.OfflineDQMCSCTF.verbose = cms.untracked.bool(True)
