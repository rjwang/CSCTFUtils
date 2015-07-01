import FWCore.ParameterSet.Config as cms


from CSCTFUtils.OfflineDQMCSCTF.offline_dqm_csctf_cfi import *


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
                                        SkipEvent = cms.untracked.vstring('ProductNotFound')
                                        )


##testing
globalTagValue = 'GR_H_V58'
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.connect = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'
process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
process.GlobalTag.globaltag = globalTagValue+'::All'
es_prefer_GlobalTag = cms.ESPrefer('GlobalTag')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
 '/store/data/Run2015A/SingleMu/RAW/v1/000/246/963/00000/06A0F6A7-5B0A-E511-BEDE-02163E0141FA.root',
# '/store/data/Run2015A/SingleMu/RAW/v1/000/246/963/00000/0C948DA8-5B0A-E511-B90C-02163E013979.root',
# '/store/data/Run2015A/SingleMu/RAW/v1/000/246/963/00000/28B6F4A8-5B0A-E511-80FC-02163E014614.root',
# '/store/data/Run2015A/SingleMu/RAW/v1/000/246/963/00000/34E885B1-5B0A-E511-B87A-02163E01478F.root',
# '/store/data/Run2015A/SingleMu/RAW/v1/000/246/963/00000/38C5F5E0-800A-E511-9C3A-02163E013675.root',
# '/store/data/Run2015A/SingleMu/RAW/v1/000/246/963/00000/3A9F3BE2-800A-E511-9B71-02163E011B58.root',
# '/store/data/Run2015A/SingleMu/RAW/v1/000/246/963/00000/3ADEACA7-5B0A-E511-B61B-02163E0138EB.root',
# '/store/data/Run2015A/SingleMu/RAW/v1/000/246/963/00000/3CB7831D-660A-E511-B891-02163E01387E.root',
# '/store/data/Run2015A/SingleMu/RAW/v1/000/246/963/00000/4EBA06B8-5B0A-E511-910C-02163E014129.root',
# '/store/data/Run2015A/SingleMu/RAW/v1/000/246/963/00000/5024B5E1-800A-E511-B833-02163E011850.root',
# '/store/data/Run2015A/SingleMu/RAW/v1/000/246/963/00000/50AEBBA8-5B0A-E511-989A-02163E0142BB.root',
# '/store/data/Run2015A/SingleMu/RAW/v1/000/246/963/00000/5A3E22A8-5B0A-E511-AD7E-02163E014718.root',
# '/store/data/Run2015A/SingleMu/RAW/v1/000/246/963/00000/6AC8CA1E-660A-E511-8CA1-02163E012123.root',
# '/store/data/Run2015A/SingleMu/RAW/v1/000/246/963/00000/74706CA8-5B0A-E511-BB45-02163E013776.root',
# '/store/data/Run2015A/SingleMu/RAW/v1/000/246/963/00000/7A6FC4E4-800A-E511-A701-02163E0138B3.root',
# '/store/data/Run2015A/SingleMu/RAW/v1/000/246/963/00000/8000E8E3-800A-E511-9654-02163E014698.root',
# '/store/data/Run2015A/SingleMu/RAW/v1/000/246/963/00000/984367AC-5B0A-E511-AE0E-02163E0135BC.root',
# '/store/data/Run2015A/SingleMu/RAW/v1/000/246/963/00000/AA60BC45-660A-E511-9ED0-02163E0137E1.root',
# '/store/data/Run2015A/SingleMu/RAW/v1/000/246/963/00000/B014A5A8-5B0A-E511-80DD-02163E014402.root',
# '/store/data/Run2015A/SingleMu/RAW/v1/000/246/963/00000/B29C0F29-550A-E511-BA6E-02163E014503.root',
# '/store/data/Run2015A/SingleMu/RAW/v1/000/246/963/00000/C473BEA8-5B0A-E511-91BB-02163E0144F1.root',
# '/store/data/Run2015A/SingleMu/RAW/v1/000/246/963/00000/C67C0DE3-800A-E511-A3D7-02163E012462.root',
# '/store/data/Run2015A/SingleMu/RAW/v1/000/246/963/00000/CCFF57E2-800A-E511-A2B2-02163E01365F.root',
# '/store/data/Run2015A/SingleMu/RAW/v1/000/246/963/00000/D01A70A8-5B0A-E511-B1D1-02163E0142F0.root',
# '/store/data/Run2015A/SingleMu/RAW/v1/000/246/963/00000/F4CA8DA7-5B0A-E511-B47F-02163E014175.root',

    )
)

process.TFileService = cms.Service("TFileService",
#                                   fileName = cms.string("/tmp/rewang/csctf_Run246963_ME-11_match_ME+11.root")
				   fileName = cms.string("/tmp/rewang/csctf_Run246963_nochange.root")
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
#
process.RawToDigi.remove("siPixelDigis")
process.RawToDigi.remove("siStripDigis")
process.RawToDigi.remove("scalersRawToDigi")
process.RawToDigi.remove("castorDigis")



process.csctfDigis.producer = cms.InputTag("rawDataCollector")



process.p = cms.Path(process.RawToDigi * process.OfflineDQMCSCTF)


### local run
process.OfflineDQMCSCTF.verbose = cms.untracked.bool(True)
