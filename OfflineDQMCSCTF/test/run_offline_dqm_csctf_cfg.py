import FWCore.ParameterSet.Config as cms


from CSCTFUtils.OfflineDQMCSCTF.offline_dqm_csctf_cfi import *


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
                                        SkipEvent = cms.untracked.vstring('ProductNotFound')
                                        )


##testing
globalTagValue = 'GR_H_V26'
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.connect = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'
process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
process.GlobalTag.globaltag = globalTagValue+'::All'
es_prefer_GlobalTag = cms.ESPrefer('GlobalTag')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/user/rewang/csctf/local_run_226129_csctf.root'
	#'/store/data/Commissioning2015/Cosmics/RAW/v1/000/237/956/00000/'
	#'/store/data/Commissioning2015/Cosmics/RAW/v1/000/237/956/00000/005F1842-E4CB-E411-BA7B-02163E011D6E.root',
	#'/store/data/Commissioning2015/Cosmics/RAW/v1/000/237/956/00000/00928605-DBCA-E411-A754-02163E0122AB.root',
	#'/store/data/Commissioning2015/Cosmics/RAW/v1/000/237/956/00000/00EB7808-E1CB-E411-A176-02163E01270C.root',
	#'/store/data/Commissioning2015/Cosmics/RAW/v1/000/237/956/00000/0266BBD2-1BCB-E411-BDA3-02163E0123B0.root',
	#'/store/data/Commissioning2015/Cosmics/RAW/v1/000/237/956/00000/02804D14-E8CA-E411-8DC9-02163E01257A.root',
	#'/store/data/Commissioning2015/Cosmics/RAW/v1/000/237/956/00000/02DEEA48-E4CB-E411-BE14-02163E012059.root',
	#'/store/data/Commissioning2015/Cosmics/RAW/v1/000/237/956/00000/04118245-E4CB-E411-A91F-02163E01283E.root',
	#'/store/data/Commissioning2015/Cosmics/RAW/v1/000/237/956/00000/0607F057-E4CB-E411-A4BC-02163E011CDD.root',
	#'/store/data/Commissioning2015/Cosmics/RAW/v1/000/237/956/00000/06220D5A-00CB-E411-8F09-02163E012785.root',
	#'/store/data/Commissioning2015/Cosmics/RAW/v1/000/237/956/00000/08699835-F6CA-E411-B762-02163E011801.root',
	#'/store/data/Commissioning2015/Cosmics/RAW/v1/000/237/956/00000/0A19AAF0-EFCA-E411-AF04-02163E012379.root',
	'/store/data/Run2015A/Commissioning/RAW/v1/000/246/923/00000/0EDA384D-EC09-E511-95CE-02163E0143E9.root',
    )
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("csctf_Run246923.root")
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
