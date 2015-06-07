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

	#'/store/data/Run2015A/Commissioning/RAW/v1/000/246/923/00000/0E4F1045-070A-E511-BEBE-02163E011ABC.root',
	#'/store/data/Run2015A/Commissioning/RAW/v1/000/246/923/00000/0EDA384D-EC09-E511-95CE-02163E0143E9.root',
	#'/store/data/Run2015A/Commissioning/RAW/v1/000/246/923/00000/7C9D5DEF-070A-E511-AEC9-02163E012925.root',
	#'/store/data/Run2015A/Commissioning/RAW/v1/000/246/923/00000/A66CD04C-070A-E511-A37A-02163E011ACE.root',
	#'/store/data/Run2015A/Commissioning/RAW/v1/000/246/923/00000/BAA28E49-070A-E511-93A4-02163E0142B3.root',

	#'/store/data/Run2015A/Commissioning/RAW/v1/000/246/926/00000/6A496D1A-F709-E511-B2B9-02163E0145D2.root',
	#'/store/data/Run2015A/Commissioning/RAW/v1/000/246/926/00000/C44F8E28-F909-E511-940D-02163E014295.root',

	#'/store/data/Run2015A/Commissioning/RAW/v1/000/246/930/00000/0C3D5BEF-F809-E511-83F9-02163E01462A.root',
	#'/store/data/Run2015A/Commissioning/RAW/v1/000/246/930/00000/6AA54584-F909-E511-8DEF-02163E011DE4.root',

	'/store/data/Run2015A/Commissioning/RAW/v1/000/246/963/00000/DA999EA5-5B0A-E511-A03E-02163E0145FF.root',
	'/store/data/Run2015A/Commissioning/RAW/v1/000/246/963/00000/DCB89DA3-5B0A-E511-8A63-02163E014192.root',
	'/store/data/Run2015A/Commissioning/RAW/v1/000/246/963/00000/DCF64CA7-5B0A-E511-8B0E-02163E0141F1.root',
	'/store/data/Run2015A/Commissioning/RAW/v1/000/246/963/00000/E4BCCC1C-660A-E511-BEBD-02163E011B7C.root',
	'/store/data/Run2015A/Commissioning/RAW/v1/000/246/963/00000/E63ADDA5-5B0A-E511-8059-02163E014649.root',
	'/store/data/Run2015A/Commissioning/RAW/v1/000/246/963/00000/EA3BD4A4-5B0A-E511-B256-02163E011A81.root',
	'/store/data/Run2015A/Commissioning/RAW/v1/000/246/963/00000/EAF602E1-800A-E511-A457-02163E0137FE.root',
	'/store/data/Run2015A/Commissioning/RAW/v1/000/246/963/00000/EC43DFA8-5B0A-E511-ACBD-02163E0119DC.root',
	'/store/data/Run2015A/Commissioning/RAW/v1/000/246/963/00000/ECAE73A4-5B0A-E511-9487-02163E0141D2.root',
	'/store/data/Run2015A/Commissioning/RAW/v1/000/246/963/00000/ECE490A6-5B0A-E511-97FE-02163E0138BC.root',
	'/store/data/Run2015A/Commissioning/RAW/v1/000/246/963/00000/EED08BA4-5B0A-E511-8304-02163E0146EE.root',
	'/store/data/Run2015A/Commissioning/RAW/v1/000/246/963/00000/EEF4FF1F-660A-E511-9162-02163E012BAA.root',
	'/store/data/Run2015A/Commissioning/RAW/v1/000/246/963/00000/F26B03DC-530A-E511-A2CC-02163E011D8F.root',
	'/store/data/Run2015A/Commissioning/RAW/v1/000/246/963/00000/F6AD2709-5D0A-E511-A268-02163E014281.root',
	'/store/data/Run2015A/Commissioning/RAW/v1/000/246/963/00000/F6E7A2A5-5B0A-E511-A109-02163E01427A.root',
	'/store/data/Run2015A/Commissioning/RAW/v1/000/246/963/00000/FAE8B7A4-5B0A-E511-AFAA-02163E013775.root',
	'/store/data/Run2015A/Commissioning/RAW/v1/000/246/963/00000/FC19C2AD-5B0A-E511-B4D3-02163E0139DE.root',

    )
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("/tmp/rewang/csctf_Run246963.root")
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
