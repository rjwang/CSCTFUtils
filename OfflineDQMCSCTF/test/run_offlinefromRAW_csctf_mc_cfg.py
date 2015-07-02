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
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = cms.untracked.vstring(

'file:/afs/cern.ch/user/r/rewang/mnt/CSCTF/B0T/MC_SingleMu_FlatPt20_100_GENSIM_DIGI_L1_0.root',
'file:/afs/cern.ch/user/r/rewang/mnt/CSCTF/B0T/MC_SingleMu_FlatPt20_100_GENSIM_DIGI_L1_1.root',
'file:/afs/cern.ch/user/r/rewang/mnt/CSCTF/B0T/MC_SingleMu_FlatPt20_100_GENSIM_DIGI_L1_2.root',
'file:/afs/cern.ch/user/r/rewang/mnt/CSCTF/B0T/MC_SingleMu_FlatPt20_100_GENSIM_DIGI_L1_3.root',
'file:/afs/cern.ch/user/r/rewang/mnt/CSCTF/B0T/MC_SingleMu_FlatPt20_100_GENSIM_DIGI_L1_4.root',
'file:/afs/cern.ch/user/r/rewang/mnt/CSCTF/B0T/MC_SingleMu_FlatPt20_100_GENSIM_DIGI_L1_5.root',
'file:/afs/cern.ch/user/r/rewang/mnt/CSCTF/B0T/MC_SingleMu_FlatPt20_100_GENSIM_DIGI_L1_6.root',
'file:/afs/cern.ch/user/r/rewang/mnt/CSCTF/B0T/MC_SingleMu_FlatPt20_100_GENSIM_DIGI_L1_7.root',
'file:/afs/cern.ch/user/r/rewang/mnt/CSCTF/B0T/MC_SingleMu_FlatPt20_100_GENSIM_DIGI_L1_8.root',
'file:/afs/cern.ch/user/r/rewang/mnt/CSCTF/B0T/MC_SingleMu_FlatPt20_100_GENSIM_DIGI_L1_9.root',
'file:/afs/cern.ch/user/r/rewang/mnt/CSCTF/B0T/MC_SingleMu_FlatPt20_100_GENSIM_DIGI_L1_10.root',
'file:/afs/cern.ch/user/r/rewang/mnt/CSCTF/B0T/MC_SingleMu_FlatPt20_100_GENSIM_DIGI_L1_11.root',
'file:/afs/cern.ch/user/r/rewang/mnt/CSCTF/B0T/MC_SingleMu_FlatPt20_100_GENSIM_DIGI_L1_12.root',
'file:/afs/cern.ch/user/r/rewang/mnt/CSCTF/B0T/MC_SingleMu_FlatPt20_100_GENSIM_DIGI_L1_13.root',
'file:/afs/cern.ch/user/r/rewang/mnt/CSCTF/B0T/MC_SingleMu_FlatPt20_100_GENSIM_DIGI_L1_14.root',
'file:/afs/cern.ch/user/r/rewang/mnt/CSCTF/B0T/MC_SingleMu_FlatPt20_100_GENSIM_DIGI_L1_15.root',
'file:/afs/cern.ch/user/r/rewang/mnt/CSCTF/B0T/MC_SingleMu_FlatPt20_100_GENSIM_DIGI_L1_16.root',
'file:/afs/cern.ch/user/r/rewang/mnt/CSCTF/B0T/MC_SingleMu_FlatPt20_100_GENSIM_DIGI_L1_17.root',
'file:/afs/cern.ch/user/r/rewang/mnt/CSCTF/B0T/MC_SingleMu_FlatPt20_100_GENSIM_DIGI_L1_18.root',
'file:/afs/cern.ch/user/r/rewang/mnt/CSCTF/B0T/MC_SingleMu_FlatPt20_100_GENSIM_DIGI_L1_19.root',
'file:/afs/cern.ch/user/r/rewang/mnt/CSCTF/B0T/MC_SingleMu_FlatPt20_100_GENSIM_DIGI_L1_20.root',
'file:/afs/cern.ch/user/r/rewang/mnt/CSCTF/B0T/MC_SingleMu_FlatPt20_100_GENSIM_DIGI_L1_21.root',
'file:/afs/cern.ch/user/r/rewang/mnt/CSCTF/B0T/MC_SingleMu_FlatPt20_100_GENSIM_DIGI_L1_22.root',
'file:/afs/cern.ch/user/r/rewang/mnt/CSCTF/B0T/MC_SingleMu_FlatPt20_100_GENSIM_DIGI_L1_23.root',
'file:/afs/cern.ch/user/r/rewang/mnt/CSCTF/B0T/MC_SingleMu_FlatPt20_100_GENSIM_DIGI_L1_24.root',




    )
)

process.TFileService = cms.Service("TFileService",
#                                   fileName = cms.string("/tmp/rewang/csctf_Run246963_ME-11_match_ME+11.root")
				   fileName = cms.string("Input_MC_v1.root")
                                  )


#process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

#-------------------------------------
# sequences needed for L1 trigger DQM
#

# standard unpacking sequence
#process.load("Configuration.StandardSequences.RawToDigi_Data_cff")

#process.RawToDigi = cms.Sequence(process.csctfDigis
#                    )
#
#process.RawToDigi.remove("siPixelDigis")
#process.RawToDigi.remove("siStripDigis")
#process.RawToDigi.remove("scalersRawToDigi")
#process.RawToDigi.remove("castorDigis")



#process.csctfDigis.producer = cms.InputTag("rawDataCollector")

process.OfflineDQMCSCTF.lctProducer = cms.InputTag("simCscTriggerPrimitiveDigis")

process.p = cms.Path(process.OfflineDQMCSCTF)


### local run
#process.OfflineDQMCSCTF.verbose = cms.untracked.bool(True)
