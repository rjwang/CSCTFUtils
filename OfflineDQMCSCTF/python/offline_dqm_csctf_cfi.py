import FWCore.ParameterSet.Config as cms

process = cms.Process("OfflineDQMCSCTF")

process.OfflineDQMCSCTF = cms.EDAnalyzer('OfflineDQMCSCTF',
	lctProducer = cms.InputTag("csctfDigis"),
	dtag = cms.string('csctf'),
	GenParticles = cms.untracked.InputTag("genParticles"),
	verbose = cms.untracked.bool(False),
	gangedME1a = cms.untracked.bool(False)
)
