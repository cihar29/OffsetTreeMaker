# PYTHON configuration file for class: OffsetTreeMaker
# Author: C. Harrington
# Date:  19 - January - 2015

import FWCore.ParameterSet.Config as cms

process = cms.Process("Ana")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

readFiles = cms.untracked.vstring()
process.source = cms.Source ("PoolSource", fileNames = readFiles)
readFiles.extend( [

  '/store/relval/CMSSW_8_0_0/ZeroBias/RECO/80X_dataRun2_relval_v0_RelVal_zb2015D-v1/10000/0041FE0D-46DA-E511-B85B-0025905A6122.root'

  # '/store/data/Run2015D/ZeroBias/AOD/16Dec2015-v1/110000/EA449360-1EAF-E511-8D54-00266CFAE748.root'

  # '/store/mc/RunIIFall15DR76/SingleNeutrino/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/004252D9-C094-E511-928F-AC162DA8C2B0.root'

] );

isMC = cms.bool(False)

if isMC:
  OutputName = "_MC"
  era = "Fall15_25nsV2_MC"
  jecLevels = "ak4PFCHSL1FastL2L3Corrector"

  #process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff" )
  #from Configuration.AlCa.GlobalTag import GlobalTag
  #process.GlobalTag = GlobalTag( process.GlobalTag, '76X_mcRun2_asymptotic_v12' )

else:
  OutputName = "_Data"
  era = "Fall15_25nsV2_DATA"
  jecLevels = "ak4PFCHSL1FastL2L3ResidualCorrector"

  process.load( "Configuration.Geometry.GeometryIdeal_cff" )
  process.load( "Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff" )
  process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff" )
  from Configuration.AlCa.GlobalTag import GlobalTag
  process.GlobalTag = GlobalTag( process.GlobalTag, '80X_dataRun2_Prompt_v1' )

  # ZeroBias Trigger
  process.HLTZeroBias =cms.EDFilter("HLTHighLevel",
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    HLTPaths = cms.vstring('HLT_ZeroBias_*'),
    eventSetupPathsKey = cms.string(''),
    andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw = cms.bool(False)
  )

  #Beam Halo
  process.load('RecoMET.METFilters.CSCTightHaloFilter_cfi')

  #HCAL HBHE
  process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
  process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
  process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
    inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResultRun2Tight'),
    reverseDecision = cms.bool(False)
  )

  #Bad EE Supercrystal filter
  #process.load('RecoMET.METFilters.eeBadScFilter_cfi')

#Jet Corrections
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
    connect = cms.string('sqlite_file:'+era+'.db'),
    toGet =  cms.VPSet(
        cms.PSet(
            record = cms.string("JetCorrectionsRecord"),
            tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PFchs"),
            label= cms.untracked.string("AK4PFchs")
        )
    )
)
process.es_prefer_jec = cms.ESPrefer("PoolDBESSource","jec")

process.load('JetMETCorrections.Configuration.JetCorrectors_cff')

process.ak4PFJetsCHSl1l2l3 = cms.EDProducer('CorrectedPFJetProducer',
    src         = cms.InputTag("ak4PFJetsCHS"),
    correctors  = cms.VInputTag(jecLevels)
)

process.pf = cms.EDAnalyzer("OffsetTreeMaker",
    numSkip = cms.int32(101),
    RootFileName = cms.string("Offset" + OutputName + ".root"),
    isMC = isMC,
    reweight = cms.bool(False),
    trackTag = cms.InputTag("generalTracks"),
    pfTag = cms.InputTag("particleFlow"),
    pvTag = cms.InputTag("offlinePrimaryVertices"),
    muTag = cms.InputTag("addPileupInfo"),
    rhoTag = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoC0Tag = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
    rhoCCTag = cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp"),
    pfJetTag = cms.InputTag("ak4PFJetsCHS"),
    corrPfJetTag = cms.InputTag("ak4PFJetsCHSl1l2l3")
)

process.myseq = cms.Sequence( process.ak4PFJetsCHSl1l2l3 * process.pf )

if isMC :
  process.p = cms.Path( process.myseq )
else:
  process.p = cms.Path( process.HLTZeroBias * 
                        process.CSCTightHaloFilter *
                        process.HBHENoiseFilterResultProducer *
                        process.ApplyBaselineHBHENoiseFilter *
                        #process.eeBadScFilter *
                        process.myseq )
