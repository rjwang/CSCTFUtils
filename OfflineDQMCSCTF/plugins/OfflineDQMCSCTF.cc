// -*- C++ -*-
//
// Package:    CSCTFUtils/OfflineDQMCSCTF
// Class:      OfflineDQMCSCTF
//
/**\class OfflineDQMCSCTF OfflineDQMCSCTF.cc CSCTFUtils/OfflineDQMCSCTF/plugins/OfflineDQMCSCTF.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Renjie Wang
//         Created:  Sun, 22 Mar 2015 18:18:25 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


// includes to fetch all reguired data products from the edm::Event
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/L1CSCTrackFinder/interface/L1CSCTrackCollection.h"
#include "DataFormats/L1CSCTrackFinder/interface/L1CSCStatusDigiCollection.h"
#include "L1Trigger/CSCCommonTrigger/interface/CSCTriggerGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "DataFormats/L1CSCTrackFinder/interface/CSCTriggerContainer.h"
#include "DataFormats/L1CSCTrackFinder/interface/TrackStub.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuRegionalCand.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTCand.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTExtendedCand.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTReadoutCollection.h"

// Sector Receiver LUT class to transform wire/strip numbers to eta/phi observables
#include "L1Trigger/CSCTrackFinder/interface/CSCSectorReceiverLUT.h"
#include "CondFormats/L1TObjects/interface/L1MuTriggerScales.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerScalesRcd.h"
#include "CondFormats/L1TObjects/interface/L1MuTriggerPtScale.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerPtScaleRcd.h"

#include "DataFormats/L1CSCTrackFinder/interface/L1CSCStatusDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/L1CSCTrackFinder/interface/L1CSCTrackCollection.h"
#include "DataFormats/L1CSCTrackFinder/interface/CSCTriggerContainer.h"
#include "DataFormats/L1CSCTrackFinder/interface/TrackStub.h"

#include "TMath.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <unistd.h>

#include "CSCTFUtils/OfflineDQMCSCTF/interface/DataEvtSummaryHandler.h"
//
// class declaration
//

class OfflineDQMCSCTF : public edm::EDAnalyzer {
public:
    explicit OfflineDQMCSCTF(const edm::ParameterSet&);
    ~OfflineDQMCSCTF();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    bool gangedME1a_;
    bool verbose_;
    CSCSectorReceiverLUT *srLUTs_[5];


    DataEvtSummaryHandler summaryHandler_;


//    edm::InputTag lctProducer;
    edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection> corrlctsToken_;

    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
OfflineDQMCSCTF::OfflineDQMCSCTF(const edm::ParameterSet& iConfig) :
//    lctProducer( iConfig.getParameter< edm::InputTag >("lctProducer") )
    corrlctsToken_(	consumes<CSCCorrelatedLCTDigiCollection>(iConfig.getParameter< edm::InputTag >("lctProducer")) )
{
    //now do what ever initialization is needed
    verbose_ = iConfig.getUntrackedParameter<bool>("verbose", false);

    //now do what ever initialization is needed
    edm::Service<TFileService> fs;
    summaryHandler_.initTree(  fs->make<TTree>("data","Event Summary") );
    TFileDirectory baseDir=fs->mkdir(iConfig.getParameter<std::string>("dtag"));


    // instantiate standard on-fly SR LUTs from CSC TF emulator package
    bzero(srLUTs_,sizeof(srLUTs_));
    int endcap=1, sector=1; // assume SR LUTs are all same for every sector in either of endcaps
    bool TMB07=true; // specific TMB firmware
    // Create a dummy pset for SR LUTs
    edm::ParameterSet srLUTset;
    srLUTset.addUntrackedParameter<bool>("ReadLUTs", false);
    srLUTset.addUntrackedParameter<bool>("Binary",   false);
    srLUTset.addUntrackedParameter<std::string>("LUTPath", "./");
    for(int station=1,fpga=0; station<=4 && fpga<5; station++) {
        if(station==1)
            for(int subSector=0; subSector<2 && fpga<5; subSector++)
                srLUTs_[fpga++] = new CSCSectorReceiverLUT(endcap, sector, subSector+1, station, srLUTset, TMB07);
        else
            srLUTs_[fpga++] = new CSCSectorReceiverLUT(endcap, sector, 0, station, srLUTset, TMB07);
    }


    gangedME1a_ = iConfig.getUntrackedParameter<bool>("gangedME1a", false);

}


OfflineDQMCSCTF::~OfflineDQMCSCTF()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    for(int i=0; i<5; i++)
        delete srLUTs_[i]; //free the array of pointers

}


//
// member functions
//

// ------------ method called for each event  ------------
void
OfflineDQMCSCTF::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;

    summaryHandler_.resetStruct();
    //event summary to be filled
    DataEvtSummary_t &ev=summaryHandler_.getEvent();

    edm::ESHandle<CSCGeometry> cscGeometry;
    iSetup.get<MuonGeometryRecord>().get(cscGeometry);


    edm::Handle<CSCCorrelatedLCTDigiCollection> corrlcts;
    iEvent.getByToken(corrlctsToken_, corrlcts);

    //std::cout << "\n ======= lctProducer ======= " << std::endl;
    std::vector<int> saveEndcaps;
    std::vector<int> saveStations;
    std::vector<int> saveRings;
    std::vector<int> saveCscId;
    std::vector<int> saveSector;

    std::vector<double> saveGblPhi;
    std::vector<double> saveGblEta;
    std::vector<double> saveGblZ;

    for(CSCCorrelatedLCTDigiCollection::DigiRangeIterator csc=corrlcts.product()->begin(); csc!=corrlcts.product()->end(); csc++) {
        CSCCorrelatedLCTDigiCollection::Range range1 = corrlcts.product()->get((*csc).first);
        for(CSCCorrelatedLCTDigiCollection::const_iterator lct=range1.first; lct!=range1.second; lct++) {
            int endcap    = (*csc).first.endcap()-1; //0: +Z, 1: -Z
            int station   = (*csc).first.station()-1; //0,1,2,3
            int sector    = (*csc).first.triggerSector()-1; //0,1,2,3,4,5
            int subSector = CSCTriggerNumbering::triggerSubSectorFromLabels((*csc).first); //0,1,2
            int ring      = (*csc).first.ring(); // 1,2,3
            int cscId     = (*csc).first.triggerCscId()-1;
            int fpga      = ( subSector ? subSector-1 : station+1 );
            int strip     = lct -> getStrip();
            int keyWire   = lct -> getKeyWG();
            int bx        = lct -> getBX();
            int quality   = lct -> getQuality();
            int bend      = lct -> getBend();

            int endcapAssignment = 1;
            int shift = 1;
            float sectorArg = sector;
            //float sectorArg = j;

            if( endcap == 1 ) {
                endcapAssignment = -1;
                shift = 2;
                //sectorArg = sector - 6;
            }

            int signedStation = (station + shift)* endcapAssignment;
            if( (station == 0) && (endcap == 0)) signedStation = subSector - 1;
            if( (station == 0) && (endcap == 1)) signedStation = (-1)*subSector;

            float chamberArg1 = cscId * 0.1 + sectorArg;
            //float chamberArg1 = i*0.1 + sectorArg;
            //std::cout << "First" << i << " " << sectorArg << " " << chamberArg1 << std::endl;

            float chamberArg11 = chamberArg1;
            if(sectorArg == 1) chamberArg1 = chamberArg11 - 0.1;
            if(sectorArg == 2) chamberArg1 = chamberArg11 - 0.2;
            if(sectorArg == 3) chamberArg1 = chamberArg11 - 0.3;
            if(sectorArg == 4) chamberArg1 = chamberArg11 - 0.4;
            if(sectorArg == 5) chamberArg1 = chamberArg11 - 0.5;

            if(verbose_) std::cout << "cscId, station, sector, endcap, sectorArg, chamber Arg: " << cscId << ", " << station << ", " <<sector << ", " << endcap << ", " << chamberArg1 << ", " << signedStation << std::endl;

            //csctfChamberOccupancies->Fill(chamberArg1, signedStation);
            //int bunchX = ( (lct->getBX()) - 6 );

            //int timingSectorArg = 3*(sector) + (lct->getMPCLink());
            //if( endcap == 1) timingSectorArg = 3*(sector + 6) + (lct->getMPCLink());
            //std::cout << "Sector, MPCLink, TSA, endcap: " << sector << ", " << lct->getMPCLink() << ", " << timingSectorArg << ", " << endcap << std::endl;

            //csctfbx->Fill(timingSectorArg, bunchX );
            //std::cout << "LCT'S, encap: " << endcap << ", station: " << station << ", sector: " << sector << ", subSector: " << subSector << ", cscId: " << cscId << std:: endl;
            //End JAG

            // Check if Det Id is within pysical range:
            if( endcap<0||endcap>1 || sector<0||sector>6 || station<0||station>3 || cscId<0||cscId>8 || fpga<0||fpga>4) {
                cout << "L1CSCTF: CSC TP are out of range: " <<"  endcap: "<<(endcap+1)<<"  station: "<<(station+1) <<"  sector: "<<(sector+1)<<"  subSector: "<<subSector<<"  fpga: "<<fpga<<"  cscId: "<<(cscId+1) << endl;
                edm::LogError("L1CSCTF: CSC TP are out of range: ")<<"  endcap: "<<(endcap+1)<<"  station: "<<(station+1) <<"  sector: "<<(sector+1)<<"  subSector: "<<subSector<<"  fpga: "<<fpga<<"  cscId: "<<(cscId+1);
                continue;
            }


            //unsigned short int testlocalPhi=0;






            lclphidat lclPhi;
            try {
                lclPhi = srLUTs_[fpga]->localPhi(lct->getStrip(), lct->getPattern(), lct->getQuality(), lct->getBend(), gangedME1a_);
            } catch(cms::Exception &) {
                bzero(&lclPhi,sizeof(lclPhi));
            }

            gblphidat gblPhi;
            try {
                gblPhi = srLUTs_[fpga]->globalPhiME(lclPhi.phi_local, lct->getKeyWG(), cscId+1, gangedME1a_);
            } catch(cms::Exception &) {
                bzero(&gblPhi,sizeof(gblPhi));
            }

            gbletadat gblEta;
            try {
                gblEta = srLUTs_[fpga]->globalEtaME(lclPhi.phi_bend_local, lclPhi.phi_local, lct->getKeyWG(), cscId+1, gangedME1a_);
            } catch(cms::Exception &) {
                bzero(&gblEta,sizeof(gblEta));
            }


            //TrackStub
            csctf::TrackStub theStub((*lct), (*csc).first);
            theStub.setPhiPacked(gblPhi.global_phi);
            theStub.setEtaPacked(gblEta.global_eta);



            // SR LUT gives packed eta and phi values -> normilize them to 1 by scale them to 'max' and shift by 'min'
            //float etaG = gblEta.global_eta/127. * 1.5 + 0.9;
            //float etaG = gblEta.global_eta/128. * 1.6 + 0.9;
            float etaG = theStub.etaValue();
            //float phi_temp = (gblPhi.global_phi + ( sector + (endcap?0:6) )*4096 + station*4096*12) * 1./(4*4096*12);
            //float phiG = fmod( phi_temp + (( sector + (endcap?0:6) )*M_PI*60/180.) + (M_PI*14/180.), 2.*M_PI);

            //theStub.phiValue() = theStub.phiPacked()*62*M_PI/(180.*4096.)
            float phiG = fmod( theStub.phiValue()+15.0*M_PI/180+(sector)*60.0*M_PI/180, 2.*M_PI );
            //csctfoccupancies->Fill( gblEta.global_eta/127. * 1.5 + 0.9, (gblPhi.global_phi + ( sector + (endcap?0:6) )*4096 + station*4096*12) * 1./(4*4096*12) );

            //cout << gblPhi.global_phi * thePhiBinning;

            //cout << "thePhiBinning: " << CSCTFConstants::SECTOR_RAD/(1<<CSCBitWidths::kGlobalPhiDataBitWidth) << endl;
            //cout << "test: " << theStub.phiPacked()*62*M_PI/(180.*4096.) << endl;
            //thePhiBinning = CSCTFConstants::SECTOR_RAD/(1<<CSCBitWidths::kGlobalPhiDataBitWidth);
            //SECTOR_RAD = 62*M_PI/180.
            //SECTOR_RAD = CSCTFConstants::SECTOR_DEG*CSCTFConstants::RAD_PER_DEGREE;


            if (verbose_) {
                cout << "\n ===== Current LCT Values ====== \n";
                cout << " Endcap    = " << endcap           << endl;
                cout << " Station   = " << station          << endl;
                cout << " Sector    = " << sector           << endl;
                cout << " SubSector = " << subSector                << endl;
                cout << " Ring      = " << ring             << endl;
                cout << " CscId     = " << cscId            << endl;
                cout << " Fpga      = " << fpga             << endl;
                cout << " Strip     = " << strip            << endl;
                cout << " Bx        = " << bx                       << " -----> " << theStub.BX()           << endl;
                cout << " Quality   = " << quality          << endl;
                cout << " Bend      = " << bend             << endl;
                cout << " KeyWire   = " << keyWire          << endl;
                cout << " Phi Bit   = " << gblPhi.global_phi        << " -----> " << theStub.phiPacked()    << endl;
                cout << " Eta Bit   = " << gblEta.global_eta        << " -----> " << theStub.etaPacked()    << endl;
                cout << " Phi Gbl   = " << phiG             << endl;
                cout << " Eta Gbl   = " << etaG             << " -----> " << theStub.etaValue()     << endl;
                cout << " Phi Loc   = " << lclPhi.phi_local         << " -----> " << theStub.phiValue() << " ---> " << theStub.phiPacked()*62*M_PI/(180.*4096.) << endl;
            }


            // quality cuts
            if(quality<10) continue;


            CSCDetId _cscid_(endcap+1, station+1, ring, sector+1);
            double globalZ(0);
            globalZ = cscGeometry->idToDet(_cscid_)->position().z();

            saveEndcaps.push_back(endcap);
            saveStations.push_back(station);
            saveRings.push_back(ring);
	    saveCscId.push_back(cscId);
	    saveSector.push_back(sector);

            //saveGblPhi.push_back(phiG);
            saveGblPhi.push_back(TMath::ACos(fabs(cos(phiG))) ); // for resolution study
            saveGblEta.push_back(etaG);
            saveGblZ.push_back(globalZ);



            //double xx = 1.*TMath::Cos(phiG);
            //double yy = 1.*TMath::Sin(phiG);
            //double zz = 1.*sinh(etaG);


            /////////////////////////////////////////

            /*
                        std::cout << " endcap: " << endcap
                                  << " station: " << station
                                  << " ring: " << ring
                                  << " CscId: " << cscId
                                  << " sector: " << sector
                                  << " subSector: " << subSector
                                  << " strip: " << strip
                                  << " WG: " << keyWire
                                  << " GlobalPhi: " << phiG
                                  << " GlobalEta: " << etaG
                                  << " xx: " << xx
                                  << " yy: " << yy
                                  << " zz: " << zz
                                  << std::endl;
            */

        } // lct != range1.scond
    } // csc!=corrlcts.product()->end()


    //resolution of LUTs
    ev.nlcts = 0;
    if(saveGblPhi.size()>4) {

        bool hasSameSector(true);
        for(size_t i=0; i<saveSector.size(); i++) {
            if(i==saveSector.size()-1) break;
            hasSameSector &= (saveSector[i]==saveSector[i+1]);
        }

        if(hasSameSector) {

            std::cout << "\n ======= LUT Resolution ======= " << std::endl;
            for(size_t j=0; j<saveGblPhi.size(); j++) {
		std::cout << " Endcap: "    << saveEndcaps[j]
			  << " Station: "   << saveStations[j]
			  << " Ring: " 	    << saveRings[j]
			  << " CSCID: "     << saveCscId[j]
			  << " Sector: "    << saveSector[j]
                	  << " GlobalPhi: " << saveGblPhi[j]
                          << " GlobalEta: " << saveGblEta[j]
			  << " GlobalZ: "   << saveGblZ[j]
                          << std::endl;

                ev.lct_gblphi[ev.nlcts] = saveGblPhi[j];
                ev.lct_gbleta[ev.nlcts] = saveGblEta[j];
		ev.lct_gblZ[ev.nlcts] = saveGblZ[j];
                ev.lct_endcap[ev.nlcts] = saveEndcaps[j];
                ev.lct_station[ev.nlcts]= saveStations[j];
                ev.lct_ring[ev.nlcts]   = saveRings[j];

                ev.nlcts++;
            }
        }


    }


    summaryHandler_.fillTree();

}


// ------------ method called once each job just before starting event loop  ------------
void
OfflineDQMCSCTF::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
OfflineDQMCSCTF::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
/*
void
OfflineDQMCSCTF::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
OfflineDQMCSCTF::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
OfflineDQMCSCTF::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
OfflineDQMCSCTF::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
OfflineDQMCSCTF::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(OfflineDQMCSCTF);
