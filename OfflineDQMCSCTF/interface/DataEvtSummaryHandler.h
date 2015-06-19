#ifndef dataevtsummaryhandler_h
#define dataevtsummaryhandler_h

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <string.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <set>
#include <cmath>

#include "Math/LorentzVector.h"
#include "TMath.h"
#include "TVector2.h"
#include "TTree.h"
#include "TLorentzVector.h"

#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
typedef std::vector<LorentzVector> LorentzVectorCollection;

#define MAXLUTS 50

struct DataEvtSummary_t {

    Int_t run,lumi,event;

    //gen level event
    Int_t nlcts_m;
    Float_t lct_m_gblphi[MAXLUTS], lct_m_lclphi[MAXLUTS], lct_m_pkdphi[MAXLUTS];
    Float_t lct_m_gbleta[MAXLUTS], lct_m_pkdeta[MAXLUTS];
    Float_t lct_m_gblZ[MAXLUTS];
    Int_t lct_m_station[MAXLUTS], lct_m_ring[MAXLUTS], lct_m_endcap[MAXLUTS], lct_m_sector[MAXLUTS], lct_m_bptx[MAXLUTS];
    Int_t lct_m_cscid[MAXLUTS];
    Int_t lct_m_strip[MAXLUTS], lct_m_keywire[MAXLUTS];

    Int_t nlcts_p;
    Float_t lct_p_gblphi[MAXLUTS], lct_p_lclphi[MAXLUTS], lct_p_pkdphi[MAXLUTS];
    Float_t lct_p_gbleta[MAXLUTS], lct_p_pkdeta[MAXLUTS];
    Float_t lct_p_gblZ[MAXLUTS];
    Int_t lct_p_station[MAXLUTS], lct_p_ring[MAXLUTS], lct_p_endcap[MAXLUTS], lct_p_sector[MAXLUTS], lct_p_bptx[MAXLUTS];
    Int_t lct_p_cscid[MAXLUTS];
    Int_t lct_p_strip[MAXLUTS], lct_p_keywire[MAXLUTS];


};

class DataEvtSummaryHandler {
public:
    //
    DataEvtSummaryHandler();
    ~DataEvtSummaryHandler();

    //current event
    DataEvtSummary_t evSummary_;
    DataEvtSummary_t &getEvent() {
        return evSummary_;
    }

    //write mode
    bool initTree(TTree *t);
    void fillTree();

    //read mode
    //bool attachToTree(TTree *t);
    int getEntries() { return (t_ ? t_->GetEntriesFast() : 0); }
    void getEntry(int ientry) {
    	resetStruct();
    	if(t_) t_->GetEntry(ientry);
    }

    void resetStruct();

private:
    //the tree
    TTree *t_;
};

#endif
