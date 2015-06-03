#include "CSCTFUtils/OfflineDQMCSCTF/interface/DataEvtSummaryHandler.h"

using namespace std;

//
DataEvtSummaryHandler::DataEvtSummaryHandler()
{
}

//
bool DataEvtSummaryHandler::initTree(TTree *t)
{
    if(t==0) return false;
    t_ = t;


    t_->Branch("nlcts", 	&evSummary_.nlcts,   "nlcts/I");
    t_->Branch("lct_gblphi",         evSummary_.lct_gblphi,           "lct_gblphi[nlcts]/F");
    t_->Branch("lct_gbleta",         evSummary_.lct_gbleta,           "lct_gbleta[nlcts]/F");
    t_->Branch("lct_gblZ",         evSummary_.lct_gblZ,           "lct_gblZ[nlcts]/F");



    t_->Branch("lct_station",         evSummary_.lct_station,           "lct_station[nlcts]/I");
    t_->Branch("lct_ring",         evSummary_.lct_ring,           "lct_ring[nlcts]/I");
    t_->Branch("lct_endcap",         evSummary_.lct_endcap,           "lct_endcap[nlcts]/I");






    return true;
}

//
bool DataEvtSummaryHandler::attachToTree(TTree *t)
{
    if(t==0) return false;
    t_ = t;


    t_->SetBranchAddress("nlcts", 	&evSummary_.nlcts);
    t_->SetBranchAddress("lct_gblphi",           evSummary_.lct_gblphi);
    t_->SetBranchAddress("lct_gbleta",           evSummary_.lct_gbleta);
    t_->SetBranchAddress("lct_gblZ",           evSummary_.lct_gblZ);


    t_->SetBranchAddress("lct_station",           evSummary_.lct_station);
    t_->SetBranchAddress("lct_ring",           evSummary_.lct_ring);
    t_->SetBranchAddress("lct_endcap",           evSummary_.lct_endcap);




    return true;
}


//
void DataEvtSummaryHandler::resetStruct()
{
    evSummary_.nlcts=0;
}

//
void DataEvtSummaryHandler::fillTree()
{
    if(t_) t_->Fill();
}

//
DataEvtSummaryHandler::~DataEvtSummaryHandler()
{
}
