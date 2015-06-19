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


    t_->Branch("nlcts_m", 	&evSummary_.nlcts_m,   "nlcts_m/I");

    t_->Branch("lct_m_gblphi",         evSummary_.lct_m_gblphi,           "lct_m_gblphi[nlcts_m]/F");
    t_->Branch("lct_m_gbleta",         evSummary_.lct_m_gbleta,           "lct_m_gbleta[nlcts_m]/F");
    t_->Branch("lct_m_gblZ",         evSummary_.lct_m_gblZ,           "lct_m_gblZ[nlcts_m]/F");
    t_->Branch("lct_m_station",         evSummary_.lct_m_station,           "lct_m_station[nlcts_m]/I");
    t_->Branch("lct_m_ring",         evSummary_.lct_m_ring,           "lct_m_ring[nlcts_m]/I");
    t_->Branch("lct_m_cscid",         evSummary_.lct_m_cscid,           "lct_m_cscid[nlcts_m]/I");
    t_->Branch("lct_m_endcap",         evSummary_.lct_m_endcap,           "lct_m_endcap[nlcts_m]/I");
    t_->Branch("lct_m_sector",         evSummary_.lct_m_sector,           "lct_m_sector[nlcts_m]/I");
    t_->Branch("lct_m_bptx",         evSummary_.lct_m_bptx,           "lct_m_bptx[nlcts_m]/I");

    t_->Branch("nlcts_p",       &evSummary_.nlcts_p,   "nlcts_p/I");
    t_->Branch("lct_p_gblphi",         evSummary_.lct_p_gblphi,           "lct_p_gblphi[nlcts_p]/F");
    t_->Branch("lct_p_gbleta",         evSummary_.lct_p_gbleta,           "lct_p_gbleta[nlcts_p]/F");
    t_->Branch("lct_p_gblZ",         evSummary_.lct_p_gblZ,           "lct_p_gblZ[nlcts_p]/F");
    t_->Branch("lct_p_station",         evSummary_.lct_p_station,           "lct_p_station[nlcts_p]/I");
    t_->Branch("lct_p_ring",         evSummary_.lct_p_ring,           "lct_p_ring[nlcts_p]/I");
    t_->Branch("lct_p_cscid",         evSummary_.lct_p_cscid,           "lct_p_cscid[nlcts_p]/I");
    t_->Branch("lct_p_endcap",         evSummary_.lct_p_endcap,           "lct_p_endcap[nlcts_p]/I");
    t_->Branch("lct_p_sector",         evSummary_.lct_p_sector,           "lct_p_sector[nlcts_p]/I");
    t_->Branch("lct_p_bptx",         evSummary_.lct_p_bptx,           "lct_p_bptx[nlcts_p]/I");







    return true;
}

//
bool DataEvtSummaryHandler::attachToTree(TTree *t)
{
    if(t==0) return false;
    t_ = t;


    t_->SetBranchAddress("nlcts_m", 	&evSummary_.nlcts_m);

    t_->SetBranchAddress("lct_m_gblphi",           evSummary_.lct_m_gblphi);
    t_->SetBranchAddress("lct_m_gbleta",           evSummary_.lct_m_gbleta);
    t_->SetBranchAddress("lct_m_gblZ",           evSummary_.lct_m_gblZ);
    t_->SetBranchAddress("lct_m_station",           evSummary_.lct_m_station);
    t_->SetBranchAddress("lct_m_ring",           evSummary_.lct_m_ring);
    t_->SetBranchAddress("lct_m_cscid",           evSummary_.lct_m_cscid);
    t_->SetBranchAddress("lct_m_endcap",           evSummary_.lct_m_endcap);
    t_->SetBranchAddress("lct_m_sector",           evSummary_.lct_m_sector);
    t_->SetBranchAddress("lct_m_bptx",           evSummary_.lct_m_bptx);



    t_->SetBranchAddress("nlcts_p",       &evSummary_.nlcts_p);
    t_->SetBranchAddress("lct_p_gblphi",           evSummary_.lct_p_gblphi);
    t_->SetBranchAddress("lct_p_gbleta",           evSummary_.lct_p_gbleta);
    t_->SetBranchAddress("lct_p_gblZ",           evSummary_.lct_p_gblZ);
    t_->SetBranchAddress("lct_p_station",           evSummary_.lct_p_station);
    t_->SetBranchAddress("lct_p_ring",           evSummary_.lct_p_ring);
    t_->SetBranchAddress("lct_p_cscid",           evSummary_.lct_p_cscid);
    t_->SetBranchAddress("lct_p_endcap",           evSummary_.lct_p_endcap);
    t_->SetBranchAddress("lct_p_sector",           evSummary_.lct_p_sector);
    t_->SetBranchAddress("lct_p_bptx",           evSummary_.lct_p_bptx);






    return true;
}


//
void DataEvtSummaryHandler::resetStruct()
{
    evSummary_.nlcts_m=0;
    evSummary_.nlcts_p=0;
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
