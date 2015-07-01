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


    //event info
    t_->Branch("run",        	&evSummary_.run,            "run/I");
    t_->Branch("lumi",       	&evSummary_.lumi,           "lumi/I");
    t_->Branch("event",      	&evSummary_.event,          "event/I");

    //mc truth
    t_->Branch("nmcparticles",  &evSummary_.nmcparticles,   "nmcparticles/I");
    t_->Branch("mc_vx",         evSummary_.mc_vx,           "mc_vx[nmcparticles]/F");
    t_->Branch("mc_vy",         evSummary_.mc_vy,           "mc_vy[nmcparticles]/F");
    t_->Branch("mc_vz",         evSummary_.mc_vz,           "mc_vz[nmcparticles]/F");
    t_->Branch("mc_px",         evSummary_.mc_px,           "mc_px[nmcparticles]/F");
    t_->Branch("mc_py",         evSummary_.mc_py,           "mc_py[nmcparticles]/F");
    t_->Branch("mc_pz",         evSummary_.mc_pz,           "mc_pz[nmcparticles]/F");
    t_->Branch("mc_en",         evSummary_.mc_en,           "mc_en[nmcparticles]/F");
    t_->Branch("mc_id",         evSummary_.mc_id,           "mc_id[nmcparticles]/I");
    t_->Branch("mc_status",     evSummary_.mc_status,       "mc_status[nmcparticles]/I");


    t_->Branch("nlcts_m", 	&evSummary_.nlcts_m,   "nlcts_m/I");

    t_->Branch("lct_m_gblphi",         evSummary_.lct_m_gblphi,           "lct_m_gblphi[nlcts_m]/F");
    t_->Branch("lct_m_lclphi",         evSummary_.lct_m_lclphi,           "lct_m_lclphi[nlcts_m]/F");
    t_->Branch("lct_m_pkdphi",         evSummary_.lct_m_pkdphi,           "lct_m_pkdphi[nlcts_m]/F");
    t_->Branch("lct_m_gbleta",         evSummary_.lct_m_gbleta,           "lct_m_gbleta[nlcts_m]/F");
    t_->Branch("lct_m_pkdeta",         evSummary_.lct_m_pkdeta,           "lct_m_pkdeta[nlcts_m]/F");
    t_->Branch("lct_m_gblZ",         evSummary_.lct_m_gblZ,           "lct_m_gblZ[nlcts_m]/F");
    t_->Branch("lct_m_station",         evSummary_.lct_m_station,           "lct_m_station[nlcts_m]/I");
    t_->Branch("lct_m_ring",         evSummary_.lct_m_ring,           "lct_m_ring[nlcts_m]/I");
    t_->Branch("lct_m_cscid",         evSummary_.lct_m_cscid,           "lct_m_cscid[nlcts_m]/I");
    t_->Branch("lct_m_endcap",         evSummary_.lct_m_endcap,           "lct_m_endcap[nlcts_m]/I");
    t_->Branch("lct_m_sector",         evSummary_.lct_m_sector,           "lct_m_sector[nlcts_m]/I");
    t_->Branch("lct_m_subsector",         evSummary_.lct_m_subsector,           "lct_m_subsector[nlcts_m]/I");
    t_->Branch("lct_m_bptx",         evSummary_.lct_m_bptx,           "lct_m_bptx[nlcts_m]/I");
    t_->Branch("lct_m_strip",         evSummary_.lct_m_strip,           "lct_m_strip[nlcts_m]/I");
    t_->Branch("lct_m_keywire",         evSummary_.lct_m_keywire,           "lct_m_keywire[nlcts_m]/I");


    t_->Branch("nlcts_p",       &evSummary_.nlcts_p,   "nlcts_p/I");

    t_->Branch("lct_p_gblphi",         evSummary_.lct_p_gblphi,           "lct_p_gblphi[nlcts_p]/F");
    t_->Branch("lct_p_lclphi",         evSummary_.lct_p_lclphi,           "lct_p_lclphi[nlcts_p]/F");
    t_->Branch("lct_p_pkdphi",         evSummary_.lct_p_pkdphi,           "lct_p_pkdphi[nlcts_p]/F");
    t_->Branch("lct_p_gbleta",         evSummary_.lct_p_gbleta,           "lct_p_gbleta[nlcts_p]/F");
    t_->Branch("lct_p_pkdeta",         evSummary_.lct_p_pkdeta,           "lct_p_pkdeta[nlcts_p]/F");
    t_->Branch("lct_p_gblZ",         evSummary_.lct_p_gblZ,           "lct_p_gblZ[nlcts_p]/F");
    t_->Branch("lct_p_station",         evSummary_.lct_p_station,           "lct_p_station[nlcts_p]/I");
    t_->Branch("lct_p_ring",         evSummary_.lct_p_ring,           "lct_p_ring[nlcts_p]/I");
    t_->Branch("lct_p_cscid",         evSummary_.lct_p_cscid,           "lct_p_cscid[nlcts_p]/I");
    t_->Branch("lct_p_endcap",         evSummary_.lct_p_endcap,           "lct_p_endcap[nlcts_p]/I");
    t_->Branch("lct_p_sector",         evSummary_.lct_p_sector,           "lct_p_sector[nlcts_p]/I");
    t_->Branch("lct_p_subsector",         evSummary_.lct_p_subsector,           "lct_p_subsector[nlcts_p]/I");
    t_->Branch("lct_p_bptx",         evSummary_.lct_p_bptx,           "lct_p_bptx[nlcts_p]/I");
    t_->Branch("lct_p_strip",         evSummary_.lct_p_strip,           "lct_p_strip[nlcts_p]/I");
    t_->Branch("lct_p_keywire",         evSummary_.lct_p_keywire,           "lct_p_keywire[nlcts_p]/I");



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
