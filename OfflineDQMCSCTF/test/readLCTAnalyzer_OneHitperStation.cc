#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TH1.h>
#include <TVirtualFitter.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>

#include "TFile.h"
#include "TTree.h"
#include "TObject.h"
#include <utility>      // std::pair
#include <iostream>     // std::cout

#include <vector>
#include <cmath>
using namespace ROOT::Math;
using namespace std;

#define MAXLUTS 50

struct CSCTF_t {
    int status;
    float dphi,deta;
    float rawphi,raweta,rawZ;
    int endcap,station,ring,sector,cscid;
    float calphi,caleta;
    float p0,p1,p2,p3,distance;
    float dr0;
    double chi2;
};

double distance2(float x,float y,float z, double *par);

double deltaPhi(double phi1, double phi2)
{
    double result = phi1 - phi2;
    while (result > M_PI) result -= 2*M_PI;
    while (result <= -M_PI) result += 2*M_PI;
    return result;
}

//double
CSCTF_t
getDeltaPhiEta(double *par, double phi, double eta, double z)
{
    CSCTF_t myfit;

    // delta phi
    double x = par[0] + par[1]*z;
    double y = par[2] + par[3]*z;

    double cal_phi = TMath::ATan( fabs(y/x) );
    if(x<0 && y>0) 	cal_phi = M_PI-cal_phi;
    else if(x<0 && y<0) cal_phi = M_PI+cal_phi;
    else if(x>0 && y<0) cal_phi = 2*M_PI - cal_phi;

    double det_phi = deltaPhi(cal_phi,phi);


    // delta Eta
    double r = sqrt(x*x+y*y);
    double cal_theta = TMath::ATan( fabs(r/z) );
    double cal_eta = (-1.)*log( TMath::Tan( cal_theta/2. ) );

    if(z<0) cal_eta = -1.*fabs(cal_eta);
    else cal_eta = fabs(cal_eta);

    double det_eta = cal_eta-eta;

    myfit.deta=det_eta;
    myfit.caleta=cal_eta;

    myfit.dphi=det_phi;
    myfit.calphi=cal_phi;

    return myfit;
}

double getDistance(double *par, double phi, double eta, double z)
{
    double theta = 2*TMath::ATan(exp(-1.*eta));
    double r = fabs(z*TMath::Tan(theta));
    double y = r*TMath::Sin(phi);
    double x = r*TMath::Cos(phi);

    if(phi>=0                  && phi<M_PI/2. )       {
        x=fabs(x);
        y=fabs(y);
    } else if(phi>=M_PI/2.       && phi<M_PI    )       {
        x=-1.*fabs(x);
        y=fabs(y);
    } else if(phi>=M_PI          && phi<3.*M_PI/2.)     {
        x=-1.*fabs(x);
        y=-1.*fabs(y);
    } else if(phi>=3.*M_PI/2.    && phi<2.*M_PI )       {
        x=fabs(x);
        y=-1.*fabs(y);
    }

    double distance = distance2(x,y,z,par);
    return sqrt(distance);
}

// define the parameteric line equation
void line(double t, double *p, double &x, double &y, double &z)
{
    // a parameteric line is define from 6 parameters but 4 are independent
    // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
    // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1;
    x = p[0] + p[1]*t;
    y = p[2] + p[3]*t;
    z = t;
}

// calculate distance line-point
double distance2(float x,float y,float z, double *par)
{
    // distance line point is D= | (xp-x0) cross  ux |
    // where ux is direction of line and x0 is a point in the line (like t = 0)
    XYZVector xp(x,y,z);
    XYZVector x0(par[0], par[2], 0. );
    XYZVector x1(par[0] + par[1], par[2] + par[3], 1. );
    XYZVector u = (x1-x0).Unit();
    double d2 = ((xp-x0).Cross(u)).Mag2();
    return d2;
}



// function to be minimized
void SumDistance2(int &, double *, double & sum, double * par, int )
{
    // the TGraph must be a global variable
    TGraph2D * gr = dynamic_cast<TGraph2D*>( (TVirtualFitter::GetFitter())->GetObjectFit() );
    assert(gr != 0);
    double * x = gr->GetX();
    double * y = gr->GetY();
    double * z = gr->GetZ();
    int npoints = gr->GetN();
    sum = 0;
    for (int i  = 0; i < npoints; ++i) {
        double d = distance2(x[i],y[i],z[i],par);
        sum += d;
    }
}

std::vector<CSCTF_t>
dofit(std::vector<double> allphiG, std::vector<double> alletaG, std::vector<double> allzG,
      std::vector<int> allendcap, std::vector<int> allsector, std::vector<int> allstation, std::vector<int> allring, std::vector<int> allcscid)
{
    std::vector<CSCTF_t> LUTrels;
    for(size_t lct=0; lct<allphiG.size(); lct++) {

        std::vector<double> tofit_phiG;
        std::vector<double> tofit_etaG;
        std::vector<double> tofit_zG;

        for(size_t t=0; t<allphiG.size(); t++) {
            if(t==lct) continue;
	    // remove the segments with same station of the target hit
	    if(allstation[t] == allstation[lct]) {
		cout << "--> remove same stations" << endl;
		continue; // remove the case has same station
	    }

            tofit_phiG.push_back(allphiG[t]);
            tofit_etaG.push_back(alletaG[t]);
            tofit_zG.push_back(allzG[t]);
        }

        TGraph2D * gr = new TGraph2D();

        for(size_t n=0; n<tofit_phiG.size(); n++) {
            double phiG = tofit_phiG[n];
            double etaG = tofit_etaG[n];
            double zG = tofit_zG[n];

            double theta = 2.*TMath::ATan(exp(-1.*etaG));

            double z = zG;
	    double r = fabs(z*TMath::Tan(theta));
	    double y = r*TMath::Sin(phiG);
	    double x = r*TMath::Cos(phiG);

            if(phiG>=0 			&& phiG<M_PI/2.	) 	{
                x=fabs(x);
                y=fabs(y);
            } else if(phiG>=M_PI/2. 	&& phiG<M_PI	) 	{
                x=-1.*fabs(x);
                y=fabs(y);
            } else if(phiG>=M_PI 	&& phiG<3.*M_PI/2.) 	{
                x=-1.*fabs(x);
                y=-1.*fabs(y);
            } else if(phiG>=3.*M_PI/2. 	&& phiG<2.*M_PI	) 	{
                x=fabs(x);
                y=-1.*fabs(y);
            }

            gr->SetPoint(n,x,y,z);
        }

        TVirtualFitter *min = TVirtualFitter::Fitter(0,4);
        min->SetObjectFit(gr);
        min->SetFCN(SumDistance2);

        Double_t arglist[10];
        arglist[0] = 0;//3;
        min->ExecuteCommand("SET PRINT",arglist,1);

        double pStart[4] = {1,1,1,1};
        min->SetParameter(0,"x0",pStart[0],0.01,0,0);
        min->SetParameter(1,"Ax",pStart[1],0.01,0,0);
        min->SetParameter(2,"y0",pStart[2],0.01,0,0);
        min->SetParameter(3,"Ay",pStart[3],0.01,0,0);

        arglist[0] = 2000; // number of function calls
        arglist[1] = 0.001; // tolerance
        int status = min->ExecuteCommand("MIGRAD",arglist,2);

        //if (minos) min->ExecuteCommand("MINOS",arglist,0);
        int nvpar,nparx;
        double amin,edm, errdef;
        min->GetStats(amin,edm,errdef,nvpar,nparx);
        min->PrintResults(1,amin);

	cout << "amin: " << amin << endl;

        // get fit parameters
        double parFit[4];
        for (int j = 0; j < 4; ++j) {
            parFit[j] = min->GetParameter(j);
            //cout  << "par: " << parFit[j] << endl;
        }

        CSCTF_t DeltaPhiEta = getDeltaPhiEta( parFit,allphiG[lct],alletaG[lct],allzG[lct] );
        double distance = getDistance( parFit,allphiG[lct],alletaG[lct],allzG[lct] );
	float dr0 = distance2(0.,0.,0.,parFit);

        CSCTF_t myfit;
        myfit.status = status;
        myfit.dphi = DeltaPhiEta.dphi;
        myfit.deta = DeltaPhiEta.deta;
        myfit.rawphi = allphiG[lct];
        myfit.raweta = alletaG[lct];
        myfit.endcap = allendcap[lct];
	myfit.sector = allsector[lct];
        myfit.station = allstation[lct];
        myfit.ring = allring[lct];
	myfit.cscid = allcscid[lct];
        myfit.rawZ = allzG[lct];
        myfit.calphi = DeltaPhiEta.calphi;
        myfit.caleta = DeltaPhiEta.caleta;
        myfit.p0 = parFit[0];
        myfit.p1 = parFit[1];
        myfit.p2 = parFit[2];
        myfit.p3 = parFit[3];
        myfit.distance = distance;
	myfit.dr0 = dr0;
        myfit.chi2 = amin;

        LUTrels.push_back(myfit);
    }

    return LUTrels;
}

void readLCTAnalyzer_OneHitperStation(TString _input_="./csctf_Run246926.root", TString _output_="output_test.root", int entry0=0)
{

    TFile *InFile = TFile::Open(_input_, "READ");
    TTree *t_ = (TTree *) InFile->Get("OfflineDQMCSCTF/data");

    TString outputLOG = _output_;
    outputLOG.ReplaceAll(".root",".log");
    ofstream outfile;
    outfile.open ("/tmp/rewang/"+outputLOG);



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



    t_->SetBranchAddress("nlcts_m", &nlcts_m);
    t_->SetBranchAddress("lct_m_gblphi", lct_m_gblphi);
    t_->SetBranchAddress("lct_m_lclphi", lct_m_lclphi);
    t_->SetBranchAddress("lct_m_pkdphi", lct_m_pkdphi);
    t_->SetBranchAddress("lct_m_gbleta", lct_m_gbleta);
    t_->SetBranchAddress("lct_m_pkdeta", lct_m_pkdeta);
    t_->SetBranchAddress("lct_m_gblZ", lct_m_gblZ);
    t_->SetBranchAddress("lct_m_station", lct_m_station);
    t_->SetBranchAddress("lct_m_ring", lct_m_ring);
    t_->SetBranchAddress("lct_m_cscid", lct_m_cscid);
    t_->SetBranchAddress("lct_m_endcap", lct_m_endcap);
    t_->SetBranchAddress("lct_m_sector", lct_m_sector);
    t_->SetBranchAddress("lct_m_bptx", lct_m_bptx);
    t_->SetBranchAddress("lct_m_strip", lct_m_strip);
    t_->SetBranchAddress("lct_m_keywire", lct_m_keywire);

    t_->SetBranchAddress("nlcts_p", &nlcts_p);
    t_->SetBranchAddress("lct_p_gblphi", lct_p_gblphi);
    t_->SetBranchAddress("lct_p_lclphi", lct_p_lclphi);
    t_->SetBranchAddress("lct_p_pkdphi", lct_p_pkdphi);
    t_->SetBranchAddress("lct_p_gbleta", lct_p_gbleta);
    t_->SetBranchAddress("lct_p_pkdeta", lct_p_pkdeta);
    t_->SetBranchAddress("lct_p_gblZ", lct_p_gblZ);
    t_->SetBranchAddress("lct_p_station", lct_p_station);
    t_->SetBranchAddress("lct_p_ring", lct_p_ring);
    t_->SetBranchAddress("lct_p_cscid", lct_p_cscid);
    t_->SetBranchAddress("lct_p_endcap", lct_p_endcap);
    t_->SetBranchAddress("lct_p_sector", lct_p_sector);
    t_->SetBranchAddress("lct_p_bptx", lct_p_bptx);
    t_->SetBranchAddress("lct_p_strip", lct_p_strip);
    t_->SetBranchAddress("lct_p_keywire", lct_p_keywire);



    //hists
    TH1F *h_deltaPhi = new TH1F("h_deltaPhi",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",500,-.8,.8);
    TH1F *h_deltaEta = new TH1F("h_deltaEta",";#it{#eta}(Fit) - #it{#eta}(LUT);Events",500,-0.8,0.8);

    TH2F *h_deltaEtavsPhi_MEp = new TH2F("h_deltaEtavsPhi_MEp","; #it{#phi}(LUT) [rad], all ME+; #it{#eta}(Fit) - #it{#eta}(LUT);Events",500,0,2.*M_PI,500,-1.,1.);
    TH2F *h_deltaEtavsPhi_MEm = new TH2F("h_deltaEtavsPhi_MEm","; #it{#phi}(LUT) [rad], all ME-; #it{#eta}(Fit) - #it{#eta}(LUT);Events",500,0,2.*M_PI,500,-1.,1.);

    TH2F *h_deltaPhivsEta_MEp = new TH2F("h_deltaPhivsEta_MEp","; #it{#eta}(LUT), all ME+; #it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",500,0.8,2.5,500,-1.,1.);
    TH2F *h_deltaPhivsEta_MEm = new TH2F("h_deltaPhivsEta_MEm","; #it{#eta}(LUT), all ME-; #it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",500,-2.5,-0.8,500,-1.,1.);

    /////////////////////

    TH2F *h_dEtavs10degPhi_MEp = new TH2F("h_dEtavs10degPhi_MEp","; #it{#phi}(LUT) [degree], all ME+; #it{#eta}(Fit) - #it{#eta}(LUT);Events",288,0,360,500,-0.2,0.2);
    TH2F *h_dPhivs10degPhi_MEp = new TH2F("h_dPhivs10degPhi_MEp","; #it{#phi}(LUT) [degree], all ME+; #it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",288,0,360,500,-0.2,0.2);

    TH2F *h_dEtavs10degPhi_MEm = new TH2F("h_dEtavs10degPhi_MEm","; #it{#phi}(LUT) [degree], all ME-; #it{#eta}(Fit) - #it{#eta}(LUT);Events",288,0,360,500,-0.2,0.2);
    TH2F *h_dPhivs10degPhi_MEm = new TH2F("h_dPhivs10degPhi_MEm","; #it{#phi}(LUT) [degree], all ME-; #it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",288,0,360,500,-0.2,0.2);

    /////////////////////

    TH2F *h_dEtavs10degPhi_MEp11 = new TH2F("h_dEtavs10degPhi_MEp11","; #it{#phi}(Fit) [degree], ME+1/1; #it{#eta}(Fit) - #it{#eta}(LUT);Events",288,0,360,100,-0.2,0.2);
    TH2F *h_dEtavs10degPhi_MEp11a = new TH2F("h_dEtavs10degPhi_MEp11a","; #it{#phi}(Fit) [degree], ME+1/1a; #it{#eta}(Fit) - #it{#eta}(LUT);Events",288,0,360,100,-0.2,0.2);
    TH2F *h_dEtavs10degPhi_MEp11b = new TH2F("h_dEtavs10degPhi_MEp11b","; #it{#phi}(Fit) [degree], ME+1/1b; #it{#eta}(Fit) - #it{#eta}(LUT);Events",288,0,360,100,-0.2,0.2);
    TH2F *h_dEtavs10degPhi_MEp12 = new TH2F("h_dEtavs10degPhi_MEp12","; #it{#phi}(Fit) [degree], ME+1/2; #it{#eta}(Fit) - #it{#eta}(LUT);Events",288,0,360,100,-0.2,0.2);
    TH2F *h_dEtavs10degPhi_MEp13 = new TH2F("h_dEtavs10degPhi_MEp13","; #it{#phi}(Fit) [degree], ME+1/3; #it{#eta}(Fit) - #it{#eta}(LUT);Events",288,0,360,100,-0.2,0.2);
    TH2F *h_dEtavs10degPhi_MEp21 = new TH2F("h_dEtavs10degPhi_MEp21","; #it{#phi}(Fit) [degree], ME+2/1; #it{#eta}(Fit) - #it{#eta}(LUT);Events",288,0,360,100,-0.2,0.2);
    TH2F *h_dEtavs10degPhi_MEp22 = new TH2F("h_dEtavs10degPhi_MEp22","; #it{#phi}(Fit) [degree], ME+2/2; #it{#eta}(Fit) - #it{#eta}(LUT);Events",288,0,360,100,-0.2,0.2);
    TH2F *h_dEtavs10degPhi_MEp31 = new TH2F("h_dEtavs10degPhi_MEp31","; #it{#phi}(Fit) [degree], ME+3/1; #it{#eta}(Fit) - #it{#eta}(LUT);Events",288,0,360,100,-0.2,0.2);
    TH2F *h_dEtavs10degPhi_MEp32 = new TH2F("h_dEtavs10degPhi_MEp32","; #it{#phi}(Fit) [degree], ME+3/2; #it{#eta}(Fit) - #it{#eta}(LUT);Events",288,0,360,100,-0.2,0.2);
    TH2F *h_dEtavs10degPhi_MEp41 = new TH2F("h_dEtavs10degPhi_MEp41","; #it{#phi}(Fit) [degree], ME+4/1; #it{#eta}(Fit) - #it{#eta}(LUT);Events",288,0,360,100,-0.2,0.2);
    TH2F *h_dEtavs10degPhi_MEp42 = new TH2F("h_dEtavs10degPhi_MEp42","; #it{#phi}(Fit) [degree], ME+4/2; #it{#eta}(Fit) - #it{#eta}(LUT);Events",288,0,360,100,-0.2,0.2);

    TH2F *h_dEtavs10degPhi_MEm11 = new TH2F("h_dEtavs10degPhi_MEm11","; #it{#phi}(Fit) [degree], ME-1/1; #it{#eta}(Fit) - #it{#eta}(LUT);Events",288,0,360,100,-0.2,0.2);
    TH2F *h_dEtavs10degPhi_MEm11a = new TH2F("h_dEtavs10degPhi_MEm11a","; #it{#phi}(Fit) [degree], ME-1/1a; #it{#eta}(Fit) - #it{#eta}(LUT);Events",288,0,360,100,-0.2,0.2);
    TH2F *h_dEtavs10degPhi_MEm11b = new TH2F("h_dEtavs10degPhi_MEm11b","; #it{#phi}(Fit) [degree], ME-1/1b; #it{#eta}(Fit) - #it{#eta}(LUT);Events",288,0,360,100,-0.2,0.2);
    TH2F *h_dEtavs10degPhi_MEm12 = new TH2F("h_dEtavs10degPhi_MEm12","; #it{#phi}(Fit) [degree], ME-1/2; #it{#eta}(Fit) - #it{#eta}(LUT);Events",288,0,360,100,-0.2,0.2);
    TH2F *h_dEtavs10degPhi_MEm13 = new TH2F("h_dEtavs10degPhi_MEm13","; #it{#phi}(Fit) [degree], ME-1/3; #it{#eta}(Fit) - #it{#eta}(LUT);Events",288,0,360,100,-0.2,0.2);
    TH2F *h_dEtavs10degPhi_MEm21 = new TH2F("h_dEtavs10degPhi_MEm21","; #it{#phi}(Fit) [degree], ME-2/1; #it{#eta}(Fit) - #it{#eta}(LUT);Events",288,0,360,100,-0.2,0.2);
    TH2F *h_dEtavs10degPhi_MEm22 = new TH2F("h_dEtavs10degPhi_MEm22","; #it{#phi}(Fit) [degree], ME-2/2; #it{#eta}(Fit) - #it{#eta}(LUT);Events",288,0,360,100,-0.2,0.2);
    TH2F *h_dEtavs10degPhi_MEm31 = new TH2F("h_dEtavs10degPhi_MEm31","; #it{#phi}(Fit) [degree], ME-3/1; #it{#eta}(Fit) - #it{#eta}(LUT);Events",288,0,360,100,-0.2,0.2);
    TH2F *h_dEtavs10degPhi_MEm32 = new TH2F("h_dEtavs10degPhi_MEm32","; #it{#phi}(Fit) [degree], ME-3/2; #it{#eta}(Fit) - #it{#eta}(LUT);Events",288,0,360,100,-0.2,0.2);
    TH2F *h_dEtavs10degPhi_MEm41 = new TH2F("h_dEtavs10degPhi_MEm41","; #it{#phi}(Fit) [degree], ME-4/1; #it{#eta}(Fit) - #it{#eta}(LUT);Events",288,0,360,100,-0.2,0.2);
    TH2F *h_dEtavs10degPhi_MEm42 = new TH2F("h_dEtavs10degPhi_MEm42","; #it{#phi}(Fit) [degree], ME-4/2; #it{#eta}(Fit) - #it{#eta}(LUT);Events",288,0,360,100,-0.2,0.2);

    /////////////////////

    TH2F *h_dEtavsFitEta_MEp11 = new TH2F("h_dEtavsFitEta_MEp11","; #it{#eta}(Fit), ME+1/1; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,0.9,2.5,100,-0.2,0.2);
    TH2F *h_dEtavsFitEta_MEp11a = new TH2F("h_dEtavsFitEta_MEp11a","; #it{#eta}(Fit), ME+1/1a; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,0.9,2.5,100,-0.2,0.2);
    TH2F *h_dEtavsFitEta_MEp11b = new TH2F("h_dEtavsFitEta_MEp11b","; #it{#eta}(Fit), ME+1/1b; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,0.9,2.5,100,-0.2,0.2);
    TH2F *h_dEtavsFitEta_MEp12 = new TH2F("h_dEtavsFitEta_MEp12","; #it{#eta}(Fit), ME+1/2; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,0.9,2.5,100,-0.2,0.2);
    TH2F *h_dEtavsFitEta_MEp13 = new TH2F("h_dEtavsFitEta_MEp13","; #it{#eta}(Fit), ME+1/3; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,0.9,2.5,100,-0.2,0.2);
    TH2F *h_dEtavsFitEta_MEp21 = new TH2F("h_dEtavsFitEta_MEp21","; #it{#eta}(Fit), ME+2/1; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,0.9,2.5,100,-0.2,0.2);
    TH2F *h_dEtavsFitEta_MEp22 = new TH2F("h_dEtavsFitEta_MEp22","; #it{#eta}(Fit), ME+2/2; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,0.9,2.5,100,-0.2,0.2);
    TH2F *h_dEtavsFitEta_MEp31 = new TH2F("h_dEtavsFitEta_MEp31","; #it{#eta}(Fit), ME+3/1; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,0.9,2.5,100,-0.2,0.2);
    TH2F *h_dEtavsFitEta_MEp32 = new TH2F("h_dEtavsFitEta_MEp32","; #it{#eta}(Fit), ME+3/2; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,0.9,2.5,100,-0.2,0.2);
    TH2F *h_dEtavsFitEta_MEp41 = new TH2F("h_dEtavsFitEta_MEp41","; #it{#eta}(Fit), ME+4/1; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,0.9,2.5,100,-0.2,0.2);
    TH2F *h_dEtavsFitEta_MEp42 = new TH2F("h_dEtavsFitEta_MEp42","; #it{#eta}(Fit), ME+4/2; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,0.9,2.5,100,-0.2,0.2);

    TH2F *h_dEtavsFitEta_MEm11 = new TH2F("h_dEtavsFitEta_MEm11","; #it{#eta}(Fit), ME-1/1; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,-2.5,-0.9,100,-0.2,0.2);
    TH2F *h_dEtavsFitEta_MEm11a = new TH2F("h_dEtavsFitEta_MEm11a","; #it{#eta}(Fit), ME-1/1a; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,-2.5,-0.9,100,-0.2,0.2);
    TH2F *h_dEtavsFitEta_MEm11b = new TH2F("h_dEtavsFitEta_MEm11b","; #it{#eta}(Fit), ME-1/1b; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,-2.5,-0.9,100,-0.2,0.2);
    TH2F *h_dEtavsFitEta_MEm12 = new TH2F("h_dEtavsFitEta_MEm12","; #it{#eta}(Fit), ME-1/2; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,-2.5,-0.9,100,-0.2,0.2);
    TH2F *h_dEtavsFitEta_MEm13 = new TH2F("h_dEtavsFitEta_MEm13","; #it{#eta}(Fit), ME-1/3; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,-2.5,-0.9,100,-0.2,0.2);
    TH2F *h_dEtavsFitEta_MEm21 = new TH2F("h_dEtavsFitEta_MEm21","; #it{#eta}(Fit), ME-2/1; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,-2.5,-0.9,100,-0.2,0.2);
    TH2F *h_dEtavsFitEta_MEm22 = new TH2F("h_dEtavsFitEta_MEm22","; #it{#eta}(Fit), ME-2/2; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,-2.5,-0.9,100,-0.2,0.2);
    TH2F *h_dEtavsFitEta_MEm31 = new TH2F("h_dEtavsFitEta_MEm31","; #it{#eta}(Fit), ME-3/1; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,-2.5,-0.9,100,-0.2,0.2);
    TH2F *h_dEtavsFitEta_MEm32 = new TH2F("h_dEtavsFitEta_MEm32","; #it{#eta}(Fit), ME-3/2; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,-2.5,-0.9,100,-0.2,0.2);
    TH2F *h_dEtavsFitEta_MEm41 = new TH2F("h_dEtavsFitEta_MEm41","; #it{#eta}(Fit), ME-4/1; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,-2.5,-0.9,100,-0.2,0.2);
    TH2F *h_dEtavsFitEta_MEm42 = new TH2F("h_dEtavsFitEta_MEm42","; #it{#eta}(Fit), ME-4/2; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,-2.5,-0.9,100,-0.2,0.2);


    TH2F *h_dEtavsLUTEta_MEp11 = new TH2F("h_dEtavsLUTEta_MEp11","; #it{#eta}(LUT), ME+1/1; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,0.9,2.5,100,-0.2,0.2);
    TH2F *h_dEtavsLUTEta_MEp11a = new TH2F("h_dEtavsLUTEta_MEp11a","; #it{#eta}(LUT), ME+1/1a; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,0.9,2.5,100,-0.2,0.2);
    TH2F *h_dEtavsLUTEta_MEp11b = new TH2F("h_dEtavsLUTEta_MEp11b","; #it{#eta}(LUT), ME+1/1b; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,0.9,2.5,100,-0.2,0.2);
    TH2F *h_dEtavsLUTEta_MEp12 = new TH2F("h_dEtavsLUTEta_MEp12","; #it{#eta}(LUT), ME+1/2; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,0.9,2.5,100,-0.2,0.2);
    TH2F *h_dEtavsLUTEta_MEp13 = new TH2F("h_dEtavsLUTEta_MEp13","; #it{#eta}(LUT), ME+1/3; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,0.9,2.5,100,-0.2,0.2);
    TH2F *h_dEtavsLUTEta_MEp21 = new TH2F("h_dEtavsLUTEta_MEp21","; #it{#eta}(LUT), ME+2/1; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,0.9,2.5,100,-0.2,0.2);
    TH2F *h_dEtavsLUTEta_MEp22 = new TH2F("h_dEtavsLUTEta_MEp22","; #it{#eta}(LUT), ME+2/2; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,0.9,2.5,100,-0.2,0.2);
    TH2F *h_dEtavsLUTEta_MEp31 = new TH2F("h_dEtavsLUTEta_MEp31","; #it{#eta}(LUT), ME+3/1; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,0.9,2.5,100,-0.2,0.2);
    TH2F *h_dEtavsLUTEta_MEp32 = new TH2F("h_dEtavsLUTEta_MEp32","; #it{#eta}(LUT), ME+3/2; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,0.9,2.5,100,-0.2,0.2);
    TH2F *h_dEtavsLUTEta_MEp41 = new TH2F("h_dEtavsLUTEta_MEp41","; #it{#eta}(LUT), ME+4/1; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,0.9,2.5,100,-0.2,0.2);
    TH2F *h_dEtavsLUTEta_MEp42 = new TH2F("h_dEtavsLUTEta_MEp42","; #it{#eta}(LUT), ME+4/2; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,0.9,2.5,100,-0.2,0.2);

    TH2F *h_dEtavsLUTEta_MEm11 = new TH2F("h_dEtavsLUTEta_MEm11","; #it{#eta}(LUT), ME-1/1; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,-2.5,-0.9,100,-0.2,0.2);
    TH2F *h_dEtavsLUTEta_MEm11a = new TH2F("h_dEtavsLUTEta_MEm11a","; #it{#eta}(LUT), ME-1/1a; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,-2.5,-0.9,100,-0.2,0.2);
    TH2F *h_dEtavsLUTEta_MEm11b = new TH2F("h_dEtavsLUTEta_MEm11b","; #it{#eta}(LUT), ME-1/1b; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,-2.5,-0.9,100,-0.2,0.2);
    TH2F *h_dEtavsLUTEta_MEm12 = new TH2F("h_dEtavsLUTEta_MEm12","; #it{#eta}(LUT), ME-1/2; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,-2.5,-0.9,100,-0.2,0.2);
    TH2F *h_dEtavsLUTEta_MEm13 = new TH2F("h_dEtavsLUTEta_MEm13","; #it{#eta}(LUT), ME-1/3; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,-2.5,-0.9,100,-0.2,0.2);
    TH2F *h_dEtavsLUTEta_MEm21 = new TH2F("h_dEtavsLUTEta_MEm21","; #it{#eta}(LUT), ME-2/1; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,-2.5,-0.9,100,-0.2,0.2);
    TH2F *h_dEtavsLUTEta_MEm22 = new TH2F("h_dEtavsLUTEta_MEm22","; #it{#eta}(LUT), ME-2/2; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,-2.5,-0.9,100,-0.2,0.2);
    TH2F *h_dEtavsLUTEta_MEm31 = new TH2F("h_dEtavsLUTEta_MEm31","; #it{#eta}(LUT), ME-3/1; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,-2.5,-0.9,100,-0.2,0.2);
    TH2F *h_dEtavsLUTEta_MEm32 = new TH2F("h_dEtavsLUTEta_MEm32","; #it{#eta}(LUT), ME-3/2; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,-2.5,-0.9,100,-0.2,0.2);
    TH2F *h_dEtavsLUTEta_MEm41 = new TH2F("h_dEtavsLUTEta_MEm41","; #it{#eta}(LUT), ME-4/1; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,-2.5,-0.9,100,-0.2,0.2);
    TH2F *h_dEtavsLUTEta_MEm42 = new TH2F("h_dEtavsLUTEta_MEm42","; #it{#eta}(LUT), ME-4/2; #it{#eta}(Fit) - #it{#eta}(LUT);Events",100,-2.5,-0.9,100,-0.2,0.2);


    ////////////////////

    TH2F *h_dPhivs10degPhi_MEp11 = new TH2F("h_dPhivs10degPhi_MEp11","; #it{#phi}(Fit) [degree], ME+1/1; #it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",288,0,360,500,-0.2,0.2);
    TH2F *h_dPhivs10degPhi_MEp11a = new TH2F("h_dPhivs10degPhi_MEp11a","; #it{#phi}(Fit) [degree], ME+1/1a; #it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",288,0,360,500,-0.2,0.2);
    TH2F *h_dPhivs10degPhi_MEp11b = new TH2F("h_dPhivs10degPhi_MEp11b","; #it{#phi}(Fit) [degree], ME+1/1b; #it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",288,0,360,500,-0.2,0.2);
    TH2F *h_dPhivs10degPhi_MEp12 = new TH2F("h_dPhivs10degPhi_MEp12","; #it{#phi}(Fit) [degree], ME+1/2; #it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",288,0,360,500,-0.2,0.2);
    TH2F *h_dPhivs10degPhi_MEp13 = new TH2F("h_dPhivs10degPhi_MEp13","; #it{#phi}(Fit) [degree], ME+1/3; #it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",288,0,360,500,-0.2,0.2);
    TH2F *h_dPhivs10degPhi_MEp21 = new TH2F("h_dPhivs10degPhi_MEp21","; #it{#phi}(Fit) [degree], ME+2/1; #it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",288,0,360,500,-0.2,0.2);
    TH2F *h_dPhivs10degPhi_MEp22 = new TH2F("h_dPhivs10degPhi_MEp22","; #it{#phi}(Fit) [degree], ME+2/2; #it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",288,0,360,500,-0.2,0.2);
    TH2F *h_dPhivs10degPhi_MEp31 = new TH2F("h_dPhivs10degPhi_MEp31","; #it{#phi}(Fit) [degree], ME+3/1; #it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",288,0,360,500,-0.2,0.2);
    TH2F *h_dPhivs10degPhi_MEp32 = new TH2F("h_dPhivs10degPhi_MEp32","; #it{#phi}(Fit) [degree], ME+3/2; #it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",288,0,360,500,-0.2,0.2);
    TH2F *h_dPhivs10degPhi_MEp41 = new TH2F("h_dPhivs10degPhi_MEp41","; #it{#phi}(Fit) [degree], ME+4/1; #it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",288,0,360,500,-0.2,0.2);
    TH2F *h_dPhivs10degPhi_MEp42 = new TH2F("h_dPhivs10degPhi_MEp42","; #it{#phi}(Fit) [degree], ME+4/2; #it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",288,0,360,500,-0.2,0.2);

    TH2F *h_dPhivs10degPhi_MEm11 = new TH2F("h_dPhivs10degPhi_MEm11","; #it{#phi}(Fit) [degree], ME-1/1; #it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",288,0,360,500,-0.2,0.2);
    TH2F *h_dPhivs10degPhi_MEm11a = new TH2F("h_dPhivs10degPhi_MEm11a","; #it{#phi}(Fit) [degree], ME-1/1a; #it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",288,0,360,500,-0.2,0.2);
    TH2F *h_dPhivs10degPhi_MEm11b = new TH2F("h_dPhivs10degPhi_MEm11b","; #it{#phi}(Fit) [degree], ME-1/1b; #it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",288,0,360,500,-0.2,0.2);
    TH2F *h_dPhivs10degPhi_MEm12 = new TH2F("h_dPhivs10degPhi_MEm12","; #it{#phi}(Fit) [degree], ME-1/2; #it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",288,0,360,500,-0.2,0.2);
    TH2F *h_dPhivs10degPhi_MEm13 = new TH2F("h_dPhivs10degPhi_MEm13","; #it{#phi}(Fit) [degree], ME-1/3; #it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",288,0,360,500,-0.2,0.2);
    TH2F *h_dPhivs10degPhi_MEm21 = new TH2F("h_dPhivs10degPhi_MEm21","; #it{#phi}(Fit) [degree], ME-2/1; #it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",288,0,360,500,-0.2,0.2);
    TH2F *h_dPhivs10degPhi_MEm22 = new TH2F("h_dPhivs10degPhi_MEm22","; #it{#phi}(Fit) [degree], ME-2/2; #it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",288,0,360,500,-0.2,0.2);
    TH2F *h_dPhivs10degPhi_MEm31 = new TH2F("h_dPhivs10degPhi_MEm31","; #it{#phi}(Fit) [degree], ME-3/1; #it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",288,0,360,500,-0.2,0.2);
    TH2F *h_dPhivs10degPhi_MEm32 = new TH2F("h_dPhivs10degPhi_MEm32","; #it{#phi}(Fit) [degree], ME-3/2; #it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",288,0,360,500,-0.2,0.2);
    TH2F *h_dPhivs10degPhi_MEm41 = new TH2F("h_dPhivs10degPhi_MEm41","; #it{#phi}(Fit) [degree], ME-4/1; #it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",288,0,360,500,-0.2,0.2);
    TH2F *h_dPhivs10degPhi_MEm42 = new TH2F("h_dPhivs10degPhi_MEm42","; #it{#phi}(Fit) [degree], ME-4/2; #it{#phi}(Fit) - #it{#phi}(LUT) [rad];Events",288,0,360,500,-0.2,0.2);

    ///////////////









    TH1F* h_deltaPhiMEp11 = new TH1F("h_deltaPhiMEp11",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME+1/1;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEp11a = new TH1F("h_deltaPhiMEp11a",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME+1/1a;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEp11b = new TH1F("h_deltaPhiMEp11b",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME+1/1b;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEp12 = new TH1F("h_deltaPhiMEp12",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME+1/2;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEp13 = new TH1F("h_deltaPhiMEp13",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME+1/3;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEp21 = new TH1F("h_deltaPhiMEp21",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME+2/1;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEp22 = new TH1F("h_deltaPhiMEp22",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME+2/2;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEp31 = new TH1F("h_deltaPhiMEp31",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME+3/1;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEp32 = new TH1F("h_deltaPhiMEp32",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME+3/2;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEp41 = new TH1F("h_deltaPhiMEp41",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME+4/1;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEp42 = new TH1F("h_deltaPhiMEp42",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME+4/2;Events",500,-.8,.8);

    TH1F* h_deltaPhiMEm11 = new TH1F("h_deltaPhiMEm11",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME-1/1;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEm11a = new TH1F("h_deltaPhiMEm11a",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME-1/1a;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEm11b = new TH1F("h_deltaPhiMEm11b",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME-1/1b;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEm12 = new TH1F("h_deltaPhiMEm12",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME-1/2;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEm13 = new TH1F("h_deltaPhiMEm13",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME-1/3;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEm21 = new TH1F("h_deltaPhiMEm21",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME-2/1;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEm22 = new TH1F("h_deltaPhiMEm22",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME-2/2;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEm31 = new TH1F("h_deltaPhiMEm31",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME-3/1;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEm32 = new TH1F("h_deltaPhiMEm32",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME-3/2;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEm41 = new TH1F("h_deltaPhiMEm41",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME-4/1;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEm42 = new TH1F("h_deltaPhiMEm42",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME-4/2;Events",500,-.8,.8);

    TH1F* h_deltaEtaMEp11 = new TH1F("h_deltaEtaMEp11",";#it{#eta}(Fit) - #it{#eta}(LUT), ME+1/1;Events",500,-0.8,0.8);
    TH1F* h_deltaEtaMEp11a = new TH1F("h_deltaEtaMEp11a",";#it{#eta}(Fit) - #it{#eta}(LUT), ME+1/1a;Events",500,-0.8,0.8);
    TH1F* h_deltaEtaMEp11b = new TH1F("h_deltaEtaMEp11b",";#it{#eta}(Fit) - #it{#eta}(LUT), ME+1/1b;Events",500,-0.8,0.8);
    TH1F* h_deltaEtaMEp12 = new TH1F("h_deltaEtaMEp12",";#it{#eta}(Fit) - #it{#eta}(LUT), ME+1/2;Events",500,-0.8,0.8);
    TH1F* h_deltaEtaMEp13 = new TH1F("h_deltaEtaMEp13",";#it{#eta}(Fit) - #it{#eta}(LUT), ME+1/3;Events",500,-0.8,0.8);
    TH1F* h_deltaEtaMEp21 = new TH1F("h_deltaEtaMEp21",";#it{#eta}(Fit) - #it{#eta}(LUT), ME+2/1;Events",500,-0.8,0.8);
    TH1F* h_deltaEtaMEp22 = new TH1F("h_deltaEtaMEp22",";#it{#eta}(Fit) - #it{#eta}(LUT), ME+2/2;Events",500,-0.8,0.8);
    TH1F* h_deltaEtaMEp31 = new TH1F("h_deltaEtaMEp31",";#it{#eta}(Fit) - #it{#eta}(LUT), ME+3/1;Events",500,-0.8,0.8);
    TH1F* h_deltaEtaMEp32 = new TH1F("h_deltaEtaMEp32",";#it{#eta}(Fit) - #it{#eta}(LUT), ME+3/2;Events",500,-0.8,0.8);
    TH1F* h_deltaEtaMEp41 = new TH1F("h_deltaEtaMEp41",";#it{#eta}(Fit) - #it{#eta}(LUT), ME+4/1;Events",500,-0.8,0.8);
    TH1F* h_deltaEtaMEp42 = new TH1F("h_deltaEtaMEp42",";#it{#eta}(Fit) - #it{#eta}(LUT), ME+4/2;Events",500,-0.8,0.8);

    TH1F* h_deltaEtaMEm11 = new TH1F("h_deltaEtaMEm11",";#it{#eta}(Fit) - #it{#eta}(LUT), ME-1/1;Events",500,-0.8,0.8);
    TH1F* h_deltaEtaMEm11a = new TH1F("h_deltaEtaMEm11a",";#it{#eta}(Fit) - #it{#eta}(LUT), ME-1/1a;Events",500,-0.8,0.8);
    TH1F* h_deltaEtaMEm11b = new TH1F("h_deltaEtaMEm11b",";#it{#eta}(Fit) - #it{#eta}(LUT), ME-1/1b;Events",500,-0.8,0.8);
    TH1F* h_deltaEtaMEm12 = new TH1F("h_deltaEtaMEm12",";#it{#eta}(Fit) - #it{#eta}(LUT), ME-1/2;Events",500,-0.8,0.8);
    TH1F* h_deltaEtaMEm13 = new TH1F("h_deltaEtaMEm13",";#it{#eta}(Fit) - #it{#eta}(LUT), ME-1/3;Events",500,-0.8,0.8);
    TH1F* h_deltaEtaMEm21 = new TH1F("h_deltaEtaMEm21",";#it{#eta}(Fit) - #it{#eta}(LUT), ME-2/1;Events",500,-0.8,0.8);
    TH1F* h_deltaEtaMEm22 = new TH1F("h_deltaEtaMEm22",";#it{#eta}(Fit) - #it{#eta}(LUT), ME-2/2;Events",500,-0.8,0.8);
    TH1F* h_deltaEtaMEm31 = new TH1F("h_deltaEtaMEm31",";#it{#eta}(Fit) - #it{#eta}(LUT), ME-3/1;Events",500,-0.8,0.8);
    TH1F* h_deltaEtaMEm32 = new TH1F("h_deltaEtaMEm32",";#it{#eta}(Fit) - #it{#eta}(LUT), ME-3/2;Events",500,-0.8,0.8);
    TH1F* h_deltaEtaMEm41 = new TH1F("h_deltaEtaMEm41",";#it{#eta}(Fit) - #it{#eta}(LUT), ME-4/1;Events",500,-0.8,0.8);
    TH1F* h_deltaEtaMEm42 = new TH1F("h_deltaEtaMEm42",";#it{#eta}(Fit) - #it{#eta}(LUT), ME-4/2;Events",500,-0.8,0.8);



    TH1F* h_deltaPhiMEpSec0 = new TH1F("h_deltaPhiMEpSec0",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME+, Sector 0;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEpSec1 = new TH1F("h_deltaPhiMEpSec1",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME+, Sector 1;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEpSec2 = new TH1F("h_deltaPhiMEpSec2",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME+, Sector 2;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEpSec3 = new TH1F("h_deltaPhiMEpSec3",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME+, Sector 3;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEpSec4 = new TH1F("h_deltaPhiMEpSec4",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME+, Sector 4;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEpSec5 = new TH1F("h_deltaPhiMEpSec5",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME+, Sector 5;Events",500,-.8,.8);

    TH1F* h_deltaPhiMEmSec0 = new TH1F("h_deltaPhiMEmSec0",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME-, Sector 0;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEmSec1 = new TH1F("h_deltaPhiMEmSec1",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME-, Sector 1;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEmSec2 = new TH1F("h_deltaPhiMEmSec2",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME-, Sector 2;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEmSec3 = new TH1F("h_deltaPhiMEmSec3",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME-, Sector 3;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEmSec4 = new TH1F("h_deltaPhiMEmSec4",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME-, Sector 4;Events",500,-.8,.8);
    TH1F* h_deltaPhiMEmSec5 = new TH1F("h_deltaPhiMEmSec5",";#it{#phi}(Fit) - #it{#phi}(LUT) [rad], ME-, Sector 5;Events",500,-.8,.8);

    TH1F* h_deltaEtaMEpSec0 = new TH1F("h_deltaEtaMEpSec0",";#it{#eta}(Fit) - #it{#eta}(LUT), ME+, Sector 0;Events",500,-.8,.8);
    TH1F* h_deltaEtaMEpSec1 = new TH1F("h_deltaEtaMEpSec1",";#it{#eta}(Fit) - #it{#eta}(LUT), ME+, Sector 1;Events",500,-.8,.8);
    TH1F* h_deltaEtaMEpSec2 = new TH1F("h_deltaEtaMEpSec2",";#it{#eta}(Fit) - #it{#eta}(LUT), ME+, Sector 2;Events",500,-.8,.8);
    TH1F* h_deltaEtaMEpSec3 = new TH1F("h_deltaEtaMEpSec3",";#it{#eta}(Fit) - #it{#eta}(LUT), ME+, Sector 3;Events",500,-.8,.8);
    TH1F* h_deltaEtaMEpSec4 = new TH1F("h_deltaEtaMEpSec4",";#it{#eta}(Fit) - #it{#eta}(LUT), ME+, Sector 4;Events",500,-.8,.8);
    TH1F* h_deltaEtaMEpSec5 = new TH1F("h_deltaEtaMEpSec5",";#it{#eta}(Fit) - #it{#eta}(LUT), ME+, Sector 5;Events",500,-.8,.8);

    TH1F* h_deltaEtaMEmSec0 = new TH1F("h_deltaEtaMEmSec0",";#it{#eta}(Fit) - #it{#eta}(LUT), ME-, Sector 0;Events",500,-.8,.8);
    TH1F* h_deltaEtaMEmSec1 = new TH1F("h_deltaEtaMEmSec1",";#it{#eta}(Fit) - #it{#eta}(LUT), ME-, Sector 1;Events",500,-.8,.8);
    TH1F* h_deltaEtaMEmSec2 = new TH1F("h_deltaEtaMEmSec2",";#it{#eta}(Fit) - #it{#eta}(LUT), ME-, Sector 2;Events",500,-.8,.8);
    TH1F* h_deltaEtaMEmSec3 = new TH1F("h_deltaEtaMEmSec3",";#it{#eta}(Fit) - #it{#eta}(LUT), ME-, Sector 3;Events",500,-.8,.8);
    TH1F* h_deltaEtaMEmSec4 = new TH1F("h_deltaEtaMEmSec4",";#it{#eta}(Fit) - #it{#eta}(LUT), ME-, Sector 4;Events",500,-.8,.8);
    TH1F* h_deltaEtaMEmSec5 = new TH1F("h_deltaEtaMEmSec5",";#it{#eta}(Fit) - #it{#eta}(LUT), ME-, Sector 5;Events",500,-.8,.8);






    TH1F* h_m_dr0 = new TH1F("h_m_dr0",";Distance to IP [mm], all ME-;Events",500,0,2000);
    TH1F* h_p_dr0 = new TH1F("h_p_dr0",";Distance to IP [mm], all ME+;Events",500,0,2000);

    TH1F* h_m_chi2 = new TH1F("h_m_chi2",";Chi2 [mm], all ME-;Events",500,0,1);
    TH1F* h_p_chi2 = new TH1F("h_p_chi2",";Chi2 [mm], all ME+;Events",500,0,1);

    TH1F* h_m_p0 = new TH1F("h_m_p0",";p_{0}, Minus EndCap;Events",500,-2000,2000);
    TH1F* h_m_p1 = new TH1F("h_m_p1",";p_{1}, Minus EndCap;Events",500,-2000,2000);
    TH1F* h_m_p2 = new TH1F("h_m_p2",";p_{2}, Minus EndCap;Events",500,-2000,2000);
    TH1F* h_m_p3 = new TH1F("h_m_p3",";p_{3}, Minus EndCap;Events",500,-2000,2000);
    TH1F *h_m_distance = new TH1F("h_m_distance",";Distance, Minus EndCap;Events",500,0,2000.);

    TH1F* h_p_p0 = new TH1F("h_p_p0",";p_{0}, Plus EndCap;Events",500,-2000,2000);
    TH1F* h_p_p1 = new TH1F("h_p_p1",";p_{1}, Plus EndCap;Events",500,-2000,2000);
    TH1F* h_p_p2 = new TH1F("h_p_p2",";p_{2}, Plus EndCap;Events",500,-2000,2000);
    TH1F* h_p_p3 = new TH1F("h_p_p3",";p_{3}, Plus EndCap;Events",500,-2000,2000);
    TH1F *h_p_distance = new TH1F("h_p_distance",";Distance, Plus EndCap;Events",500,0,2000.);





    //int evStart = 0;
    //int evEnd = -1;
    //printf("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");
    //printf("Scanning the ntuple :");
    // all entries and fill the histograms
    Int_t totnentries = (Int_t)t_->GetEntries();
    int nentries=50000;
    int entry1 = entry0+nentries;
    if(entry0>totnentries) return;
    if(entry1>totnentries) entry1=totnentries;
    //int treeStep = nentries/50;
    for (Int_t i=entry0; i<entry1; i++) {
        t_->GetEntry(i);

        //if((i-evStart)%treeStep==0) {
        //    printf(".");
        //    fflush(stdout);
        //}

        for(int iendcap=0; iendcap<2; iendcap++) {

            int nlcts=0;
            if(iendcap==0) {
		nlcts = nlcts_p;
	    }
            else {
		nlcts = nlcts_m;
	    }

            if(nlcts<1) continue;

            std::vector<double> allphiG;
            std::vector<double> alletaG;
            std::vector<double> allzG;
            std::vector<int> allendcap, allstation, allring, allsector, allbptx, allcscid;
	    std::vector<int> diffstation;

            outfile << "==================================" << endl;
            for(int n=0; n<nlcts; n++) {

                double phiG = lct_p_gblphi[n];
                double etaG = lct_p_gbleta[n];
                double zG = lct_p_gblZ[n];
                int endcap = lct_p_endcap[n];
		int sector = lct_p_sector[n];
		int bptx = lct_p_bptx[n];
                int ring = lct_p_ring[n];
                int station = lct_p_station[n];
		int cscid = lct_p_cscid[n];
		int strip = lct_p_strip[n];
		int keywire = lct_p_keywire[n];

                if(iendcap==1) {
                    phiG = lct_m_gblphi[n];
                    etaG = lct_m_gbleta[n];
                    zG = lct_m_gblZ[n];
                    endcap = lct_m_endcap[n];
		    sector = lct_m_sector[n];
		    bptx = lct_m_bptx[n];
                    ring = lct_m_ring[n];
                    station = lct_m_station[n];
		    cscid = lct_m_cscid[n];
		    strip = lct_m_strip[n];
		    keywire = lct_m_keywire[n];
                }

		if(bptx!=6) continue;


		if(station==0 && ring==1){
			//wire 12-15 is shared: 11-14 in software
			if(strip>=127 && keywire>14) continue;
			if(strip<127 && keywire<11) continue;
		}


		auto it = std::find(diffstation.begin(), diffstation.end(), station);
		if (it == diffstation.end()) {
			diffstation.push_back(station);
		}

                //outfile << "phiG: "<< phiG  << " etaG: " << etaG << endl;
                allphiG.push_back(phiG);
                alletaG.push_back(etaG);
                allzG.push_back(zG);
                allendcap.push_back(endcap);
		allsector.push_back(sector);
		allbptx.push_back(bptx);
                allstation.push_back(station);
                allring.push_back(ring);
		allcscid.push_back(cscid);
            }
            //outfile << "------------------" << endl;

	    bool ifOneHitStat(diffstation.size()==allstation.size());
	    if(!ifOneHitStat) continue;
	    if(diffstation.size()<2) continue;

            //if(allphiG.size()<5) continue;
            std::vector<CSCTF_t> deltaPhiEtas = dofit(allphiG,alletaG,allzG,allendcap,allsector,allstation,allring,allcscid);


            for(size_t j=0; j<deltaPhiEtas.size(); j++) {
                int status = deltaPhiEtas[j].status;
                if(status!=0) continue;
                float d_phi = deltaPhiEtas[j].dphi;
                float d_eta = deltaPhiEtas[j].deta;

                float rawphi = deltaPhiEtas[j].rawphi;
                float raweta = deltaPhiEtas[j].raweta;
                float rawZ = deltaPhiEtas[j].rawZ;

                float calphi = deltaPhiEtas[j].calphi;
                float caleta = deltaPhiEtas[j].caleta;

                int station = deltaPhiEtas[j].station;
                int endcap = deltaPhiEtas[j].endcap;
		int sector = deltaPhiEtas[j].sector;
                int ring   = deltaPhiEtas[j].ring;
		int cscid = deltaPhiEtas[j].cscid;

                float p0 = deltaPhiEtas[j].p0;
                float p1 = deltaPhiEtas[j].p1;
                float p2 = deltaPhiEtas[j].p2;
                float p3 = deltaPhiEtas[j].p3;
                float distance = deltaPhiEtas[j].distance;
		float dr0 = deltaPhiEtas[j].dr0;
		float chi2 = deltaPhiEtas[j].chi2;

		//if(chi2>0.1) continue;
                //if( fabs(p0)>200 || fabs(p1)>200 || fabs(p2)>200 || fabs(p3)>200 || distance>50 || chi2>0.1 || dr0 > 500) continue;

                //if( fabs(p0)>200 || fabs(p1)>200 || fabs(p2)>200 || fabs(p3)>200 || chi2>0.5 ||  dr0 > 500) continue;
		if( chi2>0.05 || dr0>500) continue;




                if(iendcap==1) {
                    h_m_p0 ->Fill(p0);
                    h_m_p1 ->Fill(p1);
                    h_m_p2 ->Fill(p2);
                    h_m_p3 ->Fill(p3);
                    h_m_distance->Fill(distance);
		    h_m_dr0 -> Fill(dr0);
		    h_m_chi2 -> Fill(chi2);
                } else if(iendcap==0) {
                    h_p_p0 ->Fill(p0);
                    h_p_p1 ->Fill(p1);
                    h_p_p2 ->Fill(p2);
                    h_p_p3 ->Fill(p3);
                    h_p_distance->Fill(distance);
		    h_p_dr0 -> Fill(dr0);
		    h_p_chi2 -> Fill(chi2);
                }

		outfile << " endcap: " << endcap
			<< " station: " << station
			<< " ring: " << ring
			<< " sector: " << sector
			<< " cscid: " << cscid
			<< endl;


                //if(fabs(d_phi)>2)
                outfile << "status: " << std::setprecision(1) << status << "\t"
                        << " rawphi: " << std::setprecision(5) << rawphi << "\t"
                        << " calphi: " << std::setprecision(5) << calphi << "\t"
                        << " d_phi: " << std::setprecision(5) << d_phi << "\t"
                        << " raweta: " << std::setprecision(5) << raweta << "\t"
                        << " caleta: " << std::setprecision(5) << caleta << "\t"
                        << " d_eta: " << std::setprecision(5) << d_eta << "\t"
                        << " rawZ: " << std::setprecision(5) << rawZ << "\t"
                        << " p0: " << std::setprecision(5) << p0
                        << " p1: " << std::setprecision(5) << p1
                        << " p2: " << std::setprecision(5) << p2
                        << " p3: " << std::setprecision(5) << p3
                        << " distance: " << std::setprecision(5) << distance
                        ;


                h_deltaPhi->Fill(d_phi);
                h_deltaEta->Fill(d_eta);

		if(endcap==0 && sector==0){ // plus endcap
			h_deltaPhiMEpSec0 -> Fill(d_phi);
			h_deltaEtaMEpSec0 -> Fill(d_eta);
		}
		if(endcap==0 && sector==1){ // plus endcap
			h_deltaPhiMEpSec1 -> Fill(d_phi);
			h_deltaEtaMEpSec1 -> Fill(d_eta);
		}
		if(endcap==0 && sector==2){ // plus endcap
			h_deltaPhiMEpSec2 -> Fill(d_phi);
			h_deltaEtaMEpSec2 -> Fill(d_eta);
		}
		if(endcap==0 && sector==3){ // plus endcap
			h_deltaPhiMEpSec3 -> Fill(d_phi);
			h_deltaEtaMEpSec3 -> Fill(d_eta);
		}
		if(endcap==0 && sector==4){ // plus endcap
			h_deltaPhiMEpSec4 -> Fill(d_phi);
			h_deltaEtaMEpSec4 -> Fill(d_eta);
		}
		if(endcap==0 && sector==5){ // plus endcap
			h_deltaPhiMEpSec5 -> Fill(d_phi);
			h_deltaEtaMEpSec5 -> Fill(d_eta);
		}


		if(endcap==1 && sector==0){ // minus endcap
			h_deltaPhiMEmSec0 -> Fill(d_phi);
			h_deltaEtaMEmSec0 -> Fill(d_eta);
		}
		if(endcap==1 && sector==1){ // minus endcap
			h_deltaPhiMEmSec1 -> Fill(d_phi);
			h_deltaEtaMEmSec1 -> Fill(d_eta);
		}
		if(endcap==1 && sector==2){ // minus endcap
			h_deltaPhiMEmSec2 -> Fill(d_phi);
			h_deltaEtaMEmSec2 -> Fill(d_eta);
		}
		if(endcap==1 && sector==3){ // minus endcap
			h_deltaPhiMEmSec3 -> Fill(d_phi);
			h_deltaEtaMEmSec3 -> Fill(d_eta);
		}
		if(endcap==1 && sector==4){ // minus endcap
			h_deltaPhiMEmSec4 -> Fill(d_phi);
			h_deltaEtaMEmSec4 -> Fill(d_eta);
		}
		if(endcap==1 && sector==5){ // minus endcap
			h_deltaPhiMEmSec5 -> Fill(d_phi);
			h_deltaEtaMEmSec5 -> Fill(d_eta);
		}



		if(endcap==0){
			h_deltaEtavsPhi_MEp->Fill(rawphi,d_eta);
			h_deltaPhivsEta_MEp->Fill(raweta,d_phi);

			h_dEtavs10degPhi_MEp -> Fill(rawphi*180./M_PI,d_eta);
			h_dPhivs10degPhi_MEp -> Fill(rawphi*180./M_PI,d_phi);
		}
		if(endcap==1){
			h_deltaEtavsPhi_MEm->Fill(rawphi,d_eta);
			h_deltaPhivsEta_MEm->Fill(raweta,d_phi);

                        h_dEtavs10degPhi_MEm -> Fill(rawphi*180./M_PI,d_eta);
                        h_dPhivs10degPhi_MEm -> Fill(rawphi*180./M_PI,d_phi);
		}




                //ME1/1
                if (station == 0 && ring == 1) {
                    if(endcap == 0) {
                        outfile << " ME_P_1_1" << endl;
                        h_deltaPhiMEp11  -> Fill(d_phi);
                        h_deltaEtaMEp11  -> Fill(d_eta);
			//shift phi scale by 5 degree, ME1/1/1 starts from 355 degree, ME1/1/2 starts from 5 degree
			double newphi = calphi*180./M_PI;
			newphi += 5;
			if(newphi>=360) newphi -= 360;

			h_dEtavs10degPhi_MEp11 -> Fill(newphi,d_eta);
			h_dPhivs10degPhi_MEp11 -> Fill(newphi,d_phi);
			h_dEtavsFitEta_MEp11 -> Fill(caleta,d_eta);
			h_dEtavsLUTEta_MEp11 -> Fill(raweta,d_eta);
			if(fabs(raweta)>2.05) {
				h_deltaPhiMEp11a  -> Fill(d_phi);
				h_deltaEtaMEp11a  -> Fill(d_eta);
                        	h_dEtavs10degPhi_MEp11a -> Fill(newphi,d_eta);
                        	h_dPhivs10degPhi_MEp11a -> Fill(newphi,d_phi);
				h_dEtavsFitEta_MEp11a -> Fill(caleta,d_eta);
				h_dEtavsLUTEta_MEp11a -> Fill(raweta,d_eta);
			}
			else{
				h_deltaPhiMEp11b  -> Fill(d_phi);
				h_deltaEtaMEp11b  -> Fill(d_eta);
                        	h_dEtavs10degPhi_MEp11b -> Fill(newphi,d_eta);
                        	h_dPhivs10degPhi_MEp11b -> Fill(newphi,d_phi);
				h_dEtavsFitEta_MEp11b -> Fill(caleta,d_eta);
				h_dEtavsLUTEta_MEp11b -> Fill(raweta,d_eta);
			}
                    }
                    if(endcap == 1) {
                        outfile << " ME_M_1_1" << endl;
                        h_deltaPhiMEm11 -> Fill(d_phi);
                        h_deltaEtaMEm11 -> Fill(d_eta);
                        //shift phi scale by 5 degree, ME1/1/1 starts from 355 degree, ME1/1/2 starts from 5 degree
                        double newphi = calphi*180./M_PI;
                        newphi += 5;
                        if(newphi>=360) newphi -= 360;
                        h_dEtavs10degPhi_MEm11 -> Fill(newphi,d_eta);
                        h_dPhivs10degPhi_MEm11 -> Fill(newphi,d_phi);
			h_dEtavsFitEta_MEm11 -> Fill(caleta,d_eta);
			h_dEtavsLUTEta_MEm11 -> Fill(raweta,d_eta);
                        if(fabs(raweta)>2.05) {
                                h_deltaPhiMEm11a  -> Fill(d_phi);
                                h_deltaEtaMEm11a  -> Fill(d_eta);
                                h_dEtavs10degPhi_MEm11a -> Fill(newphi,d_eta);
                                h_dPhivs10degPhi_MEm11a -> Fill(newphi,d_phi);
				h_dEtavsFitEta_MEm11a -> Fill(caleta,d_eta);
				h_dEtavsLUTEta_MEm11a -> Fill(raweta,d_eta);
			}
                        else{
                                h_deltaPhiMEm11b  -> Fill(d_phi);
                                h_deltaEtaMEm11b  -> Fill(d_eta);
                                h_dEtavs10degPhi_MEm11b -> Fill(newphi,d_eta);
                                h_dPhivs10degPhi_MEm11b -> Fill(newphi,d_phi);
				h_dEtavsFitEta_MEm11b -> Fill(caleta,d_eta);
				h_dEtavsLUTEta_MEm11b -> Fill(raweta,d_eta);
                        }
                    }
                }
                //ME1/2
                if (station == 0 && ring == 2) {
                    if(endcap == 0) {
                        outfile << " ME_P_1_2" << endl;
                        h_deltaPhiMEp12  -> Fill(d_phi);
                        h_deltaEtaMEp12  -> Fill(d_eta);
                        h_dEtavs10degPhi_MEp12 -> Fill(calphi*180./M_PI,d_eta);
                        h_dPhivs10degPhi_MEp12 -> Fill(calphi*180./M_PI,d_phi);
			h_dEtavsFitEta_MEp12 -> Fill(caleta,d_eta);
			h_dEtavsLUTEta_MEp12 -> Fill(raweta,d_eta);
                    }
                    if(endcap == 1) {
                        outfile << " ME_M_1_2" << endl;
                        h_deltaPhiMEm12 -> Fill(d_phi);
                        h_deltaEtaMEm12 -> Fill(d_eta);
                        h_dEtavs10degPhi_MEm12 -> Fill(calphi*180./M_PI,d_eta);
                        h_dPhivs10degPhi_MEm12 -> Fill(calphi*180./M_PI,d_phi);
			h_dEtavsFitEta_MEm12 -> Fill(caleta,d_eta);
			h_dEtavsLUTEta_MEm12 -> Fill(raweta,d_eta);
                    }
                }
                //ME1/3
                if (station == 0 && ring == 3) {
                    if(endcap == 0) {
                        outfile << " ME_P_1_3" << endl;
                        h_deltaPhiMEp13  -> Fill(d_phi);
                        h_deltaEtaMEp13  -> Fill(d_eta);
                        h_dEtavs10degPhi_MEp13 -> Fill(calphi*180./M_PI,d_eta);
                        h_dPhivs10degPhi_MEp13 -> Fill(calphi*180./M_PI,d_phi);
			h_dEtavsFitEta_MEp13 -> Fill(caleta,d_eta);
			h_dEtavsLUTEta_MEp13 -> Fill(raweta,d_eta);
                    }
                    if(endcap == 1) {
                        outfile << " ME_M_1_3" << endl;
                        h_deltaPhiMEm13 -> Fill(d_phi);
                        h_deltaEtaMEm13 -> Fill(d_eta);
                        h_dEtavs10degPhi_MEm13 -> Fill(calphi*180./M_PI,d_eta);
                        h_dPhivs10degPhi_MEm13 -> Fill(calphi*180./M_PI,d_phi);
			h_dEtavsFitEta_MEm13 -> Fill(caleta,d_eta);
			h_dEtavsLUTEta_MEm13 -> Fill(raweta,d_eta);
                    }
                }

                //ME2/1
                if (station == 1 && ring == 1) {
                    if(endcap == 0) {
                        outfile << " ME_P_2_1" << endl;
                        h_deltaPhiMEp21  -> Fill(d_phi);
                        h_deltaEtaMEp21  -> Fill(d_eta);
                        h_dEtavs10degPhi_MEp21 -> Fill(calphi*180./M_PI,d_eta);
                        h_dPhivs10degPhi_MEp21 -> Fill(calphi*180./M_PI,d_phi);
			h_dEtavsFitEta_MEp21 -> Fill(caleta,d_eta);
			h_dEtavsLUTEta_MEp21 -> Fill(raweta,d_eta);
                    }
                    if(endcap == 1) {
                        outfile << " ME_M_2_1" << endl;
                        h_deltaPhiMEm21 -> Fill(d_phi);
                        h_deltaEtaMEm21 -> Fill(d_eta);
                        h_dEtavs10degPhi_MEm21 -> Fill(calphi*180./M_PI,d_eta);
                        h_dPhivs10degPhi_MEm21 -> Fill(calphi*180./M_PI,d_phi);
			h_dEtavsFitEta_MEm21 -> Fill(caleta,d_eta);
			h_dEtavsLUTEta_MEm21 -> Fill(raweta,d_eta);
                    }
                }
                //ME2/2
                if (station == 1 && ring == 2) {
                    if(endcap == 0) {
                        outfile << " ME_P_2_2" << endl;
                        h_deltaPhiMEp22  -> Fill(d_phi);
                        h_deltaEtaMEp22  -> Fill(d_eta);
                        h_dEtavs10degPhi_MEp22 -> Fill(calphi*180./M_PI,d_eta);
                        h_dPhivs10degPhi_MEp22 -> Fill(calphi*180./M_PI,d_phi);
			h_dEtavsFitEta_MEp22 -> Fill(caleta,d_eta);
			h_dEtavsLUTEta_MEp22 -> Fill(raweta,d_eta);
                    }
                    if(endcap == 1) {
                        outfile << " ME_M_2_2" << endl;
                        h_deltaPhiMEm22 -> Fill(d_phi);
                        h_deltaEtaMEm22 -> Fill(d_eta);
                        h_dEtavs10degPhi_MEm22 -> Fill(calphi*180./M_PI,d_eta);
                        h_dPhivs10degPhi_MEm22 -> Fill(calphi*180./M_PI,d_phi);
			h_dEtavsFitEta_MEm22 -> Fill(caleta,d_eta);
			h_dEtavsLUTEta_MEm22 -> Fill(raweta,d_eta);
                    }
                }
                //ME3/1
                if (station == 2 && ring == 1) {
                    if(endcap == 0) {
                        outfile << " ME_P_3_1" << endl;
                        h_deltaPhiMEp31  -> Fill(d_phi);
                        h_deltaEtaMEp31  -> Fill(d_eta);
                        h_dEtavs10degPhi_MEp31 -> Fill(calphi*180./M_PI,d_eta);
                        h_dPhivs10degPhi_MEp31 -> Fill(calphi*180./M_PI,d_phi);
			h_dEtavsFitEta_MEp31 -> Fill(caleta,d_eta);
			h_dEtavsLUTEta_MEp31 -> Fill(raweta,d_eta);
                    }
                    if(endcap == 1) {
                        outfile << " ME_M_3_1" << endl;
                        h_deltaPhiMEm31 -> Fill(d_phi);
                        h_deltaEtaMEm31 -> Fill(d_eta);
                        h_dEtavs10degPhi_MEm31 -> Fill(calphi*180./M_PI,d_eta);
                        h_dPhivs10degPhi_MEm31 -> Fill(calphi*180./M_PI,d_phi);
			h_dEtavsFitEta_MEm31 -> Fill(caleta,d_eta);
			h_dEtavsLUTEta_MEm31 -> Fill(raweta,d_eta);
                    }
                }
                //ME3/2
                if (station == 2 && ring == 2) {
                    if(endcap == 0) {
                        outfile << " ME_P_3_2" << endl;
                        h_deltaPhiMEp32  -> Fill(d_phi);
                        h_deltaEtaMEp32  -> Fill(d_eta);
                        h_dEtavs10degPhi_MEp32 -> Fill(calphi*180./M_PI,d_eta);
                        h_dPhivs10degPhi_MEp32 -> Fill(calphi*180./M_PI,d_phi);
			h_dEtavsFitEta_MEp32 -> Fill(caleta,d_eta);
			h_dEtavsLUTEta_MEp32 -> Fill(raweta,d_eta);
                    }
                    if(endcap == 1) {
                        outfile << " ME_M_3_2" << endl;
                        h_deltaPhiMEm32 -> Fill(d_phi);
                        h_deltaEtaMEm32 -> Fill(d_eta);
                        h_dEtavs10degPhi_MEm32 -> Fill(calphi*180./M_PI,d_eta);
                        h_dPhivs10degPhi_MEm32 -> Fill(calphi*180./M_PI,d_phi);
			h_dEtavsFitEta_MEm32 -> Fill(caleta,d_eta);
			h_dEtavsLUTEta_MEm32 -> Fill(raweta,d_eta);
                    }
                }
                //ME4/1
                if (station == 3 && ring == 1) {
                    if(endcap == 0) {
                        outfile << " ME_P_4_1" << endl;
                        h_deltaPhiMEp41  -> Fill(d_phi);
                        h_deltaEtaMEp41  -> Fill(d_eta);
                        h_dEtavs10degPhi_MEp41 -> Fill(calphi*180./M_PI,d_eta);
                        h_dPhivs10degPhi_MEp41 -> Fill(calphi*180./M_PI,d_phi);
			h_dEtavsFitEta_MEp41 -> Fill(caleta,d_eta);
			h_dEtavsLUTEta_MEp41 -> Fill(raweta,d_eta);
                    }
                    if(endcap == 1) {
                        outfile << " ME_M_4_1" << endl;
                        h_deltaPhiMEm41 -> Fill(d_phi);
                        h_deltaEtaMEm41 -> Fill(d_eta);
                        h_dEtavs10degPhi_MEm41 -> Fill(calphi*180./M_PI,d_eta);
                        h_dPhivs10degPhi_MEm41 -> Fill(calphi*180./M_PI,d_phi);
			h_dEtavsFitEta_MEm41 -> Fill(caleta,d_eta);
			h_dEtavsLUTEta_MEm41 -> Fill(raweta,d_eta);
                    }
                }
                //ME4/2
                if (station == 3 && ring == 2) {
                    if(endcap == 0) {
                        outfile << " ME_P_4_2" << endl;
                        h_deltaPhiMEp42  -> Fill(d_phi);
                        h_deltaEtaMEp42  -> Fill(d_eta);
                        h_dEtavs10degPhi_MEp42 -> Fill(calphi*180./M_PI,d_eta);
                        h_dPhivs10degPhi_MEp42 -> Fill(calphi*180./M_PI,d_phi);
			h_dEtavsFitEta_MEp42 -> Fill(caleta,d_eta);
			h_dEtavsLUTEta_MEp42 -> Fill(raweta,d_eta);
                    }
                    if(endcap == 1) {
                        outfile << " ME_M_4_2" << endl;
                        h_deltaPhiMEm42 -> Fill(d_phi);
                        h_deltaEtaMEm42 -> Fill(d_eta);
                        h_dEtavs10degPhi_MEm42 -> Fill(calphi*180./M_PI,d_eta);
                        h_dPhivs10degPhi_MEm42 -> Fill(calphi*180./M_PI,d_phi);
			h_dEtavsFitEta_MEm42 -> Fill(caleta,d_eta);
			h_dEtavsLUTEta_MEm42 -> Fill(raweta,d_eta);
                    }
                }


            }
            outfile << "==================================\n";
            outfile << "\n\n";



        }

        //return;
    } //n entries
    printf("\n");


    TFile *OutFile = TFile::Open(_output_, "RECREATE");
    h_deltaPhi->Write();
    h_deltaEta->Write();


    h_deltaEtavsPhi_MEp->Write();
    h_deltaEtavsPhi_MEm->Write();
    h_deltaPhivsEta_MEp->Write();
    h_deltaPhivsEta_MEm->Write();



    h_deltaPhiMEp11->Write();
    h_deltaPhiMEp11a->Write();
    h_deltaPhiMEp11b->Write();
    h_deltaPhiMEp12->Write();
    h_deltaPhiMEp13->Write();
    h_deltaPhiMEp21->Write();
    h_deltaPhiMEp22->Write();
    h_deltaPhiMEp31->Write();
    h_deltaPhiMEp32->Write();
    h_deltaPhiMEp41->Write();
    h_deltaPhiMEp42->Write();

    h_deltaPhiMEm11->Write();
    h_deltaPhiMEm11a->Write();
    h_deltaPhiMEm11b->Write();
    h_deltaPhiMEm12->Write();
    h_deltaPhiMEm13->Write();
    h_deltaPhiMEm21->Write();
    h_deltaPhiMEm22->Write();
    h_deltaPhiMEm31->Write();
    h_deltaPhiMEm32->Write();
    h_deltaPhiMEm41->Write();
    h_deltaPhiMEm42->Write();


    h_deltaEtaMEp11->Write();
    h_deltaEtaMEp11a->Write();
    h_deltaEtaMEp11b->Write();
    h_deltaEtaMEp12->Write();
    h_deltaEtaMEp13->Write();
    h_deltaEtaMEp21->Write();
    h_deltaEtaMEp22->Write();
    h_deltaEtaMEp31->Write();
    h_deltaEtaMEp32->Write();
    h_deltaEtaMEp41->Write();
    h_deltaEtaMEp42->Write();

    h_deltaEtaMEm11->Write();
    h_deltaEtaMEm11a->Write();
    h_deltaEtaMEm11b->Write();
    h_deltaEtaMEm12->Write();
    h_deltaEtaMEm13->Write();
    h_deltaEtaMEm21->Write();
    h_deltaEtaMEm22->Write();
    h_deltaEtaMEm31->Write();
    h_deltaEtaMEm32->Write();
    h_deltaEtaMEm41->Write();
    h_deltaEtaMEm42->Write();

    h_m_p0->Write();
    h_m_p1->Write();
    h_m_p2->Write();
    h_m_p3->Write();
    h_m_distance->Write();

    h_p_p0->Write();
    h_p_p1->Write();
    h_p_p2->Write();
    h_p_p3->Write();
    h_p_distance->Write();

    h_deltaPhiMEpSec0->Write();
    h_deltaPhiMEpSec1->Write();
    h_deltaPhiMEpSec2->Write();
    h_deltaPhiMEpSec3->Write();
    h_deltaPhiMEpSec4->Write();
    h_deltaPhiMEpSec5->Write();

    h_deltaPhiMEmSec0->Write();
    h_deltaPhiMEmSec1->Write();
    h_deltaPhiMEmSec2->Write();
    h_deltaPhiMEmSec3->Write();
    h_deltaPhiMEmSec4->Write();
    h_deltaPhiMEmSec5->Write();


    h_deltaEtaMEpSec0->Write();
    h_deltaEtaMEpSec1->Write();
    h_deltaEtaMEpSec2->Write();
    h_deltaEtaMEpSec3->Write();
    h_deltaEtaMEpSec4->Write();
    h_deltaEtaMEpSec5->Write();

    h_deltaEtaMEmSec0->Write();
    h_deltaEtaMEmSec1->Write();
    h_deltaEtaMEmSec2->Write();
    h_deltaEtaMEmSec3->Write();
    h_deltaEtaMEmSec4->Write();
    h_deltaEtaMEmSec5->Write();

    h_m_dr0->Write();
    h_p_dr0->Write();

    h_dEtavs10degPhi_MEp->Write();
    h_dPhivs10degPhi_MEp->Write();
    h_dEtavs10degPhi_MEm->Write();
    h_dPhivs10degPhi_MEm->Write();



    h_dEtavs10degPhi_MEp11->Write();
    h_dEtavs10degPhi_MEp11a->Write();
    h_dEtavs10degPhi_MEp11b->Write();
    h_dEtavs10degPhi_MEp12->Write();
    h_dEtavs10degPhi_MEp13->Write();
    h_dEtavs10degPhi_MEp21->Write();
    h_dEtavs10degPhi_MEp22->Write();
    h_dEtavs10degPhi_MEp31->Write();
    h_dEtavs10degPhi_MEp32->Write();
    h_dEtavs10degPhi_MEp41->Write();
    h_dEtavs10degPhi_MEp42->Write();

    h_dEtavs10degPhi_MEm11->Write();
    h_dEtavs10degPhi_MEm11a->Write();
    h_dEtavs10degPhi_MEm11b->Write();
    h_dEtavs10degPhi_MEm12->Write();
    h_dEtavs10degPhi_MEm13->Write();
    h_dEtavs10degPhi_MEm21->Write();
    h_dEtavs10degPhi_MEm22->Write();
    h_dEtavs10degPhi_MEm31->Write();
    h_dEtavs10degPhi_MEm32->Write();
    h_dEtavs10degPhi_MEm41->Write();
    h_dEtavs10degPhi_MEm42->Write();


    h_dEtavsFitEta_MEp11->Write();
    h_dEtavsFitEta_MEp11a->Write();
    h_dEtavsFitEta_MEp11b->Write();
    h_dEtavsFitEta_MEp12->Write();
    h_dEtavsFitEta_MEp13->Write();
    h_dEtavsFitEta_MEp21->Write();
    h_dEtavsFitEta_MEp22->Write();
    h_dEtavsFitEta_MEp31->Write();
    h_dEtavsFitEta_MEp32->Write();
    h_dEtavsFitEta_MEp41->Write();
    h_dEtavsFitEta_MEp42->Write();
    h_dEtavsFitEta_MEm11->Write();
    h_dEtavsFitEta_MEm11a->Write();
    h_dEtavsFitEta_MEm11b->Write();
    h_dEtavsFitEta_MEm12->Write();
    h_dEtavsFitEta_MEm13->Write();
    h_dEtavsFitEta_MEm21->Write();
    h_dEtavsFitEta_MEm22->Write();
    h_dEtavsFitEta_MEm31->Write();
    h_dEtavsFitEta_MEm32->Write();
    h_dEtavsFitEta_MEm41->Write();
    h_dEtavsFitEta_MEm42->Write();


    h_dEtavsLUTEta_MEp11->Write();
    h_dEtavsLUTEta_MEp11a->Write();
    h_dEtavsLUTEta_MEp11b->Write();
    h_dEtavsLUTEta_MEp12->Write();
    h_dEtavsLUTEta_MEp13->Write();
    h_dEtavsLUTEta_MEp21->Write();
    h_dEtavsLUTEta_MEp22->Write();
    h_dEtavsLUTEta_MEp31->Write();
    h_dEtavsLUTEta_MEp32->Write();
    h_dEtavsLUTEta_MEp41->Write();
    h_dEtavsLUTEta_MEp42->Write();
    h_dEtavsLUTEta_MEm11->Write();
    h_dEtavsLUTEta_MEm11a->Write();
    h_dEtavsLUTEta_MEm11b->Write();
    h_dEtavsLUTEta_MEm12->Write();
    h_dEtavsLUTEta_MEm13->Write();
    h_dEtavsLUTEta_MEm21->Write();
    h_dEtavsLUTEta_MEm22->Write();
    h_dEtavsLUTEta_MEm31->Write();
    h_dEtavsLUTEta_MEm32->Write();
    h_dEtavsLUTEta_MEm41->Write();
    h_dEtavsLUTEta_MEm42->Write();





    h_dPhivs10degPhi_MEp11->Write();
    h_dPhivs10degPhi_MEp11a->Write();
    h_dPhivs10degPhi_MEp11b->Write();
    h_dPhivs10degPhi_MEp12->Write();
    h_dPhivs10degPhi_MEp13->Write();
    h_dPhivs10degPhi_MEp21->Write();
    h_dPhivs10degPhi_MEp22->Write();
    h_dPhivs10degPhi_MEp31->Write();
    h_dPhivs10degPhi_MEp32->Write();
    h_dPhivs10degPhi_MEp41->Write();
    h_dPhivs10degPhi_MEp42->Write();

    h_dPhivs10degPhi_MEm11->Write();
    h_dPhivs10degPhi_MEm11a->Write();
    h_dPhivs10degPhi_MEm11b->Write();
    h_dPhivs10degPhi_MEm12->Write();
    h_dPhivs10degPhi_MEm13->Write();
    h_dPhivs10degPhi_MEm21->Write();
    h_dPhivs10degPhi_MEm22->Write();
    h_dPhivs10degPhi_MEm31->Write();
    h_dPhivs10degPhi_MEm32->Write();
    h_dPhivs10degPhi_MEm41->Write();
    h_dPhivs10degPhi_MEm42->Write();



    h_m_chi2->Write();
    h_p_chi2->Write();



    outfile.close();

    OutFile->Close();
    InFile->Close();
}

