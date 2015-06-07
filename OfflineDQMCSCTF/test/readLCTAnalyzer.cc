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
    float dphi;
    float deta;
    float rawphi,raweta,rawZ;
    int endcap,station,ring;
    float calphi,caleta;
    float p0,p1,p2,p3,dR;
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
    double cal_theta = TMath::ATan( fabs(y/z) );
    if(z<0) cal_theta = M_PI - cal_theta;
    double cal_eta = (-1.)*TMath::Log10(TMath::Tan(cal_theta/2.))/TMath::LogE();
    if(z<0) cal_eta = -1.*fabs(cal_eta);

    double det_eta = cal_eta-eta;

    myfit.deta=det_eta;
    myfit.caleta=cal_eta;

    myfit.dphi=det_phi;
    myfit.calphi=cal_phi;

    return myfit;
}

double getDetR(double *par, double phi, double eta, double z)
{
    double theta = 2*TMath::ATan(exp(-1.*eta));
    double y = z*TMath::Tan(theta);
    double x = y/TMath::Tan(phi);

    if(phi>=0                  && phi<M_PI/2. )       { x=fabs(x);            y=fabs(y); }
    else if(phi>=M_PI/2.       && phi<M_PI    )       { x=-1.*fabs(x);        y=fabs(y); }
    else if(phi>=M_PI          && phi<3.*M_PI/2.)     { x=-1.*fabs(x);        y=-1.*fabs(y); }
    else if(phi>=3.*M_PI/2.    && phi<2.*M_PI )       { x=fabs(x);            y=-1.*fabs(y); }

    double deltaR = distance2(x,y,z,par);
    return sqrt(deltaR);
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
		std::vector<int> allendcap, std::vector<int> allstation, std::vector<int> allring)
{
    std::vector<CSCTF_t> LUTrels;
    for(size_t lct=0; lct<allphiG.size(); lct++) {

        std::vector<double> tofit_phiG;
        std::vector<double> tofit_etaG;
	std::vector<double> tofit_zG;

        for(size_t t=0; t<allphiG.size(); t++) {
            if(t==lct) continue;
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
	    double y = z*TMath::Tan(theta);
	    double x = y/TMath::Tan(phiG);

	    if(phiG>=0 			&& phiG<M_PI/2.	) 	{ x=fabs(x); 		y=fabs(y); }
	    else if(phiG>=M_PI/2. 	&& phiG<M_PI	) 	{ x=-1.*fabs(x); 	y=fabs(y); }
	    else if(phiG>=M_PI 		&& phiG<3.*M_PI/2.) 	{ x=-1.*fabs(x); 	y=-1.*fabs(y); }
	    else if(phiG>=3.*M_PI/2. 	&& phiG<2.*M_PI	) 	{ x=fabs(x); 		y=-1.*fabs(y); }

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

        // get fit parameters
        double parFit[4];
        for (int j = 0; j < 4; ++j) {
            parFit[j] = min->GetParameter(j);
            //cout  << "par: " << parFit[j] << endl;
        }

        CSCTF_t DeltaPhiEta = getDeltaPhiEta( parFit,allphiG[lct],alletaG[lct],allzG[lct] );
        double delta_R = getDetR( parFit,allphiG[lct],alletaG[lct],allzG[lct] );

	CSCTF_t myfit;
	myfit.status = status;
	myfit.dphi = DeltaPhiEta.dphi;
	myfit.deta = DeltaPhiEta.deta;
	myfit.rawphi = allphiG[lct];
	myfit.raweta = alletaG[lct];
	myfit.endcap = allendcap[lct];
	myfit.station = allstation[lct];
	myfit.ring = allring[lct];
	myfit.rawZ = allzG[lct];
	myfit.calphi = DeltaPhiEta.calphi;
        myfit.caleta = DeltaPhiEta.caleta;
	myfit.p0 = parFit[0];
	myfit.p1 = parFit[1];
	myfit.p2 = parFit[2];
	myfit.p3 = parFit[3];
        myfit.dR = delta_R;

	LUTrels.push_back(myfit);
    }

    return LUTrels;
}

void readLCTAnalyzer(TString _input_="./csctf_Run246926.root", TString _output_="output_test.root", bool isplusEndCap=true)
{

    TFile *InFile = TFile::Open(_input_, "READ");
    TTree *t_ = (TTree *) InFile->Get("OfflineDQMCSCTF/data");

    ofstream outfile;
    outfile.open ("output_csctf.log");

    //gen level event
    Int_t nlcts;
    Float_t lct_gblphi[MAXLUTS], lct_gbleta[MAXLUTS], lct_gblZ[MAXLUTS];
    Int_t lct_station[MAXLUTS], lct_ring[MAXLUTS], lct_endcap[MAXLUTS];


    //mc truth
    t_->SetBranchAddress("nlcts", &nlcts);
    t_->SetBranchAddress("lct_gblphi", lct_gblphi);
    t_->SetBranchAddress("lct_gbleta", lct_gbleta);
    t_->SetBranchAddress("lct_gblZ", lct_gblZ);
    t_->SetBranchAddress("lct_station", lct_station);
    t_->SetBranchAddress("lct_ring", lct_ring);
    t_->SetBranchAddress("lct_endcap", lct_endcap);




    //hists
    TH1F *h_deltaPhi = new TH1F("d_phi",";#it{#phi}(LUT) - #it{#phi}(Fit);Events",500,-7.,7.);
    TH1F *h_dR = new TH1F("h_dR",";#Delta R;Events",500,0,2000.);
    TH1F *h_deltaEta = new TH1F("d_eta",";#it{#eta}(LUT) - #it{#eta}(Fit);Events",500,-2.5,2.5);


    TH1F* h_deltaPhiMEp11 = new TH1F("h_deltaPhiMEp11",";#it{#phi}(LUT) - #it{#phi}(Fit);Events",500,-7.,7.);
    TH1F* h_deltaPhiMEp12 = new TH1F("h_deltaPhiMEp12",";#it{#phi}(LUT) - #it{#phi}(Fit);Events",500,-7.,7.);
    TH1F* h_deltaPhiMEp13 = new TH1F("h_deltaPhiMEp13",";#it{#phi}(LUT) - #it{#phi}(Fit);Events",500,-7.,7.);
    TH1F* h_deltaPhiMEp21 = new TH1F("h_deltaPhiMEp21",";#it{#phi}(LUT) - #it{#phi}(Fit);Events",500,-7.,7.);
    TH1F* h_deltaPhiMEp22 = new TH1F("h_deltaPhiMEp22",";#it{#phi}(LUT) - #it{#phi}(Fit);Events",500,-7.,7.);
    TH1F* h_deltaPhiMEp31 = new TH1F("h_deltaPhiMEp31",";#it{#phi}(LUT) - #it{#phi}(Fit);Events",500,-7.,7.);
    TH1F* h_deltaPhiMEp32 = new TH1F("h_deltaPhiMEp32",";#it{#phi}(LUT) - #it{#phi}(Fit);Events",500,-7.,7.);
    TH1F* h_deltaPhiMEp41 = new TH1F("h_deltaPhiMEp41",";#it{#phi}(LUT) - #it{#phi}(Fit);Events",500,-7.,7.);
    TH1F* h_deltaPhiMEp42 = new TH1F("h_deltaPhiMEp42",";#it{#phi}(LUT) - #it{#phi}(Fit);Events",500,-7.,7.);

    TH1F* h_deltaPhiMEm11 = new TH1F("h_deltaPhiMEm11",";#it{#phi}(LUT) - #it{#phi}(Fit);Events",500,-7.,7.);
    TH1F* h_deltaPhiMEm12 = new TH1F("h_deltaPhiMEm12",";#it{#phi}(LUT) - #it{#phi}(Fit);Events",500,-7.,7.);
    TH1F* h_deltaPhiMEm13 = new TH1F("h_deltaPhiMEm13",";#it{#phi}(LUT) - #it{#phi}(Fit);Events",500,-7.,7.);
    TH1F* h_deltaPhiMEm21 = new TH1F("h_deltaPhiMEm21",";#it{#phi}(LUT) - #it{#phi}(Fit);Events",500,-7.,7.);
    TH1F* h_deltaPhiMEm22 = new TH1F("h_deltaPhiMEm22",";#it{#phi}(LUT) - #it{#phi}(Fit);Events",500,-7.,7.);
    TH1F* h_deltaPhiMEm31 = new TH1F("h_deltaPhiMEm31",";#it{#phi}(LUT) - #it{#phi}(Fit);Events",500,-7.,7.);
    TH1F* h_deltaPhiMEm32 = new TH1F("h_deltaPhiMEm32",";#it{#phi}(LUT) - #it{#phi}(Fit);Events",500,-7.,7.);
    TH1F* h_deltaPhiMEm41 = new TH1F("h_deltaPhiMEm41",";#it{#phi}(LUT) - #it{#phi}(Fit);Events",500,-7.,7.);
    TH1F* h_deltaPhiMEm42 = new TH1F("h_deltaPhiMEm42",";#it{#phi}(LUT) - #it{#phi}(Fit);Events",500,-7.,7.);

    TH1F* h_deltaEtaMEp11 = new TH1F("h_deltaEtaMEp11",";#it{#eta}(LUT) - #it{#eta}(Fit);Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEp12 = new TH1F("h_deltaEtaMEp12",";#it{#eta}(LUT) - #it{#eta}(Fit);Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEp13 = new TH1F("h_deltaEtaMEp13",";#it{#eta}(LUT) - #it{#eta}(Fit);Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEp21 = new TH1F("h_deltaEtaMEp21",";#it{#eta}(LUT) - #it{#eta}(Fit);Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEp22 = new TH1F("h_deltaEtaMEp22",";#it{#eta}(LUT) - #it{#eta}(Fit);Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEp31 = new TH1F("h_deltaEtaMEp31",";#it{#eta}(LUT) - #it{#eta}(Fit);Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEp32 = new TH1F("h_deltaEtaMEp32",";#it{#eta}(LUT) - #it{#eta}(Fit);Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEp41 = new TH1F("h_deltaEtaMEp41",";#it{#eta}(LUT) - #it{#eta}(Fit);Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEp42 = new TH1F("h_deltaEtaMEp42",";#it{#eta}(LUT) - #it{#eta}(Fit);Events",500,-2.5,2.5);

    TH1F* h_deltaEtaMEm11 = new TH1F("h_deltaEtaMEm11",";#it{#eta}(LUT) - #it{#eta}(Fit);Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEm12 = new TH1F("h_deltaEtaMEm12",";#it{#eta}(LUT) - #it{#eta}(Fit);Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEm13 = new TH1F("h_deltaEtaMEm13",";#it{#eta}(LUT) - #it{#eta}(Fit);Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEm21 = new TH1F("h_deltaEtaMEm21",";#it{#eta}(LUT) - #it{#eta}(Fit);Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEm22 = new TH1F("h_deltaEtaMEm22",";#it{#eta}(LUT) - #it{#eta}(Fit);Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEm31 = new TH1F("h_deltaEtaMEm31",";#it{#eta}(LUT) - #it{#eta}(Fit);Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEm32 = new TH1F("h_deltaEtaMEm32",";#it{#eta}(LUT) - #it{#eta}(Fit);Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEm41 = new TH1F("h_deltaEtaMEm41",";#it{#eta}(LUT) - #it{#eta}(Fit);Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEm42 = new TH1F("h_deltaEtaMEm42",";#it{#eta}(LUT) - #it{#eta}(Fit);Events",500,-2.5,2.5);

    TH1F* h_p0 = new TH1F("h_p0",";p_{0};Events",500,-2000,2000);
    TH1F* h_p1 = new TH1F("h_p1",";p_{1};Events",500,-2000,2000);
    TH1F* h_p2 = new TH1F("h_p2",";p_{2};Events",500,-2000,2000);
    TH1F* h_p3 = new TH1F("h_p3",";p_{3};Events",500,-2000,2000);




    //int evStart = 0;
    //int evEnd = -1;
    //printf("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");
    //printf("Scanning the ntuple :");
    // all entries and fill the histograms
    Int_t nentries = (Int_t)t_->GetEntries();
    //int treeStep = nentries/50;
    for (Int_t i=0; i<nentries; i++) {
        t_->GetEntry(i);

        //if((i-evStart)%treeStep==0) {
        //    printf(".");
        //    fflush(stdout);
        //}

        if(nlcts<1) continue;

        std::vector<double> allphiG;
        std::vector<double> alletaG;
        std::vector<double> allzG;
        std::vector<int> allendcap, allstation, allring;

	outfile << "==================================" << endl;
        for(int n=0; n<nlcts; n++) {
            double phiG = lct_gblphi[n];
            double etaG = lct_gbleta[n];
	    double zG = lct_gblZ[n];
            double endcap = lct_endcap[n];
            double ring = lct_ring[n];
            double station = lct_station[n];

	    if(!isplusEndCap && endcap==0) continue; // only minus endcap side
	    if(isplusEndCap && endcap==1) continue; // only plus endcap side

	    //outfile << "phiG: "<< phiG  << " etaG: " << etaG << endl;
            allphiG.push_back(phiG);
            alletaG.push_back(etaG);
	    allzG.push_back(zG);
            allendcap.push_back(endcap);
            allstation.push_back(station);
            allring.push_back(ring);
        }
        //outfile << "------------------" << endl;

	if(allphiG.size()<5) continue;
	std::vector<CSCTF_t> deltaPhiEtas = dofit(allphiG,alletaG,allzG,allendcap,allstation,allring);


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
            int ring   = deltaPhiEtas[j].ring;

	    float p0 = deltaPhiEtas[j].p0;
	    float p1 = deltaPhiEtas[j].p1;
	    float p2 = deltaPhiEtas[j].p2;
	    float p3 = deltaPhiEtas[j].p3;
	    float dR = deltaPhiEtas[j].dR;
	    h_p0 ->Fill(p0);
	    h_p1 ->Fill(p1);
	    h_p2 ->Fill(p2);
	    h_p3 ->Fill(p3);
	    h_dR->Fill(dR);



	    if( fabs(p0)>200 || fabs(p1)>200 || fabs(p2)>200 || fabs(p3)>200 || dR>50) continue;

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
			;


            h_deltaPhi->Fill(d_phi);
            h_deltaEta->Fill(d_eta);


            //ME1/1
            if (station == 0 && ring == 1) {
                if(endcap == 0) {
		    outfile << " ME_P_1_1" << endl;
                    h_deltaPhiMEp11  -> Fill(d_phi);
                    h_deltaEtaMEp11  -> Fill(d_eta);
                }
                if(endcap == 1) {
		    outfile << " ME_M_1_1" << endl;
                    h_deltaPhiMEm11 -> Fill(d_phi);
                    h_deltaEtaMEm11 -> Fill(d_eta);
                }
            }
            //ME1/2
            if (station == 0 && ring == 2) {
                if(endcap == 0) {
		    outfile << " ME_P_1_2" << endl;
                    h_deltaPhiMEp12  -> Fill(d_phi);
                    h_deltaEtaMEp12  -> Fill(d_eta);
                }
                if(endcap == 1) {
		    outfile << " ME_M_1_2" << endl;
                    h_deltaPhiMEm12 -> Fill(d_phi);
                    h_deltaEtaMEm12 -> Fill(d_eta);
                }
            }
            //ME1/3
            if (station == 0 && ring == 3) {
                if(endcap == 0) {
		    outfile << " ME_P_1_3" << endl;
                    h_deltaPhiMEp13  -> Fill(d_phi);
                    h_deltaEtaMEp13  -> Fill(d_eta);
                }
                if(endcap == 1) {
		    outfile << " ME_M_1_3" << endl;
                    h_deltaPhiMEm13 -> Fill(d_phi);
                    h_deltaEtaMEm13 -> Fill(d_eta);
                }
            }

            //ME2/1
            if (station == 1 && ring == 1) {
                if(endcap == 0) {
		    outfile << " ME_P_2_1" << endl;
                    h_deltaPhiMEp21  -> Fill(d_phi);
                    h_deltaEtaMEp21  -> Fill(d_eta);
                }
                if(endcap == 1) {
		    outfile << " ME_M_2_1" << endl;
                    h_deltaPhiMEm21 -> Fill(d_phi);
                    h_deltaEtaMEm21 -> Fill(d_eta);
                }
            }
            //ME2/2
            if (station == 1 && ring == 2) {
                if(endcap == 0) {
		    outfile << " ME_P_2_2" << endl;
                    h_deltaPhiMEp22  -> Fill(d_phi);
                    h_deltaEtaMEp22  -> Fill(d_eta);
                }
                if(endcap == 1) {
		    outfile << " ME_M_2_2" << endl;
                    h_deltaPhiMEm22 -> Fill(d_phi);
                    h_deltaEtaMEm22 -> Fill(d_eta);
                }
            }
            //ME3/1
            if (station == 2 && ring == 1) {
                if(endcap == 0) {
		    outfile << " ME_P_3_1" << endl;
                    h_deltaPhiMEp31  -> Fill(d_phi);
                    h_deltaEtaMEp31  -> Fill(d_eta);
                }
                if(endcap == 1) {
		    outfile << " ME_M_3_1" << endl;
                    h_deltaPhiMEm31 -> Fill(d_phi);
                    h_deltaEtaMEm31 -> Fill(d_eta);
                }
            }
            //ME3/2
            if (station == 2 && ring == 2) {
                if(endcap == 0) {
		    outfile << " ME_P_3_2" << endl;
                    h_deltaPhiMEp32  -> Fill(d_phi);
                    h_deltaEtaMEp32  -> Fill(d_eta);
                }
                if(endcap == 1) {
		    outfile << " ME_M_3_2" << endl;
                    h_deltaPhiMEm32 -> Fill(d_phi);
                    h_deltaEtaMEm32 -> Fill(d_eta);
                }
            }
            //ME4/1
            if (station == 3 && ring == 1) {
                if(endcap == 0) {
		    outfile << " ME_P_4_1" << endl;
                    h_deltaPhiMEp41  -> Fill(d_phi);
                    h_deltaEtaMEp41  -> Fill(d_eta);
                }
                if(endcap == 1) {
		    outfile << " ME_M_4_1" << endl;
                    h_deltaPhiMEm41 -> Fill(d_phi);
                    h_deltaEtaMEm41 -> Fill(d_eta);
                }
            }
            //ME4/2
            if (station == 3 && ring == 2) {
                if(endcap == 0) {
		    outfile << " ME_P_4_2" << endl;
                    h_deltaPhiMEp42  -> Fill(d_phi);
                    h_deltaEtaMEp42  -> Fill(d_eta);
                }
                if(endcap == 1) {
		    outfile << " ME_M_4_2" << endl;
                    h_deltaPhiMEm42 -> Fill(d_phi);
                    h_deltaEtaMEm42 -> Fill(d_eta);
                }
            }


        }
        outfile << "==================================\n";
	outfile << "\n\n";




        //return;
    }
    printf("\n");


    TFile *OutFile = TFile::Open(_output_, "RECREATE");
    h_deltaPhi->Write();
    h_dR->Write();
    h_deltaEta->Write();

    h_deltaPhiMEp11->Write();
    h_deltaPhiMEp12->Write();
    h_deltaPhiMEp13->Write();
    h_deltaPhiMEp21->Write();
    h_deltaPhiMEp22->Write();
    h_deltaPhiMEp31->Write();
    h_deltaPhiMEp32->Write();
    h_deltaPhiMEp41->Write();
    h_deltaPhiMEp42->Write();

    h_deltaPhiMEm11->Write();
    h_deltaPhiMEm12->Write();
    h_deltaPhiMEm13->Write();
    h_deltaPhiMEm21->Write();
    h_deltaPhiMEm22->Write();
    h_deltaPhiMEm31->Write();
    h_deltaPhiMEm32->Write();
    h_deltaPhiMEm41->Write();
    h_deltaPhiMEm42->Write();


    h_deltaEtaMEp11->Write();
    h_deltaEtaMEp12->Write();
    h_deltaEtaMEp13->Write();
    h_deltaEtaMEp21->Write();
    h_deltaEtaMEp22->Write();
    h_deltaEtaMEp31->Write();
    h_deltaEtaMEp32->Write();
    h_deltaEtaMEp41->Write();
    h_deltaEtaMEp42->Write();

    h_deltaEtaMEm11->Write();
    h_deltaEtaMEm12->Write();
    h_deltaEtaMEm13->Write();
    h_deltaEtaMEm21->Write();
    h_deltaEtaMEm22->Write();
    h_deltaEtaMEm31->Write();
    h_deltaEtaMEm32->Write();
    h_deltaEtaMEm41->Write();
    h_deltaEtaMEm42->Write();

    h_p0->Write();
    h_p1->Write();
    h_p2->Write();
    h_p3->Write();

    outfile.close();

    OutFile->Close();
    InFile->Close();
}

