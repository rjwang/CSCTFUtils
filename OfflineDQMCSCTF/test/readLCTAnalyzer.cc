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

double distance2(float x,float y,float z, double *par);

double deltaPhi(double phi1, double phi2)
{
    double result = phi1 - phi2;
    while (result > M_PI) result -= 2*M_PI;
    while (result <= -M_PI) result += 2*M_PI;
    return result;
}

double getDeltaPhi(double *par, double phi, double eta, double z)
{
    double theta = 2*TMath::ATan(exp(-1.*eta));
    double y = z*TMath::Tan(theta);
    double t = (y - par[2])/par[3];
    double x = par[0] + par[1]*t;

    double cal_phi = TMath::ATan( y/x );
    if(cal_phi<0) cal_phi *= -1;
    double det_phi = deltaPhi(cal_phi,phi);

    cout << "calphi: " << cal_phi << " phiG: " << phi << " Det: " << det_phi << endl;

    return det_phi;
}

double getDetR(double *par, double phi, double eta, double z)
{
    double theta = 2*TMath::ATan(exp(-1.*eta));
    double y = z*TMath::Tan(theta);
    double x = y/TMath::Tan(phi);
    double deltaR = distance2(x,y,z,par);
    return sqrt(deltaR);
}

double getDeltaEta(double *par, double phi, double eta, double z)
{
    double tanphi = TMath::Tan(phi);
    double t = (par[0]*tanphi-par[2])/(par[3]-par[1]*tanphi);
    double y = par[2] + par[3]*t;

    double theta = TMath::ATan( y/z );
    double cal_eta = (-1.)*TMath::Log10(fabs(TMath::Tan(theta/2.)))/TMath::LogE();
    if(theta<0) cal_eta *= -1;

    cout << "theta: " << theta << " caleta: " << cal_eta << " etaG: " << eta << " DetEta: " << cal_eta-eta << endl;

    return cal_eta-eta;
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

std::vector<std::pair<double,double> > dofit(std::vector<double> allphiG, std::vector<double> alletaG, std::vector<double> allzG)
{
    std::vector<std::pair<double,double> > LUTrels;
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

	    double theta = 2*TMath::ATan(exp(-1.*etaG));
	    //double r = zG/TMath::Cos(theta);

	    //double x = r*TMath::Sin(theta)*TMath::Cos(phiG);
	    //double y = r*TMath::Sin(theta)*TMath::Sin(phiG);
	    double z = zG;
	    double y = z*TMath::Tan(theta);
	    double x = y/TMath::Tan(phiG);

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
        min->ExecuteCommand("MIGRAD",arglist,2);

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

        double delta_phi = getDeltaPhi( parFit,allphiG[lct],alletaG[lct],allzG[lct] );
        double delta_eta = getDeltaEta( parFit,allphiG[lct],alletaG[lct],allzG[lct] );

        //double delta_R = getDetR( parFit,allphiG[lct],alletaG[lct],allzG[lct] );

        cout << "delta phi: " << delta_phi << " delta_eta: " << delta_eta << endl;

        std::pair<double,double> deltaLUT(delta_phi,delta_eta);
	//std::pair<double,double> deltaLUT(delta_R,delta_eta);
        LUTrels.push_back(deltaLUT);

    }

    return LUTrels;
}

void readLCTAnalyzer(TString _input_="./csctf_test.root", TString _output_="output_test.root")
{

    TFile *InFile = TFile::Open(_input_, "READ");
    TTree *t_ = (TTree *) InFile->Get("OfflineDQMCSCTF/data");


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
    TH1F *h_deltaPhi = new TH1F("d_phi",";#phi - #phi_{fit};Events",500,-5.,5.);
    TH1F *h_deltaR = new TH1F("d_R",";#Delta R;Events",500,-5000.,5000.);
    TH1F *h_deltaEta = new TH1F("d_eta",";#eta - #eta_{fit};Events",500,-2.5,2.5);


    TH1F* h_deltaPhiMEp11 = new TH1F("h_deltaPhiMEp11",";#phi - #phi_{fit};Events",500,-5.,5.);
    TH1F* h_deltaPhiMEp12 = new TH1F("h_deltaPhiMEp12",";#phi - #phi_{fit};Events",500,-5.,5.);
    TH1F* h_deltaPhiMEp13 = new TH1F("h_deltaPhiMEp13",";#phi - #phi_{fit};Events",500,-5.,5.);
    TH1F* h_deltaPhiMEp21 = new TH1F("h_deltaPhiMEp21",";#phi - #phi_{fit};Events",500,-5.,5.);
    TH1F* h_deltaPhiMEp22 = new TH1F("h_deltaPhiMEp22",";#phi - #phi_{fit};Events",500,-5.,5.);
    TH1F* h_deltaPhiMEp31 = new TH1F("h_deltaPhiMEp31",";#phi - #phi_{fit};Events",500,-5.,5.);
    TH1F* h_deltaPhiMEp32 = new TH1F("h_deltaPhiMEp32",";#phi - #phi_{fit};Events",500,-5.,5.);
    TH1F* h_deltaPhiMEp41 = new TH1F("h_deltaPhiMEp41",";#phi - #phi_{fit};Events",500,-5.,5.);
    TH1F* h_deltaPhiMEp42 = new TH1F("h_deltaPhiMEp42",";#phi - #phi_{fit};Events",500,-5.,5.);

    TH1F* h_deltaPhiMEm11 = new TH1F("h_deltaPhiMEm11",";#phi - #phi_{fit};Events",500,-5.,5.);
    TH1F* h_deltaPhiMEm12 = new TH1F("h_deltaPhiMEm12",";#phi - #phi_{fit};Events",500,-5.,5.);
    TH1F* h_deltaPhiMEm13 = new TH1F("h_deltaPhiMEm13",";#phi - #phi_{fit};Events",500,-5.,5.);
    TH1F* h_deltaPhiMEm21 = new TH1F("h_deltaPhiMEm21",";#phi - #phi_{fit};Events",500,-5.,5.);
    TH1F* h_deltaPhiMEm22 = new TH1F("h_deltaPhiMEm22",";#phi - #phi_{fit};Events",500,-5.,5.);
    TH1F* h_deltaPhiMEm31 = new TH1F("h_deltaPhiMEm31",";#phi - #phi_{fit};Events",500,-5.,5.);
    TH1F* h_deltaPhiMEm32 = new TH1F("h_deltaPhiMEm32",";#phi - #phi_{fit};Events",500,-5.,5.);
    TH1F* h_deltaPhiMEm41 = new TH1F("h_deltaPhiMEm41",";#phi - #phi_{fit};Events",500,-5.,5.);
    TH1F* h_deltaPhiMEm42 = new TH1F("h_deltaPhiMEm42",";#phi - #phi_{fit};Events",500,-5.,5.);

    TH1F* h_deltaEtaMEp11 = new TH1F("h_deltaEtaMEp11",";#eta - #eta_{fit};Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEp12 = new TH1F("h_deltaEtaMEp12",";#eta - #eta_{fit};Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEp13 = new TH1F("h_deltaEtaMEp13",";#eta - #eta_{fit};Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEp21 = new TH1F("h_deltaEtaMEp21",";#eta - #eta_{fit};Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEp22 = new TH1F("h_deltaEtaMEp22",";#eta - #eta_{fit};Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEp31 = new TH1F("h_deltaEtaMEp31",";#eta - #eta_{fit};Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEp32 = new TH1F("h_deltaEtaMEp32",";#eta - #eta_{fit};Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEp41 = new TH1F("h_deltaEtaMEp41",";#eta - #eta_{fit};Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEp42 = new TH1F("h_deltaEtaMEp42",";#eta - #eta_{fit};Events",500,-2.5,2.5);

    TH1F* h_deltaEtaMEm11 = new TH1F("h_deltaEtaMEm11",";#eta - #eta_{fit};Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEm12 = new TH1F("h_deltaEtaMEm12",";#eta - #eta_{fit};Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEm13 = new TH1F("h_deltaEtaMEm13",";#eta - #eta_{fit};Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEm21 = new TH1F("h_deltaEtaMEm21",";#eta - #eta_{fit};Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEm22 = new TH1F("h_deltaEtaMEm22",";#eta - #eta_{fit};Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEm31 = new TH1F("h_deltaEtaMEm31",";#eta - #eta_{fit};Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEm32 = new TH1F("h_deltaEtaMEm32",";#eta - #eta_{fit};Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEm41 = new TH1F("h_deltaEtaMEm41",";#eta - #eta_{fit};Events",500,-2.5,2.5);
    TH1F* h_deltaEtaMEm42 = new TH1F("h_deltaEtaMEm42",";#eta - #eta_{fit};Events",500,-2.5,2.5);


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
        std::vector<double> allendcap, allstation, allring;

        for(int n=0; n<nlcts; n++) {
            double phiG = lct_gblphi[n];
            double etaG = lct_gbleta[n];
	    double zG = lct_gblZ[n];
            double endcap = lct_endcap[n];
            double ring = lct_ring[n];
            double station = lct_station[n];

	    if(zG<0) etaG *= -1.;
            allphiG.push_back(phiG);
            alletaG.push_back(etaG);
	    allzG.push_back(zG);
            allendcap.push_back(endcap);
            allstation.push_back(station);
            allring.push_back(ring);
        }

        std::vector<std::pair<double,double> > deltaPhiEtas = dofit(allphiG,alletaG,allzG);


        for(size_t j=0; j<deltaPhiEtas.size(); j++) {
            double d_phi = deltaPhiEtas[j].first;
            double d_eta = deltaPhiEtas[j].second;
            cout << "d_phi: " << d_phi << " d_eta: " << d_eta << endl;
            h_deltaPhi->Fill(d_phi);
	    //h_deltaR->Fill(d_phi);
            h_deltaEta->Fill(d_eta);

            int station = allstation[j];
            int endcap = allendcap[j];
            int ring   = allring[j];

            //ME1/1
            if (station == 0 && ring == 1) {
                if(endcap == 0) {
                    h_deltaPhiMEp11  -> Fill(d_phi);
                    h_deltaEtaMEp11  -> Fill(d_eta);
                }
                if(endcap == 1) {
                    h_deltaPhiMEm11 -> Fill(d_phi);
                    h_deltaEtaMEm11 -> Fill(d_eta);
                }
            }
            //ME1/2
            if (station == 0 && ring == 2) {
                if(endcap == 0) {
                    h_deltaPhiMEp12  -> Fill(d_phi);
                    h_deltaEtaMEp12  -> Fill(d_eta);
                }
                if(endcap == 1) {
                    h_deltaPhiMEm12 -> Fill(d_phi);
                    h_deltaEtaMEm12 -> Fill(d_eta);
                }
            }
            //ME1/3
            if (station == 0 && ring == 3) {
                if(endcap == 0) {
                    h_deltaPhiMEp13  -> Fill(d_phi);
                    h_deltaEtaMEp13  -> Fill(d_eta);
                }
                if(endcap == 1) {
                    h_deltaPhiMEm13 -> Fill(d_phi);
                    h_deltaEtaMEm13 -> Fill(d_eta);
                }
            }

            //ME2/1
            if (station == 1 && ring == 1) {
                if(endcap == 0) {
                    h_deltaPhiMEp21  -> Fill(d_phi);
                    h_deltaEtaMEp21  -> Fill(d_eta);
                }
                if(endcap == 1) {
                    h_deltaPhiMEm21 -> Fill(d_phi);
                    h_deltaEtaMEm21 -> Fill(d_eta);
                }
            }
            //ME2/2
            if (station == 1 && ring == 2) {
                if(endcap == 0) {
                    h_deltaPhiMEp22  -> Fill(d_phi);
                    h_deltaEtaMEp22  -> Fill(d_eta);
                }
                if(endcap == 1) {
                    h_deltaPhiMEm22 -> Fill(d_phi);
                    h_deltaEtaMEm22 -> Fill(d_eta);
                }
            }
            //ME3/1
            if (station == 2 && ring == 1) {
                if(endcap == 0) {
                    h_deltaPhiMEp31  -> Fill(d_phi);
                    h_deltaEtaMEp31  -> Fill(d_eta);
                }
                if(endcap == 1) {
                    h_deltaPhiMEm31 -> Fill(d_phi);
                    h_deltaEtaMEm31 -> Fill(d_eta);
                }
            }
            //ME3/2
            if (station == 2 && ring == 2) {
                if(endcap == 0) {
                    h_deltaPhiMEp32  -> Fill(d_phi);
                    h_deltaEtaMEp32  -> Fill(d_eta);
                }
                if(endcap == 1) {
                    h_deltaPhiMEm32 -> Fill(d_phi);
                    h_deltaEtaMEm32 -> Fill(d_eta);
                }
            }
            //ME4/1
            if (station == 3 && ring == 1) {
                if(endcap == 0) {
                    h_deltaPhiMEp41  -> Fill(d_phi);
                    h_deltaEtaMEp41  -> Fill(d_eta);
                }
                if(endcap == 1) {
                    h_deltaPhiMEm41 -> Fill(d_phi);
                    h_deltaEtaMEm41 -> Fill(d_eta);
                }
            }
            //ME4/2
            if (station == 3 && ring == 2) {
                if(endcap == 0) {
                    h_deltaPhiMEp42  -> Fill(d_phi);
                    h_deltaEtaMEp42  -> Fill(d_eta);
                }
                if(endcap == 1) {
                    h_deltaPhiMEm42 -> Fill(d_phi);
                    h_deltaEtaMEm42 -> Fill(d_eta);
                }
            }


        }



        //return;
    }
    printf("\n");


    TFile *OutFile = TFile::Open(_output_, "RECREATE");
    h_deltaPhi->Write();
    h_deltaR->Write();
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



    OutFile->Close();
    InFile->Close();
}

