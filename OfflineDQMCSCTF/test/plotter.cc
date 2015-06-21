//
//
#include "/afs/cern.ch/user/r/rewang/MainResults/AnaCode/plotPost.h"





void plotter(TString Input="output_test.root")
{

    TFile *infile = TFile::Open(Input, "READ");

    vector<TString> All1Dhists;

    All1Dhists.push_back("h_deltaPhi");
    All1Dhists.push_back("h_deltaEta");

    All1Dhists.push_back("h_deltaPhiMEp11");
    All1Dhists.push_back("h_deltaPhiMEp11a");
    All1Dhists.push_back("h_deltaPhiMEp11b");
    All1Dhists.push_back("h_deltaPhiMEp12");
    All1Dhists.push_back("h_deltaPhiMEp13");
    All1Dhists.push_back("h_deltaPhiMEp21");
    All1Dhists.push_back("h_deltaPhiMEp22");
    All1Dhists.push_back("h_deltaPhiMEp31");
    All1Dhists.push_back("h_deltaPhiMEp32");
    All1Dhists.push_back("h_deltaPhiMEp41");
    All1Dhists.push_back("h_deltaPhiMEp42");


    All1Dhists.push_back("h_deltaPhiMEm11");
    All1Dhists.push_back("h_deltaPhiMEm11a");
    All1Dhists.push_back("h_deltaPhiMEm11b");
    All1Dhists.push_back("h_deltaPhiMEm12");
    All1Dhists.push_back("h_deltaPhiMEm13");
    All1Dhists.push_back("h_deltaPhiMEm21");
    All1Dhists.push_back("h_deltaPhiMEm22");
    All1Dhists.push_back("h_deltaPhiMEm31");
    All1Dhists.push_back("h_deltaPhiMEm32");
    All1Dhists.push_back("h_deltaPhiMEm41");
    All1Dhists.push_back("h_deltaPhiMEm42");

    All1Dhists.push_back("h_deltaEtaMEp11");
    All1Dhists.push_back("h_deltaEtaMEp11a");
    All1Dhists.push_back("h_deltaEtaMEp11b");
    All1Dhists.push_back("h_deltaEtaMEp12");
    All1Dhists.push_back("h_deltaEtaMEp13");
    All1Dhists.push_back("h_deltaEtaMEp21");
    All1Dhists.push_back("h_deltaEtaMEp22");
    All1Dhists.push_back("h_deltaEtaMEp31");
    All1Dhists.push_back("h_deltaEtaMEp32");
    All1Dhists.push_back("h_deltaEtaMEp41");
    All1Dhists.push_back("h_deltaEtaMEp42");


    All1Dhists.push_back("h_deltaEtaMEm11");
    All1Dhists.push_back("h_deltaEtaMEm11a");
    All1Dhists.push_back("h_deltaEtaMEm11b");
    All1Dhists.push_back("h_deltaEtaMEm12");
    All1Dhists.push_back("h_deltaEtaMEm13");
    All1Dhists.push_back("h_deltaEtaMEm21");
    All1Dhists.push_back("h_deltaEtaMEm22");
    All1Dhists.push_back("h_deltaEtaMEm31");
    All1Dhists.push_back("h_deltaEtaMEm32");
    All1Dhists.push_back("h_deltaEtaMEm41");
    All1Dhists.push_back("h_deltaEtaMEm42");

    All1Dhists.push_back("h_m_p0");
    All1Dhists.push_back("h_m_p1");
    All1Dhists.push_back("h_m_p2");
    All1Dhists.push_back("h_m_p3");
    All1Dhists.push_back("h_m_distance");

    All1Dhists.push_back("h_p_p0");
    All1Dhists.push_back("h_p_p1");
    All1Dhists.push_back("h_p_p2");
    All1Dhists.push_back("h_p_p3");
    All1Dhists.push_back("h_p_distance");



    All1Dhists.push_back("h_deltaPhiMEpSec0");
    All1Dhists.push_back("h_deltaPhiMEpSec1");
    All1Dhists.push_back("h_deltaPhiMEpSec2");
    All1Dhists.push_back("h_deltaPhiMEpSec3");
    All1Dhists.push_back("h_deltaPhiMEpSec4");
    All1Dhists.push_back("h_deltaPhiMEpSec5");

    All1Dhists.push_back("h_deltaPhiMEmSec0");
    All1Dhists.push_back("h_deltaPhiMEmSec1");
    All1Dhists.push_back("h_deltaPhiMEmSec2");
    All1Dhists.push_back("h_deltaPhiMEmSec3");
    All1Dhists.push_back("h_deltaPhiMEmSec4");
    All1Dhists.push_back("h_deltaPhiMEmSec5");

    All1Dhists.push_back("h_deltaEtaMEpSec0");
    All1Dhists.push_back("h_deltaEtaMEpSec1");
    All1Dhists.push_back("h_deltaEtaMEpSec2");
    All1Dhists.push_back("h_deltaEtaMEpSec3");
    All1Dhists.push_back("h_deltaEtaMEpSec4");
    All1Dhists.push_back("h_deltaEtaMEpSec5");

    All1Dhists.push_back("h_deltaEtaMEmSec0");
    All1Dhists.push_back("h_deltaEtaMEmSec1");
    All1Dhists.push_back("h_deltaEtaMEmSec2");
    All1Dhists.push_back("h_deltaEtaMEmSec3");
    All1Dhists.push_back("h_deltaEtaMEmSec4");
    All1Dhists.push_back("h_deltaEtaMEmSec5");

    All1Dhists.push_back("h_p_dr0");
    All1Dhists.push_back("h_m_dr0");
    All1Dhists.push_back("h_m_chi2");
    All1Dhists.push_back("h_p_chi2");




    vector<TString> All2Dhists;
    All2Dhists.push_back("h_deltaEtavsPhi_MEp");
    All2Dhists.push_back("h_deltaEtavsPhi_MEm");
    All2Dhists.push_back("h_deltaPhivsEta_MEp");
    All2Dhists.push_back("h_deltaPhivsEta_MEm");


    All2Dhists.push_back("h_dPhivs10degPhi_MEp");
    All2Dhists.push_back("h_dPhivs10degPhi_MEm");
    All2Dhists.push_back("h_dEtavs10degPhi_MEp");
    All2Dhists.push_back("h_dEtavs10degPhi_MEm");


    All2Dhists.push_back("h_dEtavs10degPhi_MEp11");
    All2Dhists.push_back("h_dEtavs10degPhi_MEp11a");
    All2Dhists.push_back("h_dEtavs10degPhi_MEp11b");
    All2Dhists.push_back("h_dEtavs10degPhi_MEp12");
    All2Dhists.push_back("h_dEtavs10degPhi_MEp13");
    All2Dhists.push_back("h_dEtavs10degPhi_MEp21");
    All2Dhists.push_back("h_dEtavs10degPhi_MEp22");
    All2Dhists.push_back("h_dEtavs10degPhi_MEp31");
    All2Dhists.push_back("h_dEtavs10degPhi_MEp32");
    All2Dhists.push_back("h_dEtavs10degPhi_MEp41");
    All2Dhists.push_back("h_dEtavs10degPhi_MEp42");
    All2Dhists.push_back("h_dEtavs10degPhi_MEm11");
    All2Dhists.push_back("h_dEtavs10degPhi_MEm11a");
    All2Dhists.push_back("h_dEtavs10degPhi_MEm11b");
    All2Dhists.push_back("h_dEtavs10degPhi_MEm12");
    All2Dhists.push_back("h_dEtavs10degPhi_MEm13");
    All2Dhists.push_back("h_dEtavs10degPhi_MEm21");
    All2Dhists.push_back("h_dEtavs10degPhi_MEm22");
    All2Dhists.push_back("h_dEtavs10degPhi_MEm31");
    All2Dhists.push_back("h_dEtavs10degPhi_MEm32");
    All2Dhists.push_back("h_dEtavs10degPhi_MEm41");
    All2Dhists.push_back("h_dEtavs10degPhi_MEm42");


    All2Dhists.push_back("h_dPhivs10degPhi_MEp11");
    All2Dhists.push_back("h_dPhivs10degPhi_MEp11a");
    All2Dhists.push_back("h_dPhivs10degPhi_MEp11b");
    All2Dhists.push_back("h_dPhivs10degPhi_MEp12");
    All2Dhists.push_back("h_dPhivs10degPhi_MEp13");
    All2Dhists.push_back("h_dPhivs10degPhi_MEp21");
    All2Dhists.push_back("h_dPhivs10degPhi_MEp22");
    All2Dhists.push_back("h_dPhivs10degPhi_MEp31");
    All2Dhists.push_back("h_dPhivs10degPhi_MEp32");
    All2Dhists.push_back("h_dPhivs10degPhi_MEp41");
    All2Dhists.push_back("h_dPhivs10degPhi_MEp42");
    All2Dhists.push_back("h_dPhivs10degPhi_MEm11");
    All2Dhists.push_back("h_dPhivs10degPhi_MEm11a");
    All2Dhists.push_back("h_dPhivs10degPhi_MEm11b");
    All2Dhists.push_back("h_dPhivs10degPhi_MEm12");
    All2Dhists.push_back("h_dPhivs10degPhi_MEm13");
    All2Dhists.push_back("h_dPhivs10degPhi_MEm21");
    All2Dhists.push_back("h_dPhivs10degPhi_MEm22");
    All2Dhists.push_back("h_dPhivs10degPhi_MEm31");
    All2Dhists.push_back("h_dPhivs10degPhi_MEm32");
    All2Dhists.push_back("h_dPhivs10degPhi_MEm41");
    All2Dhists.push_back("h_dPhivs10degPhi_MEm42");




    All2Dhists.push_back("h_dEtavsFitEta_MEp11");
    All2Dhists.push_back("h_dEtavsFitEta_MEp11a");
    All2Dhists.push_back("h_dEtavsFitEta_MEp11b");
    All2Dhists.push_back("h_dEtavsFitEta_MEp12");
    All2Dhists.push_back("h_dEtavsFitEta_MEp13");
    All2Dhists.push_back("h_dEtavsFitEta_MEp21");
    All2Dhists.push_back("h_dEtavsFitEta_MEp22");
    All2Dhists.push_back("h_dEtavsFitEta_MEp31");
    All2Dhists.push_back("h_dEtavsFitEta_MEp32");
    All2Dhists.push_back("h_dEtavsFitEta_MEp41");
    All2Dhists.push_back("h_dEtavsFitEta_MEp42");
    All2Dhists.push_back("h_dEtavsFitEta_MEm11");
    All2Dhists.push_back("h_dEtavsFitEta_MEm11a");
    All2Dhists.push_back("h_dEtavsFitEta_MEm11b");
    All2Dhists.push_back("h_dEtavsFitEta_MEm12");
    All2Dhists.push_back("h_dEtavsFitEta_MEm13");
    All2Dhists.push_back("h_dEtavsFitEta_MEm21");
    All2Dhists.push_back("h_dEtavsFitEta_MEm22");
    All2Dhists.push_back("h_dEtavsFitEta_MEm31");
    All2Dhists.push_back("h_dEtavsFitEta_MEm32");
    All2Dhists.push_back("h_dEtavsFitEta_MEm41");
    All2Dhists.push_back("h_dEtavsFitEta_MEm42");

    All2Dhists.push_back("h_dEtavsLUTEta_MEp11");
    All2Dhists.push_back("h_dEtavsLUTEta_MEp11a");
    All2Dhists.push_back("h_dEtavsLUTEta_MEp11b");
    All2Dhists.push_back("h_dEtavsLUTEta_MEp12");
    All2Dhists.push_back("h_dEtavsLUTEta_MEp13");
    All2Dhists.push_back("h_dEtavsLUTEta_MEp21");
    All2Dhists.push_back("h_dEtavsLUTEta_MEp22");
    All2Dhists.push_back("h_dEtavsLUTEta_MEp31");
    All2Dhists.push_back("h_dEtavsLUTEta_MEp32");
    All2Dhists.push_back("h_dEtavsLUTEta_MEp41");
    All2Dhists.push_back("h_dEtavsLUTEta_MEp42");
    All2Dhists.push_back("h_dEtavsLUTEta_MEm11");
    All2Dhists.push_back("h_dEtavsLUTEta_MEm11a");
    All2Dhists.push_back("h_dEtavsLUTEta_MEm11b");
    All2Dhists.push_back("h_dEtavsLUTEta_MEm12");
    All2Dhists.push_back("h_dEtavsLUTEta_MEm13");
    All2Dhists.push_back("h_dEtavsLUTEta_MEm21");
    All2Dhists.push_back("h_dEtavsLUTEta_MEm22");
    All2Dhists.push_back("h_dEtavsLUTEta_MEm31");
    All2Dhists.push_back("h_dEtavsLUTEta_MEm32");
    All2Dhists.push_back("h_dEtavsLUTEta_MEm41");
    All2Dhists.push_back("h_dEtavsLUTEta_MEm42");








/*
    All1Dhists.clear();
    All1Dhists.push_back("h_deltaEtaMEmSec5");

    All2Dhists.clear();
    All2Dhists.push_back("h_dPhivs10degPhi_MEp");
*/

    //gStyle->SetOptStat(0);




    gStyle->SetOptStat("emruoi");

    for(size_t ihist=0; ihist<All1Dhists.size(); ihist++) {

        TH1F* hist_2 = (TH1F*) infile->Get(All1Dhists[ihist]);
        if(hist_2==NULL) {
            cout << "hist is NULL!" << endl;
            continue;
        }

        TCanvas *c = new TCanvas("c", "c", 700, 550);
        //TCanvas *c = new TCanvas("c", "c", 900, 550);
        TPad* t1 = new TPad("t1","t1", 0.0, 0.0, 1.0, 1.00);
        t1->Draw();
        t1->cd();
        //t1->SetBottomMargin(0.3);
        t1->SetRightMargin(0.03);
        t1->SetLogy(1);
        //c->Divide(1,2);


        hist_2->Draw("");

        if(All1Dhists[ihist].Contains("Eta")|| All1Dhists[ihist].Contains("Phi")) {
            hist_2->Fit("gaus","+","",-0.1,0.1);
        }


        if(hist_2->GetFunction("gaus"))
            hist_2->GetFunction("gaus")->SetLineWidth(2);

        //hx_->GetXaxis()->SetTitle("Vertices");
        hist_2->GetYaxis()->SetTitle("Events");
        hist_2->SetFillColor(kYellow);
        if(All1Dhists[ihist].Contains("MEp") || All1Dhists[ihist].Contains("_p_")) hist_2->SetFillColor(kGreen-7);
        if(All1Dhists[ihist].Contains("MEm") || All1Dhists[ihist].Contains("_m_")) hist_2->SetFillColor(kOrange-3);


        c->Modified();
        c->Update();
        TPaveStats *stats =  (TPaveStats*) hist_2->GetListOfFunctions()->FindObject("stats");
        stats->SetFillStyle(0);
        stats->SetName("");
        stats->SetX1NDC(.75);
        stats->SetY1NDC(.60);
        stats->SetX2NDC(.95);
        stats->SetY2NDC(.93);
        stats->SetTextColor(2);

        c->Update();

        c->SaveAs(All1Dhists[ihist]+".png");
        c->SaveAs(All1Dhists[ihist]+".pdf");

        delete c;
        delete hist_2;
    }




    for(size_t ihist=0; ihist<All2Dhists.size(); ihist++) {

        TH2F* hist_temp = (TH2F*) infile->Get(All2Dhists[ihist]);
        if(hist_temp==NULL) {
            cout << "hist is NULL!" << endl;
            continue;
        }

        TCanvas *c = new TCanvas("c", "c", 700, 550);
        TPad* t1 = new TPad("t1","t1", 0.0, 0.0, 1.0, 1.0);
        t1->Draw();
        t1->cd();
        //t1->SetBottomMargin(0.3);
        t1->SetRightMargin(0.12);
        //c->Divide(1,2);
        hist_temp->Draw("col Z");
        //hx_->GetXaxis()->SetTitle("Vertices");
        //hx_->GetYaxis()->SetTitle("<E_{Y}^{miss}>");


        c->Modified();
        c->Update();
        TPaveStats *stats =  (TPaveStats*) hist_temp->GetListOfFunctions()->FindObject("stats");
	stats->SetFillStyle(0);
        stats->SetName("");
        //stats->SetFillColor();
        stats->SetX1NDC(.70);
        stats->SetY1NDC(.67);
        stats->SetX2NDC(.85);
        stats->SetY2NDC(.93);
        stats->SetTextColor(2);




        c->SaveAs(All2Dhists[ihist]+".png");
        c->SaveAs(All2Dhists[ihist]+".pdf");

        delete c;
        delete hist_temp;

    }


    infile->Close();

}
