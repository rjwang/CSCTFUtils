//
//
#include "/afs/cern.ch/user/r/rewang/MainResults/AnaCode/plotPost.h"





void plotter()
{

    TFile *infile = TFile::Open("output_test.root", "READ");

    vector<TString> All1Dhists;

    All1Dhists.push_back("d_phi");
    All1Dhists.push_back("d_eta");

    All1Dhists.push_back("h_deltaPhiMEp11");
    All1Dhists.push_back("h_deltaPhiMEp12");
    All1Dhists.push_back("h_deltaPhiMEp13");
    All1Dhists.push_back("h_deltaPhiMEp21");
    All1Dhists.push_back("h_deltaPhiMEp22");
    All1Dhists.push_back("h_deltaPhiMEp31");
    All1Dhists.push_back("h_deltaPhiMEp32");
    All1Dhists.push_back("h_deltaPhiMEp41");
    All1Dhists.push_back("h_deltaPhiMEp42");


    All1Dhists.push_back("h_deltaPhiMEm11");
    All1Dhists.push_back("h_deltaPhiMEm12");
    All1Dhists.push_back("h_deltaPhiMEm13");
    All1Dhists.push_back("h_deltaPhiMEm21");
    All1Dhists.push_back("h_deltaPhiMEm22");
    All1Dhists.push_back("h_deltaPhiMEm31");
    All1Dhists.push_back("h_deltaPhiMEm32");
    All1Dhists.push_back("h_deltaPhiMEm41");
    All1Dhists.push_back("h_deltaPhiMEm42");

    All1Dhists.push_back("h_deltaEtaMEp11");
    All1Dhists.push_back("h_deltaEtaMEp12");
    All1Dhists.push_back("h_deltaEtaMEp13");
    All1Dhists.push_back("h_deltaEtaMEp21");
    All1Dhists.push_back("h_deltaEtaMEp22");
    All1Dhists.push_back("h_deltaEtaMEp31");
    All1Dhists.push_back("h_deltaEtaMEp32");
    All1Dhists.push_back("h_deltaEtaMEp41");
    All1Dhists.push_back("h_deltaEtaMEp42");


    All1Dhists.push_back("h_deltaEtaMEm11");
    All1Dhists.push_back("h_deltaEtaMEm12");
    All1Dhists.push_back("h_deltaEtaMEm13");
    All1Dhists.push_back("h_deltaEtaMEm21");
    All1Dhists.push_back("h_deltaEtaMEm22");
    All1Dhists.push_back("h_deltaEtaMEm31");
    All1Dhists.push_back("h_deltaEtaMEm32");
    All1Dhists.push_back("h_deltaEtaMEm41");
    All1Dhists.push_back("h_deltaEtaMEm42");

    All1Dhists.push_back("h_p0");
    All1Dhists.push_back("h_p1");
    All1Dhists.push_back("h_p2");
    All1Dhists.push_back("h_p3");
    All1Dhists.push_back("h_dR");














    //gStyle->SetOptStat(0);




    gStyle->SetOptStat("nemruoi");

    for(size_t ihist=0; ihist<All1Dhists.size(); ihist++) {

                TH1F* hist_2 = (TH1F*) infile->Get(All1Dhists[ihist]);
                if(hist_2==NULL) {
                    cout << "hist is NULL!" << endl;
                    continue;
                }

                TCanvas *c = new TCanvas("c", "c", 700, 550);
                TPad* t1 = new TPad("t1","t1", 0.0, 0.0, 1.0, 1.00);
                t1->Draw();
                t1->cd();
                //t1->SetBottomMargin(0.3);
                t1->SetRightMargin(0.05);
		t1->SetLogy(1);
                //c->Divide(1,2);
                hist_2->Draw("hist");
                //hx_->GetXaxis()->SetTitle("Vertices");
                hist_2->GetYaxis()->SetTitle("Events");
		hist_2->SetFillColor(kYellow);

		c->SaveAs(All1Dhists[ihist]+".png");
                c->SaveAs(All1Dhists[ihist]+".pdf");

                delete c;
		delete hist_2;
    }





    infile->Close();

}
