#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TString.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraph.h"

#include <iostream>
#include <map>
#include <vector>

using namespace std;

#include "/Users/wztbxl/Documents/macros/function.C"

const int n_Modules = 4;
const int n_Layers = 2;
const int n_Groups = 3;
TString LayerName[n_Layers] = {"v","h"};

int CheckReconstruct(int nFiles = 100)
{
    gStyle->SetOptStat(0);


    TH2D* nRc_vs_MC = new TH2D("nRc_vs_MC","nRc_vs_MC;MC Hits;RC Clusters",200,0,200,200,0,200);
    TH1D* hRC_rate = new TH1D("hRC_rate","hRC_rate;nRC/nMC;Ratio",1000,0,10);
    TH1D* hRC_Correct_rate = new TH1D("hRC_Correct_rate","hRC_Correct_rate;hRC_Correct/hMC;Ratio",200,0,2);

    TCanvas* c1 = new TCanvas("c1","c1",1600,1200);
    c1->SetLogz(1);
    TLegend* leg = new TLegend(0.5,0.91,0.9,0.99);
    leg->SetLineColor(0);
    leg->SetNColumns(2);

    for ( int i = 0; i < nFiles; i++)
    {
        double nCorrect_Hits = 0;

        TFile* fRC = new TFile(Form("out/Evts_10/Cluster_output_Evts10_%i_v1.root",i));
        TFile* fMC = new TFile(Form("../stgc-cluster-sim/output/output_%i.root",i));

        TGraph* RCHits = (TGraph*)fRC->Get("Hits_map_RC");
        setGraph(RCHits,8,1,2,2,1);
        TGraph* MCHits = (TGraph*)fMC->Get("Hits_map");
        setGraph(MCHits,29,2,1,1,1);
        MCHits->SetTitle("");
        TH2D* hDIG2 = (TH2D*)fMC->Get("hGEN0");

        leg->AddEntry(RCHits,"RC Cluster","p");
        leg->AddEntry(MCHits,"MC Cluster","p");

        hDIG2->Draw("colz");
        MCHits->Draw("psame");
        RCHits->Draw("psame");
        leg->Draw("same");

        c1->SaveAs(Form("./plots/CheckRC/%i.png",i));
        c1->SaveAs(Form("./plots/CheckRC/%i.pdf",i));
        c1->Clear();
        leg->Clear();

        nRc_vs_MC->Fill(MCHits->GetN(),RCHits->GetN());
        hRC_rate->Fill((double)RCHits->GetN()/MCHits->GetN());
        for (int j = 0; j < RCHits->GetN(); j++)
        {
            double xRC = 0;
            double yRC = 0;
            double xMC = 0;
            double yMC = 0;
            int is_Correct = 0;

            RCHits->GetPoint(j,xRC,yRC);
            for (int k = 0; k < MCHits->GetN(); k++)
            {
                MCHits->GetPoint(k,xMC,yMC);
                double R = 0;
                R = pow((xRC-xMC),2)+pow((yRC-yMC),2);
                R = sqrt(R);
                if (R < 5) nCorrect_Hits++;
                // cout << "R = " << R << endl;
            }
        }
        cout << "n_Correct_RC = " << nCorrect_Hits << endl;
        cout << endl; 
        hRC_Correct_rate->Fill(nCorrect_Hits/MCHits->GetN());
    }

    TFile* outFile = new TFile("CheckRC.root","recreate");
    nRc_vs_MC->Write();
    hRC_rate->Write();
    hRC_Correct_rate->Write();

    outFile->Write();
    outFile->Close();

    return 1;
}