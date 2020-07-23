#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TString.h"
#include "TStyle.h"
#include "TLegend.h"

#include <iostream>
#include <map>
#include <vector>

using namespace std;

const int n_Modules = 4;
const int n_Layers = 2;
const int n_Groups = 3;
TString LayerName[n_Layers] = {"v","h"};

int DrawDerivativePlot(int nEvts = 10)
{
    gStyle->SetOptStat(0);
    
    TFile* f = new TFile("Cluster_output.root");

    int File_idx = 0;
    TString hisname,hisname1;
    TCanvas* c1 = new TCanvas("c1","c1",1600,1200);
    TLegend* leg = new TLegend(0.6,0.91,0.89,0.95);
    leg->SetLineColor(0);
    leg->SetNColumns(2);
    for (int i = 0; i < n_Modules ; i++)
    {
        for ( int j = 0; j < n_Layers; j++ )
        {
            for ( int k = 0; k < n_Groups; k++ )
            {
                hisname = Form("hDIG3L%d%sG%d",i+1,LayerName[j].Data(),k);
                hisname1 = Form("hDIG3L%d%sG%d_hDerivative",i+1,LayerName[j].Data(),k);
                TH1D* his = (TH1D*)f->Get(hisname);
                TH1D* his1 = (TH1D*)f->Get(hisname1);
                leg->AddEntry(his,"ADC","pl");
                leg->AddEntry(his1,"1th derivative","pl");
                his->SetLineColor(1);
                his1->SetLineColor(2);
                his->GetYaxis()->SetRangeUser(-500,500);
                his->Draw("hist");
                his1->Draw("histsame");
                leg->Draw("same");
                TString savename = Form("./plots/Derivative/hDIG3L%d%sG%d.png",i+1,LayerName[j].Data(),k);
                c1->SaveAs(savename.Data());
                c1->Clear();
                leg->Clear();
            }
        }
    }

    return 1;
}