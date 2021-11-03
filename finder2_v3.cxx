#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TString.h"
#include "TGraph.h"
#include "TLine.h"
#include "TVector2.h"

#include <iostream>
#include <map>
#include <vector>
#include <algorithm>

#include "../Public_input/Parameters.h"
#include "sTGC_DataDst.cxx"

using namespace std;

sTGC_DataDst* sTGC_Evt = new sTGC_DataDst();

const int Kernel_size = 1;
const int n_Modules = 4;
const int n_Layers = 2;
const int n_Groups = 3;
const int n_Diag_Groups = 2;
const int nHits = 100;
map <string, TH1D*> ADC_Dis_1D; // 1D ADC distribution
map <string, TH1D*> ADC_Dis_1D_Diag; // 1D ADC distribution
map <string, TH1D*> ADC_Der_1D; // Derivative of the 1D ADC distribution
map <string, TH1D*> ADC_Der_2nd_1D; // 2nd order Derivative of the 1D ADC distribution
map<string, TLine *> Edges; // 
vector <int> minimum;// bin number of maxmium and minimum in "signal region" of each layer 
vector <int> maximum;
vector <int> start_point;// bin number the start point and end point of every "signal region" of each layer
vector <int> end_point;
map <string, vector<int> > minimum_DiffLayer; //record the information of every layer
map <string, vector<int> > maximum_DiffLayer;
map <string, vector<int> > start_point_DiffLayer;
map <string, vector<int> > end_point_DiffLayer;
map <int, vector<double> > Diag_Position;
multimap < int ,int > ADC_per_Channels;
vector <int> ADC_Channels;
vector <double> cluster_Posi;
// Array to save the position
double Xhit_position[n_Modules][n_Groups][nHits] = {0.0};
double Yhit_position[n_Modules][n_Groups][nHits] = {0.0};
double Diag_hit_position[n_Modules][n_Groups][nHits] = {0.0};

TString LayerName[2] = {"v","h"};

const float pent_base = 537.0;    // in mm
const float pent_nib = 179.0;     // in mm
const float pent_shift = 101.6;   // in mm
const int n_clusters_to_gen = 1; // Number of clusters in this "event"
const float noise_prob = 0.2;     // 0 -1, 1 is max noise
// const float noise_level = 1; // in ADC units (pre integration)
const float noise_level = 1.e-10;  // in ADC units (pre integration)
const float strip_pitch = 3.2;     // digitize strip pitch
const size_t sat_above = 1024;     // ADC
const float cluster_max_adc = 120; // ADC (pre integration)
int test_edge = 1;
int Edge_Evt = 1;

TH1D * hDerivative( TH1D * hin, int kernel_size = 1 )
{
	// Step 2
	TH1D * hAv = (TH1D*)hin->Clone( "hAv" );
	hAv->Reset();
	for ( int i = 1; i < hin->GetNbinsX()+1; i++ ){
		float bin_content = 0;
		for ( int j = i; j < i + kernel_size; j++ ){
            if (j > hin->GetNbinsX()) break;
			// check to make sure j is in bounds (not past last bin)
			bin_content += hin->GetBinContent( j );
		}
		bin_content = bin_content / (float)kernel_size;
		
		hAv->SetBinContent( i, bin_content );
	
	}
	
	
	TH1D * hDerivative = (TH1D*)hin->Clone( "hDerivative" );
    TString name = hin->GetName();
    name = name+"_hDerivative";
    hDerivative->SetName(name);
	hDerivative->Reset();
	// now Loop over hAv and compute derivative between neighboring bins
	for ( int i = 1; i < hAv->GetNbinsX(); i++ ){
		float bc1 = hAv->GetBinContent( i );
		float bc2 = hAv->GetBinContent( i+1 );
        float x1 = hAv->GetBinCenter(i);
        float x2 = hAv->GetBinCenter(i+1);
		// compute derivative // NOTE: need to get bin center for x1 and x2
		float d = (bc2 - bc1) / (x2 - x1);
		hDerivative->SetBinContent( i, d );
	}
	
	return hDerivative;

}
TH1D* BKG_Substract(TH1D* hist, int BKG_level) // using a constant background to do the prelimary study
{
    TH1D* hAfter_Sub = (TH1D*)hist->Clone();
    for (int i = 1; i < hist->GetNbinsX()+1; i++)
    {
        double bc = hist->GetBinContent(i);
        double bc_sub = bc-BKG_level;
        if ( bc-BKG_level < 0 )  bc_sub = 0;
        hAfter_Sub->SetBinContent(i,bc_sub);
    }

    return hAfter_Sub;

}

bool find_key_point( string name, vector <int> &vec_start, vector <int> &vec_end, vector <int> &vec_min, vector <int> &vec_max)
{
    int is_signal = 0;
    TH1D* hist_start = new TH1D("hist_start","hist_start",200,0,200);
    hist_start->SetFillColor(1);
    TH1D* hist_end = new TH1D("hist_end","hist_end",200,0,200);
    hist_end->SetFillColor(2);
    TH1D* hist_min = new TH1D("hist_min","hist_min",200,0,200);
    hist_min->SetFillColor(4);
    TH1D* hist_max = new TH1D("hist_max","hist_max",200,0,200);
    hist_max->SetFillColor(6);
    for (int i = 1; i < ADC_Der_1D[name]->GetNbinsX()-1; i++)
    {
        double bc1 = ADC_Der_1D[name]->GetBinContent(i);
        double bc2 = ADC_Der_1D[name]->GetBinContent(i+1);
        if (bc1 != 0 && i == 1 && bc2 !=0) 
        {
            is_signal = 1;
            vec_start.push_back(i);
            hist_start->SetBinContent(i,200);
        }
        else
        {
            if (bc1 == 0 && bc2 > 0) 
            {
                double Dis_bc1 = ADC_Dis_1D[name]->GetBinContent(i+2);
                double Dis_bc2 = ADC_Dis_1D[name]->GetBinContent(i+3);
                double Dis_bc3 = ADC_Dis_1D[name]->GetBinContent(i+4);
                // cout << " bin number of Dis_bc1 = " << i+2 << " bin number of Dis_bc2 = " << i+3 << " bin number of Dis_bc3 = " << i+4 << endl;
                // cout << " Dis_bc1 = " << Dis_bc1 << " Dis_bc2 = " << Dis_bc2 << " Dis_bc3 = " << Dis_bc3 << endl;
                if ( Dis_bc1 > 0 && Dis_bc2 > 0 )
                {
                    is_signal = 1;
                    vec_start.push_back(i+1);
                    hist_start->SetBinContent(i+1,200);
                }
            } // start piont
            if (bc1 < 0 && bc2 == 0) 
            {
                if (is_signal == 1)
                {
                    vec_end.push_back(i);
                    hist_end->SetBinContent(i,200);
                    is_signal = 0;
                }
            } // end point
            if (bc1 > 0 && bc2 < 0) 
            {
                if (is_signal == 1)
                {
                    vec_max.push_back(i);hist_max->SetBinContent(i,1000);
                }
            } // maximum
            if (bc1 < 0 && bc2 > 0) 
            {
                if (is_signal == 1)
                {
                    vec_min.push_back(i);hist_min->SetBinContent(i,1000);
                }
            } // minimum
            if (bc1 == 0 && bc2 == 0) continue;
        }
    }
    TCanvas* c1 = new TCanvas("c1","c1",1600,1200);
    ADC_Dis_1D[name]->DrawClone("hist");
    hist_start->Print();
    hist_start->DrawClone("same");
    hist_end->DrawClone("same");
    hist_min->DrawClone("same");
    hist_max->DrawClone("same");
    TString savename = "./plots/"+name+".png";
    c1->SaveAs(savename.Data());

    delete hist_start;
    delete hist_end;
    delete hist_min;
    delete hist_max;
    hist_start = NULL;
    hist_end = NULL;
    hist_min = NULL;
    hist_max = NULL;
    delete c1 ;

    return kTRUE;

}

int Get_Cluster_Position( TH1D* hGroup, vector < double > &cluster_posi, vector <int> &vec_ADC_per_Channels , int cluster_width = 4 , int Diagnoal_tag = 0) // cluster_width is the strip number to search, defult is 4, Diagnoal_tag = 1 is the diagnoal strip
{
    vec_ADC_per_Channels.clear();
    for ( int i = 1; i < hGroup->GetNbinsX(); i++)
    {
        vec_ADC_per_Channels.push_back(hGroup->GetBinContent(i));
    }

    sort(vec_ADC_per_Channels.rbegin(),vec_ADC_per_Channels.rend());
    TString name = hGroup->GetName();
    string name1 = name.Data();
    string Module ;
    string Direction;
    string Group;
    int module;
    int direction;
    int group;
    int Channel_number;
    if ( name1.size() < 11 )
    {
        Module = name1.substr(6,1);
        Direction = name1.substr(7,1);
        Group = name1.substr(9,1);
        module = stoi(Module);
        group = stoi(Group);
    }
    if ( name1.size() > 11)
    {
        Module = name1.substr(6,1);
        Direction = name1.substr(12,1);
        module = stoi(Module);
        direction = stoi(Direction);
    }

    while ( hGroup->GetMaximum() > 0 )
    {
        int i = 0; //right
        int j = 0; //left
        int bin = hGroup->GetMaximumBin();
        double sumADC = 0;
        double strip_ADC = 0;
        double bc = hGroup->GetBinContent(bin);
        double bc_right = bc;
        double bc_left = bc;
        int left_tag = 1;
        int right_tag = 1;
        double Max_BCID = 0;
        double BCID = 0;
        // while( right_tag == 1 || left_tag == 1 ) // find one cluster, change this from && to ||, in logical, both right and left should be 0
        while( right_tag == 1 && left_tag == 1 ) // find one cluster 
        {
            if (i == 0) {sumADC+=bc; strip_ADC = strip_ADC+bc*((bin*3.2)-1.6); hGroup->SetBinContent(bin,0); i++; j++; }
            if ( hGroup->GetBinContent(bin+i) < bc_right && hGroup->GetBinContent(bin+i) > 0 ) 
            {
                sumADC+=bc_right; 
                strip_ADC = strip_ADC+bc_right*(((bin+i)*3.2)-1.6); 
                bc_right = hGroup->GetBinContent(bin+i);
                hGroup->SetBinContent(bin+i,0);
                i++;
            }
            else
            {
                right_tag = 0;
            }
            
            if ( hGroup->GetBinContent(bin-j) < bc_left && hGroup->GetBinContent(bin-j) > 0 ) 
            {
                sumADC+=bc_left; 
                strip_ADC = strip_ADC+bc_left*(((bin-j)*3.2)-1.6); 
                bc_left = hGroup->GetBinContent(bin-j);
                hGroup->SetBinContent(bin-j,0);
                j++;
            }
            else
            {
                left_tag = 0;
            }
        }
        
        int width = i+j-1;
        if (width >= cluster_width)
        {
            // cout << "i = " << i <<" j = " << j << endl;
            double position = strip_ADC/sumADC;
            if ( Module == "3" && direction == "v" ) position = position + 0.8;
            if ( Module == "4" && direction == "v" ) position = position + 0.8;
            cluster_posi.push_back(position);
        }
        else 
        {
            if ( (bin == 166 || bin == 167 || bin == 0 || bin == 1) && Diagnoal_tag == 0 )
            {
                // cout << "i = " << i <<" j = " << j << endl;
                double position = bin*3.2+1.6;
                if ( Module == "3" && direction == "v" ) position = position - 0.8;
                if ( Module == "4" && direction == "v" ) position = position + 0.8;
                if ( bin == 166 || bin == 167 )  cluster_posi.push_back(bin*3.2+1.6);
                if ( bin == 0 || bin == 1 )  cluster_posi.push_back(bin*3.2+1.6);
            }
            if ( (bin == 151 || bin == 150 || bin == 0 || bin == 1) && Diagnoal_tag == 1 )
            {
                double position = bin*3.2+1.6;
                if ( bin == 166 || bin == 167 )  cluster_posi.push_back(bin*3.2+1.6);
                if ( bin == 0 || bin == 1 )  cluster_posi.push_back(bin*3.2+1.6);
            }
        }
    }

    return 1;
}

int stripGroup( int sx, int sy ){
    if ( (sx < 55 && sy < 150) || sy >= 150 ) return 0;
    if ( (sx < 110 && sy < 95) || sy >= 95 ) return 1;
    return 2;
}

int strip_to_position_group(int hGroup, int vGroup, int sx, int sy) // hGroup's read is y coordinate, sy，vGroup's read out is x coordiante, if return is negative, that means it is not a real hit
{
    int position = -999;
    if (hGroup == 0 && vGroup == 0 && sx < 55 & sy < 55) position = 1;
    if (hGroup == 1 && vGroup == 0 && sx >=55 && sx < 110 && sy < 55 ) position = 2;
    if ( (hGroup == 0 && vGroup == 1 && sx < 55 && sy >= 55 && sy < 110) ) position = 4;
    if ( (hGroup == 2 && vGroup == 0 && sx >= 110 && sx < 150 && sy < 55) || (hGroup == 2 && vGroup == 0 && sx >= 150 && sx < 167 && sy < 60) ) position = 3;
    if ( (hGroup == 1 && vGroup == 1 && sx >= 55 && sx < 95 && sy < 110 && sy >= 55) || (hGroup == 1 && vGroup == 1 && sx >= 95 && sx < 115 && sy < 115 && sy >= 55) ) position = 5;
    if ( (hGroup == 2 && vGroup == 1 && sx >= 110 && sx < 167 && sy < 95 && sy >= 55) ) position = 6;
    if ( (hGroup == 0 && vGroup == 2 && sx < 55 && sy < 167 && sy >= 110) ) position = 7;
    if ( (hGroup == 1 && vGroup == 2 && sx >=55 && sx < 95 && sy >= 110 && sy < 150) ) position = 8;

    return position;
}

// int xy_to_position_group(int hGroup, int vGroup, double x, double y) // hGroup's read is y coordinate, sy，vGroup's read out is x coordiante, if return is negative, that means it is not a real hit
// {
//     int position = -999;
//     if (hGroup == 0 && vGroup == 0 && x < 55*3.2 & y < 55*3.2) position = 1;
//     if (hGroup == 1 && vGroup == 0 && x >=55*3.2 && x < 110*3.2 && y < 55*3.2 ) position = 2;
//     if ( (hGroup == 0 && vGroup == 1 && x < 55*3.2 && y >= 55*3.2 && y < 110*3.2) ) position = 4;
//     if ( (hGroup == 2 && vGroup == 0 && x >= 110*3.2 && x < 150*3.2 && y < 55*3.2) || (hGroup == 2 && vGroup == 0 && x >= 150*3.2 && x < 167*3.2 && y < 60*3.2) ) position = 3;
//     if ( (hGroup == 1 && vGroup == 1 && x >= 55*3.2 && x < 95*3.2 && y < 110*3.2 && y >= 55*3.2) || (hGroup == 1 && vGroup == 1 && x >= 95*3.2 && x < 115*3.2 && y < 115*3.2 && y >= 55*3.2) ) position = 5;
//     if ( (hGroup == 2 && vGroup == 1 && x >= 110*3.2 && x < 167*3.2 && y < 95*3.2 && y >= 55*3.2) ) position = 6;
//     if ( (hGroup == 0 && vGroup == 2 && x < 55*3.2 && y < 150*3.2 && y >= 110*3.2) || ( hGroup == 0 && vGroup == 2 && x < 60*3.2 && y < 167*3.2 && y >= 150*3.2) ) position = 7;
//     if ( (hGroup == 1 && vGroup == 2 && x >=55*3.2 && x < 95*3.2 && y >= 110*3.2 && y < 150*3.2) ) position = 8;

//     return position;
// }

int xy_to_position_group(int hGroup, int vGroup, double x, double y) // hGroup's read is y coordinate, sy，vGroup's read out is x coordiante, if return is negative, that means it is not a real hit
{
    int position = -999;
    if (hGroup == 0 && vGroup == 0 && x < 55*3.2 & y < 55*3.2) position = 1;
    if (hGroup == 1 && vGroup == 0 && x >=55*3.2 && x < 110*3.2 && y < 55*3.2 ) position = 2;
    if ( (hGroup == 0 && vGroup == 1 && x < 55*3.2 && y >= 55*3.2 && y < 110*3.2) ) position = 4;
    if ( (hGroup == 2 && vGroup == 0 && x >= 110*3.2 && x < 150*3.2 && y < 55*3.2) || (hGroup == 2 && vGroup == 0 && x >= 150*3.2 && x < 167*3.2 && y < 60*3.2) ) position = 3;
    if ( (hGroup == 1 && vGroup == 1 && x >= 55*3.2 && x < 95*3.2 && y < 110*3.2 && y >= 55*3.2) || (hGroup == 1 && vGroup == 1 && x >= 95*3.2 && x < 128*3.2 && y < 128*3.2 && y >= 55*3.2) ) position = 5;
    if ( (hGroup == 2 && vGroup == 1 && x >= 110*3.2 && x < 167*3.2 && y < 95*3.2 && y >= 55*3.2) ) position = 6;
    if ( (hGroup == 0 && vGroup == 2 && x < 55*3.2 && y < 150*3.2 && y >= 110*3.2) || ( hGroup == 0 && vGroup == 2 && x < 60*3.2 && y < 167*3.2 && y >= 150*3.2) ) position = 7;
    if ( (hGroup == 1 && vGroup == 2 && x >=55*3.2 && x < 95*3.2 && y >= 110*3.2 && y < 150*3.2) ) position = 8;

    return position;
}

double cluster_position(TString name, int first_bin, int last_bin)
{
    double sumADC = 0;
    double strip_ADC = 0;
    for (int i = first_bin; i <= last_bin; i++)
    {
        double bc = ADC_Dis_1D[name.Data()]->GetBinContent(i);
        sumADC += bc;
        strip_ADC = strip_ADC+bc*((i*3.2)-1.6);// because we used center of strip to calculate the position
    }
    return strip_ADC/sumADC;
}
 
double coordinate_convertion( int module, double x, double y, double &globalx, double &globaly) // convert local coordinate to global coordinate, For direction, 0 is local x, 1 is local y
{
    if (module != 1 && module !=2 && module != 3 && module != 4) 
    {
        cout << " module number should be 1 or 2 or 3 or 4 !!!!!" << endl;
        return 0;
    }
    if ( module == 1 ) {globalx = x; globaly = y;}
    if ( module == 2 ) {globalx = -x; globaly = y;}
    if ( module == 3 ) {globalx = -x-101.6; globaly = -y;}
    if ( module == 4 ) {globalx = x+101.6; globaly = -y;}

    return  1;
}

void initEdges() //function using to init the strip edges
{
    // define the edge of strip group
    Edges["M1Edge1"] = new TLine(base_point1, base_point1, base_point2, base_point1);
    Edges["M1Edge2"] = new TLine(base_point1, base_point1, base_point1, base_point2);
    Edges["M1Edge3"] = new TLine(base_point2, base_point1, base_point2, base_point3);
    Edges["M1Edge4"] = new TLine(base_point1, base_point2, base_point3, base_point2);
    Edges["M1Edge5"] = new TLine(base_point3, base_point2, base_point2, base_point3);
    Edges["M2Edge1"] = new TLine(-base_point1, base_point1, -base_point2, base_point1);
    Edges["M2Edge2"] = new TLine(-base_point1, base_point1, -base_point1, base_point2);
    Edges["M2Edge3"] = new TLine(-base_point2, base_point1, -base_point2, base_point3);
    Edges["M2Edge4"] = new TLine(-base_point1, base_point2, -base_point3, base_point2);
    Edges["M2Edge5"] = new TLine(-base_point3, base_point2, -base_point2, base_point3);
    Edges["M3Edge1"] = new TLine(-base_point1 - 101.6, -base_point1, -base_point2 - 101.6, -base_point1);
    Edges["M3Edge2"] = new TLine(-base_point1 - 101.6, -base_point1, -base_point1 - 101.6, -base_point2);
    Edges["M3Edge3"] = new TLine(-base_point2 - 101.6, -base_point1, -base_point2 - 101.6, -base_point3);
    Edges["M3Edge4"] = new TLine(-base_point1 - 101.6, -base_point2, -base_point3 - 101.6, -base_point2);
    Edges["M3Edge5"] = new TLine(-base_point3 - 101.6, -base_point2, -base_point2 - 101.6, -base_point3);
    Edges["M4Edge1"] = new TLine(base_point1 + 101.6, -base_point1, base_point2 + 101.6, -base_point1);
    Edges["M4Edge2"] = new TLine(base_point1 + 101.6, -base_point1, base_point1 + 101.6, -base_point2);
    Edges["M4Edge3"] = new TLine(base_point2 + 101.6, -base_point1, base_point2 + 101.6, -base_point3);
    Edges["M4Edge4"] = new TLine(base_point1 + 101.6, -base_point2, base_point3 + 101.6, -base_point2);
    Edges["M4Edge5"] = new TLine(base_point3 + 101.6, -base_point2, base_point2 + 101.6, -base_point3);

    Edges["M1hG0Edge1"] = new TLine(x2, py1, x2, y5);
    Edges["M1hG0Edge2"] = new TLine(x2, y5, x7, y5);
    Edges["M1hG1Edge1"] = new TLine(x3, py1, x3, y6);
    Edges["M1hG1Edge2"] = new TLine(x3, y6, x8, y6);
    Edges["M1vG0Edge1"] = new TLine(x1, y2, x5, y2);
    Edges["M1vG0Edge2"] = new TLine(x5, y2, x5, y7);
    Edges["M1vG1Edge1"] = new TLine(x1, y3, x6, y3);
    Edges["M1vG1Edge2"] = new TLine(x6, y3, x6, y8);
    Edges["M2hG0Edge1"] = new TLine(-x2, py1, -x2, y5);
    Edges["M2hG0Edge2"] = new TLine(-x2, y5, -x7, y5);
    Edges["M2hG1Edge1"] = new TLine(-x3, py1, -x3, y6);
    Edges["M2hG1Edge2"] = new TLine(-x3, y6, -x8, y6);
    Edges["M2vG0Edge1"] = new TLine(-x1, y2, -x5, y2);
    Edges["M2vG0Edge2"] = new TLine(-x5, y2, -x5, y7);
    Edges["M2vG1Edge1"] = new TLine(-x1, y3, -x6, y3);
    Edges["M2vG1Edge2"] = new TLine(-x6, y3, -x6, y8);
    Edges["M3hG0Edge1"] = new TLine(-x2 - 101.6, -py1, -x2 - 101.6, -y5);
    Edges["M3hG0Edge2"] = new TLine(-x2 - 101.6, -y5, -x7 - 101.6, -y5);
    Edges["M3hG1Edge1"] = new TLine(-x3 - 101.6, -py1, -x3 - 101.6, -y6);
    Edges["M3hG1Edge2"] = new TLine(-x3 - 101.6, -y6, -x8 - 101.6, -y6);
    Edges["M3vG0Edge1"] = new TLine(-x1 - 101.6, -y2, -x5 - 101.6, -y2);
    Edges["M3vG0Edge2"] = new TLine(-x5 - 101.6, -y2, -x5 - 101.6, -y8);
    Edges["M3vG1Edge1"] = new TLine(-x1 - 101.6, -y3, -x6 - 101.6, -y3);
    Edges["M3vG1Edge2"] = new TLine(-x6 - 101.6, -y3, -x6 - 101.6, -y8);
    Edges["M4hG0Edge1"] = new TLine(x2 + 101.6, -py1, x2 + 101.6, -y5);
    Edges["M4hG0Edge2"] = new TLine(x2 + 101.6, -y5, x7 + 101.6, -y5);
    Edges["M4hG1Edge1"] = new TLine(x3 + 101.6, -py1, x3 + 101.6, -y6);
    Edges["M4hG1Edge2"] = new TLine(x3 + 101.6, -y6, x8 + 101.6, -y6);
    Edges["M4vG0Edge1"] = new TLine(x1 + 101.6, -y2, x5 + 101.6, -y2);
    Edges["M4vG0Edge2"] = new TLine(x5 + 101.6, -y2, x5 + 101.6, -y7);
    Edges["M4vG1Edge1"] = new TLine(x1 + 101.6, -y3, x6 + 101.6, -y3);
    Edges["M4vG1Edge2"] = new TLine(x6 + 101.6, -y3, x6 + 101.6, -y8);
}

bool in_bounds( float x, float y ){
    if ( x>0 && x < pent_nib && y < pent_base && y > 0 )
        return true;
    if ( x>0 && x < pent_base && y < pent_nib && y > 0 )
        return true;

    float dy = pent_base - pent_nib;
    float dx = pent_nib - pent_base;
    float yedge = (dy/dx) * ( x - pent_nib ) + pent_base;
    // LOG_F( INFO, "yedge = %f @ x = %f", yedge, x );

    if ( y < yedge && y < pent_base && x < pent_base && x>0 && y>0 ) return true;

    return false;
}

int in_bounds_quad( float x, float y ){
    if ( in_bounds( x, y ) ){ // top right
        return 1;
    } else if ( in_bounds( -x, y ) ){ //top left
        return 2;
    } else if ( in_bounds( -x - pent_shift, -y ) ){ // bottom left
        return 3;
    } else if ( in_bounds( x - pent_shift, -y ) ){ // bottom right
        return 4;
    }

    return 0;
}

double Distance_Point2Line_2D( TLine* line, double x, double y )
{
  double distance;
  TVector2 A; TVector2 B; TVector2 M;
  TVector2 AB; TVector2 BM;
  A.Set( line->GetX1(), line->GetY1() );
  B.Set( line->GetX2(), line->GetY2() );
  M.Set( x, y );
  AB = A-B;
  BM = M-B;
  distance = abs( AB ^ BM ) / AB.Mod();
  return distance;
}

int is_valued_edge( TLine* line, double x, double y )
{
  int tag = -1 ;

  if ( line->GetX1() == line->GetX2() )
  {
    double y1 = 0; double y2 = 0;
    if ( line->GetY1() >= line->GetY2() ) { y1 = line->GetY1(); y2 = line->GetY2(); }
    else { y1 = line->GetY2(); y2 = line->GetY1(); }

    if ( y <= y1 && y >= y2 ) tag = 1;
    else tag = 0;
  }

  if ( line->GetY1() == line->GetY2() )
  {
    double x1 = 0; double x2 = 0;
    if ( line->GetX1() >= line->GetX2() ) { x1 = line->GetX1(); x2 = line->GetX2(); }
    else { x1 = line->GetX2(); x2 = line->GetX1(); }

    if ( x <= x1 && x >= x2 ) tag = 1;
    else tag = 0;
  }

  if ( (line->GetX1() != line->GetX2() ) && ( line->GetY1() != line->GetY2() ) ) tag = 1;

  return tag;
}

int Rejected_edge(double x, double y, int nStrips)
{
    int Edge_tag = 0;

    //define the edge of sTGC detector
    double Det_edge_x[6] = {base_point1, base_point2, base_point2, base_point3, base_point1, base_point1};
    double Det_edge_y[6] = {base_point1, base_point1, base_point3, base_point2, base_point2, base_point1};

    int q = in_bounds_quad(x, y);
    if ( q <= 0 ) return -1;

    //define the detector edge and strip edge
    for ( int i = 1; i <= nDetEdges; i++ )
    {
        TString name = Form( "M%dEdge%d", q, i );
        string linename = name.Data();
        // double Distance = Edges[ linename ]->DistancetoPrimitive( x, y );
        // Edges[ linename ]->Print();
        double Distance = Distance_Point2Line_2D( (TLine*)Edges[ linename], x, y );
        // cout << name << " Distance = " << Distance << endl;
        if ( Distance < nStrips * 3.2 ) Edge_tag = 1;
        for ( int j = 0; j < 2; j++ )
        {
            for ( int k = 1; k <= 2; k++ )
            {
                int val_edge = -1;
                name = Form( "M%dhG%dEdge%d", q, j, k );
                linename = name.Data();
                val_edge = is_valued_edge( (TLine*)Edges[ linename], x, y );
                Distance = Distance_Point2Line_2D( (TLine*)Edges[ linename], x, y );
                if ( Distance < nStrips * 3.2 && val_edge == 1 ) Edge_tag = 1;
                name = Form( "M%dvG%dEdge%d", q, j, k );
                linename = name.Data();
                val_edge = is_valued_edge( (TLine*)Edges[ linename], x, y );
                Distance = Distance_Point2Line_2D( (TLine*)Edges[ linename], x, y );
                if ( Distance < nStrips * 3.2 && val_edge == 1 ) Edge_tag = 1;
            }
        }
    }

    return Edge_tag;

}

TGraph* Get_RealHits()
{
    TGraph* RC_Hits = new TGraph();
    int n_Hits = 0;
    for (int i = 0; i < n_Modules; i++) // too many loops // module loop
    {
        for (int j = 0; j < n_Groups; j++) // x group loop 
        {
            for (int k = 0; k < nHits; k++) // x hits loop
            {
                for (int l = 0; l < n_Groups; l++) // y group loop
                {
                    for (int n = 0; n < nHits; n++) // y hits loop
                    {
                        double x = Xhit_position[i][j][k];
                        double y = Yhit_position[i][l][n];
                        if (x < 1.e-3 || x > 536) x = 0;
                        if (y < 1.e-3 || y > 536) y = 0;
                        if (x == 0 || y == 0) continue;
                        int RealHit_flag = 0;
                        if ( x > y ) //Diagnoal Group 1, rotate π/4, use Y information
                        {
                            cout << " now in moduel " << i+1 << " x group "<< j << " y group " << l << endl;
                            cout <<  " (x,y) = " << "(" << x << "," << y << ")" <<endl;
                            TVector2 vec(x,y);
                            TVector2 vec_rot = vec.Rotate(TMath::Pi()/4);
                            cout << " rotate (x,y) = (" << vec_rot.X() << "," << vec_rot.Y() << ")" <<endl;  
                            for( auto n : Diag_Position[(i+1)*10+1] )
                            {
                                cout << "n = " << n << endl;
                                if ( abs(n-vec_rot.Y()) < 1 ) 
                                {
                                    RealHit_flag = 1;
                                }
                            }
                        }
                        if ( x < y ) //Diagnoal Group 2, rotate π/4, use X information
                        {
                            cout << " now in moduel " << i+1 << " x group "<< j << " y group " << l << endl;
                            cout <<  " (x,y) = " << "(" << x << "," << y << ")" <<endl;
                            TVector2 vec(x,y);
                            TVector2 vec_rot = vec.Rotate(-TMath::Pi()/4);
                            cout << " rotate (x,y) = (" << vec_rot.X() << "," << vec_rot.Y() << ")" <<endl;  
                            for( auto n : Diag_Position[(i+1)*10+2] )
                            {
                                cout << "n = " << n << endl;
                                if ( abs(n-vec_rot.X()) < 1 ) 
                                {
                                    RealHit_flag = 1;
                                }
                            }
                        }
                        // (x,y) need to pass the group selection and can match one of diagnoal hit 
                        // if ( RealHit_flag == 0 ) continue;
                        // if (xy_to_position_group(l,j,x,y) < 0 ) continue;
                        if (xy_to_position_group(l,j,x,y) < 0 || RealHit_flag == 0 ) continue;
                        cout << xy_to_position_group(l,j,x,y) << endl;
                        double Gx,Gy;
                        coordinate_convertion(i+1,x,y,Gx,Gy);
                        cout << "local strip x " << (int)(abs(x))/3.2 << " local strip y = " << (int)y/3.2 << endl;
                        cout << "Find hits in moduel " << i+1 << " (x,y) = " << "(" << Gx << "," << Gy << ")" <<endl;
                        RC_Hits->SetPoint(n_Hits,Gx,Gy);
                        n_Hits++;
                    }
                }
            }
        }
    }
    return RC_Hits;
}


int finder2_v3( TString inputName = "../stgc-cluster-sim/output/1DEff_test/Evts_10_0.root", TString outputName = "Cluster_output_v3_test.root")
{
    gRandom = new TRandom3();
    TH1D* h1D_Efficiency = new TH1D("1D_efficiency","1D_efficiency;efficiency;counts",110,0,1.1);
    TH1D* h1D_Efficiency_total = new TH1D("1D_efficiency_total","1D_efficiency_total;efficiency;counts",110,0,1.1);
    TH2D* hX_Missed_Hits = new TH2D("hX_Missed_Hits","hX_Missed_Hits",400,-640,640,400,-640,640);
    TH2D* hY_Missed_Hits = new TH2D("hY_Missed_Hits","hY_Missed_Hits",400,-640,640,400,-640,640);

    //open the input file
    TFile* inputFile = new TFile( inputName );
    if ( !inputFile->IsOpen() )
    {
       cout << "Failed to open the input file, Please Check it !!" << endl;
        return 0;
    }
    
    // read the histograms
    TString hisname;
    int File_idx = 0;
    for (int i = 0; i < n_Modules ; i++)
    {
        for ( int j = 0; j < n_Layers; j++ )
        {
            for ( int k = 0; k < n_Groups; k++ )
            {
                //read and save the histograms in the map
                hisname = Form("hDIG3L%d%sG%d",i+1,LayerName[j].Data(),k);
                ADC_Dis_1D[hisname.Data()] = (TH1D*)inputFile->Get(hisname);
                ADC_Dis_1D[hisname.Data()]->Print();
                ADC_Dis_1D[hisname.Data()] = (TH1D*)BKG_Substract(ADC_Dis_1D[hisname.Data()], 1); // test sample wo noise
                // ADC_Dis_1D[hisname.Data()] = (TH1D*)BKG_Substract(ADC_Dis_1D[hisname.Data()], 60);
                ADC_Der_1D[hisname.Data()] = hDerivative(ADC_Dis_1D[hisname.Data()],1);
                ADC_Der_2nd_1D[hisname.Data()] = hDerivative(ADC_Der_1D[hisname.Data()],1);

                File_idx++;
            }
        }
    }

    //read the diagnoal histograms 
    File_idx = 0;
    for (int i = 0; i < n_Modules; i++ )
    {
       for ( int j = 0; j < n_Layers; j++ )
       {
           for (int k = 0; k < n_Diag_Groups; k++)
           {
                // hDIG3L%dDiagG%d;//group 1 or 2
                // read and save the diagonal histograms
                hisname = Form("hDIG3L%dDiagG%d",i+1,k+1);
                ADC_Dis_1D_Diag[hisname.Data()] = (TH1D*)inputFile->Get(hisname);
                ADC_Dis_1D_Diag[hisname.Data()] = (TH1D*)BKG_Substract(ADC_Dis_1D_Diag[hisname.Data()], 1); // test sample wo noise
                File_idx++;
           }
       } 
    }

    int n_correct_cluster_total = 0;
    // loop all the histograms and get the 1D position, for the X and Y direction
    for ( auto nh: ADC_Dis_1D )
    {
        if ( nullptr != nh.second )
        {
            cout << endl;
            TString name = nh.first;
            cout << name << endl;
            string name1 = name.Data();
            cluster_Posi.clear();
            ADC_Channels.clear();
            Get_Cluster_Position(nh.second, cluster_Posi, ADC_Channels, 4,0); // function to get the 1D position
            string Module, group, direction;
            // hDIG3L%d%sG%d get the 1D direction
            Module = name1.substr(6,1);
            group = name1.substr(9,1);
            direction = name1.substr(7,1);
            int Module1, group1;
            Module1 = stoi(Module);
            group1 = stoi(group);
            // number of RC hist in this histograms
            cout << "number of RC hits = " << cluster_Posi.size() << endl;
            for ( int i = 0; i < cluster_Posi.size(); i++ )
            {
                double Coord = cluster_Posi.at(i);
                //save the hits from X and Y
                if (direction == "v")
                {
                    // cout << "Now in " << direction << endl;
                    Xhit_position[Module1-1][group1][i] = Coord;
                    cout << "in " << i << "th Hits, RC position = " <<Coord << endl;
                }
                if (direction == "h")
                {
                    // cout << "Now in " << direction << endl;
                    Yhit_position[Module1-1][group1][i] = Coord;
                    cout << "in " << i << "th Hits, RC position = " <<Coord << endl;
                    cout << Coord <<endl;
                }
            }

            // check the resolution of 1D clusters
            int n_correct_cluster = 0;
            // check vertiacal direction
            if ( direction == "v" )
            {
                TString graph_name = Form("gDIG1L%dvG%d",Module1,group1);
                string gname = graph_name.Data();
                TGraph* gr = (TGraph*)inputFile->Get(graph_name);
                cout << " number of MC hits = " << gr->GetN() << endl;
                for ( int i = 0; i < gr->GetN(); i++ )
                {
                    double MCx,MCy;
                    gr->GetPoint(i,MCx,MCy);
                    cout << "MCx = " << MCx << " MCy = " << MCy << endl;
                    int ReConstruct_tag = 0;
                    for( int j = 0; j < cluster_Posi.size(); j++)
                    {
                        double RCx = cluster_Posi.at(j);
                        double RCy = 0;
                        coordinate_convertion(Module1,RCx,RCy,RCx,RCy);//convert the local coordinate to global coordinate 
                        double R = abs(RCx-MCx);
                        cout << "RCx = " << RCx << " RCy = " << RCy << " R = " << R << endl;
                        if ( R < 3 * 0.1 ) //3 sigma with 1.4*3.2
                        {
                            n_correct_cluster++;
                            ReConstruct_tag = 1;
                        }
                        if ( ReConstruct_tag == 0 ) hX_Missed_Hits->Fill(MCx,MCy);

                    }
                }
                h1D_Efficiency->Fill((double)n_correct_cluster/(double)gr->GetN());
                cout << "n_correct_cluster = " << n_correct_cluster << endl;
                cout << "total clusters = " << gr->GetN() << endl;
                n_correct_cluster_total = n_correct_cluster_total+n_correct_cluster;
            }
            // horizontal direction
            if ( direction == "h" )
            {
                TString graph_name = Form("gDIG1L%dhG%d",Module1,group1);
                string gname = graph_name.Data();
                TGraph* gr = (TGraph*)inputFile->Get(graph_name);
                cout << " number of MC hits = " << gr->GetN() << endl;
                for ( int i = 0; i < gr->GetN(); i++ )
                {
                    double MCx,MCy;
                    gr->GetPoint(i,MCx,MCy);
                    cout << "MCx = " << MCx << " MCy = " << MCy << endl;
                    int ReConstruct_tag = 0;
                    for( int j = 0; j < cluster_Posi.size(); j++)
                    {
                        double RCy = cluster_Posi.at(j);
                        double RCx = 0;
                        coordinate_convertion(Module1,RCx,RCy,RCx,RCy);
                        double R = abs(RCy-MCy);
                        cout << "RCx = " << RCx << " RCy = " << RCy << " R = " << R << endl;
                        if ( R < 3*0.1 ) //3 sigma with 1.4*3.2
                        {
                            n_correct_cluster++;
                            ReConstruct_tag = 1;
                        }
                        if ( ReConstruct_tag == 0 ) hY_Missed_Hits->Fill(MCx,MCy);

                    }
                }
                n_correct_cluster_total = n_correct_cluster_total+n_correct_cluster;
                cout << "n_correct_cluster = " << n_correct_cluster << endl;
                cout << "total clusters = " << gr->GetN() << endl;
                h1D_Efficiency->Fill((double)n_correct_cluster/(double)gr->GetN());
            }

        }
    }
    // get the diagonal cluster information
    cout << "Diagnoal information started !!!!!!!!!! " << endl;
    for ( auto nh: ADC_Dis_1D_Diag )
    {
        if ( nullptr != nh.second )
        {
            cout << endl;
            TString name = nh.first;
            cout << name << endl;
            string name1 = name.Data();
            cluster_Posi.clear();
            ADC_Channels.clear();
            Get_Cluster_Position(nh.second, cluster_Posi, ADC_Channels, 4,1); // function to get the 1D position
            string Module, group, direction;
            // hDIG3L%dDiagG%d0 get the 1D direction
            Module = name1.substr(6,1);
            group = name1.substr(12,1);
            int Module1, group1;
            Module1 = stoi(Module);
            group1 = stoi(group);
            // number of RC hist in this histograms
            cout << "number of RC hits = " << cluster_Posi.size() << endl;
            for ( int i = 0; i < cluster_Posi.size(); i++ )
            {
                double Coord = cluster_Posi.at(i);
                //save the hits from X and Y
                if (group1 == 1)
                {
                    // cout << "Now in " << direction << endl;
                    Diag_hit_position[Module1-1][group1][i] = Coord;
                    Diag_Position[Module1*10+group1].push_back(Coord);
                    cout << "in " << i << "th Hits, RC position = " <<Coord << endl;
                }
                if (group1 == 2)
                {
                    // cout << "Now in " << direction << endl;
                    Diag_hit_position[Module1-1][group1][i] = Coord;
                    Diag_Position[Module1*10+group1].push_back(Coord);
                    cout << "in " << i << "th Hits, RC position = " <<Coord << endl;
                }
            }
            // How to check the diagnoal resolution?

        }
    }
    TGraph* MC_map = (TGraph*)inputFile->Get("Hits_map");
    h1D_Efficiency_total->Fill((double)n_correct_cluster_total/MC_map->GetN()*2);
    MC_map->Delete();

    int n_Hits = 0;
    /*
    TGraph* hits_map = new TGraph(); // the rc map without the combine of the edge hits
    hits_map->SetNameTitle("Hits_map_RC","Hits_map_RC");
    cout << " combining 1D to 2D" <<endl;
    // combin the points from 1D to 2D
    for (int i = 0; i < n_Modules; i++) // too many loops // module loop
    {
        for (int j = 0; j < n_Groups; j++) // x group loop 
        {
            for (int k = 0; k < nHits; k++) // x hits loop
            {
                for (int l = 0; l < n_Groups; l++) // y group loop
                {
                    for (int n = 0; n < nHits; n++) // y hits loop
                    {
                        double x = Xhit_position[i][j][k];
                        double y = Yhit_position[i][l][n];
                        if (x < 1.e-3 || x > 536) x = 0;
                        if (y < 1.e-3 || y > 536) y = 0;
                        if (x == 0 || y == 0) continue;
                        cout << " now in moduel " << i+1 << " x group "<< j << " y group " << l << endl;
                        cout <<  " (x,y) = " << "(" << x << "," << y << ")" <<endl;
                        // if xy pass the group selection, how about the diagnoal?
                        if (xy_to_position_group(l,j,x,y) < 0) continue;
                        cout << xy_to_position_group(l,j,x,y) << endl;
                        double Gx,Gy;
                        coordinate_convertion(i+1,x,y,Gx,Gy);
                        cout << "local strip x " << (int)(abs(x))/3.2 << " local strip y = " << (int)y/3.2 << endl;
                        cout << "Find hits in moduel " << i+1 << " (x,y) = " << "(" << Gx << "," << Gy << ")" <<endl;
                        hits_map->SetPoint(n_Hits,Gx,Gy);
                        n_Hits++;
                    }
                    
                }

            }
        }
    */
    TGraph* hits_map = Get_RealHits();
    hits_map->SetNameTitle("Hits_map_RC","Hits_map_RC");
    cout << " combining 1D to 2D" <<endl;
   

    TGraph* hits_map2 = new TGraph;
    hits_map2->SetNameTitle("Hits_map2_RC","Hits_map2_RC");
    vector <int> ClosePoints;
    vector <int> RemovedPoints;

    TGraph* hits_map_Diag1 = new TGraph();
    hits_map_Diag1->SetNameTitle("hits_map_Diag1","hits_map_Diag1");
    TGraph* hits_map_Diag2 = new TGraph();
    hits_map_Diag2->SetNameTitle("hits_map_Diag2","hits_map_Diag2");
    
    //save the rotated information to check the diagnoal strip resolution
    int nDiag1_Hits = 0;
    int nDiag2_Hits = 0;
    for (int i = 0; i < hits_map->GetN(); i++)
    {
        double x = 0, y = 0;
        hits_map->GetPoint(i,x,y);
        if ( x > y )
        {
            TVector2 vec(x,y);
            TVector2 vec_rot = vec.Rotate(TMath::Pi()/4);
            hits_map_Diag1->SetPoint(nDiag1_Hits,vec_rot.X(),vec_rot.Y());
            nDiag1_Hits++;
        }
        if ( x < y )
        {
            TVector2 vec(x,y);
            TVector2 vec_rot = vec.Rotate(-TMath::Pi()/4);
            hits_map_Diag2->SetPoint(nDiag2_Hits,vec_rot.X(),vec_rot.Y());
            nDiag2_Hits++;
        }
    }
    
    
    for (int i = 0; i < hits_map->GetN(); i++)
    {
        int is_Removed_point = 0;
        double x = 0; double y = 0; 
        hits_map->GetPoint(i,x,y);
        ClosePoints.clear();
        // int k = 0;
        cout << "number of RemovedPoints = " << RemovedPoints.size() << endl;
        if ( RemovedPoints.size() >= 1)
        {
            for( int j = 0;  j < RemovedPoints.size(); j++)
            {
                // if ()
                cout << j << endl;
                if ( i == RemovedPoints.at(j) ) is_Removed_point = 1;
            }
        }
        
        if ( is_Removed_point == 1) continue;

        for( int j = i+1; j < hits_map->GetN(); j++ )
        {
            double x_test = 0; double y_test = 0;
            hits_map->GetPoint(j,x_test,y_test);
            cout << " Distance = " << sqrt( (x-x_test)*(x-x_test) + (y-y_test)*(y-y_test) ) << endl;
            if ( sqrt( (x-x_test)*(x-x_test) + (y-y_test)*(y-y_test) ) < 1.56 )
            {
                ClosePoints.push_back(j);
                RemovedPoints.push_back(j);
            }
        }
        
        if ( ClosePoints.size() == 0 ) hits_map2->SetPoint(i,x,y);
        if ( ClosePoints.size() >= 1 )
        {
            for ( int j = 0 ; j < ClosePoints.size(); j ++ )
            {
                double x_test = 0; double y_test = 0;
                hits_map->GetPoint( ClosePoints.at(j) , x_test , y_test );
                x = ( x + x_test )/ 2;
                y = ( y + y_test )/ 2;
            }
            hits_map2->SetPoint(i,x,y);
        }
        cout << "debug " << endl;
    }

    TFile* outFile = new TFile( outputName.Data(), "recreate" );
    for ( auto nh : ADC_Dis_1D ){
        if ( nullptr != nh.second )
        {
            TString name = nh.first;
            name = name;
            nh.second->Write();
            // cout << name << endl; 
            ADC_Der_1D[name.Data()]->Write();
            ADC_Der_2nd_1D[name.Data()]->Write();
            }
    }
    h1D_Efficiency->Write();
    h1D_Efficiency_total->Write();
    hY_Missed_Hits->Write();
    hX_Missed_Hits->Write();
    hits_map->Write();
    hits_map2->Write();
    hits_map_Diag1->Write();
    hits_map_Diag2->Write();
    outFile->Write();
    outFile->Close();

    return 1;

}