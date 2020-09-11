#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TString.h"
#include "TGraph.h"

#include <iostream>
#include <map>
#include <vector>
#include <algorithm>

using namespace std;

const int Kernel_size = 1;
const int n_Modules = 4;
const int n_Layers = 2;
const int n_Groups = 3;
const int nHits = 100;
map <string, TH1D*> ADC_Dis_1D; // 1D ADC distribution
map <string, TH1D*> ADC_Der_1D; // Derivative of the 1D ADC distribution
map <string, TH1D*> ADC_Der_2nd_1D; // 2nd order Derivative of the 1D ADC distribution
vector <int> minimum;// bin number of maxmium and minimum in "signal region" of each layer 
vector <int> maximum;
vector <int> start_point;// bin number the start point and end point of every "signal region" of each layer
vector <int> end_point;
map <string, vector<int> > minimum_DiffLayer; //record the information of every layer
map <string, vector<int> > maximum_DiffLayer;
map <string, vector<int> > start_point_DiffLayer;
map <string, vector<int> > end_point_DiffLayer;
multimap < int ,int > ADC_per_Channels;
vector <int> ADC_Channels;
vector <double> cluster_Posi;
double Xhit_position[n_Modules][n_Groups][nHits] = {0.0};
double Yhit_position[n_Modules][n_Groups][nHits] = {0.0};

TString LayerName[2] = {"v","h"};

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

int Get_Cluster_Position( TH1D* hGroup, vector < double > &cluster_posi, vector <int> &vec_ADC_per_Channels , int cluster_width = 4 ) // cluster_width is the strip number to search, defult is 4
{
    vec_ADC_per_Channels.clear();
    for ( int i = 1; i < hGroup->GetNbinsX(); i++)
    {
        vec_ADC_per_Channels.push_back(hGroup->GetBinContent(i));
    }

    sort(vec_ADC_per_Channels.rbegin(),vec_ADC_per_Channels.rend());

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
        while( right_tag == 1 && left_tag == 1 ) // find one cluster
        {
            if (i == 0) {sumADC+=bc; strip_ADC = strip_ADC+bc*((bin*3.2)-1.35); hGroup->SetBinContent(bin,0); i++; j++; }
            if ( hGroup->GetBinContent(bin+i) < bc_right && hGroup->GetBinContent(bin+i) > 0 ) 
            {
                sumADC+=bc_right; 
                strip_ADC = strip_ADC+bc_right*(((bin+i)*3.2)-1.35); 
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
                strip_ADC = strip_ADC+bc_left*(((bin-j)*3.2)-1.35); 
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
            cout << "i = " << i <<" j = " << j << endl;
            cluster_posi.push_back(strip_ADC/sumADC);
        }
        else
        {
            if ( bin == 166 || bin == 167 || bin == 0 || bin == 1)
            {
                cout << "i = " << i <<" j = " << j << endl;
                if ( bin == 166 || bin == 167 )  cluster_posi.push_back(bin*3.2+1.35);
                if ( bin == 0 || bin == 1 )  cluster_posi.push_back(bin*3.2+1.35);
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
        strip_ADC = strip_ADC+bc*((i*3.2)-1.35);// because we used center of strip to calculate the position
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


int finder2_v3( TString inputName = "../stgc-cluster-sim/output_10Evts_minimumADC.root", TString outputName = "Cluster_output_v3.root")
{
    gRandom = new TRandom3();

    TFile* inputFile = new TFile( inputName );
    if ( !inputFile->IsOpen() )
    {
       cout << "Failed to open the input file, Please Check it !!" << endl;
        return 0;
    }
    
    TString hisname;
    int File_idx = 0;
    for (int i = 0; i < n_Modules ; i++)
    {
        for ( int j = 0; j < n_Layers; j++ )
        {
            for ( int k = 0; k < n_Groups; k++ )
            {
                hisname = Form("hDIG3L%d%sG%d",i+1,LayerName[j].Data(),k);
                ADC_Dis_1D[hisname.Data()] = (TH1D*)inputFile->Get(hisname);
                ADC_Dis_1D[hisname.Data()]->Print();
                ADC_Dis_1D[hisname.Data()] = (TH1D*)BKG_Substract(ADC_Dis_1D[hisname.Data()], 0); // test sample wo noise
                // ADC_Dis_1D[hisname.Data()] = (TH1D*)BKG_Substract(ADC_Dis_1D[hisname.Data()], 60);
                ADC_Der_1D[hisname.Data()] = hDerivative(ADC_Dis_1D[hisname.Data()],1);
                ADC_Der_2nd_1D[hisname.Data()] = hDerivative(ADC_Der_1D[hisname.Data()],1);

                File_idx++;
            }
        }
    }

    for ( auto nh: ADC_Dis_1D )
    {
        if ( nullptr != nh.second )
        {
            TString name = nh.first;
            cout << name << endl;
            string name1 = name.Data();
            cluster_Posi.clear();
            ADC_Channels.clear();
            Get_Cluster_Position(nh.second, cluster_Posi, ADC_Channels, 4);
            string Module, group, direction;
            // hDIG3L%d%sG%d
            Module = name1.substr(6,1);
            group = name1.substr(9,1);
            direction = name1.substr(7,1);
            int Module1, group1;
            Module1 = stoi(Module);
            group1 = stoi(group);
            for ( int i = 0; i < cluster_Posi.size(); i++ )
            {
                double Coord = cluster_Posi.at(i);
                if (direction == "v")
                {
                    Xhit_position[Module1-1][group1][i] = Coord;
                    cout << Coord << endl;
                }
                if (direction == "h")
                {
                    Yhit_position[Module1-1][group1][i] = Coord;
                    cout << Coord <<endl;
                }
            }
        }
    }

    int n_Hits = 0;
    TGraph* hits_map = new TGraph();
    hits_map->SetNameTitle("Hits_map_RC","Hits_map_RC");
    cout << " combining 1D to 2D" <<endl;
    for (int i = 0; i < n_Modules; i++) // too many loops
    {
        for (int j = 0; j < n_Groups; j++) // x group
        {
            for (int k = 0; k < nHits; k++)
            {
                for (int l = 0; l < n_Groups; l++) // y group
                {
                    for (int n = 0; n < nHits; n++)
                    {
                        double x = Xhit_position[i][j][k];
                        double y = Yhit_position[i][l][n];
                        if (x < 1.e-3 || x > 536) x = 0;
                        if (y < 1.e-3 || y > 536) y = 0;
                        if (x == 0 || y == 0) continue;
                        cout << " now in moduel " << i+1 << " x group "<< j << " y group " << l << endl;
                        cout <<  " (x,y) = " << "(" << x << "," << y << ")" <<endl;
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
    }
    TGraph* hits_map2 = new TGraph;
    hits_map2->SetNameTitle("Hits_map2_RC","Hits_map2_RC");
    vector <int> ClosePoints;
    vector <int> RemovedPoints;
    
    
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
    hits_map->Write();
    hits_map2->Write();
    outFile->Write();
    outFile->Close();

    return 1;

}