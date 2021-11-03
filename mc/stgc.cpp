//rule of strip number 
// 0/0/0/0/000
// Disk(1-4)/module(1-4)/layer(1-4)/Row(1-3)/strip()
// now the rule of layer is lower diagnoal is 1, x is 2, y is 3, upper diagnoal is 4; 

#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TLine.h"
#include "TVector2.h"
#include "TTree.h"

#include <iostream>
#include <map>
#include <ctime>
#include <algorithm>

#define LOGURU_IMPLEMENTATION 1
#include "StsTGCTreeStructure.h"
#include "loguru.h"
#include "../Public_input/Parameters.h"

StsTGCData mEvtData;
TTree* mEvtTree;

const float pent_base = 537.0;    // in mm
const float pent_nib = 179.0;     // in mm
const float pent_base_sqrt = pent_base/TMath::Sqrt(2);
const float pent_nib_sqrt = pent_nib/TMath::Sqrt(2);
const float pent_shift = 101.6;   // in mm
const int n_clusters_to_gen = 10; // Number of clusters in this "event"
const float noise_prob = 0.2;     // 0 -1, 1 is max noise
// const float noise_level = 1; // in ADC units (pre integration)
const float noise_level = 1.e-10;  // in ADC units (pre integration)
const float strip_pitch = 3.2;     // digitize strip pitch
const size_t sat_above = 1024;     // ADC
const float cluster_max_adc = 120; // ADC (pre integration)
int test_edge = 1;
int Edge_Evt = 1;
const double pi = TMath::Pi();

int debug = 0; // 0 do not print out the debug information. 1 print the debug information
int nChannel = 0;
int nDiag_1_Hits = 0;
int nDiag_2_Hits = 0;


TF1 *fClusterProfile = nullptr;
TF1 *fBreitWigner = nullptr;

using namespace std;
map<string, TH1 *> hist;
map<string, TGraph *> graph;
map<string, TLine *> Edges;
TGraph *hits_map = new TGraph();
TH1D* hDelta_BCID = new TH1D("hDelta_BCID",";#Delta_{BCID};nCounts",10, -5, 5);
TGraph* Diag_hits1_map = new TGraph();
TGraph* Diag_hits2_map = new TGraph();

Double_t lorentzianPeak(Double_t *x, Double_t *par)
{
    return (0.5 * par[0] * par[1] / TMath::Pi()) / TMath::Max(1.e-10, (x[0] - par[2]) * (x[0] - par[2]) + .25 * par[1] * par[1]);
}

bool in_bounds(float x, float y)
{
    if (x > 0 && x < pent_nib && y < pent_base && y > 0)
        return true;
    if (x > 0 && x < pent_base && y < pent_nib && y > 0)
        return true;

    float dy = pent_base - pent_nib;
    float dx = pent_nib - pent_base;
    float yedge = (dy / dx) * (x - pent_nib) + pent_base;
    // LOG_F( INFO, "yedge = %f @ x = %f", yedge, x );

    if (y < yedge && y < pent_base && x < pent_base && x > 0 && y > 0)
        return true;

    return false;
}

int in_bounds_diagonal( float x, float y )
{
    //lower part 
    if (x > 0 && x < pent_nib && y < pent_base && y > 0 && x > y)
        return 1;
    if (x > 0 && x < pent_base && y < pent_nib && y > 0 && x > y)
        return 1;

    float dy = pent_base - pent_nib;
    float dx = pent_nib - pent_base;
    float yedge = (dy / dx) * (x - pent_nib) + pent_base;
    // LOG_F( INFO, "yedge = %f @ x = %f", yedge, x );

    if (y < yedge && y < pent_base && x < pent_base && x > 0 && y > 0 && x > y)
        return 1;

    // upper part
    if (x > 0 && x < pent_nib && y < pent_base && y > 0 && x < y)
        return 2;
    if (x > 0 && x < pent_base && y < pent_nib && y > 0 && x < y)
        return 2;

    dy = pent_base - pent_nib;
    dx = pent_nib - pent_base;
    yedge = (dy / dx) * (x - pent_nib) + pent_base;
    // LOG_F( INFO, "yedge = %f @ x = %f", yedge, x );

    if (y < yedge && y < pent_base && x < pent_base && x > 0 && y > 0 && x < y)
        return 2;
    
    return 0;
}

int in_bounds_quad(float x, float y)
{
    if (in_bounds(x, y))
    { // top right
        return 1;
    }
    else if (in_bounds(-x, y))
    { //top left
        return 2;
    }
    else if (in_bounds(-x - pent_shift, -y))
    { // bottom left
        return 3;
    }
    else if (in_bounds(x - pent_shift, -y))
    { // bottom right
        return 4;
    }

    return 0;
}
int diagonal_part(float x, float y)
{
    double x_local = -1, y_local = -1;
    if ( in_bounds_quad(x,y) == 1 )// top right
        {x_local = x; y_local = y; return in_bounds_diagonal(x_local,y_local);}
    if ( in_bounds_quad(x,y) == 2 )// top left
        {x_local = -x; y_local =y; return in_bounds_diagonal(x_local,y_local);}
    if ( in_bounds_quad(x,y) == 3 )// bottom left
        {x_local = -x - pent_shift; y_local = -y; return in_bounds_diagonal(x_local,y_local);}
    if ( in_bounds_quad(x,y) == 4 )// bottom right
        {x_local = x - pent_shift; y_local = -y; return in_bounds_diagonal(x_local,y_local);}

    return 0;
}

int map_to_strip(float x, float y, int &sx, int &sy)
{
    int q = in_bounds_quad(x, y);
    if (q == 1 || q == 2)
    {
        sx = x / strip_pitch;
        sy = y / strip_pitch;
    }
    else if (q == 3)
    {
        sx = (x) / strip_pitch;
        sy = y / strip_pitch;
    }
    else if (q == 4)
    {
        sx = (x) / strip_pitch;
        sy = y / strip_pitch;
    }
    return q;
}

int map_to_local_strip(float x, float y, int &sx, int &sy)
{
    int q = in_bounds_quad(x, y);
    if (q == 1 || q == 2)
    {
        sx = abs(x / strip_pitch);
        sy = abs(y / strip_pitch);
    }
    else if (q == 3)
    {
        sx = abs((x + pent_shift) / strip_pitch);
        sy = abs(y / strip_pitch);
    }
    else if (q == 4)
    {
        sx = abs((x - pent_shift) / strip_pitch);
        sy = abs(y / strip_pitch);
    }
    return q;
}

int map_to_local_strip_diagonal(float x, float y)
{
    int d_p = diagonal_part(x, y);
    int q = in_bounds_quad(x, y);
    int s;
    double x_rot, y_rot;
    

    if ( q == 1 )
    {
        TVector2 vec(x,y);
        if (d_p == 1) 
        {
            TVector2 vec_rot = vec.Rotate(TMath::Pi()/4);
            s = vec_rot.Y()/strip_pitch;
        }
        if (d_p == 2)
        {
            TVector2 vec_rot = vec.Rotate(-TMath::Pi()/4);
            s = vec_rot.X()/strip_pitch;
        }
    }
    if ( q == 2 )
    {
        TVector2 vec(-x,y);
        if (d_p == 1) 
        {
            TVector2 vec_rot = vec.Rotate(TMath::Pi()/4);
            s = vec_rot.Y()/strip_pitch;
        }
        if (d_p == 2)
        {
            TVector2 vec_rot = vec.Rotate(-TMath::Pi()/4);
            s = vec_rot.X()/strip_pitch;
        }
    }
    if ( q == 3 )
    {
        TVector2 vec(-(x+pent_shift),-y);
        if (d_p == 1) 
        {
            TVector2 vec_rot = vec.Rotate(TMath::Pi()/4);
            s = vec_rot.Y()/strip_pitch;
        }
        if (d_p == 2)
        {
            TVector2 vec_rot = vec.Rotate(-TMath::Pi()/4);
            s = vec_rot.X()/strip_pitch;
        }
    }
    if ( q == 4 )
    {
        TVector2 vec(x-pent_shift,-y);
        if (d_p == 1) 
        {
            TVector2 vec_rot = vec.Rotate(TMath::Pi()/4);
            s = vec_rot.Y()/strip_pitch;
        }
        if (d_p == 2)
        {
            TVector2 vec_rot = vec.Rotate(-TMath::Pi()/4);
            s = vec_rot.X()/strip_pitch;
        }
    }
    return s;
}

int stripGroup(int sx, int sy)
{
    if ((sx < 55 && sy < 150) || sy >= 150)
        return 0;
    if ((sx < 110 && sy < 95) || sy >= 95)
        return 1;
    return 2;
}

void fill_cluster_GEN(float x, float y)
{
    TH2 *hGENx = (TH2 *)hist["hGEN0x"];
    TH2 *hGENy = (TH2 *)hist["hGEN0y"];
    TH2 *hGEN_diagonalx = (TH2 *)hist["hGEN0_diagonalx"];
    TH2 *hGEN_diagonaly = (TH2 *)hist["hGEN0_diagonaly"];
    TH2 *hDIGx = (TH2 *)hist["hDIG1x"];
    TH2 *hDIGy = (TH2 *)hist["hDIG1y"];
    TH2 *hDIG_diagonalx = (TH2 *)hist["hDIG1_diagonalx"];
    TH2 *hDIG_diagonaly = (TH2 *)hist["hDIG1_diagonaly"];
    TH2 *hDIGxRecorder = (TH2 *)hist["hDIG1xRecorder"];
    TH2 *hDIGyRecorder = (TH2 *)hist["hDIG1yRecorder"];
    // int bx = hGENx->GetXaxis()->FindBin( x );
    // int by = hGENx->GetYaxis()->FindBin( y );

    TF1* ftime = new TF1("ftime","-gaus+3",-10,10);
    ftime->SetParameters(3,0,2);

    float maxADC = 0;
    while (maxADC < noise_level * 7)
    {
        maxADC = gRandom->Rndm() * cluster_max_adc;
    }
    // maxADC = 100;

    // float maxADC = gRandom->Rndm() * cluster_max_adc;
    LOG_F(INFO, "Cluster at (%f, %f) w maxADC = %f", x, y, maxADC);

    // int wSpan = 20;
    // for ( int i = bx - wSpan; i < bx + wSpan; i++ ){
    //     for ( int j = by - wSpan; j < by + wSpan; j++ ){
    //         float bcx = hGEN->GetXaxis()->GetBinCenter( i );
    //         float bcy = hGEN->GetYaxis()->GetBinCenter( j );

    //         int q = in_bounds_quad( bcx, bcy );
    //         if ( q <= 0) continue;

    //         float r = sqrt( pow( bcx - x, 2 ) + pow( bcy - y, 2 ) );
    //         float v = fClusterProfile->Eval( r ) * maxADC;
    //         // float v = fBreitWigner->Eval(r) * maxADC;
    //         // LOG_F( INFO, "---- Cluster ADC = %f @ r = %f", v, r );
    //         hGEN->Fill( bcx, bcy, v );

    //         int sx = -1, sy = -1;
    //         map_to_strip( bcx, bcy, sx, sy );
    //         hDIG->Fill( sx, sy, v );

    //         map_to_local_strip( bcx, bcy, sx, sy );
    //         string hname = TString::Format( "hDIG1L%d", q ).Data();
    //         int xG = stripGroup( sx, sy );
    //         int yG = stripGroup( sy, sx ); // not a typo, just reuse function with rotation
    //         if ( hist.count( hname ) > 0 ){
    //             ((TH2*)hist[ hname ]) -> Fill( sx, sy, v );

    //             hname = TString::Format( "hDIG1L%dhG%d", q, xG ).Data();
    //             ((TH2*)hist[ hname ]) -> Fill( sx, sy, v );
    //             if ( i == bx && j == by )
    //             {
    //                 string gname = TString::Format( "gDIG1L%dhG%d", q, xG ).Data();
    //                 ((TGraph*)graph[ gname ])->SetPoint( graph[ gname ]->GetN(),x,y );
    //             }

    //             hname = TString::Format( "hDIG1L%dvG%d", q, yG ).Data();
    //             ((TH2*)hist[ hname ]) -> Fill( sx, sy, v );
    //             if ( i == bx && j == by )
    //             {
    //                 string gname = TString::Format( "gDIG1L%dvG%d", q, yG ).Data();
    //                 ((TGraph*)graph[ gname ])->SetPoint( graph[ gname ]->GetN(),x,y );
    //             }

    //         }

    //     }
    // }
    //try to seperate x and y direction sample

    int bx = 0, by = 0;
    map_to_strip(x, y, bx, by);

    int wSpan = 10;
    int BCID = gRandom->Rndm() * 1023;
    for (int i = bx - wSpan; i < bx + wSpan; i++) // sample in X direction
    {
        float bcx = i * 3.2 + 1.6;
        float bLowEdge = i * 3.2;
        float bHighEdge = (i + 1) * 3.2;
        
        // float bcx = hGENx->GetXaxis()->GetBinCenter( i );
        int q = in_bounds_quad(bcx, y); // avoid pass the lower or upper limit
        if (q <= 0)
            continue;

        float r = bcx - x;
        float rL = bLowEdge - x;
        float rH = bHighEdge - x;// will use in intergation method 
        float v = fClusterProfile->Eval(r) * maxADC;
        // float v = fClusterProfile->Integral( rL, rH ) * maxADC;
        // cout << i <<" th bin " << endl;
        // cout << "bin center = " << bcx << " r = " << r << " v = " << v << endl;

        hGENx->Fill( bcx, y, v );


        int sx = -1, sy = -1;
        map_to_strip(bcx, y, sx, sy);
        hDIGx->Fill(sx, sy, v);
        hDIGxRecorder->Fill(sx, sy, r);

        map_to_local_strip(bcx, y, sx, sy);
        // cout << "sx = " << sx << " sy = " << sy << endl;
        string hname = TString::Format("hDIG1L%d", q).Data();
        int xG = stripGroup(sy, sx); // not a typo, just reuse function with rotation

        //sample time information
        int delta_BCID;
        delta_BCID = -3;
        // cout << "origin = "<< BCID << endl; 
        // cout <<" old = " << mEvtData.BCID.back() << endl;
        int delta_bin = abs(i-bx);
        delta_BCID = (int)ftime->Eval(delta_bin);
        int smear = gRandom->Uniform(-1.5,1.5);
        delta_BCID = delta_BCID+smear;
        if ( delta_BCID > 0 ) delta_BCID = delta_BCID-1;
        //check this channel have or did not have signal before
        int is_usedCh = std::count(mEvtData.Channel.begin(), mEvtData.Channel.end(), 1.e7+q*1.e6+2*1.e5+xG*1.e4+i);
        if( !is_usedCh )
        {
            mEvtData.Channel.push_back(1.e7+q*1.e6+2*1.e5+xG*1.e4+i);
            mEvtData.ADC.push_back((int)v);
            mEvtData.BCID.push_back(BCID+delta_BCID);
        }
        hDelta_BCID->Fill(delta_BCID);


        if (hist.count(hname) > 0)
        {
            ((TH2 *)hist[hname])->Fill(sx, sy, v);

            hname = TString::Format("hDIG1L%dvG%d", q, xG).Data();
            ((TH2 *)hist[hname])->Fill(sx, sy, v);
            if (i == bx)
            {
                string gname = TString::Format("gDIG1L%dvG%d", q, xG).Data();
                ((TGraph *)graph[gname])->SetPoint(graph[gname]->GetN(), x, y);
            }
        }
        
        nChannel++;
    }
    for (int j = by - wSpan; j < by + wSpan; j++)
    {
        float bcy = j * 3.2 + 1.6;
        float bLowEdge = j * 3.2;
        float bHighEdge = (j + 1) * 3.2;
        // float bcy = hGENy->GetYaxis()->GetBinCenter( j );
        int q = in_bounds_quad(x, bcy);
        if (q <= 0)
            continue;

        float r = bcy - y;
        float rL = bLowEdge - y;
        float rH = bHighEdge - y;
        float v = fClusterProfile->Eval(r) * maxADC;
        // float v = fClusterProfile->Integral( rL, rH ) * maxADC;

        hGENy->Fill( x, bcy, v );

        int sx = -1, sy = -1;
        map_to_strip(x, bcy, sx, sy);
        hDIGy->Fill(sx, sy, v);
        hDIGyRecorder->Fill(sx, sy, r);

        map_to_local_strip(x, bcy, sx, sy);
        string hname = TString::Format("hDIG1L%d", q).Data();
        int yG = stripGroup(sx, sy);
        
        //sample time information
        int delta_BCID;
        delta_BCID = -3;
        int delta_bin = abs(j-bx);
        delta_BCID = (int)ftime->Eval(delta_bin);
        int smear = gRandom->Uniform(-1.5,1.5);
        delta_BCID = delta_BCID+smear;
        if ( delta_BCID > 0 ) delta_BCID = delta_BCID-1;
        int is_usedCh = std::count(mEvtData.Channel.begin(), mEvtData.Channel.end(),  1.e7+q*1.e6+3*1.e5+yG*1.e4+j);
        if ( !is_usedCh )
        {
            mEvtData.Channel.push_back(1.e7+q*1.e6+3*1.e5+yG*1.e4+j);
            mEvtData.ADC.push_back((int)v);
            mEvtData.BCID.push_back(BCID+delta_BCID);
        }
        hDelta_BCID->Fill(delta_BCID);

        if (hist.count(hname) > 0)
        {
            ((TH2 *)hist[hname])->Fill(sx, sy, v);

            hname = TString::Format("hDIG1L%dhG%d", q, yG).Data();
            ((TH2 *)hist[hname])->Fill(sx, sy, v);
            if (j == by)
            {
                string gname = TString::Format("gDIG1L%dhG%d", q, yG).Data();
                ((TGraph *)graph[gname])->SetPoint(graph[gname]->GetN(), x, y);
            }
        }
    }
    // the diagonal histograms, It good to have 2D plots, but how to draw it?
    // the calculation is after rotation 
    int diag_bin = map_to_local_strip_diagonal(x,y);
    // get (x,y) belong to which part and which module
    int q = in_bounds_quad(x, y);
    int d_p = diagonal_part(x, y);
    //lower part, y is the diagnoal
    cout << "diag_bin = " << diag_bin << endl;
    cout << "d_p = " << d_p << endl;
    if ( d_p == 1 )
    {
        TVector2 vec;
        if ( q == 1 ) vec.Set(x,y);
        if ( q == 2 ) vec.Set(-x,y);
        if ( q == 3 ) vec.Set(-(x+pent_shift),-y);
        if ( q == 4 ) vec.Set(x-pent_shift,-y);
        TVector2 vec_rot = vec.Rotate(pi/4);
        Diag_hits1_map->SetPoint(nDiag_1_Hits,vec_rot.X(),vec_rot.Y());
        nDiag_1_Hits++;
        for( int i = diag_bin - wSpan; i < diag_bin + wSpan; i++)
        {
            float bc = i * 3.2 + 1.6;// bin center at rotated frame.
            float bLowEdge = i * 3.2;
            float bHighEdge = (i + 1) * 3.2;
            // float bcy = hGENy->GetYaxis()->GetBinCenter( j );
            if ( i < 1 || i > 151 ) continue;

            float r = bc - vec_rot.Y();
            float rL = bLowEdge - vec_rot.Y();
            float rH = bHighEdge - vec_rot.Y();
            float v = fClusterProfile->Eval(r) * maxADC;

            int sx = vec_rot.X()/3.2;
            string hname = TString::Format("hDIG1L%dDiagG%d", q, d_p).Data();
            if (debug == 1 ) cout << "Filling histogram : " << hname << endl;

            //sample the time information
            int delta_BCID;
            delta_BCID = -3;
            int delta_bin = abs(i-bx);
            delta_BCID = (int)ftime->Eval(delta_bin);
            int smear = gRandom->Uniform(-1.5,1.5);
            delta_BCID = delta_BCID+smear;
            if ( delta_BCID > 0 ) delta_BCID = delta_BCID-1;
            int is_usedCh = std::count(mEvtData.Channel.begin(),mEvtData.Channel.end(), 1.e7+q*1.e6+1*1.e5+1*1.e4+i);
            if( !is_usedCh )
            {
                mEvtData.Channel.push_back(1.e7+q*1.e6+1*1.e5+1*1.e4+i);
                mEvtData.ADC.push_back((int)v);
                mEvtData.BCID.push_back(BCID+delta_BCID);
            }
            hDelta_BCID->Fill(delta_BCID);

            if (hist.count(hname) > 0)
            {
                if (debug == 1 ) 
                { 
                    cout << "r = " << r << " : (bc-y) (" << bc << "-" <<  vec_rot.Y() << ")" << endl;
                    cout << " x strip = " << sx << " y strip = " << i << " v = " << v << endl;
                }
                ((TH2 *)hist[hname])->Fill(sx, i, v);
            }
        }
    }

    if ( d_p == 2 )
    {
        TVector2 vec;
        if ( q == 1 ) vec.Set(x,y);
        if ( q == 2 ) vec.Set(-x,y);
        if ( q == 3 ) vec.Set(-(x+pent_shift),-y);
        if ( q == 4 ) vec.Set(x-pent_shift,-y);
        TVector2 vec_rot = vec.Rotate(-pi/4);
        Diag_hits2_map->SetPoint(nDiag_2_Hits,vec_rot.X(),vec_rot.Y());
        nDiag_2_Hits++;
        for( int i = diag_bin - wSpan; i < diag_bin + wSpan; i++)
        {
            float bc = i * 3.2 + 1.6;// bin center at rotated frame.
            float bLowEdge = i * 3.2;
            float bHighEdge = (i + 1) * 3.2;
            // float bcy = hGENy->GetYaxis()->GetBinCenter( j );
            if ( i < 1 || i > 151 ) continue;

            float r = bc - vec_rot.X();
            float rL = bLowEdge - vec_rot.X();
            float rH = bHighEdge - vec_rot.X();
            float v = fClusterProfile->Eval(r) * maxADC;

            int sy = vec_rot.Y()/3.2;
            string hname = TString::Format("hDIG1L%dDiagG%d", q, d_p).Data();
            if (debug == 1 ) cout << "Filling histogram : " << hname << endl;

            //sample the time information
            int delta_BCID;
            delta_BCID = -3;
            int delta_bin = abs(i-bx);
            delta_BCID = (int)ftime->Eval(delta_bin);
            int smear = gRandom->Uniform(-1.5,1.5);
            delta_BCID = delta_BCID+smear;
            if ( delta_BCID > 0 ) delta_BCID = delta_BCID-1;
            int is_usedCh = std::count(mEvtData.Channel.begin(),mEvtData.Channel.end(), 1.e7+q*1.e6+4*1.e5+1*1.e4+i);
            if( !is_usedCh )
            {
                mEvtData.Channel.push_back(1.e7+q*1.e6+4*1.e5+1*1.e4+i);
                mEvtData.ADC.push_back((int)v);
                mEvtData.BCID.push_back(BCID+delta_BCID);
            }

            if (hist.count(hname) > 0)
            {
                if (debug == 1 ) 
                { 
                    cout << "r = " << r << " : (bc-x) (" << bc << "-" <<  vec_rot.X() << ")" << endl;
                    cout << " x strip = " << i << " y strip = " << sy << " v = " << v << endl;
                }
                ((TH2 *)hist[hname])->Fill(i, sy, v);
            }
        }
    }

}

void saturate(TH1 *h)
{
    for (int i = 1; i < h->GetNbinsX(); i++)
    {
        if (h->GetBinContent(i) < sat_above)
            continue;

        h->SetBinContent(i, (int)(sat_above - 1));
    }
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

void BookTree()
{
    cout << " Booking the event tree " << endl;
    mEvtTree = new TTree("sTGC_DataDst","sTGC_DataDst");
    mEvtTree->SetAutoSave(100000000); // 100 MB

    // Event information
    mEvtTree->Branch("mEvtID", &mEvtData.mEventId, "mEvtID/I");
    mEvtTree->Branch("nFireChannel", &mEvtData.nFireChannel, "nFireChannel/I");

    // Channel information
    mEvtTree->Branch("Channel", &mEvtData.Channel);
    mEvtTree->Branch("BCID", &mEvtData.BCID);
    mEvtTree->Branch("ADC", &mEvtData.ADC);

}

int main(int argc, char **argv)
{

    loguru::init(argc, argv);
    // loguru::g_stderr_verbosity = 8;
    //add argv as input
    if (argc != 3)
        cout << "2 arguements is need , 1 is outfile name, 2 is logname" << endl;
    loguru::add_file(argv[2], loguru::Truncate, loguru::Verbosity_MAX);
    // loguru::add_file("everything.log", loguru::Truncate, loguru::Verbosity_MAX);

    LOG_SCOPE_FUNCTION(INFO);

    gRandom = new TRandom3();
    gRandom->SetSeed((int)time(NULL));
    initEdges();
    BookTree();
    memset(&mEvtData, 0, sizeof(mEvtData));

    hist["hquad1"] = new TH2F("hquad1", ";x;y", 600, 0, 600, 600, 0, 600);
    hist["hstgc"] = new TH2F("hstgc", ";x;y", 1400, -700, 700, 1200, -600, 600);

    hist["hGEN0x"] = new TH2F("hGEN0x", ";x;y", 1400, -700, 700, 1200, -600, 600);
    hist["hGEN0y"] = new TH2F("hGEN0y", ";x;y", 1400, -700, 700, 1200, -600, 600);
    hist["hGEN0_diagonalx"] = new TH2F("hGEN0_diagonalx", ";x;y", 1400, -700, 700, 1200, -600, 600);
    hist["hGEN0_diagonaly"] = new TH2F("hGEN0_diagonaly", ";x;y", 1400, -700, 700, 1200, -600, 600);
    hist["hGEN1"] = new TH2F("hGEN1", ";x;y", 1400, -700, 700, 1200, -600, 600);
    hist["hGEN2"] = new TH2F("hGEN2", ";x;y", 1400, -700, 700, 1200, -600, 600);
    hist["hDIG0"] = new TH2F("hDIG0", ";strip x; strip y", 400, -200, 200, 400, -200, 200);
    hist["hDIG1x"] = new TH2F("hDIG1x", ";strip x; strip y", 400, -200, 200, 400, -200, 200);
    hist["hDIG1y"] = new TH2F("hDIG1y", ";strip x; strip y", 400, -200, 200, 400, -200, 200);
    hist["hDIG1_diagonalx"] = new TH2F("hDIG1_diagonalx", ";strip x; strip y", 400, -200, 200, 400, -200, 200);
    hist["hDIG1_diagonaly"] = new TH2F("hDIG1_diagonaly", ";strip x; strip y", 400, -200, 200, 400, -200, 200);
    hist["hDIG1xRecorder"] = new TH2F("hDIG1xRecorder", ";strip x; strip y", 400, -200, 200, 400, -200, 200);
    hist["hDIG1yRecorder"] = new TH2F("hDIG1yRecorder", ";strip x; strip y", 400, -200, 200, 400, -200, 200);
    hist["hDIG2"] = new TH2F("hDIG2", ";strip x; strip y", 400, -200, 200, 400, -200, 200);
    hist["hGEN2_Diagnoal_1"] = new TH2F("hGEN2_Diagnoal_1", ";x;y", 1400, -700, 700, 1200, -600, 600);
    hist["hGEN2_Diagnoal_2"] = new TH2F("hGEN2_Diagnoal_2", ";x;y", 1400, -700, 700, 1200, -600, 600);
    hist["hGEN2_Diagnoal_rotate_1"] = new TH2F("hGEN2_Diagnoal_rotate_1", ";x;y", 1400, -700, 700, 1200, -600, 600);
    hist["hGEN2_Diagnoal_rotate_2"] = new TH2F("hGEN2_Diagnoal_rotate_2", ";x;y", 1400, -700, 700, 1200, -600, 600);
    // add the histogram 
    for (int i = 0; i < 3; i++)
    {
        for (int iMod = 1; iMod < 5; iMod++)
        {
            string n = TString::Format("hDIG%dL%d", i, iMod).Data();
            LOG_F(INFO, "Making : %s", n.c_str());
            hist[n] = new TH2F(n.c_str(), TString::Format("Module %d;local strip x; local strip y", iMod), 200, 0, 200, 200, 0, 200);
            for (int ixG = 0; ixG < 3; ixG++)
            {
                string n = TString::Format("hDIG%dL%dhG%d", i, iMod, ixG).Data();
                hist[n] = new TH2F(n.c_str(), TString::Format("Module %d stripGroup %d;local strip x; local strip y", iMod, ixG), 200, 0, 200, 200, 0, 200);
                if (i == 1)
                {
                    string gname = TString::Format("gDIG%dL%dhG%d", i, iMod, ixG).Data();
                    graph[gname] = new TGraph();
                    graph[gname]->SetNameTitle(gname.data(), gname.data());
                }
            }
            for (int iyG = 0; iyG < 3; iyG++)
            {
                string n = TString::Format("hDIG%dL%dvG%d", i, iMod, iyG).Data();
                hist[n] = new TH2F(n.c_str(), TString::Format("Module %d yGroup %d;local strip x; local strip y", iMod, iyG), 200, 0, 200, 200, 0, 200);
                if (i == 1)
                {
                    string gname = TString::Format("gDIG%dL%dvG%d", i, iMod, iyG).Data();
                    graph[gname] = new TGraph();
                    graph[gname]->SetNameTitle(gname.data(), gname.data());
                }
            }
            // diagnoal only have 2 group
            for ( int iDiagG = 0; iDiagG < 2; iDiagG++)
            {
                string n = TString::Format("hDIG%dL%dDiagG%d", i, iMod, iDiagG+1).Data();
                if ( iDiagG == 0 )
                    hist[n] = new TH2F(n.c_str(), TString::Format("Module %d Diagnoal1_stripGroup %d;local strip x; local strip diagonal", iMod, iDiagG+1), 200, 0, 200, 200, 0, 200);
                if ( iDiagG == 1 )
                    hist[n] = new TH2F(n.c_str(), TString::Format("Module %d Diagnoal1_stripGroup %d;local strip diagonal;local strip y;", iMod, iDiagG+1), 200, 0, 200, 200, 0, 200);
            }
        }
    }

    fClusterProfile = new TF1("fClusterProfile", "gaus");
    fClusterProfile->SetParameters(1.0, 0, 1.4 * 3.2);
    // fClusterProfile->SetParameters( 1.0, 0, 1.4 * 3.2 );
    fClusterProfile->SetRange(-500, 500);
    fClusterProfile->SetNpx(1.e6);
    fBreitWigner = new TF1("fBreitWinger", lorentzianPeak, -500, 500, 3);
    fBreitWigner->SetParameters(-2.02447e+06, -1.89650e+00, -5.47860e-01);
    fBreitWigner->SetNpx(1.e6);

    // Just show that the bounds mapping is working
    for (int i = 0; i < 1000000; i++)
    {

        float rx = gRandom->Rndm() * 1400 - 700;
        float ry = gRandom->Rndm() * 1400 - 700;

        if (in_bounds(rx, ry))
            hist["hquad1"]->Fill(rx, ry);
        if (in_bounds_quad(rx, ry) > 0)
        {
            hist["hstgc"]->Fill(rx, ry);

            // ****************************************************************
            // fill noise
            // ****************************************************************
            if (gRandom->Rndm() < noise_prob)
            {
                float nv = gRandom->Rndm() * noise_level;
                ((TH2 *)hist["hGEN2"])->Fill(rx, ry, nv);
                if (in_bounds_diagonal(rx,ry) == 1)
                    {
                        ((TH2 *)hist["hGEN2_Diagnoal_1"])->Fill(rx, ry, nv);
                        TVector2 vec(rx,ry);
                        TVector2 vec_rot = vec.Rotate(TMath::Pi()/4);
                        ((TH2 *)hist["hGEN2_Diagnoal_rotate_1"])->Fill(vec_rot.X(),vec_rot.Y(),nv);
                    }
                if (in_bounds_diagonal(rx,ry) == 2)
                    {
                        ((TH2 *)hist["hGEN2_Diagnoal_2"])->Fill(rx, ry, nv);
                        TVector2 vec(rx,ry);
                        TVector2 vec_rot = vec.Rotate(-TMath::Pi()/4);
                        ((TH2 *)hist["hGEN2_Diagnoal_rotate_2"])->Fill(vec_rot.X(),vec_rot.Y(),nv);
                    }
                int sx = -1, sy = -1;
                map_to_strip(rx, ry, sx, sy);
                ((TH2 *)hist["hDIG0"])->Fill(sx, sy, nv);

                int q = map_to_local_strip(rx, ry, sx, sy);
                string hname = TString::Format("hDIG0L%d", q).Data();

                int xG = stripGroup(sx, sy);
                int yG = stripGroup(sy, sx); // not a typo, just reuse function with rotation

                if (hist.count(hname) > 0)
                {
                    ((TH2 *)hist[hname])->Fill(sx, sy, nv);

                    hname = TString::Format("hDIG0L%dhG%d", q, xG).Data();
                    ((TH2 *)hist[hname])->Fill(sx, sy, nv);

                    hname = TString::Format("hDIG0L%dvG%d", q, yG).Data();
                    ((TH2 *)hist[hname])->Fill(sx, sy, nv);
                }

                int d_p = diagonal_part(rx, ry);
                hname = TString::Format("hDIG0L%dDiagG%d", q, d_p).Data();
                // cout << "histo name is " << hname << "nCounts" << hist.count(hname) << endl;
                if ( d_p == 1)
                {
                    TVector2 vec;
                    if ( q == 1 ) vec.Set(rx,ry);
                    if ( q == 2 ) vec.Set(-rx,ry);
                    if ( q == 3 ) vec.Set(-(rx+pent_shift),-ry);
                    if ( q == 4 ) vec.Set(rx-pent_shift,-ry);
                    TVector2 vec_rot = vec.Rotate(pi/4);

                    if (hist.count(hname) > 0)
                    {
                        sy = map_to_local_strip_diagonal(rx,ry);
                        sx = vec_rot.X()/3.2;
                        ((TH2 *)hist[hname])->Fill(sx, sy, nv);
                    }
                }
                if ( d_p == 2)
                {
                    TVector2 vec;
                    if ( q == 1 ) vec.Set(rx,ry);
                    if ( q == 2 ) vec.Set(-rx,ry);
                    if ( q == 3 ) vec.Set(-(rx+pent_shift),-ry);
                    if ( q == 4 ) vec.Set(rx-pent_shift,-ry);
                    TVector2 vec_rot = vec.Rotate(-pi/4);

                    if (hist.count(hname) > 0)
                    {
                        sx = map_to_local_strip_diagonal(rx,ry);
                        sy = vec_rot.Y()/3.2;
                        ((TH2 *)hist[hname])->Fill(sx, sy, nv);
                    }
                }
            }
        }
    }

    int nClusters = 0;
    nChannel = 0;
    mEvtData.mEventId = 0;
    while (nClusters < n_clusters_to_gen)
    {
        
        float rx = gRandom->Rndm() * 1400 - 700;
        float ry = gRandom->Rndm() * 1400 - 700;
        // float rx = -296.280823;
        // float ry = -496.532684;
        int Edge_tag = Rejected_edge( rx, ry, 1 );
        if ( Edge_tag == -1) {cout << " cluster position Wrong !!!" << endl; continue;}
        if ( Edge_tag == Edge_Evt && test_edge == 1 ) continue;

        if (in_bounds_quad(rx, ry) > 0)
        {
            ((TH2 *)hist["hGEN0x"])->Fill(rx, ry, 1);
            ((TH2 *)hist["hGEN0y"])->Fill(rx, ry, 1);
            // diagnoal debug part
            if (debug == 1)
            {
                int sx, sy;
                int q = map_to_local_strip(rx, ry, sx, sy);
                int d_p = diagonal_part(rx, ry);
                TVector2 vec;
                if ( q == 1 ) vec.Set(rx,ry);
                if ( q == 2 ) vec.Set(-rx,ry);
                if ( q == 3 ) vec.Set(-(rx+pent_shift),-ry);
                if ( q == 4 ) vec.Set(rx-pent_shift,-ry);

                cout << "rx = " << rx << " ry = " << ry << endl;
                cout << " the diagnoal group is " << d_p << endl;
                TVector2 vec_rot;
                if ( d_p == 1)
                    { vec_rot = vec.Rotate(pi/4); cout << "Rotate : pi/4 " << endl; LOG_F(INFO, "Rotate : pi/4  " ); }
                if ( d_p == 2)
                    { vec_rot = vec.Rotate(-pi/4); cout << "Rotate : -pi/4 " << endl; LOG_F(INFO, "Rotate : -pi/4  " );}
                cout << "the rotated coordinate is （ " << vec_rot.X() << " , " << vec_rot.Y() << " )" << endl;
                LOG_F(INFO, "Rotated cluster at (%f, %f) ",  vec_rot.X(), vec_rot.Y() );

            }
            if (debug == 1) cout << "before Digitizing " << endl;
            fill_cluster_GEN(rx, ry);
            if (debug == 1) cout << "after Digitizing " << endl;
            hits_map->SetPoint(nClusters, rx, ry);
            nClusters++;
        }

    }
    mEvtData.nFireChannel = nChannel;
    mEvtTree->Fill();


    hist["hGEN2"]->Add(hist["hGEN1"]);

    //  hDIG2 = hDIG0 + hDIG1 (noise + signal)
    hist["hDIG2"]->Add(hist["hDIG0"]);
    hist["hDIG2"]->Add(hist["hDIG1x"]);
    hist["hDIG2"]->Add(hist["hDIG1y"]);

    // add for each module and strip group
    // also project them to make the 1D
    for (int iMod = 1; iMod <= 4; iMod++)
    {
        string n = TString::Format("L%d", iMod).Data();
        hist[("hDIG2" + n)]->Add(hist[("hDIG0" + n)]);
        hist[("hDIG2" + n)]->Add(hist[("hDIG1" + n)]);

        for (int iG = 0; iG < 3; iG++)
        {
            string n = TString::Format("L%dhG%d", iMod, iG).Data();
            hist[("hDIG2" + n)]->Add(hist[("hDIG0" + n)]);
            hist[("hDIG2" + n)]->Add(hist[("hDIG1" + n)]);
            // DIIG3 is 1D and DIG4 is saturated
            hist[("hDIG3" + n)] = ((TH2 *)hist[("hDIG2" + n)])->ProjectionY(("hDIG3" + n).c_str());

            hist[("hDIG4" + n)] = (TH1 *)hist[("hDIG3" + n)]->Clone(("hDIG4" + n).c_str());
            saturate(hist[("hDIG4" + n)]);

            n = TString::Format("L%dvG%d", iMod, iG).Data();
            hist[("hDIG2" + n)]->Add(hist[("hDIG0" + n)]);
            hist[("hDIG2" + n)]->Add(hist[("hDIG1" + n)]);
            hist[("hDIG3" + n)] = ((TH2 *)hist[("hDIG2" + n)])->ProjectionX(("hDIG3" + n).c_str());

            hist[("hDIG4" + n)] = (TH1 *)hist[("hDIG3" + n)]->Clone(("hDIG4" + n).c_str());
            saturate(hist[("hDIG4" + n)]);
        }
        for( int iG = 0; iG < 2; iG++)
        {
            string n = TString::Format("L%dDiagG%d", iMod, iG+1).Data();
            hist[ ( "hDIG2" + n ) ]->Add(hist[("hDIG0" + n)]);
            hist[ ( "hDIG2" + n ) ]->Add(hist[("hDIG1" + n)]);
            if ( iG == 0)
                hist[ ( "hDIG3" + n ) ] = ((TH2 *)hist[("hDIG2" + n)])->ProjectionY(("hDIG3" + n).c_str());
            if ( iG == 1)
                hist[ ( "hDIG3" + n ) ] = ((TH2 *)hist[("hDIG2" + n)])->ProjectionX                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  (("hDIG3" + n).c_str());

        }
    }


    hist["hDIG1_diagonaly"]->Print("base");
    // TFile * fout = new TFile( "output.root", "RECREATE" );
    TFile *fout = new TFile(argv[1], "RECREATE");
    fout->cd();

    for (auto nh : hist)
    {
        if (nullptr != nh.second)
            nh.second->Write();
    }
    for (auto nh : graph)
    {
        if (nullptr != nh.second)
            nh.second->Write();
    }

    hist["hDIG1_diagonaly"]->Write();
    fClusterProfile->Write();
    fBreitWigner->Write();
    hits_map->SetNameTitle("Hits_map", "Hits_map");
    hits_map->Write();
    mEvtTree->Write();
    hDelta_BCID->Write();
    Diag_hits1_map->SetNameTitle("Diag_hits1_map", "Diag_hits1_map");
    Diag_hits1_map->Write();
    Diag_hits2_map->SetNameTitle("Diag_hits2_map", "Diag_hits2_map");
    Diag_hits2_map->Write();

    fout->Write();
    fout->Close();
}