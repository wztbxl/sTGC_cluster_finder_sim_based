

#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TGraph.h"

#include <iostream>
#include <map>
#include <ctime>

#define LOGURU_IMPLEMENTATION 1
#include "loguru.h"


const float pent_base = 537.0; // in mm
const float pent_nib = 179.0; // in mm
const float pent_shift = 101.6; // in mm
const int n_clusters_to_gen = 10; // Number of clusters in this "event"
const float noise_prob = 0.2; // 0 -1, 1 is max noise
// const float noise_level = 1; // in ADC units (pre integration)
const float noise_level = 1.e-10; // in ADC units (pre integration)
const float strip_pitch = 3.2; // digitize strip pitch
const size_t sat_above = 1024; // ADC
const float cluster_max_adc = 120; // ADC (pre integration)

TF1 * fClusterProfile = nullptr;
TF1 * fBreitWigner = nullptr;

using namespace std;
map<string, TH1*> hist;
map<string, TGraph*> graph;
TGraph* hits_map = new TGraph();

Double_t lorentzianPeak(Double_t *x, Double_t *par) 
{
  return (0.5*par[0]*par[1]/TMath::Pi()) / TMath::Max(1.e-10,(x[0]-par[2])*(x[0]-par[2])+ .25*par[1]*par[1]);
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

int map_to_strip( float x, float y, int &sx, int &sy ){
    int q = in_bounds_quad( x, y );
    if ( q == 1 || q == 2 ){
        sx = x / strip_pitch;
        sy = y / strip_pitch;
    } else if ( q == 3 ){
        sx = (x ) / strip_pitch;
        sy = y / strip_pitch;
    } else if ( q == 4 ){
        sx = (x ) / strip_pitch;
        sy = y / strip_pitch;
    }
    return q;
}

int map_to_local_strip( float x, float y, int &sx, int &sy ){
    int q = in_bounds_quad( x, y );
    if ( q == 1 || q == 2 ){
        sx = abs(x / strip_pitch);
        sy = abs(y / strip_pitch);
    } else if ( q == 3 ){
        sx = abs((x + pent_shift ) / strip_pitch);
        sy = abs(y / strip_pitch);
    } else if ( q == 4 ){
        sx = abs((x - pent_shift) / strip_pitch);
        sy = abs(y / strip_pitch);
    }
    return q;
}

int stripGroup( int sx, int sy ){
    if ( (sx < 55 && sy < 150) || sy >= 150 ) return 0;
    if ( (sx < 110 && sy < 95) || sy >= 95 ) return 1;
    return 2;
}

void fill_cluster_GEN( float x, float y ){
    TH2 * hGENx = (TH2*)hist[ "hGEN0x" ];
    TH2 * hGENy = (TH2*)hist[ "hGEN0y" ];
    TH2 * hDIGx = (TH2*)hist[ "hDIG1x" ];
    TH2 * hDIGy = (TH2*)hist[ "hDIG1y" ];
    TH2 * hDIGxRecorder = (TH2*)hist[ "hDIG1xRecorder" ];
    TH2 * hDIGyRecorder = (TH2*)hist[ "hDIG1yRecorder" ];
    // int bx = hGENx->GetXaxis()->FindBin( x );
    // int by = hGENx->GetYaxis()->FindBin( y );

    float maxADC = 0;
    while ( maxADC < noise_level*7 )
    {
        maxADC = gRandom->Rndm() * cluster_max_adc;
    }
    // maxADC = 100;
    
    // float maxADC = gRandom->Rndm() * cluster_max_adc;
    LOG_F( INFO, "Cluster at (%f, %f) w maxADC = %f", x, y, maxADC );

    
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
    map_to_strip( x, y, bx, by );

    int wSpan = 10;
    for ( int i = bx - wSpan; i < bx + wSpan; i++ ) // sample in X direction
    {
        float bcx = i * 3.2 + 1.6;
        float bLowEdge = i * 3.2;
        float bHighEdge = (i+1) * 3.2;
        // float bcx = hGENx->GetXaxis()->GetBinCenter( i );
        int q = in_bounds_quad( bcx, y );
        if ( q <= 0 ) continue;

        float r = bcx - x;
        float rL = bLowEdge - x;
        float rH = bHighEdge - x;
        float v = fClusterProfile->Eval( r ) * maxADC;
        // float v = fClusterProfile->Integral( rL, rH ) * maxADC;
        // cout << i <<" th bin " << endl; 
        // cout << "bin center = " << bcx << " r = " << r << " v = " << v << endl;

        // hGENx->Fill( bcx, y, v );

        int sx = -1, sy = -1;
        map_to_strip( bcx, y, sx, sy );
        hDIGx->Fill( sx, sy, v );
        hDIGxRecorder->Fill( sx, sy, r );

        map_to_local_strip( bcx, y, sx, sy);
        // cout << "sx = " << sx << " sy = " << sy << endl;
        string hname = TString::Format( "hDIG1L%d", q ).Data();
        int xG = stripGroup( sy, sx ); // not a typo, just reuse function with rotation
        if ( hist.count( hname ) > 0 )
        {
            ((TH2*)hist[ hname ]) -> Fill( sx, sy, v );

            hname = TString::Format( "hDIG1L%dvG%d", q, xG ).Data();
            ((TH2*)hist[ hname ]) -> Fill( sx, sy, v );
            if ( i == bx )
            {
                string gname = TString::Format( "gDIG1L%dvG%d", q, xG ).Data();
                ((TGraph*)graph[ gname ])->SetPoint( graph[ gname ]->GetN(),x,y );
            }
        }
    }
    for ( int j = by - wSpan; j < by + wSpan; j++ )
    {
        float bcy = j * 3.2 + 1.6;
        float bLowEdge = j * 3.2;
        float bHighEdge = (j+1) * 3.2;
        // float bcy = hGENy->GetYaxis()->GetBinCenter( j );
        int q = in_bounds_quad( x, bcy );
        if ( q <= 0 ) continue;

        float r = bcy - y;
        float rL = bLowEdge - y;
        float rH = bHighEdge - y;
        float v = fClusterProfile->Eval( r ) * maxADC;
        // float v = fClusterProfile->Integral( rL, rH ) * maxADC;

        // hGENy->Fill( x, bcy, v );

        int sx = -1, sy = -1;
        map_to_strip( x, bcy, sx, sy );
        hDIGy->Fill( sx, sy, v );
        hDIGyRecorder->Fill( sx, sy, r ); 

        map_to_local_strip( x, bcy, sx, sy);
        string hname = TString::Format( "hDIG1L%d", q ).Data();
        int yG = stripGroup( sx, sy );

        if ( hist.count( hname ) > 0 )
        {
            ((TH2*)hist[ hname ]) -> Fill( sx, sy, v );

            hname = TString::Format( "hDIG1L%dhG%d", q, yG ).Data();
            ((TH2*)hist[ hname ]) -> Fill( sx, sy, v );
            if ( j == by )
            {
                string gname = TString::Format( "gDIG1L%dhG%d", q, yG ).Data();
                ((TGraph*)graph[ gname ])->SetPoint( graph[ gname ]->GetN(),x,y );
            }
        }
    }
    
        
}

void saturate( TH1 * h ){
    for ( int i = 1; i < h->GetNbinsX(); i++ ){
        if ( h->GetBinContent( i ) < sat_above ) continue;

        h->SetBinContent( i, (int)(sat_above - 1) );
    }
}

int main( int argc, char** argv ){

    loguru::init(argc, argv);
    // loguru::g_stderr_verbosity = 8;
    //add argv as input
    if (argc != 3) cout << "2 arguements is need , 1 is outfile name, 2 is logname" << endl;
    loguru::add_file(argv[2], loguru::Truncate, loguru::Verbosity_MAX);
    // loguru::add_file("everything.log", loguru::Truncate, loguru::Verbosity_MAX);

    LOG_SCOPE_FUNCTION( INFO );

    gRandom = new TRandom3();
    gRandom->SetSeed((int) time(NULL) );

    

    hist["hquad1"] = new TH2F( "hquad1", ";x;y", 600, 0, 600, 600, 0, 600 );
    hist["hstgc"] = new TH2F( "hstgc", ";x;y", 1400, -700, 700, 1200, -600, 600 );

    hist["hGEN0x"] = new TH2F( "hGEN0x", ";x;y", 1400, -700, 700, 1200, -600, 600 );
    hist["hGEN0y"] = new TH2F( "hGEN0y", ";x;y", 1400, -700, 700, 1200, -600, 600 );
    hist["hGEN1"] = new TH2F( "hGEN1", ";x;y", 1400, -700, 700, 1200, -600, 600 );
    hist["hGEN2"] = new TH2F( "hGEN2", ";x;y", 1400, -700, 700, 1200, -600, 600 );
    hist["hDIG0"] = new TH2F( "hDIG0", ";strip x; strip y", 400, -200, 200, 400, -200, 200 );
    hist["hDIG1x"] = new TH2F( "hDIG1x", ";strip x; strip y", 400, -200, 200, 400, -200, 200 );
    hist["hDIG1y"] = new TH2F( "hDIG1y", ";strip x; strip y", 400, -200, 200, 400, -200, 200 );
    hist["hDIG1xRecorder"] = new TH2F( "hDIG1xRecorder", ";strip x; strip y", 400, -200, 200, 400, -200, 200 );
    hist["hDIG1yRecorder"] = new TH2F( "hDIG1yRecorder", ";strip x; strip y", 400, -200, 200, 400, -200, 200 );
    hist["hDIG2"] = new TH2F( "hDIG2", ";strip x; strip y", 400, -200, 200, 400, -200, 200 );
    for ( int i = 0; i < 3; i++ ){
        for ( int iMod = 1; iMod < 5; iMod++ ){
            string n = TString::Format( "hDIG%dL%d", i, iMod ).Data();
            LOG_F( INFO, "Making : %s", n.c_str() );
            hist[n] = new TH2F( n.c_str(), TString::Format("Module %d;local strip x; local strip y", iMod), 200, 0, 200, 200, 0, 200 );
            for ( int ixG = 0; ixG < 3; ixG++ ){
                string n = TString::Format( "hDIG%dL%dhG%d", i, iMod, ixG ).Data();
                hist[n] = new TH2F( n.c_str(), TString::Format("Module %d stripGroup %d;local strip x; local strip y", iMod, ixG), 200, 0, 200, 200, 0, 200 );
                if ( i == 1)
                {
                    string gname = TString::Format( "gDIG%dL%dhG%d", i, iMod, ixG ).Data();
                    graph[gname] = new TGraph();
                    graph[gname]->SetNameTitle(gname.data(),gname.data());
                    
                }

            }
            for ( int iyG = 0; iyG < 3; iyG++ ){
                string n = TString::Format( "hDIG%dL%dvG%d", i, iMod, iyG ).Data();
                hist[n] = new TH2F( n.c_str(), TString::Format("Module %d yGroup %d;local strip x; local strip y", iMod, iyG), 200, 0, 200, 200, 0, 200 );
                if ( i == 1 )
                {
                    string gname = TString::Format( "gDIG%dL%dvG%d", i, iMod, iyG ).Data();
                    graph[gname] = new TGraph();
                    graph[gname]->SetNameTitle(gname.data(),gname.data());
                    
                }
            }
        }
    }

    fClusterProfile = new TF1( "fClusterProfile", "gaus" );
    fClusterProfile->SetParameters( 1.0, 0, 1.4 * 3.2 );
    // fClusterProfile->SetParameters( 1.0, 0, 1.4 * 3.2 );
    fClusterProfile->SetRange( -500, 500 );
    fClusterProfile->SetNpx(1.e6);
    fBreitWigner = new TF1("fBreitWinger",lorentzianPeak,-500,500,3);
    fBreitWigner->SetParameters(-2.02447e+06,-1.89650e+00,-5.47860e-01);
    fBreitWigner->SetNpx(1.e6);



    // Just show that the bounds mapping is working
    for ( int i = 0; i <1000000; i++ ){

        float rx = gRandom->Rndm() * 1400 - 700;
        float ry = gRandom->Rndm() * 1400 - 700;

        if ( in_bounds( rx, ry ) )
            hist["hquad1"]->Fill( rx, ry );
        if ( in_bounds_quad( rx, ry ) > 0 ){
            hist["hstgc"]->Fill( rx, ry );

            // ****************************************************************
            // fill noise
            // ****************************************************************
            if ( gRandom->Rndm() < noise_prob ){
                float nv = gRandom->Rndm() * noise_level;
                ((TH2*)hist["hGEN2"])->Fill( rx, ry,  nv );
                int sx = -1, sy = -1;
                map_to_strip( rx, ry, sx, sy );
                ((TH2*)hist["hDIG0"])->Fill( sx, sy, nv );

                int q = map_to_local_strip( rx, ry, sx, sy );
                string hname = TString::Format( "hDIG0L%d", q ).Data();

                int xG = stripGroup( sx, sy );
                int yG = stripGroup( sy, sx ); // not a typo, just reuse function with rotation

                if ( hist.count( hname ) > 0 ){
                    ((TH2*)hist[ hname ]) -> Fill( sx, sy, nv );

                    hname = TString::Format( "hDIG0L%dhG%d", q, xG ).Data();
                    ((TH2*)hist[ hname ]) -> Fill( sx, sy, nv );

                    hname = TString::Format( "hDIG0L%dvG%d", q, yG ).Data();
                    ((TH2*)hist[ hname ]) -> Fill( sx, sy, nv );
                }
            }
        }

    }


    int nClusters = 0;
    while( nClusters < n_clusters_to_gen ){
        float rx = gRandom->Rndm() * 1400 - 700;
        float ry = gRandom->Rndm() * 1400 - 700;
        // float rx = 90.799255;
        // float ry = 194.756744;
        if ( in_bounds_quad( rx, ry ) > 0 ){
            ((TH2*)hist["hGEN0x"])->Fill( rx, ry, 1 );
            ((TH2*)hist["hGEN0y"])->Fill( rx, ry, 1 );
            fill_cluster_GEN( rx, ry );
            hits_map->SetPoint(nClusters,rx,ry);
            nClusters++;
        }
    }

    hist["hGEN2"]->Add( hist["hGEN1"] );

    //  hDIG2 = hDIG0 + hDIG1 (noise + signal)
    hist["hDIG2"]->Add( hist[ "hDIG0" ] );
    hist["hDIG2"]->Add( hist[ "hDIG1x" ] );
    hist["hDIG2"]->Add( hist[ "hDIG1y" ] );

    // add for each module and strip group
    // also project them to make the 1D
    for ( int iMod = 1; iMod <= 4; iMod ++ ){
        string n = TString::Format( "L%d", iMod ).Data();
        hist[ ("hDIG2" + n) ] -> Add( hist[("hDIG0" + n)] );
        hist[ ("hDIG2" + n) ] -> Add( hist[("hDIG1" + n)] );

        for ( int iG = 0; iG < 3; iG++ ){
            string n = TString::Format( "L%dhG%d", iMod, iG ).Data();
            hist[ ("hDIG2" + n) ] -> Add( hist[("hDIG0" + n)] );
            hist[ ("hDIG2" + n) ] -> Add( hist[("hDIG1" + n)] );
            // DIIG3 is 1D and DIG4 is saturated
            hist[ ("hDIG3" + n) ] = ((TH2*)hist[ ("hDIG2" + n) ])->ProjectionY( ("hDIG3" + n).c_str() );

            hist[ ("hDIG4" + n) ] = (TH1*)hist[ ("hDIG3" + n) ]->Clone( ("hDIG4" + n).c_str() );
            saturate( hist[ ("hDIG4" + n) ] );

            n = TString::Format( "L%dvG%d", iMod, iG ).Data();
            hist[ ("hDIG2" + n) ] -> Add( hist[("hDIG0" + n)] );
            hist[ ("hDIG2" + n) ] -> Add( hist[("hDIG1" + n)] );
            hist[ ("hDIG3" + n) ] = ((TH2*)hist[ ("hDIG2" + n) ])->ProjectionX( ("hDIG3" + n).c_str() );

            hist[ ("hDIG4" + n) ] = (TH1*)hist[ ("hDIG3" + n) ]->Clone( ("hDIG4" + n).c_str() );
            saturate( hist[ ("hDIG4" + n) ] );
        }
    }


    // TFile * fout = new TFile( "output.root", "RECREATE" );
    TFile * fout = new TFile( argv[1], "RECREATE" );
    fout->cd();

    for ( auto nh : hist ){
        if ( nullptr != nh.second )
            nh.second->Write();
    }
    for ( auto nh : graph ){
        if ( nullptr != nh.second )
            nh.second->Write();
    }



    fClusterProfile->Write();
    fBreitWigner->Write();
    hits_map->SetNameTitle("Hits_map","Hits_map");
    hits_map->Write();

    fout->Write();
    fout->Close();

}