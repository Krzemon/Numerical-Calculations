#include <iostream>
#include "photon_diffusion.h"

using namespace std;

int main()
{
    int nlayers=3,
    nx=200,
    ny=200;
    double x_max=0.2,
    x_source=0.1,
    dx_source=0.0,
    x_detect=0.15,
    dx_detect=0.01,
    rx=0.8,
    ry=0.6;
    int N=200000;

    PHOTON_DIFFUSION_2D ob ;
    ob.xmax =x_max;
    ob.x_source =x_source;
    ob.dx_source =dx_source;
    ob.x_detect =x_detect;
    ob.dx_detect =dx_detect;
    ob.nx =nx;
    ob.ny =ny;
    ob.rx0 =rx;
    ob.ry0 =ry;
    ob.nlayers = nlayers;
    ob.layers_data [1][0]=1.; // absorption
    ob.layers_data [1][1]=10.; // scattering
    ob.layers_data [1][2]=0.02; // width
    ob.layers_data [1][3]=0.75; // g_anizo
    ob.layers_data [1][4]=1.; // n_refraction
    ob.layers_data [2][0]=10.; // absorption
    ob.layers_data [2][1]=210.; // scattering
    ob.layers_data [2][2]=0.02; // width
    ob.layers_data [2][3]=0.075; // g_anizo
    ob.layers_data [2][4]=1.5; // n_refraction
    ob.layers_data [3][0]=10.; // absorption
    ob.layers_data [3][1]=90.; // scattering
    ob.layers_data [3][2]=0.02; // width
    ob.layers_data [3][3]=0.95; // g_anizo
    ob.layers_data [3][4]=1.0;
    ob.init ();
    ob.write_all_paths =0;
    ob.write_source_detection_paths =0;
    for ( int k =0; k <N ;k++){
        ob.single_path();
    }
    double suma_abs=0;
    double suma_ref =0;
    double suma_tr=0;

    FILE *fp;
    fp=fopen("../data/absorpcja.dat","w");
    fprintf(fp,"\n");

    FILE *fr;
    fr=fopen("../data/reflektancja.dat","w");
    fprintf(fr,"\n");

    FILE *ft;
    ft=fopen("../data/transmitancja.dat","w");
    fprintf(ft,"\n");

    for(int i=0; i<=nx; i++)
    {
        for(int j=0; j<=ny; j++)
        {
            suma_abs = suma_abs + ob.absorption[i][j];
            fprintf(fp,"%15.5E \t",ob.absorption[i][j]);

        }
        fprintf(fp,"\n");
        suma_ref = suma_ref +ob.reflectance[i];
        suma_tr = suma_tr +ob.transmittance[i];
        fprintf(fr,"%15.5E \n",ob.reflectance[i]);
        fprintf(ft,"%15.5E \n",ob.transmittance[i]);
    }
    cout << "W sumie " << (suma_abs + suma_ref + suma_tr)/(double)N << " absorpcja " << (suma_abs)/(double)N << " reflektancja " << (suma_ref )/(double)N << " transmitancja " << (suma_tr)/(double)N << endl;


    fclose(fp);
    fclose(fr);
    fclose(ft);



    return 0;
}
