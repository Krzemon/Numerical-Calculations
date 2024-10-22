#include <iostream>
#include <cmath>
#include <ctime>
#include "u_xt_exact.h"

using namespace std;

int main()
{
    double x_start = 0.0;
    double t_start = 10.0*pow(10,-9);//22.0*pow(10,-9);
    int npaths = 10000;

    double L= 0.25*pow(10,-6), C1=100*pow(10,-12), R=12.5, G=0.5*pow(10,-3), Rg = 75, Rl = 12.5, sigma = 0.75*pow(10,-9), v =pow(10,9), l=2.0, t_0=7.5*pow(10,-9);

    double c = 1/sqrt(L*C1);
    double mi=G/C1;
    double lambda = 0.5*(R/L - G/C1);

    double R0 = sqrt(L/C1);
    double zeta = R0/(R0+Rg);
    double Gamma_g = (Rg - R0)/(Rg + R0);
    double Gamma_l = (Rl - R0)/(Rl + R0);

    double fxt, bxt;
    int liczba = 500;
    double uxt[liczba];

    double x, x_write[liczba]={0};
    double u_dokladne[liczba] ={0};


    double lenght = l/(liczba-1);

    srand(time(NULL));
    FILE *f;
    f = fopen("../data/lab9_t50_n5.dat", "w");

    for(int j=0; j<liczba; j++){
        x_start = lenght*j;
        for (int i = 1; i < 3; i++) {
            double suma = 0;
            for (int n = 0; n < npaths; n++) {
                double t = t_start;
                x = x_start;
                double eta = 1;
                int sign = pow((-1),i);
                while (t > 0) {
                    double U1 = rand() / (RAND_MAX + 1.0);
                    double s = (-1) * log(U1) / (lambda + mi);
                    //cout << s << "\t" << c << "\t" << s*c << endl;
                    if (sign == -1) {
                            //cout << n << endl;
                        if ((x - c*s) > 0) {
                            eta = eta * lambda / (lambda + mi);
                            //cout << mi << endl;
                        } else {
                            s = x / c;
                            double Vg = sin(2*M_PI*pow(10,9)*(t-s))*exp(-(t-s-t_0)*(t-s-t_0)/(2*sigma*sigma));
                            suma = suma + eta * zeta * Vg;
                            eta = eta * Gamma_g;
                            //cout << n << endl;
                        }
                        x = x - c*s;
                        t = t - s;
                    } else if (sign == 1) {
                        if ((x + c*s) < l) {
                            eta = eta * lambda / (lambda + mi);
                        } else {
                            s = (l - x) / c;
                            eta = eta * Gamma_l;
                            //cout << n << endl;
                        }
                        x = x + c*s;
                        t = t - s;
                    }
                    sign = -sign;
                    //cout << x  << "\t" << s << "\t" << c*s << endl;
                }
            }

            if (i == 1) {
                fxt = suma / npaths;
                //cout << fxt << endl;
            }
            if (i == 2) {
                bxt = suma / npaths;
                //cout << bxt[n] << endl;
            }
            x_write[j]=lenght*j;
            uxt[j] = bxt + fxt;

            u_dokladne[j] = u_xt_exact ( x_write[j] , t_start , t_0 , pow(10,9) , sigma, R , G , L , C1 , Rg , Rl, l , 1000 , 100 );
        }
    }

    for(int n=0; n<500;n++){
       // cout << bxt[n] << "\t" << x[n] << endl;
        fprintf(f, "%.16f \t %.16f \t %.16f \n", x_write[n], uxt[n], u_dokladne[n]);
    }
    cout << u_xt_exact ( x_write[0] , t_start , t_0 , pow(10,9) , sigma, R , G , L , C1 , Rg , Rl, l , 1000 , 100 ) << endl;

   return 0;
}
