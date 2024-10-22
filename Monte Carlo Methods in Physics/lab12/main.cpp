using namespace std;
#include "dsmc_2d.h"

int main (){
    DSMC_2D ob ;
    ob.read ("../i.dat"); // wczytujemy dane z pliku wejściowego
    ob.init(); // automatyczna inicjalizacja położeń i prędkości
    // ob.write_position_velocity("../data/rv.dat"); // zapis ustawień początkowych
    ob.nthreads = 20; // obliczenia na 20 rdzeniach :)
    ob.icol = 1; // cząstki zderzają się
    // ob.evolution(0.0 ,200); // wykonujemy 20 tysięcy kroków ( tmax - nieznany )
    // ob.hist_velocity_all("../data/hist2.dat" ,5.0 ,50); // zapis histogramu prędkosci do pliku
    // ob.write_position_velocity("../data/rv.dat"); // zapis położeń i prędkości końcowych do pliku


    ob.evolution(0.0 ,1000); // 5000 iteracji
    ob.hist_velocity_all("../zad_4_data_05/hist_1k.dat" ,5.0 ,50);
    // ob.write_position_velocity("../data/rv_2k.dat");
    ob.write_nptv("../zad_4_data_05/nptv_1k.dat", 1);
    
return 0;
}