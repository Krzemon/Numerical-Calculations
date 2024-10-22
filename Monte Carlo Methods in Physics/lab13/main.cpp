using namespace std;
#include "dsmc_2d.h"


int main (){
    DSMC_2D ob ;
    ob.read ("../i.dat"); // wczytujemy dane z pliku wejściowego
    // ob.read ("../pos_vel_start.dat");
    
    ob.init(); // automatyczna inicjalizacja położeń i prędkości
    // ob.write_position_velocity("../rv_left.dat");
    
    ob.nthreads = 16; // obliczenia na jednym rdzeniu
    ob.icol = 1; // cząstki zderzają się

    ob.evolution(0.0 ,2000); // wykonujemy 2 tysięcy kroków ( tmax - nieznany 

    // ob.hist_velocity_all("../data/hist2_test.dat" ,5.0 ,50); // zapis histogramu prędkosci do pliku
    ob.write_position_velocity("../data/rv.dat"); // zapis położeń i prędkości końcowych do pliku

return 0;
}