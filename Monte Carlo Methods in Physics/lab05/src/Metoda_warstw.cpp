#include "Metoda_warstw.h"


double uniform() {
    return static_cast<double>(std::rand()) / RAND_MAX;
}


void integral(std::function<double(double x)> fun, int K[4], double a, double b){
    double sum = 0, sum_var = 0, temp = 0;
    int kk = 0;
    int N = static_cast<int>(std::pow(10, K[3]));
    
    for (int i = 0;  i < N; ++i){
        // losowana zmienna z przedzialu całkowania
        double U1 = uniform();
        double x = a + (b-a)*U1; 
        temp  = fun(x) * (b-a);
        sum += temp;
        sum_var  += std::pow(sum, 2);

        if( i == static_cast<int>(std::pow(10, K[kk++]))){
            std::cout << "dla liczby iteracji " << std::pow(10, i) << " całka wynosi:" << sum / i << std::endl;
            std::cout << "wariancja wynosi: " << sum_var / i << std::endl;
            std::cout << "sigma wynosi: " << (sum_var - sum) / i << std::endl;

        }
    }
}
