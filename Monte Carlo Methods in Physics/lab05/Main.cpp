#include <iostream>
#include "Metoda_warstw.h"

double function_1(double x){
    return 1 + tanh(x);
}

double function_2(double x){
    return 1 / (1 + std::pow(x, 2));
}

double function_3(double x){
    return std::pow(std::cos(M_PI* x), 10);
}


int main(){
    std::cout << std::endl;

    int K[] = {2, 3, 4, 5};
    std::cout << "Całka oznaczona z funkcji f(x) = 1 + tanh(x) w granicach [-3,3]\n";

    integral(function_1, K, -3, 3); 

    // std::cout << "Całka oznaczona z funkcji f(x) = 1 / (1 + std::pow(x, 2)) w granicach [0,10]\n";
    // for (int i = 0; i < 4; ++i){
    //     integral(function_2, K, 0, 10);
    // }

    // std::cout << "Całka oznaczona z funkcji f(x) = std::pow(std::cos(M_PI* x), 10) w granicach [0,1]\n";
    // for (int i = 0; i < 4; ++i){
    //     integral(function_3, K, 0, 1);
    // }

    return 0;
}