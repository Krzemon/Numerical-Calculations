#include <iostream>
#include "Brenner.h"

int main() {

    int n = 60;
    std::vector<double> atoms(n);

    // zadanie 1 
    double pot = Potental(atoms, n);


    return 0;
}