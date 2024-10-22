#include "Includes.h"

double randNormal() {
    static std::default_random_engine generator;
    static std::normal_distribution<double> distribution(0.0, 1.0);
    return distribution(generator);
}