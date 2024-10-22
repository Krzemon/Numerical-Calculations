#include "geometry.h"

const int DISPLAY_PRECISION = 7;

int main() {
    const double RA = 1.0;
    const double RB = std::sqrt(2.0) * RA;
    double x_A = RA + RB;
    
    // 1. Obliczanie wspólnej powierzchnii kół za pomocą metody Monte Carlo
    std::vector<double> X_RA(10000), Y_RA(10000);
    std::vector<double> X_RB(10000), Y_RB(10000);

    sferycznie_konturowany(X_RA, Y_RA, RA);
    sferycznie_konturowany(X_RB, Y_RB, RB);

    circle_Norm(X_RA, Y_RA);
    circle_Norm(X_RB, Y_RB);

    translation(X_RA, Y_RA, x_A, 0.0);

    double overlapArea = calculateOverlapArea(X_RA, Y_RA, RA, RB, x_A, 10000); // Obliczanie wspólnej powierzchnii
    std::cout << "Overlap Area: " << overlapArea << std::endl;

    // Zapisywanie punktów do plików
    drawPoints(X_RA, Y_RA, "../data/points_A.dat");
    drawPoints(X_RB, Y_RB, "../data/points_B.dat");

    // ------------------------------------------------

    const int newSize = 1000000; // 10^6

    X_RA.resize(newSize);
    Y_RA.resize(newSize);
    X_RB.resize(newSize);
    Y_RB.resize(newSize);

    // obliczyć wartości S̄_{α,β} i σ_{S̄_{α,β}} dla n =10^k , k = 2, 3, 4, 5, 6
    std::ofstream file_A("../data/A.dat", std::ios::out);
    std::ofstream file_B("../data/B.dat", std::ios::out);
    std::ofstream file_C("../data/C.dat", std::ios::out);
    std::ofstream file_D("../data/D.dat", std::ios::out);

    double mean, secondMoment, variance, standardDeviation;
    // ------------------------------------------------
    // α = A, xA = RB + 0.5 · RA
    x_A = RB + 0.5 * RA;

    sferycznie_konturowany(X_RA, Y_RA, RA);
    circle_Norm(X_RA, Y_RA);
    translation(X_RA, Y_RA, x_A, 0.0);

    
    for (int k = 2; k <= 6; ++k) {
            overlapArea = calculateOverlapArea(X_RA, Y_RA, RA, RB, x_A, std::pow(10, k));
            mean = overlapArea;
            secondMoment = overlapArea * calculateCircleArea(RA);
            variance = calculateVariance(mean, secondMoment, std::pow(10, k));
            standardDeviation = calculateStandardDeviation(mean, secondMoment, std::pow(10, k));

            file_A << std::fixed << std::setprecision(DISPLAY_PRECISION);
            file_A << k << '\t' << overlapArea << '\t' << standardDeviation << std::endl;
    }

    // ------------------------------------------------
    // α = A, xA = 0
    x_A = 0;

    sferycznie_konturowany(X_RA, Y_RA, RA);
    circle_Norm(X_RA, Y_RA);
    translation(X_RA, Y_RA, x_A, 0.0);

    
    for (int k = 2; k <= 6; ++k) {
            overlapArea = calculateOverlapArea(X_RA, Y_RA, RA, RB, x_A, std::pow(10, k));
            mean = overlapArea;
            secondMoment = overlapArea * calculateCircleArea(RA);
            variance = calculateVariance(mean, secondMoment, std::pow(10, k));
            standardDeviation = calculateStandardDeviation(mean, secondMoment, std::pow(10, k));

            file_B << std::fixed << std::setprecision(DISPLAY_PRECISION);
            file_B << k << '\t' << overlapArea << '\t' << standardDeviation << std::endl;
    }

    // ------------------------------------------------
    // α = B, xA = RB + 0.5 · RA
    x_A = RB + 0.5 * RA;

    sferycznie_konturowany(X_RB, Y_RB, RB);
    circle_Norm(X_RB, Y_RB);
    translation(X_RA, Y_RA, -x_A, 0.0);

    
    for (int k = 2; k <= 6; ++k) {
            overlapArea = calculateOverlapArea(X_RB, Y_RB, RA, RB, -x_A, std::pow(10, k));
            mean = overlapArea;
            secondMoment = overlapArea * calculateCircleArea(RB);
            variance = calculateVariance(mean, secondMoment, std::pow(10, k));
            standardDeviation = calculateStandardDeviation(mean, secondMoment, std::pow(10, k));

            file_C << std::fixed << std::setprecision(DISPLAY_PRECISION);
            file_C << k << '\t' << overlapArea << '\t' << standardDeviation << std::endl;
    }

    // ------------------------------------------------
    // α = B, xA = 0
    x_A = 0;

    sferycznie_konturowany(X_RB, Y_RB, RB);
    circle_Norm(X_RB, Y_RB);
    translation(X_RA, Y_RA, -x_A, 0.0);

    
    for (int k = 2; k <= 6; ++k) {
            overlapArea = calculateOverlapArea(X_RB, Y_RB, RA, RB, -x_A, std::pow(10, k));
            mean = overlapArea;
            secondMoment = overlapArea * calculateCircleArea(RB);
            variance = calculateVariance(mean, secondMoment, std::pow(10, k));
            standardDeviation = calculateStandardDeviation(mean, secondMoment, std::pow(10, k));

            file_D << std::fixed << std::setprecision(DISPLAY_PRECISION);
            file_D << k << '\t' << overlapArea << '\t' << standardDeviation << std::endl;
    }

    return 0;
}
