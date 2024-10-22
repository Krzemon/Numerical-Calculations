#include "rozklady_2D.h"

void Box_Muller(std::vector<double> &X, std::vector<double> &Y){
    for(int i = 0; i < n; ++i){
        double U1 = uniform();
        double U2 = uniform();
        X[i] = std::sqrt(-2 * std::log(1 - U1)) * std::cos(2 * M_PI * U2);
        Y[i] = std::sqrt(-2 * std::log(1 - U1)) * std::sin(2 * M_PI * U2);
    }
}

void Circle_Norm(std::vector<double> &X, std::vector<double> &Y){
    for(int i = 0; i < n; ++i){
        double U1 = uniform();

        X[i] *= std::sqrt(U1);
        Y[i] *= std::sqrt(U1);
    }
}



void Variance(std::vector<double> &X, std::vector<double> &Y, double cov[2][2]){
    double x;
    double EX = 0.0;
    double EY = 0.0;
    double EXsquare = 0.0;
    double EYsquare = 0.0;
    double EXY = 0.0;

    for (int i = 0; i < n ; ++i) {
        EX += X[i];
        EY += Y[i];
        EXsquare += X[i] * X[i];
        EYsquare += Y[i] * Y[i];
        EXY += X[i] * Y[i];
    }
        double EX_mean = EX / n;
        double EY_mean = EY / n;

        cov[0][0] = (EXsquare - EX_mean * EX_mean) / n;
        cov[1][1] = (EYsquare - EY_mean * EY_mean) / n;
        cov[0][1] = (EXY - EX_mean * EY_mean) / n;
        cov[1][0] = cov[0][1];
}

void Print_Cov(double cov[2][2]){
    std::cout << "Covariance matrix: " << std::endl;
    std::cout << cov[0][0] << " " << cov[0][1] << std::endl;
    std::cout << cov[1][0] << " " << cov[1][1] << std::endl;

    double rho = cov[0][1] / std::sqrt(cov[0][0] * cov[1][1]);
    std::cout << "Correlation coefficient: " << rho << std::endl;

}


void Transformance(double scale1, double scale2, double angle, double vec1[2], double vec2[2],
                   std::vector<double> &X, std::vector<double> &Y,
                   std::vector<double> &X_result, std::vector<double> &Y_result){

    // transformacja obrotu i skalowania

    vec1[0] = scale1 * std::cos(angle);
    vec1[1] = scale1 * std::sin(angle);
    
    vec2[0] = -scale2 * std::sin(angle);
    vec2[1] = scale1 * std::cos(angle);

    for(int i = 0; i < n; ++i){
        X_result[i] = vec1[0] * X[i] + vec2[0] * Y[i];
        Y_result[i] = vec1[1] * X[i] + vec2[1] * Y[i];
    }

}
