#include "rozklady_2D.h"

int main() {

    // zadanie 1
    std::vector<double> X(n);
    std::vector<double> Y(n);

    Box_Muller(X, Y);

    std::ofstream Box_Muller_File("../data/Box_Muller.dat", std::ios::out);

    for (int i = 0; i < n; i++){
        Box_Muller_File << X[i] << " " << Y[i] << std::endl;
    }

    Box_Muller_File.close();

    // zadanie 2
    std::vector<double> X_norm(n);
    std::vector<double> Y_norm(n);
    for(int i = 0; i < n; ++i){
        X_norm[i] = X[i] / std::sqrt(X[i] * X[i] + Y[i] * Y[i]);
        Y_norm[i] = Y[i] / std::sqrt(X[i] * X[i] + Y[i] * Y[i]);
    }
    // on circle     

    std::ofstream Circle_On("../data/Circle_On.dat", std::ios::out);

    for (int i = 0; i < n; i++){
        Circle_On << X_norm[i] << " " << Y_norm[i] << std::endl;
    }

    Circle_On.close();

    // in circle   
    Circle_Norm(X_norm, Y_norm);

    std::ofstream Circle_Norm_File("../data/Circle_Norm.dat", std::ios::out);

    for (int i = 0; i < n; i++){
        Circle_Norm_File << X_norm[i] << " " << Y_norm[i] << std::endl;
    }

    Circle_Norm_File.close();

    // zadanie 3
    double alpha = M_PI / 4;
    double b1 = 1;
    double b2 = 0.2;

    std::vector<double> X_elipsa(n);
    std::vector<double> Y_elipsa(n);

    double r1[2];
    double r2[2];

    Transformance(b1, b2, alpha, r1, r2, X_norm, Y_norm, X_elipsa, Y_elipsa);

    std::ofstream Ellipse_File("../data/Ellipse.dat", std::ios::out);

    for (int i = 0; i < n; i++){
        Ellipse_File << X_elipsa[i] << " " << Y_elipsa[i] << std::endl;
    }
    
    Ellipse_File.close();

    double cov[2][2] = {0.0};
    Variance(X_elipsa, Y_elipsa, cov);

    Print_Cov(cov);

    std::cout << "\n ";
    std::cout << "-------------------";
    std::cout << "\n " << std::endl;

    // zadanie 4

    std::vector<double> X_norm_2(n);
    std::vector<double> Y_norm_2(n);

    std::vector<double> X_norm_rotate(n);
    std::vector<double> Y_norm_rotate(n);

    double rn1[2];
    double rn2[2];

    Box_Muller(X_norm_2, Y_norm_2);

    Transformance(b1, b2, alpha, rn1, rn2, X_norm_2, Y_norm_2, X_norm_rotate, Y_norm_rotate);

    std::ofstream NormRotateFile("../data/NormRotate.dat", std::ios::out);

    for (int i = 0; i < n; i++){
        NormRotateFile << X_norm_rotate[i] << " " << Y_norm_rotate[i] << std::endl;
    }

    NormRotateFile.close();


    double norm_cov[2][2] = {0.0};

    Variance(X_norm_rotate, Y_norm_rotate, norm_cov);

    Print_Cov(norm_cov);


    return 0;
}