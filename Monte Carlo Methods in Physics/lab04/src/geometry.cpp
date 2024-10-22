#include "geometry.h"

double uniform() {
    return static_cast<double>(std::rand()) / RAND_MAX;
}

void sferycznie_konturowany(std::vector<double> &X, std::vector<double> &Y, double R) {
    int n = X.size();
    for(int i = 0; i < n; ++i){
        double U1 = uniform();
        double U2 = uniform();
        // double R = std::sqrt(-2 * std::log(U1));
        X[i] = R * std::sin(2 * M_PI * U2);
        Y[i] = R * std::cos(2 * M_PI * U2);
    }
}

void circle_Norm(std::vector<double> &X, std::vector<double> &Y){
    int n = X.size();
    for(int i = 0; i < n; ++i){
        double U3 = uniform();
        X[i] = X[i] * std::sqrt(U3);
        Y[i] = Y[i] * std::sqrt(U3);
    }
}

void translation(std::vector<double> &X, std::vector<double> &Y, double x_A, double y_A){
    int n = X.size();
    for(int i = 0; i < n; ++i){
        X[i] += x_A;
        Y[i] += y_A;
    }
}

double calculateOverlapArea(const std::vector<double> &X, const std::vector<double> &Y, double R_alpha, double R_beta, double x_A, int n) {
    double overlapAreaSum = 0.0;
    double overlapAreaSquaredSum = 0.0;

    for (int i = 0; i < n; ++i) {
        double distance = std::sqrt((X[i] + x_A) * (X[i] + x_A) + Y[i] * Y[i]);
        double overlapArea = 0.0;

        if (distance <= R_alpha && distance <= R_beta) {
            overlapArea = M_PI * R_alpha * R_alpha;
        } else if (distance <= R_alpha && distance > R_beta) {
            double angle1 = std::acos(R_beta / distance);
            double angle2 = std::acos((distance * distance + R_beta * R_beta - R_alpha * R_alpha) / (2 * distance * R_beta));
            overlapArea = angle1 * R_beta * R_beta + angle2 * R_alpha * R_alpha - R_alpha * std::sin(angle2) * distance;
        } else if (distance > R_alpha && distance <= R_beta) {
            double angle1 = std::acos(R_alpha / distance);
            double angle2 = std::acos((distance * distance + R_alpha * R_alpha - R_beta * R_beta) / (2 * distance * R_alpha));
            overlapArea = angle1 * R_alpha * R_alpha + angle2 * R_beta * R_beta - R_beta * std::sin(angle2) * distance;
        }

        overlapAreaSum += overlapArea;
        overlapAreaSquaredSum += overlapArea * overlapArea;
    }

    double meanOverlapArea = overlapAreaSum / n;
    double secondMoment = overlapAreaSquaredSum / n;
    return meanOverlapArea;
}


void drawPoints(const std::vector<double> &X, const std::vector<double> &Y, const std::string &filename) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Cannot open output file: " << filename << std::endl;
        return;
    }

    for (size_t i = 0; i < X.size(); ++i) {
        outFile << X[i] << " " << Y[i] << std::endl;
    }

    outFile.close();
}


double calculateVariance(double mean, double secondMoment, int N) {
    return (secondMoment - mean * mean) / N;
}

double calculateStandardDeviation(double mean, double secondMoment, int N) {
    double variance = secondMoment - mean * mean;
    return std::sqrt(variance / N);
}

double calculateCircleArea(double R_alpha) {
    return M_PI * R_alpha * R_alpha;
}

