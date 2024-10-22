#include "wiener.h"

const double D = 1.0;
const int N_max = 10000;
const double dt = 0.1;
const double t_max = 100.0;
const double sigma_dt = std::sqrt(2 * D * dt);

double x[N_max] = {0};
double y[N_max] = {0};

void saveCoefficients(double D_xx, double D_xy, double D_yy, int iteration) {
    std::ofstream file;
    file.open("../data/coefficients.dat", std::ios::app);
    file << iteration * dt << " " << D_xx << " " << D_xy << " " << D_yy << std::endl;
    file.close();
}

void saveData(const std::string& filename, const double x1[], const double y1[]) {
    std::ofstream file;
    file.open("../data/" + filename + ".dat", std::ios::out | std::ios::trunc); // Otwarcie pliku do nadpisania
    for (int i = 0; i < N_max; ++i) {
        file << x[i] << " " << y[i] << std::endl;
    }
    file.close();
}

void SimStep() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0, sigma_dt);

    for (int i = 0; i < N_max; ++i) {
        x[i] += d(gen);
        y[i] += d(gen);
    }
}

double mean(const double x[], int n){
    double sum = 0;
    for(int i = 0; i < n; ++i){
        sum += x[i];
    }
    return sum / n;
}

double mean_square(const double x[], int n){
    double sum = 0;
    for(int i = 0; i < n; ++i){
        sum += x[i] * x[i];
    }
    return sum / n;
}

void wiener() {
    int n_t_max = t_max / dt;

    double D_xx[n_t_max-1], D_xy[n_t_max-1], D_yy[n_t_max-1];

    std::ofstream file("../data/coefficients.dat");
    file.close();

    for (int i = 1; i < n_t_max; ++i) {
        SimStep();
        double t = i * dt;
        double sum_x = 0, sum_y = 0, sum_xx = 0, sum_xy = 0, sum_yy = 0;

        if(t == 0.1)
            saveData("coordinate_0-1", x, y);
        else if(t == 1.0)
            saveData("coordinate_1-0", x, y);
        else if(t == 5.0)
            saveData("coordinate_5-0", x, y);


        for (int j = 0; j < N_max; ++j) {
            sum_x += x[j];
            sum_y += y[j];
            sum_xx += x[j] * x[j];
            sum_xy += x[j] * y[j];
            sum_yy += y[j] * y[j];
        }
    
        sum_x /= N_max;
        sum_y /= N_max;
        sum_xx /= N_max;
        sum_xy /= N_max;
        sum_yy /= N_max;

        D_xx[i-1] = (sum_xx - sum_x * sum_x / N_max) / (2 * t);
        D_xy[i-1] = (sum_xy - sum_x * sum_y / N_max) / (2 * t);
        D_yy[i-1] = (sum_yy - sum_y * sum_y / N_max) / (2 * t);

        saveCoefficients(D_xx[i-1], D_xy[i-1], D_yy[i-1], i);
    }


    double D_xx_avg = mean(D_xx, n_t_max);
    double D_xy_avg = mean(D_xy, n_t_max);
    double D_yy_avg = mean(D_yy, n_t_max);
    double D2_xx_avg = mean_square(D_xx, n_t_max);
    double D2_xy_avg = mean_square(D_xy, n_t_max);
    double D2_yy_avg = mean_square(D_yy, n_t_max);
    double sigma_D_xx = sqrt((D2_xx_avg - D_xx_avg * D_xx_avg) / n_t_max);
    double sigma_D_xy = sqrt((D2_xy_avg - D_xy_avg * D_xy_avg) / n_t_max);
    double sigma_D_yy = sqrt((D2_yy_avg - D_yy_avg * D_yy_avg) / n_t_max);

    std::cout << "D_xx = " << D_xx_avg << " +/- " << sigma_D_xx << std::endl;
    std::cout << "D_xy = " << D_xy_avg << " +/- " << sigma_D_xy << std::endl;
    std::cout << "D_yy = " << D_yy_avg << " +/- " << sigma_D_yy << std::endl;

}