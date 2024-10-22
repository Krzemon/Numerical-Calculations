#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

size_t N = 301;
const double L = 100.0 / 0.05292;
const double m = 0.067;
const double h_bar = 1.0;
const double dx = L / (N - 1);

const double to_meV = 27211.6;
//const size_t num_energies = 1000;
const size_t num_energies = 5000;
const double E_min = 0.0;
const double E_max = 35.0;
const double dE = (E_max - E_min) / num_energies;
//const double dE_atom = dE / to_meV;

// Funkcja do obliczania wartości funkcji falowej w punkcie i+1
double Psi(const double psi_i, const double psi_im1, const double E, const double V) {
    return (-2 * m / (h_bar * h_bar) * (E - V) * dx * dx * psi_i + 2 * psi_i - psi_im1);
}

// Funkcja do obliczania normy z funkcji falowej
double Integral(const std::vector<double>& psi){
    double sum = 0.0;
    for (int j = 0; j < N; ++j) {
        sum += psi[j] * psi[j];
    }
    sum *= dx;
    return std::sqrt(sum);
}

// Funkcja do obliczania miejsc zerowych
std::vector<double> ZeroPointEnergy(const std::vector<std::vector<double>>& psi, double lewe[], double prawe[], int ile){
    std::vector<double> miejscaZerowe(ile, 0.0);

    for (int i=0; i<ile; ++i){
        while (prawe[i]-lewe[i]>1e-6){
            double srodek = (lewe[i]+prawe[i])/2;
            if (psi[srodek][N-1]*psi[lewe[i]][N-1]<0){
                prawe[i] = srodek;
            }
            else{
                lewe[i] = srodek;
            }
        }
        miejscaZerowe[i] = (lewe[i] + prawe[i]) / 2;
    }
    return miejscaZerowe;
}


// Funkcja do obliczenia psi_3 i C
void CalculatePsiAndC(std::vector<std::vector<std::vector<double>>>& psi_3, std::vector<std::vector<double>>& C,
                      std::vector<double>& W, const std::vector<double>& E_3, double E_min, double dE, double to_meV,
                      double W_min, double dW) {
    const int W_size = psi_3.size();
    const int num_energies = psi_3[0].size();

    for (int k = 0; k < W_size; ++k) {
        W[k] = W_min + k * dW;
        for (int j = 0; j < num_energies; ++j) {
            psi_3[k][j][1] = 1.0;
            for (int i = 1; i < N; ++i) {
                double V = 0;
                if (i == N / 2)
                    V = -W[k] / to_meV;
                psi_3[k][j][i+1] = Psi(psi_3[k][j][i], psi_3[k][j][i-1], (E_min + dE * num_energies) / to_meV, V);
            }
        }
    }

    for (int k = 0; k < W_size; ++k) {
        for (int j = 0; j < num_energies; ++j) {
            C[k][j] = Integral(psi_3[k][j]);
        }
    }

    // Normalizacja psi_3
    for (int k = 0; k < W_size; ++k) {
        for (int j = 0; j < num_energies; ++j) {
            for (int i = 0; i < N; ++i) {
                psi_3[k][j][i] = psi_3[k][j][i] / C[k][j];
            }
        }
    }

    // Zapisz dane do pliku
    std::ofstream dataFile("../data/psi_W_7.txt");
    for (int k = 0; k < W_size; ++k) {
        for (int j = 0; j < num_energies; ++j) {
            for (int i = 0; i <= N; ++i) {
                dataFile << psi_3[k][j][i] << " ";
            }
            dataFile << "\n";
        }
    }
    dataFile.close();
}



int main() {
    std::vector<std::vector<double>> psi(num_energies, std::vector<double>(N));
    const double V_0 = 0;
    for (int i = 0; i < num_energies; ++i) {
        psi[i][1] = 1.0;
    }

    for (int i = 0; i < num_energies; ++i) {
        for (int j = 1; j < N-1; ++j) {
            psi[i][j+1] = Psi(psi[i][j], psi[i][j-1], (E_min + dE * i)/to_meV, V_0);
        }
    }

    //ZADANIE 1
    std::ofstream psi_file("../data/psi_file.txt", std::ios::out);
    for (int i = 0; i < num_energies; ++i) {
        for (int j = 0; j < N; ++j) {
            psi_file << psi[i][j] << " ";
        }
        psi_file << std::endl;
    }
    psi_file.close();

    size_t zero_energie_index = 80;         //dla 1000 to jest punkt 16
    std::vector<double> psi_0(N);   //blisko 0
    std::vector<double> psi_p(N);   //+5%
    std::vector<double> psi_m(N);   //-5%
    psi_0[1] = 1.0;
    for (int j = 1; j < N-1; ++j) {
        psi_0[j+1] = Psi(psi_0[j], psi_0[j-1], (E_min + dE * zero_energie_index)/to_meV, V_0);
    }
    psi_p[1] = 1.0;
    for (int j = 1; j < N-1; ++j) {
        psi_p[j+1] = Psi(psi_p[j], psi_p[j-1], 1.05*(E_min + dE * zero_energie_index)/to_meV, V_0);
    }
    psi_m[1] = 1.0;
    for (int j = 1; j < N-1; ++j) {
        psi_m[j+1] = Psi(psi_m[j], psi_m[j-1], 0.95*(E_min + dE * zero_energie_index)/to_meV, V_0);
    }

    double norm_0 = Integral(psi_0);
    double norm_p = Integral(psi_p);
    double norm_m = Integral(psi_m);

    for (int i = 0; i < N; ++i) {
        psi_0[i] /= norm_0;
        psi_p[i] /= norm_p;
        psi_m[i] /= norm_m;
    }

    std::ofstream energy_0("../data/0_energy_0.txt", std::ios::out);
    for (int i = 0; i < N; ++i) {
        energy_0 << psi_0[i] << std::endl;
    }
    energy_0.close();

    std::ofstream energy_p("../data/0_energy_p.txt", std::ios::out);
    for (int i = 0; i < N; ++i) {
        energy_p << psi_p[i] << std::endl;
    }
    energy_p.close();

    std::ofstream energy_m("../data/0_energy_m.txt", std::ios::out);
    for (int i = 0; i < N; ++i) {
        energy_m << psi_m[i] << std::endl;
    }
    energy_m.close();


    //ZADANIE 2
    // double zero_l[] = {10, 50, 130, 240, 390, 570, 770}; //indeksy stanow wlasnych energii
    // double zero_p[] = {30, 80, 160, 270, 410, 590, 790}; // num_E = 1000
    
    double zero_l[] = {50, 250, 650, 1200, 1950, 2850, 3850}; //indeksy stanow wlasnych energii
    double zero_p[] = {150, 400, 800, 1350, 2050, 2950, 3950}; // num_E = 5000


    int ile = 7;
    for (int i = 0; i < ile; ++i)
    {
        zero_l[i] *= dE;
        zero_p[i] *= dE;
    }

    std::ofstream zero_point;
    if (N == 101) {
        zero_point.open("../zero_point_N_101.txt");
    } else if (N == 301) {
        zero_point.open("../zero_point_N_301.txt");
    } else {
        zero_point.open("../zero_point_N_inne.txt"); 
    }

    // Wypisanie zakresu energii do pliku
    zero_point << "Zakres energii: " << E_min << " - " << E_max << std::endl;
    std::cout  << "Zakres energii: " << E_min << " - " << E_max << std::endl;

    std::vector<double> miejscaZerowe = ZeroPointEnergy(psi, zero_l, zero_p, ile);

    zero_point << "Miejsca zerowe metoda bisekcji: " << std::endl;
    std::cout  << "Miejsca zerowe metoda bisekcji: " << std::endl;
    for (int i = 0; i < 7; ++i) {
        zero_point << miejscaZerowe[i] << ", ";
        std::cout  << miejscaZerowe[i] << ", ";    
    }
    zero_point << std::endl;
    std::cout  << std::endl;

    double zero_eq[7]{}; // Miejsca zerowe z równania
    for (int i = 0; i < 7; ++i) {
        zero_eq[i] = to_meV * (i+1)*(i+1)*M_PI*M_PI*h_bar*h_bar/(2*m*L*L);
    }
    zero_point << "Miejsca zerowe wyniki dokładne: " << std::endl;
    std::cout  << "Miejsca zerowe wyniki dokładne: " << std::endl;
    for (int i = 0; i < 7; ++i) { 
        zero_point << zero_eq[i] << ", ";
        std::cout  << zero_eq[i] << ", ";
    }
    zero_point << std::endl;
    std::cout  << std::endl;

    zero_point.close();



    //ZADANIE 3
    //if i == N/2 then V = -W
    //wyznaczyc 7 miejsc zerowych w funkcji W w przedziale 0-1eV 
    


    const int W_size = 1000;
    const double W_min = 0.0;
    const double W_max = 1000.0;
    const double dW = 1.00;
    
    std::vector<std::vector<std::vector<double>>> psi_3(W_size, std::vector<std::vector<double>>(num_energies, std::vector<double>(N + 1, 0.0)));
    std::vector<std::vector<double>> C(W_size, std::vector<double>(num_energies, 0.0));
    std::vector<double> W(W_size, 0.0);
    std::vector<double> E_3(num_energies, 0.0);

    //liczy i zapisuje do pliku psi_3
    CalculatePsiAndC(psi_3, C, W, E_3, E_min, dE, to_meV, W_min, dW);

    return 0;
}
