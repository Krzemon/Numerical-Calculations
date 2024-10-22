#include "LJ_Potential.h"
#include "Particle.h"
#include "Obstacle.h"
#include "Point.h"

int main() {
    double Lx = 10.0;
    double Ly = 10.0;
    double Lx_2 = 2 * Lx;
    double Ly_2 = 2 * Ly;

    int npartsx = 20;
    int npartsy = 20;
    int npartsx_2 = 2 * npartsx;
    int nparts = npartsx * npartsy + 4 * (npartsx_2 - 1);

    double Xh = 0.2;
    double Xc = 0.1;
    double ixh = Xh * 2 * Lx;
    double ixc = (1 - Xc) * 2 * Lx;
    double FX = 40.0;
    double Fc = 0.95;

    double m = 1.0;
    double mp = 1.0e20;
    int nsteps = 100;
    double dt = 0.01;

    double sigma = 2.0 * Lx / (npartsx - 1) ;
    double epsilon = 5.0;
    double rcut = 2.5 * sigma;

    // Inicjalizacja ljPotential
    LJ_Potential &ljPotential = LJ_Potential::getInstance(sigma, epsilon, rcut);

    std::vector<Particle> particles;
    double x_min = -10.0;
    double y_min = -10.0;
    double y_max = 10.0;
    double x_max = 10.0;
    x_min += 1e-4;
    y_min += 1e-4;
    x_max -= 1e-4;
    y_max -= 1e-4;

/******************** Obsługa punktów brzegu rury *************************/

    for (int ix = 1; ix < npartsx_2 - 1; ++ix) {
        double x = x_min + sigma * ix / 2;
        particles.push_back(Particle(x, y_min, mp, 1));
        particles.push_back(Particle(x + sigma / 4, y_min + sigma / 5, mp, 1));
    }
    particles.push_back(Particle(x_max, y_min, mp, 1));
    particles.push_back(Particle(x_max + sigma / 4, y_min + sigma / 5, mp, 1));

    for (int ix = 1; ix < npartsx_2 - 1; ++ix) {
        double x = x_min + sigma * ix / 2;
        particles.push_back(Particle(x, y_max, mp, 1));
        particles.push_back(Particle(x - sigma / 4, y_max - sigma / 5, mp, 1));
    }
    particles.push_back(Particle(x_max, y_max, mp, 1));
    particles.push_back(Particle(x_max - sigma / 4, y_max - sigma / 5, mp, 1));

    int vec_size = static_cast<int>(particles.size());
    std::cout << "Liczba punktów rury: "<<vec_size << std::endl;
/******************** Inicjalizacja przeszkody *************************/ 
    // Create vector of particles
    std::vector<Particle> obstacleParticles;

    // Define the obstacle shape
    std::vector<Point> obstacleShape = {
        Point(-2, -3),
        Point(-2, 2),
        Point(-3, 0),
        Point(-1, 0)
    };

    // Populate obstacleParticles using obstacleShape
    for (const auto& obstaclePoint : obstacleShape) {
        obstacleParticles.push_back(Particle(obstaclePoint.getX(), obstaclePoint.getY(), mp, 1));
    }

    // Convert vector of particles to vector of points
    std::vector<Point> obstaclePoints;
    for (const auto& particle : obstacleParticles) {
        obstaclePoints.push_back(Point(particle.xpos, particle.ypos));
    }

    // Create obstacle
    Obstacle obstacle(obstaclePoints, mp, sigma, epsilon, rcut);


/******************** Obsługa punktów cząstek przepływających dla przeszkody *************************/ 
    int flow_particles = 0;
    for (int iy = 0; iy < npartsy; ++iy) {
        for (int ix = 0; ix < npartsx; ++ix) {
            double x = x_min + (static_cast<double>(ix) / (npartsx - 2)) * (x_max - x_min);
            double y = y_min + (static_cast<double>(iy) / (npartsy - 2)) * (y_max - y_min);

            // Losowa prędkość z rozkładu normalnego
            double vx = 2 * randNormal();
            double vy = 2 * randNormal();

            if (!obstacle.isInsideObstacle(Point(x, y)) && !obstacle.isNearObstacle(Point(x, y), sigma / 5)) {
                particles.push_back(Particle(x, y, m, vx, vy));
            }
        }

        if (iy < npartsy - 2) {
            double vx = 2 * randNormal();
            double vy = 2 * randNormal();
            double x = x_min + sigma;
            double y = y_min + sigma + iy * sigma;

            // Sprawdź czy punkt znajduje się wewnątrz przeszkody lub jest zbyt blisko
            if (!obstacle.isInsideObstacle(Point(x, y)) && obstacle.isNearObstacle(Point(x, y), sigma / 5)) {
                particles.push_back(Particle(x, y, m, vx, vy));
            }
        }
    }
    int vec_size_1 = static_cast<int>(particles.size());
    int vec_size_2 = vec_size_1 - vec_size;     // moglbym uzyc tego ale gdzies mi uciekaja czastki z polozeniami ktore sa NAN,
                                                // nie ruszam w analizie bo mam zrzut pamieci
    for (int i = 0; i < nparts; ++i) {
        if (!std::isnan(particles[i].xpos) && !std::isnan(particles[i].ypos))
            flow_particles++;
    }
    std::cout << "Liczba cząstek przepływających: "<< flow_particles - vec_size<< std::endl; //odejmuje czastki rury
/******************** Usunięcie plików przed zapisem *************************/ 

    std::string positionFilePath = "data/position.txt";
    std::ofstream positionFile(positionFilePath, std::ios::trunc);  // Otwarcie pliku w trybie "trunc" (czyli usunięcie istniejącej zawartości)

    if (!positionFile.is_open()) {
        std::cerr << "Błąd podczas otwierania pliku position.txt do zapisu." << std::endl;
    } else {
        positionFile.close();
    }

    std::string velocityFilePath = "data/velocity.txt";
    std::ofstream velocityFile(velocityFilePath, std::ios::trunc);  // Otwarcie pliku w trybie "trunc" (czyli usunięcie istniejącej zawartości)

    if (!velocityFile.is_open()) {
        std::cerr << "Błąd podczas otwierania pliku velocity.txt do zapisu." << std::endl;
    } else {
        velocityFile.close();
    }

    std::string temperatureFilePath = "data/temperature.txt";
    std::ofstream temperatureFile(temperatureFilePath, std::ios::trunc);  // Otwarcie pliku w trybie "trunc" (czyli usunięcie istniejącej zawartości)

    if (!temperatureFile.is_open()) {
        std::cerr << "Błąd podczas otwierania pliku temperature.txt do zapisu." << std::endl;
    } else {
        temperatureFile.close();
    }

    std::string particle_velocityFilePath = "data/particle_velocity.txt";
    std::ofstream particle_velocity(particle_velocityFilePath, std::ios::trunc);  // Otwarcie pliku w trybie "trunc" (czyli usunięcie istniejącej zawartości)

    if (!particle_velocity.is_open()) {
        std::cerr << "Błąd podczas otwierania pliku particle_velocity.txt do zapisu." << std::endl;
    } else {
        particle_velocity.close();
    }


/******************** Pętla czasowa symulacji *************************/ 

    for (int n = 1; n <= nsteps; ++n) {
        // Resetowanie sił dla każdej cząstki
        for (int i = 0; i < nparts; ++i) {
            particles[i].xfor = 0.0;
            particles[i].yfor = 0.0;
        }

        // Obliczanie sił między cząstkami
        for (int i = 0; i < nparts; ++i) {
            for (int j = 0; j < nparts; ++j) {
                if (i != j) {
                    double dx = particles[i].xpos - particles[j].xpos;
                    double dy = particles[i].ypos - particles[j].ypos;
                    double d = sqrt(dx * dx + dy * dy);

                    if (d < rcut) {
                        double fx = ljPotential.calculateForce(d) * dx / d;
                        double fy = ljPotential.calculateForce(d) * dy / d;
                        particles[i].xfor += fx;
                        particles[j].xfor -= fx;
                        particles[i].yfor += fy;
                        particles[j].yfor -= fy;
                    }
                }
            }
        }

        // Dodawanie siły do cząstek przed barierą
        for (int i = 0; i < nparts; ++i) {
            if (particles[i].xpos < ixh && i > vec_size) {
                particles[i].xfor += FX;
            }
        }

        obstacle.reactToObstacle(particles);
        
        for (int i = 0; i < nparts; ++i) {
            particles[i].integrate(dt, Fc, ixh);

            // Warunki brzegowe
            if (particles[i].xpos > x_max) {
                particles[i].xpos -= Lx_2;
            }
            if (particles[i].ypos > y_max) {
                particles[i].ypos -= Ly_2;
            }
            if (particles[i].xpos < x_min) {
                particles[i].xpos += Lx_2;
            }
            if (particles[i].ypos < y_min) {
                particles[i].ypos += Ly_2;
            }
        }

        // Dodawanie danych położeń do pliku
        std::ofstream positionFile(positionFilePath, std::ios::app);
        if (positionFile.is_open()) {
            for (int i = 0; i < nparts; i++) {
                // Sprawdzanie czy cząstka jest wewnątrz przeszkody
                if (obstacle.isInsideObstacle(Point(particles[i].xpos, particles[i].ypos))) {
                    continue; // Pomijanie cząstek wewnątrz przeszkody
                }

                positionFile << n << "," << particles[i].xpos << "," << particles[i].ypos << "\n";
            }
            positionFile << "\n\n";
            positionFile.close();
        } else {
            std::cerr << "Błąd podczas zapisywania danych do pliku position.txt." << std::endl;
        }

        // Dodawanie danych prędkości do pliku // tylko cząstki przepływające, realnie płynie 307 cząstek
        std::ofstream velocityFile(velocityFilePath, std::ios::app);
        if (velocityFile.is_open()) {
            for (int i = 0; i < nparts; i++) {
                // Sprawdzanie czy cząstka jest wewnątrz przeszkody
                if (obstacle.isInsideObstacle(Point(particles[i].xpos, particles[i].ypos))) {
                    continue; // Pomijanie cząstek wewnątrz przeszkody
                }

                if (i > 156 && !std::isnan(particles[i].xpos) && !std::isnan(particles[i].ypos)) {
                    double velocity = particles[i].xvel * particles[i].xvel + particles[i].yvel * particles[i].yvel;
                    velocityFile << n << "," << particles[i].xvel << "," << particles[i].yvel << "," << velocity << "\n";
                }
            }
            velocityFile << "\n\n";
            velocityFile.close();
        } else {
            std::cerr << "Błąd podczas zapisywania danych do pliku velocity.txt." << std::endl;
        }

        int who = 200;    // zakres od 157 do 556, inaczej od vec_size do vec_size_1

        // Dodawanie danych prędkości jednej cząstki do pliku
        std::ofstream particle_velocity(particle_velocityFilePath, std::ios::app);
        if (particle_velocity.is_open() && who >= vec_size && who < vec_size_1) {
            // Sprawdzanie czy cząstka jest wewnątrz przeszkody
            if (!obstacle.isInsideObstacle(Point(particles[who].xpos, particles[who].ypos))) {
                particle_velocity << n << "," << particles[who].xvel << "," << particles[who].yvel << "\n";
            } else {
                std::cerr << "Cząstka znajduje się wewnątrz przeszkody. Pomijanie zapisu do pliku particle_velocity.txt." << std::endl;
            }
            particle_velocity.close();
        } else {
            std::cerr << "Błąd podczas zapisywania danych do pliku particle_velocity.txt." << std::endl;
        }

        double kineticEnergy = 0.0;
        double k_B = 1;             //1.38064852e-23;

        for (int i = vec_size; i < vec_size_1; ++i) {
            // Sprawdzanie czy cząstka jest wewnątrz przeszkody
            if (!obstacle.isInsideObstacle(Point(particles[i].xpos, particles[i].ypos))) {
                if (!std::isnan(particles[i].xvel) && !std::isnan(particles[i].yvel)) {
                    kineticEnergy += 0.5 * m * (particles[i].xvel * particles[i].xvel + particles[i].yvel * particles[i].yvel);
                }
            } else {
                std::cerr << "Cząstka znajduje się wewnątrz przeszkody. Pomijanie w obliczeniach kinetycznej energii." << std::endl;
            }
        }

        double temperature = (1.0 / (k_B * (2 * nparts - 2))) * kineticEnergy;

        // Dodawanie danych temperatury do pliku
        std::ofstream temperatureFile(temperatureFilePath, std::ios::app);
        if (temperatureFile.is_open()) {
            temperatureFile << n << "," << temperature << "\n";
            temperatureFile.close();
        } else {
            std::cerr << "Błąd podczas zapisywania danych do pliku temperature.txt." << std::endl;
        }
    }

    return 0;
}


