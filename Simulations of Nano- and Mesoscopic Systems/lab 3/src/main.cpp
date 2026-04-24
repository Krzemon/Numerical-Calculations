#include "config.hpp"
#include "params.hpp"
#include "DoubleQuantumDot.hpp"
#include "TimeEvolution.hpp"
#include <iomanip>
#include <omp.h>
#include <vector>
#include <mutex>
 
 
int main()
{
    stream_config(std::cout, 9, false);
    ex_counter(); // 1
    {
        DoubleQuantumDot dqd {};
        auto f4 = prepareDataFile("energies_4states.dat");
        f4.out << "# F[kV/cm]  E0[meV]  E1[meV]  E2[meV]  E3[meV]\n";
 
        const int NF = 81; // ilosc krokow przemiatania potencjału
        for (int k = 0; k < NF; ++k) {
            double F_kVcm = -2.0 + 4.0 * k / (NF - 1); // [-2,2]
            double F_au = F_kVcm_to_au(F_kVcm);
 
            Eigen::MatrixXd eigvecs;
            Eigen::VectorXd vals = dqd.solveEigen(F_au, eigvecs);

            stream_config(f4.out);    
            f4.out << F_kVcm
                   << "  " << au_to_meV(vals(0))
                   << "  " << au_to_meV(vals(1))
                   << "  " << au_to_meV(vals(2))
                   << "  " << au_to_meV(vals(3))
                   << "\n";
        }
        dqd << std::cout;
        std::cout << " [saved]: " << f4.path << '\n';
    }
 
    stream_config(std::cout, 9, false);
    ex_counter(); // 2
    {
        DoubleQuantumDot dqd {};
        Eigen::MatrixXd eigvecs;
        Eigen::VectorXd eigvals = dqd.solveEigen(0.0, eigvecs);
 
        Eigen::VectorXd psi0 = eigvecs.col(0);
        Eigen::VectorXd psi1 = eigvecs.col(1);
 
        double dE_au = eigvals(1) - eigvals(0);
        std::cout << "dE = " << au_to_meV(dE_au) << " meV\n";
 
        auto f = prepareDataFile("wavefunctions.dat");
        f.out << "# x[nm]  psi0[a.u.]  psi1[a.u.]\n";
        stream_config(f.out);
        const auto& x = dqd.x(); // węzły
        
        for (int i = 0; i < dqd.get_params().n; ++i)
            f.out << au_to_nm(x(i)) << "  " << psi0(i) << "  " << psi1(i) << "\n";
        std::cout << " [saved]: " << f.path << '\n';
    }

    stream_config(std::cout, 9, false);
    ex_counter(); // 3
    {
        DoubleQuantumDot dqd {};
        EvolutionParams ep;
        Eigen::MatrixXd eigvecs;
        Eigen::VectorXd eigvals = dqd.solveEigen(0.0, eigvecs);
        Eigen::VectorXd psi0 = eigvecs.col(0);
        Eigen::VectorXd psi1 = eigvecs.col(1);
        double dE_au  = eigvals(1) - eigvals(0);
        double dE_meV = au_to_meV(dE_au);
 
        // a: omega = dE (rezonansowa)
        {
            ep.omega_meV = dE_meV;
            TimeEvolution te(dqd, psi0, psi1, ep);
            te.run();
            // std::cout << te;
            te << std::cout;
 
            auto f = prepareDataFile("evolution_resonance.dat");
            te.write(f.out);
            std::cout << " [saved]: " << f.path << '\n';
        }
        // b: omega = 0.95 * dE
        {
            ep.omega_meV = 0.95 * dE_meV;
            TimeEvolution te(dqd, psi0, psi1, ep);
            te.run();
            // std::cout << te;
            te << std::cout;
 
            auto f = prepareDataFile("evolution_95percent.dat");
            te.write(f.out);
            std::cout << " [saved]: " << f.path << '\n';
        }
    }
 
    stream_config(std::cout, 9, false);
    ex_counter(); // 4
    {
        DoubleQuantumDot dqd {};

        auto fmap = prepareDataFile("map_p1_t_F.dat");
        fmap.out << "# F[kV/cm]  t[ns]  p1\n";

        const int NF = 101;  // F w [0, 0.2] kV/cm, gesta siatka

        Eigen::MatrixXd eigvecs0;
        Eigen::VectorXd eigvals0 = dqd.solveEigen(0.0, eigvecs0);
        const double    dE_meV_global = au_to_meV(eigvals0(1) - eigvals0(0));
        const Eigen::VectorXd psi0_ref = eigvecs0.col(0);
        const Eigen::VectorXd psi1_ref = eigvecs0.col(1);

        // Wyniki z kazdego watku trzymamy osobno, potem scalamy w kolejnosci
        std::vector<std::vector<std::string>> results(NF);

        #pragma omp parallel for schedule(dynamic)
        for (int k = 0; k < NF; ++k)
        {
            EvolutionParams ep;
            ep.F_kVcm    = 0.2 * k / (NF - 1);
            ep.omega_meV = dE_meV_global;

            TimeEvolution te(dqd, psi0_ref, psi1_ref, ep);
            te.run();

            const auto& res = te.get_result();
            std::vector<std::string>& lines = results[k];
            for (std::size_t i = 0; i < res.t.size(); ++i) {
                lines.push_back(
                    std::to_string(ep.F_kVcm) + "  " +
                    std::to_string(au_to_ns(res.t[i])) + "  " +
                    std::to_string(res.p1[i]) + "\n"
                );
            }

            #pragma omp critical
            std::cerr << "  [map] F=" << ep.F_kVcm << " kV/cm done (thread " << omp_get_thread_num() << ")\n";
        }

        // Zapis w kolejnosci k=0..NF-1
        stream_config(fmap.out);
        for (int k = 0; k < NF; ++k) {
            for (const auto& line : results[k])
                fmap.out << line;
            fmap.out << "\n";
        }

        std::cout << " [saved]: " << fmap.path << '\n';
    }


    std::cout << std::endl;
    return 0;
}