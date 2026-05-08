#include "config.hpp"
#include "params.hpp"
#include "SingleBarrier.hpp"
#include "RTDSolver.hpp"
#include "QPCSolver.hpp"
#include <omp.h>

int main()
{
    stream_config(std::cout, 9, false);

    ex_counter(); // 1
    {
        SingleBarrier sb(5.0, 0.27, 0.063, calc_m_AlGaAs(0.3), 500);
        sb.scan_and_save(0.0, 1.0, 0.0001, "zad1_single_barrier.dat");
    }

    ex_counter(); // 2
    {
        SingleBarrier sb(5.0, 0.27, 0.063, calc_m_AlGaAs(0.3), 500);
        sb.scan_var_mass_and_save(0.0, 1.0, 0.0001, "zad2_single_barrier_var_mass.dat");
    }

    ex_counter(); // 3a
    {
        RTDSolver rtd(5.0, 3.0, 0.27, 0.063, calc_m_AlGaAs(0.3), 530);
        rtd.scan_TE_and_save(0.0, 1.0, 0.0001, "zad3_rtd_TE.dat");
    }

    ex_counter(); // 3b
    {
        RTDSolver rtd(5.0, 3.0, 0.27, 0.063, calc_m_AlGaAs(0.3), 530);
        rtd.compute_IV_and_save(0.0, 0.5, 0.001, 0.087, "zad3_rtd_IV.dat");
    }

    ex_counter(); // 4a
    {
        QPCSolver qpc;
        qpc.save_En("zad4_En.dat");
    }

    ex_counter(); // 4b
    {
        QPCSolver qpc;
        qpc.compute_GE_and_save(0.0, 0.2, 0.001, "zad4_GE.dat");
    }

    ex_counter(); // 4c
    {
        QPCSolver qpc;
        qpc.compute_GVg_and_save(0.0, 25.0, 0.01, "zad4_GVg.dat");
    }

    std::cout << std::endl;
    return 0;
}
