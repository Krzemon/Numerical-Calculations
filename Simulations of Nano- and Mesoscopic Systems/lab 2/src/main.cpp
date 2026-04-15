#include "config.hpp"
#include "ImaginaryTime.hpp"
#include "HartreeFock.hpp"

int main()
{
    stream_config();

    ex_counter();
    {
        ImaginaryTime im_time;
        im_time.run();
        std::cout << im_time << std::endl;

        auto f = prepareDataFile("energies_k.dat");
        im_time.write_energies(f.out);
        std::cout << " [saved]: " << f.path << '\n';
    }
 
    ex_counter();
    {
        auto fIT = prepareDataFile("E_vs_a_IT.dat");
        fIT.out << "# a [nm]  E_IT [eV]\n";
 
        for (double a = 30.0; a <= 60.0; a += 5.0) {
            ImaginaryTime it(a);
            it.run();
            fIT.out << a << '\t' << it.energy() << '\n';
            std::cout << it;
        }
        std::cout << " [saved]: " << fIT.path << '\n';
    }
 
    ex_counter();
    {
        {
            ImaginaryTime it30(30.0);
            it30.run();
            auto f = prepareDataFile("density_a30.dat");
            it30.write_density(f.out);
            std::cout << " [saved]: " << f.path << '\n';
        }
        std::cout << "--\n";
        {
            ImaginaryTime it60(60.0);
            it60.run();
            auto f = prepareDataFile("density_a60.dat");
            it60.write_density(f.out);
            std::cout << " [saved]: " << f.path << '\n';
        }
    }
 
    ex_counter();
    {
        auto fHF = prepareDataFile("E_vs_a_HF.dat");
        fHF.out << "# a [nm]  E_HF [eV]\n";
 
        for (double a = 30.0; a <= 60.0; a += 5.0) {
            HartreeFock hf(a);
            hf.run();
            fHF.out << a << '\t' << hf.energy() << '\n';
        }
        std::cout << " [saved]: " << fHF.path << '\n';
    }

    std::cout << std::endl;
    return 0;
}
