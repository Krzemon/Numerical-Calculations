#!/usr/bin/env python3
import sys
import matplotlib.pyplot as plt
from pathlib import Path

def main():
    if len(sys.argv) != 3:
        print("Uzycie: python3 plot_zad2.py <plik_danych> <plik_wykresu_u>")
        sys.exit(1)

    data_file = Path(sys.argv[1])
    plot_file = Path(sys.argv[2])

    if not data_file.exists():
        print(f"Blad: plik danych '{data_file}' nie istnieje.")
        sys.exit(1)

    try:
        with open(data_file) as f:
            content = f.read().strip()
    except Exception as e:
        print(f"Blad przy wczytywaniu pliku: {e}")
        sys.exit(1)

    blocks = content.split("\n\n")

    N_list, x_list, u_d_list, u_n_list = [], [], [], []

    for block in blocks:
        lines = [ln.strip() for ln in block.splitlines() if ln.strip()]
        try:
            header = lines[0]
            N_val = int(header.split('=')[1].strip())
        except Exception as e:
            print(f"Blad przy odczycie N w bloku:\n{block}\n{e}")
            continue
        try:
            x = list(map(float, lines[1].split()))
            u_d = list(map(float, lines[2].split()))
            u_n = list(map(float, lines[3].split()))
        except Exception as e:
            print(f"Blad przy konwersji danych:\n{block}\n{e}")
            continue

        N_list.append(N_val)
        x_list.append(x)
        u_d_list.append(u_d)
        u_n_list.append(u_n)

    plt.figure(figsize=(10, 6))
    for idx, N in enumerate(N_list):
        plt.plot(x_list[idx], u_n_list[idx], linestyle='--', label=f'u_n(x), N = {N}', linewidth=4)
    plt.plot(x_list[idx], u_d_list[idx], linestyle='-', label=f'u_d(x)', color='black', linewidth=2)
    plt.title('Metoda kolokacji — rozwiązanie numeryczne u_n(x) vs dokładne u_d(x)')
    plt.xlabel('x')
    plt.ylabel('u(x)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(plot_file, dpi=200)
    plt.close()
    print(f'Zapisano wykres szukanej funkcji do: "{plot_file.resolve()}"')

if __name__ == "__main__":
    main()