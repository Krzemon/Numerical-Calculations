#!/usr/bin/env python3
import sys
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from io import BytesIO

def main():
    if len(sys.argv) != 3:
        print("Użycie: python3 create_gif.py <plik_danych> <plik_gifa>")
        return

    data_file = sys.argv[1]
    gif_file = sys.argv[2]

    frames = []
    xs, ys, us = [], [], []
    frame_number = None

    with open(data_file) as f:
        for line in f:
            line = line.strip()

            # Nowa sekcja (klatka)
            if line.startswith("# klatka"):
                # Jeżeli poprzednia sekcja istnieje → generuj klatkę
                if xs:
                    frames.append(create_frame(xs, ys, us, frame_number))
                    xs, ys, us = [], [], []  # reset

                # Wyciągamy numer klatki
                parts = line.split("=")[1].split("/")
                frame_number = int(parts[0].strip())
                continue

            # Inna linia komentarza pomijamy
            if line.startswith("#") or len(line) == 0:
                continue

            # Wczytanie x y u
            parts = line.split()
            if len(parts) >= 3:
                x, y, u = map(float, parts[:3])
                xs.append(x)
                ys.append(y)
                us.append(u)

    # Ostatnia sekcja
    if xs:
        frames.append(create_frame(xs, ys, us, frame_number))

    if not frames:
        print("Brak danych do GIF-a!")
        return

    # Zapis GIF-a
    frames[0].save(
        gif_file,
        save_all=True,
        append_images=frames[1:],
        duration=100,
        loop=0
    )
    print("Utworzono GIF:", gif_file)


def create_frame(xs, ys, us, frame_number):
    """Tworzy klatkę z mapą rozwiązania u(x,y)."""
    plt.figure(figsize=(6, 5))
    sc = plt.scatter(xs, ys, c=us, s=8, cmap='cividis')
    plt.colorbar(sc, label="u(x, y)")
    plt.title(f"Klatka {frame_number}")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.xlim(min(xs), max(xs))
    plt.ylim(min(ys), max(ys))
    plt.gca().set_aspect("equal", adjustable="box")
    plt.tight_layout()

    buf = BytesIO()
    plt.savefig(buf, format='png', dpi=200)
    plt.close()

    buf.seek(0)
    return Image.open(buf)


if __name__ == "__main__":
    main()