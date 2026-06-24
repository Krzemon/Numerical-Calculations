"""
Przykład nanodrutu w KWANT - obliczenie konduktancji G(E)
Na podstawie laboratorium: Pakiet KWANT - transport elektronowy
"""

import kwant
import numpy as np
import matplotlib.pyplot as plt

# ── Parametry układu ──────────────────────────────────────────────────────────
a   = 5e-9          # krok siatki [m], dx = 5 nm
t   = 1.0           # energia przeskoku (jednostki umowne; t = ℏ²/2m*dx²)
W   = 10            # szerokość nanodrutu [węzły] → 2W = 50 nm
L   = 30            # długość nanodrutu [węzły]

# ── Budowa układu (scattering region) ────────────────────────────────────────
lat = kwant.lattice.square(a, norbs=1)
syst = kwant.Builder()

# Wnętrze nanodrutu: energia nawęzłowa = 4t (człon kinetyczny)
syst[(lat(i, j) for i in range(L) for j in range(W))] = 4 * t

# Przeskoki w kierunku x i y: -t
syst[lat.neighbors()] = -t

# ── Kontakty (leads) ──────────────────────────────────────────────────────────
sym_left  = kwant.TranslationalSymmetry((-a, 0))
sym_right = kwant.TranslationalSymmetry(( a, 0))

lead_left  = kwant.Builder(sym_left)
lead_right = kwant.Builder(sym_right)

for lead, sym in [(lead_left, sym_left), (lead_right, sym_right)]:
    lead[(lat(0, j) for j in range(W))] = 4 * t
    lead[lat.neighbors()] = -t

syst.attach_lead(lead_left)
syst.attach_lead(lead_right)

# ── Finalizacja i obliczenia ───────────────────────────────────────────────────
syst = syst.finalized()

energies = np.linspace(0.05, 3.95, 200)   # zakres energii [w jednostkach t]
conductance = []

for E in energies:
    smat = kwant.smatrix(syst, E)
    conductance.append(smat.transmission(1, 0))   # G w jednostkach e²/h

# ── Wykres ────────────────────────────────────────────────────────────────────
plt.figure(figsize=(8, 4))
plt.plot(energies, conductance, color='steelblue', lw=1.5)
plt.xlabel('Energia [t]')
plt.ylabel('Konduktancja [e²/h]')
plt.title(f'Konduktancja nanodrutu (W={W}·a, L={L}·a)')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('conductance.png', dpi=150)
plt.show()
print("Zapisano wykres → conductance.png")