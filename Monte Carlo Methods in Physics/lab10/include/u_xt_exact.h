#ifndef U_XT_EXACT_H
#define U_XT_EXACT_H

#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include<vector>
#include<complex>
using namespace std;

/*****************************************************************************
 *  liczymy rozwiazanie dokladne dla impulsu generowanego przez zrodlo
 *  oznaczenia:
 * 			x -polozenie, t - aktualny czas, t0 - maksimum sygnalu zrodla, Omega=2*Pi*f - (f to czestotliwosc zrodla)
 *                sigma - szerokosc impulsu,  R - opor linii, G - konduktancja linii, C-pojemnosc linii, Rg - opor zrodla
 *                Rl-opor cewki w x=l,length-dlugosc linii, number_nodes - liczba wezlow w calkowaniu po czestosci (omedze),
 *			n_sum_terms - maksymalna liczba wyrazow uwzgledniana w sumowaniu wkladow ui (i=0,1,...,n_sum_terms)
 *
 *
 *****************************************************************************/


double u_xt_exact(double x, double t, double t0,  double freq, double sigma, double R, double G, double L, double C, double Rg, double Rl, double length, int number_nodes, int n_sum_terms);

#endif // U_XT_EXACT_H

