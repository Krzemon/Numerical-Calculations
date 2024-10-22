#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cstdlib>
#include <functional>

#include <bits/stdc++.h> // tanh(x)

double uniform();

// liczba iteracji, lewy przedział, prawy przedział
void integral(std::function<double(double x)> fun, int K[4], double a, double b);