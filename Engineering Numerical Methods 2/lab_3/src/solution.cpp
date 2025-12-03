#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <filesystem>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "functions.hpp"
#include "config.hpp"

void calc_zad_1(std::ofstream& out);
void calc_zad_2(std::ofstream& out, int M);
void calc_zad_3(std::ofstream& out, int M);

void general_eigen(const std::vector<std::vector<double>>& A,
                   const std::vector<std::vector<double>>& B,
                   int N,
                   std::vector<double>& evals,
                   std::vector<std::vector<double>>& evecs)
{
    gsl_matrix* a = gsl_matrix_alloc(N,N);
    gsl_matrix* b = gsl_matrix_alloc(N,N);
    gsl_vector* eval = gsl_vector_alloc(N);
    gsl_matrix* evec = gsl_matrix_alloc(N,N);

    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++){
            gsl_matrix_set(a,i,j,A[i][j]);
            gsl_matrix_set(b,i,j,B[i][j]);
        }

    gsl_eigen_gensymmv_workspace* w = gsl_eigen_gensymmv_alloc(N);
    gsl_eigen_gensymmv(a,b,eval,evec,w);
    gsl_eigen_gensymmv_free(w);
    gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_VAL_ASC);

    evals.assign(N,0.0);
    evecs.assign(N,std::vector<double>(N,0.0));
    for(int k=0;k<N;k++){
        evals[k] = gsl_vector_get(eval,k);
        for(int i=0;i<N;i++)
            evecs[k][i] = gsl_matrix_get(evec,i,k);
    }

    gsl_matrix_free(a);
    gsl_matrix_free(b);
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
}

void calc_zad_1(std::ofstream& out) {
    int M = 5;
    auto xnodes = generate_nodes(M,1.0);
    std::vector<std::vector<double>> S,O;
    assemble_global_matrices(M,xnodes,S,O);

    int N = 2*M+1;
    // zapis macierzy S do pliku
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++) out << std::setprecision(12) << S[i][j] << " ";
        out << "\n";
    }
}

void calc_zad_2(std::ofstream& out, int M, const std::string& psi_prefix) {
    int N = 2*M+1;
    double a0=0.4, a1=2.0, da=0.05;

    std::vector<double> xnodes;
    std::vector<std::vector<double>> S,O;
    std::vector<double> evals;
    std::vector<std::vector<double>> evecs;

    // ----------------- wartości własne -----------------
    out << "# alpha";
    int max_ev_print=5; // tylko mu = 0..4
    for(int mu=0; mu<max_ev_print; mu++) out << " E" << mu;
    out << "\n";

    for(double alpha=a0; alpha<=a1+1e-12; alpha+=da){
        xnodes = generate_nodes(M, alpha);
        assemble_global_matrices(M, xnodes, S, O);
        general_eigen(S,O,N,evals,evecs);

        out << std::setprecision(8) << alpha;
        int how = std::min((int)evals.size(), max_ev_print);
        for(int m=0;m<how;m++) out << " " << std::setprecision(12) << evals[m];
        out << "\n";
    }

    // ----------------- funkcje własne dla alpha=1.4 -----------------
    double alpha_special = 1.4;
    xnodes = generate_nodes(M, alpha_special);
    assemble_global_matrices(M, xnodes, S, O);
    general_eigen(S,O,N,evals,evecs);

    // zapis funkcji własnych do pliku w gęstej siatce
    fs::path u_path;
    std::ofstream fu;
    std::string filename = psi_prefix + std::to_string(M) + "_alpha_1_4.dat";
    prepareDataFile(u_path, fu, filename);
    fu.close();

    write_u_of_x_file(u_path, M, xnodes, evecs, max_ev_print);
}

void calc_zad_3(std::ofstream& out, int M, const std::string& psi_prefix) {
    calc_zad_2(out, M, psi_prefix);
}