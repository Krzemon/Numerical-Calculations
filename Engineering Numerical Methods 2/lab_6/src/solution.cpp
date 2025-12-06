#include "solution.hpp"
#include "config.hpp"

void calc_zad_1(std::ofstream &out)
{
    make_grid(20, -5, 5, 20, -5, 5);
    out << std::fixed << std::setprecision(6);
    out << "# index x y nx=" << g_nx << " ny=" << g_ny << "\n";
    for (const auto &node : g_nodes)
    {
        out << node.idx << " " << node.x << " " << node.y << "\n";
    }

    out << "\n# elementy i j k\n";
    for (std::size_t m = 0; m < g_elem_nodes.size(); ++m)
    {
        const NodeSet &e = g_elem_nodes[m];
        out << m << " "
            << e[0]->idx << " "
            << e[1]->idx << " "
            << e[2]->idx << "\n";
    }
}

void save_global_E(std::ofstream &out)
{
    out << std::fixed << std::setprecision(6);
    for (int i = 0; i < g_N; ++i)
    {
        for (int j = 0; j < g_N; ++j)
        {
            out << g_E[i * g_N + j] << " ";
        }
        out << "\n";
    }
}

void save_global_O(std::ofstream &out)
{
    out << std::fixed << std::setprecision(6);
    for (int i = 0; i < g_N; ++i)
    {
        for (int j = 0; j < g_N; ++j)
        {
            out << g_O[i * g_N + j] << " ";
        }
        out << "\n";
    }
}

void calc_zad_2(std::ofstream &E_out, std::ofstream &O_out, std::ofstream &evals_out)
{
    save_global_E(E_out);
    save_global_O(O_out);
    solve_generalized_eigen_lapack(g_E.data(), g_O.data(), g_N, evals, evecs);

    int K = std::min(10, g_N);

    evals_out << std::setprecision(15);
    for (int k = 0; k < K; k++)
    {
        double lambda = evals[k];
        if (lambda < 0 && lambda > -1e-12)
            lambda = 0.0;
        evals_out << lambda << "\n";
    }

    const int Nx = 201;
    const int Ny = 201;

    for (int k = 0; k < K; ++k)
    {
        std::string fname = "mode_" + std::to_string(k) + ".dat";
        auto fh = prepareDataFile(fname);
        std::ofstream &fout = fh.out;

        fout << std::setprecision(12);
        fout << "# lambda " << evals[k] << "\n";
        fout << "# x y phi\n";

        double dx = (g_xmax - g_xmin) / (Nx - 1);
        double dy = (g_ymax - g_ymin) / (Ny - 1);

        for (int iy = 0; iy < Ny; ++iy)
        {
            double y = g_ymin + iy * dy;
            for (int ix = 0; ix < Nx; ++ix)
            {
                double x = g_xmin + ix * dx;

                double uk = 0.0;
                bool found = false;

                // iterujemy po wszystkich elementach (dwa trojkaty na kwadrat)
                for (int m = 0; m < g_M; m += 2)
                {
                    double x_left = g_elem_nodes[m][0]->x;
                    double x_right = g_elem_nodes[m][1]->x;
                    double y_bottom = g_elem_nodes[m][0]->y;
                    double y_top = g_elem_nodes[m][2]->y;

                    if (x < x_left || x > x_right || y < y_bottom || y > y_top)
                        continue;

                    // wybór trójkąta w kwadracie
                    int elem_idx = ((x - x_left) + (y - y_bottom) > (x_right - x_left)) ? m + 1 : m;
                    const NodeSet &tri = g_elem_nodes[elem_idx];

                    // mały układ liniowy do współrzędnych lokalnych
                    double x0 = tri[0]->x, y0 = tri[0]->y;
                    double x1 = tri[1]->x, y1 = tri[1]->y;
                    double x2 = tri[2]->x, y2 = tri[2]->y;

                    double A[2][2] = {{x1 - x0, x2 - x0}, {y1 - y0, y2 - y0}};
                    double b[2] = {x - x0, y - y0};

                    // rozwiąż 2x2 układ: A * xi = b
                    double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
                    if (std::fabs(det) < 1e-14)
                        continue;
                    double xi1 = (b[0] * A[1][1] - b[1] * A[0][1]) / det;
                    double xi2 = (A[0][0] * b[1] - A[1][0] * b[0]) / det;

                    // funkcje kształtu dla trójkąta liniowego
                    double phi0 = 1.0 - xi1 - xi2;
                    double phi1 = xi1;
                    double phi2 = xi2;

                    uk = phi0 * evecs[tri[0]->idx][k] +
                         phi1 * evecs[tri[1]->idx][k] +
                         phi2 * evecs[tri[2]->idx][k];

                    found = true;
                    break;
                }

                if (!found)
                    uk = 0.0; // punkt poza elementami
                fout << x << " " << y << " " << uk << "\n";
            }
        }
        std::cout << "Zapisano: " << fh.path << "\n";
    }
}

void calc_zad_3()
{
    int targets[] = {1, 2}; // drugi i trzeci wektor wlasny
    c_2.resize(g_N);
    c_3.resize(g_N);

    for (int idx = 0; idx < 2; ++idx)
    {
        std::vector<double> *out = (idx == 0 ? &c_2 : &c_3);
        int k = targets[idx];
        double norm_sq = 0.0;
        for (int i = 0; i < g_N; ++i)
        {
            double sum = 0.0;
            for (int j = 0; j < g_N; ++j)
            {
                sum += evecs[j][k] * g_O[i * g_N + j];
            }
            norm_sq += evecs[i][k] * sum;
        }
        double inv_norm = 1.0 / std::sqrt(norm_sq);

        for (int i = 0; i < g_N; ++i)
        {
            (*out)[i] = evecs[i][k] * inv_norm;
        }
    }
}

void calc_zad_4()
{
    double omega3 = std::sqrt(evals[2]);
    g_y0.resize(g_N);
    g_v0.resize(g_N);
    for (int i = 0; i < g_N; ++i)
    {
        g_y0[i] = c_2[i];
        g_v0[i] = omega3 * c_3[i];
    }
}

void calc_zad_5(std::ofstream &out, bool gif)
{
    const double beta = 0.25;
    const double omega2 = std::sqrt(evals[1]);
    const double Tmax = 2.0 * M_PI / omega2;
    const double dt = Tmax / 10000.0;
    const int Nsteps = static_cast<int>(Tmax / dt);

    double snapshot_times[5] = {0.0, 0.25 * Tmax, 0.5 * Tmax, 0.75 * Tmax, Tmax};
    int snapshot_counter = 0;
    int gif_counter = 0;

    std::vector<double> O_local = g_O;
    std::vector<double> E_local = g_E;

    gsl_vector *y = gsl_vector_alloc(g_N);
    gsl_vector *v = gsl_vector_alloc(g_N);
    gsl_vector *y_new = gsl_vector_alloc(g_N);
    gsl_vector *v_new = gsl_vector_alloc(g_N);
    // gsl_vector *y_old = gsl_vector_alloc(g_N);
    gsl_vector *tmp_vec = gsl_vector_alloc(g_N);
    gsl_vector *tmp_vec2 = gsl_vector_alloc(g_N);
    gsl_vector *c2 = gsl_vector_alloc(g_N);
    gsl_vector *c3 = gsl_vector_alloc(g_N);

    for (int i = 0; i < g_N; ++i)
    {
        gsl_vector_set(y, i, g_y0[i]);
        gsl_vector_set(v, i, g_v0[i]);
        gsl_vector_set(c2, i, c_2[i]);
        gsl_vector_set(c3, i, c_3[i]);
    }
    // gsl_vector_memcpy(y_old, y);

    // A = O + beta*dt^2*E
    std::vector<double> A_row(g_N * g_N, 0.0);
    for (int i = 0; i < g_N; ++i)
        for (int j = 0; j < g_N; ++j)
            A_row[i * g_N + j] = O_local[i * g_N + j] + beta * dt * dt * E_local[i * g_N + j];

    std::vector<double> A_col(g_N * g_N);
    row_to_col_major(A_row.data(), A_col.data(), g_N);
    cholesky_factor(A_col.data(), g_N);

    // O_LLT do v
    std::vector<double> O_col(g_N * g_N);
    row_to_col_major(O_local.data(), O_col.data(), g_N);
    cholesky_factor(O_col.data(), g_N);

    out << "# t yk_O_c2 yk_O_c3 yk_O_yk yk_E_yk\n";
    out << std::fixed << std::setprecision(10);

    auto write_outputs = [&](double t, gsl_vector *y_vec)
    {
        gsl_matrix_view Oview = gsl_matrix_view_array(O_local.data(), g_N, g_N);
        gsl_matrix_view Eview = gsl_matrix_view_array(E_local.data(), g_N, g_N);

        double s1, s2, s3, s4;

        // y^T O c2
        gsl_blas_dgemv(CblasNoTrans, 1.0, &Oview.matrix, c2, 0.0, tmp_vec);
        gsl_blas_ddot(y_vec, tmp_vec, &s1);

        // y^T O c3
        gsl_blas_dgemv(CblasNoTrans, 1.0, &Oview.matrix, c3, 0.0, tmp_vec);
        gsl_blas_ddot(y_vec, tmp_vec, &s2);

        // y^T O y
        gsl_vector *temp = gsl_vector_calloc(g_N);
        gsl_blas_dgemv(CblasNoTrans, 1.0, &Oview.matrix, y_vec, 0.0, temp);
        gsl_blas_ddot(y_vec, temp, &s3);
        gsl_vector_free(temp);

        // y^T E y
        temp = gsl_vector_calloc(g_N);
        gsl_blas_dgemv(CblasNoTrans, 1.0, &Eview.matrix, y_vec, 0.0, temp);
        gsl_blas_ddot(y_vec, temp, &s4);
        gsl_vector_free(temp);

        out << t << " " << s1 << " " << s2 << " " << s3 << " " << s4 << "\n";
    };

    write_outputs(0.0, y);

    const double coef2 = ((2.0 * beta - 1.0) / 2.0) * dt * dt;
    gsl_matrix_view Oview = gsl_matrix_view_array(O_local.data(), g_N, g_N);
    gsl_matrix_view Eview = gsl_matrix_view_array(E_local.data(), g_N, g_N);

    for (int k = 0; k < Nsteps; ++k)
    {
        // RHS = O*y + dt*O*v + coef2*E*y
        gsl_blas_dgemv(CblasNoTrans, 1.0, &Oview.matrix, y, 0.0, tmp_vec);
        gsl_blas_dgemv(CblasNoTrans, 1.0, &Oview.matrix, v, 0.0, tmp_vec2);
        gsl_blas_daxpy(dt, tmp_vec2, tmp_vec);
        gsl_blas_dgemv(CblasNoTrans, 1.0, &Eview.matrix, y, 0.0, tmp_vec2);
        gsl_blas_daxpy(coef2, tmp_vec2, tmp_vec);

        std::vector<double> rhs(g_N);
        for (int i = 0; i < g_N; ++i)
            rhs[i] = gsl_vector_get(tmp_vec, i);
        solve_cholesky(A_col.data(), rhs.data(), g_N);
        for (int i = 0; i < g_N; ++i)
            gsl_vector_set(y_new, i, rhs[i]);

        // O*v_{k+1} = O*v_k - dt/2*(E*y_k + E*y_{k+1})
        gsl_blas_dgemv(CblasNoTrans, 1.0, &Eview.matrix, y, 0.0, tmp_vec);
        gsl_blas_dgemv(CblasNoTrans, 1.0, &Eview.matrix, y_new, 0.0, tmp_vec2);
        gsl_vector_add(tmp_vec, tmp_vec2);
        gsl_blas_dgemv(CblasNoTrans, 1.0, &Oview.matrix, v, 0.0, v_new);
        gsl_blas_daxpy(-dt / 2.0, tmp_vec, v_new);

        solve_cholesky(O_col.data(), v_new->data, g_N);


        // gsl_vector_memcpy(y_old, y);
        gsl_vector_memcpy(y, y_new);
        gsl_vector_memcpy(v, v_new);

        // // warunki brzegowe
        // for (int i = 0; i < g_N; ++i)
        // {
        //     const Node &node = g_nodes[i];
        //     if (std::abs(node.x - g_xmin) < 1e-10 || std::abs(node.x - g_xmax) < 1e-10 ||
        //         std::abs(node.y - g_ymin) < 1e-10 || std::abs(node.y - g_ymax) < 1e-10)
        //     {
        //         gsl_vector_set(y, i, g_y0[i]);
        //         gsl_vector_set(v, i, g_v0[i]);
        //     }
        // }

        if (snapshot_counter < 5 && ((k + 1) * dt >= snapshot_times[snapshot_counter] - 1e-10))
        {
            g_y_snapshots[snapshot_counter].resize(g_N);
            for (int i = 0; i < g_N; ++i)
                g_y_snapshots[snapshot_counter][i] = gsl_vector_get(y, i);
            snapshot_counter++;
        }

        if ((k + 1) % 100 == 0 || (k + 1 == Nsteps))
        {
            write_outputs((k + 1) * dt, y);

            if (gif)
            {
                gif_snapshots[gif_counter].resize(g_N);
                for (int i = 0; i < g_N; ++i)
                    gif_snapshots[gif_counter][i] = gsl_vector_get(y, i);
                gif_counter++;
            }
        }
    }

    gsl_vector_free(y);
    gsl_vector_free(v);
    gsl_vector_free(y_new);
    gsl_vector_free(v_new);
    gsl_vector_free(tmp_vec);
    gsl_vector_free(tmp_vec2);
    // gsl_vector_free(y_old);
    gsl_vector_free(c2);
    gsl_vector_free(c3);
}

void calc_zad_6(std::ofstream &out_1, std::ofstream &out_2, std::ofstream &out_3, std::ofstream &out_4, std::ofstream &out_5)
{
    const int Nx = 200;
    const int Ny = 200;
    const int Nsnap = 5;
    std::ofstream *outs[Nsnap] = {&out_1, &out_2, &out_3, &out_4, &out_5};

    for (int snap = 0; snap < Nsnap; ++snap)
    {
        std::ofstream &fout = *outs[snap];
        fout << std::fixed << std::setprecision(12);
        fout << "# klatka = " << snap + 1 << "/" << Nsnap << "\n";
        fout << "# x y u\n";

        double dx = (g_xmax - g_xmin) / (Nx - 1);
        double dy = (g_ymax - g_ymin) / (Ny - 1);

        for (int iy = 0; iy < Ny; ++iy)
        {
            double y = g_ymin + iy * dy;
            for (int ix = 0; ix < Nx; ++ix)
            {
                double x = g_xmin + ix * dx;

                double u = 0.0;
                bool found = false;

                for (int m = 0; m < g_M; m += 2)
                {
                    double x_left = g_elem_nodes[m][0]->x;
                    double x_right = g_elem_nodes[m][1]->x;
                    double y_bottom = g_elem_nodes[m][0]->y;
                    double y_top = g_elem_nodes[m][2]->y;

                    if (x < x_left || x > x_right || y < y_bottom || y > y_top)
                        continue;

                    int elem_idx = ((x - x_left) + (y - y_bottom) > (x_right - x_left)) ? m + 1 : m;
                    const NodeSet &tri = g_elem_nodes[elem_idx];

                    double x0 = tri[0]->x, y0 = tri[0]->y;
                    double x1 = tri[1]->x, y1 = tri[1]->y;
                    double x2 = tri[2]->x, y2 = tri[2]->y;

                    double A[2][2] = {{x1 - x0, x2 - x0}, {y1 - y0, y2 - y0}};
                    double b[2] = {x - x0, y - y0};

                    double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
                    if (std::fabs(det) < 1e-14)
                        continue;
                    double xi1 = (b[0] * A[1][1] - b[1] * A[0][1]) / det;
                    double xi2 = (A[0][0] * b[1] - A[1][0] * b[0]) / det;

                    double phi0 = 1.0 - xi1 - xi2;
                    double phi1 = xi1;
                    double phi2 = xi2;

                    u = phi0 * g_y_snapshots[snap][tri[0]->idx] +
                        phi1 * g_y_snapshots[snap][tri[1]->idx] +
                        phi2 * g_y_snapshots[snap][tri[2]->idx];

                    found = true;
                    break;
                }

                if (!found)
                    u = 0.0;
                fout << x << " " << y << " " << u << "\n";
            }
        }
        std::cout << "Zapisano mapę dla klatki " << snap + 1 << "/" << Nsnap << "\n";
    }
}

void write_gif_snapshots(std::ofstream &out)
{
    const int Nx = 200;
    const int Ny = 200;
    const int Nsnap = 100; // liczba snapshotów
    double dx = (g_xmax - g_xmin) / (Nx - 1);
    double dy = (g_ymax - g_ymin) / (Ny - 1);

    out << std::fixed << std::setprecision(12);

    for (int snap = 0; snap < Nsnap; ++snap)
    {
        out << "# klatka = " << snap + 1 << "/" << Nsnap << "\n";
        out << "# x y u\n";

        for (int iy = 0; iy < Ny; ++iy)
        {
            double y = g_ymin + iy * dy;
            for (int ix = 0; ix < Nx; ++ix)
            {
                double x = g_xmin + ix * dx;

                double u = 0.0;
                bool found = false;

                for (int m = 0; m < g_M; m += 2)
                {
                    double x_left = g_elem_nodes[m][0]->x;
                    double x_right = g_elem_nodes[m][1]->x;
                    double y_bottom = g_elem_nodes[m][0]->y;
                    double y_top = g_elem_nodes[m][2]->y;

                    if (x < x_left || x > x_right || y < y_bottom || y > y_top)
                        continue;

                    int elem_idx = ((x - x_left) + (y - y_bottom) > (x_right - x_left)) ? m + 1 : m;
                    const NodeSet &tri = g_elem_nodes[elem_idx];

                    double x0 = tri[0]->x, y0 = tri[0]->y;
                    double x1 = tri[1]->x, y1 = tri[1]->y;
                    double x2 = tri[2]->x, y2 = tri[2]->y;

                    double A[2][2] = {{x1 - x0, x2 - x0}, {y1 - y0, y2 - y0}};
                    double b[2] = {x - x0, y - y0};

                    double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
                    if (std::fabs(det) < 1e-14)
                        continue;

                    double xi1 = (b[0] * A[1][1] - b[1] * A[0][1]) / det;
                    double xi2 = (A[0][0] * b[1] - A[1][0] * b[0]) / det;

                    double phi0 = 1.0 - xi1 - xi2;
                    double phi1 = xi1;
                    double phi2 = xi2;

                    u = phi0 * gif_snapshots[snap][tri[0]->idx] +
                        phi1 * gif_snapshots[snap][tri[1]->idx] +
                        phi2 * gif_snapshots[snap][tri[2]->idx];

                    found = true;
                    break;
                }

                if (!found)
                    u = 0.0;
                out << x << " " << y << " " << u << "\n";
            }
        }
    }
    std::cout << "Zapisano dane do tworzenia gifu\n";
}