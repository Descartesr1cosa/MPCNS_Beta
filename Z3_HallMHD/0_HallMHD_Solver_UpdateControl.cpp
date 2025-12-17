#include "HallMHD_Solver.h"

bool HallMHDSolver::UpdateControlAndOutput()
{
    // 更新控制
    control_.Update();

    if (fmod(control_.Nstep, control_.Resouput_step) == 0)
        PrintMinMaxDiagnostics_();

    if (control_.if_outfile)
        output_.output_field();
    // output_.output_plt_cell_field(output_.var_defaut_plt_name); //For debug
    // output_.output_plt_field(); // output_.output_plt_cell_field(output_.var_defaut_plt_name); // output_field();
    if (control_.if_stop)
    {
        // output_.output_plt_cell_field(output_.var_defaut_plt_name); //For debug
        output_.output_field();
        return true;
    }
    return false;
}

// 统计 inner 域的 min/max，并做 MPI 全局归约
void HallMHDSolver::PrintMinMaxDiagnostics_()
{
    const int myid = par_->GetInt("myid");

    double rho_min_l = std::numeric_limits<double>::infinity();
    double rho_max_l = -std::numeric_limits<double>::infinity();
    double p_min_l = std::numeric_limits<double>::infinity();
    double p_max_l = -std::numeric_limits<double>::infinity();

    double bx_min_l = std::numeric_limits<double>::infinity();
    double bx_max_l = -std::numeric_limits<double>::infinity();
    double by_min_l = std::numeric_limits<double>::infinity();
    double by_max_l = -std::numeric_limits<double>::infinity();
    double bz_min_l = std::numeric_limits<double>::infinity();
    double bz_max_l = -std::numeric_limits<double>::infinity();

    double bmag_min_l = std::numeric_limits<double>::infinity();
    double bmag_max_l = -std::numeric_limits<double>::infinity();

    double divb_absmax_l = 0.0;

    const int nblock = fld_->num_blocks();
    for (int ib = 0; ib < nblock; ++ib)
    {
        auto &U = fld_->field(fid_U, ib);
        auto &PV = fld_->field(fid_PV, ib);
        auto &Bcel = fld_->field(fid_Bcell, ib);
        auto &divB = fld_->field("divB", ib);

        const Int3 lo = U.inner_lo();
        const Int3 hi = U.inner_hi();

        for (int k = lo.k; k < hi.k; ++k)
            for (int j = lo.j; j < hi.j; ++j)
                for (int i = lo.i; i < hi.i; ++i)
                {
                    const double rho = U(i, j, k, 0);
                    const double p = PV(i, j, k, 3);

                    const double bx = Bcel(i, j, k, 0);
                    const double by = Bcel(i, j, k, 1);
                    const double bz = Bcel(i, j, k, 2);
                    const double bmag = std::sqrt(bx * bx + by * by + bz * bz);

                    const double divb = divB(i, j, k, 0);

                    rho_min_l = std::min(rho_min_l, rho);
                    rho_max_l = std::max(rho_max_l, rho);
                    p_min_l = std::min(p_min_l, p);
                    p_max_l = std::max(p_max_l, p);

                    bx_min_l = std::min(bx_min_l, bx);
                    bx_max_l = std::max(bx_max_l, bx);
                    by_min_l = std::min(by_min_l, by);
                    by_max_l = std::max(by_max_l, by);
                    bz_min_l = std::min(bz_min_l, bz);
                    bz_max_l = std::max(bz_max_l, bz);

                    bmag_min_l = std::min(bmag_min_l, bmag);
                    bmag_max_l = std::max(bmag_max_l, bmag);

                    divb_absmax_l = std::max(divb_absmax_l, std::abs(divb));
                }
    }

    // --- MPI global reduction ---
    double rho_min_g, rho_max_g, p_min_g, p_max_g;
    double bx_min_g, bx_max_g, by_min_g, by_max_g, bz_min_g, bz_max_g;
    double bmag_min_g, bmag_max_g;
    double divb_absmax_g;

    PARALLEL::mpi_barrier();
    MPI_Allreduce(&rho_min_l, &rho_min_g, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&rho_max_l, &rho_max_g, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);
    MPI_Allreduce(&p_min_l, &p_min_g, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
    MPI_Allreduce(&p_max_l, &p_max_g, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);

    MPI_Allreduce(&bx_min_l, &bx_min_g, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
    MPI_Allreduce(&bx_max_l, &bx_max_g, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);
    MPI_Allreduce(&by_min_l, &by_min_g, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
    MPI_Allreduce(&by_max_l, &by_max_g, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);
    MPI_Allreduce(&bz_min_l, &bz_min_g, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
    MPI_Allreduce(&bz_max_l, &bz_max_g, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);

    MPI_Allreduce(&bmag_min_l, &bmag_min_g, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
    MPI_Allreduce(&bmag_max_l, &bmag_max_g, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);

    MPI_Allreduce(&divb_absmax_l, &divb_absmax_g, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);

    // --- print on rank 0 only ---
    if (myid == 0)
    {
        std::printf(
            "\t\trho      =[%.6e, %.6e]\n"
            "\t\tp        =[%.6e, %.6e]\n"
            "\t\tBmag     =[%.6e, %.6e]\n"
            "\t\tBx       =[%.6e, %.6e]\tBy  =[%.6e, %.6e]\tBz  =[%.6e, %.6e]\n"
            "\t\tmax|divB|=%.6e\n\n",
            rho_min_g, rho_max_g, p_min_g, p_max_g,
            bmag_min_g, bmag_max_g,
            bx_min_g, bx_max_g, by_min_g, by_max_g, bz_min_g, bz_max_g,
            divb_absmax_g);
    }
}