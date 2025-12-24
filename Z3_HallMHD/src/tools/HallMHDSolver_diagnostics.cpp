
#include "0_basic/MPI_WRAPPER.h"

#include "HallMHD_Solver.h"

void HallMHDSolver::PrintMinMaxDiagnostics_()
{
    const int myid = par_->GetInt("myid");

    double rho_min_l = std::numeric_limits<double>::infinity();
    double rho_max_l = -std::numeric_limits<double>::infinity();
    double p_min_l = std::numeric_limits<double>::infinity();
    double p_max_l = -std::numeric_limits<double>::infinity();
    double V_min_l = std::numeric_limits<double>::infinity();
    double V_max_l = -std::numeric_limits<double>::infinity();

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
        auto &U = fld_->field(fid_.fid_U, ib);
        auto &PV = fld_->field(fid_.fid_PV, ib);
        auto &Bcel = fld_->field(fid_.fid_Bcell, ib);
        auto &divB = fld_->field(fid_.fid_divB, ib);
        auto &jac = fld_->field(fid_.fid_Jac, ib);

        const Int3 lo = U.inner_lo();
        const Int3 hi = U.inner_hi();

        for (int k = lo.k; k < hi.k; ++k)
            for (int j = lo.j; j < hi.j; ++j)
                for (int i = lo.i; i < hi.i; ++i)
                {
                    const double rho = U(i, j, k, 0);
                    const double p = PV(i, j, k, 3);
                    const double u = PV(i, j, k, 0);
                    const double v = PV(i, j, k, 1);
                    const double w = PV(i, j, k, 2);

                    const double bx = Bcel(i, j, k, 0);
                    const double by = Bcel(i, j, k, 1);
                    const double bz = Bcel(i, j, k, 2);
                    const double bmag = std::sqrt(bx * bx + by * by + bz * bz);

                    const double divb = divB(i, j, k, 0) * jac(i, j, k, 0);

                    rho_min_l = std::min(rho_min_l, rho);
                    rho_max_l = std::max(rho_max_l, rho);
                    p_min_l = std::min(p_min_l, p);
                    p_max_l = std::max(p_max_l, p);
                    V_min_l = std::min(V_min_l, sqrt(u * u + v * v + w * w));
                    V_max_l = std::max(V_max_l, sqrt(u * u + v * v + w * w));

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
    double rho_min_g, rho_max_g, p_min_g, p_max_g, V_min_g, V_max_g;
    double bx_min_g, bx_max_g, by_min_g, by_max_g, bz_min_g, bz_max_g;
    double bmag_min_g, bmag_max_g;
    double divb_absmax_g;

    PARALLEL::mpi_barrier();

    double mins_l[7] = {rho_min_l, p_min_l, bx_min_l, by_min_l, bz_min_l, bmag_min_l, V_min_l};
    double mins_g[7];
    double maxs_l[8] = {rho_max_l, p_max_l, bx_max_l, by_max_l, bz_max_l, bmag_max_l, divb_absmax_l, V_max_l};
    double maxs_g[8];

    PARALLEL::mpi_min(mins_l, mins_g, 7);
    PARALLEL::mpi_max(maxs_l, maxs_g, 8);

    // 解包
    rho_min_g = mins_g[0];
    p_min_g = mins_g[1];
    bx_min_g = mins_g[2];
    by_min_g = mins_g[3];
    bz_min_g = mins_g[4];
    bmag_min_g = mins_g[5];
    V_min_g = mins_g[6];

    rho_max_g = maxs_g[0];
    p_max_g = maxs_g[1];
    bx_max_g = maxs_g[2];
    by_max_g = maxs_g[3];
    bz_max_g = maxs_g[4];
    bmag_max_g = maxs_g[5];
    divb_absmax_g = maxs_g[6];
    V_max_g = maxs_g[7];

    // --- print on rank 0 only ---
    if (myid == 0)
    {
        std::printf(
            "\t\trho   =[%.6e, %.6e]\tp            =[%.6e, %.6e]\tVel =[%.6e, %.6e]\n"
            "\t\tBmag  =[%.6e, %.6e]\tmax|divB_PHI|=[%26.6e]\n"
            "\t\tBx    =[%.6e, %.6e]\tBy           =[%.6e, %.6e]\tBz  =[%.6e, %.6e]\n",
            rho_min_g, rho_max_g, p_min_g, p_max_g, V_min_g, V_max_g,
            bmag_min_g, bmag_max_g, divb_absmax_g,
            bx_min_g, bx_max_g, by_min_g, by_max_g, bz_min_g, bz_max_g);
        std::cout << std::endl
                  << std::flush;
    }
}